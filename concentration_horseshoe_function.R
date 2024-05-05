# horseshoe concentration model
concentration_horseshoe <- function(Y, burn = 1000, iter = 5000, lambda = 1,
                                    approximate = FALSE, fixed_threshold = 0,
                                    auto.threshold = FALSE, t = 5, 
                                    adapt_p0 = 0, adapt_p1 = -4.6*10^(-4)) {
  N <- nrow(Y)
  p <- ncol(Y)
  S <- t(Y) %*% Y
  xi <- 1
  eta <- matrix(1, nrow = p, ncol = p)
  diag(eta) <- 0
  estimate_Omega <- diag(1, nrow = p, ncol = p)
  estimate_Sigma <- S / N
  m_eff <- rep(p-1, p)
  if (fixed_threshold == 0) {
    fixed_threshold <- 1/(p^2)
  }
  Omega_hat <- matrix(0, nrow = p, ncol = p)
  Sigma_hat <- matrix(0, nrow = p, ncol = p)
  Eta_hat <- matrix(0, nrow = p, ncol = p)
  concentration_graph <- matrix(0, nrow = p, ncol = p)
  OmegaSamples <- array(0, dim = c(p, p, iter + burn))
  SigmaSamples <- array(0, dim = c(p, p, iter + burn))
  XiSamples <- rep(0, iter + burn)
  EtaSamples <- array(0, dim = c(p, p, iter + burn))
  ActiveMatrixSamples <- array(0, dim = c(p, p, iter + burn))
  
  # sampling start
  for (i in 1:(iter + burn)) {
    #if (i %% 100 == 0) {
    #  cat("iter : ", i, "\n")
    #}
    # active_matrix renewal
    active_matrix <- matrix(0, nrow = p, ncol = p)
    # update Omega
    for (j in 1:p) {
      v12 <- xi * eta[, j]
      if (approximate == TRUE) {
        if (auto.threshold == TRUE) {
          if (i %% t == 0) {
            u_i <- stats::runif(1, 0, 1)
            p_i <- exp(adapt_p0 + adapt_p1 * i)
            if (u_i < p_i) {
              m_eff[j] <- sum(1/((S[j, j]+lambda)^(-1) * xi * eta[-j, j] *
                                   diag(estimate_Omega[-j, -j]) + 1))
            }
          }
          threshold <- sort(eta[-j, j])[ceiling(m_eff[j])]
          active_matrix[-j, j] <- ifelse(eta[-j, j] <= threshold, 1, active_matrix[-j, j])
          active_matrix[j, -j] <- active_matrix[-j, j]
          k <- which(active_matrix[, j] == 1)
          h <- which(active_matrix[-j, j] == 1)
        } else {
          active_matrix[-j, j] <- ifelse(v12[-j] < 1/fixed_threshold, 1, 
                                         active_matrix[-j, j])
          active_matrix[j, -j] <- active_matrix[-j, j]
          k <- which(active_matrix[, j] == 1)
          h <- which(active_matrix[-j, j] == 1)
        }
      } else {
        k <- c(1:p)[-j]
        h <- 1:(p-1)
      }
      active_num <- length(k)
      sig11 <- estimate_Sigma[-j, -j]
      sig12 <- estimate_Sigma[-j, j]
      sig22 <- estimate_Sigma[j, j]
      inv_omega11 <- sig11 - sig12 %*% t(sig12) / sig22
      # u sampling
      if(active_num != 0) {
        inv_C <- (S[j, j] + lambda) * inv_omega11[h, h] + 
          diag(v12[k], nrow = active_num, ncol = active_num)
        inv_C <- (inv_C + t(inv_C))/2
        C_chol <- chol(inv_C)
        ########################################################################
        # 에러 검출용 코드 #####################################################
        ########################################################################
        #C_chol <- tryCatch(chol(inv_C),
        #                   error = function(e) {
        #                     list(now_iter = i,
        #                          now_j = j,
        #                          now_meff = m_eff[j],
        #                          now_omega = estimate_Omega,
        #                          now_active = active_matrix[, j],
        #                          now_active_matrix = active_matrix,
        #                          all_omega = OmegaSamples,
        #                          all_active_matrix = ActiveMatrixSamples,
        #                          all_xi = XiSamples,
        #                          all_eta = EtaSamples)
        #                   })
        #if (is.list(C_chol))
        #  return(C_chol)
        ########################################################################
        ########################################################################
        x <- solve(t(C_chol), S[k, j])
        mu <- -solve(C_chol, x)
        f <- rnorm(active_num, mean = 0, sd = 1)
        f_star <- solve(C_chol, f)
        u <- mu + f_star
        estimate_Omega[k, j] <- u
        estimate_Omega[j, k] <- u
        # inactive u sampling
        if (p-active_num-1 != 0) {
          #n_f <- rnorm(p-active_num-1, 0, 1)
          #n_u <- -(S[-c(k, j), j]/v12[-c(k, j)]) + n_f/sqrt(v12[-c(k, j)])
          #estimate_Omega[-c(k, j), j] <- n_u
          #estimate_Omega[j, -c(k, j)] <- n_u
          estimate_Omega[-c(k, j), j] <- 0.0001
          estimate_Omega[j, -c(k, j)] <- 0.0001
        }
      } else { 
        # 모든 성분이 inactive인 경우
        #f <- rnorm(p-1, 0, 1)
        #u <- -(S[-j, j]/v12[-j]) + f/sqrt(v12[-j])
        #estimate_Omega[-j, j] <- u
        #estimate_Omega[j, -j] <- u
        estimate_Omega[-j, j] <- 0.0001
        estimate_Omega[j, -j] <- 0.0001
      }
      u_star <- estimate_Omega[-j, j]
      # v sampling
      v <- rgamma(1, shape = N/2 + 1, rate = (S[j, j] + lambda)/2)
      estimate_Omega[j, j] <- v + t(u_star) %*% inv_omega11 %*% u_star
      # sigma update
      z <- inv_omega11 %*% u_star
      estimate_Sigma[-j, -j] <- inv_omega11 + z %*% t(z) / v
      estimate_Sigma[-j, j] <- -z/v
      estimate_Sigma[j, -j] <- -z/v
      estimate_Sigma[j, j] <- 1/v
      # eta update
      upsi <- stats::runif(p-1, 0, 1/(1 + eta[-j, j]))
      tempps <- (xi * u_star^2)/2
      ub <- (1 - upsi)/upsi
      Fub <- 1 - exp(-tempps * ub)
      Fub[Fub < (1e-04)] <- 1e-04
      up <- stats::runif(p-1, 0, Fub)
      new_eta <- -log(1 - up)/tempps
      new_eta <- ifelse(new_eta <= 2.220446e-16, 2.220446e-16, new_eta)
      eta[-j, j] <- new_eta
      eta[j, -j] <- new_eta
      
      # eta update : rejection sampler version
      #eps <- (xi * u_star^2)/2
      #new_eta <- rejection_sampler(eps, a = 0.2, b = 10)
      #new_eta <- ifelse(new_eta <= 2.220446e-16, 2.220446e-16, new_eta)
      #eta[-j, j] <- new_eta
      #eta[j, -j] <- new_eta
    }
    
    # xi update : halfCauchy
    tempt <- sum(estimate_Omega^2 * eta)/4
    utau <- stats::runif(1, 0, 1/(1 + xi))
    ubt <- (1 - utau)/utau
    Fubt <- stats::pgamma(ubt, (p*(p-1)+2)/4, scale = 1/tempt)
    Fubt <- max(Fubt, 1e-08)
    ut <- stats::runif(1, 0, Fubt)
    xi <- stats::qgamma(ut, (p*(p-1)+2)/4, scale = 1/tempt)
    
    #xi update : fixed
    #xi <- 100
    
    # xi update : Metropolis algorithm
    #tempt <- sum(estimate_Omega^2 * eta)/4
    #new_xi <- exp(stats::rnorm(1, mean = log(xi), sd = 0.8))
    #dxi <- stats::dgamma(c(xi, new_xi), (p*(p-1)+2)/4, scale = 1/tempt)
    #dxi <- ifelse(dxi <= 2.220446e-16, 2.220446e-16, dxi)
    #accept_prob <- ((1+xi)*new_xi*dxi[2])/((1+new_xi)*xi*dxi[1])
    #if (runif(1, 0, 1) < accept_prob) {
    #  xi <- new_xi
    #}
    
    OmegaSamples[, , i] <- estimate_Omega
    SigmaSamples[, , i] <- estimate_Sigma
    EtaSamples[, , i] <- eta
    XiSamples[i] <- xi
    ActiveMatrixSamples[, , i] <- active_matrix
    if (i > burn) {
      Omega_hat <- Omega_hat + estimate_Omega
      Sigma_hat <- Sigma_hat + estimate_Sigma
      Eta_hat <- Eta_hat + eta
      concentration_graph <- concentration_graph + active_matrix
    }
  }
  
  model_info <- list(Y = Y, burn = burn, iter = iter, pi = pi, lambda = lambda)
  Omega_hat <- Omega_hat/iter
  Sigma_hat <- Sigma_hat/iter
  Eta_hat <- Eta_hat/iter
  Xi_hat <- mean(XiSamples[(burn+1):(iter+burn)])
  concentration_graph <- ifelse(concentration_graph/iter > 0.5, 1, 0)
  
  result <- list(info = model_info, OmegaSamples = OmegaSamples, 
                 SigmaSamples = SigmaSamples, EtaSamples = EtaSamples, 
                 XiSamples = XiSamples, 
                 ActiveMatrixSamples = ActiveMatrixSamples, 
                 OmegaHat = Omega_hat, SigmaHat = Sigma_hat, EtaHat = Eta_hat, 
                 XiHat = Xi_hat, Graph = concentration_graph)
}