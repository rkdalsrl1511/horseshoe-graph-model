concentration_sssl <- function(Y, burn = 1000, iter = 5000, pi, lambda = 1, 
                               v0 = 0.02, h = 50) {
  
  N <- nrow(Y)
  p <- ncol(Y)
  S <- t(Y) %*% Y
  v1 <- h * v0
  V <- matrix(v1, nrow = p, ncol = p)
  diag(V) <- 0
  estimate_Omega <- matrix(0, nrow = p, ncol = p)
  estimate_Sigma <- S / N
  Omega_hat <- matrix(0, nrow = p, ncol = p)
  Sigma_hat <- matrix(0, nrow = p, ncol = p)
  concentraion_graph <- matrix(0, nrow = p, ncol = p)
  OmegaSamples <- array(0, dim = c(p, p, iter))
  SigmaSamples <- array(0, dim = c(p, p, iter))
  VSamples <- array(0, dim = c(p, p, iter))
  
  # sampling start
  for (i in 1:(iter + burn)) {
    if (i %% 1000 == 0) {
      cat("iter : ", i, "\n")
    }
    # update Omega
    for (j in 1:p) {
      # inverse omega11 matrix
      sig11 <- estimate_Sigma[-j, -j]
      sig12 <- estimate_Sigma[-j, j]
      sig22 <- estimate_Sigma[j, j]
      v12 <- V[-j, j]^2
      inv_omega11 <- sig11 - sig12 %*% t(sig12) / sig22
      
      # u sampling
      inv_C <- (S[j, j] + lambda) * inv_omega11 + diag(v12^(-1))
      C_chol <- chol(inv_C)
      x <- solve(t(C_chol), S[-j, j])
      mu <- -solve(C_chol, x)
      f <- rnorm(p-1, mean = 0, sd = 1)
      f_star <- solve(C_chol, f)
      u <- mu + f_star
      
      # v sampling
      v <- rgamma(1, shape = N/2 + 1, rate = (S[j, j] + lambda)/2)
      
      # omega update
      estimate_Omega[-j, j] <- u
      estimate_Omega[j, -j] <- u
      estimate_Omega[j, j] <- v + t(u) %*% inv_omega11 %*% u
      
      # sigma update
      z <- inv_omega11 %*% u
      estimate_Sigma[-j, -j] <- inv_omega11 + z %*% t(z) / v
      estimate_Sigma[-j, j] <- -z/v
      estimate_Sigma[j, -j] <- -z/v
      estimate_Sigma[j, j] <- 1/v
      
      # V matrix update
      w1 <- -log(v0) - (0.5*u^2)/(v0^2) + log(1 - pi)
      w2 <- -log(v1) - (0.5*u^2)/(v1^2) + log(pi)
      w_max <- ifelse(w1 >= w2, w1, w2)
      w <- exp(w2 - w_max) / (exp(w1 - w_max) + exp(w2 - w_max))
      v_star <- ifelse(runif(p-1, 0, 1) < w, v1, v0)
      V[-j, j] <- v_star
      V[j, -j] <- v_star
    }
    if (i > burn) {
      OmegaSamples[, , i-burn] <- estimate_Omega
      SigmaSamples[, , i-burn] <- estimate_Sigma
      VSamples[, , i-burn] <- V
      Omega_hat <- Omega_hat + estimate_Omega
      Sigma_hat <- Sigma_hat + estimate_Sigma
      concentraion_graph <- concentraion_graph + V
      
    }
  }
  
  model_info <- list(Y = Y, burn = burn, iter = iter, pi = pi, lambda = lambda, 
                     v0 = v0, v1 = v1)
  
  Omega_hat <- Omega_hat/iter
  Sigma_hat <- Sigma_hat/iter
  concentraion_graph <- ifelse(concentraion_graph/iter > (v0 + v1)/2, 1, 0)
  
  result <- list(info = model_info, OmegaSamples = OmegaSamples, 
                 SigmaSamples = SigmaSamples, VSamples = VSamples, 
                 OmegaHat = Omega_hat, SigmaHat = Sigma_hat,
                 Graph = concentraion_graph)
  
}