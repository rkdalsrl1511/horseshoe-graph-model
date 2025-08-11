
# ghs_like_mcmc.R
# Fully Bayesian Graphical Horseshoe-Like (GHS-LIKE) MCMC sampler in R
# --------------------------------------------------------------------
# This implementation follows Algorithm 2 ("GHS-LIKE-MCMC") and the posterior
# updates for tau^2 and xi given in eq. (3.5) of the paper:
# - Column-wise Gibbs updates for (beta, gamma) using blockwise formulas
# - Local scales t^2 and auxiliary m updates using the slash-normal / Pareto mixture
# - Global scale tau^2 with Makalic–Schmidt half-Cauchy augmentation
#
# INPUTS
#   S        : p x p sample covariance (i.e., X'X / n)
#   n        : sample size (so that nS = n * S)
#   burnin   : number of burn-in iterations
#   nmc      : number of saved posterior samples
#   thin     : thinning factor (save every 'thin' iterations after burn-in)
#   seed     : (optional) RNG seed
#   save_omega: if TRUE, returns 3D array of Omegas; if FALSE, returns last draw and running mean
#   t_cap    : numeric stability cap for t^2 entries (default 1e15, entries above are left unchanged)
#
# RETURNS
#   A list with posterior draws / summaries, including last Omega, last Sigma, tau2 path,
#   and either sampled Omegas (if save_omega=TRUE) or running mean of Omega.
#
# NOTE
#   Shapes/rates use R's rgamma parameterization: shape, rate.
#   Inv-Gamma(α, β) means 1 / Gamma(shape=α, rate=β).
#
# --------------------------------------------------------------------

rinvgamma <- function(n, shape, rate) {
  1 / rgamma(n, shape = shape, rate = rate)
}

ghs_like_mcmc <- function(Y, n, burnin, nmc, thin = 1, seed = NULL,
                          save_omega = FALSE, t_cap = 1e15, verbose = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  S <- t(Y) %*% Y
  p <- ncol(S)

  # Initialize
  Omega <- diag(p)
  Sigma <- diag(p)  # Omega^{-1}
  T2    <- matrix(1.0, p, p)  # stores t^2, symmetric, zero diagonal
  Maux  <- matrix(0.5, p, p)  # stores m in (0,1)
  diag(T2) <- 0; diag(Maux) <- 0

  tau2 <- 1.0
  xi   <- 1.0

  n_total <- burnin + nmc * thin
  save_idx <- 0L
  Omegas <- array(NA_real_, dim = c(p, p, nmc))
  tau2_path <- numeric(nmc)

  if (verbose) cat(sprintf("Starting GHS-LIKE-MCMC: p=%d, burnin=%d, saved=%d, thin=%d\n", p, burnin, nmc, thin))

  for (iter in 1:n_total) {
    # Sweep over columns i = 1..p
    for (i in 1:p) {
      # 1) Sample gamma ~ Gamma(n/2 + 1, 2/s_ii)
      sii <- S[i, i]
      S_col <- S[-i, i, drop = FALSE]
      gamma_i <- rgamma(1, shape = n/2 + 1, rate = max(sii, .Machine$double.eps) / 2)

      # 2) Compute Omega^{-1}_{(-i)(-i)} via Sigma blocks
      Sigma_11 <- Sigma[-i, -i, drop = FALSE]
      sigma_i  <- Sigma[-i, i, drop = FALSE]
      sigma_ii <- Sigma[i, i]
      OmegaInv_11 <- Sigma_11 - (sigma_i %*% t(sigma_i)) / sigma_ii

      # 3) Construct C = [ s_ii * OmegaInv_11 + (diag(tau^2 / t2))^{-1} ]^{-1}
      t2_vec <- as.numeric(T2[-i, i])
      # numerical cap condition for stability (those above cap won't be updated)
      upd_mask <- (t2_vec <= t_cap)
      Dinv <- diag(t2_vec / tau2, nrow = p - 1, ncol = p - 1)
      A <- sii * OmegaInv_11 + Dinv
      chol_A <- chol(A)

      # 4) Sample beta ~ N( -C S_{(-i)i}, C )
      mu_beta <- solve(t(chol_A), S_col)
      mu_beta <- -solve(chol_A, mu_beta)
      z <- rnorm(p - 1)
      beta <- as.numeric(mu_beta + solve(chol_A, z))

      # 5) Variable transform: omega_{(-i)i} = beta; omega_{ii} = gamma + beta' OmegaInv_11 beta
      Omega[-i, i] <- beta
      Omega[i, -i] <- beta
      Omega[i, i]  <- as.numeric(gamma_i + t(beta) %*% OmegaInv_11 %*% beta)

      # 6) Local scales: sample m in (0,1), then t^2
      rate_m <- t2_vec / 2
      u <- runif(length(rate_m))
      one_minus_e <- 1 - exp(-pmin(rate_m, 700))
      m_new <- -log(1 - u * one_minus_e) / pmax(rate_m, .Machine$double.eps)
      m_new <- pmin(m_new, 1 - 1e-12)
      
      # rate for t^2: m/2 + omega^2/(2 tau^2)
      rate_t <- m_new / 2 + (beta^2) / (2 * tau2)
      t2_prop <- rgamma(p - 1, shape = 1.5, rate = rate_t)

      # apply stability mask: do not update entries with previous t2 > t_cap
      t2_vec[upd_mask] <- t2_prop[upd_mask]
      T2[-i, i] <- t2_vec
      T2[i, -i] <- t2_vec
      Maux[-i, i] <- m_new
      Maux[i, -i] <- m_new

      # 7) Update Sigma blocks (Sherman–Morrison-style formulas)
      tmp <- OmegaInv_11 %*% beta
      Sigma[-i, -i] <- OmegaInv_11 + (tmp %*% t(tmp)) / gamma_i
      Sigma[-i, i]  <- - tmp / gamma_i
      Sigma[i, -i]  <- - tmp / gamma_i
      Sigma[i, i]   <- 1 / gamma_i
    } # end for i

    # 8) Sample global tau^2 and xi (Makalic–Schmidt)
    #    tau^2 | xi, ... ~ InvGamma( (p(p-1)/2 + 1)/2,  1/xi + sum_{i<j} t2_ij * omega_ij^2 / 2 )
    #    xi | tau ~ InvGamma( 1 + 1/tau^2, 1 )
    upper_idx <- upper.tri(T2, diag = FALSE)
    sum_t2w2 <- sum(T2[upper_idx] * (Omega[upper_idx]^2))
    shape_tau <- (p * (p - 1) / 2 + 1) / 2
    rate_tau  <- (1 / xi) + 0.5 * sum_t2w2
    tau2 <- rinvgamma(1, shape = shape_tau, rate = rate_tau)

    xi <- rinvgamma(1, shape = 1 + 1 / tau2, rate = 1)

    # Save after burn-in with thinning
    if (iter > burnin && ((iter - burnin) %% thin == 0)) {
      save_idx <- save_idx + 1L
      tau2_path[save_idx] <- tau2
      Omegas[, , save_idx] <- Omega
      
      if (verbose && (save_idx %% max(1, floor(nmc / 10)) == 0)) {
        cat(sprintf("Saved %d / %d posterior draws...\n", save_idx, nmc))
      }
    }
  } # end for iter

  out <- list(
    Omega_last = Omega,
    Sigma_last = Sigma,
    tau2_path  = tau2_path,
    tau2_last  = tau2,
    xi_last    = xi,
    T2_last    = T2,
    M_last     = Maux,
    burnin     = burnin,
    nmc        = nmc,
    thin       = thin
  )
  
  out$Omegas <- Omegas
  class(out) <- "ghs_like_mcmc"
  out
}

# Simple helper to compute posterior mean precision matrix if save_omega=FALSE
posterior_mean_precision <- function(fit) {
  if (!inherits(fit, "ghs_like_mcmc")) stop("fit must be object from ghs_like_mcmc()")
  if (!is.null(fit$Omega_mean)) return(fit$Omega_mean)
  apply(fit$Omegas, c(1, 2), mean)
}

# Omegas: M x p x p (저장된 사후 표본)  ← ghs_like_mcmc(..., save_omega=TRUE)로 받아야 합니다.
edges_from_ci <- function(Omegas, level = 0.50) {
  stopifnot(length(dim(Omegas)) == 3)
  M <- dim(Omegas)[3]; p <- dim(Omegas)[2]
  alpha <- 1 - level
  qlo <- alpha/2; qhi <- 1 - alpha/2
  
  # 사후 분위수 행렬 (p x p)
  lower <- apply(Omegas, 1:2, quantile, probs = qlo, names = FALSE)
  upper <- apply(Omegas, 1:2, quantile, probs = qhi, names = FALSE)
  
  adj <- (lower > 0) | (upper < 0)  # 0을 배제하면 엣지
  diag(adj) <- FALSE
  adj
}
