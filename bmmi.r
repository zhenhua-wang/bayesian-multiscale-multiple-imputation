library(Matrix)
library(MASS)
library(dlm)

bmmi <- function(num_iter, y, y_agg, miss, miss_agg,
                 a, R, tau, kappa, alpha, beta) {
  op <- options(digits = 7)
  on.exit(options(op))
  k <- dim(y)[1]
  T <- dim(y)[2]
  num_years <- T / 4
  y_rep <- NULL
  y_agg_rep <- NULL

  ## starting values
  theta <- matrix(0, k, 4*num_years)
  sigma2 <- apply(y, 1, var)
  xi <- replicate(k, 1/rgamma(1, alpha, beta))

  for (iter in 1:num_iter) {
    ## step 1 sample the latent process
    ## https://www.jarad.me/courses/stat615/slides/DLMs/DLMs.pdf
    for (j in 1:k) {
      dlm_mod <- dlmModPoly(1, dV = sigma2[j], dW = xi[j] * sigma2[j],
        m0 = rnorm(1, a, R))
      dlm_filter <- dlmFilter(y[j, ], dlm_mod)
      theta[j, ] <- as.vector(dlmBSample(dlm_filter))[-1]
    }

    ## step 2 sample the hyperparameters
    for (j in 1:k) {
      alpha_star_j <- alpha[j] + (T - 1)/2
      beta_star_j <-
        zapsmall(beta[j] +
                   sum(0.5 * (theta[j, -1] - theta[j, -T])^2 / sigma2[j]))
      xi[j] <- zapsmall(1/rgamma(1, alpha_star_j, beta_star_j))
    }
    for (j in 1:k) {
      tau_star_j <- tau[j] + (2*T - 1)/2
      kappa_star_j <-
        zapsmall(kappa[j] + sum((y[j, ] - theta[j, ])^2) +
                   sum(0.5 * (theta[j, -1] - theta[j, -T])^2 / xi[j]))
      sigma2[j] <- zapsmall(1/rgamma(1, tau_star_j, kappa_star_j))
    }

    ## step 3 sample z
    y_trans <- matrix(0, num_years, 4*k)
    theta_trans <- matrix(0, num_years, 4*k)
    miss_trans <- matrix(FALSE, num_years, 4*k)
    for (t_prime in 1:num_years) {
      y_trans[t_prime, ] <- as.vector(
        t(y[, (4*(t_prime-1)+1):(4*t_prime)]))
      theta_trans[t_prime, ] <- as.vector(
        t(theta[, (4*(t_prime-1)+1):(4*t_prime)]))
      miss_trans[t_prime, ] <- as.vector(
        t(miss[, (4*(t_prime-1)+1):(4*t_prime)]))
    }
    z <- cbind(y_trans, y_agg)
    miss_z <- cbind(miss_trans, miss_agg)

    H <- rbind(
      diag(k * 4),
      do.call(cbind, replicate(k, diag(4), simplify=FALSE)),
      kronecker(diag(k), t(rep(1, 4))))
    mu <- theta_trans %*% t(H)
    V <- diag(unlist(lapply(sigma2, rep, times = 4)))
    SIGMA <- H %*% V %*% t(H)

    for (t_prime in 1:num_years) {
      if (any(miss_z[t_prime, ])) {
        mu_tm <- mu[t_prime, ][miss_z[t_prime, ]]
        SIGMA_mm <- SIGMA[miss_z[t_prime, ], miss_z[t_prime, ]]
        SIGMA_mo <- SIGMA[miss_z[t_prime, ], !miss_z[t_prime, ]]
        SIGMA_om <- SIGMA[!miss_z[t_prime, ], miss_z[t_prime, ]]
        SIGMA_oo <- SIGMA[!miss_z[t_prime, ], !miss_z[t_prime, ]]
        SIGMA_oo_decomp <- eigen(SIGMA_oo)
        P_SIGMA_oo <- zapsmall(SIGMA_oo_decomp$vectors)
        D_SIGMA_oo <- zapsmall(SIGMA_oo_decomp$values)
        D_SIGMA_oo_star_idx <- D_SIGMA_oo > 0
        D_SIGMA_oo_star_inv <- matrix(0, sum(D_SIGMA_oo_star_idx), sum(D_SIGMA_oo_star_idx))
        diag(D_SIGMA_oo_star_inv) <- zapsmall(1 / D_SIGMA_oo[D_SIGMA_oo_star_idx])
        P_SIGMA_oo_star <- P_SIGMA_oo[, D_SIGMA_oo_star_idx]
        SIGMA_oo_psedoinv <- zapsmall(
          P_SIGMA_oo_star %*% D_SIGMA_oo_star_inv %*% t(P_SIGMA_oo_star))
        gamma_tm <- zapsmall(
          mu_tm + SIGMA_mo %*% SIGMA_oo_psedoinv %*%
            (z[t_prime, ][!miss_z[t_prime, ]] - mu[t_prime, ][!miss_z[t_prime, ]]))
        Omega <- zapsmall(
          SIGMA_mm - SIGMA_mo %*% SIGMA_oo_psedoinv %*% SIGMA_om)
        Omega <- (Omega + t(Omega)) / 2 # solve numerical issue
        Omega_decomp <- eigen(Omega)
        P_omega <- zapsmall(Omega_decomp$vectors)
        D_omega <- zapsmall(Omega_decomp$values)
        D_omega_star_idx <- D_omega > 0
        D_omega_star_rank <- sum(D_omega_star_idx)
        if (D_omega_star_rank == 0) next
        D_omega_star_sqrt <- matrix(0, D_omega_star_rank, D_omega_star_rank)
        diag(D_omega_star_sqrt) <- zapsmall(sqrt(D_omega[D_omega_star_idx]))
        P_omega_star <- P_omega[, D_omega_star_idx]
        u <- mvrnorm(1, rep(0, D_omega_star_rank), diag(D_omega_star_rank))
        z[t_prime, ][miss_z[t_prime, ]] <- zapsmall(
          gamma_tm + P_omega_star %*% D_omega_star_sqrt %*% u)
      }
    }

    ## step 4 impute missing values
    y_impu <- matrix(0, k, 4*num_years)
    for (j in 1:k) {
      y_impu[j, ] <- as.vector(t(z[, (4*(j-1)+1):(4*j)]))
    }
    y_agg_impu <- z[, (4*k+1):dim(z)[2]]
    y[miss] <- y_impu[miss]
    y_agg[miss_agg] <- y_agg_impu[miss_agg]

    ## y_rep <- rbind(y_rep, y)
    ## y_agg_rep <- rbind(y_agg_rep, y_agg)
    cat(iter, "\r")
  }
  return(list(y = y, y_agg = y_agg, theta = theta))
}

get_qcew_data <- function(series) {
  num_years <- (dim(series)[2] - 2) / (4 + 1)
  k <- dim(series)[1] - 1
  T <- num_years * 4
  data <- series[, -1:-2]
  y <- matrix(0, k, 4*num_years)
  y_agg <- matrix(0, num_years, 4 + k)
  for (t_prime in 1:num_years) {
    ts_idx <- (4*(t_prime-1)+1):(4*t_prime)
    y_agg[t_prime, 1:4] <- as.matrix(data[1, ts_idx + (t_prime - 1)*1])
    y[, ts_idx] <- as.matrix(data[-1, ts_idx + (t_prime - 1)*1])
    for (j in 1:k) {
      y_agg[t_prime, (4+1):(4+k)] <- data[-1, 4*t_prime + t_prime*1]
    }
  }
  miss <- y == 0
  miss_agg <- y_agg == 0
  y[miss] <- NA
  y_agg[miss_agg] <- NA
  return(list(y = y, y_agg = y_agg,
    miss = miss, miss_agg = miss_agg, k = k, T = T))
}
