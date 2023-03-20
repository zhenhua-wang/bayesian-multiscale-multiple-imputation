library(Matrix)
library(MASS)
library(dlm)

bmmi <- function(num_iter, y, y_agg, miss, miss_agg,
                 a, R, tau, kappa, alpha, beta,
                 positive_threhold = 1e-10) {
  k <- dim(y)[1]
  T <- dim(y)[2]
  num_years <- T / 4

  ## starting values
  theta0 <- rnorm(1, a, R)
  theta <- matrix(0, k, 4*num_years)
  sigma2 <- apply(y, 1, var)
  xi <- c(1, 1, 1)

  for (iter in 1:num_iter) {
    ## step 1 sample the latent process
    ## https://www.jarad.me/courses/stat615/slides/DLMs/DLMs.pdf
    for (j in 1:k) {
      dlm_mod <- dlmModPoly(order = 1, dV = sigma2[j], dW = xi[j] * sigma2[j], m0 = theta0)
      dlm_filter <- dlmFilter(y[j, ], dlm_mod)
      theta[j, ] <- as.vector(dlmBSample(dlm_filter))[-1]
    }

    ## step 2 sample the hyperparameters
    alpha_star <- alpha + (T - 1)/2
    beta_star <- c()
    for (j in 1:k) {
      beta_star[j] <- beta[j] +
        sum(0.5 * (theta[j, -1] - theta[j, -T])^2 / sigma2[j])
    }
    xi <- 1/rgamma(k, alpha_star, beta_star)
    tau_star <- tau + (2*T - 1)/2
    kappa_star <- c()
    for (j in 1:k) {
      kappa_star[j] <- kappa[j] + sum((y[j, ] - theta[j, ])^2) +
        sum(0.5 * (theta[j, -1] - theta[j, -T])^2 / xi[j])
    }
    sigma2 <- 1/rgamma(k, tau_star, kappa_star)

    ## step 3 sample z
    y_trans <- matrix(0, num_years, 4*k)
    theta_trans <- matrix(0, num_years, 4*k)
    miss_trans <- matrix(FALSE, num_years, 4*k)
    for (t_prime in 1:num_years) {
      y_trans[t_prime, ] <- as.vector(t(y[, (4*(t_prime-1)+1):(4*t_prime)]))
      theta_trans[t_prime, ] <- as.vector(t(theta[, (4*(t_prime-1)+1):(4*t_prime)]))
      miss_trans[t_prime, ] <- as.vector(t(miss[, (4*(t_prime-1)+1):(4*t_prime)]))
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
      mu_tm <- mu[t_prime, ][miss_z[t_prime, ]]
      SIGMA_mm <- SIGMA[miss_z[t_prime, ], miss_z[t_prime, ]]
      SIGMA_mo <- SIGMA[miss_z[t_prime, ], !miss_z[t_prime, ]]
      SIGMA_om <- SIGMA[!miss_z[t_prime, ], miss_z[t_prime, ]]
      SIGMA_oo <- SIGMA[!miss_z[t_prime, ], !miss_z[t_prime, ]]
      SIGMA_oo_decomp <- eigen(SIGMA_oo)
      P <- SIGMA_oo_decomp$vectors
      D_star_idx <- SIGMA_oo_decomp$values > positive_threhold
      D_star_inv <- matrix(0, sum(D_star_idx), sum(D_star_idx))
      diag(D_star_inv) <- 1 / SIGMA_oo_decomp$values[D_star_idx]
      P_star <- P[,  D_star_idx]
      SIGMA_oo_psedoinv <- P_star %*% D_star_inv %*% t(P_star)
      gamma_tm <- mu_tm - SIGMA_mo %*% SIGMA_oo_psedoinv %*%
        (z[t_prime, ][!miss_z[t_prime, ]] - mu[t_prime, ][!miss_z[t_prime, ]])
      Omega <- SIGMA_mm - SIGMA_mo %*% SIGMA_oo_psedoinv %*% SIGMA_om
      Omega_decomp <- eigen(Omega)
      P_omega <- Omega_decomp$vectors
      D_omega_star_idx <- Omega_decomp$values > positive_threhold
      D_omega_star_sqrt <- matrix(0, sum(D_omega_star_idx), sum(D_omega_star_idx))
      diag(D_omega_star_sqrt) <- sqrt(Omega_decomp$values[D_omega_star_idx])
      P_omega_star <- P_omega[, D_omega_star_idx]
      Omega_rank <- rankMatrix(Omega)[1]
      u <- mvrnorm(1, rep(0, sum(D_omega_star_idx)), diag(sum(D_omega_star_idx))) #TODO
      z[t_prime, ][miss_z[t_prime, ]] <- gamma_tm + P_omega_star %*% D_omega_star_sqrt %*% u
    }

    ## step 4 impute missing values
    y_impu <- matrix(0, k, 4*num_years)
    for (j in 1:k) {
      y_impu[j, ] <- as.vector(t(z[, (4*(j-1)+1):(4*j)]))
    }
    y_agg_impu <- z[, (4*k+1):dim(z)[2]]
    y[miss] <- y_impu[miss]
    y_agg[miss_agg] <- y_agg_impu[miss_agg]
  }
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
  return(list(y = y, y_agg = y_agg, miss = miss, miss_agg = miss_agg, k = k, T = T))
}
