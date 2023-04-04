library(Matrix)
library(MASS)
library(dlm)

pinv <- function(mat) {
  ## assuming mat is symmetric
  mat_decomp <- eigen(mat)
  P_mat <- zapsmall(mat_decomp$vectors)
  D_mat <- zapsmall(mat_decomp$values)
  D_mat_star_idx <- D_mat > 0
  D_mat_star_inv <- matrix(0, sum(D_mat_star_idx), sum(D_mat_star_idx))
  diag(D_mat_star_inv) <- 1 / D_mat[D_mat_star_idx]
  P_mat_star <- P_mat[, D_mat_star_idx]
  mat_pinv <- P_mat_star %*% D_mat_star_inv %*% t(P_mat_star)
  return(mat_pinv)
}

rnorm_singular <- function(Mu, Var) {
  ## assuming omega is symmetric
  Var_decomp <- eigen(Var)
  P_Var <- zapsmall(Var_decomp$vectors)
  D_Var <- zapsmall(Var_decomp$values)
  D_Var_star_idx <- D_Var > 0
  D_Var_star_rank <- sum(D_Var_star_idx)
  D_Var_star_sqrt <- matrix(0, D_Var_star_rank, D_Var_star_rank)
  diag(D_Var_star_sqrt) <- sqrt(D_Var[D_Var_star_idx])
  P_Var_star <- P_Var[, D_Var_star_idx]
  u <- mvrnorm(1, rep(0, D_Var_star_rank), diag(D_Var_star_rank))
  sample <- Mu + P_Var_star %*% D_Var_star_sqrt %*% u
  return(sample)
}

y_to_ytrans <- function(y, num_years, k) {
  y_trans <- NULL
  for (t_prime in 1:num_years) {
    y_trans <- rbind(y_trans, as.vector(
      t(y[, (4*(t_prime-1)+1):(4*t_prime)])))
  }
  return(y_trans)
}

ytrans_to_y <- function(y_trans, num_years, k) {
  y <- NULL
  for (j in 1:k) {
    y <- rbind(y, as.vector(t(y_trans[, (4*(j-1)+1):(4*j)])))
  }
  return(y)
}

bmmi <- function(num_iter, num_burning, y, y_agg, miss, miss_agg,
                 a, R, tau, kappa, alpha, beta) {
  k <- dim(y)[1]
  N <- dim(y)[2]
  num_years <- N / 4
  y_rep <- NULL
  y_agg_rep <- NULL
  sigma2_rep <- NULL
  xi_rep <- NULL
  theta_rep <- NULL

  ## starting values
  theta <- matrix(0, k, 4*num_years)
  sigma2 <- apply(y, 1, var)
  xi <- replicate(k, 1/rgamma(1, alpha, beta))

  for (iter in 1:num_iter) {
    ## step 1 sample the latent process
    ## https://www.jarad.me/courses/stat615/slides/DLMs/DLMs.pdf
    for (j in 1:k) {
      dlm_mod <- dlmModPoly(1, dV = sigma2[j], dW = xi[j] * sigma2[j],
        m0 = a, C0 = R)
      dlm_filter <- dlmFilter(y[j, ], dlm_mod)
      theta[j, ] <- as.vector(dlmBSample(dlm_filter))[-1]
    }

    ## step 2 sample the hyperparameters
    for (j in 1:k) {
      ## sample xi
      alpha_star_j <- alpha[j] + (N - 1)/2
      beta_star_j <- beta[j] +
        sum(0.5 * (theta[j, -1] - theta[j, -N])^2 / sigma2[j])
      xi[j] <- 1/rgamma(1, alpha_star_j, beta_star_j)
      ## sample sigma2
      tau_star_j <- tau[j] + (2*N - 1)/2
      kappa_star_j <- kappa[j] + sum((y[j, ] - theta[j, ])^2) +
        sum(0.5 * (theta[j, -1] - theta[j, -N])^2 / xi[j])
      sigma2[j] <- 1/rgamma(1, tau_star_j, kappa_star_j)
    }

    ## step 3 sample z
    y_trans <- y_to_ytrans(y, num_years, k)
    theta_trans <- y_to_ytrans(theta, num_years, k)
    miss_trans <- y_to_ytrans(miss, num_years, k)
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
        mis_idx <- miss_z[t_prime, ]
        obs_idx <- !miss_z[t_prime, ]
        mu_tm <- mu[t_prime, ][mis_idx]
        ## sample from posterior
        SIGMA_oo_pinv <- pinv(SIGMA[obs_idx, obs_idx])
        gamma_tm <-
          mu_tm + SIGMA[mis_idx, obs_idx] %*% SIGMA_oo_pinv %*%
          (z[t_prime, ][obs_idx] - mu[t_prime, ][obs_idx])
        Omega <-
          SIGMA[mis_idx, mis_idx] -
          SIGMA[mis_idx, obs_idx] %*% SIGMA_oo_pinv %*%
          SIGMA[obs_idx, mis_idx]
        Omega <- (Omega + t(Omega)) / 2 # solve numerical issue
        z[t_prime, ][mis_idx] <- rnorm_singular(gamma_tm, Omega)
      }
    }

    ## step 4 impute missing values
    y_impu <- ytrans_to_y(z[, 1:(4*k)], num_years, k)
    y_agg_impu <- z[, (4*k+1):dim(z)[2]]
    y[miss] <- y_impu[miss]
    y_agg[miss_agg] <- y_agg_impu[miss_agg]

    if (iter > num_burning) {
      y_rep <- rbind(y_rep, y)
      y_agg_rep <- rbind(y_agg_rep, y_agg)
      sigma2_rep <- rbind(sigma2_rep, sigma2)
      xi_rep <- rbind(xi_rep, xi)
      theta_rep <- rbind(theta_rep, theta)
    }
    cat(iter, "\r")
  }
  return(list(y = y, y_agg = y_agg, theta = theta,
    y_rep = y_rep, y_agg_rep = y_agg_rep,
    sigma2_rep = sigma2_rep, xi_rep = xi_rep,
    theta_rep = theta_rep))
}

get_qcew_data <- function(series) {
  num_years <- (dim(series)[2] - 2) / (4 + 1)
  k <- dim(series)[1] - 1
  N <- num_years * 4
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
    miss = miss, miss_agg = miss_agg, k = k, N = N))
}

plot_series <- function(y, miss) {
  k <- dim(y)[1]
  N <- dim(y)[2]
  x <- 1:N
  par(mfrow = c(1, 3))
  for (j in 1:k) {
    mis_j <- miss[j, ]
    obs_j <- !miss[j, ]
    plot(x, y[j, ], type = "l")
    points(x[obs_j], y[j, obs_j], col = "black")
    points(x[mis_j], y[j, mis_j], col = "blue", pch = 17, cex = 2)
  }
}


plot_mcmc_series <- function(y_mcmc, miss, theta_mcmc, return_table = FALSE) {
  k <- dim(y_mcmc)[1]
  N <- dim(y_mcmc)[3]
  ## get mean
  y_mean <- apply(y_mcmc, c(1, 3), mean)
  theta_mean <- apply(theta_mcmc, c(1, 3), mean)
  ## get 95 CI
  y_lower <- apply(y_mcmc, c(1, 3), quantile, probs = 0.025)
  y_upper <- apply(y_mcmc, c(1, 3), quantile, probs = 0.975)
  theta_lower <- apply(theta_mcmc, c(1, 3), quantile, probs = 0.025)
  theta_upper <- apply(theta_mcmc, c(1, 3), quantile, probs = 0.975)
  x <- 1:N
  ## par(mfrow = c(3, 1))
  for (j in 1:k) {
    mis_j <- miss[j, ]
    obs_j <- !miss[j, ]
    ylim_upper <- max(y_mean[j, ], y_upper[j, mis_j], theta_upper[j, ])
    ylim_lower <- min(y_mean[j, ], y_lower[j, mis_j], theta_lower[j, ])
    plot(x, theta_mean[j, ], type = "l", ylim = c(ylim_lower, ylim_upper))
    abline(h=0)
    lines(x, theta_lower[j, ], type = "l", col = "red", lty=2)
    lines(x, theta_upper[j, ], type = "l", col = "red", lty=2)
    points(x[obs_j], y_mean[j, obs_j], col = "black")
    points(x[mis_j], y_mean[j, mis_j], col = "blue", pch = 17, cex = 1)
    arrows(x[mis_j], y_lower[j, mis_j], x[mis_j], y_upper[j, mis_j],
      length = 0.05, angle = 90, code = 3, lty = 2)
  }
  ## return CI
  if (return_table) {
    CI_table <- data.frame(
      y_mean = as.vector(t(y_mean)),
      y_lower = as.vector(t(y_lower)),
      y_upper = as.vector(t(y_upper)),
      serie_id = rep(seq(1, k), each = N))
    return(CI_table)
  }
}

impute_series <- function(data_list) {
  k <- data_list$k
  N <- data_list$N
  num_years <- N / 4
  y <- data_list$y / 1000
  y_agg <- data_list$y_agg / 1000
  miss <- data_list$miss
  miss_agg <- data_list$miss_agg
  num_iter <- 10000
  num_burning <- 5000
  a <- mean(y, na.rm = TRUE)
  R <- 10^10
  tau <- rep(0.01, k)
  kappa <- rep(0.01, k)
  alpha <- rep(3, k)
  beta <- rep(0.1, k)
  ## initial impute
  ## time series
  for (j in 1:k) {
    y[j, ][miss[j, ]] <- mean(y[j, ], na.rm=TRUE)
  }
  ## quarterly total
  for (t in 1:num_years) {
    for (q in 1:4) {
      miss_curr <- miss_agg[t, q]
      if (miss_curr) {
        y_agg[t, q] <- sum(y[, (t-1)*4 + q])
      }
    }
  }
  ## annual total
  for (t in 1:num_years) {
    for (j in 1:k) {
      miss_curr <- miss_agg[t, 4+j]
      if (miss_curr) {
        y_agg[t, 4+j] <- sum(y[j, ((t-1)*4+1):(t*4)])
      }
    }
  }
  ## impute
  imputed <- bmmi(num_iter, num_burning, y, y_agg, miss, miss_agg, a, R, tau, kappa, alpha, beta)
  ## plot
  y_imputed_mcmc <- array(0, dim = c(k, (num_iter-num_burning), N))
  theta_imputed_mcmc <- array(0, dim = c(k, (num_iter-num_burning), N))
  for (j in 1:k) {
    for (i in 1:(num_iter-num_burning)) {
      y_imputed_mcmc[j, i, ] <- imputed$y_rep[((i-1)*k + 1):(i*k), ][j, ]
      theta_imputed_mcmc[j, i, ] <- imputed$theta_rep[((i-1)*k + 1):(i*k), ][j, ]
    }
  }
  return(plot_mcmc_series(y_imputed_mcmc, miss, theta_imputed_mcmc, return_table = TRUE))
}
