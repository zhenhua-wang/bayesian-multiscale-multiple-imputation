library(Matrix)
library(MASS)
library(dlm)
library(VGAM)
source("./bmmi.r")

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

bmdp <- function(num_iter, num_burning, y,
                 a, R, tau, kappa, alpha, beta,
                 epsilon, DELTA, hidden = FALSE) {
  k <- dim(y)[1]
  N <- dim(y)[2]
  num_years <- N / 4
  b_dp <- DELTA / epsilon
  y_true_rep <- NULL
  y_true_agg_rep <- NULL
  sigma2_rep <- NULL
  xi_rep <- NULL
  theta_rep <- NULL

  ## intermediate values
  H <- rbind(
    diag(k * 4),
    do.call(cbind, replicate(k, diag(4), simplify=FALSE)),
    kronecker(diag(k), t(rep(1, 4))))
  y_trans <- y_to_ytrans(y, num_years, k)
  y_agg <- matrix(NA, num_years, 4 + k)
  y_agg <- compute_y_agg_from_y(y, y_agg)
  z <- cbind(y_trans, y_agg)

  ## starting values
  theta <- matrix(0, k, 4*num_years)
  sigma2 <- apply(y, 1, var)
  xi <- replicate(k, 1/rgamma(1, alpha, beta))
  y_true <- y
  if (hidden) {
    y_true[1, ] <- laplace_mechanism(y_true[1, ], epsilon, DELTA)
    y_true[2, ] <- laplace_mechanism(y_true[2, ], epsilon, DELTA)
    y_true[3, ] <- laplace_mechanism(y_true[3, ], epsilon, DELTA)
    y_true <- y_true
  }
  y_agg_true <- matrix(NA, num_years, 4 + k)
  y_agg_true <- compute_y_agg_from_y(y_true, y_agg_true)
  y_true_trans <- y_to_ytrans(y_true, num_years, k)
  z_true <- cbind(y_true_trans, y_agg_true)

  for (iter in 1:num_iter) {
    ## step 1 sample true total value
    theta_trans <- y_to_ytrans(theta, num_years, k)
    V <- diag(unlist(lapply(sigma2, rep, times = 4)))
    if (hidden) {
      for(year in 1:num_years) {
        ## sample latent variable x_t'_i
        x_year <- c()
        for (i in 1:dim(y_true_trans)[2]) {
          temp <- abs(2 * (y_trans[year, i] - y_true_trans[year, i]))
          if (temp == 0) {
            x_year <- c(x_year, 1 / rgamma(1, 1 / 2, 1 / 8))
          } else {
            mu_x <- (b_dp / 2) / temp
            lambda_x <- 1 / 4
            x_year <- c(x_year, rinv.gaussian(1, mu_x, lambda_x))
          }
        }
        ## sample z_t'
        Q <- pinv(diag(x_year))
        ## mu <- H %*% theta_trans[year, ]
        ## SIGMA <- H %*% V %*% t(H)
        ## Omega_z <- pinv(pinv(Q) + pinv(SIGMA))
        Omega_z <- pinv(pinv(Q) + pinv(V))
        Omega_z <- (Omega_z + t(Omega_z)) / 2
        ## Theta_z <- Omega_z %*% (pinv(Q) %*% z[year, ] + pinv(SIGMA) %*% mu)
        Theta_z <- Omega_z %*% (pinv(Q) %*% y_trans[year, ] + pinv(V) %*% theta_trans[year, ])
        y_true_trans[year, ] <- rnorm_singular(Theta_z, Omega_z)
      }
      ## update y_true and y_agg_true
      ## y_true <- ytrans_to_y(z_true[, 1:(4*k)], num_years, k)
      ## y_agg_true <- z_true[, (4*k+1):dim(z_true)[2]]
      y_true <- ytrans_to_y(y_true_trans, num_years, k)
      y_agg_true <- compute_y_agg_from_y(y_true, y_agg_true)
    }

    ## step 2 sample the latent process
    for (j in 1:k) {
      dlm_mod <- dlmModPoly(1, dV = sigma2[j], dW = xi[j] * sigma2[j],
        m0 = a, C0 = R)
      dlm_filter <- dlmFilter(y_true[j, ], dlm_mod)
      theta[j, ] <- as.vector(dlmBSample(dlm_filter))[-1]
    }

    ## step 3 sample the hyperparameters
    for (j in 1:k) {
      ## sample xi
      alpha_star_j <- alpha[j] + (N - 1)/2
      beta_star_j <- beta[j] +
        sum(0.5 * (theta[j, -1] - theta[j, -N])^2 / sigma2[j])
      xi[j] <- 1/rgamma(1, alpha_star_j, beta_star_j)
      ## sample sigma2
      tau_star_j <- tau[j] + (2*N - 1)/2
      kappa_star_j <- kappa[j] + sum(0.5 * (y_true[j, ] - theta[j, ])^2) +
        sum(0.5 * (theta[j, -1] - theta[j, -N])^2 / xi[j])
      sigma2[j] <- 1/rgamma(1, tau_star_j, kappa_star_j)
    }

    if (iter > num_burning) {
      y_true_rep <- rbind(y_true_rep, y_true)
      y_true_agg_rep <- rbind(y_true_agg_rep, y_agg_true)
      sigma2_rep <- rbind(sigma2_rep, sigma2)
      xi_rep <- rbind(xi_rep, xi)
      theta_rep <- rbind(theta_rep, theta)
    }
    cat(iter, "\r")
  }
  return(list(y = y, y_agg = y_agg, theta = theta,
    y_true_rep = y_true_rep, y_true_agg_rep = y_true_agg_rep,
    sigma2_rep = sigma2_rep, xi_rep = xi_rep,
    theta_rep = theta_rep))
}

gaussian_mechanism <- function(z, epsilon, delta, DELTA) {
  sigma_N <- DELTA * (c + sqrt(c^2 + epsilon)) / (epsilon * sqrt(2))
  N <- rnorm(length(z), 0, sigma_N)
  z + N
}

laplace_mechanism <- function(z, epsilon, DELTA) {
  b_dp <- DELTA / epsilon
  rlaplace(length(z), z, b_dp)
}

compute_y_agg_from_y <- function(y, y_agg,
                                 miss_agg = matrix(1,
                                   nrow(y_agg),
                                   ncol(y_agg))) {
  num_years <- dim(y_agg)[1]
  k <- dim(y)[1]
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
  return(y_agg)
}

release_series <- function(data_list, epsilon, delta, DELTA, num_iter = 10000, num_burning = 5000, DP = FALSE, hidden = FALSE) {
  k <- data_list$k
  N <- data_list$N
  num_years <- N / 4
  if (DP) {
    y <- data_list$y_dp / 1000
  } else {
    y <- data_list$y / 1000}
  y_agg <- data_list$y_agg / 1000
  miss <- data_list$miss
  miss_agg <- data_list$miss_agg
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
  imputed <- bmdp(num_iter, num_burning, y, a, R, tau, kappa, alpha, beta, epsilon, DELTA, hidden)
  return(imputed)
}
