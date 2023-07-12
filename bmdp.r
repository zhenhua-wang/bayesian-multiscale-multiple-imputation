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

bmdp <- function(num_iter, num_burning, y, y_agg, miss, miss_agg,
                 a, R, tau, kappa, alpha, beta,
                 epsilon, DELTA, hidden = FALSE) {
  k <- dim(y)[1]
  N <- dim(y)[2]
  num_years <- N / 4
  y_true_rep <- NULL
  y_agg_rep <- NULL
  sigma2_rep <- NULL
  xi_rep <- NULL
  theta_rep <- NULL

  ## starting values
  theta <- matrix(0, k, 4*num_years)
  sigma2 <- apply(y, 1, var)
  xi <- replicate(k, 1/rgamma(1, alpha, beta))
  b_dp <- DELTA / epsilon
  y_true <- y

  for (iter in 1:num_iter) {
    ## step 1 sample true total value
    if (hidden) {
      for(j in 1:k) {
        for (t in 1:N) {
          current_x <- y_true[j, k]
          proposed_x <- rlaplace(1, current_x, b_dp)
          log.r <-
            dnorm(y[j, k], proposed_x, b_dp, log = TRUE) -
            dnorm(y[j, k], current_x, b_dp, log = TRUE) +
            dnorm(proposed_x, theta[j, k], sqrt(sigma2[j]), log = TRUE) -
            dnorm(current_x, theta[j, k], sqrt(sigma2[j]), log = TRUE)
          if (log(runif(1)) < log.r) {
            y_true[j, k] <- proposed_x
          } else {
            y_true[j, k] <- current_x
          }
        }
      }
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
      y_agg_rep <- rbind(y_agg_rep, y_agg)
      sigma2_rep <- rbind(sigma2_rep, sigma2)
      xi_rep <- rbind(xi_rep, xi)
      theta_rep <- rbind(theta_rep, theta)
    }
    cat(iter, "\r")
  }
  return(list(y = y, y_agg = y_agg, theta = theta,
    y_true_rep = y_true_rep, y_agg_rep = y_agg_rep,
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
  imputed <- bmdp(num_iter, num_burning, y, y_agg, miss, miss_agg, a, R, tau, kappa, alpha, beta, epsilon, DELTA, hidden)
  ## plot
  y_imputed_mcmc <- array(0, dim = c(k, (num_iter-num_burning), N))
  theta_imputed_mcmc <- array(0, dim = c(k, (num_iter-num_burning), N))
  for (j in 1:k) {
    for (i in 1:(num_iter-num_burning)) {
      y_imputed_mcmc[j, i, ] <- imputed$y_true_rep[((i-1)*k + 1):(i*k), ][j, ]
      theta_imputed_mcmc[j, i, ] <- imputed$theta_rep[((i-1)*k + 1):(i*k), ][j, ]
    }
  }
  return(plot_mcmc_series(y_imputed_mcmc, miss, theta_imputed_mcmc))
}
