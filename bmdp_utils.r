require(Matrix)
require(MASS)
require(dlm)
source('./bmmi.r')

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
