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
