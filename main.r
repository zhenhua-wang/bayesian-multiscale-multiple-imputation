library(Matrix)
library(MASS)
library(dlm)
source("./bmmi.r")

load("./data/series4")
data_list <- get_qcew_data(series4)
k <- data_list$k
T <- data_list$T
num_years <- T / 4
y <- data_list$y
y_agg <- data_list$y_agg
miss <- data_list$miss
miss_agg <- data_list$miss_agg

num_iter <- 10000
num_burning <- 5000
a <- 0
R <- 1
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
  for (j in (4+1):(4+k)) {
    miss_curr <- miss_agg[t, j]
    if (miss_curr) {
      y_agg[t, j] <- sum(y[j-4, ((t-1)*4+1):(t*4)])
    }
  }
}

## run
imputed <- bmmi(num_iter, y, y_agg, miss, miss_agg, a, R, tau, kappa, alpha, beta)
y_imputed <- imputed$y
y_agg_imputed <- imputed$y_agg
theta_imputed <- imputed$theta
