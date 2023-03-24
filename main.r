library(Matrix)
library(MASS)
library(dlm)
source("./bmmi.r")

load("./data/series3")
data_list <- get_qcew_data(series3)
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

plot_series(y, miss)

## run
imputed <- bmmi(num_iter, num_burning, y, y_agg, miss, miss_agg, a, R, tau, kappa, alpha, beta)
y_imputed <- imputed$y
y_agg_imputed <- imputed$y_agg
theta_imputed <- imputed$theta
sigma2_imputed <- imputed$sigma2_rep
xi_imputed <- imputed$xi_rep

png(filename="./mcmc.png")
par(mfrow = c(2, 3))
plot(1:(num_iter-num_burning), imputed$sigma2_rep[, 1], xlab = "", ylab = "sigma2[1]")
plot(1:(num_iter-num_burning), imputed$sigma2_rep[, 2], xlab = "", ylab = "sigma2[2]")
plot(1:(num_iter-num_burning), imputed$sigma2_rep[, 3], xlab = "", ylab = "sigma2[3]")
plot(1:(num_iter-num_burning), imputed$xi_rep[, 1], xlab = "", ylab = "x1[1]")
plot(1:(num_iter-num_burning), imputed$xi_rep[, 2], xlab = "", ylab = "x1[2]")
plot(1:(num_iter-num_burning), imputed$xi_rep[, 3], xlab = "", ylab = "x1[3]")
dev.off()

y_imputed_mcmc <- array(0, dim = c(k, (num_iter-num_burning), N))
for (j in 1:k) {
  for (i in 1:(num_iter-num_burning)) {
    y_imputed_mcmc[j, i, ] <- imputed$y_rep[((i-1)*k + 1):(i*k), ][j, ]
  }
}

plot_mcmc_series(y_imputed_mcmc, miss)
