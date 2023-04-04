library(Matrix)
library(MASS)
library(dlm)
source("./bmmi.r")

load("./data/series1")
data_list <- get_qcew_data(series1)
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
theta_imputed_mcmc <- array(0, dim = c(k, (num_iter-num_burning), N))
for (j in 1:k) {
  for (i in 1:(num_iter-num_burning)) {
    y_imputed_mcmc[j, i, ] <- imputed$y_rep[((i-1)*k + 1):(i*k), ][j, ]
    theta_imputed_mcmc[j, i, ] <- imputed$theta_rep[((i-1)*k + 1):(i*k), ][j, ]
  }
}

y_mcmc <- y_imputed_mcmc
k <- dim(y_mcmc)[1]
N <- dim(y_mcmc)[3]
## get mean
y_mean <- apply(y_mcmc, c(1, 3), mean)
## get 95 CI
y_lower <- apply(y_mcmc, c(1, 3), quantile, probs = 0.025)
y_upper <- apply(y_mcmc, c(1, 3), quantile, probs = 0.975)
x <- 1:N
par(mfrow = c(3, 1))
for (j in 1:k) {
  mis_j <- miss[j, ]
  obs_j <- !miss[j, ]
  ylim_upper <- max(y_mean[j, ], y_upper[j, mis_j])
  ylim_lower <- min(y_mean[j, ], y_lower[j, mis_j])
  plot(x, y_mean[j, ], type = "l", ylim = c(ylim_lower, ylim_upper))
  points(x[obs_j], y_mean[j, obs_j], col = "black")
  points(x[mis_j], y_mean[j, mis_j], col = "blue", pch = 17, cex = 1)
  arrows(x[mis_j], y_lower[j, mis_j], x[mis_j], y_upper[j, mis_j],
    length = 0.05, angle = 90, code = 3, lty = 2)
}

y_CI <- plot_mcmc_series(y_imputed_mcmc, miss, theta_imputed_mcmc, TRUE)
y_CI$y <- as.vector(t(y_mean))
y_CI[as.vector(t(!miss)), c("y_mean", "y_lower", "y_upper")] <- NA
y_CI[y_CI$serie_id == 1, ]
