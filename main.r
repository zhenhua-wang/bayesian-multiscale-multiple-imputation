library(Matrix)
library(MASS)
library(dlm)
source("./bmmi.r")

load("./data/series1")
data_list <- get_qcew_data(series1)
k <- data_list$k
T <- data_list$T
num_years <- T / 4
y <- data_list$y
y_agg <- data_list$y_agg
miss <- data_list$miss
miss_agg <- data_list$miss_agg

num_iter <- 1
a <- 0
R <- 10^10
tau <- rep(0.01, k)
kappa <- rep(0.01, k)
alpha <- rep(3, k)
beta <- rep(0.1, k)

## initial impute
for (j in 1:k) {
  y[j, ][miss[j, ]] <- mean(y[j, ], na.rm=TRUE)
}
for (t in 1:num_years) {
  y_agg[t, ][miss_agg[t, ]] <- mean(y_agg[t, ], na.rm=TRUE)
}

## run
bmmi(num_iter, y, y_agg, miss, miss_agg, a, R, tau, kappa, alpha, beta)
