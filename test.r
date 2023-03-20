library(dlm)
load("./data/series1")
source("./bmmi.r")
data_list <- get_qcew_data(series1)
k <- data_list$k
T <- data_list$T
y <- data_list$y
for (j in 1:k) {
  y[j, ][is.na(y[j, ])] <- mean(y[j, ], na.rm = TRUE)
}
## y <- Nile

drawIGpost <- function(x, a=0, b=0) {
  return(1/rgamma(1, a+length(x)/2, b+sum(x^2)/2))
}

num_iters <- 3000
tau <- 0.01
kappa <- 0.01
alpha <- 3
beta <- 0.1
theta <- array(0, dim = c(num_iters, dim(y)[1], dim(y)[2]))
sigma2 <- replicate(3, 1/rgamma(1, alpha, beta))
xi <- replicate(3, 1/rgamma(1, tau, kappa))
for (i in 1:num_iters) {
  cat(i,"\r")
  for (j in 1:k) {
    # Sample states
    mod <- dlmModPoly(1, dV = sigma2[j], dW = xi[j], rnorm(1, 0, 10^10))
    filt <- dlmFilter(y[j, ], mod)
    theta_curr <- dlmBSample(filt)
    theta[i, j, ] <- theta_curr[-1]
    # Sample V and W
    alpha_star <- alpha + (T - 1)/2
    beta_star <- beta +
      sum(0.5 * (theta_curr[-1] - theta_curr[-T])^2 / sigma2[j])
    xi[j] <- 1/rgamma(1, alpha_star, beta_star)
    tau_star <- tau + (2*T - 1)/2
    kappa_star <- kappa + sum((y[j, ] - theta_curr[-1])^2) +
      sum(0.5 * (theta_curr[-1] - theta_curr[-T])^2 / xi[j])
    sigma2[j] <- 1/rgamma(1, tau_star, kappa_star)
  }
}
