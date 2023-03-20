library(dlm)
load("./data/series1")
source("./bmmi.r")
data_list <- get_qcew_data(series1)
k <- data_list$k
N <- data_list$T
y <- data_list$y
for (j in 1:k) {
  y[j, ][is.na(y[j, ])] <- mean(y[j, ], na.rm = TRUE)
}

drawIGpost <- function(x, a=0, b=0) {
  return(1/rgamma(1, a+length(x)/2, b+sum(x^2)/2))
}

alpha <- 3
beta <- 0.1
tau <- 0.01
kappa <- 0.01
n.reps <- 50000
sigma2.reps <- matrix(0, n.reps, k)
xi.reps <- matrix(0, n.reps, k)
theta.reps <- array(0, dim = c(n.reps, j, dim(y)[2]))
sigma2 <- replicate(k, 1/rgamma(1, tau, kappa))
xi <- replicate(k, 1/rgamma(1, alpha, beta))

for (i in 1:n.reps) {
  cat(i,"\r")
  for (j in 1:k) {
    # Sample states
    mod <- dlmModPoly(1, dV = sigma2[j], dW = xi[j] * sigma2[j])
    filt <- dlmFilter(y[j, ], mod)
    theta <- dlmBSample(filt)[-1]
    # Sample V and W
    alpha_star <- alpha + (N - 1)/2
    beta_star <- beta +
      sum(0.5 * (theta[-1] - theta[-N])^2 / sigma2[j])
    xi <- 1/rgamma(k, alpha_star, beta_star)
    tau_star <- tau + (2*N - 1)/2
    kappa_star <- kappa + sum((y[j, ] - theta)^2) +
      sum(0.5 * (theta[-1] - theta[-N])^2 / xi[j])
    sigma2 <- 1/rgamma(k, tau_star, kappa_star)
    # Save iterations
    sigma2.reps[i, ] <- sigma2
    xi.reps[i, ] <- xi
    theta.reps[i, j, ] <- theta
  }
}

plot(1:n.reps, sigma2.reps[, 1])
plot(1:n.reps, xi.reps[, 1])
