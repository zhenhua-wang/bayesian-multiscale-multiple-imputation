library(dlm)
load("./data/series1")
source("./bmmi.r")
data_list <- get_qcew_data(series1)
k <- data_list$k
y <- data_list$y
for (j in 1:k) {
  y[j, ][is.na(y[j, ])] <- mean(y[j, ], na.rm = TRUE)
}
## y <- Nile

drawIGpost <- function(x, a=0, b=0) {
  return(1/rgamma(1, a+length(x)/2, b+sum(x^2)/2))
}

num_iters <- 5000
theta <- array(0, dim = c(num_iters, dim(y)[1], dim(y)[2]))
V <- replicate(3, 1/rgamma(1, 3, 0.1))
W <- replicate(3, 1/rgamma(1, 0.01, 0.01))
for (i in 1:num_iters) {
  cat(i,"\r")
  for (j in 1:k) {
    # Sample states
    mod <- dlmModPoly(1, dV = V[j], dW = W[j], rnorm(1, 0, 10^10))
    filt <- dlmFilter(y[j, ], mod)
    theta_curr <- dlmBSample(filt)
    # Sample V and W
    V[j] <- drawIGpost(y[j, ]-theta_curr[-1])
    W[j] <- drawIGpost(theta_curr[-1] - theta_curr[-length(theta_curr)])
    theta[i, j, ] <- theta_curr[-1]
  }
}
