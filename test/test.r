library(dlm)
load("./data/series3")
source("./bmmi.r")
data_list <- get_qcew_data(series3)
k <- data_list$k
T <- data_list$T
y <- data_list$y
for (j in 1:k) {
  y[j, ][is.na(y[j, ])] <- mean(y[j, ], na.rm = TRUE)
}

drawIGpost <- function(x, a=0, b=0) {
  return(1/rgamma(1, a+length(x)/2, b+sum(x^2)/2))
}

n.reps <- 30000
V.reps <- c()
W.reps <- c()
theta.reps <- array(0, dim = c(n.reps, j, dim(y)[2]))
V <- apply(y, 1, var)
W <- replicate(3, 1/rgamma(1, 0.01, 0.01))

for (i in 1:n.reps) {
  cat(i,"\r")
  for (j in 1:k) {
    # Sample states
    mod <- dlmModPoly(1, dV = V[j], dW = W[j])
    filt <- dlmFilter(y[j, ], mod)
    theta <- dlmBSample(filt)[-1]
    # Sample V and W
    V[j] <- drawIGpost(y[j, ]-theta)
    W[j] <- drawIGpost(theta[-1]-theta[-length(y[, ])])
    # Save iterations
    V.reps[i] <- V
    W.reps[i] <- W
    theta.reps[i, j, ] <- theta
  }
}

plot(1:n.reps, V.reps)
plot(1:n.reps, W.reps)
