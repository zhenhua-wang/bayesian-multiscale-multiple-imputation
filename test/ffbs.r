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
y <- y[1, ]

drawIGpost <- function(x, a=0, b=0) {
  return(1/rgamma(1, a+length(x)/2, b+sum(x^2)/2))
}

n.reps <- 10000
V.reps <- c()
W.reps <- c()
theta.reps <- matrix(0, n.reps, length(y))
V <- 1/rgamma(1, 0.01, 0.01)
W <- 1/rgamma(1, 0.01, 0.01)

for (i in 1:n.reps) {
  cat(i,"\r")
  # Sample states
  mod <- dlmModPoly(1, dV = V, dW = W)
  filt <- dlmFilter(y, mod)
  theta <- dlmBSample(filt)
  # Sample V and W
  V <- drawIGpost(y-theta[-1])
  W <- drawIGpost(theta[-1]-theta[-length(y)])
  # Save iterations
  V.reps[i] <- V
  W.reps[i] <- W
  theta.reps[i,] <- theta[-1]
}

plot(1:n.reps, V.reps)
plot(1:n.reps, W.reps)
