ffbs <- function(y, F, G, V, W) {
  # Initialization
  n <- length(y)
  m <- rep(0, n)
  C <- rep(0, n)
  alpha <- rep(0, n)
  beta <- rep(0, n)
  C[n] <- 1 / (G + V)
  m[n] <- y[n] - F
  alpha[n] <- m[n] / (G + V)

  # Forward filtering step
  for (t in 2:n) {
    C[t-1] <- 1 / (1 / C[t] + 1 / W)
    m[t-1] <- (m[t] / C[t] + alpha[t-1]) * C[t-1]
    alpha[t-1] <- F * m[t-1]
    beta[t] <- (alpha[t-1] - m[t]) / C[t]
    C[t-1] <- C[t-1] * (1 - F^2)
  }

  # Backward smoothing step
  m_smooth <- m
  C_smooth <- C
  for (t in (n-1):1) {
    m_smooth[t] <- m[t] + C[t] * (F * m_smooth[t+1] - alpha[t])
    C_smooth[t] <- C[t] + C[t] * (F^2 * C_smooth[t+1] - C[t] - beta[t+1]) / C_smooth[t+1]
  }

  # Return results
  return(list(m=m_smooth, C=C_smooth))
}

# Generate some data
set.seed(123)
n <- 100
y <- rnorm(n, mean=0, sd=1)
F <- 0.9
G <- 1
V <- 1
W <- 0.1

# Run the FFBS algorithm
results <- ffbs(y, F, G, V, W)

# Plot the results
plot(results$m, type="l", main="Smoothed Means")
plot(results$C, type="l", main="Smoothed Variances")
