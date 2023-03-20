library(Matrix)
library(MASS)
library(dlm)
source("./bmmi.r")

num_iter <- 1
k <- 1
T <- 4
num_years <- T / 4
a <- 0
R <- 10^10
tau <- rep(0.01, k)
kappa <- rep(0.01, k)
alpha <- rep(3, k)
beta <- rep(0.1, k)

## data
y <- matrix(c(1.5, 2, 5, 1, 9), num_years, k * 4)
H <- rbind(
  diag(k * 4),
  do.call(cbind, replicate(k, diag(4), simplify=FALSE)),
  kronecker(diag(k), t(rep(1, 4))))
y_all <- y %*% t(H)
y <- matrix(0, k, 4*num_years)
for (j in 1:k) {
  y[j, ] <- as.vector(t(y_all[, (4*(j-1)+1):(4*j)]))
}
y_agg <- y_all[, (4*k+1):dim(y_all)[2]] # (k, num_years + 1)

miss <- matrix(as.logical(rbinom(y, 1, 0.3)), k, 4*num_years)
miss_agg <- matrix(as.logical(rbinom(y_agg, 1, 0.3)), num_years, num_years + 1)
y[miss] <- NA
y_agg[miss_agg] <- NA
for (j in 1:k) {
  y[j, ][miss[j, ]] <- mean(y[j, ], na.rm=TRUE)
}
for (t in 1:num_years) {
  y_agg[t, ][miss_agg[t, ]] <- mean(y_agg[t, ], na.rm=TRUE)
}

## run
bmmi(num_iter, y, y_agg, miss, miss_agg, a, R, tau, kappa, alpha, beta)


## test
library(monomvn)
y_trans <- matrix(0, num_years, 4*k)
theta_trans <- matrix(0, num_years, 4*k)
miss_trans <- matrix(FALSE, num_years, 4*k)
for (t_prime in 1:num_years) {
  y_trans[t_prime, ] <- as.vector(t(y[, (4*(t_prime-1)+1):(4*t_prime)]))
  theta_trans[t_prime, ] <- as.vector(t(theta[, (4*(t_prime-1)+1):(4*t_prime)]))
  miss_trans[t_prime, ] <- as.vector(t(miss[, (4*(t_prime-1)+1):(4*t_prime)]))
}
z <- cbind(y_trans, y_agg)
miss_z <- cbind(miss_trans, miss_agg)
H <- rbind(
  diag(k * 4),
  do.call(cbind, replicate(k, diag(4), simplify=FALSE)),
  kronecker(diag(k), t(rep(1, 4))))
mu <- theta_trans %*% t(H)
V <- diag(unlist(lapply(sigma2, rep, times = 4)))
SIGMA <- solve(rwish(3, matrix(replicate(9*9, 1/rgamma(1, 1, 1)), 9, 9)))

t_prime = 1
mu_tm <- mu[t_prime, ][miss_z[t_prime, ]]
SIGMA_mm <- SIGMA[miss_z[t_prime, ], miss_z[t_prime, ]]
SIGMA_mo <- SIGMA[miss_z[t_prime, ], !miss_z[t_prime, ]]
SIGMA_om <- SIGMA[!miss_z[t_prime, ], miss_z[t_prime, ]]
SIGMA_oo <- SIGMA[!miss_z[t_prime, ], !miss_z[t_prime, ]]
SIGMA_oo_decomp <- eigen(SIGMA_oo)
P <- SIGMA_oo_decomp$vectors
D_star_idx <- SIGMA_oo_decomp$values > positive_threhold
D_star_inv <- matrix(0, sum(D_star_idx), sum(D_star_idx))
diag(D_star_inv) <- 1 / SIGMA_oo_decomp$values[D_star_idx]
P_star <- P[,  D_star_idx]
SIGMA_oo_psedoinv <- P_star %*% D_star_inv %*% t(P_star)
gamma_tm <- mu_tm - SIGMA_mo %*% SIGMA_oo_psedoinv %*%
  (z[t_prime, ][!miss_z[t_prime, ]] - mu[t_prime, ][!miss_z[t_prime, ]])
Omega <- SIGMA_mm - SIGMA_mo %*% SIGMA_oo_psedoinv %*% SIGMA_om
Omega_decomp <- eigen(Omega)
P_omega <- Omega_decomp$vectors
D_omega_star_idx <- Omega_decomp$values > positive_threhold
D_omega_star_sqrt <- matrix(0, sum(D_omega_star_idx), sum(D_omega_star_idx))
diag(D_omega_star_sqrt) <- sqrt(Omega_decomp$values[D_omega_star_idx])
P_omega_star <- P_omega[, D_omega_star_idx]
Omega_rank <- rankMatrix(Omega)[1]
u <- mvrnorm(1, rep(0, sum(D_omega_star_idx)), diag(sum(D_omega_star_idx)))
z[t_prime, ][miss_z[t_prime, ]] <- gamma_tm + P_omega_star %*% D_omega_star_sqrt %*% u
