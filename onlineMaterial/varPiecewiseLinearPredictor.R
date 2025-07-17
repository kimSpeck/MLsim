# Below is an R snippet demonstrating how to plug in μ and σ, compute α, then compute E[X+] and Var(X+) accordingly.

# Example parameters:
mu    <- 0
sigma <- 1

# Compute alpha = -mu / sigma
alpha <- -mu / sigma

# Standard normal PDF and CDF at alpha:
phiZ <- dnorm(-mu / sigma)  # PDF φ(z)
PhiZ <- pnorm(-mu / sigma)  # CDF Φ(z)

# 1) Compute E(X+):
E <- mu * (1 - PhiZ) + sigma * phiZ

# 2) Compute E[(X+)^2]:
E2 <- (mu^2 + sigma^2) * (1 - PhiZ) + mu * sigma * phiZ

# 3) Final variance:
Var <- E2 - E^2

# Print results:
cat("E[X*] =", E_plus, "\n")
cat("Var(X*) =", Var_plus, "\n")
