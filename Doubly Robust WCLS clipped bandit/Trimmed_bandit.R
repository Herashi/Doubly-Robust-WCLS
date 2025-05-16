set.seed(42)
T <- 1000

# True parameters for data generation
true_beta <- c(intercept = 0, A = -0.2, S = 0.8, AS = 0.2)
phi <- sqrt(0.5)       # AR(1) coefficient
sigma <- 0.5      # Error SD
eta <- 2          # Exploration temperature

# Bandit initialization
beta_hat <- matrix(0, ncol = 1, nrow = 4)  # Coefficients: [Intercept, A, S, A:S]
P <- diag(100, 4)                          # Prior precision matrix
epsilon_prev <- 0                           # AR error state

# Storage
history <- data.frame(
  S = numeric(T),
  A = numeric(T),
  Y = numeric(T),
  p = numeric(T)
)

for (t in 1:T) {
  # Generate context
  S_t <- sample(c(-1, 1), 1)
  
  # Calculate action probability (trimmed)
  if(t == 1) {
    delta <- 0
  } else {
    delta <- beta_hat[2] + beta_hat[4] * S_t  # A + A:S effect
  }
  p_t <- plogis(eta * delta)
  p_t <- pmax(pmin(p_t, 0.99), 0.01)
  
  # Select action
  A_t <- rbinom(1, 1, p_t)
  
  # Generate outcome with AR(1) errors
  epsilon_t <- phi * epsilon_prev + rnorm(1, 0, sigma)
  Y_t <- true_beta %*% c(1, A_t, S_t, A_t * S_t) + epsilon_t
  
  # Update RLS estimates
  x_t <- matrix(c(1, A_t, S_t, A_t * S_t), ncol = 1)
  y_hat <- crossprod(x_t, beta_hat)[1,1]
  K <- (P %*% x_t) / as.numeric(1 + crossprod(x_t, P) %*% x_t)
  beta_hat <<- beta_hat +  K  %*% (Y_t - y_hat)
  P <<- P - tcrossprod(K, x_t) %*% P
  
  # Store results and update state
  history[t,] <- list(S_t, A_t, Y_t, p_t)
  epsilon_prev <- epsilon_t
}

# Results inspection
par(mfrow = c(2,2))
plot(history$p, type = "l", main = "Action Probabilities", ylim = c(0,1))
acf(history$Y, main = "Outcome Autocorrelation")
plot(history$Y, type = "l", main = "Outcome Series")
cat("Estimated coefficients:\n")
print(setNames(round(beta_hat, 3), c("Intercept", "A", "S", "A:S")))