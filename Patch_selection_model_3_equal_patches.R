# Function to input parameters
input_parameters <- function() {
  capacity <- 10
  x_critical <- 3
  alpha <- 1
  horizon <- 20
  n_patch <- 3
  
  beta <- c(0, 0.004, 0.02)
  lambda <- c(1, 0.4, 0.6)
  y <- c(0, 3, 5)
  
  list(capacity = capacity, x_critical = x_critical, alpha = alpha, horizon = horizon, n_patch = n_patch,
       beta = beta, lambda = lambda, y = y)
}

# Function to initialize F1
initialize_f1 <- function(capacity, x_critical) {
  f1 <- rep(0, capacity + 1)
  f1[(x_critical + 2):(capacity + 1)] <- 1
  f1
}

# Function to solve DPE
solve_dpe <- function(f1, params) {
  capacity <- params$capacity
  x_critical <- params$x_critical
  alpha <- params$alpha
  n_patch <- params$n_patch
  beta <- params$beta
  lambda <- params$lambda
  y <- params$y
  
  f0 <- rep(0, capacity + 1)
  pstar <- rep(NA, capacity + 1)
  
  for (x in (x_critical + 1):capacity) {
    rhs <- rep(0, n_patch)
    
    for (i in 1:n_patch) {
      x_prime <- x - alpha + y[i]
      if (x_prime > capacity) x_prime <- capacity
      if (x_prime < x_critical) x_prime <- x_critical
      x2 <- x - alpha
      if (x2 < x_critical) x2 <- x_critical
      term1 <- lambda[i] * f1[x_prime + 1]
      term2 <- (1 - lambda[i]) * f1[x2 + 1]
      rhs[i] <- (1 - beta[i]) * (term1 + term2)
    }
    
    vmax <- max(rhs)
    imax <- which.max(rhs)
    
    f0[x + 1] <- vmax
    pstar[x + 1] <- imax
  }
  
  list(f0 = f0, pstar = pstar)
}

# Function to print results
print_results <- function(t, f1, pstar, capacity, x_critical) {
  cat("TIME: ", t, "\n")
  cat(" X F(X,t,T) Opt Patch\n")
  cat(" ----------------------------------------------------------- \n")
  for (x in x_critical:capacity) {
    cat(x, f1[x + 1], pstar[x + 1], "\n")
  }
  cat(" ----------------------------------------------------------- \n")
}

# Main program
params <- input_parameters()
f1 <- initialize_f1(params$capacity, params$x_critical)

for (t in (params$horizon - 1):1) {
  result <- solve_dpe(f1, params)
  f1 <- result$f0
  pstar <- result$pstar
  print_results(t, f1, pstar, params$capacity, params$x_critical)
}
