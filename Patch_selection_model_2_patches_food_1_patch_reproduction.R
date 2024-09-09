# Define the parameters
x_crit <- 0      # Critical level of reserves (animal dies if reserves fall to this level)
x_max <- 30      # Maximum possible level of reserves
t_max <- 20      # Length of the season (days)
x_rep <- 4       # Level of reserves needed for reproduction

# Parameters for patches
a <- c(1, 1, 1)    # Energetic cost of foraging in patches 1, 2, and 3
m <- c(0.01, 0.05, 0.02)  # Chance of mortality if patch 1, 2, or 3 is visited
p <- c(0.2, 0.5, 0)       # Probability that food is found in patches 1, 2, and 3
y <- c(2, 4, 0)    # Energetic value of food in patches 1, 2, and 3
c <- 4             # Maximum reproductive output in the reproductive patch (patch 3) in a single day

# Parameters for g(x, T)
A <- 60            # Asymptotic value of g(x, T) when x becomes very large
x_0 <- 0.25 * x_max # x_0 is a parameter used in the function g(x, T)

# Define the function g(x, T)
g <- function(x) {
  # Calculates the terminal fitness value at time T based on reserve level x
  A * (x - x_crit) / (x - x_crit + x_0)
}

# Initialize matrices to store fitness values (F) and optimal patches (opt_patch)
F <- matrix(0, nrow = x_max + 1, ncol = t_max)
opt_patch <- matrix(0, nrow = x_max + 1, ncol = t_max)

# Set the terminal fitness values using g(x, T)
for (x in (x_crit + 1):x_max) {
  # For each reserve level from x_crit + 1 to x_max, calculate the fitness at time T
  F[x + 1, t_max] <- g(x)
}

# Backward iteration to calculate F(x, t) and opt_patch
for (t in (t_max - 1):1) {
  for (x in x_crit:x_max) {
    if (x == x_crit) {
      # If the reserve level is at the critical level, fitness is 0 (the animal is dead)
      F[x + 1, t] <- 0
      opt_patch[x + 1, t] <- 0
    } else {
      # Calculate the fitness value of visiting patch 1
      if (x - a[1] >= x_crit) {
        V1 <- (1 - m[1]) * (p[1] * F[min(x - a[1] + y[1], x_max) + 1, t + 1] +
                              (1 - p[1]) * F[max(x - a[1], x_crit) + 1, t + 1])
      } else {
        V1 <- 0
      }
      
      # Calculate the fitness value of visiting patch 2
      if (x - a[2] >= x_crit) {
        V2 <- (1 - m[2]) * (p[2] * F[min(x - a[2] + y[2], x_max) + 1, t + 1] +
                              (1 - p[2]) * F[max(x - a[2], x_crit) + 1, t + 1])
      } else {
        V2 <- 0
      }
      
      # Calculate the fitness value of visiting patch 3 (reproductive patch)
      if (x <= x_rep) {
        V3 <- (1 - m[3]) * F[max(x - a[3], x_crit) + 1, t + 1]
      } else if (x_rep < x & x <= x_rep + c) {
        V3 <- (x - x_rep) + (1 - m[3]) * F[max(x_rep - a[3], x_crit) + 1, t + 1]
      } else {
        V3 <- c + (1 - m[3]) * F[max(x - c - a[3], x_crit) + 1, t + 1]
      }
      
      # Find the maximum fitness value among the three patches
      V_values <- c(V1, V2, V3)
      F[x + 1, t] <- max(V_values)
      # Store the index of the patch that gives the maximum fitness value
      opt_patch[x + 1, t] <- which.max(V_values)
    }
  }
}

# Convert matrices to data frames for better readability
F_df <- as.data.frame(F)
opt_patch_df <- as.data.frame(opt_patch)

# Add reserve level as row names
row.names(F_df) <- row.names(opt_patch_df) <- 0:x_max

# Print the results
print("Maximized fitness values (F):")
print(F_df)
print("Optimal patches (opt_patch):")
print(opt_patch_df)
