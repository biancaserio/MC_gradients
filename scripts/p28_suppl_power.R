# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  

# Power analysis for Serio, Yilmaz et al. (2024)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  

### power analysis for whole model

# Load the pwr package
library(pwr)


# Define the effect size f2
f2 <- 0.1 # Effect size

# Perform power analysis
result <- pwr.f2.test(u = 2,            # Number of predictors (estradiol + progesterone; testosterone + cortisol; estradiol + testosterone)
                      f2 = f2,          # Effect size f2
                      sig.level = 0.05, # Significance level (uncorrected)
                      power = 0.8)      # Desired power

# Extract the required sample size (n)
required_n <- ceiling(result$v + result$u + 1)

# Print the result
cat("Required number of observations (n) for the desired power:", required_n, "\n")



### power analysis for the effect of 1 hormone in the model (via simulations)

# Load the packages
library(lme4)
library(simr)



# Set seed for reproducibility
set.seed(123)

# Define the number of observations
n <- 250  # Number of observations

# Generate simulated data
data <- data.frame(
  estradiol = scale(rnorm(n)),
  progesterone = scale(rnorm(n))
)

# Define the effect size for estradiol
effect_size <- 0.2  # Desired effect size

# Create a response variable with the specified effect size plus some noise
data$gradient_val_F <- effect_size * data$estradiol + effect_size * data$progesterone + rnorm(n)  # Adding random noise

# Fit the linear model
model <- lm(gradient_val_F ~ estradiol + progesterone, data = data)

# Perform power analysis
power <- powerSim(model, nsim = 10000, test = fixed("estradiol", method = "t"))

# Print the power result
print(power)

