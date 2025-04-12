###############################################################################
# Simulation data
###############################################################################
# 1. Poisson Emission Distribution
set.seed(929)

# Parameters
n <- 100
states <- 1:3
pi <- c(0.5, 0.3, 0.2)  # Initial state probabilities
A <- matrix(c(0.5, 0.3, 0.2,
              0.3, 0.5, 0.2,
              0.2, 0.3, 0.5), nrow = 3, byrow = TRUE)  # Transition matrix
lambda <- c(3, 6, 9)  # Poisson rates for each state

# Containers
Z <- numeric(n)              # hidden states
X <- numeric(n)              # observations

# Simulate first state and observation
Z[1] <- sample(states, 1, prob = pi)
X[1] <- rpois(1, lambda[Z[1]])

# Simulate the sequence
for (t in 2:n) {
  Z[t] <- sample(states, 1, prob = A[Z[t - 1], ])
  X[t] <- rpois(1, lambda[Z[t]])
}

# Combine into a data frame
hmm_data_poisson <- data.frame(time = 1:n, state = Z, observation = X)
###############################################################################
# 2. Gaussian Emission Distributions
set.seed(929)

# Parameters
n <- 150
states <- 1:3
pi <- c(0.5, 0.3, 0.2)
A <- matrix(c(0.4, 0.3, 0.3,
              0.3, 0.5, 0.2,
              0.4, 0.3, 0.3), nrow = 3, byrow = TRUE)
mu <- c(0, 5, 10)
sigma <- c(1, 2, 3)

# Containers
Z <- numeric(n)              # hidden states
X <- numeric(n)              # observations

# Simulate first state and observation
Z[1] <- sample(states, 1, prob = pi)
X[1] <- rnorm(1, mean = mu[Z[1]], sd = sigma[Z[1]])

# Simulate the sequence
for (t in 2:n) {
  Z[t] <- sample(states, 1, prob = A[Z[t - 1], ])
  X[t] <- rnorm(1, mean = mu[Z[t]], sd = sigma[Z[t]])
}

# Combine into a data frame
hmm_data_gaussian <- data.frame(time = 1:n, state = Z, observation = X)

###############################################################################
# Real-world data sets
###############################################################################
# 1. NVIDIA stock price
# Load the required package
library(quantmod)

# Get NVDA stock data from Yahoo Finance
NVDA_table <- getSymbols("NVDA", src = "yahoo", return.class = "xts", auto.assign = FALSE)

# Filter to a specific date range
NVDA_table <- stats::window(NVDA_table, start = as.Date("2024-01-02"), end = as.Date("2025-04-08"))

# Extract the closing prices
NVDA <- NVDA_table$NVDA.Close

nvda_vec <- as.numeric(NVDA)
nvda_ts <- ts(nvda_vec, start = c(2024, 1), frequency = length(nvda_vec))
###############################################################################
# 2. TESLA stock price

# Get TSLA stock data from Yahoo Finance
TSLA_table <- getSymbols("TSLA", src = "yahoo", return.class = "xts", auto.assign = FALSE)

# Filter to a specific date range
TSLA_table <- stats::window(TSLA_table, start = as.Date("2024-01-02"), end = as.Date("2025-04-08"))

# Extract the closing prices
TSLA <- TSLA_table$TSLA.Close

tesla_vec <- as.numeric(TSLA)
tesla_ts <- ts(tesla_vec, start = c(2024, 1), frequency = length(tesla_vec))
################################################################################
# 3. The Toronto temperature data
tem_2024 <- read.csv("data/2024.csv")
tem_2024 <- tem_2024 %>% 
  select(Date.Time, Year, Month, Day, Max.Temp...C., Min.Temp...C., 
         Mean.Temp...C., Total.Rain..mm., Total.Precip..mm., 
         Spd.of.Max.Gust..km.h.)

tem_2024$Date.Time <- as.Date(tem_2024$Date.Time)
april_oct <- na.omit(tem_2024) %>%
  filter(Date.Time >= as.Date("2024-04-01") & Date.Time <= as.Date("2024-10-31"))
plot(x = april_oct$Date.Time, 
     y = april_oct$Mean.Temp...C., 
     type = "l",
     main = "Daily Mean Temperature (Apr–Oct 2024)", 
     xlab = "Date", 
     ylab = "Mean Temp (°C)")

################################################################################
# 4. Earthquake data
dat <- read.table("http://www.hmms-for-time-series.de/second/data/earthquakes.txt")
# Or set your own path

x <- dat[, 2]
d <- dat[, 1]
n <- length(x)
################################################################################
# 5. Toronto Bicycle Theft
library(here)
library(tidyverse)
bicycle_thefts <- read.csv(here("./data/bicycle_thefts.csv"))

bicycle_data <- bicycle_thefts %>%
  group_by(OCC_YEAR, OCC_MONTH) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    OCC_MONTH = factor(OCC_MONTH, levels = month.name),
    MONTH_NUM = match(as.character(OCC_MONTH), month.name)
  ) %>%
  arrange(OCC_YEAR, MONTH_NUM) %>% 
  filter(OCC_YEAR >= 2014)

