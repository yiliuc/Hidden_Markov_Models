###############################################################################
# Simulation studies
###############################################################################
set.seed(929)

# Parameters
n <- 10000
states <- 1:3
pi <- c(0.5, 0.3, 0.2)
A <- matrix(c(0.7, 0.2, 0.1,
              0.3, 0.5, 0.2,
              0.2, 0.3, 0.5), nrow = 3, byrow = TRUE)
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
hmm_data <- data.frame(time = 1:n, state = Z, observation = X)
###############################################################################
library(HiddenMarkov)

# Number of states
nstates <- 3

# Initialize the transition probabilities
Pi <- matrix(c(0.8, 0.1, 0.1,
               0.1, 0.8, 0.1,
               0.1, 0.1, 0.8), nrow=nstates, byrow=TRUE)

# Initialize the initial distribution
delta <- rep(1/nstates, nstates)

# Initial emission parameters
mu <- c(1, 2, 3)
sigma <- c(1, 1, 1)

model <- dthmm(hmm_data$observation, Pi=Pi, delta=delta, distn="norm",
               pm=list(mean=mu, sd=sigma))
fitted_model <- BaumWelch(model)
summary(fitted_model)
###############################################################################
# Application to real-world data sets
###############################################################################
# Load the required package
library(quantmod)

# Get NVDA stock data from Yahoo Finance
NVDA_table <- getSymbols("NVDA", src = "yahoo", return.class = "xts", auto.assign = FALSE)

# Filter to a specific date range
NVDA_table <- stats::window(NVDA_table, start = as.Date("2024-01-02"), end = as.Date("2025-04-08"))

# Extract the closing prices
NVDA <- NVDA_table$NVDA.Close

nvda_vec <- as.numeric(NVDA)
nvda_ts <- ts(nvda_vec, start = c(2024, 1), frequency = 252)
################################################################################
# The temperature data
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
################################################################################
################################################################################



################################################################################
# Textbook, using Poisson
################################################################################
pois.HMM.pn2pw <- function(m, lambda, gamma, delta = NULL, stationary = TRUE) {
  tlambda <- log(lambda)
  if (m == 1) return(tlambda)
  
  foo <- log(gamma / diag(gamma))
  tgamma <- as.vector(foo[!diag(m)])
  
  if (stationary) {
    tdelta <- NULL
  } else {
    tdelta <- log(delta[-1] / delta[1])
  }
  
  parvect <- c(tlambda, tgamma, tdelta)
  return(parvect)
}

pois.HMM.pw2pn <- function(m, parvect, stationary = TRUE) {
  lambda <- exp(parvect[1:m])
  gamma <- diag(m)
  
  if (m == 1) {
    return(list(lambda = lambda, gamma = gamma, delta = 1))
  }
  
  gamma[!gamma] <- exp(parvect[(m + 1):(m * m)])
  gamma <- gamma / apply(gamma, 1, sum)
  
  if (stationary) {
    delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  } else {
    foo <- c(1, exp(parvect[(m * m + 1):(m * m + m - 1)]))
    delta <- foo / sum(foo)
  }
  
  return(list(lambda = lambda, gamma = gamma, delta = delta))
}

pois.HMM.mllk <- function(parvect, x, m, stationary = TRUE, ...) {
  if (m == 1) return(-sum(dpois(x, exp(parvect), log = TRUE)))
  
  n <- length(x)
  pn <- pois.HMM.pw2pn(m, parvect, stationary = stationary)
  
  foo <- pn$delta * dpois(x[1], pn$lambda)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  
  for (i in 2:n) {
    if (!is.na(x[i])) {
      P <- dpois(x[i], pn$lambda)
    } else {
      P <- rep(1, m)
    }
    
    foo <- foo %*% pn$gamma * P
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
  }
  
  mllk <- -lscale
  return(mllk)
}


pois.HMM.mle <- function(x, m, lambda0, gamma0, delta0 = NULL, stationary = TRUE, ...) {
  parvect0 <- pois.HMM.pn2pw(m, lambda0, gamma0, delta0, stationary = stationary)
  
  mod <- nlm(pois.HMM.mllk, parvect0, x = x, m = m, stationary = stationary)
  
  pn <- pois.HMM.pw2pn(m = m, mod$estimate, stationary = stationary)
  
  mllk <- mod$minimum
  np <- length(parvect0)
  AIC <- 2 * (mllk + np)
  n <- sum(!is.na(x))
  BIC <- 2 * mllk + np * log(n)
  
  list(
    m = m,
    lambda = pn$lambda,
    gamma = pn$gamma,
    delta = pn$delta,
    code = mod$code,
    mllk = mllk,
    AIC = AIC,
    BIC = BIC
  )
}

pois.HMM.forecast <- function(xf, h = 1, x, mod)
{
  n    <- length(x)
  nxf  <- length(xf)
  dxf  <- matrix(0, nrow = h, ncol = nxf)
  foo  <- mod$delta * dpois(x[1], mod$lambda)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  
  for (i in 2:n)
  {
    foo <- foo %*% mod$gamma * dpois(x[i], mod$lambda)
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
  }
  
  for (i in 1:h)
  {
    foo <- foo %*% mod$gamma
    for (j in 1:mod$m)
      dxf[i, ] <- dxf[i, ] + foo[j] * dpois(xf, mod$lambda[j])
  }
  
  return(dxf)
}

################################################################################
# Read earthquake data
dat <- read.table("http://www.hmms-for-time-series.de/second/data/earthquakes.txt")
# Or set your own path

x <- dat[, 2]
d <- dat[, 1]
n <- length(x)

# ============================ fit 2-state HMM ============================
m <- 2
lambda0 <- c(15, 25) # intial parameters
gamma0 <- matrix(
  c(0.9, 0.1,
    0.1, 0.9),
  m, m, byrow = TRUE
)

mod2s <- pois.HMM.mle(x, m, lambda0, gamma0, stationary = TRUE)

delta0 <- c(1, 1) / 2 # Initial distribution
mod2h <- pois.HMM.mle(x, m, lambda0, gamma0, delta = delta0, stationary = FALSE)

mod2s
mod2h

# ============================ fit 3-state HMM ============================
m <- 3
lambda0 <- c(10, 20, 30)
gamma0 <- matrix(
  c(0.8, 0.1, 0.1,
    0.1, 0.8, 0.1,
    0.1, 0.1, 0.8),
  m, m, byrow = TRUE
)

mod3s <- pois.HMM.mle(x, m, lambda0, gamma0, stationary = TRUE)

delta0 <- c(1, 1, 1) / 3
mod3h <- pois.HMM.mle(x, m, lambda0, gamma0, delta = delta0, stationary = FALSE)

mod3s
mod3h
################################################################################
#=== Use it for 1 - step - ahead and plot the forecast d i s t r i b u t i o n .
h <- 1
xf <- 0:50
forecasts <- pois.HMM.forecast(xf, h, x, mod3s)
fc <- forecasts[1,]

par(mfrow = c(1, 1), las = 1)
plot(xf, fc, #fc,  
     type = "h",
     main = paste("Earthquake series: forecast distribution for", d[n] + 1),
     xlim = c(0, max(xf)), ylim = c(0, 0.12),
     xlab = "count", ylab = "probability", lwd = 3)

#=== Forecast 1 -4 steps ahead and plot these .
h <- 4
xf <- 0:45
forecasts <- pois.HMM.forecast(xf, h, x, mod3s)

par(mfrow = c(2, 2), las = 1)
for (i in 1:4) {
  fc <- forecasts[i,]
  plot(xf, fc, type = "h",
       main = paste("Forecast distribution for", d[n] + i),
       xlim = c(0, max(xf)), ylim = c(0, 0.12),
       xlab = "count", ylab = "probability", lwd = 3)
}

#=== Compute the marginal d i s t r i b u t i o n ( called " dstat " below )
# for mod3h .
#=== This is also the long - term forecast .
m <- 3
lambda <- mod3h$lambda
delta <- solve(t(diag(m) - mod3h$gamma + 1), rep(1, m))

dstat <- numeric(length(xf))
for (j in 1:m) {
  dstat <- dstat + delta[j] * dpois(xf, lambda[j])
}

#=== Compare the 30 - year - ahead forecast with the long - term forecast .
h <- 30
xf <- 0:45
forecasts <- pois.HMM.forecast(xf, h, x, mod3h)
fc <- forecasts[h,]

par(mfrow = c(1, 1), las = 1)
plot(xf, fc, type = "h",
     main = paste("Forecast distribution for", d[n] + h),
     xlim = c(0, max(xf)), ylim = c(0, 0.12),
     xlab = "count", ylab = "probability", lwd = 3)
lines(xf, dstat, col = "gray", lwd = 3)

#=== Compare the 50 - year - ahead forecast with the long - term forecast .
h <- 50
xf <- 0:45
forecasts <- pois.HMM.forecast(xf, h, x, mod3h)
fc <- forecasts[h,]

par(mfrow = c(1, 1), las = 1)
plot(xf, fc, type = "h",
     main = paste("Forecast distribution for", d[n] + h),
     xlim = c(0, max(xf)), ylim = c(0, 0.12),
     xlab = "count", ylab = "probability", lwd = 3)
lines(xf, dstat, col = "gray", lwd = 3)




################################################################################
################################################################################
################################################################################



################################################################################
# Gaussian
################################################################################
normal.HMM.pn2pw <- function(m, mu, sigma, gamma, delta = NULL, stationary = TRUE) {
  tmu <- mu                   # Means are unconstrained
  tsigma <- log(sigma)        # Log-transform standard deviations
  
  # Transition matrix transformation (same as Poisson HMM)
  foo <- log(gamma / diag(gamma))
  tgamma <- as.vector(foo[!diag(m)])
  
  if (stationary) {
    tdelta <- NULL
  } else {
    tdelta <- log(delta[-1] / delta[1])
  }
  
  parvect <- c(tmu, tsigma, tgamma, tdelta)
  return(parvect)
}

normal.HMM.pw2pn <- function(m, parvect, stationary = TRUE) {
  mu <- parvect[1:m]
  sigma <- exp(parvect[(m + 1):(2 * m)])
  
  gamma <- diag(m)
  gamma[!gamma] <- exp(parvect[(2 * m + 1):(2 * m + m * (m - 1))])
  gamma <- gamma / rowSums(gamma)
  
  if (stationary) {
    delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  } else {
    foo <- c(1, exp(parvect[(2 * m + m * (m - 1) + 1):(2 * m + m * (m - 1) + m - 1)]))
    delta <- foo / sum(foo)
  }
  
  return(list(mu = mu, sigma = sigma, gamma = gamma, delta = delta))
}

normal.HMM.mllk <- function(parvect, x, m, stationary = TRUE, ...) {
  if (m == 1) return(-sum(dnorm(x, mean = parvect[1], sd = exp(parvect[2]), log = TRUE)))
  
  n <- length(x)
  pn <- normal.HMM.pw2pn(m, parvect, stationary = stationary)
  
  # Initial probabilities
  foo <- pn$delta * dnorm(x[1], mean = pn$mu, sd = pn$sigma)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo # Initial probability
  
  for (i in 2:n) {
    if (!is.na(x[i])) {
      P <- dnorm(x[i], mean = pn$mu, sd = pn$sigma)
    } else {
      P <- rep(1, m)
    }
    
    foo <- foo %*% pn$gamma * P
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
  }
  
  mllk <- -lscale
  return(mllk)
}

normal.HMM.mle <- function(x, m, mu0, sigma0, gamma0, delta0 = NULL, stationary = TRUE, ...) {
  # Step 1: transform natural → working parameters
  parvect0 <- normal.HMM.pn2pw(m, mu0, sigma0, gamma0, delta0, stationary = stationary)
  
  # Step 2: minimize negative log-likelihood
  mod <- nlm(normal.HMM.mllk, parvect0, x = x, m = m, stationary = stationary)
  
  # Step 3: transform working → natural parameters
  pn <- normal.HMM.pw2pn(m = m, parvect = mod$estimate, stationary = stationary)
  
  # Step 4: model diagnostics
  mllk <- mod$minimum
  np <- length(parvect0)
  n <- sum(!is.na(x))
  
  AIC <- 2 * (mllk + np)
  BIC <- 2 * mllk + np * log(n)
  
  # Step 5: return results
  return(list(
    m = m,
    mu = pn$mu,
    sigma = pn$sigma,
    gamma = pn$gamma,
    delta = pn$delta,
    code = mod$code,
    mllk = mllk,
    AIC = AIC,
    BIC = BIC
  ))
}

normal.HMM.forecast <- function(xf, h = 1, x, mod) {
  n    <- length(x)
  nxf  <- length(xf)
  dxf  <- matrix(0, nrow = h, ncol = nxf)
  
  # Initial state probabilities after first observation
  foo <- mod$delta * dnorm(x[1], mean = mod$mu, sd = mod$sigma)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  
  # Filtering through the observed sequence
  for (i in 2:n) {
    foo <- foo %*% mod$gamma * dnorm(x[i], mean = mod$mu, sd = mod$sigma)
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
  }
  
  # Forecast h steps ahead
  for (i in 1:h) {
    foo <- foo %*% mod$gamma
    for (j in 1:mod$m) {
      dxf[i, ] <- dxf[i, ] + foo[j] * dnorm(xf, mean = mod$mu[j], sd = mod$sigma[j])
    }
  }
  
  return(dxf)
}

################################################################################
# ============================ fit 2-state Gaussian HMM ============================
m <- 2
mu0 <- c(10, 20)             # Initial means
sigma0 <- c(2, 2)            # Initial standard deviations
gamma0 <- matrix(
  c(0.9, 0.1,
    0.1, 0.9),
  m, m, byrow = TRUE
)

# Stationary model
mod2s <- normal.HMM.mle(x, m, mu0, sigma0, gamma0, stationary = TRUE)

# Non-stationary model
delta0 <- c(1, 1) / 2
mod2h <- normal.HMM.mle(x, m, mu0, sigma0, gamma0, delta0, stationary = FALSE)

# ============================ fit 3-state Gaussian HMM ============================
m <- 3
mu0 <- c(5, 15, 25)
sigma0 <- c(1.5, 2.5, 2)
gamma0 <- matrix(
  c(0.8, 0.1, 0.1,
    0.1, 0.8, 0.1,
    0.1, 0.1, 0.8),
  m, m, byrow = TRUE
)

# Stationary model
mod3s <- normal.HMM.mle(x, m, mu0, sigma0, gamma0, stationary = TRUE)

# Non-stationary model
delta0 <- c(1, 1, 1) / 3
mod3h <- normal.HMM.mle(x, m, mu0, sigma0, gamma0, delta0, stationary = FALSE)

# View the models
mod2s
mod2h
mod3s
mod3h




m <- 3

# Initialize the transition probabilities
gamma0 <- matrix(c(0.8, 0.1, 0.1,
               0.1, 0.8, 0.1,
               0.1, 0.1, 0.8), nrow=nstates, byrow=TRUE)

# Initialize the initial distribution
delta0 <- rep(1/nstates, nstates)

# Initial emission parameters
mu0 <- c(1, 2, 3)
sigma0 <- c(1, 1, 1)

x <- hmm_data[, 3]
