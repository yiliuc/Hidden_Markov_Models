source("codes.R")
source("data_cleaning.R")
################################################################################
# Gaussian HMM implementations, Textbook example
################################################################################
# Read earthquake data
dat <- read.table("http://www.hmms-for-time-series.de/second/data/earthquakes.txt")
# Or set your own path

x <- dat[, 2]
d <- dat[, 1]
n <- length(x)

################################################################################
# ============================ fit 2-state Gaussian HMM ============================
m <- 2
mu0 <- c(15, 25)             # Initial means
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
mu0 <- c(10, 20, 30)
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
mod2hq
mod3s
mod3h

#=== Use it for 1 - step - ahead and plot the forecast d i s t r i b u t i o n .
h <- 1
xf <- 0:50
forecasts <- normal.HMM.forecast(xf, h, x, mod3s)
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
forecasts <- normal.HMM.forecast(xf, h, x, mod3s)

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
mu <- mod3h$mu
sigma <- mod3h$sigma

# Stationary distribution
delta <- solve(t(diag(m) - mod3h$gamma + 1), rep(1, m))

# Mixture density over grid xf
dstat <- numeric(length(xf))
for (j in 1:m) {
  dstat <- dstat + delta[j] * dnorm(xf, mean = mu[j], sd = sigma[j])
}


#=== Compare the 30 - year - ahead forecast with the long - term forecast .
h <- 30
xf <- 0:50
forecasts <- normal.HMM.forecast(xf, h, x, mod3h)
fc <- forecasts[h,]

par(mfrow = c(1, 1), las = 1)
plot(xf, fc, type = "h",
     main = paste("Forecast distribution for", d[n] + h),
     xlim = c(0, max(xf)), ylim = c(0, 0.12),
     xlab = "count", ylab = "probability", lwd = 3)
lines(xf, dstat, col = "gray", lwd = 3)

#=== Compare the 50 - year - ahead forecast with the long - term forecast .
h <- 50
xf <- 0:50
forecasts <- normal.HMM.forecast(xf, h, x, mod3h)
fc <- forecasts[h,]

par(mfrow = c(1, 1), las = 1)
plot(xf, fc, type = "h",
     main = paste("Forecast distribution for", d[n] + h),
     xlim = c(0, max(xf)), ylim = c(0, 0.12),
     xlab = "count", ylab = "probability", lwd = 3)
lines(xf, dstat, col = "gray", lwd = 3)



m <- 3

# Initialize the transition probabilities
gamma0 <- matrix(c(0.8, 0.1, 0.1,
                   0.1, 0.8, 0.1,
                   0.1, 0.1, 0.8), nrow=3, byrow=TRUE)

# Initialize the initial distribution
delta0 <- rep(1/3, 3)

# Initial emission parameters
mu0 <- c(1, 2, 3)
sigma0 <- c(1, 1, 1)

x <- hmm_data[, 3]


################################################################################
# NVIDIA example
################################################################################
library(forecast)
diff <- na.omit(diff(log(nvda_ts)))
x <- (as.numeric(diff))^2
ggtsdisplay(diff)

# ============================ fit 2-state Gaussian HMM ============================
m <- 2
mu0 <- c(-0.05, 0)             # Initial means
sigma0 <- c(0.2, 0.2)            # Initial standard deviations
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
mu0 <- c(10, 20, 30)
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

#=== Use it for 1 - step - ahead and plot the forecast d i s t r i b u t i o n .
h <- 1
xf <- seq(-0.1, 0.1, by = 0.01)
forecasts <- normal.HMM.forecast(xf, h, x, mod3h)
fc <- forecasts[1,]

par(mfrow = c(1, 1), las = 1)
plot(xf, fc, #fc,  
     type = "h",
     main = paste("Earthquake series: forecast distribution for", d[n] + 1),
     xlim = c(0, max(xf)), ylim = c(0, 20),
     xlab = "count", ylab = "probability", lwd = 3)

#=== Forecast 1 -4 steps ahead and plot these .
h <- 4
xf <- seq(-0.1, 0.1, by = 0.01)
forecasts <- normal.HMM.forecast(xf, h, x, mod3h)

par(mfrow = c(2, 2), las = 1)
for (i in 1:4) {
  fc <- forecasts[i,]
  plot(xf, fc, type = "h",
       main = paste("Forecast distribution for", d[n] + i),
       xlim = c(min(xf), max(xf)), ylim = c(0, 20),
       xlab = "count", ylab = "probability", lwd = 3)
}

#=== Compare the 50 - year - ahead forecast with the long - term forecast .
h <- 50
xf <- seq(-0.1, 0.1, by = 0.01)
forecasts <- normal.HMM.forecast(xf, h, x, mod3h)
fc <- forecasts[h,]

par(mfrow = c(1, 1), las = 1)
plot(xf, fc, type = "h",
     main = paste("Forecast distribution for", d[n] + h),
     xlim = c(min(xf), max(xf)), ylim = c(0, 20),
     xlab = "count", ylab = "probability", lwd = 3)




################################################################################
# NVIDIA example: EM + Forecast
################################################################################

################################################################################
library(HiddenMarkov)
#########################################
# Number of states
nstates <- 3

# Initialize the transition probabilities
Pi <- matrix(c(0.8, 0.1, 0.1,
               0.1, 0.8, 0.1,
               0.1, 0.1, 0.8), nrow=nstates, byrow=TRUE)

# Initialize the initial distribution
delta <- rep(1/nstates, nstates)

# Initial emission parameters
mu <- c(-0.05, 0, 0.05)
sigma <- c(0.1, 0.1, 0.1)

x <- (as.numeric(diff))^2

model <- dthmm(x, Pi=Pi, delta=delta, distn="norm",
               pm=list(mean=mu, sd=sigma))
fitted_model <- BaumWelch(model)
summary(fitted_model)
mod_new <- convert.HiddenMarkov.output(fitted_model, x = x)

#########################################
# Number of states
nstates <- 2

# Initialize the transition probabilities
Pi <- matrix(c(0.8, 0.1,
               0.1, 0.8), nrow=nstates, byrow=TRUE)

# Initialize the initial distribution
delta <- rep(1/nstates, nstates)

# Initial emission parameters
mu <- c(-0.05, 0.05)
sigma <- c(0.1, 0.1)

x <- (as.numeric(diff))^2

model <- dthmm(x, Pi=Pi, delta=delta, distn="norm",
               pm=list(mean=mu, sd=sigma))
fitted_model <- BaumWelch(model)
summary(fitted_model)


mod_new <- convert.HiddenMarkov.output(fitted_model, x = x)

#=== Use it for 1 - step - ahead and plot the forecast d i s t r i b u t i o n .

h <- 1
xf <- seq(-0.1, 0.1, length.out = 10000)  # finer grid
forecasts <- normal.HMM.forecast(xf, h, x, mod_new)
fc <- forecasts[h, ]

par(mfrow = c(1, 1), las = 1)
plot(xf, fc,
     type = "l",
     main = bquote("Forecast density of" ~ X[T+.(h)]),
     xlim = c(-0.005, 0.005),
     ylim = c(0, max(fc) * 1.05),
     xlab = expression(x),
     ylab = "Density",
     lwd = 2)

#=== Use it for 1:4 - step - ahead and plot the forecast d i s t r i b u t i o n .
h <- 4
xf <- seq(-0.1, 0.1, length.out = 10000)  # finer grid
forecasts <- normal.HMM.forecast(xf, h, x, mod_new)

par(mfrow = c(2, 2), las = 1, mar = c(4, 4, 2, 1))  # tidy layout

for (i in 1:h) {
  fc <- forecasts[i, ]
  
  plot(xf, fc,
       type = "l",
       main = bquote("Forecast density for " ~ X[T+.(i)]),
       xlim = c(-0.005, 0.005),
       ylim = c(0, max(fc) * 1.05),  # auto-scaled per panel
       xlab = expression(x),
       ylab = "Density",
       lwd = 2,
       col = "steelblue")
  
  # Optional: vertical line for reference (e.g., last observed value)
  abline(v = tail(x, 1), col = "red", lty = 2)
}

h <- 50
xf <- seq(-0.1, 0.1, length.out = 10000)  # finer grid
forecasts <- normal.HMM.forecast(xf, h, x, mod_new)
fc <- forecasts[h, ]
sum
par(mfrow = c(1, 1), las = 1)
plot(xf, fc,
     type = "l",
     main = bquote("Forecast density of" ~ X[T+.(h)]),
     xlim = c(-0.005, 0.005),
     ylim = c(0, max(fc) * 1.05),
     xlab = expression(x),
     ylab = "Density",
     lwd = 2)

################################################################################
# Confidence interval
h <- 50
xf <- seq(-0.1, 0.1, length.out = 10000)  # finer grid
forecasts <- normal.HMM.forecast(xf, h, x, mod_new)
# Setup
dx <- diff(xf)[1]  # grid spacing
H <- h             # total forecast steps
results <- data.frame(
  step = 1:H,
  mean = NA,
  lower_95 = NA,
  upper_95 = NA
)

# Loop over all forecast steps
for (i in 1:H) {
  fc <- forecasts[i, ]
  fc <- fc / sum(fc * dx)  # normalize
  
  # Mean
  mu <- sum(xf * fc * dx)
  
  # CDF and quantile function
  cdf <- cumsum(fc * dx)
  qfun <- approxfun(cdf, xf)
  
  # 95% confidence interval
  lwr <- qfun(0.025)
  upr <- qfun(0.975)
  
  # Store results
  results$mean[i] <- mu
  results$lower_95[i] <- lwr
  results$upper_95[i] <- upr
}

# Display the result
print(results)

