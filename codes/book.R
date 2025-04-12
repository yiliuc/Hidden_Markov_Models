################################################################################
# Codes of HMM from textbook, using Poisson
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

pois.HMM.lforward <- function(x, mod) {
  n <- length(x)
  lalpha <- matrix(NA, mod$m, n)
  
  foo <- mod$delta * dpois(x[1], mod$lambda)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  lalpha[, 1] <- lscale + log(foo)
  
  for (i in 2:n) {
    foo <- foo %*% mod$gamma * dpois(x[i], mod$lambda)
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
    lalpha[, i] <- log(foo) + lscale
  }
  
  return(lalpha)
}

pois.HMM.lbackward <- function(x, mod) {
  n <- length(x)
  m <- mod$m
  lbeta <- matrix(NA, m, n)
  
  lbeta[, n] <- rep(0, m)
  foo <- rep(1 / m, m)
  lscale <- log(m)
  
  for (i in (n - 1):1) {
    foo <- mod$gamma %*% (dpois(x[i + 1], mod$lambda) * foo)
    lbeta[, i] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo / sumfoo
    lscale <- lscale + log(sumfoo)
  }
  
  return(lbeta)
}

pois.HMM.state_probs <- function(x, mod) {
  n <- length(x)
  la <- pois.HMM.lforward(x, mod)
  lb <- pois.HMM.lbackward(x, mod)
  c <- max(la[, n])
  llk <- c + log(sum(exp(la[, n] - c)))
  stateprobs <- matrix(NA, ncol = n, nrow = mod$m)
  for (i in 1:n) {
    stateprobs[, i] <- exp(la[, i] + lb[, i] - llk)
  }
  return(stateprobs)
}

# Note that state output ‘statepreds’ is a matrix even if h=1.
pois.HMM.state_prediction <- function(h = 1, x, mod) {
  n <- length(x)
  la <- pois.HMM.lforward(x, mod)
  c <- max(la[, n])
  llk <- c + log(sum(exp(la[, n] - c)))
  
  statepreds <- matrix(NA, ncol = h, nrow = mod$m)
  foo <- exp(la[, n] - llk)
  
  for (i in 1:h) {
    foo <- foo %*% mod$gamma
    statepreds[, i] <- foo
  }
  
  return(statepreds)
}

pois.HMM.local_decoding <- function(x, mod) {
  n <- length(x)
  stateprobs <- pois.HMM.state_probs(x, mod)
  ild <- rep(NA, n)
  for (i in 1:n) {
    ild[i] <- which.max(stateprobs[, i])
  }
  ild
}

pois.HMM.viterbi <- function(x, mod) {
  n <- length(x)
  xi <- matrix(0, n, mod$m)
  
  # Initialization
  foo <- mod$delta * dpois(x[1], mod$lambda)
  xi[1, ] <- foo / sum(foo)
  
  # Forward recursion
  for (i in 2:n) {
    foo <- apply(xi[i - 1, ] * mod$gamma, 2, max) * dpois(x[i], mod$lambda)
    xi[i, ] <- foo / sum(foo)
  }
  
  # Backtracking
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ])
  
  for (i in (n - 1):1) {
    iv[i] <- which.max(mod$gamma[, iv[i + 1]] * xi[i, ])
  }
  
  return(iv)
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
xf <- 0:45
forecasts <- pois.HMM.forecast(xf, h, x, mod3s)
fc <- forecasts[1,]

par(mfrow = c(1, 1), las = 1)
plot(xf, fc, #fc,  
     type = "h",
     main = paste("Earthquake series: forecast distribution for", d[n] + 1),
     xlim = c(0, max(xf)), ylim = c(0, 0.12),
     xlab = "count", ylab = "probability", lwd = 3)
# lines(xf, dstat, col = "gray", lwd = 3)

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
state_probs <- pois.HMM.state_probs(x, mod3h)
