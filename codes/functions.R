################################################################################
# Gaussian HMM functions
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
  } else if (stationary == FALSE) {
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
  alpha <- pn$delta * dnorm(x[1], mean = pn$mu, sd = pn$sigma)
  sum_alpha <- sum(alpha)
  lscale <- log(sum_alpha)
  alpha <- alpha / sum_alpha # Initial probability
  
  for (i in 2:n) {
    if (!is.na(x[i])) {
      P <- dnorm(x[i], mean = pn$mu, sd = pn$sigma)
    } else {
      P <- rep(1, m)
    }
    if (any(!is.finite(P))) {
      cat("Invalid P at t =", i, "; x[i] =", x[i], "\n")}
    
    alpha <- alpha %*% pn$gamma * P
    sum_alpha <- sum(alpha)
    # print(sum_alpha)
    
    # if (!is.finite(sum_alpha) || sum_alpha == 0) {
    #   cat("Invalid sum_alpha at t =", i, "; sum_alpha =", sum_alpha, "\n")}
      
    lscale <- lscale + log(sum_alpha)
    alpha <- alpha / sum_alpha
  }
  
  neg_log_likelihood <- -lscale
  return(neg_log_likelihood)
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

normal.HMM.lforward <- function(x, mod) {
  n <- length(x)
  lalpha <- matrix(NA, mod$m, n)
  
  # Initial step
  foo <- mod$delta * dnorm(x[1], mean = mod$mean, sd = mod$sd)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  lalpha[, 1] <- lscale + log(foo)
  
  # Forward loop
  for (i in 2:n) {
    foo <- foo %*% mod$gamma * dnorm(x[i], mean = mod$mean, sd = mod$sd)
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
    lalpha[, i] <- log(foo) + lscale
  }
  
  return(lalpha)
}

normal.HMM.lbackward <- function(x, mod) {
  n <- length(x)
  m <- mod$m
  lbeta <- matrix(NA, m, n)
  
  lbeta[, n] <- rep(0, m)
  foo <- rep(1 / m, m)
  lscale <- log(m)
  
  for (i in (n - 1):1) {
    foo <- mod$gamma %*% (dnorm(x[i + 1], mean = mod$mean, sd = mod$sd) * foo)
    lbeta[, i] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo / sumfoo
    lscale <- lscale + log(sumfoo)
  }
  
  return(lbeta)
}

normal.HMM.state_probs <- function(x, mod) {
  n <- length(x)
  la <- normal.HMM.lforward(x, mod)
  lb <- normal.HMM.lbackward(x, mod)
  c <- max(la[, n])
  llk <- c + log(sum(exp(la[, n] - c)))
  stateprobs <- matrix(NA, ncol = n, nrow = mod$m)
  for (i in 1:n) {
    stateprobs[, i] <- exp(la[, i] + lb[, i] - llk)
  }
  return(stateprobs)
}

# Note that state output ‘statepreds’ is a matrix even if h=1.
normal.HMM.state_prediction <- function(h = 1, x, mod) {
  n <- length(x)
  la <- normal.HMM.lforward(x, mod)
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

normal.HMM.local_decoding <- function(x, mod) {
  n <- length(x)
  stateprobs <- normal.HMM.state_probs(x, mod)
  ild <- rep(NA, n)
  for (i in 1:n) {
    ild[i] <- which.max(stateprobs[, i])
  }
  ild
}

normal.HMM.viterbi <- function(x, mod) {
  n <- length(x)
  xi <- matrix(0, n, mod$m)
  
  # Initialization (t = 1)
  foo <- mod$delta * dnorm(x[1], mean = mod$mean, sd = mod$sd)
  xi[1, ] <- foo / sum(foo)
  
  # Forward recursion (t = 2 to n)
  for (i in 2:n) {
    foo <- apply(xi[i - 1, ] * mod$gamma, 2, max) * dnorm(x[i], mean = mod$mean, sd = mod$sd)
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



normal.HMM.forecast <- function(xf, h = 1, x, mod) {
  n    <- length(x)
  nxf  <- length(xf)
  dxf  <- matrix(0, nrow = h, ncol = nxf)
  
  # Initial state probabilities after first observation
  foo <- mod$delta * dnorm(x[1], mean = mod$mean, sd = mod$sd)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  
  # Filtering through the observed sequence
  for (i in 2:n) {
    foo <- foo %*% mod$gamma * dnorm(x[i], mean = mod$mean, sd = mod$sd)
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
  }
  
  # Forecast h steps ahead
  for (i in 1:h) {
    foo <- foo %*% mod$gamma
    for (j in 1:mod$m) {
      dxf[i, ] <- dxf[i, ] + foo[j] * dnorm(xf, mean = mod$mean[j], sd = mod$sd[j])
    }
  }
  
  return(dxf)
}

################################################################################
# Transfer the format of outputs from HiddenMarkov to normal.HMM.mle
convert.HiddenMarkov.output <- function(hmm_result, x = NULL) {
  distn <- hmm_result$distn
  gamma <- hmm_result$Pi
  delta <- hmm_result$delta
  m <- nrow(gamma)
  
  # Initialize output variables
  mean <- sd <- lambda <- NULL
  code <- mllk <- AIC <- BIC <- NA
  
  # Extract emission parameters based on distribution
  if (distn == "norm") {
    mean <- hmm_result$pm$mean
    sd <- hmm_result$pm$sd
  } else if (distn == "pois") {
    lambda <- hmm_result$pm$lambda
  } else {
    stop("Unsupported distribution: ", distn)
  }
  
  # Compute log-likelihood and diagnostics if x is provided
  if (!is.null(x)) {
    if (distn == "norm") {
      parvect <- normal.HMM.pn2pw(m, mean, sd, gamma, delta, stationary = FALSE)
      mllk <- normal.HMM.mllk(parvect, x = x, m = m, stationary = FALSE)
    } else if (distn == "pois") {
      parvect <- pois.HMM.pn2pw(m, lambda, gamma, delta, stationary = FALSE)
      mllk <- pois.HMM.mllk(parvect, x = x, m = m, stationary = FALSE)
    }
    np <- length(parvect)
    n <- sum(!is.na(x))
    AIC <- 2 * (mllk + np)
    BIC <- 2 * mllk + np * log(n)
  }
  
  # Output
  return(list(
    m = m,
    distn = distn,
    mean = mean,
    sd = sd,
    lambda = lambda,
    gamma = gamma,
    delta = delta,
    code = code,
    mllk = mllk,
    AIC = AIC,
    BIC = BIC
  ))
}
