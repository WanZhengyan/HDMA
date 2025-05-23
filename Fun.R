# Generate data for linear regression
Create_data_lin <- function(n, p, Sigma, beta, sd) {
  x <- mvrnorm(n, rep(0, p), Sigma)
  mu = x %*% beta
  y <- rnorm(n, mean = 0, sd = sd) + mu
  Data <- as.matrix(cbind(x, y))
  return(Data)
}

# Generate data for logistic regression
Create_data_log <- function(n, p, Sigma, beta) {
  x <- mvrnorm(n, rep(0, p), Sigma)
  linear_base = x %*% beta
  y0_prob <- 1 / (1 + exp(-linear_base))
  y0 <- rbinom(n = n, size = 1.0, prob = y0_prob)
  Data <- as.matrix(cbind(x, y0))
  return(list(Data = Data, y0_prob = y0_prob))
}

# Generate regression coefficients
Create_beta_lin <- function(type, p) {
  if (type == 1) {
    beta <- c(rep(1, 5), rep(0.2, 10), rep(1, 5), rep(0, p - 20))
  } else if (type == 2) {
    beta <- 5 * c(seq(1, p, 1)^(-2))
  } else if (type == 3) {
    beta <- 5 * c(exp(-0.3 * seq(1, p, 1)))
  }
  return(beta)
}

Create_beta_log <- function(type, p) {
  if (type == 1) {
    beta <- c(rep(3, 5), rep(1, 10), rep(-0.2, 5), rep(0, p - 20))
  } else if (type == 2) {
    beta <- 5 * c(rep(1, 5), c(seq(1, p - 5, 1)^(-4)))
  } else if (type == 3) {
    beta <- 5 * c(rep(1, 5), c(exp(-0.5 * seq(1, p - 5, 1))))
  }
  return(beta)
}
# Generate covariance matrix for regressors
covariance <- function(p, rho, type) {
  if (type == 1) {
    M = matrix(rep(1:p, p), ncol = p, byrow = F)
    Sigma = rho^(abs(M - t(M)))
    return(Sigma)
  } else if (type == 2) {
    return(toeplitz(c(1, rho, rep(0, p - 2))))
  }
}

# Color for the figures
Gen_morandi <- function(n) {
  hues <- seq(0, 1, length.out = n + 1)[1:n]
  saturation <- runif(n, 0.1, 0.3)
  value <- runif(n, 0.7, 0.9)
  hsv(h = hues, s = saturation, v = value)
}

### Loss functions
l2_loss <- function(X, Y, beta) {
  mu <- X %*% beta
  loss <- (Y - mu)^2
  return(mean(loss) / 2)
}

logistic_loss <- function(X, Y, beta) {
  mu <- X %*% beta
  loss <- log(1 + exp(mu)) - Y * mu
  return(mean(loss))
}

quantile_loss <- function(X, Y, beta, tau) {
  mu <- X %*% beta
  loss <- 0.5 * abs(Y - mu) + (tau - 0.5) * (Y - mu)
  return(mean(loss))
}

gradl2_loss <- function(X, Y, beta) {
  n = length(X[, 1])
  mu <- X %*% beta
  return((mu - Y) / n)
}

gradlogistic_loss <- function(X, Y, beta) {
  n = length(X[, 1])
  mu <- X %*% beta
  return(((1 / (1 + exp(
    -mu
  ))) - Y) / n)
}


Create_beta_inference <- function(p) {
  beta <- c(2, 0.5, 1, rep(0, p - 3))
  return(beta)
}


L_p <- function(Y, mu, family) {
  mu = as.vector(mu)
  if (family == "gaussian") {
    return(mu - Y)
  } else if (family == "binomial") {
    return(exp(mu) / (1 + exp(mu)) - Y)
  }
  
}

L_pp <- function(Y, mu, family) {
  mu = as.vector(mu)
  if (family == "gaussian") {
    return(rep(1, length(Y)))
  } else if (family == "binomial") {
    return(exp(mu) / ((1 + exp(mu))^2))
  }
  
}

# Some function for FGMA
dSCAD <- function(t, a = 3.7) {
  t = abs(t)
  return(sapply(t, function(x) {
    ifelse(x <= 1, 1, max(0, a - x) / (a - 1))
  }))
}

dMCP <- function(t, a = 3) {
  t = abs(t)
  return(sapply(t, function(x) {
    max(0, 1 - (x / a))
  }))
}
