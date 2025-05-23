pflasso <- function(X = Data_X,
                    Y = Data_Y,
                    intercept = T) {
  n = length(X[, 1])
  p = length(X[1, ])
  lasso_fit <- cv.glmnet(
    x = X,
    y = Y,
    family = "gaussian",
    alpha = 1,
    intercept = intercept
  )
  beta_lasso <- as.vector(coef(lasso_fit, s = "lambda.min"))
  idx <- which(beta_lasso[-1] != 0)
  if (intercept == T) {
    X <- scale(X, scale = F)
    center = as.vector(attr(X, "scaled:center"))
    mean_Y = mean(Y)
    Y = Y - mean_Y
    beta <- rep(0, p + 1)
    if (length(idx) > 0) {
      Qlasso = l2_loss(X, Y, beta_lasso[-1])
      K = min(n - 1, length(idx))
      ordered_idx <- order(beta_lasso[-1], decreasing = T)
      idx_j = ordered_idx[1:K]
      ols_fit = as.vector(coef(glm(Y ~ X[, idx_j] - 1, family = gaussian())))
      gamma1 = (l2_loss(X[, idx_j], Y, ols_fit) - Qlasso) / 2
      for (k in 1:K) {
        idx_j = ordered_idx[1:k]
        ols_fit = as.vector(coef(glm(Y ~ X[, idx_j] - 1, family = gaussian())))
        beta[(1 + idx_j)] = ols_fit
        if ((l2_loss(X, Y, beta[-1]) - Qlasso) <= gamma1) {
          break
        }
      }
    }
    beta[1] = mean_Y - sum(beta[-1] * center)
    return(beta)
  } else{
    beta <- rep(0, p)
    if (length(idx) > 0) {
      Qlasso = l2_loss(X, Y, beta_lasso[-1])
      K = min(n - 1, length(idx))
      ordered_idx <- order(beta_lasso[-1], decreasing = T)
      idx_j = ordered_idx[1:K]
      ols_fit = as.vector(coef(glm(Y ~ X[, idx_j] - 1, family = gaussian())))
      gamma1 = (l2_loss(X[, idx_j], Y, ols_fit) - Qlasso) / 2
      for (k in 1:K) {
        idx_j = ordered_idx[1:k]
        ols_fit = as.vector(coef(glm(Y ~ X[, idx_j] - 1, family = gaussian())))
        beta[idx_j] = ols_fit
        if ((l2_loss(X, Y, beta) - Qlasso) <= gamma1) {
          break
        }
      }
    }
    return(beta)
  }
}

plasso <- function(X = Data_X,
                   Y = Data_Y,
                   family = "binomial",
                   intercept = T,
                   tau = NULL) {
  n = length(X[, 1])
  p = length(X[1, ])
  if (family == "binomial") {
    lasso_fit <- cv.glmnet(
      x = X,
      y = Y,
      family = "binomial",
      alpha = 1,
      intercept = intercept
    )
    beta_lasso <- as.vector(coef(lasso_fit, s = "lambda.min"))
    idx <- which(beta_lasso[-1] != 0)
    if (length(idx) > 1) {
      beta_fit <- as.vector(coef(
        cv.glmnet(
          x = X[, idx],
          y = Y,
          family = "binomial",
          alpha = 1,
          intercept = intercept
        )
      ))
      if (intercept == T) {
        beta = rep(0, p + 1)
        beta[c(1, (1 + idx))] = beta_fit
      } else{
        beta = rep(0, p)
        beta[idx] = beta_fit[-1]
      }
    } else if (length(idx) == 1) {
      if (intercept == T) {
        beta = rep(0, p + 1)
        beta_fit <- as.vector(coef(glm(Y ~ X[, idx], family = binomial())))
        beta[1] = beta_fit[1]
        beta[1 + idx] = beta_fit[-1]
      } else{
        beta = rep(0, p)
        beta_fit <- as.vector(coef(glm(Y ~ X[, idx] - 1, family = binomial())))
        beta[idx] = beta_fit
      }
    } else{
      beta = c(log(sum(Y == 1) / sum(Y == 0)), rep(0, p))
    }
    return(beta)
  } else if (family == "quantile") {
    fit <- rq.pen.cv(X, Y, penalty = "LASSO", tau = tau)
    ordered_index <- order(abs(as.vector(coef(fit, lambda = "min"))[-1]), decreasing =
                             T)
    p0 = min(sum(as.vector(coef(fit, lambda = "min"))[-1] != 0), n / log(n))
    beta = rep(0, p + 1)
    if (p0 > 0) {
      idx = ordered_index[1:p0]
      beta[c(1, (1 + idx))] = as.vector(coef(rq(Y ~ X[, idx], tau = tau)))
    } else{
      beta[1] = as.vector(coef(rq(Y ~ 1, tau = tau)))
    }
    return(beta)
    
  }
  
  
}
