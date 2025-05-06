### projection onto simplex
proj_simplex <- function(v) {
  s = 0
  rho = 0
  U = v
  while (length(U) > 0) {
    G = U[which(U >= U[1])]
    L = U[which(U < U[1])]
    if ((s + sum(G)) - (rho + length(G)) * U[1] < 1) {
      s = s + sum(G)
      rho = rho + length(G)
      U = L
    } else{
      U = U[-1]
    }
  }
  v = v - ((s - 1) / rho)
  v[v < 0] = 0
  return(v)
}

### Fast greedy model averaging
FGMA <- function(pred_matrix,
                 Y,
                 K,
                 family,
                 L0 = 1e-6,
                 inflation = 2,
                 iter.max,
                 eps) {
  if (family == "gaussian") {
    loss <- l2_loss
    grad <- gradl2_loss
  } else if (family == "binomial") {
    loss <- logistic_loss
    grad <- gradlogistic_loss
  }
  n = length(pred_matrix[, 1])
  Aold = 1
  L = L0
  J <- 1
  Q <- Inf
  for (j in 1:K) {
    lam <- rep(0, K)
    lam[j] <- lam[j] + 1
    if (loss(X = pred_matrix, Y = Y, beta = lam) < Q) {
      J <- j
      Q <- loss(X = pred_matrix, Y = Y, beta = lam)
    }
  }
  wold <- rep(0, K)
  wold[J] <- 1
  z = wold
  obj_values = c(Q)
  for (k in 1:iter.max) {
    Q1 = 1
    Q2 = 0
    while (Q1 > Q2) {
      w = proj_simplex(z - as.vector(t(pred_matrix) %*% grad(
        X = pred_matrix, Y = Y, beta = z
      ) / L))
      Q1 = loss(X = pred_matrix, Y = Y, beta = w)
      Q21 = loss(X = pred_matrix, Y = Y, beta = z)
      Q22 = sum((w - z) * as.vector(t(pred_matrix) %*% grad(
        X = pred_matrix, Y = Y, beta = z
      )))
      Q23 = sum((w - z) ^ 2) * L / 2
      Q2 = Q21 + Q22 + Q23
      L = inflation * L
    }
    obj_values <- c(obj_values, Q1)
    L = L / inflation
    Anew = (1 + sqrt(1 + 4 * Aold ^ 2)) / 2
    z = w + (w - wold) * (Aold - 1) / Anew
    Aold = Anew
    wold = w
    
    if (min(as.vector(t(pred_matrix) %*% grad(
      X = pred_matrix, Y = Y, beta = w
    )) - sum(as.vector(
      t(pred_matrix) %*% grad(X = pred_matrix, Y = Y, beta = w)
    ) * w)) > -eps) {
      break
    }
  }
  return(list(weights = w, obj_values = obj_values))
  
}

### A function for finding a lambdamin for which the solution of the penalized logistic regression exists
find_lambda_min <- function(X,
                            Y,
                            lambdalist,
                            alpha,
                            standardize,
                            intercept,
                            penalty,
                            LLAiter) {
  if (penalty == "Lasso") {
    return(0)
  }
  p = length(X[1, ])
  n = length(X[, 1])
  if (penalty == "MCP") {
    dPen = dMCP
  } else if (penalty == "SCAD"){
    dPen = dSCAD
  }
  loss_vec <- rep(0, length(lambdalist))
  if (intercept == T) {
    beta_matrix <- matrix(0, nrow = (1 + p), ncol = 0)
  } else{
    beta_matrix <- matrix(0, nrow = p, ncol = 0)
  }
  for (lambda in lambdalist) {
    for (s in 1:(LLAiter + 1)) {
      if (s == 1) {
        beta1 <- rep(0, p + 1)
      }
      pf = dPen(beta1[-1] / lambda)
      if (sum(pf) == 0) {
        pf = rep(1, p)
        lambdaiter = lambda
      } else{
        lambdaiter = lambda * sum(pf) / p
        pf = p * pf / sum(pf)
      }
      beta1 <- as.vector(coef(
        glmnet(
          x = X,
          y = Y,
          family = "binomial",
          lambda = lambdaiter,
          penalty.factor = pf,
          alpha = alpha,
          standardize = standardize,
          intercept = intercept
        )
      ))
    }
    if (intercept == F) {
      beta_matrix <- cbind(beta_matrix, beta1[-1])
    } else{
      beta_matrix <- cbind(beta_matrix, beta1)
    }
  }
  if (intercept == F) {
    loss_vec <- apply(beta_matrix, 2, function(v) {
      min(1e6, logistic_loss(X, Y, v))
    })
  } else{
    loss_vec <- apply(beta_matrix, 2, function(v) {
      min(1e6, logistic_loss(cbind(1, X), Y, v))
    })
  }
  if (sum(loss_vec[-length(loss_vec)] / loss_vec[-1] + loss_vec[-1] /
          loss_vec[-length(loss_vec)] > 1e2) > 0) {
    lambdamin = lambdalist[which(loss_vec[-length(loss_vec)] / loss_vec[-1] +
                                   loss_vec[-1] / loss_vec[-length(loss_vec)] > 1e2)[1]]
  } else{
    lambdamin = 0
  }
  return(lambdamin)
}

### A function for feature ranking
Ranking <- function(X = X,
                    Y = Y,
                    Kne = 4,
                    d2 = 10,
                    tau = tau,
                    nfolds = nfolds,
                    family = "gaussian",
                    standardize = standardize,
                    intercept = F,
                    alpha = 1,
                    penalty = penalty,
                    LLAiter = LLAiter) {
  p = length(X[1, ])
  n = length(X[, 1])
  if (family == "gaussian") {
    family1 = gaussian()
  } else if (family == "binomial") {
    family1 = binomial()
  }
  if (penalty == "MCP") {
    dPen = dMCP
  } else{
    dPen = dSCAD
  }
  if (family == "gaussian") {
    loss <- l2_loss
  } else if (family == "binomial") {
    loss <- logistic_loss
  } else if (family == "quantile") {
    loss <- function(...) {
      quantile_loss(..., tau = tau)
    }
  }
  # Calculate initial estimator
  if (family == "quantile") {
    fit <- rq.pen.cv(X, Y, penalty = "LASSO", tau = tau)
    importance = abs(as.vector(coef(fit, lambda = "min")))[-1]
    lambda <- fit$gtr$lambda
  } else{
    fit <- cv.glmnet(
      x = X,
      y = Y,
      family = family,
      nfolds = nfolds,
      standardize = standardize,
      intercept = intercept,
      alpha = alpha
    )
    lambda = fit$lambda.min
    importance = abs(as.vector(coef(fit))[-1])
    lambdamin = 0
    if (penalty != "Lasso") {
      n1 = floor(n / nfolds)
      lambdalist = fit$lambda[which(fit$lambda >= fit$lambda.min)]
      if (family == "binomial") {
        lambdamin = find_lambda_min(
          X = X,
          Y = Y,
          lambdalist = lambdalist,
          alpha = alpha,
          standardize = standardize,
          intercept = intercept,
          penalty = penalty,
          LLAiter = LLAiter
        )
      }
      loss_vec <- rep(0, length(lambdalist))
      for (m in 1:nfolds) {
        if (m == nfolds) {
          idx_i <- ((m - 1) * n1 + 1):n
        } else{
          idx_i <- ((m - 1) * n1 + 1):(m * n1)
        }
        if (intercept == T) {
          beta_matrix <- matrix(0, nrow = (1 + length(X[1, ])), ncol = 0)
        } else{
          beta_matrix <- matrix(0, nrow = length(X[1, ]), ncol = 0)
        }
        for (lambda in lambdalist) {
          for (s in 1:(LLAiter + 1)) {
            if (s == 1) {
              beta1 <- rep(0, p + 1)
            }
            pf = dPen(beta1[-1] / lambda)
            if (sum(pf) == 0) {
              pf = rep(1, p)
              lambdaiter = lambda
            } else{
              lambdaiter = lambda * sum(pf) / p
              pf = p * pf / sum(pf)
            }
            beta1 <- as.vector(coef(
              glmnet(
                x = X[-idx_i, ],
                y = Y[-idx_i],
                family = family,
                lambda = lambdaiter,
                penalty.factor = pf,
                alpha = alpha,
                standardize = standardize,
                intercept = intercept
              )
            ))
          }
          if (intercept == F) {
            beta_matrix <- cbind(beta_matrix, beta1[-1])
          } else{
            beta_matrix <- cbind(beta_matrix, beta1)
          }
        }
        if (intercept == F) {
          loss_vec <- loss_vec + apply(beta_matrix, 2, function(v) {
            loss(X[idx_i, ], Y[idx_i], v)
          })
        } else{
          loss_vec <- loss_vec + apply(beta_matrix, 2, function(v) {
            loss(cbind(1, X[idx_i, ]), Y[idx_i], v)
          })
        }
      }
      lambda <- max(lambdalist[which.min(loss_vec)], lambdamin)
      for (s in 1:(LLAiter + 1)) {
        if (s == 1) {
          beta1 <- rep(0, p + 1)
        }
        pf = dPen(beta1[-1] / lambda)
        if (sum(pf) == 0) {
          pf = rep(1, p)
          lambdaiter = lambda
        } else{
          lambdaiter = lambda * sum(pf) / p
          pf = p * pf / sum(pf)
        }
        beta1 <- as.vector(coef(
          glmnet(
            x = X,
            y = Y,
            family = family,
            lambda = lambdaiter,
            penalty.factor = pf,
            alpha = alpha,
            standardize = standardize,
            intercept = intercept
          )
        ))
      }
      importance = abs(beta1[-1])
    }
  }
  zero <- which(importance == 0)
  # which(importance >0)
  # Ranking for the remaining predictors by marginal p-values
  if (family == "quantile") {
    xi_hat = as.numeric(quantile(Y, probs = tau))
    for (j in zero) {
      importance[j] <- mean(abs(X[, j] * as.numeric(coef(
        rq(Y ~ X[, j] - 1, tau = tau)
      )) - xi_hat))
      importance[zero] <- importance[zero] - max(importance)
    }
  } else{
    for (j in zero) {
      importance[j] <- -summary(glm(Y ~ X[, j] , family = family1))$coefficients[2, 4]
    }
  }
  ordered_index <- order(importance, decreasing = T)
  d1 = 2 * max(1, min(ceiling(sum(importance >= 0) / Kne)))
  Kne = min(Kne, floor(p / d1))
  K = floor((p - (Kne * d1)) / (d2)) + Kne
  lamcoef <- lambda * sqrt(n / log(d1 * Kne))
  return(
    list(
      ordered_index = ordered_index,
      d1 = d1,
      d2 = d2,
      Kne = Kne,
      K = K,
      importance = importance,
      lamcoef = lamcoef,
      lambdamin = lambdamin
    )
  )
  
}



### nest: "T" for nested models and "F" for nonnested models (dd is the gap between submodels);
### "mix" for nested + nonnested models (Kne is the number of nested models; d2 is the size of each nonnested model)
### d1: the gap of the nested models when nest = "mix" (default is data-driven);
### optimizer: "FGMA" for fast greedy model averaging; "GMA" for greedy model averaging (which only requires zero-order information); "CVX" for optimizer in "CVXR"; "relax" for relaxing the weight set and fitting by lasso
### d: use unpenalized estimators for the candidate models with model dimension not greater than d
### iter.max: the number of iterations in greedy model averaging algorithm; eps: tolerance in GMA\FGMA algorithm
### alpha: tuning parameter for the elastic net
### Kmax: the maximum number of candidate models
### family: "gaussian" for linear mean regression; "logistic" for logistic regression; "quantile" for linear quantile regression ("tau" is the probability we used)
### penalty: "Lasso", "MCP" or "SCAD" (using multi-step local linear approximation algorithm; LLAiter: iterations)
### threshold: hard threshold to MA estimator

HDMA <- function(X,
                 Y,
                 family = "gaussian",
                 nfolds = 5,
                 Kne = 4,
                 d1 = NULL,
                 d2 = 10,
                 d = 0,
                 Kmax = Inf,
                 nest = "mix",
                 dd = 10,
                 tau = NULL,
                 sorted = T,
                 iter.max = 200,
                 intercept = T,
                 penalty = "Lasso",
                 LLAiter = 3,
                 optimizer = "FGMA",
                 compare = F,
                 eps = 1e-5,
                 alpha = 1,
                 threshold = F) {
  p = length(X[1, ])
  n = length(X[, 1])
  n1 = floor(n / nfolds)
  if (intercept == T) {
    X <- scale(X, scale = F)
    center = as.vector(attr(X, "scaled:center"))
    beta0 = ifelse(family == "gaussian", mean(Y), log(sum(Y == 1) / sum(Y == 0)))
  }
  standardize = F
  if (penalty == "MCP") {
    dPen = dMCP
  } else{
    dPen = dSCAD
  }
  if (family == "quantile") {
    penalty = ifelse(penalty == "Lasso", "LASSO", penalty)
    optimizer = ifelse(optimizer == "FGMA", "GMA", optimizer)
  }
  ### Calculate variable importance
  Ranking_result <- Ranking(
    X = X,
    Y = Y,
    Kne = Kne,
    d2 = d2,
    tau = tau,
    family = family,
    intercept = intercept,
    standardize = standardize,
    alpha = alpha,
    penalty = penalty,
    LLAiter = LLAiter,
    nfolds = nfolds
  )
  if (sorted == T) {
    ordered_index <- Ranking_result$ordered_index
  } else{
    ordered_index <- 1:p
  }
  
  lamcoef <- Ranking_result$lamcoef
  d1 <- ifelse(is.null(d1), Ranking_result$d1, d1)
  Kne <- Ranking_result$Kne
  lambdamin <- Ranking_result$lambdamin
  K = ifelse(nest == "mix", min(Ranking_result$K, Kmax), min(floor(p / dd), Kmax))
  
  # Compute beta for each fold
  pred_matrix <- matrix(0, nrow = 0, ncol = K)
  for (m in 1:nfolds) {
    if (m == nfolds) {
      idx_i <- ((m - 1) * n1 + 1):n
    } else{
      idx_i <- ((m - 1) * n1 + 1):(m * n1)
    }
    if (intercept == T) {
      beta_matrix <- matrix(0, nrow = (1 + length(X[1, ])), ncol = 0)
    } else{
      beta_matrix <- matrix(0, nrow = length(X[1, ]), ncol = 0)
    }
    for (i in 1:K) {
      if (nest == T) {
        idx_j <- ordered_index[1:(i * dd)]
      } else if (nest == F) {
        idx_j <- ordered_index[((i - 1) * dd + 1):(i * dd)]
      } else if (nest == "mix") {
        if (i <= Kne) {
          idx_j <- ordered_index[1:(i * d1)]
        } else{
          idx_j <- ordered_index[c((Kne * d1 + (i - Kne - 1) * d2 + 1):(Kne * d1 + (i - Kne) * d2))]
        }
      }
      lambda = max(lambdamin, lamcoef * sqrt(log(length(idx_j)) / n))
      lambda = ifelse(length(idx_j) <= d, 0, lambda)
      if (family == "quantile") {
        beta1 <- as.vector(coef(
          rq.pen(
            X[-idx_i, idx_j],
            Y[-idx_i],
            tau = tau,
            penalty = penalty,
            lambda = c(lambda, lambda + 0.001)
          ),
          lambda = lambda
        ))
      } else{
        beta1 <- as.vector(coef(
          glmnet(
            x = X[-idx_i, idx_j],
            y = Y[-idx_i],
            family = family,
            lambda = lambda,
            alpha = alpha,
            standardize = standardize,
            intercept = intercept
          )
        ))
        if (penalty != "Lasso") {
          for (s in 1:LLAiter) {
            pf = dPen(beta1[-1] / lambda)
            if (sum(pf) == 0) {
              pf = rep(1, length(idx_j))
              lambdaiter = lambda
            } else{
              lambdaiter = lambda * sum(pf) / length(idx_j)
              pf = length(idx_j) * pf / sum(pf)
            }
            
            beta1 <- as.vector(coef(
              glmnet(
                x = X[-idx_i, idx_j],
                y = Y[-idx_i],
                family = family,
                lambda = lambdaiter,
                penalty.factor = pf,
                alpha = alpha,
                standardize = standardize,
                intercept = intercept
              )
            ))
          }
        }
        if (intercept == F) {
          beta1 <- beta1[-1]
        }
      }
      if (intercept == F) {
        beta <- rep(0, length(X[1, ]))
        beta[idx_j] <- beta1
        beta_matrix <- cbind(beta_matrix, beta)
      } else{
        beta <- rep(0, (length(X[1, ]) + 1))
        if (family == "quantile") {
          beta[c(1, (idx_j + 1))] <- beta1
        } else{
          beta[c(1, (idx_j + 1))] <- c(beta0, beta1[-1])
        }
        beta_matrix <- cbind(beta_matrix, beta)
        
      }
    }
    if (intercept == T) {
      pred_matrix <- rbind(pred_matrix, cbind(1, X[idx_i, ]) %*% beta_matrix)
    } else{
      pred_matrix <- rbind(pred_matrix, X[idx_i, ] %*% beta_matrix)
    }
  }
  
  # Calculate beta using all samples
  if (intercept == T) {
    beta_matrix <- matrix(0, nrow = (1 + length(X[1, ])), ncol = 0)
  } else{
    beta_matrix <- matrix(0, nrow = length(X[1, ]), ncol = 0)
  }
  for (i in 1:K) {
    if (nest == T) {
      idx_j <- ordered_index[1:(i * dd)]
    } else if (nest == F) {
      idx_j <- ordered_index[((i - 1) * dd + 1):(i * dd)]
    } else if (nest == "mix") {
      if (i <= Kne) {
        idx_j <- ordered_index[1:(i * d1)]
      } else{
        idx_j <- ordered_index[c((Kne * d1 + (i - Kne - 1) * d2 + 1):(Kne * d1 + (i - Kne) * d2))]
      }
    }
    lambda = max(lambdamin, lamcoef * sqrt(log(length(idx_j)) / n))
    lambda = ifelse(length(idx_j) <= d, 0, lambda)
    if (family == "quantile") {
      beta1 <- as.vector(coef(
        rq.pen(
          X[, idx_j],
          Y,
          tau = tau,
          penalty = penalty,
          lambda = c(lambda, lambda + 0.001)
        ),
        lambda = lambda
      ))
    } else{
      beta1 <- as.vector(coef(
        glmnet(
          x = X[, idx_j],
          y = Y,
          family = family,
          lambda = lambda,
          alpha = alpha,
          standardize = standardize,
          intercept = intercept
        )
      ))
      if (penalty != "Lasso") {
        for (s in 1:LLAiter) {
          pf = dPen(beta1[-1] / lambda)
          if (sum(pf) == 0) {
            pf = rep(1, length(idx_j))
            lambdaiter = lambda
          } else{
            lambdaiter = lambda * sum(pf) / length(idx_j)
            pf = length(idx_j) * pf / sum(pf)
          }
          beta1 <- as.vector(coef(
            glmnet(
              x = X[, idx_j],
              y = Y,
              family = family,
              lambda = lambdaiter,
              penalty.factor = pf,
              alpha = alpha,
              standardize = standardize,
              intercept = intercept
            )
          ))
        }
      }
      
      if (intercept == F) {
        beta1 <- beta1[-1]
      }
    }
    if (intercept == F) {
      beta <- rep(0, length(X[1, ]))
      beta[idx_j] <- beta1
      beta_matrix <- cbind(beta_matrix, beta)
    } else{
      beta <- rep(0, (length(X[1, ]) + 1))
      if (family == "quantile") {
        beta[c(1, (idx_j + 1))] <- beta1
      } else{
        beta[c(1, (idx_j + 1))] <- c(beta0, beta1[-1])
      }
      beta_matrix <- cbind(beta_matrix, beta)
    }
  }
  obj_values <- c()
  if (compare == T) {
    obj_values1 <- c()
    if (family == "gaussian") {
      loss <- l2_loss
    } else if (family == "binomial") {
      loss <- logistic_loss
    }
    w <- rep(0, K)
    wold <- w
    for (k in 1:iter.max) {
      alpha <- 2 / (k + 1)
      J <- 1
      Q <- Inf
      for (j in 1:K) {
        lam <- (1 - alpha) * w
        lam[j] <- lam[j] + alpha
        if (loss(X = pred_matrix,
                 Y = Y,
                 beta = lam) < Q) {
          J <- j
          Q <- loss(X = pred_matrix,
                    Y = Y,
                    beta = lam)
        }
      }
      obj_values1 <- c(obj_values1, Q)
      eJ <- rep(0, K)
      eJ[J] <- 1
      w <- w + alpha * (eJ - w)
      if ((max(abs(w - wold)) < eps) && (alpha < 0.01)) {
        break
      } else{
        wold <- w
      }
    }
    weights1 <- w
    beta_MA1 <- round(as.vector(beta_matrix %*% weights1), digits = 6)
    FGMA_fit <- FGMA(
      pred_matrix = pred_matrix,
      Y = Y,
      K = K,
      family = family,
      iter.max = iter.max,
      eps = eps
    )
    weights2 <- FGMA_fit$weights
    obj_values2 <- FGMA_fit$obj_values
    beta_MA2 <- round(as.vector(beta_matrix %*% weights2), digits = 6)
    return(
      list(
        beta_MA_GMA = beta_MA1,
        weights_GMA = weights1,
        obj_values_GMA = obj_values1,
        beta_MA_FGMA = beta_MA2,
        weights_FGMA = weights2,
        obj_values_FGMA = obj_values2
      )
    )
  }
  
  if (optimizer == "CVX") {
    w <- Variable(K)
    if (family == "gaussian") {
      objective_w_expr <- (cvxr_norm(pred_matrix %*% w - Y) ^ 2) / n
    } else if (family == "binomial") {
      objective_w_expr <- sum(logistic(pred_matrix %*% w) - Y * (pred_matrix %*% w)) / n
    } else if (family == "quantile") {
      objective_w_expr <- sum(0.5 * abs(Y - pred_matrix %*% w) + (tau - 0.5) * (Y -
                                                                                  pred_matrix %*% w))
    }
    objective_w <- Minimize(objective_w_expr)
    result <- solve(Problem(objective_w, constraints = list(w >= 0, sum(w) == 1)), solver = "ECOS")
    weights <- result$getValue(w)
    beta_MA <- round(as.vector(beta_matrix %*% weights), digits = 6)
  }
  if (optimizer == "relax") {
    if (family == "quantile") {
      w_fit = rq(Y ~ pred_matrix - 1, method = "lasso", tau = tau)
      weights <- as.vector(coef(w_fit))
      beta_MA <- round(as.vector(beta_matrix %*% weights), digits = 6)
    } else{
      w_fit <- cv.glmnet(
        x = pred_matrix,
        y = Y,
        family = family,
        alpha = 1,
        intercept = F
      )
      weights <- as.vector(coef(w_fit))[-1]
      beta_MA <- round(as.vector(beta_matrix %*% weights), digits = 6)
    }
    
  }
  
  if (optimizer == "GMA") {
    if (family == "gaussian") {
      loss <- l2_loss
    } else if (family == "binomial") {
      loss <- logistic_loss
    } else if (family == "quantile") {
      loss <- function(...) {
        quantile_loss(..., tau = tau)
      }
    }
    w <- rep(0, K)
    wold <- w
    for (k in 1:iter.max) {
      alpha <- 2 / (k + 1)
      J <- 1
      Q <- Inf
      for (j in 1:K) {
        lam <- (1 - alpha) * w
        lam[j] <- lam[j] + alpha
        if (loss(X = pred_matrix,
                 Y = Y,
                 beta = lam) < Q) {
          J <- j
          Q <- loss(X = pred_matrix,
                    Y = Y,
                    beta = lam)
        }
      }
      obj_values <- c(obj_values, Q)
      eJ <- rep(0, K)
      eJ[J] <- 1
      w <- w + alpha * (eJ - w)
      if ((max(abs(w - wold)) < eps) && (alpha < 5e-3)) {
        break
      } else{
        wold <- w
      }
    }
    weights <- w
    beta_MA <- round(as.vector(beta_matrix %*% weights), digits = 6)
  }
  
  if (optimizer == "FGMA") {
    FGMA_fit <- FGMA(
      pred_matrix = pred_matrix,
      Y = Y,
      K = K,
      family = family,
      iter.max = iter.max,
      eps = eps
    )
    weights <- FGMA_fit$weights
    obj_values <- FGMA_fit$obj_values
    beta_MA <- round(as.vector(beta_matrix %*% weights), digits = 6)
  }
  
  if (intercept == T) {
    if (threshold == T) {
      beta1 = beta_MA[-1]
      beta1[abs(beta1) < lamcoef * sqrt(log(p) / n)] = rep(0, sum(abs(beta1) <
                                                                    lamcoef * sqrt(log(p) / n)))
      beta_MA[-1] = beta1
    }
    beta_MA[1] <- beta_MA[1] - sum(beta_MA[-1] * center)
  } else{
    if (threshold == T) {
      beta_MA[abs(beta_MA) < lamcoef * sqrt(log(p) / n)] = rep(0, sum(abs(beta_MA) <
                                                                        lamcoef * sqrt(log(p) / n)))
    }
  }
  return(list(
    beta_MA = beta_MA,
    weights = weights,
    obj_values = obj_values,
    d1 = d1
  ))
}
