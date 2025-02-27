
### nest: "T" for nested models and "F" for nonnested models (dd is the gap between submodels);
### "mix" for nested + nonnested models (Kne is the number of nested models; d2 is the size of each nonnested model)
### d1: the gap of the nested models when nest = "mix" (Default is data-driven);
### optimizer: "GMA" for greedy model averaging; "CVX" for optimizer in "CVXR"; "relax" for relaxing the weight set and fitting by lasso
### dd: use unpenalized estimator for the candidate models with model dimension not greater than dd
### iter.max: the number of iterations in greedy model averaging algorithm; eps: tolerance in GMA algorithm
### weights: weights for each observation
### alpha: tuning parameter for the elastic net
### Kmax: the maximum number of candidate models
### weights: weights for each observation
MACV <- function(X = Data_X,
                 Y = Data_Y,
                 family = "gaussian",
                 nfolds = 10,
                 Kne = 2,
                 d1 = NULL,
                 d2 = 10,
                 d = 0,
                 Kmax = Inf,
                 nest = "mix",
                 dd = 10,
                 sorted = T,
                 iter.max = 500,
                 intercept = F,
                 optimizer = "GMA",
                 eps = 1e-3,
                 alpha = 1,
                 weights = NULL) {
  p = length(X[1, ])
  n = length(X[, 1])
  n1 = floor(n / nfolds)
  if(intercept==T){
    X<-scale(X,scale=F)
    center=as.vector(attr(X,"scaled:center"))
    beta0=ifelse(family=="gaussian",mean(Y),log(sum(Y==1)/sum(Y==0)))
  }
  standardize = ifelse(family=="gaussian",F,T)
    


  ### Calculate variable importance
  Ranking_result <- Ranking(
    X = X,
    Y = Y,
    Kne = Kne,
    d2 = d2,
    family = family,
    intercept = intercept,
    standardize = standardize,
    alpha = alpha,
    weights = weights
  )
  ordered_index <- Ranking_result$ordered_index
  lamcoef <- Ranking_result$lamcoef
  d1 <- ifelse(is.null(d1),Ranking_result$d1,d1)
  Kne <- Ranking_result$Kne
  K = ifelse(nest == "mix",min(Ranking_result$K, Kmax),min(floor(p / dd), Kmax))

  # Calculate beta for each fold
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
      lambda = lamcoef * sqrt(log(length(idx_j)) / n)
      lambda = ifelse(length(idx_j)<=d,0,lambda)
      if (family == "cox") {
        beta1 <- as.vector(coef(
          glmnet(
            x = X[-idx_i, idx_j],
            y = Surv(Y[-idx_i, 1], Y[-idx_i, 2]),
            family = family,
            lambda = lambda,
            alpha = alpha,
            weights = weights[-idx_i]
          )
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
            intercept = intercept,
            weights = weights[-idx_i]
          )
        ))
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
        beta[c(1, (idx_j + 1))] <- c(beta0, beta1[-1])
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
    lambda = lamcoef * sqrt(log(length(idx_j)) / n)
    lambda = ifelse(length(idx_j)<=d,0,lambda)
    if (family == "cox") {
      beta1 <- as.vector(coef(
        glmnet(
          x = X[, idx_j],
          y = Surv(Y[, 1], Y[, 2]),
          family = family,
          lambda = lambda,
          alpha = alpha,
          weights = weights
        )
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
          intercept = intercept,
          weights = weights
        )
      ))
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
      beta[c(1, (idx_j + 1))] <- c(beta0, beta1[-1])
      beta_matrix <- cbind(beta_matrix, beta)
    }
  }
  obj_values <- c()
  if (optimizer == "CVX") {
    if (family == "gaussian") {
      w <- Variable(K)
      objective_w_expr <- (cvxr_norm(pred_matrix %*% w - Y) ^ 2) / n
    } else if (family == "binomial") {
      w <- Variable(K)
      objective_w_expr <- sum(logistic(pred_matrix %*% w) - Y * (pred_matrix %*% w)) / n
    }
    objective_w <- Minimize(objective_w_expr)
    result <- solve(Problem(objective_w, constraints = list(w >= 0)))
    weights <- result$getValue(w)
    beta_MA <- round(as.vector(beta_matrix %*% weights), digits = 6)
  }
  if (optimizer == "relax") {
    if (family == "cox") {
      w_fit <- cv.glmnet(
        x = pred_matrix,
        y = Surv(Y[, 1], Y[, 2]),
        family = family,
        alpha = alpha,
      )
      weights <- abs(as.vector(coef(w_fit)))
    } else{
      w_fit <- cv.glmnet(
        x = pred_matrix,
        y = Y,
        family = family,
        alpha = 1,
        intercept = F
      )
      weights <- abs(as.vector(coef(w_fit)))[-1]
    }
    beta_MA <- round(as.vector(beta_matrix %*% weights), digits = 6)
  }
  if (optimizer == "GMA") {
    if (family == "gaussian") {
      loss <- l2_loss
    } else if (family == "binomial") {
      loss <- logistic_loss
    } else if (family == "cox") {
      loss <- cox_loss
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
        if (loss(X = pred_matrix, Y = Y, beta = lam) < Q) {
          J <- j
          Q <- loss(X = pred_matrix, Y = Y, beta = lam)
        }
      }
      obj_values <- c(obj_values, Q)
      eJ <- rep(0, K)
      eJ[J] <- 1
      w <- w + alpha * (eJ - w)
      if ((max(abs(w - wold)) < eps) && (alpha < 0.01)) {
        break
      } else{
        wold <- w
      }
    }
    weights <- w
    beta_MA <- round(as.vector(beta_matrix %*% weights), digits = 6)
  }
  if(intercept==T){
    beta_MA[1]<-beta_MA[1]-sum(beta_MA[-1]*center)
  }
  return(list(
    beta_MA = beta_MA,
    weights = weights,
    obj_values = obj_values,
    d1 = d1
  ))
}

### A function for feature ranking
Ranking <- function(X = X,
                    Y = Y,
                    Kne = 2,
                    d2 = 10,
                    family = "gaussian",
                    standardize = standardize,
                    intercept = F,
                    alpha = 1,
                    weights = NULL) {
  
  p = length(X[1, ])
  n = length(X[, 1])
  if (family == "gaussian") {
    family1 = gaussian()
  } else if (family == "binomial") {
    family1 = binomial()
  }

  # Calculate initial estimator
  if (family == "cox") {
    fit <- cv.glmnet(
      x = X,
      y = Surv(Y[, 1], Y[, 2]),
      family = family,
      standardize = F,
      alpha = alpha,
      weights = weights
    )
    importance = abs(as.vector(coef(fit)))
  } else{
    fit <- cv.glmnet(
      x = X,
      y = Y,
      family = family,
      standardize = standardize,
      intercept = intercept,
      alpha = alpha,
      weights = weights
    )
    importance = abs(as.vector(coef(fit))[-1])
  }
  lamcoef <- fit$lambda.min * sqrt(n / log(p))
  zero<-which(importance == 0)
  # which(importance >0)
  # Ranking for the remaining predictors by marginal p-values
  if (family == "cox") {
    for (j in zero) {
      importance[j] <- -summary(coxph(Surv(Y[, 1], Y[, 2]) ~ X[, j]))$coefficients[1, 4]
    }
  } else{
    for (j in zero) {
      importance[j] <- -summary(glm(Y ~ X[, j] , family = family1))$coefficients[2, 4]
    }
  }
  ordered_index <- order(importance, decreasing = T)
  d1 = max(2, min(ceiling(sum(importance > 0) / Kne )))
  Kne = min(Kne, floor(p / d1))
  K = floor((p - (Kne * d1)) / (d2)) + Kne
  return(
    list(
      ordered_index = ordered_index,
      d1 = d1,
      d2 = d2,
      Kne = Kne,
      K = K,
      importance = importance,
      lamcoef = lamcoef
    )
  )
  
}



logistic_loss<-function(X,Y,beta){
  mu<-X%*%beta
  loss<-log(1+exp(mu))-Y*mu
  return(mean(loss))
}


l2_loss<-function(X,Y,beta){
  mu<-X%*%beta
  loss<-(Y-mu)^2
  return(mean(loss)/2)
}

cox_loss<-function(X,Y,beta){
  n=length(X[,1])
  mu<-X%*%beta
  or<-order(Y[,1])
  Loss<-0
  for(i in 1:n){
    Loss<-Loss+Y[or[i],2]*(log(sum(exp(mu[or[i:n]])))-mu[or[i]])
  }
  return(as.numeric(Loss/n))
}
