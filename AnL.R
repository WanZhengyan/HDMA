
### sum1: relax the space of weight vector if sum1=F
### ds: a vector for choosing the model dimension of the candidate models
AnL<-function(X=Data_X,Y=Data_Y,ds,family="gaussian",sum1=F,intercept=T,sorted=T){
  if(intercept==T){
    X<-scale(X,scale=F)
    center=as.vector(attr(X,"scaled:center"))
    beta0=ifelse(family=="gaussian",mean(Y),log(sum(Y==1)/sum(Y==0)))
  }
  n=length(Y)
  p=length(X[1,])
  if(sorted==T){
    importance<-numeric(length=length(X[1,]))
    if(family=="gaussian"){
      family1=gaussian()
    }else{
      family1=binomial()
    }
    if(family=="gaussian"){
      importance<-as.vector(apply(X,2,function(vec){summary(glm(Y~vec,family = family1))$coefficients[2,4]}))
    }else if(family=="binomial"){
      importance<-as.vector(apply(X,2,function(vec){summary(glm(Y~vec,family = family1))$coefficients[2,4]}))
    }
    ordered_index<-order(importance, decreasing = FALSE)
  }else{
    ordered_index<-1:p
  }
  
  cv<-NULL
  weights_list<-list()
  beta_MA_list<-list()
  for(d in ds){
    K<-min(floor(sum(importance<=0.05)/d),floor(n/log(n)))
    if(sum(importance<=0.05)/d-K>0 && family=="gaussian"){
      pred_matrix<-matrix(0,nrow=0,ncol=K+1)
    }else{
      pred_matrix<-matrix(0,nrow=0,ncol=K)
    }
    for(m in 1:n){
      if(intercept==T){
        beta_matrix<-matrix(0,nrow=(1+length(X[1,])),ncol=0)
      }else{
        beta_matrix<-matrix(0,nrow=length(X[1,]),ncol=0)
      }
      for(i in 1:(K+1)){
        if(i==(K+1)){
          if(sum(importance<=0.05)/d-K>0 && family=="gaussian"){
            idx_j<-ordered_index[((i-1)*d+1):min(sum(importance<=0.05),i*d)]
          }else{
            break
          }
        }else{
          idx_j<-ordered_index[((i-1)*d+1):(i*d)]
        }
        if(intercept==F){
          beta_hat<-as.vector(coef(glm(Y[-m]~X[-m,idx_j]-1,family=family1)))
          beta<-rep(0,length(X[1,]))
          beta[idx_j]<-beta_hat
          beta_matrix<-cbind(beta_matrix,beta)
        }else{
          beta_hat<-as.vector(coef(glm(Y[-m]~X[-m,idx_j],family=family1)))
          beta<-rep(0,(length(X[1,])+1))
          beta[c(1,(idx_j+1))]<-c(beta0,beta_hat[-1])
          beta_matrix<-cbind(beta_matrix,beta)
        }
        
      }
      if(intercept==T){
        pred_matrix<-rbind(pred_matrix,cbind(1,t(X[m,]))%*%beta_matrix)
      }else{
        pred_matrix<-rbind(pred_matrix,t(X[m,])%*%beta_matrix)
      }
      
    }
    if(intercept==T){
      beta_matrix<-matrix(0,nrow=(1+length(X[1,])),ncol=0)
    }else{
      beta_matrix<-matrix(0,nrow=length(X[1,]),ncol=0)
    }
    for(i in 1:(K+1)){
      if(i==(K+1)){
        if(sum(importance<=0.05)/d-K>0 && family=="gaussian"){
          idx_j<-ordered_index[((i-1)*d+1):min(sum(importance<=0.05),i*d)]
        }else{
          break
        }
      }else{
        idx_j<-ordered_index[((i-1)*d+1):(i*d)]
      }
      if(intercept==F){
        beta_hat<-as.vector(coef(glm(Y~X[,idx_j]-1,family=family1)))
        beta<-rep(0,length(X[1,]))
        beta[idx_j]<-beta_hat
        beta_matrix<-cbind(beta_matrix,beta)
      }else{
        beta_hat<-as.vector(coef(glm(Y~X[,idx_j],family=family1)))
        beta<-rep(0,(length(X[1,])+1))
        beta[c(1,(idx_j+1))]<-c(beta0,beta_hat[-1])
        beta_matrix<-cbind(beta_matrix,beta)
      }
    }
    if(K==0){
      beta_MA<-as.vector(beta_matrix%*%c(1))
      cv<-c(cv,sum((Y-pred_matrix%*%c(1))^2))
    }else{
      if(sum(importance<=0.05)/d-K>0 && family=="gaussian"){
        w<-Variable(K+1)
      }else{
        w<-Variable(K)
      }
      if(family=="gaussian"){
        objective_w_expr<-cvxr_norm(pred_matrix%*%w-Y)^2
      }else if(family=="binomial"){
        objective_w_expr<-sum(logistic(pred_matrix%*%w)-Y*(pred_matrix%*%w))
      }
      objective_w<-Minimize(objective_w_expr)
      if(sum1==F){
        result<-solve(Problem(objective_w,constraints=list(w >= 0, w <= 1)))
      }else{
        result<-solve(Problem(objective_w,constraints=list(w >= 0, sum(w)==1)))
      }
      weights<-result$getValue(w)
      beta_MA<-as.vector(beta_matrix%*%weights)
      cv<-c(cv,result$value)
    }
    
    if(intercept==T){
      beta_MA[1]<-beta_MA[1]-sum(beta_MA[-1]*center)
    }
    
    beta_MA_list<-append(beta_MA_list,list(beta_MA))
    weights_list<-append(weights_list,list(weights))
  }
  didx_opt<-which.min(cv)
  beta_MA<-beta_MA_list[[didx_opt]]
  weights<-weights_list[[didx_opt]]

  return(list(beta_MA,weights,ds[didx_opt]))
}

