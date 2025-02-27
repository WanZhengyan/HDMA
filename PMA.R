
### lambda: tuning parameter for the penalty of model averaging
PMA<-function(X=Data_X,Y=Data_Y,lambda=2,intercept=F,sorted=F){
  p=length(X[1,])
  n=length(Y)
  colnames(X)<-1:p
  if(intercept==T){
    X<-scale(X,scale=F)
    center=as.vector(attr(X,"scaled:center"))
  }
  
  # ALASSO
  pf<-adaptive.weights(x=X,y=Y,weight.method="univariate")$weights
  trans_X<-X%*%diag((pf)^-1)
  lars_fit<-lars(x=X, y=Y, type = "lasso",normalize=F,intercept=intercept,use.Gram=F)
  solutions<-as.matrix(lars_fit$beta[2:(n+1),])
  last_enter<-as.vector(apply(solutions,2,fn_last_enter))
  
  # Prepare candidate models
  v<-c()
  BIC<-c()
  if(intercept==F){
    beta_matrix<-matrix(0,nrow=p,ncol=0)
    for(i in 1:n){
      v<-c(v,length(as.vector(which(solutions[i,]!=0))))
      lm_fit<-lm(Y~X[,as.vector(which(solutions[i,]!=0))]-1)
      beta<-rep(0,p)
      beta[as.vector(which(solutions[i,]!=0))]<-coef(lm_fit)
      beta_matrix<-cbind(beta_matrix,beta)
      BIC<-c(BIC,ifelse(sum(as.vector(which(solutions[i,]!=0)))==0,Inf,BIC(lm_fit)))
    }
  }else{
    beta_matrix<-matrix(0,nrow=p+1,ncol=0)
    for(i in 1:n){
      v<-c(v,length(as.vector(which(solutions[i,]!=0))))
      lm_fit<-lm(Y~X[,as.vector(which(solutions[i,]!=0))])
      beta<-rep(0,p+1)
      beta[c(1,as.vector(which(solutions[i,]!=0))+1)]<-coef(lm_fit)
      beta_matrix<-cbind(beta_matrix,beta)
      BIC<-c(BIC,BIC(lm_fit))
    }
  }
  
  BIC_model<-as.vector(which(solutions[which.min(BIC),]!=0))
  ordered_index<-BIC_model[order(last_enter[BIC_model], decreasing = FALSE)]
  is_in_solutions<-function(vec){
    for(i in 1:n){
      if(length(vec)==length(as.vector(which(solutions[i,]!=0)))){
        if(sum(sort(vec)!=sort(as.vector(which(solutions[i,]!=0))))==0){
          return(TRUE)
        }
      }
    }
    return(FALSE)
  }
  if(intercept==F){
    for(i in 1:(length(ordered_index)-1)){
      if(!is_in_solutions(ordered_index[1:i])){
        lm_fit<-lm(Y~X[,ordered_index[1:i]]-1)
        beta<-rep(0,p)
        beta[ordered_index[1:i]]<-coef(lm_fit)
        beta_matrix<-cbind(beta_matrix,beta)
        v<-c(v,i)
      }
    }
    lambda=lambda*sum((Y-X%*%beta_matrix[,which.min(BIC)])^2)/(n-(length(BIC_model)))
  }else{
    for(i in 1:(length(ordered_index)-1)){
      if(!is_in_solutions(ordered_index[1:i])){
        lm_fit<-lm(Y~X[,ordered_index[1:i]])
        beta<-rep(0,p+1)
        beta[c(1,ordered_index[1:i]+1)]<-coef(lm_fit)
        beta_matrix<-cbind(beta_matrix,beta)
        v<-c(v,i)
      }
    }
    lambda=lambda*sum((Y-cbind(1,X)%*%beta_matrix[,which.min(BIC)])^2)/(n-(length(BIC_model)))
  }
  

  # Optimize weights
  if(intercept==F){
    pred_matrix<-X%*%beta_matrix
  }else{
    pred_matrix<-cbind(1,X)%*%beta_matrix
  }
  K<-length(beta_matrix[1,])
  w<-Variable(K)
  objective_w_expr<-cvxr_norm(pred_matrix%*%w-Y)^2+lambda*sum(v*w)
  objective_w<-Minimize(objective_w_expr)
  result<-solve(Problem(objective_w,constraints=list(w >= 0, sum(w) == 1)))
  weights<-result$getValue(w)
  beta_MA<-round(as.vector(beta_matrix%*%weights),digits=6)
  if(intercept==T){
    beta_MA[1]<-beta_MA[1]-sum(beta_MA[-1]*center)
  }
  return(list(beta_MA,K,lambda))
}

fn_last_enter<-function(vec){
  if(sum(vec!=0)==0){
    return(Inf)
  }else{
    vec<-c(0,vec)
    return(max(which(vec[1:max(which(vec!=0))]==0)))
  }
}

