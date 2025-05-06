bootstrap_sapply<-function(B, X, W_hat, S_hat, idxlist){
  n=length(X[,1])
  p=length(X[1,])
  Z=abs(S_hat%*%W_hat%*%t(X)%*%matrix(rnorm(n*B, mean = 0, sd = 1),nrow=n,ncol=B))
  bootstraps=matrix(0,nrow=0,ncol=B)
  for(G in 1:length(idxlist)){
    bootstraps=rbind(bootstraps,apply(Z[idxlist[[G]],],2,max))
  }
  return(bootstraps)
}

Inference<-function(X,Y,family,beta_hat,beta_true,alpha=0.05,B=500,idxlist,nlambda=2){
  n=length(X[,1])
  p=length(X[1,])
  J_hat_inverse = sugm(diag(sqrt(L_pp(Y,X%*%beta_hat,family)))%*%X, verbose = F, nlambda=nlambda)
  W_hat = J_hat_inverse$icov[[J_hat_inverse$nlambda]]
  S_hat  =  diag(diag(t(diag(sqrt(L_pp(Y,X%*%beta_hat,family)))%*%X)%*%(diag(sqrt(L_pp(Y,X%*%beta_hat,family)))%*%X)))/n
  beta_tilde = beta_hat-as.vector(W_hat%*%apply(diag(L_p(Y,X%*%beta_hat,family))%*%X,2,mean))
  
  if(length(beta_true)==1){
    idxlist = list(c(1:p))
    bootstraps = bootstrap_sapply(B=B,X=diag(L_p(Y,X%*%beta_hat,family))%*%X/sqrt(n),W_hat=W_hat,S_hat=S_hat,idxlist=idxlist)
    Q = apply(bootstraps,1,quantile,1-alpha)
    upper = beta_tilde+(1/diag(S_hat))*Q/sqrt(n)
    lower = beta_tilde-(1/diag(S_hat))*Q/sqrt(n)
    return(list(upper=upper,lower=lower,beta_tilde=beta_tilde))
  }
  bootstraps = bootstrap_sapply(B=B,X=diag(L_p(Y,X%*%beta_hat,family))%*%X/sqrt(n),W_hat=W_hat,S_hat=S_hat,idxlist=idxlist)
  Q = apply(bootstraps,1,quantile,1-alpha)
  coverage = c()
  average_length = c()
  for(G in 1:length(idxlist)){
    coverage = c(coverage, all((beta_tilde[idxlist[[G]]]-(1/diag(S_hat))[idxlist[[G]]]*Q[G]/sqrt(n) < beta_true[idxlist[[G]]]) &
                                 (beta_tilde[idxlist[[G]]]+(1/diag(S_hat))[idxlist[[G]]]*Q[G]/sqrt(n) > beta_true[idxlist[[G]]])))
    average_length = c(average_length, mean((1/diag(S_hat))[idxlist[[G]]]*Q[G]/sqrt(n)*2))
  }
  return(list(coverage=coverage,average_length=average_length,beta_tilde=beta_tilde))
}

