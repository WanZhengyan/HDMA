# Generate data for linear regression
Create_data_lin<-function(n,p,Sigma,beta,sd){
  x<-mvrnorm(n,rep(0,p),Sigma)
  mu=x%*%beta
  y<-rnorm(n, mean = 0, sd = sd)+mu
  Data<-as.matrix(cbind(x,y))
  return(Data)
}

# Generate data for logistic regression
Create_data_log<-function(n,p,Sigma,beta){
  x<-mvrnorm(n,rep(0,p),Sigma)
  linear_base=x%*%beta
  y0_prob<-1/(1+exp(-linear_base))
  y0<-rbinom(n=n,size=1.0,prob=y0_prob)
  Data<-as.matrix(cbind(x,y0))
  return(list(Data=Data,y0_prob=y0_prob))
}

# Generate regression coefficients
Create_beta<-function(type,p){
  if(type==1){
    beta<-c(rep(1,5),rep(0.2,10),rep(1,5),rep(0,p-20))
  }else if(type==2){
    beta<-5*c(seq(1,p,1)^(-2))
  }else if(type==3){
    beta<-5*c(exp(-0.2*seq(1,p,1)))
  }
  return(beta)
}

# Generate covariance matrix for regressors
covariance<-function(p,rho,type){
  if(type==1){
    M = matrix(rep(1:p,p),ncol=p,byrow=F)
    Sigma = rho^(abs(M-t(M)))
    return(Sigma)
  }else if(type==2){
    return(toeplitz(c(1,0.5,rep(0,p-2))))
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
l2_loss<-function(X,Y,beta){
  mu<-X%*%beta
  loss<-(Y-mu)^2
  return(mean(loss)/2)
}

logistic_loss<-function(X,Y,beta){
  mu<-X%*%beta
  loss<-log(1+exp(mu))-Y*mu
  return(mean(loss))
}


