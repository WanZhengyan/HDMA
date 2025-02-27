alasso<-function(X=Data_X,Y=Data_Y,intercept=F){
  pf<-adaptive.weights(x=X,y=Y,weight.method="univariate")$weights
  alasso_fit<-cv.glmnet(x = X,y = Y,family = "gaussian", alpha=1,penalty.factor=pf,standardize=F,intercept=intercept)
  beta<-as.vector(coef(alasso_fit,s="lambda.min"))
  if(intercept==F){
    beta<-beta[-1]
  }
  return(beta)
}
