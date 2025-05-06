library(MASS)
library(caret)
library(glmnet)
library(ncvreg)
library(dplyr)
library(CVXR)
library(tidyr)
library(magrittr)
library(MESS)
library(flare)
library(parallel)
library(foreach)
library(doParallel)
rm(list=ls())


source("HDMA.R")
source("Fun.R")
source("Bootstrap.R")

nlist=c(100,200)
plist=c(100,200)
settings=matrix(0,ncol=3,nrow=0)

for(n in nlist){
  for(p in plist){
    for(Xtype in 1:2){
      settings=rbind(settings,c(n,p,Xtype))
    }
  }
}


packages<-c("MASS","dplyr","caret","ncvreg","glmnet","tidyr","magrittr","MESS","Matrix","flare")

simulation=NULL

t1<-proc.time()
cl<-makeCluster(8)
registerDoParallel(cl)
simulation<-foreach(setting=1:length(settings[,1]),.combine="rbind",.packages=packages) %dopar% {
  set.seed(123+setting)
  n=settings[setting,1]
  p=settings[setting,2]
  X_type=settings[setting,3]
  beta<-Create_beta_inference(p)
  Sigma = covariance(p,0.5,type=X_type)

  vec_lasso_result=matrix(NA,ncol=6,nrow=0);vec_HDMA_lasso_result=matrix(NA,ncol=6,nrow=0);
  vec_SCAD_result=matrix(NA,ncol=6,nrow=0);vec_HDMA_SCAD_result=matrix(NA,ncol=6,nrow=0);
  vec_MCP_result=matrix(NA,ncol=6,nrow=0);vec_HDMA_MCP_result=matrix(NA,ncol=6,nrow=0);
  idxlist=list(c(1:5),c(1:floor(p/5)),c(1:p))
  for(round in 1:20){
    set.seed(666+setting*round)
    Data<-Create_data_log(n,p,Sigma,beta)[[1]]
    Data_X<-Data[,1:p]
    Data_Y<-Data[,(p+1)]
    
    # lasso
    lasso_fit<-cv.glmnet(x = Data_X,y = Data_Y,family = "binomial",alpha = 1,intercept=F)
    beta_lasso<-as.vector(coef(lasso_fit,s="lambda.min"))[-1]
    inference_lasso=Inference(Data_X,Data_Y,family="binomial",beta_hat=beta_lasso,beta_true=beta,alpha=0.05,B=500,idxlist=idxlist)
    vec_lasso_result=rbind(vec_lasso_result,c(inference_lasso$coverage,inference_lasso$average_length))
    print(apply(vec_lasso_result,2,mean))
    # SCAD
    SCAD_fit<-cv.ncvreg(X = Data_X,y = Data_Y,family = "binomial", penalty="SCAD",intercept=F)
    beta_SCAD<-as.vector(coef(SCAD_fit,s="lambda.min"))[-1]
    inference_SCAD=Inference(Data_X,Data_Y,family="binomial",beta_hat=beta_SCAD,beta_true=beta,alpha=0.05,B=500,idxlist=idxlist)
    vec_SCAD_result=rbind(vec_SCAD_result,c(inference_SCAD$coverage,inference_SCAD$average_length))
    print(apply(vec_SCAD_result,2,mean))
    # MCP
    MCP_fit<-cv.ncvreg(X = Data_X,y = Data_Y,family = "binomial", penalty="MCP",intercept=F)
    beta_MCP<-as.vector(coef(MCP_fit,s="lambda.min"))[-1]
    inference_MCP=Inference(Data_X,Data_Y,family="binomial",beta_hat=beta_MCP,beta_true=beta,alpha=0.05,B=500,idxlist=idxlist)
    vec_MCP_result=rbind(vec_MCP_result,c(inference_MCP$coverage,inference_MCP$average_length))
    print(apply(vec_MCP_result,2,mean))
    # HDMA-lasso
    HDMA_lasso_fit<-HDMA(Data_X,Data_Y,family="binomial",nest="mix",intercept=F,penalty="Lasso")
    beta_HDMA_lasso<-HDMA_lasso_fit[[1]]
    inference_HDMA_lasso=Inference(Data_X,Data_Y,family="binomial",beta_hat=beta_HDMA_lasso,beta_true=beta,alpha=0.05,B=500,idxlist=idxlist)
    vec_HDMA_lasso_result=rbind(vec_HDMA_lasso_result,c(inference_HDMA_lasso$coverage,inference_HDMA_lasso$average_length))
    print(apply(vec_HDMA_lasso_result,2,mean))
    # HDMA-SCAD
    HDMA_SCAD_fit<-HDMA(Data_X,Data_Y,family="binomial",nest="mix",intercept=F,penalty="SCAD")
    beta_HDMA_SCAD<-HDMA_SCAD_fit[[1]]
    inference_HDMA_SCAD=Inference(Data_X,Data_Y,family="binomial",beta_hat=beta_HDMA_SCAD,beta_true=beta,alpha=0.05,B=500,idxlist=idxlist)
    vec_HDMA_SCAD_result=rbind(vec_HDMA_SCAD_result,c(inference_HDMA_SCAD$coverage,inference_HDMA_SCAD$average_length))
    print(apply(vec_HDMA_SCAD_result,2,mean))
    # HDMA-MCP
    HDMA_MCP_fit<-HDMA(Data_X,Data_Y,family="binomial",nest="mix",intercept=F,penalty="MCP")
    beta_HDMA_MCP<-HDMA_MCP_fit[[1]]
    inference_HDMA_MCP=Inference(Data_X,Data_Y,family="binomial",beta_hat=beta_HDMA_MCP,beta_true=beta,alpha=0.05,B=500,idxlist=idxlist)
    vec_HDMA_MCP_result=rbind(vec_HDMA_MCP_result,c(inference_HDMA_MCP$coverage,inference_HDMA_MCP$average_length))
    print(apply(vec_HDMA_MCP_result,2,mean))
    print(1)
  }
  result <- cbind(
    data.frame(
      n = n,
      p = p,
      Xtype = X_type,
      model = c("Lasso", "SCAD", "MCP", "HDMALasso", "HDMASCAD", "HDMAMCP")
    ),
    rbind(
      apply(vec_lasso_result, 2, mean),
      apply(vec_SCAD_result, 2, mean),
      apply(vec_MCP_result, 2, mean),
      apply(vec_HDMA_lasso_result, 2, mean),
      apply(vec_HDMA_SCAD_result, 2, mean),
      apply(vec_HDMA_MCP_result, 2, mean)
    )
  )
  return(result)
}
print(proc.time()-t1)
stopImplicitCluster()
simulation=cbind(simulation[,1:4],round(simulation[,-(1:4)],3))
write.csv(simulation,file="simulation_Inference_Log.csv",row.names=F)
