library(MASS)
library(ncvreg)
library(caret)
library(glmnet)
library(ggplot2)
library(dplyr)
library(CVXR)
library(tidyr)
library(magrittr)
library(MESS)
library(lars)
library(parallel)
library(foreach)
library(doParallel)
rm(list=ls())


source("PMA.R")
source("MACV.R")
source("Fun.R")
source("AnL.R")
source("ALASSO.R")


packages<-c("MASS","dplyr","caret","ncvreg","glmnet","CVXR","ggplot2","tidyr","magrittr","MESS","lars","Matrix")

nlist=c(100,200)
plist=c(1000,2000)
settings=matrix(0,ncol=4,nrow=0)
for(betatype in 1:3){
  for(n in nlist){
    for(p in plist){
      for(Xtype in 1:2){
        settings=rbind(settings,c(betatype,n,p,Xtype))
      }
    }
  }
}

simulation=NULL

t1<-proc.time()
cl<-makeCluster(6)
registerDoParallel(cl)
simulation<-foreach(setting=1:length(settings[,1]),.combine="rbind",.packages=packages) %dopar% {
# for(setting in c(3)){
  
  set.seed(666)
  i=settings[setting,1]
  n=settings[setting,2]
  p=settings[setting,3]
  X_type=settings[setting,4]
  beta<-Create_beta(i,p)
  Sigma = covariance(p,0.5,type=X_type)
  sd=0.5
  
  #generate the test set
  n_test=1000
  Data_test<-Create_data_lin(n_test,p,Sigma,beta,sd=sd)
  test_X<-Data_test[,1:p]
  test_Y<-Data_test[,(p+1)]

  vec_lasso_test=c();vec_two_stage_lasso_test=c();vec_MCP_test=c();
  vec_EN_test=c();vec_alasso_test=c();vec_SCAD_test=c();vec_HDMA_test=c();
  vec_HDMArelax_test=c();vec_AnL_test=c();vec_PMA_test=c()
  for(round in 1:100){
    set.seed(123+round)
    #generate the training set
    Data_train<-Create_data_lin(n,p,Sigma,beta,sd=sd)
    Data_X<-Data_train[,1:p]
    Data_Y<-Data_train[,(p+1)]
    
    # lasso
    lasso_fit<-cv.glmnet(x = Data_X,y = Data_Y,family = "gaussian",alpha = 1,intercept=F)
    beta_lasso<-as.vector(coef(lasso_fit,s="lambda.min"))[-1]
    vec_lasso_test<-c(vec_lasso_test,l2_loss(test_X,test_Y,beta_lasso))

    # two stage lasso
    idx<-which(beta_lasso!=0)
    two_stage_lasso_fit<-cv.glmnet(x = Data_X[,idx],y = Data_Y,family = "gaussian",alpha = 1,intercept=F)
    two_stage_lasso_beta<-rep(0,length(Data_X[1,]))
    two_stage_lasso_beta[idx]<-as.vector(coef(two_stage_lasso_fit))[-1]
    vec_two_stage_lasso_test<-c(vec_two_stage_lasso_test,l2_loss(test_X,test_Y,two_stage_lasso_beta))

    
    # MCP
    MCP_fit<-cv.ncvreg(X = Data_X,y = Data_Y,family = "gaussian", penalty="MCP",intercept=F)
    beta_MCP<-as.vector(coef(MCP_fit,s="lambda.min"))[-1]
    vec_MCP_test<-c(vec_MCP_test,l2_loss(test_X,test_Y,beta_MCP))
    
    # SCAD
    SCAD_fit<-cv.ncvreg(X = Data_X,y = Data_Y,family = "gaussian", penalty="SCAD",intercept=F)
    beta_SCAD<-as.vector(coef(SCAD_fit,s="lambda.min"))[-1]
    vec_SCAD_test<-c(vec_SCAD_test,l2_loss(test_X,test_Y,beta_SCAD))
    
    # ALASSO
    beta_alasso<-alasso(X = Data_X,Y = Data_Y)
    vec_alasso_test<-c(vec_alasso_test,l2_loss(test_X,test_Y,beta_alasso))

    # EN
    EN_fit<-cv.glmnet(x = Data_X,y = Data_Y,family = "gaussian",alpha = 0.5,intercept=F)
    beta_EN<-as.vector(coef(EN_fit,s="lambda.min"))[-1]
    vec_EN_test<-c(vec_EN_test,l2_loss(test_X,test_Y,beta_EN))

    # HDMA
    HDMA_fit<-MACV(Data_X,Data_Y,Kne=2,d2=10,family="gaussian",nest="mix",intercept=F)
    beta_HDMA<-HDMA_fit[[1]]
    vec_HDMA_test<-c(vec_HDMA_test,l2_loss(test_X,test_Y,beta_HDMA))

    # HDMArelax
    HDMArelax_fit<-MACV(Data_X,Data_Y,Kne=2,d2=10,family="gaussian",nest="mix",optimizer="relax",intercept=F)
    beta_HDMArelax<-HDMArelax_fit[[1]]
    vec_HDMArelax_test<-c(vec_HDMArelax_test,l2_loss(test_X,test_Y,beta_HDMArelax))

    # AnL
    AnL_fit<-AnL(Data_X,Data_Y,ds=c(4,5,6),family="gaussian",intercept=F)
    beta_AnL<-AnL_fit[[1]]
    vec_AnL_test<-c(vec_AnL_test,l2_loss(test_X,test_Y,beta_AnL))
    
    #PMA
    PMA_fit<-PMA(Data_X,Data_Y,lambda=log(length(Data_Y)))
    beta_PMA<-PMA_fit[[1]]
    vec_PMA_test<-c(vec_PMA_test,l2_loss(test_X,test_Y,beta_PMA))
  }

  test <- data.frame(
    LASSO = vec_lasso_test,
    TwostageLasso = vec_two_stage_lasso_test,
    ALASSO = vec_alasso_test,
    MCP = vec_MCP_test,
    SCAD = vec_SCAD_test,
    EN = vec_EN_test,
    HDMA = vec_HDMA_test,
    HDMArelax = vec_HDMArelax_test,
    AnL = vec_AnL_test,
    PMA = vec_PMA_test
  )
  
  result<-cbind(data.frame(betatype=i,n=n,p=p,Xtype=X_type,measure=c("Mean","SD","Median")),rbind(apply(test,2,mean),apply(test,2,sd),apply(test,2,median)))

  
  data_long <- pivot_longer(test, cols = names(test), names_to = "Model", values_to = "value")
  data_long$Model <- factor(data_long$Model, levels = names(test))
  pp <- ggplot(data_long, aes(x = Model, y = value, fill = Model)) +
    geom_boxplot() +
    labs(title = paste0("LinReg (p=",p,", n=",n,", Case ",setting,")"), x = NULL, y = "Prediction error") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = Gen_morandi(n = nlevels(data_long$Model)))
    ggsave(pp,filename=paste0("LinReg_p",p,"_n",n,"_Case",setting,".png"))

    return(result)
}
print(proc.time()-t1)
stopImplicitCluster()
simulation=cbind(simulation[,1:5],round(simulation[,-(1:5)],3))
write.csv(simulation,file="simulation_Lin.csv",row.names=F)
