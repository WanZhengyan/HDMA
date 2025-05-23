library(MASS)
library(dplyr)
library(caret)
library(ncvreg)
library(glmnet)
library(CVXR)
library(ggplot2)
library(tidyr)
library(Matrix)
library(MESS)
library(magrittr)
library(parallel)
library(foreach)
library(doParallel)
rm(list=ls())

source("HDMA.R")
source("Fun.R")
source("AnL.R")
source("PLasso.R")


packages<-c("MASS","dplyr","caret","ncvreg","glmnet","CVXR","ggplot2","tidyr","magrittr","Matrix","MESS")

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
cl<-makeCluster(8)
registerDoParallel(cl)
simulation<-foreach(setting=1:length(settings[,1]),.combine="rbind",.packages=packages) %dopar% {
  cat("Processing setting", setting, ": betatype =", settings[setting,1], "n =", settings[setting,2], "p =", settings[setting,3], "Xtype =", settings[setting,4], "\n")
  set.seed(123+setting)
  i=settings[setting,1]
  n=settings[setting,2]
  p=settings[setting,3]
  X_type=settings[setting,4]
  beta<-Create_beta_log(i,p)
  Sigma = covariance(p,0.5,type=X_type)

  #generate the test set
  n_test=1000
  Data_test<-Create_data_log(n_test,p,Sigma,beta)[[1]]
  test_X<-Data_test[,1:p]
  test_Y<-Data_test[,(p+1)]
  
  vec_lasso_test=c();vec_MCP_test=c();vec_plasso_test=c()
  vec_EN_test=c();vec_SCAD_test=c();vec_HDMA_Lasso_test=c();
  vec_HDMA_SCAD_test=c();vec_HDMA_MCP_test=c();vec_AnL_test=c();

  for(round in 1:100){
    
    #generate the training set
    Data_train<-Create_data_log(n,p,Sigma,beta)[[1]]
    Data_X<-Data_train[,1:p]
    Data_Y<-Data_train[,(p+1)]

    # lasso
    lasso_fit<-cv.glmnet(x = Data_X,y = Data_Y,family = "binomial",alpha = 1,intercept=F)
    beta_lasso<-as.vector(coef(lasso_fit,s="lambda.min"))[-1]
    vec_lasso_test<-c(vec_lasso_test,logistic_loss(test_X,test_Y,beta_lasso))
    
    # post lasso
    beta_plasso<-plasso(X = Data_X,Y = Data_Y,family="binomial",intercept=F)
    vec_plasso_test<-c(vec_plasso_test,logistic_loss(test_X,test_Y,beta_plasso))
    
    # MCP
    MCP_fit<-cv.ncvreg(X = Data_X,y = Data_Y,family = "binomial", penalty="MCP",warn=F)
    beta_MCP<-as.vector(coef(MCP_fit,s="lambda.min"))[-1]
    vec_MCP_test<-c(vec_MCP_test,logistic_loss(test_X,test_Y,beta_MCP))

    # SCAD
    SCAD_fit<-cv.ncvreg(X = Data_X,y = Data_Y,family = "binomial", penalty="SCAD",warn=F)
    beta_SCAD<-as.vector(coef(SCAD_fit,s="lambda.min"))[-1]
    vec_SCAD_test<-c(vec_SCAD_test,logistic_loss(test_X,test_Y,beta_SCAD))

    #EN
    EN_fit<-cv.glmnet(x = Data_X,y = Data_Y,family = "binomial",alpha = 0.5,intercept=F)
    beta_EN<-as.vector(coef(EN_fit,s="lambda.min"))[-1]
    vec_EN_test<-c(vec_EN_test,logistic_loss(test_X,test_Y,beta_EN))
    
    # HDMA-Lasso
    HDMA_Lasso_fit<-HDMA(Data_X,Data_Y,Kne=4,d2=10,family="binomial",nest="mix",penalty="Lasso",intercept=F)
    beta_HDMA_Lasso<-HDMA_Lasso_fit[[1]]
    vec_HDMA_Lasso_test<-c(vec_HDMA_Lasso_test,logistic_loss(test_X,test_Y,beta_HDMA_Lasso))

    # HDMA-SCAD
    HDMA_SCAD_fit<-HDMA(Data_X,Data_Y,Kne=4,d2=10,family="binomial",nest="mix",penalty="SCAD",intercept=F)
    beta_HDMA_SCAD<-HDMA_SCAD_fit[[1]]
    vec_HDMA_SCAD_test<-c(vec_HDMA_SCAD_test,logistic_loss(test_X,test_Y,beta_HDMA_SCAD))

    # HDMA-MCP
    HDMA_MCP_fit<-HDMA(Data_X,Data_Y,Kne=4,d2=10,family="binomial",nest="mix",penalty="MCP",intercept=F)
    beta_HDMA_MCP<-HDMA_MCP_fit[[1]]
    vec_HDMA_MCP_test<-c(vec_HDMA_MCP_test,logistic_loss(test_X,test_Y,beta_HDMA_MCP))

    # AnL
    result <- try({
      AnL_fit<-AnL(Data_X,Data_Y,ds=c(2,3,4),family="binomial",intercept=F)
    }, silent = TRUE)
    if (inherits(result, "try-error")) {
      AnL_fit<-AnL(Data_X,Data_Y,ds=1,family="binomial",intercept=F)
    } else {}
    beta_AnL<-AnL_fit[[1]]
    vec_AnL_test<-c(vec_AnL_test,logistic_loss(test_X,test_Y,beta_AnL))

  }
  
  

  test <- data.frame(
    LASSO = vec_lasso_test,
    PLasso = vec_plasso_test,
    SCAD = vec_SCAD_test,
    MCP = vec_MCP_test,
    ENet = vec_EN_test,
    AnL = vec_AnL_test,
    HDMALasso = vec_HDMA_Lasso_test,
    HDMASCAD = vec_HDMA_SCAD_test,
    HDMAMCP = vec_HDMA_MCP_test
  )
  result<-cbind(data.frame(betatype=i,n=n,p=p,Xtype=X_type,measure=c("Mean","SD","Median")),rbind(apply(test,2,mean),apply(test,2,sd),apply(test,2,median)))

  data_long <- pivot_longer(test, cols = names(test), names_to = "Model", values_to = "value")
  data_long$Model <- factor(data_long$Model, levels = names(test))
  pp <- ggplot(data_long, aes(x = Model, y = value, fill = Model)) +
    geom_boxplot() +
    labs(title = paste0("LogReg (p=",p,", n=",n,", betatype ",settings[setting,1],", Xtype ",settings[setting,4],")"), x = NULL, y = "Prediction error") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_manual(values = Gen_morandi(n = nlevels(data_long$Model)))
  
  ggsave(pp,filename=paste0("LogReg_p",p,"_n",n,"_betatype",settings[setting,1],"_Xtype",settings[setting,4],".png"),width = 8,height = 6,dpi = 300)
  result=cbind(result[,1:5],round(result[,-(1:5)],3))
  write.csv(result,file=paste0("LogReg_p",p,"_n",n,"_betatype",settings[setting,1],"_Xtype",settings[setting,4],".csv"),row.names=F)
  return(result)
}
print(proc.time()-t1)
stopImplicitCluster()
simulation=cbind(simulation[,1:5],round(simulation[,-(1:5)],3))
write.csv(simulation,file="simulation_Log.csv",row.names=F)

