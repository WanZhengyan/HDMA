library(MASS)
library(ncvreg)
library(caret)
library(glmnet)
library(rqPen)
library(quantreg)
library(ggplot2)
library(dplyr)
library(CVXR)
library(tidyr)
library(flare)
library(tidyverse)
library(MESS)
library(magrittr)
library(lars)

source("PMA.R")
source("HDMA.R")
source("Fun.R")
source("AnL.R")
source("PLASSO.R")
source("Bootstrap.R")
source("JMA_quan.R")

Data <- t(read.csv("riboflavin.csv")[,-1])
Data_X <- Data[,-1]
Data_Y <- as.vector(Data[,1])

#############################################################################
############################## Mean regression ##############################
### Compare FGMA and GMA
for(penalty in c("Lasso","MCP","SCAD")){
  set.seed(123)
  MA_fit<-HDMA(Data_X,Data_Y,Kne=4,d2=10,family="gaussian",eps=0,compare=T,penalty=penalty)
  
  ggsave(
    ggplot() +
      geom_line(
        data = data.frame(
          iterations = 0:40,
          value = MA_fit[[3]][1:41],
          group = "GMA"
        ),
        aes(x = iterations, y = value, color = group),
        linewidth = 1.5
      ) +
      geom_point(  
        data = data.frame(
          iterations = 0:40,
          value = MA_fit[[3]][1:41],
          group = "GMA"
        ),
        aes(x = iterations, y = value, color = group),
        shape = 21,      
        size = 3,        
        fill = "white",  
        stroke = 1      
      ) +
      geom_line(
        data = data.frame(
          iterations = 0:40,
          value = MA_fit[[6]][1:41],
          group = "FGMA"
        ),
        aes(x = iterations, y = value, color = group),
        linewidth = 1.5
      ) +
      geom_point(  
        data = data.frame(
          iterations = 0:40,
          value = MA_fit[[6]][1:41],
          group = "FGMA"
        ),
        aes(x = iterations, y = value, color = group),
        shape = 24,    
        size = 3,        
        fill = "white",  
        stroke = 1    
      ) +
      scale_color_manual(
        name = "Algorithm",
        values = c("GMA" = "#00008B", "FGMA" = "#8B4000"),
        labels = c("GMA" = "GMA", "FGMA" = "FGMA")
      ) +
      labs(x = "Iterations", y = "CV(w)/n") +
      theme(
        axis.text.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 13),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = c(0.95, 0.95),  
        legend.justification = c(1, 1),
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(face = "bold", size = 15)
      ),
    filename = paste0("FGMA_curve_Lin_",penalty, ".png"),
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
  )
}

### Inference
set.seed(123)
scaled_X<-scale(Data_X)
scaled_Y<-as.vector(scale(Data_Y))
gene_name <- as.vector(read.csv("riboflavin.csv")[,1])[-1]
idx<-which(as.vector(coef(cv.glmnet(x = scaled_X,y = scaled_Y,family = "gaussian",alpha=1,intercept=F),s="lambda.min"))[-1]!=0)

# lasso
lasso_fit<-cv.glmnet(x = scaled_X[,idx],y = scaled_Y,family = "gaussian",alpha = 1,intercept=F)
beta_lasso<-as.vector(coef(lasso_fit,s="lambda.min"))
inference_lasso=Inference(scaled_X[,idx],scaled_Y,family="gaussian",beta_hat=beta_lasso[-1],beta_true=0,alpha=0.05,B=500,idxlist=0)
gene_name[idx[which(inference_lasso$upper<=0|inference_lasso$lower>=0)]]

# SCAD
SCAD_fit<-cv.ncvreg(X = scaled_X[,idx],y = scaled_Y,family = "gaussian", penalty="SCAD",intercept=F)
beta_SCAD<-as.vector(coef(SCAD_fit,s="lambda.min"))
inference_SCAD=Inference(scaled_X[,idx],scaled_Y,family="gaussian",beta_hat=beta_SCAD[-1],beta_true=0,alpha=0.05,B=500,idxlist=0)
gene_name[idx[which(inference_SCAD$upper<=0|inference_SCAD$lower>=0)]]

# MCP
MCP_fit<-cv.ncvreg(X = scaled_X[,idx],y = scaled_Y,family = "gaussian", penalty="MCP",intercept=F)
beta_MCP<-as.vector(coef(MCP_fit,s="lambda.min"))
inference_MCP=Inference(scaled_X[,idx],scaled_Y,family="gaussian",beta_hat=beta_MCP[-1],beta_true=0,alpha=0.05,B=500,idxlist=0)
gene_name[idx[which(inference_MCP$upper<=0|inference_MCP$lower>=0)]]

# HDMA-lasso
HDMA_lasso_fit<-HDMA(scaled_X[,idx],scaled_Y,family="gaussian",nest="mix",intercept=F,penalty="Lasso")
beta_HDMA_lasso<-HDMA_lasso_fit[[1]]
inference_HDMA_lasso=Inference(scaled_X[,idx],scaled_Y,family="gaussian",beta_hat=beta_HDMA_lasso,beta_true=0,alpha=0.05,B=500,idxlist=0)
gene_name[idx[which(inference_HDMA_lasso$upper<=0|inference_HDMA_lasso$lower>=0)]]

# HDMA-SCAD
HDMA_SCAD_fit<-HDMA(scaled_X[,idx],scaled_Y,family="gaussian",nest="mix",intercept=F,penalty="SCAD")
beta_HDMA_SCAD<-HDMA_SCAD_fit[[1]]
inference_HDMA_SCAD=Inference(scaled_X[,idx],scaled_Y,family="gaussian",beta_hat=beta_HDMA_SCAD,beta_true=0,alpha=0.05,B=500,idxlist=0)
gene_name[idx[which(inference_HDMA_SCAD$upper<=0|inference_HDMA_SCAD$lower>=0)]]

# HDMA-MCP
HDMA_MCP_fit<-HDMA(scaled_X[,idx],scaled_Y,family="gaussian",nest="mix",intercept=F,penalty="MCP")
beta_HDMA_MCP<-HDMA_MCP_fit[[1]]
inference_HDMA_MCP=Inference(scaled_X[,idx],scaled_Y,family="gaussian",beta_hat=beta_HDMA_MCP,beta_true=0,alpha=0.05,B=500,idxlist=0)
gene_name[idx[which(inference_HDMA_MCP$upper<=0|inference_HDMA_MCP$lower>=0)]]

df_wide <- data.frame(
  variable = gene_name[idx],
  estimate_lasso = inference_lasso$beta_tilde,
  lower_lasso = inference_lasso$lower,
  upper_lasso = inference_lasso$upper,
  estimate_SCAD = inference_SCAD$beta_tilde,
  lower_SCAD = inference_SCAD$lower,
  upper_SCAD = inference_SCAD$upper,
  estimate_MCP = inference_MCP$beta_tilde,
  lower_MCP = inference_MCP$lower,
  upper_MCP = inference_MCP$upper,
  estimate_HDMA_lasso = inference_HDMA_lasso$beta_tilde,
  lower_HDMA_lasso = inference_HDMA_lasso$lower,
  upper_HDMA_lasso = inference_HDMA_lasso$upper,
  estimate_HDMA_SCAD = inference_HDMA_SCAD$beta_tilde,
  lower_HDMA_SCAD = inference_HDMA_SCAD$lower,
  upper_HDMA_SCAD = inference_HDMA_SCAD$upper,
  estimate_HDMA_MCP = inference_HDMA_MCP$beta_tilde,
  lower_HDMA_MCP = inference_HDMA_MCP$lower,
  upper_HDMA_MCP = inference_HDMA_MCP$upper
)

df_long <- df_wide %>%
  mutate(variable = factor(variable, levels = unique(variable))) %>%
  pivot_longer(
    cols = -variable,
    names_to = c(".value", "model"),
    names_pattern = "([a-zA-Z]+)_([a-zA-Z_]+)"
  ) %>%
  mutate(
    model = case_when(
      model == "lasso" ~ "LASSO",
      model == "SCAD" ~ "SCAD",
      model == "MCP" ~ "MCP",
      model == "HDMA_lasso" ~ "HDMA (LASSO)",
      model == "HDMA_SCAD" ~ "HDMA (SCAD)",
      model == "HDMA_MCP" ~ "HDMA (MCP)",
      TRUE ~ model
    ) %>% factor(levels = c("LASSO", "SCAD", "MCP", 
                            "HDMA (LASSO)", "HDMA (SCAD)", "HDMA (MCP)"))
  )

p <- ggplot(df_long, 
            aes(y = estimate, x = variable, 
                color = model, shape = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.4,
    position = position_dodge(width = 0.6),
    linewidth = 0.3
  ) +
  geom_point(
    size = 1.2,
    position = position_dodge(width = 0.6)
  ) +
  scale_color_manual(
    name = "Model",
    values = c(
      "LASSO"         = "#7bccc4", 
      "SCAD"          = "#43a2ca", 
      "MCP"           = "#0868ac", 
      "HDMA (LASSO)"  = "#88419d",  
      "HDMA (SCAD)"   = "#8c2d04",  
      "HDMA (MCP)"    = "#4d004b"   
    )
  )+
  scale_shape_manual(name="Model",values = c(15, 16, 17, 18, 19, 3)) +
  labs(x = NULL, y = "95% Simultaneous Confidence Intervals") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 45,  
      hjust = 1,
      vjust = 0.5,
      size = 8,
      margin = margin(t = -10)
    ),
    legend.position = c(0.99, 0.99),
    legend.justification = c(1, 1),
    legend.spacing.x = unit(0.1, "cm"),
    legend.key.size = unit(0.3, "cm"),
    plot.margin = margin(6, 6, 6, 6, "mm"), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  coord_cartesian(expand = FALSE) 

var_count <- length(unique(df_long$variable))
ggsave("Simultaneous_ci.pdf", plot = p,
       width = var_count * 0.6 + 5, 
       height = 15, 
       units = "cm")
    


vec_lasso_test=c();vec_pflasso_test=c();vec_MCP_test=c();
vec_EN_test=c();vec_alasso_test=c();vec_SCAD_test=c();vec_HDMA_Lasso_test=c();
vec_HDMA_SCAD_test=c();vec_HDMA_MCP_test=c();vec_AnL_test=c();vec_PMA_test=c()
vec_HDMA_OPT_test=c()
top10_weights_lasso<-matrix(0,nrow=0,ncol=10)
top10_weights_SCAD<-matrix(0,nrow=0,ncol=10)
top10_weights_MCP<-matrix(0,nrow=0,ncol=10)

for(round in 1:100){
  n=length(Data_Y)
  set.seed(123+round)
  train_part <- sample(1:n, 50)
  train_X<-Data_X[train_part,]
  train_Y<-Data_Y[train_part]
  test_X<-Data_X[-train_part,]
  test_Y<-Data_Y[-train_part]
  
  # Lasso
  lasso_fit<-cv.glmnet(x = train_X,y = train_Y,family = "gaussian",alpha = 1)
  beta_lasso<-as.vector(coef(lasso_fit,s="lambda.min"))
  vec_lasso_test<-c(vec_lasso_test,l2_loss(cbind(1,test_X),test_Y,beta_lasso))

  # Post-fit Lasso 
  beta_pflasso<-pflasso(X = train_X,Y = train_Y)
  vec_pflasso_test<-c(vec_pflasso_test,l2_loss(cbind(1,test_X),test_Y,beta_pflasso))

  # MCP
  MCP_fit<-cv.ncvreg(X = train_X,y = train_Y,family = "gaussian", penalty="MCP")
  beta_MCP<-as.vector(coef(MCP_fit,s="lambda.min"))
  vec_MCP_test<-c(vec_MCP_test,l2_loss(cbind(1,test_X),test_Y,beta_MCP))

  # SCAD
  SCAD_fit<-cv.ncvreg(X = train_X,y = train_Y,family = "gaussian", penalty="SCAD")
  beta_SCAD<-as.vector(coef(SCAD_fit,s="lambda.min"))
  vec_SCAD_test<-c(vec_SCAD_test,l2_loss(cbind(1,test_X),test_Y,beta_SCAD))

  # EN
  EN_fit<-cv.glmnet(x = train_X,y = train_Y,family = "gaussian",alpha = 0.5)
  beta_EN<-as.vector(coef(EN_fit,s="lambda.min"))
  vec_EN_test<-c(vec_EN_test,l2_loss(cbind(1,test_X),test_Y,beta_EN))

  # HDMA-Lasso
  HDMA_Lasso_fit<-HDMA(train_X,train_Y,Kne=4,d2=10,family="gaussian",nest="mix",penalty="Lasso")
  beta_HDMA_Lasso<-HDMA_Lasso_fit[[1]]
  vec_HDMA_Lasso_test<-c(vec_HDMA_Lasso_test,l2_loss(cbind(1,test_X),test_Y,beta_HDMA_Lasso))
  top10_weights_lasso<-rbind(top10_weights_lasso,HDMA_Lasso_fit$weights[1:10])
  print(mean(vec_HDMA_Lasso_test))
  # HDMA-SCAD
  HDMA_SCAD_fit<-HDMA(train_X,train_Y,Kne=4,d2=10,family="gaussian",nest="mix",penalty="SCAD")
  beta_HDMA_SCAD<-HDMA_SCAD_fit[[1]]
  vec_HDMA_SCAD_test<-c(vec_HDMA_SCAD_test,l2_loss(cbind(1,test_X),test_Y,beta_HDMA_SCAD))
  top10_weights_SCAD<-rbind(top10_weights_SCAD,HDMA_SCAD_fit$weights[1:10])
  print(mean(vec_HDMA_SCAD_test))
  # HDMA-MCP
  HDMA_MCP_fit<-HDMA(train_X,train_Y,Kne=4,d2=10,family="gaussian",nest="mix",penalty="MCP")
  beta_HDMA_MCP<-HDMA_MCP_fit[[1]]
  vec_HDMA_MCP_test<-c(vec_HDMA_MCP_test,l2_loss(cbind(1,test_X),test_Y,beta_HDMA_MCP))
  top10_weights_MCP<-rbind(top10_weights_MCP,HDMA_MCP_fit$weights[1:10])
  print(mean(vec_HDMA_MCP_test))
  # HDMA-OPT
  opt_idx<-which.min(c(HDMA_Lasso_fit$obj_values[length(HDMA_Lasso_fit$obj_values)],
                       HDMA_SCAD_fit$obj_values[length(HDMA_SCAD_fit$obj_values)],
                       HDMA_MCP_fit$obj_values[length(HDMA_MCP_fit$obj_values)]))
  vec_HDMA_OPT_test<-c(vec_HDMA_OPT_test,c(l2_loss(cbind(1,test_X),test_Y,beta_HDMA_Lasso),
                                           l2_loss(cbind(1,test_X),test_Y,beta_HDMA_SCAD),
                                           l2_loss(cbind(1,test_X),test_Y,beta_HDMA_MCP))[opt_idx])
  
  # AnL
  AnL_fit<-AnL(train_X,train_Y,ds=0.05*50*c(1:8),family="gaussian")
  beta_AnL<-AnL_fit[[1]]
  vec_AnL_test<-c(vec_AnL_test,l2_loss(cbind(1,test_X),test_Y,beta_AnL))


  # PMA
  PMA_fit<-PMA(train_X,train_Y,lambda=log(length(train_Y)))
  beta_PMA<-PMA_fit[[1]]
  vec_PMA_test<-c(vec_PMA_test,l2_loss(cbind(1,test_X),test_Y,beta_PMA))

}
test <- data.frame(
  Lasso = vec_lasso_test,
  PLasso = vec_pflasso_test,
  SCAD = vec_SCAD_test,
  MCP = vec_MCP_test,
  ENet = vec_EN_test,
  AnL = vec_AnL_test,
  PMA = vec_PMA_test,
  HDMALasso = vec_HDMA_Lasso_test,
  HDMASCAD = vec_HDMA_SCAD_test,
  HDMAMCP = vec_HDMA_MCP_test
)
colnames(test)[c(8,9,10)]<-c("HDMA(Lasso)","HDMA(SCAD)","HDMA(MCP)")
data_long <- pivot_longer(test, cols = names(test), names_to = "Model", values_to = "value")
data_long$Model <- factor(data_long$Model, levels = names(test))

pp <- ggplot(data_long, aes(x = Model, y = value, fill = Model)) + 
  geom_boxplot(
    width = 0.55,
    alpha = 0.9,         
    outlier.shape = 21,   
    outlier.size = 2.3,   
    outlier.color = "gray30",
    outlier.fill = "gray90",
    color = "gray30",    
    lwd = 0.5            
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,       
    size = 2,       
    color = "grey30",    
    fill = "#8C9B9D",   
    stroke = 0.5          
  ) +
  labs(title = "", x = NULL, y = "Prediction Error") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(face = "bold", size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", margin = margin(r = 15)),
    panel.grid.major.y = element_line(color = "gray93", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black"),
    legend.position = "none"
  ) +
  scale_fill_manual(values = Gen_morandi(n = nlevels(data_long$Model))) +
  coord_cartesian(ylim = c(0, 0.4)) +
  geom_text(
    data = aggregate(value ~ Model, data_long, mean),
    aes(label = sprintf("%.3f", value), y = value + 0.01),
    color = "black",   
    size = 4,
    fontface = "bold",
    family = "sans"
  )
ggsave(pp,filename=paste0("Real_Data_Lin.png"),width = 8,height = 6,dpi = 300)

result<-cbind(data.frame(measure=c("Mean","SD","Median")),rbind(apply(test,2,mean),apply(test,2,sd),apply(test,2,median)))
write.csv(result,file="Real_Data_Lin.csv",row.names=F)

colMeans(top10_weights_lasso)
colMeans(top10_weights_SCAD)
colMeans(top10_weights_MCP)
write.csv(top10_weights_lasso,file="Weights_Lin_lasso.csv",row.names=F)
write.csv(top10_weights_SCAD,file="Weights_Lin_SCAD.csv",row.names=F)
write.csv(top10_weights_MCP,file="Weights_Lin_MCP.csv",row.names=F)

save(test,file="Real_Data_Lin.RData")
# load("Real_Data_Lin.RData")

