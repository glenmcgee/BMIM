################################################
###        Sim A: Non-Linear Single Index    ###
################################################
## Oct 28, 2020
## Run simulation A2 on subsample of NHANES data 

## set to TRUE to run locally ## FALSE is on cluster
runLOCAL=FALSE 

## params
R <- 100000            ## no. of iterations
burn <- 0.4            ## percent burn-in
thin <- 40             ## thinning number
doLog <- FALSE         ## dont log transform the exposures
folds <- 4             ## no. of folds for CV
jump <- 0.35           ## sd of random walk for theta* bsmim
sel <- seq(burn*R+1,R,by=thin) 
set.seed(1000)

### load libraries
library(bsmim2)
library(bkmr)
library(GIGrvg)
library(refund)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(splines)
library(MASS)
library(qgcomp)


### set up cluster 
if(runLOCAL==TRUE){
  path <- "../Final Simulations/Results/" ## path for results
  n <- 400          ##  sample size
  sd <- 0.8         ##  of errors 
  mod_version <- 1  ##  model version
  iter_no <- 0      ##  iteration number
  suffix <- paste0("_simA2","_n",n,"_sd0",sub("\\.","",sd),"_iter") #paste0("_analysis","_n",n) ## output file doesnt have a suffix
  nhanes <- na.omit(read.csv("../NHANES/Data/studypop.csv"))
  # nhanes <- nhanes_raw
  # set.seed(1000+iter_no)
} else{
  path <- "Results/" ## path for results
  args <- commandArgs(trailingOnly = TRUE) ## collect arguments from shell script
  n <- as.integer(args[1])            ## get sample size via slurm
  sd <- as.integer(args[2])/10        ## get of errors via slurm
  mod_version <- as.integer(args[3])  ## get model version via slurm
  iter_no <- as.integer(args[4])      ## get iteration number
  suffix <- paste0("_simA2","_n",n,"_sd0",sub("\\.","",sd),"_iter",iter_no) ## append output file names with the iteration number, to be combined later
  print(paste0("mod",mod_version,suffix))
  nhanes <- na.omit(read.csv("studypop.csv")) ## read in complete data only 
  # set.seed(1000+iter_no) 
  # nhanes <- nhanes[sample(nrow(nhanes)),] ## randomly reorder data
}


##########################
###  Pre-process data  ###  
##########################
set.seed(0) ## use sametraining set
# nhanes <- nhanes[1:500,]
nhanes <- nhanes[sample(nrow(nhanes))[1:500],]
source("NHANES_cleandat.R")
set.seed(1000+iter_no) ## different outcomes/errors for each dataset
resample_ids <- 1:n  ## same exposures each iteration
dat <- prep_data_split(resample_ids)
# resample_ids <- sample(nrow(nhanes))[1:n]
# dat <- prep_data_split(resample_ids)
# y <- dat$y
# y_TEST <- dat$y_TEST


################################################
###         Data-Generating Mechanism        ###
################################################

# make true weight function
w1 <- c(8:1,rep(0,2),rep(-2,8)); w1 <- w1/sqrt(sum(w1^2)) ## different for A1

# exposure reponse
wx1 <- scale(X%*%w1);mn1 <- mean(X%*%w1);sd1 <- sd(X%*%w1);med1 <- c(apply(X,2,median)%*%w1-mn1)/sd1

## exposure response functions
hfun <- function(z) 1.5/(1+exp(-15*z)) -0.25*z #1.5/(1+exp(-15*z))

h <- hfun(wx1)
h_TRAIN <- h[resample_ids]
h_TEST <- h[-resample_ids]

## covariate effects 
gamma <- c(-0.43,0.00,-0.25,0.12,0.08)#,-0.06,0.23,-0.12,-0.02,0.02,0.08,-0.07,0.01,0.02,-0.04,-0.18,-0.28,-0.15)

## generate outcomes for entire dataset
yfull <- h + covariates%*%gamma + rnorm(nrow(nhanes),0,1)*sd
y <- dat$y <- yfull[resample_ids]
y_TEST <- dat$y_TEST <- yfull[-resample_ids]

#############################
###     true values       ###
#############################
## X25, X50, X75 are based on entire sample, not subset
h25 <- hfun((dat$SIM$X25[[1]]%*%w1-mn1)/sd1)
h50 <- hfun((dat$SIM$X50[[1]]%*%w1-mn1)/sd1)
h75 <- hfun((dat$SIM$X75[[1]]%*%w1-mn1)/sd1)

### true overall association
# overall_true <- diff(apply(dat$SIM$X[[1]],2,function(x) quantile(x,c(0.5,0.6)))%*%beta)
overall_true <- diff(hfun((apply(dat$SIM$X[[1]],2,function(x) quantile(x,c(0.5,0.6)))%*%w1-mn1)/sd1))



################################################
###         Fit Models                       ###
################################################
if(mod_version==1){ ## qgcomp method
  
  df <- data.frame(cbind(y,dat$SIM$X[[1]],dat$covariates))
  deg <- 1
  exp_names <- colnames(dat$SIM$X[[1]])
  fit <- qgcomp(y~PCB074+PCB099+PCB138+PCB153+PCB170+PCB180+PCB187+PCB194+PCB126+PCB169+PCB118+Dioxin1+Dioxin2+Dioxin3+Furan1+Furan2+Furan3+Furan4+V20+V21+male+bmicat2+bmicat3,expnms=exp_names, data = df,q=10)
  # pred_overall <- getOverall_qgcomp(fit,degree=deg)
  pred_TEST <- pred_QGC(fit=fit,Xnew=dat$SIM$X_TEST[[1]],Znew=dat$covariates_TEST,fitted=TRUE,exp_names=exp_names)
  pred_TRAIN <- pred_QGC(fit=fit,Xnew=dat$SIM$X[[1]],Znew=dat$covariates,fitted=TRUE,exp_names=exp_names)
  # pred25 <- pred_QGC(fit=fit,Xnew=dat$SIM$X25[[1]],Znew=0*covariates[1:nrow(dat$SIM$X25[[1]]),],fitted=FALSE,exp_names=exp_names)
  pred50 <- pred_QGC(fit=fit,Xnew=dat$SIM$X50[[1]],Znew=0*covariates[1:nrow(dat$SIM$X50[[1]]),],fitted=FALSE,exp_names=exp_names)
  # pred75 <- pred_QGC(fit=fit,Xnew=dat$SIM$X75[[1]],Znew=0*covariates[1:nrow(dat$SIM$X75[[1]]),],fitted=FALSE,exp_names=exp_names)
  # pred_ind <- NULL
  # pred_inter <- NULL
  CV <- gqcomp_crossval(fit,kfolds=folds)
  # qgcomp_list <- list(fit=fit,pred_TEST=pred_TEST,pred_overall=pred_overall,pred25=pred25,pred50=pred50,pred75=pred75,pred_ind=pred_ind,pred_inter=pred_inter,y_TEST=y_TEST,h_TEST=h_TEST,overall_true=overall_true,h25=h25,h50=h50,h75=h75,CV=CV)
  qgcomp_list <- list(fit=fit,pred_TEST=pred_TEST,pred_TRAIN=pred_TRAIN,pred50=pred50,y_TEST=y_TEST,h_TEST=h_TEST,h50=h50,CV=CV)
  save(qgcomp_list, file=paste0(path,"/Lists/qgcomp_list",suffix,".RData") )
  
}else if(mod_version==2){ ## qgcomp method--quadratic
  
  df <- data.frame(cbind(y,dat$SIM$X[[1]],dat$covariates))
  deg <- 2
  exp_names <- colnames(dat$SIM$X[[1]])
  fit <- qgcomp(y~poly(PCB074,deg,raw=T)+poly(PCB099,deg,raw=T)+poly(PCB138,deg,raw=T)+poly(PCB153,deg,raw=T)+poly(PCB170,deg,raw=T)+poly(PCB180,deg,raw=T)+poly(PCB187,deg,raw=T)+poly(PCB194,deg,raw=T)+poly(PCB126,deg,raw=T)+poly(PCB169,deg,raw=T)+poly(PCB118,deg,raw=T)+poly(Dioxin1,deg,raw=T)+poly(Dioxin2,deg,raw=T)+poly(Dioxin3,deg,raw=T)+poly(Furan1,deg,raw=T)+poly(Furan2,deg,raw=T)+poly(Furan3,deg,raw=T)+poly(Furan4,deg,raw=T)+V20+V21+male+bmicat2+bmicat3,expnms=exp_names, data = df,q=10,degree=deg)
  # pred_overall <- getOverall_qgcomp(fit,degree=deg)
  pred_TEST <- pred_QGC(fit=fit,Xnew=dat$SIM$X_TEST[[1]],Znew=dat$covariates_TEST,fitted=TRUE,exp_names=exp_names)
  pred_TRAIN <- pred_QGC(fit=fit,Xnew=dat$SIM$X[[1]],Znew=dat$covariates,fitted=TRUE,exp_names=exp_names)
  # pred25 <- pred_QGC(fit=fit,Xnew=dat$SIM$X25[[1]],Znew=0*covariates[1:nrow(dat$SIM$X25[[1]]),],fitted=FALSE,exp_names=exp_names)
  pred50 <- pred_QGC(fit=fit,Xnew=dat$SIM$X50[[1]],Znew=0*covariates[1:nrow(dat$SIM$X50[[1]]),],fitted=FALSE,exp_names=exp_names)
  # pred75 <- pred_QGC(fit=fit,Xnew=dat$SIM$X75[[1]],Znew=0*covariates[1:nrow(dat$SIM$X75[[1]]),],fitted=FALSE,exp_names=exp_names)
  # pred_ind <- NULL
  # pred_inter <- NULL
  CV <- gqcomp_crossval(fit,kfolds=folds)
  # qgcomp2_list <- list(fit=fit,pred_TEST=pred_TEST,pred_overall=pred_overall,pred25=pred25,pred50=pred50,pred75=pred75,pred_ind=pred_ind,pred_inter=pred_inter,y_TEST=y_TEST,h_TEST=h_TEST,overall_true=overall_true,h25=h25,h50=h50,h75=h75,CV=CV)
  qgcomp2_list <- list(fit=fit,pred_TEST=pred_TEST,pred_TRAIN=pred_TRAIN,pred50=pred50,y_TEST=y_TEST,h_TEST=h_TEST,h50=h50,CV=CV)
  save(qgcomp2_list, file=paste0(path,"/Lists/qgcomp2_list",suffix,".RData") )
  
}else if(mod_version==3){ ## fit 1-dim single index model (SIM) +LINEAR KERNEL
  
  fit <- bsmim2(y=y,x=dat$SIM$X,z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=FALSE,polydegree=1,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,basis.opts=NULL,draw_h=FALSE)
  # pred_overall <- predict_hnew_assoc2(fit,overall = TRUE)
  pred_TEST <- predict_hnew_X2(fit,newX=dat$SIM$X_TEST,newY=y_TEST,newZ=dat$covariates_TEST)
  pred_TRAIN <- predict_hnew_X2(fit,newX=dat$SIM$X,newY=y,newZ=dat$covariates)
  # pred25<- predict_hnew_X2(fit,newX=dat$SIM$X25)
  pred50<- predict_hnew_X2(fit,newX=dat$SIM$X50)
  # pred75<- predict_hnew_X2(fit,newX=dat$SIM$X75)
  # pred_ind <- predict_hnew_indexwise2(fit)
  # pred_inter <- NULL
  CV<- bsmim_crossval2(fit,kfolds=folds)
  # SIMlinear_list <- list(fit=fit,pred_TEST=pred_TEST,pred_overall=pred_overall,pred25=pred25,pred50=pred50,pred75=pred75,pred_ind=pred_ind,pred_inter=pred_inter,y_TEST=y_TEST,h_TEST=h_TEST,overall_true=overall_true,h25=h25,h50=h50,h75=h75,CV=CV)
  SIMlinear_list <- list(fit=fit,pred_TEST=pred_TEST,pred_TRAIN=pred_TRAIN,pred50=pred50,y_TEST=y_TEST,h_TEST=h_TEST,h50=h50,CV=CV)
  save(SIMlinear_list, file=paste0(path,"/Lists/SIMlinear_list",suffix,".RData") )
  
}else if(mod_version==4){ ## fit 1-dim single index model (SIM) 
  
  fit <- bsmim2(y=y,x=dat$SIM$X,z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,basis.opts=NULL,draw_h=FALSE)
  # pred_overall <- predict_hnew_assoc2(fit,overall = TRUE)
  pred_TEST <- predict_hnew_X2(fit,newX=dat$SIM$X_TEST,newY=y_TEST,newZ=dat$covariates_TEST)
  pred_TRAIN <- predict_hnew_X2(fit,newX=dat$SIM$X,newY=y,newZ=dat$covariates)
  # pred25<- predict_hnew_X2(fit,newX=dat$SIM$X25)
  pred50<- predict_hnew_X2(fit,newX=dat$SIM$X50)
  # pred75<- predict_hnew_X2(fit,newX=dat$SIM$X75)
  # pred_ind <- predict_hnew_indexwise2(fit)
  # pred_inter <- NULL
  CV<- bsmim_crossval2(fit,kfolds=folds)
  # SIM_list <- list(fit=fit,pred_TEST=pred_TEST,pred_overall=pred_overall,pred25=pred25,pred50=pred50,pred75=pred75,pred_ind=pred_ind,pred_inter=pred_inter,y_TEST=y_TEST,h_TEST=h_TEST,overall_true=overall_true,h25=h25,h50=h50,h75=h75,CV=CV)
  SIM_list <- list(fit=fit,pred_TEST=pred_TEST,pred_TRAIN=pred_TRAIN,pred50=pred50,y_TEST=y_TEST,h_TEST=h_TEST,h50=h50,CV=CV)
  save(SIM_list, file=paste0(path,"/Lists/SIM_list",suffix,".RData") )
  
}else if(mod_version==5){ ## fit 3-dim MIM
  
  fit <- bsmim2(y=y,x=dat$bsmim$X,z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,basis.opts=NULL,draw_h=FALSE)
  # pred_overall <- predict_hnew_assoc2(fit,overall = TRUE)
  pred_TEST <- predict_hnew_X2(fit,newX=dat$bsmim$X_TEST,newY=y_TEST,newZ=dat$covariates_TEST)
  pred_TRAIN <- predict_hnew_X2(fit,newX=dat$bsmim$X,newY=y,newZ=dat$covariates)
  # pred25<- predict_hnew_X2(fit,newX=dat$bsmim$X25)
  pred50<- predict_hnew_X2(fit,newX=dat$bsmim$X50)
  # pred75<- predict_hnew_X2(fit,newX=dat$bsmim$X75)
  # pred_ind <- predict_hnew_indexwise2(fit)
  # pred_inter <- pred_twoway(fit)
  CV<- bsmim_crossval2(fit,kfolds=folds)
  # bsmim_list <- list(fit=fit,pred_TEST=pred_TEST,pred_overall=pred_overall,pred25=pred25,pred50=pred50,pred75=pred75,pred_ind=pred_ind,pred_inter=pred_inter,y_TEST=y_TEST,h_TEST=h_TEST,overall_true=overall_true,h25=h25,h50=h50,h75=h75,CV=CV)
  bsmim_list <- list(fit=fit,pred_TEST=pred_TEST,pred_TRAIN=pred_TRAIN,pred50=pred50,y_TEST=y_TEST,h_TEST=h_TEST,h50=h50,CV=CV)
  save(bsmim_list, file=paste0(path,"/Lists/bsmim_list",suffix,".RData") )
  
}else if(mod_version==6){ ## 18-dim MIM 
  
  fit <- bsmim2(y=y,x=dat$bkmr$X,z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,basis.opts=NULL,draw_h=FALSE)
  # pred_overall <- predict_hnew_assoc2(fit,overall = TRUE)
  pred_TEST <- predict_hnew_X2(fit,newX=dat$bkmr$X_TEST,newY=y_TEST,newZ=dat$covariates_TEST)
  pred_TRAIN <- predict_hnew_X2(fit,newX=dat$bkmr$X,newY=y,newZ=dat$covariates)
  # pred25<- predict_hnew_X2(fit,newX=dat$bkmr$X25)
  pred50<- predict_hnew_X2(fit,newX=dat$bkmr$X50)
  # pred75<- predict_hnew_X2(fit,newX=dat$bkmr$X75)
  # pred_ind <- predict_hnew_indexwise2(fit)
  # pred_inter <- NULL
  CV<- bsmim_crossval2(fit,kfolds=folds)
  # bkmr_list <- list(fit=fit,pred_TEST=pred_TEST,pred_overall=pred_overall,pred25=pred25,pred50=pred50,pred75=pred75,pred_ind=pred_ind,pred_inter=pred_inter,y_TEST=y_TEST,h_TEST=h_TEST,overall_true=overall_true,h25=h25,h50=h50,h75=h75,CV=CV)
  bkmr_list <- list(fit=fit,pred_TEST=pred_TEST,pred_TRAIN=pred_TRAIN,pred50=pred50,y_TEST=y_TEST,h_TEST=h_TEST,h50=h50,CV=CV)
  save(bkmr_list, file=paste0(path,"/Lists/bkmr_list",suffix,".RData") )
  
}else if(mod_version==7){ ## fit bkmr via kmbayes +componentwise variable selection
  
  fit <-  kmbayes(y=y, Z=dat$SIM$X[[1]], X=dat$covariates, iter=R, verbose=FALSE, varsel=TRUE)
  # pred_overall <- OverallRiskSummaries(fit = fit, y = y, Z = dat$SIM$X[[1]], X = dat$covariates, qs = seq(0.05, 0.95, by = 0.05), q.fixed = 0.5, method = "exact",sel=sel); pred_overall <- ovrl_CI(pred_overall)
  pred_TEST <- get_predkm(fit,Znew=dat$SIM$X_TEST[[1]],Xnew=dat$covariates_TEST,sel=sel,fitted=TRUE,samples=TRUE)
  pred_TRAIN <- get_predkm(fit,Znew=dat$SIM$X[[1]],Xnew=dat$covariates,sel=sel,fitted=TRUE,samples=TRUE)
  # pred25 <- get_predkm(fit,dat$SIM$X25[[1]],sel=sel)
  pred50 <- get_predkm(fit,dat$SIM$X50[[1]],sel=sel)
  # pred75 <- get_predkm(fit,dat$SIM$X75[[1]],sel=sel)
  # pred_ind <- NULL
  # pred_inter <- predict_km_inter(fit)
  CV <- kmbayes_crossval(fit,sel=sel,kfolds=folds)
  # km_list <- list(fit=fit,pred_TEST=pred_TEST,pred_overall=pred_overall,pred25=pred25,pred50=pred50,pred75=pred75,pred_ind=pred_ind,pred_inter=pred_inter,y_TEST=y_TEST,h_TEST=h_TEST,overall_true=overall_true,h25=h25,h50=h50,h75=h75,CV=CV)
  km_list <- list(fit=fit,pred_TEST=pred_TEST,pred_TRAIN=pred_TRAIN,pred50=pred50,y_TEST=y_TEST,h_TEST=h_TEST,h50=h50,CV=CV)
  save(km_list, file=paste0(path,"/Lists/km_list",suffix,".RData") )
  
}else if(mod_version==8){ ## fit bkmr via kmbayes +hierarchical variable selection
  
  fit <-  kmbayes(y=y, Z=dat$SIM$X[[1]], X=dat$covariates, iter=R, verbose=FALSE, varsel=TRUE, groups=rep(1:3,times=c(8,2,8)))
  # pred_overall <- OverallRiskSummaries(fit = fit, y = y, Z = dat$SIM$X[[1]], X = dat$covariates, qs = seq(0.05, 0.95, by = 0.05), q.fixed = 0.5, method = "exact",sel=sel); pred_overall <- ovrl_CI(pred_overall)
  pred_TEST <- get_predkm(fit,Znew=dat$SIM$X_TEST[[1]],Xnew=dat$covariates_TEST,sel=sel,fitted=TRUE,samples=TRUE)
  pred_TRAIN <- get_predkm(fit,Znew=dat$SIM$X[[1]],Xnew=dat$covariates,sel=sel,fitted=TRUE,samples=TRUE)
  # pred25 <- get_predkm(fit,dat$SIM$X25[[1]],sel=sel)
  pred50 <- get_predkm(fit,dat$SIM$X50[[1]],sel=sel)
  # pred75 <- get_predkm(fit,dat$SIM$X75[[1]],sel=sel)
  # pred_ind <- NULL
  # pred_inter <- predict_km_inter(fit)
  CV <- kmbayes_crossval(fit,sel=sel,groups_vec=rep(1:3,times=c(8,2,8)),kfolds=folds)
  # kmhier_list <- list(fit=fit,pred_TEST=pred_TEST,pred_overall=pred_overall,pred25=pred25,pred50=pred50,pred75=pred75,pred_ind=pred_ind,pred_inter=pred_inter,y_TEST=y_TEST,h_TEST=h_TEST,overall_true=overall_true,h25=h25,h50=h50,h75=h75,CV=CV)
  kmhier_list <- list(fit=fit,pred_TEST=pred_TEST,pred_TRAIN=pred_TRAIN,pred50=pred50,y_TEST=y_TEST,h_TEST=h_TEST,h50=h50,CV=CV)
  save(kmhier_list, file=paste0(path,"/Lists/kmhier_list",suffix,".RData") )
  
}else if(mod_version==9){ ## fit bkmr via kmbayes +hierarchical variable selection-->one group only
  
  fit <-  kmbayes(y=y, Z=dat$SIM$X[[1]], X=dat$covariates, iter=R, verbose=FALSE, varsel=TRUE, groups=rep(1,times=c(18)))
  # pred_overall <- OverallRiskSummaries(fit = fit, y = y, Z = dat$SIM$X[[1]], X = dat$covariates, qs = seq(0.05, 0.95, by = 0.05), q.fixed = 0.5, method = "exact",sel=sel); pred_overall <- ovrl_CI(pred_overall)
  pred_TEST <- get_predkm(fit,Znew=dat$SIM$X_TEST[[1]],Xnew=dat$covariates_TEST,sel=sel,fitted=TRUE,samples=TRUE)
  pred_TRAIN <- get_predkm(fit,Znew=dat$SIM$X[[1]],Xnew=dat$covariates,sel=sel,fitted=TRUE,samples=TRUE)
  # pred25 <- get_predkm(fit,dat$SIM$X25[[1]],sel=sel)
  pred50 <- get_predkm(fit,dat$SIM$X50[[1]],sel=sel)
  # pred75 <- get_predkm(fit,dat$SIM$X75[[1]],sel=sel)
  # pred_ind <- NULL
  # pred_inter <- predict_km_inter(fit)
  CV <- kmbayes_crossval(fit,sel=sel,groups_vec=rep(1,times=c(18)),kfolds=folds)
  # kmhier_list <- list(fit=fit,pred_TEST=pred_TEST,pred_overall=pred_overall,pred25=pred25,pred50=pred50,pred75=pred75,pred_ind=pred_ind,pred_inter=pred_inter,y_TEST=y_TEST,h_TEST=h_TEST,overall_true=overall_true,h25=h25,h50=h50,h75=h75,CV=CV)
  kmhier1_list <- list(fit=fit,pred_TEST=pred_TEST,pred_TRAIN=pred_TRAIN,pred50=pred50,y_TEST=y_TEST,h_TEST=h_TEST,h50=h50,CV=CV)
  save(kmhier1_list, file=paste0(path,"/Lists/kmhier1_list",suffix,".RData") )
  
}







