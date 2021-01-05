####################################
###       NHANES analysis        ###
####################################
## Oct 8, 2020
## Analyze full NHANES data (n=1003)

## set to TRUE to run locally ## FALSE is on cluster
runLOCAL=FALSE 

## params
R <- 150000            ## no. of iterations
burn <- 0.5           ## percent burn-in
thin <- 25            ## thinning number
doLog <- FALSE        ## dont log transform the exposures
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
  mod_version <- 1  ##  model version
  iter_no <- 1 ## CV fold
  suffix <- paste0("_analysis") #paste0("_analysis","_n",n) ## output file doesnt have a suffix
  nhanes <- na.omit(read.csv("../NHANES/Data/studypop.csv"))
} else{
  path <- "Results/" ## path for results
  args <- commandArgs(trailingOnly = TRUE) ## collect arguments from shell script
  iter_no <- as.integer(args[1])  ## get iter_no via slurm (1 through 5)
  mod_version <- as.integer(args[2])  ## get model version via slurm
  suffix <- paste0("_analysis_CV",iter_no) ## append output file names with the iteration number, to be combined later
  nhanes <- na.omit(read.csv("studypop.csv")) ## read in complete data only 
}


##########################
###  Pre-process data  ###  
##########################
source("NHANES_cleandat.R")
resample_ids <- sample(nrow(nhanes))
if(iter_no==1){ ## not divisible by 5 so being careful about indexing
  CV_ids <- resample_ids[-(1:201)]
}else if(iter_no==2){
  CV_ids <- resample_ids[-(202:402)]
}else if(iter_no==3){
  CV_ids <- resample_ids[-(403:603)]
}else if(iter_no==4){
  CV_ids <- resample_ids[-(604:803)]
}else if(iter_no==5){
  CV_ids <- resample_ids[-(804:1003)]
}
dat <- prep_data_split(CV_ids)#
y <- dat$y
y_TEST <- dat$y_TEST



################################################
###         Fit Models                       ###
################################################
if(mod_version==1){ ## qgcomp method
  
  df <- data.frame(cbind(y,dat$SIM$X[[1]],dat$covariates))
  deg <- 1
  exp_names <- colnames(dat$SIM$X[[1]])
  # fit <- qgcomp(y~poly(PCB074,deg,raw=T)+poly(PCB099,deg,raw=T)+poly(PCB138,deg,raw=T)+poly(PCB153,deg,raw=T)+poly(PCB170,deg,raw=T)+poly(PCB180,deg,raw=T)+poly(PCB187,deg,raw=T)+poly(PCB194,deg,raw=T)+V10+V11+male+bmicat2+bmicat3,expnms=exp_names, data = df,q=10,degree=deg)
  fit <- qgcomp(y~PCB074+PCB099+PCB138+PCB153+PCB170+PCB180+PCB187+PCB194+PCB126+PCB169+PCB118+Dioxin1+Dioxin2+Dioxin3+Furan1+Furan2+Furan3+Furan4+V20+V21+male+bmicat2+bmicat3,expnms=exp_names, data = df,q=10)
  df_TEST <- getQuants(fit,dat$SIM$X_TEST[[1]],dat$covariates_TEST)     ## quantize data for prediction 
  pred_TEST <- predict(fit$fit,expnms=exp_names,newdata=df_TEST)
  qgcomp_list <- list(fit=fit,y=y,y_TEST=y_TEST,pred_TEST=pred_TEST)
  save(qgcomp_list, file=paste0(path,"/Lists/qgcomp_list",suffix,".RData") )
  
}else if(mod_version==2){ ## qgcomp +quadratic
  
  df <- data.frame(cbind(y,dat$SIM$X[[1]],dat$covariates))
  deg <- 2
  exp_names <- colnames(dat$SIM$X[[1]])
  # fit <- qgcomp(y~poly(PCB074,deg,raw=T)+poly(PCB099,deg,raw=T)+poly(PCB138,deg,raw=T)+poly(PCB153,deg,raw=T)+poly(PCB170,deg,raw=T)+poly(PCB180,deg,raw=T)+poly(PCB187,deg,raw=T)+poly(PCB194,deg,raw=T)+V10+V11+male+bmicat2+bmicat3,expnms=exp_names, data = df,q=10,degree=deg)
  fit <- qgcomp(y~poly(PCB074,deg,raw=T)+poly(PCB099,deg,raw=T)+poly(PCB138,deg,raw=T)+poly(PCB153,deg,raw=T)+poly(PCB170,deg,raw=T)+poly(PCB180,deg,raw=T)+poly(PCB187,deg,raw=T)+poly(PCB194,deg,raw=T)+poly(PCB126,deg,raw=T)+poly(PCB169,deg,raw=T)+poly(PCB118,deg,raw=T)+poly(Dioxin1,deg,raw=T)+poly(Dioxin2,deg,raw=T)+poly(Dioxin3,deg,raw=T)+poly(Furan1,deg,raw=T)+poly(Furan2,deg,raw=T)+poly(Furan3,deg,raw=T)+poly(Furan4,deg,raw=T)+V20+V21+male+bmicat2+bmicat3,expnms=exp_names, data = df,q=10,degree=deg)
  df_TEST <- getQuants(fit,dat$SIM$X_TEST[[1]],dat$covariates_TEST)     ## quantize data for prediction 
  pred_TEST <- predict(fit$fit,expnms=exp_names,newdata=df_TEST)
  qgcomp2_list <- list(fit=fit,y=y,y_TEST=y_TEST,pred_TEST=pred_TEST)
  save(qgcomp2_list, file=paste0(path,"/Lists/qgcomp2_list",suffix,".RData") )
  
}else if(mod_version==3){ ## fit 1-dim single index model (SIM) + LINEAR kernel

  fit <- bsmim2(y=y,x=dat$SIM$X,z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=FALSE,polydegree=1,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=0.2,basis.opts=NULL,draw_h=FALSE)
  pred_TEST <- predict_hnew_X2(fit,newX=dat$SIM$X_TEST,newY=y_TEST,newZ=dat$covariates_TEST)
  SIMlinear_list <- list(fit=fit,y=y,y_TEST=y_TEST,pred_TEST=pred_TEST)
  save(SIMlinear_list, file=paste0(path,"/Lists/SIMlinear_list",suffix,".RData") )
  
}else if(mod_version==4){ ## fit 1-dim single index model (SIM)
  
  fit <- bsmim2(y=y,x=dat$SIM$X,z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=0.2,basis.opts=NULL,draw_h=FALSE)
  pred_TEST <- predict_hnew_X2(fit,newX=dat$SIM$X_TEST,newY=y_TEST,newZ=dat$covariates_TEST)
  SIM_list <- list(fit=fit,y=y,y_TEST=y_TEST,pred_TEST=pred_TEST)
  save(SIM_list, file=paste0(path,"/Lists/SIM_list",suffix,".RData") )
  
}else if(mod_version==5){ ## fit 1-dim single index model (SIM) + inverse uniform prior
  
  fit <- bsmim2(y=y,x=dat$SIM$X,z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=FALSE,prior_theta_slab_bounds=c(0,100),stepsize_theta=0.2,basis.opts=NULL,draw_h=FALSE)
  pred_TEST <- predict_hnew_X2(fit,newX=dat$SIM$X_TEST,newY=y_TEST,newZ=dat$covariates_TEST)
  SIMinvunif_list <- list(fit=fit,y=y,y_TEST=y_TEST,pred_TEST=pred_TEST)
  save(SIMinvunif_list, file=paste0(path,"/Lists/SIMinvunif_list",suffix,".RData") )
  
}else if(mod_version==6){ ## fit 3-dim multi index model (bsmim)

  fit <- bsmim2(y=y,x=dat$bsmim$X,z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=0.2,basis.opts=NULL,draw_h=FALSE)
  pred_TEST <- predict_hnew_X2(fit,newX=dat$bsmim$X_TEST,newY=y_TEST,newZ=dat$covariates_TEST)
  bsmim_list <- list(fit=fit,y=y,y_TEST=y_TEST,pred_TEST=pred_TEST)
  save(bsmim_list, file=paste0(path,"/Lists/bsmim_list",suffix,".RData") )
  
}else if(mod_version==7){ ## fit 3-dim multi index model (bsmim) with inverse uniform prior

  fit <- bsmim2(y=y,x=dat$bsmim$X,z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=FALSE,prior_theta_slab_bounds=c(0,100),stepsize_theta=0.2,basis.opts=NULL,draw_h=FALSE)
  pred_TEST <- predict_hnew_X2(fit,newX=dat$bsmim$X_TEST,newY=y_TEST,newZ=dat$covariates_TEST)
  bsmiminvunif_list <- list(fit=fit,y=y,y_TEST=y_TEST,pred_TEST=pred_TEST)
  save(bsmiminvunif_list, file=paste0(path,"/Lists/bsmiminvunif_list",suffix,".RData") )
  
}else if(mod_version==8){ ## fit 18-dim multi index model (bkmr)

  fit <- bsmim2(y=y,x=dat$bkmr$X,z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=0.2,basis.opts=NULL,draw_h=FALSE)
  pred_TEST <- predict_hnew_X2(fit,newX=dat$bkmr$X_TEST,newY=y_TEST,newZ=dat$covariates_TEST)
  bkmr_list <- list(fit=fit,y=y,y_TEST=y_TEST,pred_TEST=pred_TEST)
  save(bkmr_list, file=paste0(path,"/Lists/bkmr_list",suffix,".RData") )
  
}else if(mod_version==9){ ## fit 18-dim multi index model (bkmr) + invunif

  fit <- bsmim2(y=y,x=dat$bkmr$X,z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=FALSE,prior_theta_slab_bounds=c(0,100),stepsize_theta=0.2,basis.opts=NULL,draw_h=FALSE)
  pred_TEST <- predict_hnew_X2(fit,newX=dat$bkmr$X_TEST,newY=y_TEST,newZ=dat$covariates_TEST)
  bkmrinvunif_list <- list(fit=fit,y=y,y_TEST=y_TEST,pred_TEST=pred_TEST)
  save(bkmrinvunif_list, file=paste0(path,"/Lists/bkmrinvunif_list",suffix,".RData") )
  
}else if(mod_version==10){ ## fit bkmr via kmbayes + componentwise variable selection
  
  fit <-  kmbayes(y=y, Z=dat$SIM$X[[1]], X=dat$covariates, iter=R, verbose=FALSE, varsel=TRUE)
  pred_TEST <- SamplePred(fit,Znew=dat$SIM$X_TEST[[1]],Xnew=dat$covariates_TEST,sel=sel);pred_TEST <- pred_CI(pred_TEST)
  km_list <- list(fit=fit,y=y,y_TEST=y_TEST,pred_TEST=pred_TEST)
  save(km_list, file=paste0(path,"/Lists/km_list",suffix,".RData") )
  
}else if(mod_version==11){ ## fit bkmr via kmbayes + low variable selection 
  
  fit <-  kmbayes(y=y, Z=dat$SIM$X[[1]], X=dat$covariates, iter=R, verbose=FALSE, varsel=TRUE, groups=rep(1:3,times=c(8,2,8)))
  pred_TEST <- SamplePred(fit,Znew=dat$SIM$X_TEST[[1]],Xnew=dat$covariates_TEST,sel=sel);pred_TEST <- pred_CI(pred_TEST)
  kmhier_list <- list(fit=fit,y=y,y_TEST=y_TEST,pred_TEST=pred_TEST)
  save(kmhier_list, file=paste0(path,"/Lists/kmhier_list",suffix,".RData") )
  
}else if(mod_version==12){ ## fit SIM +linear kernel on quantiled data
  # use qgcomp to get same cut points
  df <- data.frame(cbind(y,dat$SIM$X[[1]],dat$covariates))
  deg <- 1
  exp_names <- colnames(dat$SIM$X[[1]])
  qfit <- qgcomp(y~PCB074+PCB099+PCB138+PCB153+PCB170+PCB180+PCB187+PCB194+PCB126+PCB169+PCB118+Dioxin1+Dioxin2+Dioxin3+Furan1+Furan2+Furan3+Furan4+V20+V21+male+bmicat2+bmicat3,expnms=exp_names, data = df,q=10)
  ## fit SIM to quantile data
  df <- qfit$fit$data
  X_TEST <- list(getQuants(qfit,dat$SIM$X_TEST[[1]],dat$covariates_TEST)[,exp_names])
  fit <- bsmim2(y=y,x=list(df[,exp_names]),z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=FALSE,polydegree=1,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=0.2,basis.opts=NULL,draw_h=FALSE)
  pred_TEST <- predict_hnew_X2(fit,newX=X_TEST,newY=y_TEST,newZ=dat$covariates_TEST)
  SIMlinearquant_list <- list(fit=fit,y=y,y_TEST=y_TEST,pred_TEST=pred_TEST)
  save(SIMlinearquant_list, file=paste0(path,"/Lists/SIMlinearquant_list",suffix,".RData") )
  
}

