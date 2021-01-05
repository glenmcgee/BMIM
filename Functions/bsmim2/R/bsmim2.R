

#-------------------------------------------------------------------------------------------------
#' Bayesian Structural Multipollutant Index Model (BSMIM)
#'
#' This function fits the BSMIM. Only gaussian kernel is implemented as of now. Adapted from "bkmrdlm_multi.cpp" in the "regimes" package by Ander Wilson.
#'
#' @param y Vector of outcomes.
#' @param x A list of M matrices. Each element is an (N by L_m) matrix, with each row being an individual's L_m-vector of exposures in the mth index. 
#' @param z A matrix of covariates and confounders. This can be ommited.
#' @param group_id A vector of ids indicating cluster membership for random intercepts model.
#' @param niter Number of MCMC iterations including burnin.
#' @param nburn The number of iterations to be discarded as burnin. Must be less than niter.
#' @param nthin Thinning, every nthin-th draw from the posterior will be saved.
#' @param prior_sigma Vector of length 2 corresponding to hyperparameters \eqn{(a_{\sigma},b_{\sigma})}, the parameters for the gamma prior on \eqn{\sigma^{-2}}
#' @param prior_lambda Hyperparameter \eqn{b_{\lambda}}, the prior variance on \eqn{log(\lambda^{-1})}, (only used in the horseshoe/non-spike-slab model.
#' @param prior_lambdaB Hyperparameter \eqn{b_{\lambda_B}}, the prior variance on \eqn{log(\lambda^{-1}_B)},.
#' @param prior_lambda_shaperate Hyperparameters (\eqn{a_{\lambda}}, \eqn{b_{\lambda}}), the shape and rate parameters for gamma prior.
#' @param prior_tau Hyperparameter \eqn{\tau_{0}^2}, scale parameter for half Cauchy prior for \eqn{\tau}. (default 1)
#' @param kappa Vector of length M of scale parameters, kappa_m. If only one number is provided, we take all \eqn{\kappa_m}=kappa
#' @param basis.opts List with the entries: type = the type of basis used, either 'face' (default) or "ns" or "bs" for splines or "gam" for presmoothing the exposure with a gam following defaults from mgcv; knots = the number of knots used for method face; pve = the percent of variance explained by the PCs for method face; df = the df for ns method.
#' @param horseshoe Use the horseshoe prior for indexwise selection (default FALSE)
#' @param gaussian Use a Gaussian kernel (TRUE, default) or a polynomial kernel (FALSE)
#' @param polydegree Degree of polynomial when polynomial kernel is used.  Only applies when gaussian=FALSE.
#' @return An object of class 'bsmim'.
#' @author Glen McGee and Ander Wilson (adapted from the "regimes" package by Ander Wilson).
#' @importFrom stats model.matrix sd
#' @importFrom GIGrvg rgig
#' @export


bsmim2 <- function(y,
                   x,
                   z,
                   group_id=NULL,
                   niter,
                   nburn=round(niter/2),
                   nthin=1,
                   prior_sigma=c(0.001,0.001),
                   prior_lambda=1, ## under non-spike/slab framework (sd of lognormal)
                   prior_lambdaB=1,
                   prior_lambda_shaperate=c(1,0.1), ## under spike slab (shape and rate for gamma) ## default--> mean 10, sd 10 (same as Bobb 2018)
                   prior_tau=0.01,
                   kappa=1,
                   basis.opts.list=NULL, #list(type="face", pve=.9)
                   horseshoe=FALSE,
                   spike_slab=FALSE,
                   gauss_prior=FALSE, ## gaussian priors on thetastar ## only if spike_slab=TRUE
                   prior_theta_slab_sd=0.25, ## sd of slab
                   prior_theta_slab_bounds=c(0,100), ## bounds of unif for invunif prior on thetastar^2
                   prior_pi=c(1,1), ## prior for pi in spike & slab
                   stepsize_theta=0.1,
                   gaussian=TRUE,
                   polydegree=2,
                   draw_h=FALSE,
                   constraints=NULL, ## 0,1,2 indicating no, positive, and dirichlet constraints for weight priors
                   prior_slabpos=c(0.4,1.6), ## shape and rate for gamma prior on thetastar (under constraints=1)
                   prior_alphas=NULL,    ## M-List of alpha hyperparameters for dirichlet prior on weights (under constraints=2)
                   prior_slabrho=c(4,2)){ ## shape and rate for gamma prior on rho^{1/2} (under constraints=2)
  
  
  
  #####################
  ## other inputs
  #####################
  niter <- round(niter)
  nburn <- round(nburn)
  nthin <- round(nthin)
  if(niter<1){
    stop("niter must be a positive integer")
  }
  if(nburn<0){
    nburn <- 0;
  }
  if(nburn>=niter){
    stop("nburn must be less than niter")
  }
  if(nthin<1){
    nthin <- 1;
  }
  if((niter%%nthin)!=0){
    stop("niter must be a multiple of nburn")
  }
  if(!gaussian){
    if(missing(polydegree)){
      stop("Missing polydegree. Using polynomial kernel of order up to 2.")
      polydegree=2
    }
  }

  #####################
  ## data checks
  #####################
  # sample size
  N <- length(y)
  if(any(is.na(y))){
    stop("missing values are not allowed in y")
  }
  y <- as.matrix(y)

  # design matrix for Z
  if(!missing(z)){
    if(is.null(z)){
      stop("current version only works with at least one non-intercept covariate in z")
    }else{
      if(any(is.na(z))){
        stop("missing values are not allowed in z")
      }
      Z <- model.matrix(~z)
      Z <- Z[,qr(Z)$pivot[1:qr(Z)$rank]][,-1] ## checking multicollinearity
    }
  }else{
    stop("current version only with with at least one non-intercept covariate in z")
  }
  # check that Z has n obs
  if(nrow(Z)!=N){
    stop("number of observations in y and Z must match")
  }


  # check dimensions of X
  if(is.list(x)){
    xnames <- names(x)
    M <- length(x)
    x <- lapply(x, as.matrix)
    dims <- as.data.frame(do.call(rbind, lapply(x,dim)))
    Lm <- dims[,2]
    if(any(dims[,1]!=N)){
      stop("number of observations in y and x must match")
    }
    ## bkmr-dlm required T to be equal for all exposures
    # Tmax <- dims[1,2]
    # if(any(dims[,2]!=Tmax)){
    #   stop("number of columns of each element of x must match")
    # }
    
    # check for missing in x
    if(any(lapply(x,function(x) sum(is.na(x)))>0)){
      stop("missing values are not allowed in x")
    }
    ## current standardization is for each X_im, not for the whole matrix X_m
    # standardize X_m
    mnx <- sdx <- vector(mode="list",length=M)
    for(m in 1:M){
      mnx[[m]] <- apply(x[[m]],2,mean)
      sdx[[m]] <- apply(x[[m]],2,sd)
      x[[m]] <- t(apply(x[[m]],1,function(x) (x-mnx[[m]])/sdx[[m]])) 
      if(dim(x[[m]])[1]!=N){x[[m]] <- t(x[[m]])}
    }
    ## old version standardized all of X_m together
    # # standardize X_m
    # mnx <- sdx <- list(NA,M)
    # for(m in 1:M){
    #   mnx[m] <- mean(x[[m]])
    #   sdx[m] <- sd(x[[m]])
    #   x[[m]] <- (x[[m]]-mnx[m])/sdx[m]
    # }
  }else{
    stop("x must be a list of matrices")
  }
  
  ## check kappa: if a single number, repeate it
  if(length(kappa)>1){
    kappavec=kappa
  }else{
    kappavec <- rep(kappa,M)
  }

  ## build B matrix for variance contribution of random intercepts
  if(!is.null(group_id)){ 
    if(length(group_id)==N){
      randint <- TRUE
      Amat <- model.matrix(~as.factor(group_id)-1) ## design matrix indicating cluster membership
      Bmat <- tcrossprod(Amat)  ## block diagonal matrix indicating cluster membership
    }
  }else{
    randint <- FALSE
    Amat <- matrix(0,N,N)
    Bmat <- matrix(0,N,N)
  }

  ## check params for spike&slab prior
  if(is.null(prior_theta_slab_sd)){
    prior_theta_slab <- 0.25
  }else if(length(prior_theta_slab_sd)>1){
    prior_theta_slab <- 0.25
  }
  
  if(is.null(prior_pi)){
    prior_pi <- c(1,1)
  }else if(length(prior_pi)!=2){
    prior_pi <- c(1,1)
  }
  
  ## check params for informative priors
  if(is.null(constraints)){
    constraints <- rep(0,M)
  }else if(length(constraints)!=M){
    constraints <- rep(0,M)
  }
  if(is.null(prior_alphas)){
    prior_alphas <- vector(mode = "list", length = M)
  }
  for(mm in 1:M){
    if(length(prior_alphas[[mm]])!=Lm[mm]){
      prior_alphas[[mm]] <- rep(1,Lm[mm])
    }
  }

  
  ########################
  ## basis construction
  ########################

  if(is.null(basis.opts.list)){
    basis.opts.list <- vector(mode="list",length=M)   ## if not specified, use unstructured
  }else if(!is.null(basis.opts.list$type)){
    basis.opts.list <- rep(list(basis.opts.list),M)         ## if only one option specified, use for all indices
  }
  if(length(basis.opts.list)!=M){
    basis.opts.list <- vector(mode="list",length=M)   ## if otherwise the wrong length, use unstructured
  }
  
  ### 
  X <- B <- vector(mode="list",length=M)
  for(m in 1:M){
    basis.opts <- basis.opts.list[[m]]
    
    if(is.null(basis.opts)){
      B[[m]]$psi <- diag(ncol(x[[m]]))                ## if not specified, use unstructured
    }else if(is.null(basis.opts$type)){
      B[[m]]$psi <- diag(ncol(x[[m]]))                ## if not specified, use unstructured
    }else{
      if(toupper(basis.opts$type) %in% c("NONE","PCA","NS","BS","FACE","GAM","MEAN","AVERAGE")){
        B[[m]] <- getbasis(x[[m]],basis.opts)  
      }else{
        stop("basis type not recognized.")
      }
    }
    X[[m]] <- x[[m]]%*%B[[m]]$psi             
  }
  if(!is.null(xnames)) names(B) <- xnames
  
  # ### Just using unstructured B=I !!
  # ### commenting out Ander's code
  # X <- B <- vector(mode="list",length=M)
  # for(m in 1:M){
  #   if(toupper(basis.opts$type) %in% c("NS","BS","FACE","GAM","MEAN","AVERAGE")){
  #     ###B[[m]] <- bdlimbasis(x[[m]],basis.opts)  ## currently not actually using it 
  #     B[[m]]$psi <- diag(ncol(x[[m]])) 
  #   }else{
  #     stop("basis type not recognized.")
  #   }
  # 
  #   X[[m]] <- x[[m]]%*%B[[m]]$psi              ## currently not actually using it 
  # }
  # if(!is.null(xnames)) names(B) <- xnames



  ########################
  ## fit model
  ########################
  
  if(spike_slab==TRUE){
    if(gauss_prior==TRUE){
      if(sum(constraints==0)){
        fit <- bsmim_spikeslab_gaussprior_mcmc2(yz=cbind(y,Z),
                                                Xlist= X,
                                                a_lam=prior_lambda_shaperate[1],
                                                b_lam=prior_lambda_shaperate[2],
                                                b_lambdaB=prior_lambdaB[1],
                                                a_sig=prior_sigma[1],
                                                b_sig=prior_sigma[2],
                                                s_theta=prior_theta_slab_sd,
                                                step_theta=stepsize_theta,
                                                a_pi=prior_pi[1],
                                                b_pi=prior_pi[2],
                                                poly=(1-gaussian),
                                                d=polydegree,
                                                randint=randint,
                                                Bmat=Bmat,
                                                draw_h=draw_h,
                                                n_inner=nthin,
                                                n_outer=round(niter/nthin),
                                                n_burn=round(nburn/nthin))
        
      }else{ ### informative priors (+spikeslab varselection +gauss prior)
        fit <- bsmim_informative_mcmc2(yz=cbind(y,Z),
                                       Xlist= X,
                                       a_lam=prior_lambda_shaperate[1],
                                       b_lam=prior_lambda_shaperate[2],
                                       b_lambdaB=prior_lambdaB[1],
                                       a_sig=prior_sigma[1],
                                       b_sig=prior_sigma[2],
                                       s_theta=prior_theta_slab_sd,
                                       step_theta=stepsize_theta,
                                       a_pi=prior_pi[1],
                                       b_pi=prior_pi[2],
                                       poly=(1-gaussian),
                                       d=polydegree,
                                       randint=randint,
                                       Bmat=Bmat,
                                       draw_h=draw_h,
                                       thetaconstraint=constraints,
                                       a_slabpos=prior_slabpos[1], 
                                       b_slabpos=prior_slabpos[2], 
                                       alphas=prior_alphas,    
                                       a_rho=prior_slabrho[1],     
                                       b_rho=prior_slabrho[2],     
                                       n_inner=nthin,
                                       n_outer=round(niter/nthin),
                                       n_burn=round(nburn/nthin))
      }

    }else{
      fit <- bsmim_spikeslab_mcmc2(yz=cbind(y,Z),
                                   Xlist= X,
                                   a_lam=prior_lambda_shaperate[1],
                                   b_lam=prior_lambda_shaperate[2],
                                   b_lambdaB=prior_lambdaB[1],
                                   a_sig=prior_sigma[1],
                                   b_sig=prior_sigma[2],
                                   a_theta=prior_theta_slab_bounds[1],
                                   b_theta=prior_theta_slab_bounds[2],
                                   step_theta=stepsize_theta,
                                   a_pi=prior_pi[1],
                                   b_pi=prior_pi[2],
                                   poly=(1-gaussian),
                                   d=polydegree,
                                   randint=randint,
                                   Bmat=Bmat,
                                   draw_h=draw_h,
                                   n_inner=nthin,
                                   n_outer=round(niter/nthin),
                                   n_burn=round(nburn/nthin))
    }
  }else{
    fit <- bsmim_mcmc2(yz=cbind(y,Z),
                       Xlist= X,
                       b_lambda=prior_lambda[1],
                       b_lambdaB=prior_lambdaB[1],
                       a_sig=prior_sigma[1],
                       b_sig=prior_sigma[2],
                       tau02=prior_tau,
                       kappa=kappavec,
                       poly=(1-gaussian),
                       d=polydegree,
                       horseshoe=horseshoe,
                       randint=randint,
                       Bmat=Bmat,
                       draw_h=draw_h,
                       n_inner=nthin,
                       n_outer=round(niter/nthin),
                       n_burn=round(nburn/nthin))
  }

  ########################
  ## remove burnin
  ########################
  #remove burnin for theta, w, nu (all lists)
  fit$w <- list()
  for(m in 1:M){
    
    ## theta
    fit$theta[[m]] <- as.matrix(fit$theta[[m]][(nburn/nthin+1):(niter/nthin),])
    ## handle NAs
    fit$theta[[m]][is.na(fit$theta[[m]])] <- 0
    # flip sign to impose positivity constraint
    cmw <- rowMeans(fit$theta[[m]]%*%t(B[[m]]$psi))
    if(any(cmw<0)){
      fit$theta[[m]][which(cmw<0),] <- -fit$theta[[m]][which(cmw<0),]
    }

    # ### CHECK RESCALING:
    # ## 1) do we need this rescaling?
    # ## 2) ander parameterizes rho differently--I don't think it should be divided here... maybe im wrong
    # # rescale to divide by number of exposures
    # fit$theta[[m]] <- fit$theta[[m]]*sqrt(nrow(B[[m]]$psi))
    # fit$rho[,m] <- fit$rho[,m]*sqrt(nrow(B[[m]]$psi)) ## also ander parameterizes rho differently I think so i doubt we would divide here
  
    ## w
    fit$w[[m]] <- fit$theta[[m]]%*%t(B[[m]]$psi)
    
  }
  if(!is.null(xnames)) names(fit$theta) <- xnames
  
  if(sum(constraints)>0){
    fit$wPOS <- list()
    for(m in 1:M){
      if(constraints[m]>0){
        ## theta
        fit$thetaPOS[[m]] <- as.matrix(fit$thetaPOS[[m]][(nburn/nthin+1):(niter/nthin),])
        ## handle NAs
        fit$thetaPOS[[m]][is.na(fit$thetaPOS[[m]])] <- 0
        ##
        fit$wPOS[[m]] <- fit$thetaPOS[[m]]%*%t(B[[m]]$psi)
      }
      
    }
    if(!is.null(xnames)) names(fit$thetaPOS) <- xnames}
    


  #remove burnin for other values
  fit$lambdaInverse <- fit$lambdaInverse[(nburn/nthin+1):(niter/nthin),]
  fit$lambdaBInverse <- fit$lambdaBInverse[(nburn/nthin+1):(niter/nthin),]
  fit$rho <- fit$rho[(nburn/nthin+1):(niter/nthin),]
  fit$hsamp <- fit$hsamp[(nburn/nthin+1):(niter/nthin),]
  fit$sigma2 <- fit$sigma2[(nburn/nthin+1):(niter/nthin),]
  fit$gamma <- fit$gamma[(nburn/nthin+1):(niter/nthin),]
  if(!is.null(colnames(Z))){
    colnames(fit$gamma) <- colnames(Z)
  } 
  if(spike_slab==TRUE){
    fit$tau <- NULL
    fit$nu <- NULL
    
    for(m in 1:length(fit$theta)){
      fit$theta[[m]][is.na(fit$theta[[m]])] <- 1/sqrt(ncol(fit$theta[[m]]))
    }
    
    
    
    
  }else{
    ## remove burnin for tau and nu
    if(horseshoe==1){ ## under horseshoe prior
      fit$tau <- fit$tau[(nburn/nthin+1):(niter/nthin),]
      for(m in 1:M){ ## nu is a list
        fit$nu[[m]] <- as.matrix(fit$nu[[m]])[(nburn/nthin+1):(niter/nthin),]
      }
    }else{ ## otherwise there is no tau
      fit$tau <- NULL
      for(m in 1:M){ ## delete excess columns for nu
        fit$nu[[m]] <- as.matrix(fit$nu[[m]])[(nburn/nthin+1):(niter/nthin),1]
      }
    }
  }


  
  ########################
  ## check model fit
  ########################
  if(draw_h==FALSE){
    DIC <- NULL
    WAIC <- NULL
    LPPD <- NULL
    RMSE <- NULL
  }else{
    
    ## draw random intercepts if necessary 
    ### NOTE: currently using lppd for NEW clusters, i.e. drawing b_(K+1) from N(0,\sigma^2_B), rather than for existing clusters (this is to keep consistent with predictions in cross validation)
    if(randint==TRUE){
      bk <- mvrnorm(length(levels(as.factor(group_id))),mu=rep(0,length(fit$sigma2)),Sigma=diag(fit$lambdaBInverse*fit$sigma2))## then repeat for cluster members
      bki <- Amat%*%bk
    }else{
      bki <- matrix(0,nrow=length(y),ncol=length(fit$sigma2))
    }
        
    ## Compute DIC
    Dthetamean  <- log(2*pi) + log(mean(fit$sigma2)) + c(((y-as.matrix(apply(fit$hsamp,2,mean))-as.matrix(apply(bki,1,mean))-Z%*%apply(fit$gamma,2,mean))^2)/mean(fit$sigma2))                                 ## -2 logP(y|\hat{\theta}) where \hat{\theta} is the posterior mean
    Dbar        <- log(2*pi) + mean(log(fit$sigma2)) + apply((matrix(y,nrow=nrow(as.matrix(y)),ncol=length(fit$sigma2))-t(fit$hsamp)-bki-Z%*%t(fit$gamma))^2,1,function(x) mean(x/fit$sigma2) )  ## -2 (1/S)\sum_{s}^S logP(y|\theta^s) 
    pD  <- Dbar-Dthetamean  ## -2 (1/S)\sum_{s}^S logP(y|\theta^s) -(-2 logP(y|\hat{\theta}))
    DIC <- Dthetamean+2*pD  ## DIC formula (from BDA 3)
    DIC <- sum(DIC)         ## sum over y
    
    ## Compute WAIC2
    lppd <- log(apply( exp( -0.5*log(2*pi) -0.5*log(matrix(fit$sigma2,nrow=length(fit$sigma2),ncol=nrow(as.matrix(y)))) -0.5*t((matrix(y,nrow=nrow(as.matrix(y)),ncol=length(fit$sigma2))-t(fit$hsamp)-bki-Z%*%t(fit$gamma))^2)/fit$sigma2  ),2,mean  )) ## -2 log[(1/S)\sum_{s}^S P(y|\theta^s)] 
    pW <- apply(-0.5*log(2*pi) -0.5*log(matrix(fit$sigma2,nrow=length(fit$sigma2),ncol=nrow(as.matrix(y)))) -0.5*t((matrix(y,nrow=nrow(as.matrix(y)),ncol=length(fit$sigma2))-t(fit$hsamp)-bki-Z%*%t(fit$gamma))^2)/fit$sigma2,2,var)   ## -2 (1/S)\sum_{s}^S logP(y|\theta^s) 
    WAIC <- -2*(lppd-pW)        ## WAIC2 formula (from BDA 3)
    WAIC <- sum(WAIC) # return(list(WAIC=sum(WAIC),SE=sqrt(length(y)*var(WAIC))))
    
    ## Compute LPPD for observed data
    LPPD <- sum(lppd)
    
    ## Compute RMSE (summed over observations)
    rmse <- sqrt(apply((matrix(y,nrow=nrow(as.matrix(y)),ncol=length(fit$sigma2))-t(fit$hsamp)-bki-Z%*%t(fit$gamma))^2,1,mean))
    RMSE <- sum(rmse)
  }
    


  ########################
  ## return values
  ########################
  ## return data
  fit$y <- y
  fit$x <- x
  fit$z <- Z
  fit$xscale <- list(mean=mnx, sd=sdx)
  fit$basis <- B
  fit$randint <- randint
  fit$Amat <- Amat
  fit$Bmat <- Bmat
  
  ## return results
  fit$hsummary = data.frame(mean=fit$hmean,
                            sd = sqrt(diag(fit$hcov)),
                            lower = fit$hmean + qnorm(0.05/2) * sqrt(diag(fit$hcov)),
                            upper = fit$hmean + qnorm(1 - 0.05/2) * sqrt(diag(fit$hcov)))
  fit$DIC <- DIC
  fit$WAIC <- WAIC
  fit$LPPD <- LPPD
  fit$RMSE <- RMSE
  
  ## pass vars from original call
  fit$group_id <- group_id
  fit$niter <- niter
  fit$nburn <- nburn
  fit$nthin <- nthin
  fit$prior_sigma <- prior_sigma
  fit$prior_lambda <- prior_lambda
  fit$prior_lambda_shaperate <- prior_lambda_shaperate
  fit$prior_lambdaB <- prior_lambdaB
  fit$prior_tau <- prior_tau
  fit$kappa <- kappa
  fit$basis.opts.list <- basis.opts.list
  fit$horseshoe <- horseshoe
  fit$gaussian <- gaussian
  fit$polydegree <- polydegree
  fit$draw_h <- draw_h
  
  fit$spike_slab <- spike_slab
  fit$prior_pi <- prior_pi
  fit$gauss_prior <- gauss_prior
  fit$prior_theta_slab_sd <- prior_theta_slab_sd
  fit$prior_theta_slab_bounds <- prior_theta_slab_bounds
  
  fit$call <- match.call()
  class(fit) <- "bsmim"
  return(fit)

}
  









#-------------------------------------------------------------------------------------------------
#' K-Fold Cross Validation for BSMIM 
#'
#'
#' @param object an object of class 'bsmim' 
#' @return An list of class 'bsmim'.
#' @author Glen McGee and Ander Wilson (adapted from the "regimes" package by Ander Wilson).
#' @importFrom stats model.matrix sd
#' @importFrom GIGrvg rgig
#' @export

bsmim_crossval2 <- function(object,
                            kfolds=4,
                            correct=FALSE){
  
  if(correct==TRUE & is.null(object$LPPD)){
    print("Warning: draw_h must be set to TRUE in object in order to compute corrected accuracy. Computing uncorrected CV accuracy.")
    correct <- FALSE
  }
  
  ########################
  ## loop over k folds
  ########################
  ## split data
  id_k <- sample(rep(1:kfolds,ceiling(length(object$y)/kfolds))[1:length(object$y)],replace=F)
  if(object$randint==TRUE){
    ## sample by cluster
    id_kk <- sample(rep(1:kfolds,ceiling(length(levels(as.factor(object$group_id)))/kfolds))[1:length(levels(as.factor(object$group_id)))],replace=F)
    id_k <- id_kk[as.numeric(as.factor(object$group_id))]
  }
  
  lppd_cv <- rmse_cv <- 0  ## CV lppd = sum of lppd for test sets fit on the training sets (analogous for rmse)
  lppd_k <- rmse_k <- 0   ## average lppd for full data set, fit on each of the k training sets (analogous for rmse)
  test_out <- pred_out <- c()  ## outcome and predicted outcome
  for(kk in 1:kfolds){
    
    ## fit model to training set (k-1 folds)
    fit_k <- bsmim2(y=object$y[id_k!=kk],
                    x=lapply(object$x,function(x) as.matrix(x)[id_k!=kk,]),
                    z=object$z[id_k!=kk,],
                    group_id=object$group_id[id_k!=kk],
                    niter=object$niter,
                    nburn=object$nburn,
                    nthin=object$nthin,
                    prior_sigma=object$prior_sigma,
                    prior_lambda=object$prior_lambda,
                    prior_lambdaB=object$prior_lambdaB,
                    prior_tau=object$prior_tau,
                    kappa=object$kappa,
                    basis.opts.list=object$basis.opts.list, 
                    horseshoe=object$horseshoe,
                    gaussian=object$gaussian,
                    polydegree=object$polydegree,
                    draw_h=correct) ## if using correction, need lppd for full data
    
    ## predict on test set (kth fold)
    pred_k <- predict_hnew_X2(fit_k,
                              newX=lapply(object$x,function(x) as.matrix(x)[id_k==kk,]),
                              newY=object$y[id_k==kk],
                              newZ=object$z[id_k==kk,],
                              newAmat=object$Amat[id_k==kk,],
                              rawX=FALSE) ## using preprocessed X data
    
    ## prediction accuracy on test set
    lppd_cv <- lppd_cv + pred_k$LPPD
    rmse_cv <- rmse_cv + pred_k$RMSE
    
    ## save predicted and test outcomes
    pred_out <- c(pred_out,pred_k$pred_out$mean)
    test_out <- c(test_out,object$y[id_k==kk])
    
    ## if using correction, need prediction accuracy for full data based on just kth training set (scale by K since we take the average)
    if(correct==TRUE){                              
      lppd_k <- lppd_k + (fit_k$LPPD+pred_k$LPPD)/kfolds ## average lppd for full data set, fit on each of the k training sets
      rmse_k <- rmse_k + (fit_k$RMSE+pred_k$RMSE)/kfolds ## average rmse for full data set, fit on each of the k training sets
    }
  }
  

  
  ## compute uncorrected accuracy 
  LPPD <- lppd_cv
  RMSE <- rmse_cv
  
  if(correct==TRUE){
    ## if using correction, need prediction accuracy for full data based on full data fit
    lppd_full <- object$LPPD ## lppd from full data fit on full data
    rmse_full <- object$RMSE ## rmse from full data fit on full data
    ## compute corrected accuracy 
    LPPD_corr <- lppd_cv+(lppd_full-lppd_k)
    RMSE_corr <- rmse_cv+(rmse_full-rmse_k)
  }else{
    LPPD_corr <- NULL
    RMSE_corr <- NULL
  }
  
  
  mse_out <- mean((pred_out-test_out)^2) ## for mean outcome, i.e. fitted val
  bias_out <- mean((pred_out-test_out))
  
  
  ########################
  ## return values
  ########################
  return(list(LPPD=LPPD,
              RMSE=RMSE,
              LPPD_corr=LPPD_corr,
              RMSE_corr=RMSE_corr,
              mse_out=mse_out,
              bias_out=bias_out))
  

}



#-------------------------------------------------------------------------------------------------
#' Summarize posterior of thetas (weights)
#'
#' This function takes as input a bsmim object and summarizes the distribution of estimates of the weights, theta.
#'
#' @param y obj An object of class bsmim
#' @return A list of length M of data.frames summarizing theta estimates.
#' @author Glen McGee
#' @export
summarize_thetas <- function(obj){
  res <- list()
  for(jj in 1:length(obj$theta)){
    theta_raw <- obj$theta[[jj]]
    rho0 <- apply(theta_raw,1,function(x) sum(x!=0)==0) ## when all thetas are 0 (i.e. rho==0), theta is not well defined
    prop_rho0 <- mean(rho0)
    theta <- theta_raw[!rho0,] ## consider only well defined values # theta[apply(theta,1,function(x) sum(x!=0)==0),] <- 1/sqrt(ncol(theta))
    if(ncol(obj$theta[[jj]])==1){
      prop_theta0 <- 1
      post_mean <- 1
      theta_est <- 1
      post_lci <- 1
      post_uci <- 1
    }else{
      prop_theta0 <- apply(theta,2,function(x) mean(x==0))
      post_mean <- apply(theta,2,mean)
      theta_est <- post_mean/sqrt(sum(post_mean^2)) ## standardize posterior mean to satisfy constrain
      post_lci <- apply(theta,2,function(x) quantile(x,0.025)) ## credible interval
      post_uci <- apply(theta,2,function(x) quantile(x,0.975)) ## credible interval 
    }
        theta_df <- data.frame(PIP_RHO=1-prop_rho0,PIP=1-prop_theta0,est=theta_est,lci=post_lci,uci=post_uci) ## inference for theta conditional on rho!=0
    res[[jj]] <- theta_df
  }
  return(res)
}


