
## For the index version: we fix thetabar=mean(theta across iterations)
## then take E_i=x_i thetabar, and take quantiles of E_i across N
## do this forthe exposure of interest, and set crossM to the quantile for the whole index
## all other exposures are set to their individual median values (not index-medians)

#' Predict hnew for BSMIM (old version)
#' 
#'
#' @param object An object of class bsmim
#' @param points Number of points to predict at
#' @param qtl M-list of Lm-vectors containing quantiles to set exposure components to
#' @param newvals optional M-list of (pts by Lm)-matrices containing new values to predict at---all others get set to 0; overrides qtl and points
#'
#' @return a dataframe containing predicted values.
#' @export
#'
predict_hnew_old2 <- function(object,
                         points=20, 
                         qtls=list(NA),
                         newvals=NULL,
                         approx=FALSE){
  
  
  ########################
  ## data checks
  ########################
  ## dimensions
  N <- length(object$y)               ## number of observations, N      
  M <- length(object$x)               ## number of indices, M
  Lm <- unlist(lapply(object$x,ncol)) ## M-vector of index sizes L_m
  
  ## if no quantiles provided, or if wrong dimension, set all quantiles to medians
  if(length(qtls)!=M){
    qtls <- vector(mode = "list",length=M)
  }
  for(m in 1:M){
    if(length(qtls[[m]])!=Lm[m]){
      qtls[[m]] <- rep(0.5,Lm[m])
    }
  }
  ## if newvals is provided but given wrong dimension, set to default
  if(!is.null(newvals)){
    if(length(newvals)!=M){
      newvals <- NULL
    }else{
      for(m in 1:M){
        if(ncol(as.matrix(newvals[[m]]))!=Lm[m]){
          newvals <- NULL
        }
      }
    }
  }
  
  ######################################################
  ## Construct new grid points, quantiles and weights
  ######################################################
  ## compute weights (basically A theta^*, so that X*theta*=x(Atheta*))
  ## compute Xqlist for componentwise quantiles of every exposure
  ## compute grid of points for each index
  weights <- list()                         ## list of (S_iter by Lm) matrices with rows rho^{1/2}w (rows correspond to different iterations)
  gridpoints <- list()                      ## list of (points by Lm ) matgrices whose columns are a grid of points for hnew, for each exposure component
  Xqlist <- list()                          ## list of L_m-vectors for quantile exposure level across N observations
  psi <- list()                             ## prep psi for Rcpp function
  for(m in 1:M){   ## rescale to divide by number of data times and rho
    weights[[m]] <- sqrt(as.matrix(object$rho)[,m])*object$theta[[m]] %*% t(object$basis[[m]]$psi) # / sqrt(nrow(object$basis[[m]]$psi))         ## removing the scaling we did in for theta
    psi[[m]] <- object$basis[[m]]$psi  ## list
    
    if(is.null(newvals)){
      temp_Xq <- c() ## since empty component of list causes problems #Xqlist[[m]] <- c()
      for(l in 1:Lm[m]){
        temp_Xq <- c(temp_Xq,quantile(object$x[[m]][,l],  qtls[[m]][l]))#Xqlist[[m]] <- c(Xqlist[[m]],quantile(object$x[[m]][,Lm[m]],  qtls[[m]][Lm[m]]))                           ## get appropriate quantile of comparison for each component individually
      }
      Xqlist[[m]] <- temp_Xq
      gridpoints[[m]] <- apply(object$x[[m]],2,function(x) seq(quantile(x,0.05),quantile(x,0.95), length=points))  ## set Lmth column of mth matrix to a grid of points form the 5th to 95th percentile, with evenly spaced points (not at percentiles)
      
    }else{
      Xqlist[[m]] <- rep(0,Lm[m])
      gridpoints[[m]] <- as.matrix(newvals[[m]])  ## set Lmth column of mth matrix to a grid of new points
    }
    
  }
  
  ########################
  ## Run Rcpp Function
  ########################
  if(approx==TRUE){
    predict_hnew <- bsmim_predict_approx_old_cpp2(yz=cbind(object$y,object$z),
                                              Xlist=object$x,
                                              thetalist=object$theta,
                                              psilist=psi,
                                              rho=as.matrix(object$rho),
                                              gamma=as.matrix(object$gamma),
                                              lambdaInverse=object$lambdaInverse,
                                              lambdaBInverse=object$lambdaBInverse,
                                              sigma2=object$sigma2,
                                              weightslist=weights,
                                              gridpointslist=gridpoints,
                                              Xqlist=Xqlist,
                                              poly=(1-object$gaussian),
                                              d=object$polydegree,
                                              randint=object$randint,
                                              Bmat=object$Bmat)
  }else{
    predict_hnew <- bsmim_predict_old_cpp2(yz=cbind(object$y,object$z),
                                       Xlist=object$x,
                                       thetalist=object$theta,
                                       psilist=psi,
                                       rho=as.matrix(object$rho),
                                       gamma=as.matrix(object$gamma),
                                       lambdaInverse=object$lambdaInverse,
                                       lambdaBInverse=object$lambdaBInverse,
                                       sigma2=object$sigma2,
                                       weightslist=weights,
                                       gridpointslist=gridpoints,
                                       Xqlist=Xqlist,
                                       poly=(1-object$gaussian),
                                       d=object$polydegree,
                                       randint=object$randint,
                                       Bmat=object$Bmat)
  }
  
  
  
  
  ########################
  ## return values
  ########################
  fits <- NULL
  cc <- 1 ## column index
  for(m in 1:M){
    for(l in 1:Lm[m]){
      temp <-
        data.frame(m = m,
                   l = l,
                   # name = names(object$x)[m],
                   mean = predict_hnew$hmean[,cc],
                   sd = predict_hnew$hsd[,cc],
                   grid = gridpoints[[m]][,l]*object$xscale$sd[[m]][l]+object$xscale$mean[[m]][l] )  ## grid of points, NOT SCALED OR TRANSFORMED APPROPRIATELY
      fits <- rbind(fits,temp)
      cc <- cc+1
    }
  }
  fits$lower <- fits$mean - 1.96*fits$sd
  fits$upper <- fits$mean + 1.96*fits$sd
  
  class(fits) <- "hpred"
  return(fits)
}

#' Predict hnew for BSMIM
#' 
#'
#' @param object An object of class bsmim
#' @param points Number of points to predict at
#' @param qtl M-list of Lm-vectors containing quantiles to set exposure components to
#' @param newX optional M-list of (pts by Lm)-matrices containing new values to predict at---all others get set to 0; overrides qtl and points
#'
#' @return a dataframe containing predicted values.
#' @export
#'
predict_hnew2 <- function(object,
                          points=20, 
                          qtl_lims=c(0.05,0.95),
                          qtls=list(NA),
                          newX=NULL,
                          approx=FALSE){
  
  
  ########################
  ## data checks
  ########################
  ## dimensions
  N <- length(object$y)               ## number of observations, N      
  M <- length(object$x)               ## number of indices, M
  Lm <- unlist(lapply(object$x,ncol)) ## M-vector of index sizes L_m
  
  ## if no quantiles provided, or if wrong dimension, set all quantiles to medians
  if(length(qtls)!=M){
    qtls <- vector(mode = "list",length=M)
  }
  for(m in 1:M){
    if(length(qtls[[m]])!=Lm[m]){
      qtls[[m]] <- rep(0.5,Lm[m])
    }
  }

  
  ### check newX / set to quantiles
  if(is.null(newX)){ ### if no newX given, construct same grid as the other versions above
    
    ## Construct new grid points, quantiles and weights
    newX <- list()                      ## list of ((points x sum(Lm)) by Lm ) matrices whose columns are a grid of points for hnew
    id <- 1 ## iterator to set grid appropriately
    for(m in 1:M){   
      temp_Xq <- matrix(0,nrow=(points*sum(Lm)),ncol=Lm[m]) ## since empty component of list causes problems #Xqlist[[m]] <- c()
      for(l in 1:Lm[m]){
        temp_Xq[,l] <- quantile(object$x[[m]][,l],qtls[[m]][l])  ## set elements to median
        temp_Xq[points*(id-1)+(1:points),l] <- seq(quantile(object$x[[m]][,l],qtl_lims[1]),quantile(object$x[[m]][,l],qtl_lims[2]), length=points) ## set grid
        id <- id+1 ## increase iterator
      }
      newX[[m]] <- temp_Xq
    }
    
  }
  if(length(newX)!=M){
    stop("newX must be length M")
  }
  points <- nrow(as.matrix(newX[[1]]))
  for(m in 1:M){
    if((nrow(as.matrix(newX[[m]]))!=points) | (ncol(as.matrix(newX[[m]]))!=Lm[m])){
      stop("newX has wrong dimension")
    }else{
      newX[[m]] <- as.matrix(newX[[m]])
    }
  }
  
  
  ######################################################
  ## Construct weights (basically A theta^*, so that X*theta*=x(Atheta*))
  ######################################################
  weights <- list() ## list of (S_iter by Lm) matrices with rows rho^{1/2}w (rows correspond to different iterations)
  psi <- list()     ## prep psi for Rcpp function
  for(m in 1:M){    ## rescale to divide by number of data times and rho
    weights[[m]] <- sqrt(as.matrix(object$rho)[,m])*object$theta[[m]] %*% t(object$basis[[m]]$psi) # / sqrt(nrow(object$basis[[m]]$psi))         ## removing the scaling we did in for theta
    psi[[m]] <- object$basis[[m]]$psi  ## list
    
  }
  
  ########################
  ## Run Rcpp Function
  ########################
  if(approx==TRUE){
    predict_hnew <- bsmim_predict_approx_cpp2(yz=cbind(object$y,object$z),
                                              Xlist=object$x,
                                              thetalist=object$theta,
                                              psilist=psi,
                                              rho=as.matrix(object$rho),
                                              gamma=as.matrix(object$gamma),
                                              lambdaInverse=object$lambdaInverse,
                                              lambdaBInverse=object$lambdaBInverse,
                                              sigma2=object$sigma2,
                                              weightslist=weights,
                                              gridpointslist=newX,
                                              poly=(1-object$gaussian),
                                              d=object$polydegree,
                                              randint=object$randint,
                                              Bmat=object$Bmat)
  }else{
    predict_hnew <- bsmim_predict_cpp2(yz=cbind(object$y,object$z),
                                       Xlist=object$x,
                                       thetalist=object$theta,
                                       psilist=psi,
                                       rho=as.matrix(object$rho),
                                       gamma=as.matrix(object$gamma),
                                       lambdaInverse=object$lambdaInverse,
                                       lambdaBInverse=object$lambdaBInverse,
                                       sigma2=object$sigma2,
                                       weightslist=weights,
                                       gridpointslist=newX,
                                       poly=(1-object$gaussian),
                                       d=object$polydegree,
                                       randint=object$randint,
                                       Bmat=object$Bmat)
  }
  
  
  
  
  ########################
  ## return values
  ########################
  # fits <- NULL
  # cc <- 1 ## column index
  # for(m in 1:M){
  #   for(l in 1:Lm[m]){
  #     temp <-
  #       data.frame(m = m,
  #                  l = l,
  #                  # name = names(object$x)[m],
  #                  mean = predict_hnew$hmean[,cc],
  #                  sd = predict_hnew$hsd[,cc],
  #                  grid = gridpoints[[m]][,l]*object$xscale$sd[[m]][l]+object$xscale$mean[[m]][l] )  ## grid of points, NOT SCALED OR TRANSFORMED APPROPRIATELY
  #     fits <- rbind(fits,temp)
  #     cc <- cc+1
  #   }
  # }
  # fits$lower <- fits$mean - 1.96*fits$sd
  # fits$upper <- fits$mean + 1.96*fits$sd
  # class(fits) <- "hpred"
  # return(fits)
  fits <- data.frame(mean=predict_hnew$hmean,
                     sd=sqrt(diag(predict_hnew$hcov)))
  fits$lower <- fits$mean - 1.96*fits$sd
  fits$upper <- fits$mean + 1.96*fits$sd
  
  results <- list(fits=fits,
                  grid=newX) ## rmse for test data
  class(results) <- "hpred"
  return(results)
  

}



#' Compute exposure-associations: hnew(x*)-hnew(median(x)), setting all other exposures to median
#' 
#'
#' @param object An object of class bsmim
# #' @param points Number of points to predict at
#' @param qtl M-list of Lm-vectors containing quantiles to set exposure components to
# #' @param newX optional M-list of (pts by Lm)-matrices containing new values to predict at---all others get set to 0; overrides qtl and points
#'
#' @return a dataframe containing predicted values.
#' @export
#'
predict_hnew_assoc2 <- function(object,
                          # points=20, 
                          qtl_lims=c(0.05,0.95),
                          qtls=list(NA),
                          # newX=NULL,
                          compare_qtl=0.5,
                          overall=FALSE,
                          approx=FALSE){
  
  
  ########################
  ## data checks
  ########################
  ## dimensions
  N <- length(object$y)               ## number of observations, N      
  M <- length(object$x)               ## number of indices, M
  Lm <- unlist(lapply(object$x,ncol)) ## M-vector of index sizes L_m
  
  ## if no quantiles provided, or if wrong dimension, set all quantiles to medians
  if(length(qtls)!=M){
    qtls <- vector(mode = "list",length=M)
  }
  for(m in 1:M){
    if(length(qtls[[m]])!=Lm[m]){
      qtls[[m]] <- rep(0.5,Lm[m])
    }
  }
  
  if(length(compare_qtl)!=1){
    compare_qtl <- 0.5
  }

  
  ### check newX / set to quantiles
  # if(is.null(newX)){ ### if no newX given, construct same grid as the other versions above
    
    ## Construct new grid points, quantiles and weights
    points=19 ## (0.05,0.10,...,0.95)
    newX <- list()                      ## list of ((points x sum(Lm)) by Lm ) matrices whose columns are a grid of points for hnew
    if(overall==FALSE){
      id <- 1 ## iterator to set grid appropriately
      for(m in 1:M){   
        temp_Xq <- matrix(0,nrow=(points*sum(Lm)),ncol=Lm[m]) ## since empty component of list causes problems #Xqlist[[m]] <- c()
        for(l in 1:Lm[m]){
          temp_Xq[,l] <- quantile(object$x[[m]][,l],qtls[[m]][l])  ## set elements to median
          temp_Xq[points*(id-1)+(1:points),l] <- quantile(object$x[[m]][,l],seq(qtl_lims[1],qtl_lims[2],length=points))  ## set grid
          id <- id+1 ## increase iterator
        }
        temp_Xq <- rbind(temp_Xq,apply(object$x[[m]],2,function(x) quantile(x,compare_qtl))) ## add comparison
        newX[[m]] <- temp_Xq
      }
      
    }else{
      for(m in 1:M){   
        temp_Xq <- matrix(0,nrow=points,ncol=Lm[m]) ## since empty component of list causes problems #Xqlist[[m]] <- c()
        for(l in 1:Lm[m]){
          temp_Xq[,l] <- quantile(object$x[[m]][,l],seq(qtl_lims[1],qtl_lims[2],length=points))  ## set grid
        }
        temp_Xq <- rbind(temp_Xq,apply(object$x[[m]],2,function(x) quantile(x,compare_qtl))) ## add comparison
        newX[[m]] <- temp_Xq
      }
      
    }
    
  # }

  
  
  ######################################################
  ## Construct weights (basically A theta^*, so that X*theta*=x(Atheta*))
  ######################################################
  weights <- list() ## list of (S_iter by Lm) matrices with rows rho^{1/2}w (rows correspond to different iterations)
  psi <- list()     ## prep psi for Rcpp function
  for(m in 1:M){    ## rescale to divide by number of data times and rho
    weights[[m]] <- sqrt(as.matrix(object$rho)[,m])*object$theta[[m]] %*% t(object$basis[[m]]$psi) # / sqrt(nrow(object$basis[[m]]$psi))         ## removing the scaling we did in for theta
    psi[[m]] <- object$basis[[m]]$psi  ## list
    
  }
  
  ########################
  ## Run Rcpp Function
  ########################
  if(approx==TRUE){
    predict_hnew <- bsmim_predict_approx_cpp2(yz=cbind(object$y,object$z),
                                              Xlist=object$x,
                                              thetalist=object$theta,
                                              psilist=psi,
                                              rho=as.matrix(object$rho),
                                              gamma=as.matrix(object$gamma),
                                              lambdaInverse=object$lambdaInverse,
                                              lambdaBInverse=object$lambdaBInverse,
                                              sigma2=object$sigma2,
                                              weightslist=weights,
                                              gridpointslist=newX,
                                              poly=(1-object$gaussian),
                                              d=object$polydegree,
                                              randint=object$randint,
                                              Bmat=object$Bmat)
  }else{
    predict_hnew <- bsmim_predict_cpp2(yz=cbind(object$y,object$z),
                                       Xlist=object$x,
                                       thetalist=object$theta,
                                       psilist=psi,
                                       rho=as.matrix(object$rho),
                                       gamma=as.matrix(object$gamma),
                                       lambdaInverse=object$lambdaInverse,
                                       lambdaBInverse=object$lambdaBInverse,
                                       sigma2=object$sigma2,
                                       weightslist=weights,
                                       gridpointslist=newX,
                                       poly=(1-object$gaussian),
                                       d=object$polydegree,
                                       randint=object$randint,
                                       Bmat=object$Bmat)
  }
  
  ########################
  ## estimate contrasts
  ########################
  temp <- diag(1,length(predict_hnew$hmean))
  temp[,ncol(temp)] <- -1
  temp <- temp[-nrow(temp),]
  contrasts <- data.frame(mean=predict_hnew$hmean[-length(predict_hnew$hmean)]-predict_hnew$hmean[length(predict_hnew$hmean)],
                          sd=sqrt(diag(temp%*%predict_hnew$hcov%*%t(temp))) )
  contrasts$lower <- contrasts$mean - 1.96*contrasts$sd
  contrasts$upper <- contrasts$mean + 1.96*contrasts$sd
  ########################
  ## return values
  ########################
  # fits <- NULL
  # cc <- 1 ## column index
  # for(m in 1:M){
  #   for(l in 1:Lm[m]){
  #     temp <-
  #       data.frame(m = m,
  #                  l = l,
  #                  # name = names(object$x)[m],
  #                  mean = predict_hnew$hmean[,cc],
  #                  sd = predict_hnew$hsd[,cc],
  #                  grid = gridpoints[[m]][,l]*object$xscale$sd[[m]][l]+object$xscale$mean[[m]][l] )  ## grid of points, NOT SCALED OR TRANSFORMED APPROPRIATELY
  #     fits <- rbind(fits,temp)
  #     cc <- cc+1
  #   }
  # }
  # fits$lower <- fits$mean - 1.96*fits$sd
  # fits$upper <- fits$mean + 1.96*fits$sd
  # class(fits) <- "hpred"
  # return(fits)
  fits <- data.frame(mean=predict_hnew$hmean,
                     sd=sqrt(diag(predict_hnew$hcov)))
  fits$lower <- fits$mean - 1.96*fits$sd
  fits$upper <- fits$mean + 1.96*fits$sd
  fits <- fits[-nrow(fits),] ## remove the value used for comparisons
  
  results <- list(fits=fits,
                  contrasts=contrasts,
                  overall=overall,
                  grid=newX) ## rmse for test data
  class(results) <- "hassoc"
  return(results)
  
  
}








#' Predict hnew for BSMIM
#' 
#'
#' @param object An object of class bsmim
#' @param points Number of points to predict at
#' @param crossM Exposure to set a cross section quantile for
#' @param qtl quantile for the cross section
#' @param trueW weights for prediction for comparison to true curve; if NULL then set to posterior mean
#'
#' @return a dataframe containing predicted values.
#' @export
#'

predict_hnew_indexwise2 <- function(object,
                                   points=20, 
                                   qtl_lims=c(0.05,0.95),
                                   crossM=0, 
                                   qtl=0.5,
                                   trueW=NULL){
  
  ## dimensions
  M <- length(object$x)
  N <- length(object$y)
  
  
  ######################################################
  ## Construct new grid points, quantiles and weights
  ######################################################
  ## compute weights (basically A theta^*, so that X*theta*=x(Atheta*))
  ## compute matrix of Ebar=Xthetabar (used for index quantiles)
  ## compute grid of points for each index
  weights <- list()                         ## list of (S_iter by Lm) matrices with rows rho^{1/2}w (rows correspond to different iterations)
  Ebar <- matrix(NA,N,M)                    ## matrix of Ebar=Xthetabar, with columns corresponding to indices m
  gridpoints <- matrix(NA,points,M)         ## grid of points for hnew, for each m
  Xmedian <- list()                         ## list of L_m-vectors for median exposure level across N observations
  psi <- list()                             ## prep psi for Rcpp function
  for(m in 1:M){   ## rescale to divide by number of data times and rho
    weights[[m]] <- sqrt(as.matrix(object$rho)[,m])*object$theta[[m]] %*% t(object$basis[[m]]$psi) # / sqrt(nrow(object$basis[[m]]$psi))         ## removing the scaling we did in for theta
    if(!is.null(trueW)){
      Ebar[,m] <- apply(object$x[[m]] %*% t(trueW[[m]] %*% t(object$basis[[m]]$psi)),1,mean) # / sqrt(nrow(object$basis[[m]]$psi)))  ## compute X*theta_bar
    }else{
      theta_temp <- object$theta[[m]] ## NEW ERROR HANDLING: allow for var selection
      theta_temp[apply(theta_temp,1,function(x) sum(x!=0)==0),] <- 1/sqrt(ncol(theta_temp))
      Ebar[,m] <- apply(object$x[[m]] %*% t(theta_temp %*% t(object$basis[[m]]$psi)),1,mean) # / sqrt(nrow(object$basis[[m]]$psi)))  ## compute X*theta_bar
    }
    gridpoints[,m] <- seq(quantile(Ebar[,m],qtl_lims[1]),quantile(Ebar[,m],qtl_lims[2]), length=points)  ## set mth column to a grid of points form the 5th to 95th percentile, with evenly spaced points (not at percentiles)
    Xmedian[[m]] <- apply(object$x[[m]],2,median) 
    psi[[m]] <- object$basis[[m]]$psi  ## list
  }
  Eq <- apply(Ebar,2,function(x) quantile(x,qtl))
  
  
  ########################
  ## Run Rcpp Function
  ########################
  predict_hnew <- bsmim_predict_indexwise_cpp2(yz=cbind(object$y,object$z),
                                               Xlist=object$x,
                                               thetalist=object$theta,
                                               psilist=psi,
                                               rho=as.matrix(object$rho),
                                               gamma=as.matrix(object$gamma),
                                               lambdaInverse=object$lambdaInverse,
                                               lambdaBInverse=object$lambdaBInverse,
                                               sigma2=object$sigma2,
                                               weightslist=weights,
                                               gridpoints=gridpoints,
                                               Xmedianlist=Xmedian,
                                               Eq=Eq,
                                               crossM=crossM,
                                               poly=(1-object$gaussian),
                                               d=object$polydegree,
                                               randint=object$randint,
                                               Bmat=object$Bmat)
    

  
  ########################
  ## return values
  ########################
  fits <- NULL
  for(m in 1:M){
    temp <- 
      data.frame(m = m,
                 # name = names(object$x)[m],
                 mean = predict_hnew$hmean[,m],
                 sd = predict_hnew$hsd[,m],
                 grid = gridpoints[,m] ) ##gridpoints[,m]  )                                    ## grid of points, NOT SCALED OR TRANSFORMED APPROPRIATELY
                 # E = ((gridpoints[,m] * sum(object$basis[[m]]$psi%*%colMeans(object$theta[[m]])))-mean(object$x[[m]]%*%object$basis[[m]]$psi%*%colMeans(object$theta[[m]])))/sd(object$x[[m]]%*%object$basis[[m]]$psi%*%colMeans(object$theta[[m]])) #scale according to exposure.
    fits <- rbind(fits,temp)
  }
  
  fits$lower <- fits$mean - 1.96*fits$sd
  fits$upper <- fits$mean + 1.96*fits$sd
  if(!missing(crossM)){
    fits$cross_m <- crossM
    fits$cross <- names(object$x)[crossM]
    fits$qtl <- qtl
  }
  
  class(fits) <- "hpred_indexwise"
  return(fits)
}








#' Predict hnew for BSMIM for new exposure values
#' 
#'
#' @param object An object of class bsmim
#' @param newX M-list of (pts by Lm)-matrices containing new values to predict at---all others get set to 0; overrides qtl and points
#' @param newY N-vector of new outcomes
#' @param newZ N by P_z matrix of new covariates
#' @param newAmat N by K design matrix indicating cluster membership (only used when using random intercepts)
#' @return a list of (1) a dataframe containing mean, sd and interval for hnew, (2) matrix of samples of hnew, (3) list of exposure grid (newX)
#' @export
#'
predict_hnew_X2 <- function(object,
                            newX=NULL,
                            newY=NULL,
                            newZ=NULL,
                            newAmat=NULL,
                            rawX=TRUE,## rawX=TRUE indicates that newX needs to be preprocessed
                            qtl_lims=c(0.05,0.95),## only used if newX is NULL
                            points=20 ){    ## only used if newX is NULL
                           
  
  
  ########################
  ## data checks
  ########################
  ## dimensions
  N <- length(object$y)               ## number of observations, N      
  M <- length(object$x)               ## number of indices, M
  Lm <- unlist(lapply(object$x,ncol)) ## M-vector of index sizes L_m
  
  ### check newX
  if(is.null(newX)){ ### if no newX given, construct same grid as the other versions above

    ## Construct new grid points, quantiles and weights
    newX <- list()                      ## list of ((points x sum(Lm)) by Lm ) matrices whose columns are a grid of points for hnew
    id <- 1 ## iterator to set grid appropriately
    for(m in 1:M){   
        temp_Xq <- matrix(0,nrow=(points*sum(Lm)),ncol=Lm[m]) ## since empty component of list causes problems #Xqlist[[m]] <- c()
        for(l in 1:Lm[m]){
          temp_Xq[,l] <- quantile(object$x[[m]][,l],0.5)  ## set elements to median
          temp_Xq[points*(id-1)+(1:points),l] <- seq(quantile(object$x[[m]][,l],qtl_lims[1]),quantile(object$x[[m]][,l],qtl_lims[2]), length=points) ## set grid
          id <- id+1 ## increase iterator
        }
        newX[[m]] <- temp_Xq
    }
    
  }else{
    if(rawX==TRUE){ ## if data havent been pre-processed
      for(m in 1:M){
        newX[[m]] <- (newX[[m]]-matrix(fit$xscale$mean[[m]],ncol=Lm[m],nrow=nrow(as.matrix(newX[[m]])),byrow=TRUE))/matrix(fit$xscale$sd[[m]],ncol=Lm[m],nrow=nrow(as.matrix(newX[[m]])),byrow=TRUE)
      }
    }
  }
  if(length(newX)!=M){
    stop("newX must be length M")
  }
  points <- nrow(as.matrix(newX[[1]]))
  for(m in 1:M){
    if((nrow(as.matrix(newX[[m]]))!=points) | (ncol(as.matrix(newX[[m]]))!=Lm[m])){
      stop("newX has wrong dimension")
    }else{
      newX[[m]] <- as.matrix(newX[[m]])
    }
  }
  
  ######################################################
  ## Construct weights (basically A theta^*, so that X*theta*=x(Atheta*))
  ######################################################
  weights <- list() ## list of (S_iter by Lm) matrices with rows rho^{1/2}w (rows correspond to different iterations)
  psi <- list()     ## prep psi for Rcpp function
  for(m in 1:M){    ## rescale to divide by number of data times and rho
    weights[[m]] <- sqrt(as.matrix(object$rho)[,m])*object$theta[[m]] %*% t(object$basis[[m]]$psi) # / sqrt(nrow(object$basis[[m]]$psi))         ## removing the scaling we did in for theta
    psi[[m]] <- object$basis[[m]]$psi  ## list
    
  }
  
  
  ########################
  ## Run Rcpp Function
  ########################
  predict_hnew <- bsmim_predict_X_cpp2(yz=cbind(object$y,object$z),
                                       Xlist=object$x,
                                       thetalist=object$theta,
                                       psilist=psi,
                                       rho=as.matrix(object$rho),
                                       gamma=as.matrix(object$gamma),
                                       lambdaInverse=object$lambdaInverse,
                                       lambdaBInverse=object$lambdaBInverse,
                                       sigma2=object$sigma2,
                                       weightslist=weights,
                                       gridpointslist=newX,
                                       poly=(1-object$gaussian),
                                       d=object$polydegree,
                                       randint=object$randint,
                                       Bmat=object$Bmat)
  
  
  ########################
  ## Check model fit
  ########################
  LPPD <- RMSE <- NULL
  if(!is.null(newY) & !is.null(newZ)){              ## only if new Y and Z values are provided
    if(nrow(as.matrix(newY))==points & nrow(newZ)==points){

      ## draw random intercepts if necessary 
      ### NOTE: currently using lppd for NEW clusters, i.e. drawing b_(K+1) from N(0,\sigma^2_B), rather than for existing clusters (this is to keep consistent with predictions in cross validation)
      if(object$randint==TRUE & !is.null(newAmat)){
        bk <- mvrnorm(ncol(newAmat),mu=rep(0,length(object$sigma2)),Sigma=diag(object$lambdaBInverse*object$sigma2))
        bki <- newAmat%*%bk ## repeat for cluster members
      }else{
        bki <- matrix(0,nrow=length(newY),ncol=length(object$sigma2))
      }
      
      ## Compute LPPD for test data
      lppd <- log(apply( exp( -0.5*log(2*pi) -0.5*log(matrix(object$sigma2,nrow=length(object$sigma2),ncol=length(newY))) -0.5*t((matrix(newY,nrow=nrow(as.matrix(newY)),ncol=length(object$sigma2))-t(predict_hnew$hsamp)-bki-newZ%*%t(object$gamma))^2)/object$sigma2  ),2,mean  )) ## -2 log[(1/S)\sum_{s}^S P(y|\theta^s)] 
      LPPD <- sum(lppd)
      
      ## Compute RMSE for test data
      rmse <- sqrt(apply((matrix(newY,nrow=nrow(as.matrix(newY)),ncol=length(object$sigma2))-t(predict_hnew$hsamp)-bki-newZ%*%t(object$gamma))^2,1,mean))
      RMSE <- sum(rmse)
      
    }
  }
  
  
  ### posterior mean of predicted outcome
  pred_out <- NULL
  mean_fitted <- NULL
  if(!is.null(newZ)){
    if(nrow(newZ)==points){
      
      ## draw random intercepts if necessary 
      if(object$randint==TRUE & !is.null(newAmat)){
        bk <- mvrnorm(ncol(newAmat),mu=rep(0,length(object$sigma2)),Sigma=diag(object$lambdaBInverse*object$sigma2))
        bki <- newAmat%*%bk ## repeat for cluster members
      }else{
        bki <- matrix(0,nrow=length(newY),ncol=length(object$sigma2))
      }
      
      ## distribution of h + Z*posterior mean of gamma +bki (doesnt include an error term for predictions)
      predY <- t(predict_hnew$hsamp)+newZ%*%t(object$gamma)+bki # predict_hnew$hmean
      pred_out <- data.frame(mean=apply(predY,1,mean),
                             sd=apply(predY,1,sd),
                             lower=apply(predY,1,function(x) quantile(x,0.025)),
                             upper=apply(predY,1,function(x) quantile(x,0.975))) #
 
    }
    ## mean fitted val = hmean + Z%*%gammamean ## should be same as mean of pred_out
    mean_fitted <- predict_hnew$hmean+newZ%*%apply(object$gamma,2,mean) 
  }
  
  
  
  ########################
  ## return values
  ########################
  
  fits <- data.frame(mean=predict_hnew$hmean,
                     sd=sqrt(diag(predict_hnew$hcov)))
  fits$lower <- fits$mean - 1.96*fits$sd
  fits$upper <- fits$mean + 1.96*fits$sd
  
  
  results <- list(fits=fits,
                  hsamp=predict_hnew$hsamp,
                  grid=newX,
                  LPPD=LPPD, ## lppd for test data
                  RMSE=RMSE, ## rmse for test data
                  mean_fitted=mean_fitted, ## mean fitted val for outcomes (if newZ provided)
                  pred_out=pred_out) ##  fitted outcomes (if newZ provided)
  class(results) <- "hpred"
  return(results)
}





