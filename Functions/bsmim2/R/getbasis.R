require(refund)

#-------------------------------------------------------------------------------------------------
#' Build BDLIM basis
#'
#' \code{bdlimbasis} builds the basis for thge BDLIM model using either FACE of natural cubic splines.
#' @param X A n x p matrix of observations. Each row has the observations for one individual. The observations should be evenly spaced over a grid.
#' @param basis.opts List with the entries: type = the type of basis used, either 'face' (default) or "ns" or "bs" for splines or "gam" for presmoothing the exposure with a gam following defaults from mgcv; knots = the number of knots used for method face; pve = the percent of variance explained by the PCs for method face; df = the df for ns method.
#' @seealso \code{\link{fpca.face}}, \code{\link{ns}}, \code{\link{bs}}, and \code{\link{gam}} for details on the methods called.
#' @author Glen McGee & Ander Wilson
#' @importFrom splines ns bs
#' @importFrom refund fpca.face
#' @importFrom mgcv gam


getbasis <- function(X,basis.opts){
  
  
  #account for missing input in basis function
  if(toupper(basis.opts$type)=="NONE"){
    
    return(list(psi=diag(ncol(X)), type=basis.opts$type))
    
  }else if(toupper(basis.opts$type)=="PCA"){
    
    #correct pve if needed.
    if(is.null(basis.opts$pve) | !is.numeric(basis.opts$pve)){
      basis.opts$pve=.99
    }else{
      if(basis.opts$pve<0 | basis.opts$pve>1) basis.opts$pve <- .99
    }
    pca <- prcomp(X)
    nbases <- min(which(summary(pca)$importance[3,]>basis.opts$pve)) ## number of PCs to include
    return(list(psi=as.matrix(pca$rotation[,1:nbases]), eigenvalues=pca$sdev[1:nbases]^2, pve=basis.opts$pve,type=basis.opts$type))
    
  }else if(toupper(basis.opts$type)=="NS"){
    
    if(is.null(basis.opts$df) & !is.null(basis.opts$knots)){
      basis.opts$df <- basis.opts$knots+2
    }else if(is.null(basis.opts$df)){ 
      basis.opts$df <- 5
    }
    
    B1 <- ns(seq(1,ncol(X)),df=basis.opts$df, intercept=TRUE)
    X <-  B1 %*% qr.solve(B1,t(X)) 
    svdX <- svd(X)
    
    return(list(psi=svdX$u[,1:ncol(B1)], lambda=svdX$d[1:ncol(B1)], pve=basis.opts$pve,type=basis.opts$type))
    
  }else if(toupper(basis.opts$type)=="BS"){
    
    if(is.null(basis.opts$df)) basis.opts$df <- round(ncol(X)/5)
    B1 <- bs(seq(1,ncol(X)),df=basis.opts$df, intercept=TRUE)
    X <-  B1 %*% qr.solve(B1,t(X)) 
    svdX <- svd(X)
    
    return(list(psi=svdX$u[,1:ncol(B1)], lambda=svdX$d[1:ncol(B1)], pve=basis.opts$pve,type=basis.opts$type))
    
  }else if(toupper(basis.opts$type)=="FACE"){
    
    #correct pve if needed.
    if(is.null(basis.opts$pve) | !is.numeric(basis.opts$pve)){
      basis.opts$pve=.99
    }else{
      if(basis.opts$pve<0 | basis.opts$pve>1) pve <- .99
    }
    
    if(is.null(basis.opts$knots)) basis.opts$knots <- round(ncol(X)/3)
    fitted <- fpca.face(X, knots=basis.opts$knots,pve=basis.opts$pve)
    lambda <- fitted$evalues
    return(list(psi=fitted$efunctions, lambda=lambda, knots=basis.opts$knots,pve=basis.opts$pve,type=basis.opts$type))
    
  }else if(toupper(basis.opts$type)=="GAM"){
    
    #correct pve if needed.
    if(is.null(basis.opts$pve) | !is.numeric(basis.opts$pve)){
      basis.opts$pve=.99
    }else{
      if(basis.opts$pve<0 | basis.opts$pve>1) pve <- .99
    }
    
    
    for(i in 1:nrow(X)){
      X[i,] <- predict(gam(X[i,]~s(seq(1:ncol(X)))))
    }
    
    eigcorX <- eigen(cor(X))
    nbases <- min(which(cumsum(eigcorX$values)/sum(eigcorX$values)>basis.opts$pve))
    return(list(psi=eigcorX$vectors[,1:nbases], eigenvalues=eigcorX$values[1:nbases], pve=basis.opts$pve, df=basis.opts$df,type=basis.opts$type))
    
    
  }else{
    
    stop("basis type not recognized.")
    
  }
  
  
}