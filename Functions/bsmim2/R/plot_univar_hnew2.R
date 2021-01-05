require(ggplot2)


#' Plot component-wise response curves
#' 
#'
#' @param pred An object of class hpred
#' @param ylims An optional vector of length 2 containing the limits for the y axis; default is null
#'
#' @return a list of component-wise response-curve plots .
#' @export
#'
plot_univar_hnew2 <- function(pred,ylims=NULL,assoc=TRUE){
  
  if(class(pred)=="hassoc" & assoc==TRUE){
    df <- pred$contrasts
  }else{
    df <- pred$fits
  }
  
  if(class(pred)=="hassoc"){
    if(pred$overall==TRUE){
      df$grid <- seq(0.05,0.95,by=0.05)
      
      
      ## set ylimits
      if(length(ylims)!=2){
        miny <- min(df$lower);maxy <- max(df$upper)  
      }else if(ylims[1]>=ylims[2]){
        miny <- min(df$lower);maxy <- max(df$upper)  
      }else{
        miny <- ylims[1];maxy <- ylims[2]
      }
      
      
      p <- ggplot(df, aes(x=grid, y=mean))+
        geom_line(linetype=3)+ ## 3 is dotted
        geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.2)+
        ylim(miny,maxy)+
        ylab("Estimated exposure-response vs median")+
        xlab("Quantile")+
        ggtitle(paste("Overall"))+
        theme_bw() +
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
      
      return(p)
    }
  }
  
  
  ncomp <- 0 ## number of components
  for(mm in 1:length(pred$grid)){
    for(ll in 1:ncol(pred$grid[[mm]])){
      ncomp <- ncomp+1
    }
  }
  

  df$grid <- 0
  pts <- nrow(pred$fits)/ncomp
  nid <- 1 ## indexing
  for(mm in 1:length(pred$grid)){
    for(ll in 1:ncol(pred$grid[[mm]])){
      df$grid[(nid-1)*pts + (1:pts)] <- pred$grid[[mm]][(nid-1)*pts + (1:pts), ll]
      nid <- nid+1
    }
  }
  # ## make dataframe from object of class hpred
  # df <- data.frame(matrix(unlist(pred), nrow=length(pred$m), byrow=F))
  # colnames(df) <- names(pred)
  
  ## set ylimits
  if(length(ylims)!=2){
    miny <- min(df$lower);maxy <- max(df$upper)  
  }else if(ylims[1]>=ylims[2]){
    miny <- min(df$lower);maxy <- max(df$upper)  
  }else{
    miny <- ylims[1];maxy <- ylims[2]
  }
  
  nid <- 1
  plots <- list()
  for(mm in 1:length(pred$grid)){                 ## loop over m (indices)
    plots_m <- list()
    for(ll in 1:ncol(pred$grid[[mm]])){     ## loop over l (components)
      
      # dfplot <- df[df$m==mm & df$l==ll,]  ## subset data
      dfplot <- df[(nid-1)*pts + (1:pts),]
      nid <- nid+1
      
      p <- ggplot(dfplot, aes(x=grid, y=mean))+
        geom_line(linetype=3)+ ## 3 is dotted
        geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.2)+
        ylim(miny,maxy)+
        ylab("Estimated exposure-response (h)")+
        xlab("Exposure Component")+
        ggtitle(paste("Index",mm,", Component",ll))+
        theme_bw() +
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
      plots_m[[ll]] <- p
    }
    plots[[mm]] <- plots_m
  }
  
  return(plots)
}




#' #' Old version: Plot component-wise response curves (only difference is the format of pred$grid)
#' #' 
#' #'
#' #' @param pred An object of class hpred
#' #' @param ylims An optional vector of length 2 containing the limits for the y axis; default is null
#' #'
#' #' @return a list of component-wise response-curve plots .
#' #' @export
#' #'
#' plot_univar_hnew_old2 <- function(pred,ylims=NULL){
#'   
#'   ## make dataframe from object of class hpred
#'   df <- data.frame(matrix(unlist(pred), nrow=length(pred$m), byrow=F))
#'   colnames(df) <- names(pred)
#'   
#'   ## set ylimits
#'   if(length(ylims)!=2){
#'     miny <- min(df$lower);maxy <- max(df$upper)  
#'   }else if(ylims[1]>=ylims[2]){
#'     miny <- min(df$lower);maxy <- max(df$upper)  
#'   }else{
#'     miny <- ylims[1];maxy <- ylims[2]
#'   }
#'   
#'   plots <- list()
#'   for(mm in 1:max(df$m)){                 ## loop over m (indices)
#'     plots_m <- list()
#'     for(ll in 1:max(df$l[df$m==mm])){     ## loop over l (components)
#'       
#'       dfplot <- df[df$m==mm & df$l==ll,]  ## subset data
#'       
#'       p <- ggplot(dfplot, aes(x=grid, y=mean))+
#'         geom_line(linetype=3)+ ## 3 is dotted
#'         geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.2)+
#'         ylim(miny,maxy)+
#'         ylab("Estimated exposure-response (h)")+
#'         xlab("Exposure Component")+
#'         ggtitle(paste("Index",mm,", Component",ll))+
#'         theme_bw() +
#'         theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
#'       plots_m[[ll]] <- p
#'     }
#'     plots[[mm]] <- plots_m
#'   }
#'   
#'   return(plots)
#' }

# ### Example 
# pred_exp <- predict_hnew(test,points=50)
# pp <- plot_univar_hnew(pred_exp)
# pp














#' Plot component-wise response curves
#' 
#'
#' @param pred An object of class hpred_indexwise
#' @param ylims An optional vector of length 2 containing the limits for the y axis; default is null
#'
#' @return a list of component-wise response-curve plots .
#' @export
#'
plot_univar_hnew_indexwise2 <- function(pred,ylims=NULL){
  
  ## make dataframe from object of class hpred
  df <- data.frame(matrix(unlist(pred), nrow=length(pred$m), byrow=F))
  colnames(df) <- names(pred)
  
  
  ## set ylimits
  if(length(ylims)!=2){
    miny <- min(df$lower);maxy <- max(df$upper)  
  }else if(ylims[1]>=ylims[2]){
    miny <- min(df$lower);maxy <- max(df$upper)  
  }else{
    miny <- ylims[1];maxy <- ylims[2]
  }
  
  plots <- list()
  for(mm in 1:max(df$m)){                 ## loop over m (indices)
    
      dfplot <- df[df$m==mm,]  ## subset data
      
      p <- ggplot(dfplot, aes(x=grid, y=mean))+
        geom_line(linetype=3)+ ## 3 is dotted
        geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.2)+
        ylim(miny,maxy)+
        ylab("Estimated exposure-response (h)")+
        xlab("Exposure")+
        ggtitle(paste("Index",mm))+
        theme_bw() +
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
      
      plots[[mm]] <- p
    }

  
  return(plots)
}


# ### Example 
# pred_index <- predict_hnew_indexwise(test,points=50)
# pp <- plot_univar_hnew_indexwise(pred_index)
# pp


