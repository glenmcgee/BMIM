require(ggplot2)

#' Plot component-wise response curves
#' 
#'
#' @param object An object of class bsmim
#'
#' @return a list of trace plots.
#' @export
#'

plot_trace2 <- function(object){
  
  if(class(object)!="bsmim"){
    stop("input must be of class bsmim")
  }
  
  ## theta
  plots_theta <- list()
  for(mm in 1:length(object$theta)){                 ## loop over m (indices)
    plots_m <- list()
    for(ll in 1:ncol(object$theta[[mm]])){     ## loop over l (components)
      
      dfplot <- data.frame(theta=object$theta[[mm]][,ll],
                           iter=1:nrow(object$theta[[mm]]) )
      
      p <- ggplot(dfplot, aes(x=iter, y=theta))+
        geom_line(linetype=1)+ ## 3 is dotted
        ylab("Theta")+
        xlab("Iteration")+
        ggtitle(paste("Theta: Index",mm,", Component",ll))+
        theme_bw() +
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
      plots_m[[ll]] <- p
    }
    plots_theta[[mm]] <- plots_m
  }
  
  ## rho
  plots_rho <- list()
  for(mm in 1:ncol(as.matrix(object$rho))){                 ## loop over m (indices)

    dfplot <- data.frame(rho=as.matrix(object$rho)[,mm],
                         iter=1:nrow(as.matrix(object$rho)) )
    
    p <- ggplot(dfplot, aes(x=iter, y=rho))+
      geom_line(linetype=1)+ ## 3 is dotted
      ylab("Rho")+
      xlab("Iteration")+
      ggtitle(paste("Rho: Index",mm))+
      theme_bw() +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
    
    plots_rho[[mm]] <- p
  }
  
  ## nu
  	plots_nu <- list()
  	for(mm in 1:length(object$nu)){                 ## loop over m (indices)
      plots_m <- list()
      for(ll in 1:ncol(object$nu[[mm]])){     ## loop over l (components)
        
        dfplot <- data.frame(nu=object$nu[[mm]][,ll],
                             iter=1:nrow(object$nu[[mm]]) )
        
        p <- ggplot(dfplot, aes(x=iter, y=nu))+
          geom_line(linetype=1)+ ## 3 is dotted
          ylab("Nu")+
          xlab("Iteration")+
          ggtitle(paste("Nu: Index",mm,", Component",ll))+
          theme_bw() +
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
        plots_m[[ll]] <- p
      }
      plots_nu[[mm]] <- plots_m
    }


  
  ## gamma
  plots_gamma <- list()
  for(pp in 1:ncol(object$gamma)){                 ## loop over P covariates
    
    dfplot <- data.frame(gamma=object$gamma[,pp],
                         iter=1:nrow(object$gamma) )
    
    p <- ggplot(dfplot, aes(x=iter, y=gamma))+
      geom_line(linetype=1)+ ## 3 is dotted
      ylab("Gamma")+
      xlab("Iteration")+
      ggtitle(paste("Gamma: ",pp))+
      theme_bw() +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
    
    plots_gamma[[pp]] <- p
  }
  
  ## lambdaInverse
  dfplot <- data.frame(lambdaInverse=object$lambdaInverse,
                       iter=1:length(object$lambdaInverse) )
  plot_lambdaInverse <- ggplot(dfplot, aes(x=iter, y=lambdaInverse))+
    geom_line(linetype=1)+ ## 3 is dotted
    ylab("lambda Inverse")+
    xlab("Iteration")+
    ggtitle(paste("lambda Inverse"))+
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
    
  ## lambdaBInverse
  dfplot <- data.frame(lambdaBInverse=object$lambdaBInverse,
                       iter=1:length(object$lambdaBInverse) )
  plot_lambdaBInverse <- ggplot(dfplot, aes(x=iter, y= lambdaBInverse))+
    geom_line(linetype=1)+ ## 3 is dotted
    ylab("lambda_B Inverse")+
    xlab("Iteration")+
    ggtitle(paste("lambda_B Inverse"))+
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
  
  ## tau
    dfplot <- data.frame(tau=object$tau,
                         iter=1:length(object$tau) )
    plot_tau <- ggplot(dfplot, aes(x=iter, y=tau))+
      geom_line(linetype=1)+ ## 3 is dotted
      ylab("Tau")+
      xlab("Iteration")+
      ggtitle(paste("Tau"))+
      theme_bw() +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 

  
  
  ## sigma2
  dfplot <- data.frame(sigma2=object$sigma2,
                       iter=1:length(object$sigma2) )
  plot_sigma2 <- ggplot(dfplot, aes(x=iter, y=sigma2))+
    geom_line(linetype=1)+ ## 3 is dotted
    ylab("sigma2")+
    xlab("Iteration")+
    ggtitle(paste("sigma2"))+
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
  
  
  return(list(theta=plots_theta,
              rho=plots_rho,
              nu=plots_nu,
              tau=plot_tau,
              gamma=plots_gamma,
              lambdaInverse=plot_lambdaInverse,
              sigma2=plot_sigma2))
}





