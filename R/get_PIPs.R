
#' Get PIPs
#'
#' Get main effect and second order interaction posterior inclusion probabilities.
#'
#' @param fit Fitted model
#'
#' @return List of PIP objects.
#'
#' @export
get_PIPs <- function(fit){

  mainPIPs <- apply(fit$gamma_additive,2,mean)

  interPIPs <- apply(fit$gamma_nonadditive,2,mean)

  twowayPIPs <- matrix(NA,ncol=fit$p,nrow=fit$p)
  for(ii in (1:(fit$p-1))){
    for(jj in ((ii+1):fit$p)){
      twowayPIPs[ii,jj] <- mean(fit$gamma_nonadditive[,ii]*fit$gamma_nonadditive[,jj])
    }
  }

  return(list(mainPIPs=mainPIPs,
              interPIPs=interPIPs, ## any non-additivity
              twowayPIPs=twowayPIPs)) ## simultaneous presence in kernel

}

#' Get thetas
#'
#' @keywords internal
get_thetas <- function(fit){

  p <- fit$p
  Lm <- fit$Lm

  theta <- w <- vector(mode="list",length=p)
  for(j in 1:p){
    theta_id <- sum(Lm[0:(j-1)])+(1:Lm[j])

    theta[[j]] <- fit$theta[,theta_id] ## get theta

    w[[j]] <- fit$theta[,theta_id]%*%t(fit$Psi[[j]]) ## apply basis
    w[[j]] <- w[[j]]/apply(w[[j]],1,function(x) sqrt(sum(x^2))) ## re standardize
    w[[j]][is.na(w[[j]])] <- 0 ## deal with NAs
  }



  return(list(theta=theta,
              w=w))
}


#' Summarize weights
#'
#' Summarize thetas and ws
#'
#' @param fit Fitted model
#'
#' @return List of theta and w summaries
#'
#' @export
summarize_weights <- function(fit){

  weights <- get_thetas(fit)
  theta <- weights$theta
  w <- weights$w


  theta_sum <- lapply(theta,function(thet){

    thet <- thet[!apply(thet,1,function(x) sum(x!=0)==0),] ## theta not well defined if gamma=0

    res <- rbind(apply(thet,2,mean)/sqrt(sum(apply(thet,2,mean)^2)), ## making sure post mean is standardized
                 apply(thet,2,function(x) stats::quantile(x,0.025)),
                 apply(thet,2,function(x) stats::quantile(x,0.975)))
    rownames(res) <- c("Mean","2.5th","97.5th")
  })

  w_sum <- lapply(theta,function(thet){

    thet <- thet[!apply(thet,1,function(x) sum(x!=0)==0),] ## theta not well defined if gamma=0

    res <- rbind(apply(thet,2,mean)/sqrt(sum(apply(thet,2,mean)^2)), ## making sure post mean is standardized
                 apply(thet,2,function(x) stats::quantile(x,0.025)),
                 apply(thet,2,function(x) stats::quantile(x,0.975)))
    rownames(res) <- c("Mean","2.5th","97.5th")
  })

  return(list(theta=theta_sum,
              w=w_sum))

}

#' Plot weights
#'
#' Plot thetas and ws
#'
#' @param fit Fitted model
#'
#' @return List of theta and w plots.
#'
#' @export
plot_weights <- function(fit){

  weights <- get_thetas(fit)
  theta <- weights$theta
  w <- weights$w


  ## theta
  plots_theta <- list()
  for(j in 1:length(theta)){

    dfplot <- data.frame(theta[[j]])
    colnames(dfplot) <- 1:ncol(dfplot)
    dfplot <- tidyr::gather(dfplot,"l","theta")
    if(ncol(theta[[j]])==1){dfplot$l <- ""}
    dfplot$l <- factor(dfplot$l,levels=unique(dfplot$l)) ## order x axis

    p <- ggplot2::ggplot(dfplot, ggplot2::aes_string(x="l", y="theta"))+
      ggplot2::geom_boxplot()+##color=alpha("black",0.8)
      ggplot2::geom_hline(yintercept=0,linetype=2,alpha=0.8)+
      ggplot2::scale_y_continuous(expression(paste("Weights (",theta,")")),limits=c(-1,1)) +
      ggplot2::scale_x_discrete("Component") +
      ggplot2::ggtitle("") +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank())

    plots_theta[[j]] <- p
  }

  ## w
  plots_w <- list()
  for(j in 1:length(w)){

    dfplot <- data.frame(w[[j]])
    colnames(dfplot) <- 1:ncol(dfplot)
    dfplot <- tidyr::gather(dfplot,"l","w")
    if(ncol(w[[j]])==1){dfplot$l <- ""}
    dfplot$l <- factor(dfplot$l,levels=unique(dfplot$l)) ## order x axis

    p <- ggplot2::ggplot(dfplot, ggplot2::aes_string(x="l", y="w"))+
      ggplot2::geom_boxplot()+##color=alpha("black",0.8)
      ggplot2::geom_hline(yintercept=0,linetype=2,alpha=0.8)+
      ggplot2::scale_y_continuous(expression(paste("Weights (",w,")")),limits=c(-1,1)) +
      ggplot2::scale_x_discrete("Component") +
      ggplot2::ggtitle("") +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank())

    plots_w[[j]] <- p
  }




  return(list(theta=plots_theta,
              w=plots_w))

}


