#' Compute the relatedness summary parameters (RSPs) from a pedigree object or twice the kinship matrix
#'
#' This function computes the RSPs from the eigenvalues of twice the kinship matrix either directly or via a pedigree object.
#'
#' @param ped a pedigree object from the 'kinship2' package.  One of \code{ped} or \code{kmat} must be specified, but not both.
#' @param kmat a kinship matrix
#' @param thresholderror a scalar numeric of a permissible threshold that the mean eigenvalue is allowed to deviate from 1.
#' @param forcePD logical of whether positive definiteness of kmat should be strictly enforced.
#' @return a pedRSP object with the calculated RSPs.
#' @author Jesse D. Raffa
#' @details
#' This function calculates the eigenvalues from a pedigree object or twice the kinship matrix,
#' and using these computes the RSPs: the variance of the eigenvalue and log eigenvalues, and the geometric mean.
#' @seealso \code{\link[kinship2]{pedigree}} \code{\link[kinship2]{kinship}}
#' @export
#' @importFrom kinship2 pedigree
#' @importFrom kinship2 kinship
#' @examples
#' require(kinship2)
#' # Example from kinship2 package
#' test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
#' mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
#' dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
#' sex =c(1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 2, 2, 2))
#' tped <- with(test1, pedigree(id, dad, mom, sex))
#' computeRSPs(ped=tped)

computeRSPs <- function(ped=NULL,kmat=NULL,thresholderror=0.01,forcePD=TRUE) {
  if(is.null(ped) & is.null(kmat))
    stop("Must pass pedigree or kinship matrix")
  if(!is.null(ped) & (class(ped)!="pedigree" & class(ped)!="pedigreeList"))
    stop("ped argument not a pedigree")
  if(!is.null(kmat) & (class(kmat)!="matrix" & class(kmat)!="Matrix" & class(kmat)!="dsCMatrix"))
    stop('kmat argument not a matrix')
  if(!is.null(ped) & !is.null(kmat))
    stop("please specify either ped or kmat argument, not both")
  if(is.null(kmat))
    kmat <- 2*kinship(ped)
  if(abs(mean(diag(kmat))-1)>thresholderror) {
    message(paste("Mean eigenvalue:",mean(diag(kmat))))
    stop(paste("This should be around one.  If this is correct, set thresholderror higher than",thresholderror))
  }
  eg <- eigen(kmat)
  if(min(eg$val)<=0 & forcePD) {
    message("Kinship matrix is not positive definite. If this is correct, and you would like to proceed, set forcePD=FALSE")
    stop("Note: Most RSPs will not be finite")
  }
  lambda <- eg$val;
  Vlambda <- var(lambda)
  Vq <- var(log(lambda))
  gammapar <- exp(mean(log(lambda)) );
  ret <- list(lambda=lambda,N=length(lambda),Vlambda=Vlambda,Vq=Vq,gammapar=gammapar)
  class(ret) <- "pedRSP"
  return(ret)

}



#' @export
print.pedRSP = function(x, digits=3,...) {
  message(paste("Relatedness Summary Parameters"))
  message(paste0("N=",x$N))
  message(paste0("V_lambda=",round(x$Vlambda,digits=digits)))
  message(paste0("V_q=",round(x$Vq,digits=digits)))
  message(paste0("gamma=",round(x$gammapar,digits=digits)))
  Vqs <- mean(log(c(1.5,0.5))^2) - mean(log(c(1.5,0.5)))^2
  Nr <- ceiling(( (x$N-1)*x$Vq/Vqs + 1)/2)
  message(paste0("Effective number of sibpairs: ",Nr, " (pairs) or ", Nr*2, " individuals"))
}


#' Compute the exact theoretical expected likelihood ratio test (ELRT) and power from a pedRSP object.
#'
#' This function computes the exact theoretical expectation of the likelihood ratio test (ELRT) statistic from a pedRSP object and uses this to calculate power. This function uses all N eigenvalues to compute the ELRT.
#'
#'
#' @param pedRSP a pedRSP object computed using the computeRSPs function
#' @param truehsq a vector of at least length 1 of the true heritability.  Exactly one of truehsq and null must have length 1.
#' @param null a vector of at least length 1 of the null hypothesis heritability.
#' @param sig.level a scalar of the significance level of the hypothesis test
#' @return a pedELRT object with the approximated ELRT and power.
#' @author Jesse D. Raffa
#' @details
#' This function calculates the approximate ELRT/power from a pedRSP object using only the actual study size and the variance of the eigenvalues.
#' This approximation is generally seen to over estimate the ELRT and thus power, and should not be used in practice.
#' @seealso \code{\link[pedRSP]{computeRSPs}}
#' @export
#' @examples
#' require(kinship2)
#' test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
#' mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
#' dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
#' sex =c(1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 2, 2, 2))
#' tped <- with(test1, pedigree(id, dad, mom, sex))
#' x <- computeRSPs(ped=tped)
#' ex.elrt <- computeEXELRT(pedRSP=x,truehsq=0.3,null=seq(0.1,0.8,0.1))
#' plot(ex.elrt)

computeEXELRT <- function(pedRSP,truehsq,null,sig.level=0.05) {
  if(class(pedRSP)!='pedRSP')
    stop("pedRSP argument must be pedRSP object")
  if(length(truehsq)!=1 & length(null)!=1 )
    stop("One of truehsq and null must have length 1")
  eg <- pedRSP$lambda
  if(length(null)>1) {
    ELRT <- rep(NA,length(null))
    for(j in 1:length(null)) {
     ELRT[j] <-  1 + -sum(log(1+(eg-1)*truehsq)) - pedRSP$N + sum(log(1+(eg-1)*null[j])) + sum((1+truehsq*(eg-1))/(1+(eg-1)*null[j]))
    }

  } else {
    ELRT <- rep(NA,length(truehsq))
    for(j in 1:length(truehsq)) {
      ELRT[j] <- 1 + -sum(log(1+(eg-1)*truehsq[j])) - pedRSP$N + sum(log(1+(eg-1)*null)) + sum((1+truehsq[j]*(eg-1))/(1+(eg-1)*null))
    }
  }
  power <- 1-pchisq(qchisq(1-sig.level,df=1),1,ncp=(ELRT-1))
  if(length(truehsq)>1) {
    grid <- 'truehsq'
    h <- truehsq;
  } else {
    grid <- "null"
    h <- null;

  }
  ret <- list(type="EXELRT",obj=pedRSP,lambda=eg,N=pedRSP$N,grid=grid,h=h,ELRT=ELRT,power=power,truehsq=truehsq,null=null,sig.level=sig.level)
  class(ret) <- "pedELRT"
  return(ret)
}



#' Compute an approximation (AELRTA) to the expected likelihood ratio test (ELRT) and power from a pedRSP object.
#'
#' This function computes an approximation (AELRTA) to the expectation of the likelihood ratio test statistic from a pedRSP object and uses this to calculate power. This function uses the variance of the eigenvalues to compute the ELRT.
#'
#'
#' @param pedRSP a pedRSP object computed using the computeRSPs function
#' @param truehsq a vector of at least length 1 of the true heritability.  Exactly one of truehsq and null must have length 1.
#' @param null a vector of at least length 1 of the null hypothesis heritability.
#' @param sig.level a scalar of the significance level of the hypothesis test
#' @return a pedELRT object with the approximated ELRT and power.
#' @author Jesse D. Raffa
#' @details
#' This function calculates the approximate ELRT/power from a pedRSP object using only the actual study size and the variance of the eigenvalues.
#' This approximation is generally seen to over estimate the ELRT and thus power, and should not be used in practice.
#' @seealso \code{\link[pedRSP]{computeRSPs}}
#' @export
#' @examples
#' require(kinship2)
#' test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
#' mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
#' dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
#' sex =c(1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 2, 2, 2))
#' tped <- with(test1, pedigree(id, dad, mom, sex))
#' x <- computeRSPs(ped=tped)
#' elrta <- computeAELRTA(pedRSP=x,truehsq=0.3,null=seq(0.1,0.8,0.1))
#' plot(elrta)

computeAELRTA <- function(pedRSP,truehsq,null,sig.level=0.05) {
  if(class(pedRSP)!='pedRSP')
    stop("pedRSP argument must be pedRSP object")
  if(length(truehsq)!=1 & length(null)!=1 )
    stop("One of truehsq and null must have length 1")
  ELRT <- pedRSP$Vlambda*(pedRSP$N-1)/2*(truehsq-null)^2 +1
  power <- 1-pchisq(qchisq(1-sig.level,df=1),1,ncp=(ELRT-1))
  if(length(truehsq)>1) {
    grid <- 'truehsq'
    h <- truehsq;
  } else {
    grid <- "null"
    h <- null;

  }
  ret <- list(type="ALERTA",obj=pedRSP,Vlambda=pedRSP$Vlambda,N=pedRSP$N,grid=grid,h=h,ELRT=ELRT,power=power,truehsq=truehsq,null=null,sig.level=sig.level)
  class(ret) <- "pedELRT"
  return(ret)
}

#' Compute an approximation (AELRTB) to the expected likelihood ratio test (ELRT) and power from a pedRSP object.
#'
#' This function computes an approximation (AELRTB) to the expectation of the likelihood ratio test statistic from a pedRSP object and uses this to calculate power. This function uses the variance of the log-eigenvalues and geometric mean of eigenvalues to compute the ELRT.
#'
#'
#' @param pedRSP a pedRSP object computed using the computeRSPs function
#' @param truehsq a vector of at least length 1 of the true heritability.  Exactly one of truehsq and null must have length 1.
#' @param null a vector of at least length 1 of the null hypothesis heritability.
#' @param sig.level a scalar of the significance level of the hypothesis test
#' @return a pedELRT object with the approximated ELRT and power.
#' @author Jesse D. Raffa
#' @details
#' This function calculates the approximate ELRT/power from a pedRSP object using only the actual study size and the variance of the log-eigenvalues.
#' This approximation can both over and under estimate the exact ELRT depending on the values of truehsq and null.  See paper for more details
#' @seealso \code{\link[pedRSP]{computeRSPs}}
#' @export
#' @examples
#' require(kinship2)
#' test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
#' mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
#' dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
#' sex =c(1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 2, 2, 2))
#' tped <- with(test1, pedigree(id, dad, mom, sex))
#' x <- computeRSPs(ped=tped)
#' elrtb <- computeAELRTB(pedRSP=x,truehsq=0.3,null=seq(0.1,0.8,0.1))
#' plot(elrtb)

computeAELRTB <- function(pedRSP,truehsq,null,sig.level=0.05) {
  if(class(pedRSP)!='pedRSP')
    stop("pedRSP argument must be pedRSP object")
  if(length(truehsq)!=1 & length(null)!=1 )
    stop("One of truehsq and null must have length 1")
  ggtrue <- 1+ truehsq*(pedRSP$gammapar-1);
  ggnulls <- 1+ null*(pedRSP$gammapar-1);
  ELRT <- pedRSP$N*( - log(ggtrue/ggnulls) + ggtrue/ggnulls -1) + ((pedRSP$N-1)*pedRSP$Vq*pedRSP$gammapar/2)*( -null*(null-1)/ggnulls^2 + truehsq*(truehsq-1)/ggtrue^2 - (null*pedRSP$gammapar + null-1)*(truehsq-null)/ggnulls^3) + 1
  power <- 1-pchisq(qchisq(1-sig.level,df=1),1,ncp=(ELRT-1))
  if(length(truehsq)>1) {
    grid <- 'truehsq'
    h <- truehsq;
  } else {
    grid <- "null"
    h <- null;

  }
  ret <- list(type="ALERTB",obj=pedRSP,Vq=pedRSP$Vq,N=pedRSP$N,grid=grid,h=h,ELRT=ELRT,power=power,truehsq=truehsq,null=null,sig.level=sig.level)
  class(ret) <- "pedELRT"
  return(ret)
}

#' Compute an approximation (AELRTC) to the expected likelihood ratio test (ELRT) and power from a pedRSP object.
#'
#' This function computes an approximation (AELRTC) to the expecation of the likelihood ratio test statistic from a pedRSP object and uses this to calculate power. This function uses the variance of the log-eigenvalues to compute the ELRT.
#'
#'
#' @param pedRSP a pedRSP object computed using the computeRSPs function
#' @param truehsq a vector of at least length 1 of the true heritability.  Exactly one of truehsq and null must have length 1.
#' @param null a vector of at least length 1 of the null hypothesis heritability.
#' @param sig.level a scalar of the significance level of the hypothesis test
#' @return a pedELRT object with the approximated ELRT and power.
#' @author Jesse D. Raffa
#' @details
#' This function calculates the approximate ELRT/power from a pedRSP object using only the actual study size and the variance of the log-eigenvalues.
#' This approximation can both over and under estimate the exact ELRT depending on the values of truehsq and null.  See paper for more details
#' @seealso \code{\link[pedRSP]{computeRSPs}}
#' @export
#' @examples
#' require(kinship2)
#' test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
#' mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
#' dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
#' sex =c(1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 2, 2, 2))
#' tped <- with(test1, pedigree(id, dad, mom, sex))
#' x <- computeRSPs(ped=tped)
#' elrtc <- computeAELRTC(pedRSP=x,truehsq=0.3,null=seq(0.1,0.8,0.1))
#' plot(elrtc)

computeAELRTC <- function(pedRSP,truehsq,null,sig.level=0.05) {
  if(class(pedRSP)!='pedRSP')
    stop("pedRSP argument must be pedRSP object")
  if(length(truehsq)!=1 & length(null)!=1 )
    stop("One of truehsq and null must have length 1")
  ELRT <- pedRSP$Vq*(pedRSP$N-1)/2*(truehsq-null)^2 +1
  power <- 1-pchisq(qchisq(1-sig.level,df=1),1,ncp=(ELRT-1))
  if(length(truehsq)>1) {
    grid <- 'truehsq'
    h <- truehsq;
  } else {
    grid <- "null"
    h <- null;

  }
  ret <- list(type="ALERTC",obj=pedRSP,Vq=pedRSP$Vq,N=pedRSP$N,grid=grid,h=h,ELRT=ELRT,power=power,truehsq=truehsq,null=null,sig.level=sig.level)
  class(ret) <- "pedELRT"
  return(ret)
}


#' Compute an approximation (AELRTESS) to the expected likelihood ratio test (ELRT) and power from a pedRSP object.
#'
#' This function computes an approximation (AELRTESS) to the expecation of the likelihood ratio test statistic from a pedRSP object and uses this to calculate power. This function uses the effective number of sibpairs approximation.
#'
#'
#' @param pedRSP a pedRSP object computed using the computeRSPs function
#' @param truehsq a vector of at least length 1 of the true heritability.  Exactly one of truehsq and null must have length 1.
#' @param null a vector of at least length 1 of the null hypothesis heritability.
#' @param sig.level a scalar of the significance level of the hypothesis test
#' @return a pedELRT object with the approximated ELRT and power.
#' @author Jesse D. Raffa
#' @details
#' This function calculates the approximate ELRT/power from a pedRSP object using only the actual study size and the variance of the log-eigenvalues.
#' This approximation is usually a good approximation for intermediate levels of heritability.  When heritability is small or large, it can be overly conservative (under estimate power/ELRT).  See paper for more details
#' @seealso \code{\link[pedRSP]{computeRSPs}}
#' @export
#' @examples
#' require(kinship2)
#' test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
#' mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
#' dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
#' sex =c(1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 2, 2, 2))
#' tped <- with(test1, pedigree(id, dad, mom, sex))
#' x <- computeRSPs(ped=tped)
#' elrtess <- computeAELRTESS(pedRSP=x,truehsq=0.3,null=seq(0.1,0.8,0.1))
#' plot(elrtess)

computeAELRTESS <- function(pedRSP,truehsq,null,sig.level=0.05) {
  if(class(pedRSP)!='pedRSP')
    stop("pedRSP argument must be pedRSP object")
  if(length(truehsq)!=1 & length(null)!=1 )
    stop("One of truehsq and null must have length 1")
  Vqs <- mean(log(c(1.5,0.5))^2) - mean(log(c(1.5,0.5)))^2
  Nr <- ceiling(( (pedRSP$N-1)*pedRSP$Vq/Vqs + 1)/2)
  tmpped <- list();
  tmpped$N <- Nr*2;
  tmpped$lambda <- rep(c(1.5,0.5),Nr)
  class(tmpped) <- "pedRSP";
  ELRT <- computeEXELRT(tmpped,truehsq=truehsq,null=null,sig.level=sig.level)$ELRT;
  power <- 1-pchisq(qchisq(1-sig.level,df=1),1,ncp=(ELRT-1))
  if(length(truehsq)>1) {
    grid <- 'truehsq'
    h <- truehsq;
  } else {
    grid <- "null"
    h <- null;

  }
  ret <- list(type="ALERTESS",obj=pedRSP,Vq=pedRSP$Vq,N=pedRSP$N,grid=grid,h=h,ELRT=ELRT,power=power,truehsq=truehsq,null=null,sig.level=sig.level)
  class(ret) <- "pedELRT"
  return(ret)
}


#' Function for plotting a pedELRT object
#'
#' This function plots two graphs side-by-side of the ELRT and the power of the tests specified in the computeELRT functions
#' @param x a pedELRT object
#' @param inc.exact a logical of whether to include the exact ELRT in the plot (default=FALSE)
#' @param legend a logical of whether to include a legend in the plot (default=FALSE)
#' @param ... Additional arguments.
#' @seealso \code{\link[pedRSP]{computeAELRTA}} \code{\link[pedRSP]{computeAELRTB}} \code{\link[pedRSP]{computeAELRTC}} \code{\link[pedRSP]{computeAELRTESS}} \code{\link[pedRSP]{computeEXELRT}}
#' @author Jesse D. Raffa
#' @details
#' Plotting function for pedELRT object.  It plots the approximate or exact ELRT from what is calculated in the computeELRTX and computeEXLERT functions.
#' @examples
#' require(kinship2)
#' test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
#' mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
#' dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
#' sex =c(1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 2, 2, 2))
#' tped <- with(test1, pedigree(id, dad, mom, sex))
#' x <- computeRSPs(ped=tped)
#' elrtess <- computeAELRTESS(pedRSP=x,truehsq=0.3,null=seq(0.1,0.8,0.1))
#' plot(elrtess,inc.exact=TRUE,legend=TRUE)
#' @export
plot.pedELRT = function(x,inc.exact=FALSE,legend=FALSE,...) {
  par(mfcol=c(1,2))
  plot(x$h,x$ELRT,xlab=x$grid,ylab="ELRT",main=paste("ELRT:",x$type),type="l",ylim=c(0,max(x$ELRT)))
  if(inc.exact & x$type!="EXELRT") {
    exelrt <- computeEXELRT(x$obj,truehsq=x$truehsq,null=x$null,sig.level=x$sig.level)
    lines(exelrt$h,exelrt$ELRT,col="orange")
  }
  plot(x$h,x$power,xlab=x$grid,ylab="Power",main=paste("Power:",x$type),type="l",ylim=c(0,max(x$power)))
  if(inc.exact & x$type!="EXELRT") {
    lines(exelrt$h,exelrt$power,col="orange")
  }
  if(legend & inc.exact & x$type!="EXELRT")
    legend(x$h[which.min(x$power)]-0.02,max(x$power),lty=1,col=c("black","orange"),c(x$type,"EXELRT"))
}
