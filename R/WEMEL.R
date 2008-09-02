### utility function for the function ovWEMEL.HLR
### compute robust distances RD_i and robust weights depending on RD_i
### Copyright (C) Peter J. Rousseeuw, Andreas Christmann, 2001
### Please note, that there is NO WARRANTY, 25/JUN/2001

my.robweights <- function(x, q = 0.75) { # TV 2008-07-02: removed ask = F argument (not used)
  n <- nrow(as.matrix(x))                             #        as well as id.n = 3 argument and ... (not used) 
  p <- ncol(as.matrix(x))
  mcd.obj <- cov.mcd(x, quantile.used = floor(q*n)) # TV: adapted to MASS cov.mcd 
  # X <- mcd.obj$X #  TV: according to S-PLUS help page same as input data
  #                    --> replaced with x in file
  center <- mcd.obj$center
  cov <- mcd.obj$cov
  MD <- matrix(rep(0, n), ncol = 1)
  RD <- matrix(rep(0, n), ncol = 1)
  w <- matrix(rep(1,n), ncol = 1)
  if ( p > 0 ) {
    MD <- sqrt(mahalanobis(x, apply(x, 2, mean), n/(n - 1) * cov.wt(x)$cov))
    RD <- matrix(sqrt(mahalanobis(x, center, cov)), ncol = 1)
    M <- quantile(RD^2, probs=q, na.rm = FALSE)
    w <- M / ( apply(cbind(RD^2, rep(M,n)), 1, max) ) 
  }
  output <- list(RD=RD, MD=MD, weights=w)
  return(output)
}

### WEMEL: robust estimator in the hidden logistic regression model
### compute robust distances RD_i and robust weights depending on RD_i
### Used for S-Plus Version 4.5, Windows NT
### Copyright (C) Peter J. Rousseeuw, Andreas Christmann, 2001
### Please note, that there is NO WARRANTY, 25/JUN/2001
#---------------------------------------------------------------------
#   argument        meaning
#   --------        -------
#   x               design matrix (n,p)
#   x1              (sub-)matrix of the design matrix. The robust weights are computed
#                   w.r.t. to x1.  E.g. x1 contains only continuous variables.
#   y               response vector (n,1), values: 0 or 1
#   delta           constant, default: 0.01
#   q               quantile used for MCD and for the robust weights, default: 0.75
#   method          method to define weights, default: "MCD", alternative "PCA"
#   w               input vector of weights
#   epsilon         precision constant for the algorithm, default: 1.E-6
#   maxit           maximum number of iterations of the algorithm, default: 100
#   output          results (list)
#----------------------------------------------------------

"WEMEL" <- function(x, x1, y, delta=0.01, q=0.75, method="MCD", 
    w=rep(1,length(y)), epsilon=1.0E-6, maxit=100) {
  pihat <- max(delta, min(1-delta,mean(y)))
  delta0 <- (pihat*delta) / (1+delta)
  delta1 <- (1+pihat*delta) / (1+delta)
  ytilde <- delta0*(1-y) + delta1*y
  response <- cbind(ytilde,1-ytilde)
  if (method == "MCD") { 
    robust.weights <- my.robweights(x1, q=q)$weights 
  }
  if (method == "PCA") { 
    robust.weights <- w 
  }
  suppressWarnings(outWEMEL <- glm(response ~ x, family=binomial, weights=robust.weights,
          control=glm.control(epsilon=epsilon, maxit=maxit)))
  beta <- outWEMEL$coef
  ### construct output object
  output <- list(WEMEL=beta, outWEMEL=outWEMEL, robust.weights=robust.weights)
  class(output) <- "WEMEL"
  output$call <- match.call()
  output
}

print.WEMEL <- function(x, ...) {
  cat("WEMEL:\n")
  print(x$WEMEL, ...)
  invisible(x)
}

summary.WEMEL <- function(object, ...) {
  cat("WEMEL:\n")
  print(object$WEMEL, ...)
  cat("\nAvailable arguments:\n")
  print(names(object), ...)
  invisible(object)
}

plot.WEMEL <- function(x, which = c(1, 2), ...) {
  if (is.logical(x$outWEMEL) && is.na(x$outWEMEL)) { 
    warning("Plot not available") 
  }
  
  show <- c(FALSE, FALSE)
  show[which] <- TRUE
  
  par(mfrow = c(length(which), 1))
  
  dotList <- list(...)
  # defaults from original S-PLUS code (except for pch=23)
  dotList$xlab <- if(!hasArg(xlab)) "x'WEMEL" else dotList$xlab
  dotList$ylab <- if(!hasArg(ylab))"y" else dotList$ylab
  dotList$ylim <- if(!hasArg(ylim)) c(0, 1) else dotList$ylim
  dotList$main <- if(!hasArg(main)) "Fitted regression curve" else dotList$main
  
  # add data
  dotList$x <- x$outWEMEL$linear.predictors
  dotList$y <- x$outWEMEL$y
  
  if (is.list(x$outWEMEL)) {
  
    if (show[1]){ # Index plot of robust weights
      plot(x$robust.weights, main="Index plot of robust weights", 
           ylim=c(0,1), ylab="w")
    }
    if (show[2]){ # observed vs. predicted
      do.call("plot", dotList)
      tmp <- seq(min(x$outWEMEL$linear.predictors), max(x$outWEMEL$linear.predictors),0.05)
      lines(tmp, plogis(tmp))
    }
    invisible(x)
  }
}

### Examples
#graphsheet()
#par(mfrow=c(2,1))
### Example 1 for function MEL: data set has overlap
#set.seed(314)
#n <- 500
#beta <- matrix(c(3),ncol=1)
#x <- matrix(rnorm(n),ncol=1)
#eta <- -2 + x %*% beta
#y <- rbinom(nrow(x),1,plogis(eta))
#out <- WEMEL(x,x,y)
#print(out)
#summary(out)
#plot(out)
#
### Example 2 for function MEL: 2 bad leverage points
#x[499:500] <- c(-10,10)
#y[499:500] <- c(1,0)
#out <- WEMEL(x,x,y, delta=0.01)
#print(out)
#plot(out)

### Example 3 for function MEL: data set has no overlap
#y[eta <= -1] <- 0
#y[eta > -1] <- 1
#out <- WEMEL(x,x,y, delta=0.01)
#print(out)
#plot(out)

