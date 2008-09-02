# MEL: maximum estimated likelihood estimator
# Copyright (C) Peter J. Rousseeuw, Andreas Christmann, 2001
# Please note, that there is NO WARRANTY, 18/JUN/2001

"MEL" <- function(x, y, delta=0.01, epsilon=1.0E-6, maxit=100) {
  pihat <- max(delta, min(1-delta, mean(y)))
  delta0 <- (pihat*delta) / (1+delta)
  delta1 <- (1+pihat*delta) / (1+delta)
  ytilde <- delta0 * (1-y) + delta1 * y
  response <- cbind(ytilde, 1-ytilde)
  suppressWarnings(outMEL <- glm(response ~ x, family=binomial, 
          control=glm.control(epsilon=epsilon, maxit=maxit)))
  beta <- coef(outMEL)
  ### construct output object
  output <- list(MEL=beta, outMEL=outMEL)
  class(output) <- "MEL"
  output$call <- match.call()
  output
}

print.MEL <- function(x, ...) {
	cat("MEL:\n")
	print(x$MEL, ...)
	invisible(x)
}

summary.MEL <- function(object, ...) {
	cat("MEL:\n")
	print(object$MEL, ...)
	cat("\nAvailable arguments:\n")
	print(names(object), ...)
	invisible(object)
}

plot.MEL <- function(x, ...) {
  if (is.logical(x$outMEL) && is.na(x$outMEL)) { 
    warning("Plot not available") 
  }
  dotList <- list(...)
  # defaults from original S-PLUS code (except for pch=23)
  dotList$xlab <- if(!hasArg(xlab)) "x'MEL" else dotList$xlab
  dotList$ylab <- if(!hasArg(ylab))"y" else dotList$ylab
  dotList$ylim <- if(!hasArg(ylim)) c(0, 1) else dotList$ylim
  # add data
  dotList$x <- x$outMEL$linear.predictors
  dotList$y <- x$outMEL$y
  
  if (is.list(x$outMEL)) {
    do.call("plot", dotList)
    tmp <- seq(min(x$outMEL$linear.predictors), max(x$outMEL$linear.predictors),0.05)
    lines(tmp, plogis(tmp), col = "grey")
    invisible(x)
  }
}


