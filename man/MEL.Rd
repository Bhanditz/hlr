\name{MEL}
\alias{MEL}
\alias{plot.MEL}
\alias{print.MEL}
\alias{summary.MEL}
\title{Logistic Regression Using Maximum Estimated Likelihood}
\description{
  Logistic regression using maximum estimated likelihood (MEL) 
  in order to cope with separation issues.
}
\usage{
MEL(x, y, delta = 0.01, epsilon = 1e-06, maxit = 100)
\method{print}{MEL}(x, ...)
\method{summary}{MEL}(object, ...)
\method{plot}{MEL}(x, ...)
}
\arguments{
  \item{x}{design matrix (n, p) for function \code{MEL}; object of class
           'MEL' for the print and plot methods}
  \item{y}{response vector}
  \item{delta}{constant}
  \item{epsilon}{precision constant}
  \item{maxit}{maximum number of iterations}
  \item{object}{object of class 'MEL'}
  \item{...}{further arguments to be passed to the methods}
}
\details{
  The \code{MEL} function fits the MEL-model to the data.
  The \code{print} method displays the model coefficients.
  The \code{summary} method displays the model coefficients and
  displays the names of the components of the MEL output object.
  The \code{plot} function plots the observed response values against
    the predicted values (using the MEL linear predictor) on link scale.
    The actual logistic cdf is added as well.
}
\value{
Object of class 'MEL' with following components  
  \item{MEL}{MEL estimates of the coefficients}
  \item{outMEL}{object of class 'glm' corresponding to the
     final fit}
}
\references{Rousseeuw, P.J. and Christmann, A. (2003). 
  Robustness against separation and outliers in binary regression.
  Computational Statistics \& Data Analysis, \bold{43}, 315 -- 332.
  
  Original S-PLUS code available at
  \url{http://www.stoch.uni-bayreuth.de/de/CHRISTMANN}
}
\author{Tobias Verbeke, largely based on original S-PLUS code by 
  Peter J. Rousseeuw and Andreas Christmann}
\examples{
  par(mfrow = c(2,1))
  ### Example 1 for function MEL: data set has overlap
  set.seed(314)
  n <- 500
  beta <- matrix(c(2), ncol=1)
  x <- matrix(rnorm(n), ncol=1)
  eta <- -3 + x %*% beta
  y <- rbinom(nrow(x), 1, plogis(eta))
  out <- MEL(x, y)
  print(out)
  summary(out)
  plot(out)
  title("overlap")

  ### Example 2 for function MEL: data set has no overlap
  y[eta <= -1] <- 0
  y[eta > -1] <- 1
  out <- MEL(x,y, delta=0.01)
  out
  plot(out)
  title("no overlap")

  ### Example 3 for function MEL: all responses are equal to 0
  y <- rep(0, nrow(x))
  MEL(x,y)

  ### Example 4 for function MEL: all responses are equal to 1
  y <- rep(1, nrow(x))
  out <- MEL(x,y)
  out 
}
\keyword{models}
