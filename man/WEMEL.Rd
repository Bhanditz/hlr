\name{WEMEL}
\alias{WEMEL}
\alias{plot.WEMEL}
\alias{print.WEMEL}
\alias{summary.WEMEL}
\title{Robust Estimator in the Hidden Logistic Regression Model}
\description{
  Logistic regression using weighted maximum estimated likelihood (WEMEL) 
  in order to cope with separation issues and outliers.
}
\usage{
WEMEL(x, x1, y, delta=0.01, q=0.75, method="MCD", w=rep(1,length(y)), epsilon = 1e-06, maxit = 100)
\method{print}{WEMEL}(x, ...)
\method{summary}{WEMEL}(object, ...)
\method{plot}{WEMEL}(x, which = c(1, 2), ...)
}
\arguments{
  \item{x}{design matrix (n, p) for function \code{WEMEL}; object of class
           'WEMEL' for the print and plot methods}
  \item{x1}{(sub-)matrix of the design matrix. The robust weights are computed
            w.r.t. to \code{x1}.  E.g. \code{x1} contains only continuous 
            variables.}
  \item{y}{response vector}
  \item{delta}{constant}
  \item{q}{quantile used for MCD and for the robust weights; defaults to 0.75}
  \item{method}{method to define weights; one of "MCD" (default) or "PCA"}
  \item{w}{input vector of weights}
  \item{epsilon}{precision constant for the algorithm, default: 1.E-6}
  \item{maxit}{maximum number of iterations for the algorithm; defaults to 100}
  \item{object}{object of class 'WEMEL'}
  \item{which}{which plot should be plotted ? An index plot of the robust
  weights (\code{which=1}), the observed response values against
    the predicted values (\code{which=2}) or both (\code{which=c(1, 2)}, the
    default).} 
  \item{...}{further arguments to be passed to the methods}
}
\details{
  The \code{WEMEL} function fits the WEMEL-model to the data.
  The \code{print} method displays the model coefficients.
  The \code{summary} method displays the model coefficients and
  displays the names of the components of the WEMEL output object.
  The \code{plot} function plots either the index plot of the robust
  weights (\code{which=1}) or the observed response values against
    the predicted values (using the WEMEL linear predictor) on link 
    scale with a logistic cdf overplotted (\code{which=2}). The
    default value (\code{which=c(1, 2)}) plots both.
}
\value{
Object of class 'WEMEL' with following components  
  \item{WEMEL}{WEMEL estimates of the coefficients}
  \item{outWEMEL}{object of class 'glm' corresponding to the
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
  ### Example 1 for function WEMEL: data set has overlap
  set.seed(314)
  n <- 500
  beta <- matrix(c(3), ncol=1)
  x <- matrix(rnorm(n), ncol=1)
  eta <- -2 + x %*% beta
  y <- rbinom(nrow(x), 1, plogis(eta))
  out <- WEMEL(x, x, y)
  print(out)
  summary(out)
  plot(out)
  
  ### Example 2 for function WEMEL: 2 bad leverage points
  x[499:500] <- c(-10, 10)
  y[499:500] <- c(1, 0)
  out <- WEMEL(x, x, y, delta=0.01)
  out
  plot(out)
  
  ### Example 3 for function WEMEL: data set has no overlap
  eta <- -2 + x
  y[eta <= -1] <- 0
  y[eta > -1] <- 1
  out <- WEMEL(x, x, y, delta=0.01)
  out
  plot(out)  
}
\keyword{models}
