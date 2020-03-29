\name{seedols}
\alias{seedols}
\title{Ordinary least squares}
\description{Returns ordinary least squares estimates. And, the function results in subclass "seedols". For this function to work, either X or Y has to be one-dimensional. It is not necessary that X and Y should be predictors and response, respectively. Regardless of the position in the arguments, the one-dimensional and multi-dimensional variables become response and predictors, respectively.}

\usage{seedols(X, Y)}

\arguments{
\item{X}{numeric vector or matrix, a first set of variables}
\item{Y}{numeric vector or matrix, a second set of variables}
}
\value{
\item{coef}{the estimated coefficients for each iterative projection upto u}
\item{X}{X, the first set}
\item{Y}{Y, the second set}
}

\examples{
########  data(cookie) ########
data(cookie)
myseq<-seq(141,651,by=2)
X<-as.matrix(cookie[-c(23,61),myseq])
Y<-as.matrix(cookie[-c(23,61),701:704])

ols1 <- seedols(X[,1:4],Y[,1])
ols2 <- seedols(Y[,1],X[,1:4])
## ols1 and ols2 are the same results.
}
