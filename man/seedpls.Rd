\name{seedpls}
\alias{seedpls}
\title{Partial least squares through iterative projections}
\description{Returns partial least squares estimates through iterative projections. And, the function results in subclass "seedpls".}

\usage{seedpls(X, Y, u=5, scale=FALSE)}

\arguments{
\item{X}{numeric matrix (n * p), a set of predictors}
\item{Y}{numeric vector or matrix (n * r), responses (it can be multi-dimensional)}
\item{u}{numeric, the number of projections. The default is 5.}
\item{scale}{logical, FALSE is default. If TRUE, each predictor is standardized with mean 0 and variance 1}
}
\value{
\item{coef}{the estimated coefficients for each iterative projection upto u}
\item{u}{the maximum number of projections}
\item{X}{Predictors}
\item{Y}{Response}
\item{scale}{status of scaling predictors}
}

\examples{
########  data(cookie) ########
data(cookie)
myseq<-seq(141,651,by=2)
X<-as.matrix(cookie[-c(23,61),myseq])
Y<-as.matrix(cookie[-c(23,61),701:704])

fit.pls1 <- seedpls(X,Y[,1]) ## one-dimensional response
fit.pls2 <- seedpls(X,Y, u=6, scale=TRUE) ## four-dimensional response

}
