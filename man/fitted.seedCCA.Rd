\name{fitted.seedCCA}
\alias{fitted.seedCCA}
\title{Fitted values of ordinary and partial least squares}
\description{Returns fitted values of ordinary and partial least squares through iterative projections. It works only for subclasses "seedols" and "seedpls".}

\usage{
\method{fitted}{seedCCA}(object, u=NULL,\dots)
}


\arguments{
\item{object}{The name of an object of class "seedCCA"}
\item{u}{numeric, the number of projections. The default is NULL. This option is valid for PLS alone. The option returns the fitted values for u projections. For example, if it is specified at k, then the fitted values with k projections are returned.}
\item{\dots}{arguments passed to the fitted.method}
}

\examples{
########  data(cookie) ########
########  data(cookie) ########
data(cookie)
myseq<-seq(141,651,by=2)
X<-as.matrix(cookie[-c(23,61),myseq])
Y<-as.matrix(cookie[-c(23,61),701:704])

fit.ols1 <- seedCCA(X[,1:4], Y[,1], type="cca")
fit.pls1 <- seedCCA(X, Y[,1], type="pls")
fit.pls2 <- seedCCA(X, Y[,1], type="pls", scale=TRUE)

fitted(fit.ols1)
fitted(fit.pls1)
fitted(fit.pls1, u=4)
fitted(fit.pls2, u=4)
}
