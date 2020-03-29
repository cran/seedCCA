\name{iniCCA}
\alias{iniCCA}
\title{Initialized CCA in seeded CCA}
\description{Returns the canonical coefficient matrices from the initialized step in seeded CCA. The initialzied CCA is done only for the first set in its first argument. The \code{num.d} must be less than or equal to the dimension of the second set.}

\usage{iniCCA(X, Y, u, num.d)}

\arguments{
\item{X}{numeric matrix (n * p), the first set of variables: this set of variables alone is reduced.}
\item{Y}{numeric matrix (n * r), the second set of variables}
\item{u}{numeric, the terminiating index of the projection}
\item{num.d}{numeric, the first "num.d" eigenvectors of cov(X,Y) to replace cov(X,Y), if min(p,r) relatively bigger than n. The \code{num.d} must be less than or equal to min(p,r).}
}
\value{
\item{B}{the initialized CCAed coefficient matrix projected by the value of \code{u}}
}
\examples{
########  data(cookie) ########
data(cookie)
myseq<-seq(141,651,by=2)
X<-as.matrix(cookie[-c(23,61),myseq])
Y<-as.matrix(cookie[-c(23,61),701:704])
min.pr <- min( dim(X)[2], dim(Y)[2])
MX0 <- iniCCA(X, Y, u=4, num.d=min.pr)
ini.X <- X\%*\%MX0

########  data(nutrimouse) ########
data(nutrimouse)
Y<-as.matrix(nutrimouse$lipid)
X<-as.matrix(nutrimouse$gene)
MX0 <- iniCCA(X, Y, u=4, num.d=4)
MY0 <- iniCCA(Y, X, u=5, num.d=4)
ini.X <- X \%*\% MX0
ini.Y <- Y \%*\% MY0
}
