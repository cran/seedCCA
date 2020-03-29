\name{Pm}
\alias{Pm}
\title{Projection of a seed matrix on to the column subspace of M with respect to \code{Sx} inner-product}

\description{
The function reuturns a projection of a seed matrix on to the column subspace of M with respect to \code{Sx} inner-product.}

\usage{
Pm(M, Sx, seed)
}
\arguments{
  \item{M}{numeric matrix (p * k), a basis matrix of the column space of M}
  \item{Sx}{a inner-product matrix}
  \item{seed}{seed matrix}
}

\examples{
data(cookie)
myseq<-seq(141,651,by=2)
X<-as.matrix(cookie[-c(23,61),myseq])
Y<-as.matrix(cookie[-c(23,61),701:704])

## using cov(X,Y) as a seed matrix
seed <- cov(X,Y)
col.num <- dim(seed)[2]
M <- iniCCA(X, Y, u=4, num.d=col.num)
Sx <- cov(X)
Pm(M, Sx, seed)

## using the first 2 largest eigenvectors of cov(X,Y) as a seed matrix
seed2 <- svd(cov(X,Y))$u[,1:2]
M2 <- iniCCA(X, Y, u=4, num.d=2)
Pm(M, Sx, seed2)
}

