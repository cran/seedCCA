\name{seeding.auto.stop}
\alias{seeding.auto.stop}
\title{increments of iterative projections with automatic stopping}
\description{
Returns increments (nFu) of iterative projections of a seed matrix onto a covariance matrix upto k, which properly chosen by satisfying the terminating condition \code{eps} (\code{eps} can be selected by users).}

\usage{
seeding.auto.stop(seed, covx, n, u.max=30, eps=0.01)
}

\arguments{
  \item{seed}{numeric matrix (p * d), a seed matrix}
  \item{covx}{numeric matrix (p * p), covariance matrix of X}
  \item{n}{numeric, sample sizes}
  \item{u.max}{numeric, maximum number of projection.  The default value is equal to 30.}
  \item{eps}{numeric, a value of a condition for terminating the projection.  The default value is equal to 0.01.}
  }

\value{
  \item{nFu}{n*Fu values}
  \item{u}{the number of projection properly chosen by satisfying the terminating condition \code{eps}}
  }

\examples{
data(cookie)
myseq<-seq(141,651,by=2)
X<-as.matrix(cookie[-c(23,61),myseq])
Y<-as.matrix(cookie[-c(23,61),701:704])
seed <- cov(X,Y)
covx <- cov(X)
seeding.auto.stop(seed, covx, n=dim(X)[1])
seeding.auto.stop(seed, covx, n=dim(X)[1], u.max=20, eps=0.001)
}

