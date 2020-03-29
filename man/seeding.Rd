\name{seeding}
\alias{seeding}
\title{increments of iterative projections}
\description{
Returns increments (nFu) of iterative projections of a seed matrix onto a covariance matrix \code{u} times.)}

\usage{
seeding(seed, covx, n,  u=10)
}

\arguments{
  \item{seed}{numeric matrix (p * d), a seed matrix}
  \item{covx}{numeric matrix (p * p), covariance matrix of X}
  \item{n}{numeric, sample sizes}
  \item{u}{numeric, maximum number of projections}
}

\value{
  \item{nFu}{n*Fu values}
}

\examples{
data(cookie)
myseq<-seq(141,651,by=2)
X<-as.matrix(cookie[-c(23,61),myseq])
Y<-as.matrix(cookie[-c(23,61),701:704])
seed <- cov(X,Y)
covx <- cov(X)
seeding(seed, covx, n=dim(X)[1], u=4)
}

