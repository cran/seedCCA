\name{seeding}
\alias{seeding}
\title{Iterative projection of a seed matrix onto cov(X)}
\description{
Returns iterative projection of a seed matrix on cov(X) upto k times.}

\usage{
seeding(seed, cov.x, n, u=5)
}

\arguments{
  \item{seed}{numeric matrix (p * d), a seed matrix}
  \item{cov.x}{numeric matrix (p * p), cov(X)}
  \item{n}{sample sizes}
  \item{u}{maximum number of projection}
}

\value{
  \item{F.u_x}{n*Fu values of X}
  \item{F.u_y}{n*Fu values of Y}
}

\examples{
data(cookie)
myseq<-seq(141,651,by=2)
X<-as.matrix(cookie[-c(23,61),myseq])
Y<-as.matrix(cookie[-c(23,61),701:704])
seed <- cov(X,Y)
cov.x <- cov(X)
seeding(seed, cov.x, n=dim(X)[1], u=4)
}

