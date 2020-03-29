\name{print.seedCCA}
\alias{print.seedCCA}
\title{basic function for printing class "seedCCA"}
\description{The function controls to print class "seedCCA".
The function prints the estimated coefficents, if they exist.
For subclass "finalCCA", canonical correlations are additionally reported.
For subsclass "selectu", increments, suggested number of projections and the values of \code{type} and \code{eps} are reported.
}

\usage{
\method{print}{seedCCA}(x,\dots)
}

\arguments{
\item{x}{The name of an object of class "seedCCA"}
\item{\dots}{arguments passed to the print.method}
}

\examples{
########  data(cookie) ########
data(cookie)
myseq<-seq(141,651,by=2)
X<-as.matrix(cookie[-c(23,61),myseq])
Y<-as.matrix(cookie[-c(23,61),701:704])

fit.seed2 <- seedCCA(X,Y)
fit.seed2
print(fit.seed2)
}
