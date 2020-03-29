\name{plot.seedCCA}
\alias{plot.seedCCA}
\title{Plotting class "seedCCA" depending on the value of \code{type}}

\description{
The function is for plotting class "seedCCA". Depending on subclass defined by the value of \code{type}, its resulting plot is different.
}

\usage{
\method{plot}{seedCCA}(x, ref=90, eps=0.01, \dots)
}

\arguments{
\item{x}{The name of an object of class "seedCCA"}
\item{ref}{numeric, the value for reference line. It must be chosen between 0 and 100. It works only for subclass "finalCCA".}
\item{eps}{numeric, a value to terminate projections. It must be chosen between 0 and 1. The default is equal to 0.01. It works only for subclass "seedpls".}
\item{\dots}{arguments passed to the plot.method}
}

\details{
subclass "finalCCA": the function makes a plot for percents of cumulative canonical correlations.

subclass "seedpls": the function returns a proper number of projections and plot of the projection increment against the number of projections. A proper number of projections is indicated with a blue-color vertical bar in the plot. Only for subclass "seedpls", the output is retured. See Value part.

subclass "seedols": No plotting

subclass "selectu": the function makes a plot for increment of iterative projections by the output of subclass "selectu".
}
\value{
  \item{proper.u}{proper value of the number of projections}
  \item{nFu}{incrments (n*Fu) of the iterative projection.}
  \item{u}{the maximum number of projections from "seedpls" object}
  \item{eps}{a value for terminating the projection.  The default value is equal to 0.01.}
}

\examples{
########  data(cookie) ########
data(cookie)
myseq<-seq(141,651,by=2)
X<-as.matrix(cookie[-c(23,61),myseq])
Y<-as.matrix(cookie[-c(23,61),701:704])

fit.cca <- seedCCA(X[,1:4],Y[,1:4],type="cca")
fit.seed1 <- seedCCA(X,Y, type="seed1")
fit.pls1 <- seedCCA(X,Y[,1],type="pls")
fit.selu <- selectu(X,Y, type="seed2")

plot(fit.cca)
plot(fit.seed1, ref=95)
plot(fit.pls1)
plot(fit.pls1, eps=0.00001)
plot(fit.selu)
}
