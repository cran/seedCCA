\name{seedCCA}
\alias{seedCCA}
\title{Seeded Canonical correlation analysis}
\description{
The function \code{seedCCA} is mainly for implementing seeded canonical correlation analysis proposed by Im et al. (2015). The function conducts the following four methods, depending on the value of \code{type}. The option \code{type} has one of \code{c("cca", "seed1", "seed2", "pls")}. }

\usage{
seedCCA(X,Y,type="seed2",ux=NULL,uy=NULL,u=10,eps=0.01,cut=0.9,d=NULL,AS=TRUE,scale=FALSE)
}

\arguments{
\item{X}{numeric vector or matrix (n * p), the first set of variables}
\item{Y}{numeric vector or matrix (n * r), the second set of variables}
\item{type}{character, a choice of methods among \code{c("cca", "seed1", "seed2", "pls")}. The default is \code{"seed2"}.}
\item{ux}{numeric, maximum number of projections for X. The default is NULL. If this is not NULL, it surpasses the option \code{u} with \code{type="seed1"} and p>r. For \code{type="seed2"}, if \code{ux} and \code{uy} are not NULL, they surpass \code{u}.}
\item{uy}{numeric, maximum number of projections for Y. The default is NULL. If this is not NULL, it surpasses the option \code{u} with \code{type="seed1"} and p<r. For \code{type="seed2"}, if \code{ux} and \code{uy} are not NULL, they surpass \code{u}.}
\item{u}{numeric, maximum number of projections. The default is 10. This is used for \code{type="seed1"}, \code{type="seed2"} and {tyepe="pls"}.}
\item{eps}{numeric, the criteria to terminate iterative projections. The default is 0.01. If increment of projections is less than \code{eps}, then the iterative projection is terminated.}
\item{cut}{numeric, between 0 and 1. The default is 0.9.
If \code{d} is NULL, \code{cut} is used for automatic replacements of cov(X,Y) and cov(Y,X) with their eigenvectors, depending on the value of \code{cut}. So, if any value of \code{d} is given, \code{cut} is not effective. cov(X,Y) and cov(Y,X) are replaced with their largest eigenvectors, whose cumulative eigenvalue proportion is bigger than the value of \code{cut}. This only works for \code{type="seed2"}.}
\item{d}{numeric, the user-selected number of largest eigenvectors of cov(X, Y) and cov(Y, X). The default is NULL. This only works for \code{type="seed2"}. If any value of \code{d} is given, \code{cut} does not work.}
\item{AS}{logical, status of automatic stop of projections. The default is \code{TRUE}. If TRUE, the iterative projection is automatically stopped, when the terminaion condition \code{eps} is satisfied. If\code{AS=FALSE}, the iterative projections are stopped at the value of \code{u}. }
\item{scale}{logical. scaling predictors to have zero mean and one standard deviation. The default is \code{FALSE}. If \code{scale=TRUE}, each predictor is scaled with mean 0 and variance 1 for partial least squares. This option works only for \code{type="pls"}.}
}

\details{
Let p and r stand for the numbers of variables in the two sets and n stands for the sample size. The option of \code{type="cca"} can work only when max(p,r) < n, and \code{seedCCA} conducts standard canonical correlation analysis (Johnson and Wichern, 2007). If \code{type="cca"} is given and either p or r is equal to one, ordinary least squares (OLS) is done instead of canonical correlation analysis. If max(p,r) >= n, either \code{type="seed1"} or \code{type="seed2"} has to be chosen. This is the main purpose of \code{seedCCA}. If \code{type="seed1"}, only one set of variables, saying X with p for convenience, to have more variables than the other, saying Y with r, is initially reduced by the iterative projection approach (Cook et al. 2007). And then, the canonical correlation analysis of the initially-reduced X and the original Y is finalized. If \code{type="seed2"}, both X and Y are initially reduced. And then, the canonical correlation analysis of the two initially-reduced X and Y are finalzed. If \code{type="pls"}, partial least squares (PLS) is done. If \code{type="pls"} is given, the first set of variables in \code{seedCCA} is predictors and the second set is response. This matters The response can be multivariate. Depeding on the value of \code{type}, the resulted subclass by \code{seedCCA} are different.:

\code{type="cca"}: subclass "finalCCA" (p >2; r >2; p,r<n)

\code{type="cca"}: subclass "seedols"  (either p or r is equal to 1.)

\code{type="seed1"} and \code{type="seed2"}: subclass "finalCCA" (max(p,r)>n)

\code{type="pls"}: subclass "seedpls" (p>n and r <n)

So, \code{plot(object)} will result in different figure depending on the object.

The order of the values depending on type is follows.:

\code{type="cca"}: standard CCA (max(p,r)<n, min(p,r)>1) / "finalCCA" subclass

\code{type="cca"}: ordinary least squares (max(p,r)<n, min(p,r)=1) / "seedols" subclass

\code{type="seed1"}: seeded CCA with case1 (max(p,r)>n and p>r) / "finalCCA" subclass

\code{type="seed1"}: seeded CCA with case1 (max(p,r)>n and p<=r) / "finalCCA" subclass

\code{type="seed2"}: seeded CCA with case2 (max(p,r)>n) / "finalCCA" subclass

\code{type="pls"}: partial least squares (p>n and r<n) / "seedpls" subclass
}

\value{
\item{type="cca"}{Values with selecting \code{type="cca"}: standard CCA (max(p,r)<n, min(p,r)>1) / "finalCCA" subclass}
\item{cor}{canonical correlations}
\item{xcoef}{the estimated canonical coefficients for X}
\item{ycoef}{the estimated canonical coefficients for Y}
\item{Xscores}{the estimated canonical variates for X}
\item{Yscores}{the estimated canonical variates for Y}

\item{type="cca"}{Values with selecting \code{type="cca"}: ordinary least squares (max(p,r)<n, min(p,r)=1) / "seedols" subclass}
\item{coef}{the estimated ordinary least squares coefficients}
\item{X}{X, the first set}
\item{Y}{Y, the second set}

\item{type="seed1"}{Values with selecting \code{type="seed1"}: seeded CCA with case1 (max(p,r)>n and p>r) / "finalCCA" subclass}
\item{cor}{canonical correlations}
\item{xcoef}{the estimated canonical coefficients for X}
\item{ycoef}{the estimated canonical coefficients for Y}
\item{proper.u}{a suggested proper number of projections for X}
\item{initialMX0}{the initialized canonical coefficient matrices of X}
\item{newX}{initially-reduced X}
\item{Y}{the original Y}
\item{Xscores}{the estimated canonical variates for X}
\item{Yscores}{the estimated canonical variates for Y}


\item{type="seed1"}{Values with selecting \code{type="seed1"}: seeded CCA with case1 (max(p,r)>n and p<=r) / "finalCCA" subclass)}
\item{cor}{canonical correlations}
\item{xcoef}{the estimated canonical coefficients for X}
\item{ycoef}{the estimated canonical coefficients for Y}
\item{proper.u}{a suggested proper number of projections for Y}
\item{X}{the original X}
\item{initialMY0}{the initialized canonical coefficient matrices of Y}
\item{newY}{initially-reduced Y}
\item{Xscores}{the estimated canonical variates for X}
\item{Yscores}{the estimated canonical variates for Y}


\item{type="seed2"}{Values with selecting \code{type="seed2"}: seeded CCA with case2 (max(p,r)>n) / "finalCCA" subclass}
\item{cor}{canonical correlations}
\item{xcoef}{the estimated canonical coefficients for X}
\item{ycoef}{the estimated canonical coefficients for Y}
\item{proper.ux}{a suggested proper number of projections for X}
\item{proper.uy}{a suggested proper number of projections for Y}
\item{d}{suggested number of eigenvectors of cov(X,Y)}
\item{initialMX0}{the initialized canonical coefficient matrices of X}
\item{initialMY0}{the initialized canonical coefficient matrices of Y}
\item{newX}{initially-reduced X}
\item{newY}{initially-reduced Y}
\item{Xscores}{the estimated canonical variates for X}
\item{Yscores}{the estimated canonical variates for Y}


\item{type="pls"}{Values with selecting \code{type="pls"}:: partial least squares (p>n and r<n) / "seedpls" subclass}
\item{coef}{the estimated coefficients for each iterative projection upto u}
\item{u}{the maximum number of projections}
\item{X}{predictors}
\item{Y}{response}
\item{scale}{status of scaling predictors}
\item{cases}{the number of observations}
}
\references{
R. D. Cook, B. Li and F. Chiaromonte. Dimension reduction in regression without matrix inversion. Biometrika 2007; 94: 569-584.

Y. Im, H. Gang and JK. Yoo. High-throughput data dimension reduction via seeded canonical correlation analysis, J. Chemometrics 2015; 29: 193-199.

R. A. Johnson and D. W. Wichern. Applied Multivariate Statistical Analysis. Pearson Prentice Hall: New Jersey, USA; 6 edition. 2007; 539-574.

K. Lee and JK. Yoo. Canonical correlation analysis through linear modeling, AUST. NZ. J. STAT. 2014;  56: 59-72.

}
\examples{
######  data(cookie) ######
data(cookie)
myseq<-seq(141,651,by=2)
X<-as.matrix(cookie[-c(23,61),myseq])
Y<-as.matrix(cookie[-c(23,61),701:704])
dim(X);dim(Y)

## standard CCA
fit.cca <-seedCCA(X[,1:4], Y, type="cca")  ## standard canonical correlation analysis is done.
plot(fit.cca)

## ordinary least squares
fit.ols1 <-seedCCA(X[,1:4], Y[,1], type="cca")  ## ordinary least squares is done, because r=1.
fit.ols2 <-seedCCA(Y[,1], X[,1:4], type="cca")  ## ordinary least squares is done, because p=1.

## seeded CCA with case 1
fit.seed1 <- seedCCA(X, Y, type="seed1") ## suggested proper value of u is equal to 3.
fit.seed1.ux <- seedCCA(X, Y, ux=6, type="seed1") ## iterative projections done 6 times.
fit.seed1.uy <- seedCCA(Y, X, uy=6, type="seed1", AS=FALSE)  ## projections not done until uy=6.
plot(fit.seed1)

## partial least squares
fit.pls1 <- seedCCA(X, Y[,1], type="pls")
fit.pls.m <- seedCCA(X, Y, type="pls") ## multi-dimensional response
par(mfrow=c(1,2))
plot(fit.pls1); plot(fit.pls.m)


########  data(nutrimouse) ########
data(nutrimouse)
X<-as.matrix(nutrimouse$gene)
Y<-as.matrix(nutrimouse$lipid)
dim(X);dim(Y)

## seeded CCA with case 2
fit.seed2 <- seedCCA(X, Y, type="seed2")  ## d not specified, so cut=0.9 (default) used.
fit.seed2.99 <- seedCCA(X, Y, type="seed2", cut=0.99)  ## cut=0.99 used.
fit.seed2.d3 <- seedCCA(X, Y, type="seed2", d=3)  ## d is specified with 3.

## ux and uy specified, so proper values not suggested.
fit.seed2.uxuy <- seedCCA(X, Y, type="seed2", ux=10, uy=10)
plot(fit.seed2)
}

