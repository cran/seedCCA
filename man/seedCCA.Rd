\name{seedCCA}
\alias{seedCCA}
\title{Seeded Canonical correlation analysis}
\description{
seedCCA is to conduct CCA when the maximum of the dimensions of the two sets of variables is greater than the sample sizes.}

\usage{
seedCCA(first.set, second.set, u1=2, u2=2, case1=FALSE, num.d=4)
}

\arguments{
\item{first.set}{numeric matrix (n * p), the first set of variables}
\item{second.set}{numeric matrix (n * r), the second set of variables}
\item{u1}{the termination index of projections for the first set of variables}
\item{u2}{the termination index of projections for the second set of variables}
\item{case1}{logical If FALSE initialized CCA are conducted for both variables.}
\item{num.d}{numeric, the number of the "num.d" largest eigenvectors of cov(first.set, second.set) and cov(second.set, first.set), if case1=FALSE. The default value is equal to 4. This option does not work, if case1=TRUE.}
}

\value{
\item{initialMX0}{the initialized canonical coefficient matrices of the first set of variables}
\item{initialMY0}{the initialized canonical coefficient matrices of the second set of variables}
\item{newX}{the initially-CCAed first set of variables)}
\item{newY}{the initially-CCAed second set of variables}
\item{xcoef}{the estimated canonical coefficients for the first set of variables}
\item{ycoef}{the estimated canonical coefficients for the second set of variables}
\item{Xcanvar}{the estimated canonical variates for the first set of variables}
\item{Ycanvar}{the estimated canonical variates for the second set of varialbes}
\item{eigenvalue}{the two sets of canonical correlations}
}
\references{
Y. Im, H. Gang and JK. Yoo (2014). High-throughput data dimension reduction via seeded canonical correlation analysis, J. Chemometrics 2015; 29: 193-199.

R. A. Johnson and D. W. Wichern(2007). Applied Multivariate Statistical Analysis. Pearson Prentice Hall: New Jersey, USA; 6 edition.539-574.

R. D. Cook, B. Li and F. Chiaromonte (2007). Dimension reduction in regression without matrix inversion. Biometrika 2007; 94: 569-584.

K. Lee and JK. Yoo(2014). Canonical correlation analysis through linearmodeling. Aust. Nz. J. Stat. 2014; 56: 59-72.
}
\examples{
######  data(cookie) ######
data(cookie)
myseq<-seq(141,651,by=2)
X<-as.matrix(cookie[-c(23,61),myseq])
Y<-as.matrix(cookie[-c(23,61),701:704])
dim(X);dim(Y)
selectu(X, Y, case1=TRUE)
seedCCA(X, Y, u1=2, case1=TRUE)

########  data(nutrimouse) ########
data(nutrimouse)
X<-as.matrix(nutrimouse$gene)
Y<-as.matrix(nutrimouse$lipid)
dim(X);dim(Y)
covplot(X, Y, mind=10)
selectu(X, Y, u=10, num.d=4)
seedCCA(X, Y, u1=6, u2=5, num.d=4)
}

