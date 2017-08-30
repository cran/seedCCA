\name{covplot}
\alias{covplot}
\title{scree-ploting cov(first.set, second.set)}
\description{Returns a scree-plot of the eigenvalues of cov(first.set, second.set) to select its first d largest eigenvectors.}

\usage{covplot(first.set, second.set, mind=NULL)}

\arguments{
\item{first.set}{numeric matrix (n * p), the first set of variables}
\item{second.set}{numeric matrix (n * r), the second set of variables}
\item{mind}{numeric, the number of the eigenvalues to show their cumulative percentages. The default is NULL, and then it is equal to min(p,r)}
}

\value{
\item{eigenvalues}{the ordiered eigenvalues of cov(X,Y) }
\item{cum.percent}{the cumulative percentages of the eigenvalues}
\item{num.evecs}{a vector of the numbers of the eigenvectors which forces the cumulative percentages bigger than 0.6, 0.7, 0.8, 0.9}
}

\examples{
data(cookie)
myseq<-seq(141,651,by=2)
X<-as.matrix(cookie[-c(23,61),myseq])
Y<-as.matrix(cookie[-c(23,61),701:704])
covplot(X, Y)
covplot(X, Y, mind=4)
}
