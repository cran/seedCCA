\name{selectu}
\alias{selectu}
\title{Function that guides a selection of the terminating index when using seedCCA function
}
\description{
The usage of \code{selectu} depends on one of its arguments, \code{type}. If \code{tyep="seed1"}, the n*F_u is computed for a higher dimension one of X and Y and a proper number of prjections is reported. For example, suppose that the dimension of X is higher than Y. Then \code{selectu(X,Y, type="case1")} and \code{selectu(Y, X, u=5, type="case1")} gives the same results, and it is for X. If \code{type="seed2"}, n*F_u is computed for X and Y and proper numbers of projections for X and Y are reported. And, For \code{type="seed2"}, \code{num.d} must be specified. Its defualt value is 2. The argument \code{eps} is a terminating condition for stopping projections. The projection is stopped, when the increment is less than the value of \code{eps}. The argument \code{auto.stop=TRUE} has the function automatically stopped as soon as the increment is less than the value of \code{eps}. If not, the increments are computed until the value of \code{u.max} is reached. The function \code{selectu} results in subclass "selectu".
}

\usage{
selectu(X, Y, type="seed2", u.max=30, auto.stop=TRUE, num.d=2, eps=0.01)
}
\arguments{
  \item{X}{numeric matrix (n * p), the first set of variables}
  \item{Y}{numeric matrix (n * r), the second set of variables}
  \item{type}{character, the default is "seed2". "seed1" is for the first case of Seeded CCA (One set of variable is initially-reduced). "seed2" is for the second case of Seeded CCA (Two sets of variables are initially reduced).}
  \item{u.max}{numeric, t he maximum number of u. The default is equal to 30.}
  \item{auto.stop}{logical, The default value is TRUE. If TRUE, the iterative projection is automatically stopped, when the terminaion condition eps is satisfied. If FALSE, the iterative projections are stopped at the value of u.max.}
  \item{num.d}{numeric, the number of the "num.d" largest eigenvectors of cov(first.set, second.set), if case1=FALSE. The default value is equal to 2. This option works only for type="seed2".}
  \item{eps}{numeric, the default value is equal to 0.01. A value for terminating the projection.}

}


\details{
The order of the values depending on type is follows:

\code{type="seed1"}

\code{type="seed2"}
}

\value{
\item{type="seed1"}{Values with selecting \code{type="seed1"}:}
  \item{nFu}{incrments (n*Fu) of the iterative projection for initally reduction one set of variable.}
  \item{proper.u}{proper value of the number of projections for X}
  \item{type}{types of seeded CCA}
  \item{eps}{a value for terminating the projection.  The default value is equal to 0.01.}

\item{type="seed2"}{Values with selecting \code{type="seed2"}:}
  \item{nFu.x}{incrments (n*Fu) of the iterative projection for initially reducing X.}
  \item{nFu.y}{incrments (n*Fu) of the iterative projection for initially reducing Y.}
  \item{proper.ux}{proper value of the number of projections for X}
  \item{proper.uy}{proper value of the number of projections for Y}
  \item{type}{types of seeded CCA}
  \item{eps}{a value for terminating the projection.  The default value is equal to 0.01.}
}

\examples{
######  data(cookie) ######
data(cookie)
myseq<-seq(141,651,by=2)
X<-as.matrix(cookie[-c(23,61),myseq])
Y<-as.matrix(cookie[-c(23,61),701:704])

selectu(X, Y, type="seed1")
selectu(X, Y, type="seed1", auto.stop=FALSE)
selectu(X, Y, type="seed2", eps=0.001, num.d=3)
selectu(X, Y, type="seed2", auto.stop=FALSE)



########  data(nutrimouse) ########
data(nutrimouse)
Y<-as.matrix(nutrimouse$lipid)
X<-as.matrix(nutrimouse$gene)
selectu(X, Y, type="seed2", num.d=4)
selectu(X, Y, type="seed2", num.d=4, eps=0.001)
selectu(X, Y, type="seed2", auto.stop=FALSE, num.d=4, eps=0.001)
}

