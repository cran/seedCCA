\name{selectu}
\alias{selectu}
\title{Function that guide the selection of the terminating index when using seedCCA function
}
\description{
When you use the seedCCA function, you need to specify the terminating indexes, for the first and second sets of variables. The function provides the guidance for  selecting the indexes for the first and second sets of variables. The function returns a plot of the increment in n*Fu
against u. The plot contains a red line, which incates that the increment is equal to 0.1. If "case1=TRUE", the n*F_u and screeplot for one to have a higher dimension between the first and second argument sets are reported. For example, suppose that the dimension of the first set is higher than the second set. Then selectu(X,Y, u=5, case1=TRUE) and selectu(Y, X, u=5, case1=TRUE) provide the equivalent scree plot for X,
but the returned name of n*F_u for selectu(X, Y, u=5, case1=TRUE) is F.u1, but that for for selectu(Y, X, u=5, case1=TRUE) is F.u2.
}

\usage{
selectu(first.set, second.set, u=5, case1=FALSE, num.d=4)
}
\arguments{
  \item{first.set}{numeric matrix (n * p), the first set of variables}
  \item{second.set}{numeric matrix (n * r), the second set of variables}
  \item{u}{numeric, the maximum number of u. The default value is equal to 5.}
  \item{case1}{logical, if FALSE, initialized CCA are conducted for both variables. The default value is FALSE. }
  \item{num.d}{numeric, the number of the "num.d" largest eigenvectors of cov(first.set, second.set), if case1=FALSE. The default value is equal to 4. This option does not work, if case1=TRUE.}
}
\value{
  \item{F.u_x}{n*Fu values of X}
  \item{F.u_y}{n*Fu values of Y}
}

\examples{
######  data(cookie) ######
data(cookie)
myseq<-seq(141,651,by=2)
X<-as.matrix(cookie[-c(23,61),myseq])
Y<-as.matrix(cookie[-c(23,61),701:704])
selectu(X, Y, case1=TRUE) ## nFu1 returned
selectu(Y, X, case1=TRUE) ## the same as selectu(X, Y, case1=TRUE), but nFu2 reuturned.

selectu(X, Y, u=10, case1=TRUE)


########  data(nutrimouse) ########
data(nutrimouse)
Y<-as.matrix(nutrimouse$lipid)
X<-as.matrix(nutrimouse$gene)
selectu(X, Y, u=10, num.d=4)
}

