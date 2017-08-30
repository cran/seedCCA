selectu<-function(first.set, second.set, u=5, case1=FALSE, num.d=4){
  X <- first.set
  Y <- second.set
  X <- as.matrix(X); Y<-as.matrix(Y)
  svd.n <- num.d
  n=dim(X)[1]
  X.name<-paste(1:u, " ", "to"," ", 2:(u+1), sep="")

    if(case1==FALSE){
       seedx <- svd(cov(X,Y))$u[,(1:svd.n)]
       seedy <- svd(cov(Y,X))$u[,(1:svd.n)]
       F.u1 <- seeding(seedx, cov(X), u=u, n=dim(X)[1])
       F.u2 <- seeding(seedy, cov(Y), u=u,  n=dim(X)[1])

       par(mfrow=c(1,2))
        ## if (max(F.u1$F.u)>20){ y.max<-ceiling(max(F.u1$F.u))} else{ y.max<-20}
          plot(F.u1$F.u, type="b", ylim=c(0, 20), xlab="u1", ylab="n*Fu1", axes=F)
          axis(1, at=1:u, labels=X.name); axis(2); abline(h=0.1,col=2, lwd=2)
        ## if (max(F.u2$F.u)>20) {y.max<-ceiling(max(F.u2$F.u))} else {y.max<-20}
          plot(F.u2$F.u, type="b", ylim=c(0, 20), xlab="u2", ylab="n*Fu2", axes=F)
          axis(1, at=1:u, labels=X.name); axis(2); abline(h=0.1,col=2)
  return(list(F.u1=F.u1$F.u, F.u2=F.u2$F.u))
  }else if(dim(X)[2]>dim(Y)[2])
  {
    F.u1 <- seeding(cov(X,Y), cov(X), n=dim(X)[1], u=u)
    ## if (max(F.u1$F.u)>20) {y.max<-ceiling(max(F.u1$F.u))} else{ y.max<-20}
    plot(F.u1$F.u, type="b", ylim=c(0, 20), xlab="u1", ylab="n*Fu1", axes=F)
    axis(1, at=1:u, labels=X.name); axis(2); abline(h=0.1,col=2)
    return(list(F.u1=F.u1$F.u))
  }else
  {
    F.u2<-seeding(cov(Y,X),cov(Y), n=dim(X)[1], u=u)
    ## if (max(F.u2$F.u)>20) {y.max<-ceiling(max(F.u2$F.u))} else {y.max<-20}
    plot(F.u2$F.u, type="b", ylim=c(0, 20), xlab="u2", ylab="n*Fu2", axes=F)
    axis(1, at=1:u, labels=X.name); axis(2); abline(h=0.1,col=2)
    return(list(F.u2=F.u2$F.u))
  }
}
