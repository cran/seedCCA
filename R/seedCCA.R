#####################################################################
##  Seeded CCA
#####################################################################

seedCCA <- function(X,Y,type="seed2",ux=NULL,uy=NULL,u=10,eps=0.01,cut=0.9,d=NULL,AS=TRUE,scale=FALSE){
  X <- as.matrix(X); Y<-as.matrix(Y)
  n <- nrow(X)
  p<-dim(X)[2]
  r<-dim(Y)[2];
  M <- max(p,r)
  min.pr <- min(p,r)

  if(type=="cca"){
    if (max(p,r)>n) stop("The number of variables is bigger than sample sizes. Conduct Seeded CCA or PLS.")
    if (p==1 | r==1){
      cat("NOTE:  One of the two sets are 1-dimensional, so a linear regression via ordinary least square is fitted.", "\n")
      ans <- seedols(X,Y)}
    else{
      cat("NOTE:  The standard CCA is fitted for the two sets.", "\n")
      step2 <- finalCCA(X,Y)
      ans <- step2
      class(ans) <- c("finalCCA","seedCCA")
    }
  }

  if(type=="seed1"){
    cat("NOTE:  Seeded CCA with case 1 is fitted. The set with larger dimension is initially reduced.", "\n")
    cat("The first and second sets are denoted as X and Y, respectively. ", "\n")

    svd.n <- min(p,r)

    if(p>r & !is.null(ux)) { u.max <- ux
    }  else if(p<=r & !is.null(uy)) { u.max <- uy
    }  else {u.max <-u}

    sel.u <- selectu(X, Y, type="seed1", u.max=u.max, auto.stop=AS, num.d=svd.n, eps=eps)
    u <- sel.u$proper.u
    plot(sel.u)

    if(p>r){
      MX0<-iniCCA(X, Y, u=u, num.d=svd.n)
      new.X<-X%*%MX0
      step2<-finalCCA(new.X, Y)
      xcoef<-MX0%*%step2$xcoef
      ans<-list(cor=step2$cor, xcoef=xcoef, ycoef=step2$ycoef, proper.u=u, initialMX0=MX0, newX=new.X, Y=Y,
                Xscores=step2$Xphi, Yscores=step2$Yphi)
      class(ans) <- c("finalCCA","seedCCA")
    }else{
      MY0<-iniCCA(Y, X, u=u, num.d=svd.n)
      new.Y<-Y%*%MY0
      step2<-finalCCA(X,new.Y)
      ycoef<-MY0%*%step2$ycoef
      ans<-list(cor=step2$cor, xcoef=step2$xcoef, ycoef=ycoef, proper.u=u, X=X, initialMY0=MY0, newY=new.Y,
                Xscores=step2$Xphi, Yscores=step2$Yphi)
      class(ans) <- c("finalCCA","seedCCA")
    }
  }

  if(type=="seed2"){
    cat("NOTE:  Seeded CCA with case 2 is fitted. The two sets are initially reduced. ", "\n")
    cat("The first and second sets are denoted as X and Y, respectively. ", "\n")

    if (is.null(d)) {
      covxy <- cov(X,Y); evalues <- (svd(covxy)$d)^2
      percents <- cumsum(evalues) / sum(evalues)
      svd.n <- min(which(percents>cut))
    }else { svd.n <- d}

    if (is.null(ux) | is.null(uy)){
      sel.u <- selectu(X, Y, type="seed2", u.max=u, auto.stop=AS, num.d=svd.n, eps=eps)
      ux <- sel.u$proper.ux
      uy <- sel.u$proper.uy
      plot(sel.u)} else{
        ux<-ux ; uy<-uy }

    MX0<-iniCCA(X, Y, ux, svd.n)
    MY0<-iniCCA(Y, X, uy, svd.n)
    new.X<-X%*%MX0; new.Y<-Y%*%MY0
    step2<-finalCCA(new.X, new.Y)
    xcoef<-MX0%*%step2$xcoef
    ycoef<-MY0%*%step2$ycoef
    ans<-list(cor=step2$cor, xcoef=xcoef, ycoef=ycoef, proper.ux=ux, proper.uy=uy, d=svd.n,
              initialMX0=MX0, initialMY0=MY0, newX=new.X, newY=new.Y, Xscores=step2$Xscores, Yscores=step2$Yscores)
    class(ans) <- c("finalCCA","seedCCA")
  }

  if(type=="pls"){
    cat("Note: The first set (X) is predictors and the second set (Y) is response.", "\n")
    ans <- seedpls(X, Y, u=u, scale=scale)
  }

  return(ans)
}


#####################################################################
##  Finalized CCA
#####################################################################
finalCCA <- function(X, Y){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  cca.fit<-CCA::cc(X,Y)
  xcoef<-cca.fit$xcoef
  ycoef<-cca.fit$ycoef
  eigen<-cca.fit$cor
  MX<-X%*%xcoef; MY<-Y%*%ycoef
  ans<-list(cor=eigen, xcoef=xcoef, ycoef=ycoef, Xscores=MX, Yscores=MY )
  return(ans)
}

#####################################################################
##  Initialized CCA
#####################################################################
iniCCA<- function(X, Y, u, num.d){
  svd.n <- num.d
  dy <- dim(Y)[2] ; dx <- dim(X)[2]
  min.pr <- min(dx,dy)
  N<-nrow(X)

  if(svd.n==dy){
    s<- cov(X,Y)
    M<- seed <- s
  } else {
    s<-svd(cov(X,Y))
    M<- seed <- s$u[,(1:svd.n)]
  }

  Sx<-cov(X)
  if(u==1)
    B <- Pm(M, Sx, seed)
  if(u>1){
    M<-cov.seed<-seed
    for(i in 2:u){
      cov.seed<-Sx%*%cov.seed
      M<-cbind(M, cov.seed)
    }
    B <- Pm(M, Sx, seed)}
  return(B)
}


#####################################################################
##  partial least squares with seeded CCA
#####################################################################
seedpls <- function(X, Y, u=5, scale=FALSE){
  Y <- as.matrix(Y)
  ori.X <- X
  if (scale==FALSE) X <- as.matrix(X) else X<-scale(as.matrix(X))
  n <- dim(Y)[1]
  r <- dim(Y)[2]
  xnames <- colnames(X)
  ynames <- colnames(Y)
  if (is.null(ynames)) ynames <- paste("y", 1:r, sep="") else ynames <- ynames
  Bx <- list()
  for (i in 1:u) Bx[[i]] <- iniCCA(X, Y, u=i, num.d=r)
  names(Bx) <- paste("u=", 1:u, sep="")
  if (r==1){ for(i in 1:u) rownames(Bx[[i]]) <- xnames}
  else{for(i in 1:u) {rownames(Bx[[i]]) <- xnames; colnames(Bx[[i]])<-ynames} }

  ans<-c(list(coef=Bx), list(u=u, X=ori.X, Y=Y, scale=scale))
  class(ans) = c("seedpls","seedCCA")
  return(ans)
}

#####################################################################
#   ordinary least squares with seeded CCA
#####################################################################

seedols <- function(X, Y){

  p <- dim(as.matrix(X))[2]
  q <- dim(as.matrix(Y))[2]

  if (p==1)  coef=solve( t(scale(Y, scale = FALSE)) %*% scale(Y, scale = FALSE))%*%t(scale(Y, scale = FALSE))%*%scale(as.matrix(X),scale=FALSE)
  if (q==1)  coef=solve( t(scale(X, scale = FALSE)) %*% scale(X, scale = FALSE))%*%t(scale(X, scale = FALSE))%*%scale(as.matrix(Y),scale=FALSE)
  ans=list(coef=coef, X=X, Y=Y)
  class(ans) <- c("seedols", "seedCCA")
  return(ans)
}

########################################################################
### Function for selecting u (# of projections)
########################################################################
selectu<-function(X, Y, type="seed2", u.max=30, auto.stop=TRUE, num.d=2, eps=0.01){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  svd.n <- num.d
  n=dim(Y)[1]

  if(type=="seed2"){
    seedx <- svd(cov(X,Y))$u[,(1:svd.n)]
    seedy <- svd(cov(Y,X))$u[,(1:svd.n)]
    if (auto.stop){ Fu.x <- seeding.auto.stop(seedx, cov(X), n=n, u.max=u.max, eps=eps)
    Fu.y <- seeding.auto.stop(seedy, cov(Y), n=n, u.max=u.max, eps=eps)
    } else { Fu.x <- seeding(seedx, cov(X), n=n, u=u.max)
    Fu.y <- seeding(seedy, cov(Y), n=n, u=u.max) }
    nFu.x=Fu.x$nFu;  nFu.y=Fu.y$nFu

    if (length(which(nFu.x<eps))!=0) {ux <-min(which(nFu.x<eps))}
    else {ux <- u.max}
    if (length(which(nFu.y<eps))!=0) {uy <-min(which(nFu.y<eps))}
    else {uy <- u.max}

    ans <- list(nFu.x=Fu.x$nFu, nFu.y=Fu.y$nFu, proper.ux=ux, proper.uy=uy, type=type, eps=eps)

  }  else if (dim(X)[2] > dim(Y)[2]){
    if (auto.stop) {Fu <- seeding.auto.stop(cov(X,Y), cov(X), n=n, u.max=u.max, eps=eps)
    } else {Fu <- seeding(cov(X,Y), cov(X), n=n, u=u.max)}
    nFu=Fu$nFu
    if (length(which(nFu<eps))!=0) {u1 <-min(which(nFu<eps))
    } else {u1 <- u.max}
    ans <- list(nFu=Fu$nFu, proper.u=u1, type=type, eps=eps)

  } else { if (auto.stop) {Fu<-seeding.auto.stop(cov(Y,X), cov(Y), n=n, u.max=u.max, eps=eps)
  } else {Fu<-seeding(cov(Y,X), cov(Y), n=n, u=u.max)}
    nFu=Fu$nFu
    if (length(which(nFu<eps))!=0) {u2 <-min(which(nFu<eps))}
    else {u2 <- u.max}

    ans <- list(nFu=Fu$nFu, proper.u=u2, type=type, eps=eps)
  }
  class(ans) <- c("selectu", "seedCCA")
  return(ans)
}



########################################################################
### Seeding function
########################################################################
seeding<-function(seed, covx, n,  u=10){
  seed<-as.matrix(seed)
  Ru<-cov.seed<-seed
  m<-dim(Ru)[2];inner<-(t(Ru)%*% covx %*% Ru); rank.inner<-qr(inner)$rank
  if( m==1 ) {inner<-inner[1,1];B0<-(Ru %*% t(Ru) %*% seed)/inner}else{if(rank.inner == m){B0<-Ru %*% solve(inner, t(Ru)) %*% seed}
    else{ B0 <-Ru %*% corpcor::pseudoinverse(inner) %*% t(Ru)%*%seed}}

  Fu <- NULL

  for(i in 1:u){
    cov.seed<-covx%*%cov.seed; Ru<-cbind(Ru, cov.seed)
    m <- dim(Ru)[2]
    inner <- t(Ru) %*% covx %*% Ru
    rank.inner <- qr(inner)$rank
    if(rank.inner == m) {
      B1<-Ru %*% solve(inner, t(Ru)) %*% seed
      Delta<-B1-B0
      Fu[i]<- sum(diag(t(Delta)%*% covx %*% Delta))
      B <- B0
      B0 <- B1} else{
        B1<-Ru %*% corpcor::pseudoinverse(inner) %*% t(Ru) %*% seed
        Delta<- B1-B0
        Fu[i]<-sum(diag(t(Delta) %*% covx %*% Delta))
        B<-B0
        B0<-B1}
  }
  nFu<-Fu*n
  return(list(nFu=nFu))
}


seeding.auto.stop<-function(seed, covx, n, u.max=30, eps=0.01){
  seed<-as.matrix(seed)
  Ru<-cov.seed<-seed
  m<-dim(Ru)[2];inner<-(t(Ru)%*% covx %*% Ru); rank.inner<-qr(inner)$rank
  if( m==1 ) {inner<-inner[1,1];B0<-(Ru %*% t(Ru) %*% seed)/inner}else{if(rank.inner == m){B0<-Ru %*% solve(inner, t(Ru)) %*% seed}
    else{ B0 <-Ru %*% corpcor::pseudoinverse(inner) %*% t(Ru)%*%seed}}

  nFu <- NULL
  diff <- eps+0.1
  auto.u<-0
  i <- 1
  while(diff>eps) {
    cov.seed<-covx%*%cov.seed; Ru<-cbind(Ru, cov.seed)
    m <- dim(Ru)[2]
    inner <- t(Ru) %*% covx %*% Ru
    rank.inner <- qr(inner)$rank
    if(rank.inner == m)
    {B1<-Ru %*% solve(inner, t(Ru)) %*% seed
    Delta<-B1-B0
    diff <- n*sum(diag(t(Delta)%*% covx %*% Delta))
    nFu[i]<- diff; i<-i+1
    B <- B0
    B0 <- B1}else{
      B1<-Ru %*% corpcor::pseudoinverse(inner) %*% t(Ru) %*% seed
      Delta<- B1-B0
      diff <- n*sum(diag(t(Delta)%*% covx %*% Delta))
      nFu[i]<- diff; i<-i+1
      B<-B0
      B0<-B1}
    auto.u<-auto.u+1
    if (auto.u==u.max) {
      cat("The maximum number of iterations is reached.", "So, users must choose u bigger than",u.max,".","\n")
      break
    }
  }
  if (auto.u < u.max) {u<-auto.u}
  else {u<-u.max}

  ans=list(nFu = nFu, u=u)
  return(ans)
}
########################################################################
### coefficient ans fitted value functions for OLS and PLS
########################################################################
coef.seedCCA<- function(object, u=NULL,...){
  if ( class(object)[1]=="finalCCA" ) {
    stop("Not applicable for CCA, but applicable for OLS and PLS ")
  }

  if  ( class(object)[1]=="seedpls" ) {
    Bx <- object$coef
    umax <- object$u
    coef <- NULL

    if (is.null(u)) { {for(i in 1:umax) coef[[i]] <- round(Bx[[i]],3)}
      names(coef) <- paste("u=", 1:umax, sep="")
    }else {coef <- round(as.data.frame(t(Bx[[u]])),3)}
  }

  if  ( class(object)[1]=="seedols" ) {
    coef <- object$coef
  }
  return(coef)
}

fitted.seedCCA<- function(object, u=NULL,...){
  if ( class(object)[1]=="finalCCA" ) {
    stop("Not applicable for CCA, but applicable for OLS and PLS ")
  }

  if  ( class(object)[1]=="seedpls" ) {
    Bx <- object$coef
    y <- object$Y
    if (object$scale) {xc <- scale(object$X)} else {xc <- apply(object$X, 2, function(x){x-mean(x)})}

    fitted <- NULL

    if (is.null(u)) {u <- object$u
    for (i in 1:u) {
      newx <- xc %*%Bx[[i]]
      fitted[[i]] <- round(lm(y~newx)$fitted, 3)}
    names(fitted) <- paste("u=", 1:u, sep="")
    } else {u=u
    newx <- xc %*%Bx[[u]]
    fitted <- round(lm(y~newx)$fitted, 3)}
  }

  if  ( class(object)[1]=="seedols" ) {
    p <- dim(as.matrix(object$X))[2]
    q <- dim(as.matrix(object$Y))[2]
    coef <- object$coef
    r <- dim(coef)[1]

    if (p==r) fitted <- as.vector(scale(object$X, scale=FALSE)%*%coef)
    if (q==r) fitted <- as.vector(scale(object$Y, scale=FALSE)%*%coef)
  }
  return(fitted)
}

########################################################################
### Printing function
########################################################################

print.seedCCA <- function(x,...){
  if ( class(x)[1]=="finalCCA" ) {
    cor<-x$cor
    xcoef <- x$xcoef
    ycoef <- x$ycoef

    cat("Canonical correlation:\n")
    print(cor)

    cat("xcoef:\n")
    print(xcoef)

    cat("ycoef:\n")
    print(ycoef)

    cat("\n")
    invisible(x)
  }

  if  ( class(x)[1]=="seedpls" ) {
    cat("Coefs:\n")
    evectors<-x$coef
    u <- x$u; evs <- NULL
    for (i in 1:u) evs[[i]] <- round(evectors[[i]],4)
    names(evs) <- paste("u=", 1:u, sep="")
    print(evs)
    cat("\n")
    invisible(x)
  }

  if  ( class(x)[1]=="seedols" ) {
    cat("Coef:\n")
    evectors<-round(x$coef,4)
    print(evectors)
    cat("\n")
    invisible(x)
  }

  if  ( class(x)[1]=="selectu" ) {
    if(x$type=="seed2"){cat("$nFu.x\n")
      print(x$nFu.x)
      cat("\n")

      cat("$nFu.y\n")
      print(x$nFu.y)
      cat("\n")

      cat("$proper.ux\n")
      print(x$proper.ux)
      cat("\n")

      cat("$proper.uy\n")
      print(x$proper.uy)
      cat("\n")

      cat("$type\n")
      print(x$type)
      cat("\n")

      cat("$eps\n")
      print(x$eps)
      cat("\n")
    } else {cat("$nFu\n")
    print(x$nFu)
    cat("\n")

    cat("$proper.u\n")
    print(x$proper.u)
    cat("\n")

    cat("$type\n")
    print(x$type)
    cat("\n")

    cat("$eps\n")
    print(x$eps)
    cat("\n")
      }}
}

########################################################################
### Plotting function
########################################################################

plot.seedCCA <- function(x, ref=90, eps=0.01,...) {

  if ( class(x)[1]=="finalCCA" ) {
    ccor <- x$cor
    cum.percent <- cumsum(ccor)*100/sum(ccor)
    plot(cum.percent, type="b", ylim=c(0,100), lwd=2,
         ylab="Percent of cumulative canonical correlation", xlab="Number of canonical pairs")
    abline(h=ref, col=2, lwd=2)
    ans = NULL}

  if  ( class(x)[1]=="seedpls" ) {
    X <- x$X
    Y <- x$Y
    svd.n <- dim(Y)[2]
    n=dim(X)[1]
    u = x$u
    Fu <- seeding(cov(X,Y), cov(X), n=n, u=u)
    nFu=Fu$nFu
    if (length(which(nFu<eps))!=0) {u1 <-min(which(nFu<eps))}
    else {u1 <- u}

    cand.u <- u1+1
    X.name<-paste(1:cand.u, " ", "to"," ", 2:(cand.u+1), sep="")
    plot(nFu[1:cand.u], type="b", ylim=c(0, 20), xlab="u", ylab="n*Fu", axes=F)
    axis(1, at=1:cand.u, labels=X.name); axis(2); abline(h=eps,col=2, lwd=2)
    if ( length(which(nFu < eps))==0) {proper.u <- u
    cat("Caution: The terminating condition is NOT satisfied. The number of projections should be bigger than ",u,".","\n")
    } else {proper.u <- min(which(nFu < eps))}
    abline(v=proper.u, col=4, lwd=2)
    ans <- list(proper.u=proper.u, nFu=nFu, u=u, eps=eps)}

  if  ( class(x)[1]=="selectu" ) {
    if(x$type=="seed2"){
      ux <- x$proper.ux
      uy <- x$proper.uy

      cand.ux <- ux+1
      cand.uy <- uy+1

      X.name1<-paste(1:cand.ux, " ", "to"," ", 2:(cand.ux+1), sep="")
      X.name2<-paste(1:cand.uy, " ", "to"," ", 2:(cand.uy+1), sep="")
      eps <- x$eps

      par(mfrow=c(1,2))
      if(ux > 1){
        plot( x$nFu.x[1:cand.ux], type="b", ylim=c(0, 20), xlab="ux", ylab="nFu_x", axes=F)
        axis(1, at=1:cand.ux, labels=X.name1); axis(2); abline(h=eps,col=2, lwd=2); abline(v=ux, col=4, lwd=2)
      } else {plot(1:5, 1:5, type="n", xlab="ux", ylab="nFu_x")
        text(3, 3,"ux=1, so no plot is required.", cex=1)}

      if(uy > 1){
        plot( x$nFu.y[1:cand.uy], type="b", ylim=c(0, 20), xlab="uy", ylab="n*Fu_y", axes=F)
        axis(1, at=1:cand.uy, labels=X.name2); axis(2); abline(h=eps,col=2, lwd=2); abline(v=uy, col=4, lwd=2)
      } else {plot(1:5, 1:5, type="n", xlab="uy", ylab="nFu_y")
        text(3, 3,"uy=1, so no plot is required.", cex=1)}
    }
    else {
      if (x$proper.u==1) stop("The length of nFu is 1. There is no need for a plot")
      u <- x$proper.u
      cand.u <- u+1
      X.name<-paste(1:cand.u, " ", "to"," ", 2:(cand.u+1), sep="")
      eps <- x$eps
      plot( x$nFu[1:cand.u], type="b", ylim=c(0, 20), xlab="u", ylab="n*Fu", axes=F)
      axis(1, at=1:cand.u, labels=X.name); axis(2); abline(h=eps,col=2, lwd=2); abline(v=u, col=4, lwd=2)
    }
    ans=NULL
  }

  if  ( class(x)[1]=="seedols" ) stop("OLS, so no plot is constructed.")
  return(ans)
}

########################################################################
### Required function
########################################################################

Pm<-function(M, Sx, seed){
  M%*%corpcor::pseudoinverse(t(M)%*%Sx%*%M)%*%t(M)%*%seed
}

covplot<-function(X, Y, mind=NULL){
  temp<-cov(X,Y)
  min.pr <- min(dim(temp)[1], dim(temp)[2])
  if (is.null(mind)) mind <- min.pr

  s.value<-svd(cov(X,Y))$d^2
  ss.value<-cumsum(s.value)/sum(s.value)
  plot(s.value,main="Scree Plot of cov(X,Y)", ylab="Eigenvalues", xlab="", type="o", lwd=2)
  mtext("Number of Dimensions",side=1,line=4)

  ind <- rep(0, min.pr)
  ind[1:mind]<-1
  l<-1:length(ss.value)*ind

  upto60 <- min(which(ss.value>0.6))
  upto70 <- min(which(ss.value>0.7))
  upto80 <- min(which(ss.value>0.8))
  upto90 <- min(which(ss.value>0.9))

  text(l[ind!=0],s.value[l[ind!=0]],labels=paste(round(ss.value[l[ind!=0]],3)*100,"%"),pos=1,offset=0.2,cex=0.7,col=2)
  axis(1,at=c(1:length(ss.value[ind==1]),length(ss.value)),labels=c(round(ss.value[ind==1],3),1),line=2,cex=0.5)
  abline(v=upto60, col=2, lty=2)
  abline(v=upto70, col=3, lty=2)
  abline(v=upto80, col=4, lty=2)
  abline(v=upto90, col=5, lty=2)

  num.evecs <- c(upto60, upto70, upto80, upto90)
  names(num.evecs) <- c("up to 60%", "up to 70%", "up to 80%", "up to 90%")
  ans <- list(eigenvalue=s.value, cum.percent=ss.value, num.evecs=num.evecs)
  return(ans)
}
