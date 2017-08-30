seedCCA<-function(first.set, second.set, u1=2, u2=2, case1=FALSE, num.d=4){

  X <- first.set
  Y <- second.set
  svd.n <- num.d
  ux <- u1
  uy <- u2

  X <- as.matrix(X); Y<-as.matrix(Y)
  n <- nrow(X)
  p<-dim(X)[2]
  r<-dim(Y)[2];
  M <- max(p,r)
  min.pr <- min(p,r)

  if(M<n){
    step2 <- finalCCA(X,Y)
    ans <- step2
  }

  else if(M>=n){

    if(case1==FALSE){
      MX0<-iniCCA(X,Y,ux,svd.n)
      MY0<-iniCCA(Y,X,uy,svd.n)
      new.X<-X%*%MX0; new.Y<-Y%*%MY0
      step2<-finalCCA(new.X,new.Y)
      xcoef<-MX0%*%step2$xcoef
      ycoef<-MY0%*%step2$ycoef
      ans<-list(initialMX0=MX0, initialMY0=MY0, newX=new.X, newY=new.Y,
	xcoef=xcoef, ycoef=ycoef, Xcanvar=step2$Xphi, Ycanvar=step2$Yphi, eigenvalue=step2$eigen)
    }

    else if(case1==TRUE){
      svd.n <- min(p,r)
      if(p>r){
        MX0<-iniCCA(X, Y, ux, svd.n)
        new.X<-X%*%MX0
        step2<-finalCCA(new.X, Y)
        xcoef<-MX0%*%step2$xcoef
	ans<-list(initialMX0=MX0, newX=new.X, Y=Y,
	  xcoef=xcoef, ycoef=step2$ycoef, Xcanvar=step2$Xphi, Ycanvar=step2$Yphi, eigenvalue=step2$eigen)
      }else{
	MY0<-iniCCA(Y, X, uy, svd.n)
        new.Y<-Y%*%MY0
        step2<-finalCCA(X,new.Y)
        ycoef<-MY0%*%step2$ycoef

        ans<-list(X=X, initialMY0=MY0, newY=new.Y,
	  xcoef=step2$xcoef, ycoef=ycoef, Xcanvar=step2$Xphi, Ycanvar=step2$Yphi, eigenvalue=step2$eigen)
      }
    }
  }
  return(ans)
}
