finalCCA <- function(first.ini.set, second.ini.set){
  X <- first.ini.set
  Y <- second.ini.set
  cca.fit<-CCA::cc(X,Y)
  xcoef<-cca.fit$xcoef
  ycoef<-cca.fit$ycoef
  eigen<-cca.fit$cor
  MX<-X%*%xcoef; MY<-Y%*%ycoef
  ans<-list("Xphi"=MX, "Yphi"=MY, xcoef=xcoef, ycoef=ycoef, eigen=eigen)
  return(ans)
}
