covplot<-function(first.set, second.set, mind=NULL){
  X <- first.set
  Y <- second.set
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
