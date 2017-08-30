iniCCA<- function(first.set, second.set, u, num.d){
  X <- first.set
  Y <- second.set
  svd.n <- num.d
  dy <- dim(Y)[2] ; dx <- dim(X[2])
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
