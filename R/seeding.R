seeding<-function(seed, cov.x, n, u=5){
  R.u<-cov.seed<-seed
  m<-dim(R.u)[2];inner<-(t(R.u)%*% cov.x %*% R.u); rank.inner<-qr(inner)$rank
  if( m==1 ) {inner<-inner[1,1];B0<-(R.u %*% t(R.u) %*% seed)/inner}else{if(rank.inner == m){B0<-R.u %*% solve(inner, t(R.u)) %*% seed}
    else{ B0 <-R.u %*% pseudoinverse(inner) %*% t(R.u)%*%seed}}
  F.u <- rep(0,u)
  for(i in 1:u){
    cov.seed<-cov.x%*%cov.seed; R.u<-cbind(R.u, cov.seed)
    m <- dim(R.u)[2]
    inner <- t(R.u) %*% cov.x %*% R.u
    rank.inner <- qr(inner)$rank
    if(rank.inner == m)
    {B1<-R.u %*% solve(inner, t(R.u)) %*% seed
    Delta<-B1-B0
    F.u[i]<- sum(diag(t(Delta)%*% cov.x %*% Delta))
    B <- B0
    B0 <- B1}else{
      B1<-R.u %*% pseudoinverse(inner) %*% t(R.u) %*% seed
      Delta<- B1-B0
      F.u[i]<-sum(diag(t(Delta) %*% cov.x %*% Delta))
      B<-B0
      B0<-B1}
  }
  F.u<-F.u*n
  return(list(F.u=F.u))
}
