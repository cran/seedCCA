Pm<-function(M, Sx, seed){
  M%*%corpcor::pseudoinverse(t(M)%*%Sx%*%M)%*%t(M)%*%seed
}


