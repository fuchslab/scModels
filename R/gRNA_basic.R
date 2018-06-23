#' Gillespie Algorithm to simulate from basic model of mRNA generation and degradation
#' @param n number of simulations
#' @param r mRNA generation rate
#' @param mu mRNA degradation rate
#'
#' @export
#' @examples
#' gRNA_basic(1000,  r = 10, mu = 1)
#'
#'

gRNA_basic<-function(n,r,mu){
  t_0=0
  X_0= 0
  t_max=200
  Y<-rep(0,n)
  for (i in 1:n) {
    X<-X_0
    t_x<-t_0
    while(t_x<t_max){
      #i
      lambda_1<-r
      lambda_2<-mu*X
      lambdax=lambda_1+lambda_2

      #ii
      tau<-rexp(1,lambdax)
      tau_stern <- min(tau,t_max-t_x)
      #iii
      u<-runif(1)
      if (u<=(lambda_1/lambdax)) k=1 else  k=2

      #iv
      X=X+(k==1)-(k==2)
      #v
      t_x<-t_x+tau_stern
    }
    Y[i]<-X
  }
  return(Y)
}

