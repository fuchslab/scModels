## ---- gRNA_burst_Chr
gmRNA_burst<-function(n,r_burst, s_burst, r_deg){
    t_0=0
    X_0=0
    t_max= 200
    Y<-rep(0,n)
    for (i in 1:n) {
        X<-X_0
        t_x<-t_0
        while(t_x<t_max){
            #i
            lambda_1<-r_burst
            lambda_4<-r_deg*X
            lambdax=lambda_1+lambda_4

            #ii
            tau<-rexp(1,lambdax)
            tau_stern <- min(tau,t_max-t_x)
            #iii
            u<-runif(1)
            if (u<=(lambda_1/lambdax)) k=1 else k=4
            #iv
            if (tau<=tau_stern) {
                if (k==1) {
                    # Burst
                    X <- X + rgeom(1,1/(1+s_burst))
                }
                else if (k==4) {
                    # Degradation
                    X <- X-1
                }
            }
            #v
            t_x<-t_x+tau_stern
        }
        Y[i]<-X
    }
    return(Y)
}
