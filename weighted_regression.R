## weighted least squares regression function ##
## Weights should sum up to 1;

weighted_regression <- function (x, y, w, alpha=0.05, print_plot=FALSE) {
    n=length(x)
    w=prop.table(w)
    # w=rep(1,n)
    y_weighted=sum(y*w)/sum(w)
    x_weighted=sum(x*w)/sum(w)
    b1=sum(w*(y-y_weighted)*(x-x_weighted))/sum(w*((x-x_weighted)^2))
    b0=y_weighted-b1*x_weighted
    print(b0)
    print(b1)
    estimates=data.frame(b0=b0, b1=b1)
    
    ##calculating confidence intervals for prediction
    predicted=b0+b1*x
    s=sqrt((1/(n-2))*sum((y-predicted)^2))
    kwadratx <- (x-x_weighted)^2
    SEy=s*sqrt((1/n) + kwadratx/sum(kwadratx))
    
    t_critical=qt(1-alpha/2,n-2)
    CIup=predicted+t_critical*SEy
    CIdown=predicted-t_critical*SEy
    df_pred=data.frame(x=x, y_observed=y, y_predicted=predicted, lower_CI=CIdown, upper_CI=CIup)
    

    if(print_plot==TRUE) {
      require(ggplot2)
      p2 <- ggplot(df_pred, aes(x, y_observed))+
        geom_point() +
        geom_line(aes(x,df_pred$y_predicted), col="darkred", size=1) +
        geom_ribbon(aes(ymin=lower_CI,ymax=upper_CI),alpha=0.4) +
        labs(x="Age", y="Diversity") 
      
      plot(p2)
    }
 
    
    SSR=sum((predicted-y_weighted)^2)
    SSE=sum((y-predicted)^2)
    SST=sum((y-y_weighted)^2)
    R2=SSR/SST
    F_value=SSR*(n-2)/SSE
    F_p_val=1-pf(F_value,1,n-2)
    F_critical=qf(1-alpha, df1=1, df2=n-2) 
    estimates$Fpval=F_p_val
    print(paste("R squared ", R2))
    print(paste("F test ", F_value, "F critical ", F_critical, "p_value ", F_p_val))
    
    ##calculating confidence intervals for b1
    residuals=y-predicted
    weighted_res=sqrt(w)*residuals
    SEres=sqrt(sum(weighted_res^2)/(n-2))
    
    SEb1=SEres/sqrt(sum(w*((x-x_weighted)^2)))
    CIupb1=b1+t_critical*SEb1
    CIdownb1=b1-t_critical*SEb1
    estimates$CIupb1=CIupb1
    estimates$CIdownb1=CIdownb1
    estimates$SEb1=SEb1
    
    ##additional info for comparison of b1 coefficients
    weighted_sdx=sqrt(sum(w*((x-x_weighted)^2))/sum(w))
    weighted_sdy=sqrt(sum(w*((y-y_weighted)^2))/sum(w))
    r=sum((x-x_weighted)*(y-y_weighted))/sqrt(sum((x-x_weighted)^2)*sum((y-y_weighted)^2))
    r1=sqrt(R2)
    print(paste("weighted correlation: ", r))
    estimates$w_sdx=weighted_sdx
    estimates$w_sdy=weighted_sdy
    estimates$r=r
    estimates$n=length(x)
    
    ## comparison of b0 coefficient
    # SEb0=SEres*sqrt(1/n + (x_weighted^2)/sum(w*((x-x_weighted)^2))) #zapo?yczone z wyk?adu Asi
    SEb0 = sqrt((SEres^2)/sum(w)+(SEb1^2)*(x_weighted^2))
    CIupb0=b0+t_critical*SEb0
    CIdownb0=b0-t_critical*SEb0
    estimates$CIupb0=CIupb0
    estimates$CIdownb0=CIdownb0
    estimates$SEb0=SEb0
    
    return(list(estimates, df_pred))
}