#estimates is a data frame from weighted_regression function. 
porownanie_nachylen <- function(est1, est2, alpha=0.05, group.names = c("Male", "Female"), color = c("#F8766D", "#00BFC4")) {
    s1=est1$w_sdy*sqrt(((est1$n-1)/(est1$n-2))*(1-est1$r^2))
    s2=est2$w_sdy*sqrt(((est2$n-1)/(est2$n-2))*(1-est2$r^2))
    SE1=(s1^2)/((est1$w_sdx^2)*(est1$n-1))
    SE2=(s2^2)/((est2$w_sdx^2)*(est2$n-1))
    t=(est1$b1-est2$b1)/sqrt(SE1+SE2)
    t_critical=qt(1-alpha/2,est1$n+est2$n-4)
    print(paste('p-value: ', 2*pt(-abs(t),df=est1$n+est2$n-4)))
    
    print(paste('t-value: ', t))
    print(paste('t-critical: ', t_critical))
    
    #normal t test
    tk=(est1$b1-est2$b1)/sqrt(est1$SEb1^2 + est2$SEb1^2)
    print(paste('p-value: ', 2*pt(-abs(tk),df=est1$n+est2$n-4)))
    print(paste('t-value: ', tk))
    
    df <- data.frame(group = group.names, b0 =c(est1$b1, est2$b1), CIdown = c(est1$CIdownb1, est2$CIdownb1),
                     CIup = c(est1$CIupb1, est2$CIupb1))
    require(ggplot2)
    theme_set(theme_bw(base_size = 20))
    p4 <- ggplot(data=df, aes(x=group, y=b0, col=group)) + 
      geom_hline(aes(yintercept=0), col="red")+
      geom_errorbar(aes(ymin=CIdown, ymax=CIup), width=.2, lwd=2) +
      geom_point(color="gray33", size=2) +
      theme_bw(base_size = 14) +
      labs(title ="Slope comparison", y="b1 coeff", x="group") +
      scale_color_manual(values=color) +
      theme(axis.title =element_text(size=18), 
            strip.text= element_text(size=18),
            axis.text = element_text(size=18)) 
    # ylim(c(-0.004, 0.0025))
    print(p4)
}