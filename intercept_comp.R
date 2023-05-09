#estimates is a data frame from weighted_regression function. 
porownanie_wwolnych <- function(est1, est2, alpha=0.05, group.names = c("Male", "Female"), color = c("#F8766D", "#00BFC4")) {
  
  print(paste("-------- t-test -------"))
  tk=(est1$b0-est2$b0)/sqrt(est1$SEb0^2 + est2$SEb0^2)
  print(paste('p-value: ', 2*pt(-abs(tk),df=est1$n+est2$n-4)))
  print(paste('t-value: ', tk))
  
  df <- data.frame(group = group.names, b0 =c(est1$b0, est2$b0), CIdown = c(est1$CIdownb0, est2$CIdownb0),
                   CIup = c(est1$CIupb0, est2$CIupb0))
  require(ggplot2)
  theme_set(theme_bw(base_size = 20))
  p4 <- ggplot(data=df, aes(x=group, y=b0, col=group)) + 
    geom_errorbar(aes(ymin=CIdown, ymax=CIup), width=.2, lwd=2) +
    geom_point(color="gray33", size=2) +
    theme_bw(base_size = 14) +
    labs(title ="Intercept comparison", y="b0 coeff", x="group") +
    scale_color_manual(values=color) +
    theme(axis.title =element_text(size=18), 
          strip.text= element_text(size=18),
          axis.text = element_text(size=18)) 
  # ylim(c(-0.004, 0.0025))
  print(p4)
 
}