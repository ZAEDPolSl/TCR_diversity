bic=function(model) {
  n=length(model$residuals)
  LL=logLik(model)
  bic_value=-2*LL+log(n)*length(model$coefficients)
  res=c(bic_value=as.numeric(bic_value), LL=as.numeric(LL))
  return(res)
}

source('weighted_regression.R')
source('slopes_comp.R')
source('intercept_comp.R')


#### data - cross-sectional
status <- read.table("statusDiv.txt", header = T, row.names = 1, sep="\t")
seq_prod <- read.table("prodSeqDiv.txt", header = T, row.names = 1, sep = "\t")
seq_non <- read.table("nonSeqDiv.txt", header = T, row.names = 1, sep = "\t")


df1 <- data.frame(status = status[,1], seq_prod =seq_prod[,1], seq_non = seq_non[,1])
row.names(df1)=rownames(status1)
descr <- read.table('samples_description_Dean.txt', sep="\t", stringsAsFactors = F, header = T)
df1$AGE=descr$AGE[which(descr$patient %in% sub("xx", '', row.names(df1)))]
df1$SEX=descr$SEX[which(descr$patient %in% sub("xx", '', row.names(df1)))]


gender <- rep(0,nrow(df1))
gender[which(df1$SEX=="Female")]=1
df1$gender=gender
df1$counts <- summary$prod.total+summary$non_total



## chose i for certain type of diversity _ status/ sequence.prod / sequence.non
i="seq_prod"

## comparison on general model in men and women
women<- which(df1$SEX=="Female")
men <- which(df1$SEX=="Male")
summary(lm(df1[,i] ~ AGE*gender, data=df1, weights = df1$counts))

summary(lm(df1[,i] ~ AGE, data=df1, weights = df1$counts))
model.all <- lm(df1[,i] ~ AGE, data=df1, weights = df1$counts)

summary(lm(df1[women,i] ~ df1[women,'AGE'], data=df1, weights = df1$counts[women]))
model.women <- lm(df1[women,i] ~ df1[women,'AGE'], data=df1, weights = df1$counts[women])

summary(lm(df1[men,i] ~ df1[men,'AGE'], data=df1, weights = df1$counts[men]))
model.men <- lm(df1[men,i] ~ df1[men,'AGE'], data=df1, weights = df1$counts[men])

# weighted regression - comparison of coefficients
s1 <-weighted_regression(df1$AGE, df1[,i], prop.table(df1$counts), print_plot = T)
s1w <-weighted_regression(df1$AGE[women],df1[women,i], prop.table(df1$counts[women]), print_plot = FALSE)
s1m <-weighted_regression(df1$AGE[men],df1[men,i], prop.table(df1$counts[men]), print_plot = FALSE)
slopes_comp(s1m[[1]], s1w[[1]])
intercept_comp(s1m[[1]], s1w[[1]])


##Segmented regression
# df1 <- df1[which(df1$gender==1),] ## only for women
# df1 <- df1[which(df1$gender==0),]  ## only for men
srt <- order(df1[,i], decreasing = F)
df1 <- df1[srt,]
srt <- order(df1$AGE)
df1 <- df1[srt,]

wyniki = matrix(NA, ncol=8, nrow=nrow(df1))

## brute-force technique
for (p in 5:(nrow(df1)-5)) {
  
  if(df1$AGE[p]==df1$AGE[p+1]) next;
  
  dyoung <- df1[1:p,] 
  dold <- df1[(p+1):nrow(df1),]
  
  ## comparison of models - youngs
  m1 <- BIC(lm(dyoung[,i] ~ AGE, data=dyoung, weights = dyoung$counts))
  m2 <- BIC(lm(dyoung[,i] ~ AGE + gender, data=dyoung, weights = dyoung$counts))
  m3 <- BIC(lm(dyoung[,i] ~ AGE + gender + AGE:gender, data=dyoung, weights = dyoung$counts))
  
  if(which.min(c(m1,m2,m3))!=1) print(paste("model young,", p))
  
  ## comparison of models - old
  n1 <- BIC(lm(dold[,i] ~ AGE, data=dold, weights = dold$counts))
  n2 <- BIC(lm(dold[,i] ~ AGE + gender, data=dold, weights = dold$counts))
  n3 <- BIC(lm(dold[,i] ~ AGE + gender + AGE:gender, data=dold, weights = dold$counts))
  
  if(which.min(c(n1,n2,n3))!=1) print(paste("model old,", p))
  
  wyniki[p,1:6] = c(m1,m2,m3,n1,n2,n3)
  wyniki[p,7] = mean(c(min(m1,m2,m3), min(n1,n2,n3)))
  wyniki[p,8] = min(m1,m2,m3)+min(n1,n2,n3)
  
}

wyniki = cbind(wyniki, df1$AGE, df1$SEX, df1$gender)

# ## BIC - vs the number of donors in a group
# dfplot <- data.frame(group=rep(c("young", "old"), each=478), 
#                      setsize=rep(c(5:482),2),
#                     BIC = c(as.numeric(wyniki[5:482,1]), as.numeric(wyniki[5:482,4])), 
#                     BICn = c(as.numeric(wyniki[5:482,1])/c(5:482), as.numeric(wyniki[5:482,4])/c(482:5)))
# ggplot(dfplot, aes(setsize, BICn, col=group, group=group)) +
#   geom_point()


head(wyniki[order(wyniki[,8], decreasing = T),])
kt <- which.min(wyniki[,8])
#choice for women and status diversity - as BIC value are similar (difference less than 2)
# kt <- 6
kty <- which.min(wyniki[kt,1:3])
kto <- which.min(wyniki[kt,4:6])

dyoung <- df1[1:kt,] 
dold <- df1[(kt+1):nrow(df1),]


modely = lm(dyoung[,i] ~ AGE, data=dyoung, weights = dyoung$counts)
modelo = lm(dold[,i] ~ AGE, data=dold, weights = dold$counts)
# modelo = lm(dold[,i] ~ AGE+gender, data=dold, weights = dold$counts)
# modelo = lm(dold[,i] ~ AGE+gender+AGE:gender, data=dold, weights = dold$counts)
summary(modely)
m_young <- weighted_regression(x = dyoung$AGE, y=dyoung[,i], w = dyoung$counts)
summary(modelo)
m_old <- weighted_regression(x = dold$AGE, y=dold[,i], w = dold$counts)

clopes_comp(m_young[[1]], m_old[[1]],group.names = c("1.Young", "2.Old"), color = c("#00BFC4", "#00BFC4"))

