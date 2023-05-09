get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

##Set 2 - logitudinal - already with information about age and sex
df2 <- read.table('CD4CD8_downsampled.txt', header = T, row.names = 1, sep="\t", dec=",")
df2$gender=df2$gender-1
colnames(df2) <- c("status", "seq_prod", "seq_non", "V_prod", "V_non", "J_prod", "J_non", "counts", "AGE", "SEX", "patnum")

## chose i for certain type of diversity
i="seq_prod"

library(lme4)

#random intercept model for individual donor differences
#all slopes correlated and sex does not impact..?
fit.0<-lmer(df2[,i]~1+(1|patnum), data=df2)


#1|patnum means that every patient (patnum) has got its own intercept
#using offset can set a known a priori mean (intercept value?) 

fit.1<-lmer(df2[,i]~AGE+(1|patnum), data=df2)
summary(fit.1)
anova(fit.1)
fixef(fit.1)
ranef(fit.1)
extractAIC(fit.1)
VarCorr(fit.1)

# fit.1a<-lmer(status~SEX+(1|patnum), data=df2)
# summary(fit.1)
# anova(fit.1)
# fixef(fit.1)
# ranef(fit.1)
# extractAIC(fit.1)

#all slopes correlated with each other.
fit.2<-lmer(df2[,i]~AGE+SEX+(1|patnum), data=df2)
summary(fit.2)
anova(fit.2)
fixef(fit.2)
ranef(fit.2)
extractAIC(fit.2)
VarCorr(fit.2)

## slopes for men correlated and slopes for women correlated
fit.3<-lmer(df2[,i]~AGE+SEX+AGE:SEX+(1|patnum), data=df2)
summary(fit.3)
anova(fit.3)
fixef(fit.3)
ranef(fit.3)
extractAIC(fit.3)
VarCorr(fit.3)

# fit.3a<-lmer(status~AGE:SEX+(1|patnum), data=df2)
# fit.3b<-lmer(status~AGE + AGE:SEX+(1|patnum), data=df2)

anova(fit.0, fit.1, fit.2, fit.3)

# #individual trajectory model with random slope for time
# #i.e. we allow that men and women have different trajectories over time
# fit.4<-lmer(status~AGE+SEX+(AGE|patnum), data=df2)
# summary(fit.4)
# anova(fit.4)
# fixef(fit.4)
# ranef(fit.4)
# 
# fit.4a<-lmer(status~AGE+SEX+AGE:SEX+(AGE|patnum), data=df2)
# anova(fit.1, fit.4)


df2$gender <- as.factor(ifelse(df2$SEX=='0', 'Men', 'Female'))
df2$fitted <- fitted(fit.2)
df2$patient <- as.factor(df2$patnum)



##plot
library(ggplot2)
png("VGene_MixedEffects_model1.png", res=200, height=650, width=1200)
p <-  ggplot(df2, aes(AGE, df2[,i], group = patient, colour=gender))+
  geom_point(aes(shape = patient), size =3) +
  geom_line(aes(AGE, fitted)) +
  # facet_wrap(~cell)+
  labs(x="Age", y="V Gene Diversity", title="Diversity ~ AGE + donors") +
  theme_bw()+
  theme(axis.title =element_text(size=22),
        strip.text= element_text(size=22),
        axis.text = element_text(size=22))
plot(p)
dev.off()



##############################-------------------------------------------------------------------
## cd4 and cd8 separately

df2 <- read.table('../CD4_CD8_sep/CD4_CD8_sep.txt', header = T, row.names = 1, sep="\t", dec=",")
df2$cell <- ifelse(df2$CD4_CD8=="CD4", 0, 1)
colnames(df2) <- c("status", "seq_prod", "seq_non", "V_prod", "V_non", "J_prod", "J_non", "counts", "AGE", "SEX", "patnum", "CD4_CD8", "cell")

## chose i for certain type of diversity
i="seq_prod"

library(lme4)
#random intercept model for individual donor differences
#all slopes correlated and sex does not impact..?
fit.0<-lmer(df2[,i]~1+(1|patnum), data=df2)
fit.0a<-lmer(df2[,i]~1+(cell|patnum), data=df2)
anova(fit.0, fit.0a)

#1|patnum means that every patient (patnum) has got its own intercept

fit.1<-lmer(df2[,i]~AGE+(cell|patnum), data=df2)
summary(fit.1)
anova(fit.1)
fixef(fit.1)
ranef(fit.1)
extractAIC(fit.1)
VarCorr(fit.1)


#all slopes correlated with each other.
fit.2<-lmer(df2[,i]~AGE+SEX+(cell|patnum), data=df2)
summary(fit.2)
anova(fit.2)
fixef(fit.2)
ranef(fit.2)
extractAIC(fit.2)
VarCorr(fit.2)

## slopes for men correlated and slopes for women correlated
fit.3<-lmer(df2[,i]~AGE+SEX+AGE:SEX+(cell|patnum), data=df2)
summary(fit.3)
anova(fit.3)
fixef(fit.3)
ranef(fit.3)
extractAIC(fit.3)
VarCorr(fit.3)

fit.4<-lmer(df2[,i]~AGE + SEX+AGE:SEX + cell + (cell|patnum), data=df2)
fit.5<-lmer(df2[,i]~AGE + SEX+AGE:SEX + AGE:cell + (cell|patnum), data=df2)
fit.6<-lmer(df2[,i]~AGE + SEX+AGE:SEX + cell + AGE:cell + (cell|patnum), data=df2)
fit.7<-lmer(df2[,i]~AGE + SEX+AGE:SEX + SEX:cell + (cell|patnum), data=df2)
fit.8<-lmer(df2[,i]~AGE + SEX+AGE:SEX + cell + SEX:cell + (cell|patnum), data=df2)

anova(fit.0a, fit.1, fit.2, fit.3, fit.4, fit.5, fit.6, fit.7, fit.8)

fit.5<-lmer(df2[,i]~AGE + SEX+AGE:SEX + AGE:cell + (cell|patnum), data=df2)
fit.5a<-lmer(df2[,i]~AGE + SEX+AGE:SEX + AGE:cell + (1|patnum), data=df2)
fit.5b<-lmer(df2[,i]~AGE + SEX+AGE:SEX + AGE:cell + (1|patnum:cell), data=df2)
fit.6 <- lmer(df2[,i]~AGE + SEX+AGE:SEX + cell + AGE:cell + (cell|patnum), data=df2)
fit.6a <- lmer(df2[,i]~AGE + SEX+AGE:SEX + cell + AGE:cell + (1|patnum), data=df2)
anova(fit.5, fit.5a, fit.6, fit.6a)


df2$gender <- factor(ifelse(df2$SEX=='0', 'Men', 'Women'), levels=c("Women", "Men"))
df2$fitted <- fitted(fit.5)
df2$patient <- as.factor(df2$patnum)


##plot
library(ggplot2)
png("JGene_MixedEffects_cd4cd8_separate.png", res=200, height=800, width=1400)
p <-  ggplot(df2, aes(AGE, df2[,i], group = patient, colour=gender))+
  geom_point(aes(shape = patient), size =3) +
  geom_line(aes(AGE, fitted)) +
  facet_wrap(~CD4_CD8)+
  labs(x="Age", y="Sequence Diversity", title="Diversity ~ AGE + SEX + AGE*SEX +  AGE*cells + donors") +
  theme_bw()+
  theme(axis.title =element_text(size=22),
        strip.text= element_text(size=22),
        axis.text = element_text(size=22))
plot(p)
dev.off()
