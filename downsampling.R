## CODE FOR CALCULATING DIVERSITY BASED ON RECONSTITUTED SAMPLE REGARDING CD4/CD8 RATIO

Div <- function(d) {
  x <- (entropy(d)-(length(d)-1)/(2*sum(d)))/log(length(d))
  return(x)
}

Curate <- function(dane) {
  
  prod=which(dane[,'sequenceStatus']=="In")
  FunIn = dane[prod,]
  V3 <- which(FunIn$vGeneName=="TCRBV03-01")
  if(length(V3)>0) FunIn = FunIn[-V3,]
  
  ## data preprocessing
  dane=dane[-prod,]
  FunOut=dane[which((dane[,8] %in% badV) == FALSE),]
  ktore <- which((FunOut[,'vGeneName'] %in% badV_al2)==TRUE & FunOut[,'vGeneAllele']==2) 
  print(length(ktore))
  if(length(ktore)>0) FunOut = FunOut[-ktore,]
  ktore <- which((FunOut[,'vGeneName'] %in% badV_al3)==TRUE & FunOut[,'vGeneAllele']==3) 
  print(length(ktore))
  if(length(ktore)>0) FunOut = FunOut[-ktore,]
  
  ktore <- which((FunOut[,'jGeneName'] %in% badJ)==TRUE) 
  print(length(ktore))
  if(length(ktore)>0) FunOut = FunOut[-ktore,]
  
  j2 <- which((FunOut[,'jGeneName'] %in% badJ2) == TRUE & FunOut[, 'jGeneAllele']==2)
  print(length(j2))
  if(length(j2)>0) FunOut = FunOut[-ktore,]
  
  V3 <- which(FunOut$vGeneName=="TCRBV03-01")
  if(length(V3)>0) FunOut = FunOut[-V3,]
  
  return(list(FunIn, FunOut))
}

source('PreProcessing_human.R')


library('entropy')
library('boot')
library(ggplot2)
library(plyr)



badV = c("TCRBV01-01", "TCRBV03-02", "TCRBV05-02", "TCRBV05-03", "TCRBV05-07", "TCRBV06-07",
         "TCRBV07-01", "TCRBV07-05", "TCRBV08-01", "TCRBV08-02", "TCRBV12-01", "TCRBV12-02", 
         "TCRBV17-01", "TCRBV21-01", "TCRBV22-01", "TCRBV23-01", "TCRBV26-01", 
         "TCRBVA", "TCRBVB", "TCRBVC") ## tutaj wszystkie mozliwe allele są niefunkcjonalne!
badV_al2 = c("TCRBV07-03", "TCRBV07-04", "TCRBV16-01")
badV_al3 = c("TCRBV07-03", "TCRBV10-01", "TCRBV30-01")
badJ = c("TCRBJ02-02P")
badJ2 = c("TCRBJ02-07")

##reading in ratio for each donor
ratio <- read.table("CD4CD8ratio.txt", sep="\t", header = T, dec=",", stringsAsFactors = F, row.names = 1)

#parameter for number of iterations
k=100

#resulting files -  number of columns corresponds to number of iterations
sum_un_In <- matrix(ncol=k, nrow=nrow(ratio),
                 dimnames = list(rownames(ratio), c(1:k)))
sum_tot_In <- matrix(ncol=k, nrow=nrow(ratio),
                    dimnames = list(rownames(ratio), c(1:k)))
sum_un_Out <- matrix(ncol=k, nrow=nrow(ratio),
                    dimnames = list(rownames(ratio), c(1:k)))
sum_tot_Out <- matrix(ncol=k, nrow=nrow(ratio),
                    dimnames = list(rownames(ratio), c(1:k)))

status <- matrix(ncol=k, nrow=nrow(ratio),
                 dimnames = list(rownames(ratio), c(1:k)))
Vdiv.prod <- matrix(ncol=k, nrow=nrow(ratio),
                    dimnames = list(rownames(ratio), c(1:k)))
Jdiv.prod <- matrix(ncol=k, nrow=nrow(ratio),
                    dimnames = list(rownames(ratio), c(1:k)))
Vdiv.non <- matrix(ncol=k, nrow=nrow(ratio),
                   dimnames = list(rownames(ratio), c(1:k)))
Jdiv.non <- matrix(ncol=k, nrow=nrow(ratio),
                   dimnames = list(rownames(ratio), c(1:k)))
sequence.prod <- matrix(ncol=k, nrow=nrow(ratio),
                        dimnames = list(rownames(ratio), c(1:k)))
sequence.non <- matrix(ncol=k, nrow=nrow(ratio),
                       dimnames = list(rownames(ratio), c(1:k)))
ile.zdup <- matrix(ncol=k, nrow=nrow(ratio),
                       dimnames = list(rownames(ratio), c(1:k)))


zdup <- list()

for (i in 1:nrow(pliki)) { 
  ## input datasets  - dane4: data for CD4+; dane8 - data for CD8+
  dane4=read.table(pliki[i,1], header=TRUE, sep="\t", stringsAsFactors = FALSE)
  dane8=read.table(pliki[i,2], header=TRUE, sep="\t", stringsAsFactors = FALSE)
  
  n4 <- sum(dane4[,3]) # counts here
  n8 <- sum(dane8[,3])
  
  #up- or down-sampling regadirng the cd4/cd8 ratio
  x4 <- round(ratio[i]*n8)
  x8 <- round(n4/ratio[i])
  
  #if x4 < n4, then downsampling of CD4 sequences required
  #if x4 > n4, then x8 < n8, then downsampling of CD8 sequences required
  cum8 <- c(0, cumsum(dane8[,3]))
  cum4 <- c(0, cumsum(dane4[,3]))
  
  for(kk in 1:k) {
       ## DOWNSAMPLING 
    if(x4<n4) {
      los <- table(cut(sample(1:n4, x4, replace=F), breaks = cum4))
      tmp <- dane4
      tmp[,3] <- as.numeric(los)
      ktore0 <- which(tmp[,3]==0)
      if(length(ktore0)>0) tmp <- tmp[-ktore0,]
      dane <- rbind(dane8,tmp)
    } else {
      los <- table(cut(sample(1:n8, x8, replace=F), breaks = cum8))
      tmp <- dane8
      tmp[,3] <- as.numeric(los)
      ktore0 <- which(tmp[,3]==0)
      if(length(ktore0)>0) tmp <- tmp[-ktore0,]
      dane <- rbind(dane4,tmp)
    }
    
    
    ## filtering of duplicated nulecotide sequences
    dup <- which(duplicated(dane$nucleotide, fromLast = T)==1)
    zdup[[kk]] <- dup
    ile.zdup[i,kk]<- length(dup)
    
    if(length(dup>0)){
      #only unique seuqences
      for(dd in dup) {
        dane[dd,3] = dane[dd,3] + dane$count[which(dane$nucleotide==dane$nucleotide[dd])[2]]
      }
      dane <- dane[-which(duplicated(dane$nucleotide)==1),]
    }
    
    ## data preprocessing
    dane <- PreProcessing(dane)[[1]]
    
    FunIn <- Curate(dane)[[1]]
    FunOut <- Curate(dane)[[2]]
    
    sum_un_In[i,kk]=nrow(FunIn)
    sum_tot_In[i,kk]=sum(as.numeric(FunIn[,3]))
    sum_un_Out[i,kk]=nrow(FunOut)
    sum_tot_Out[i,kk]=sum(as.numeric(FunOut[,3]))
    
    
    counts <- c(sum(as.numeric(FunIn[,3])), sum(as.numeric(FunOut[,3])))
    status[i,kk] <- Div(counts)
    
    
    # sequence diversity
    sequence.prod[i,kk] <- Div(as.numeric(FunIn[,3]))
    sequence.non[i,kk] <- Div(as.numeric(FunOut[,3]))

    ## V gene diversity
    # total counts of genes:
    counts=as.matrix(table(FunIn[,8], FunIn[,3]))
    a=as.numeric(colnames(counts))
    counts=as.matrix(colSums(a*t(counts)))
    Vdiv.prod[i,kk] <- Div(as.numeric(counts))

    counts=as.matrix(table(FunOut[,8], FunOut[,3]))
    a=as.numeric(colnames(counts))
    counts=as.matrix(colSums(a*t(counts)))
    Vdiv.non[i,kk] <- Div(as.numeric(counts))

    ## J gene diversity
    counts=as.matrix(table(FunIn[,'jGeneName'], FunIn[,3]))
    a=as.numeric(colnames(counts))
    counts=as.matrix(colSums(a*t(counts)))
    Jdiv.prod[i,kk] <- Div(as.numeric(counts))

    counts=as.matrix(table(FunOut[,'jGeneName'], FunOut[,3]))
    a=as.numeric(colnames(counts))
    counts=as.matrix(colSums(a*t(counts)))
    Jdiv.non[i,kk] <- Div(as.numeric(counts))

    
    ##saving files
    write.table(sum_un_In, file='sum_unIn.txt', quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t", na="")
    write.table(sum_un_Out, file='sum_unOut.txt', quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t", na="")
    write.table(sum_tot_In, file='sum_totIn.txt', quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t", na="")
    write.table(sum_tot_Out, file='sum_totOut.txt', quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t", na="")
    
    write.table(status, file='status.txt', quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t", na="")
    write.table(Vdiv.prod, file='V_prod.txt', quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t", na="")
    write.table(Vdiv.non, file='V_non', quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t", na="")
    write.table(Jdiv.prod, file='J_prod.txt', quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t", na="")
    write.table(Jdiv.non, file='J_non.txt', quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t", na="")
    write.table(sequence.prod, file='sequence_prod.txt', quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t", na="")
    write.table(sequence.non, file='sequence_non.txt', quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t", na="")
  }
 
}

