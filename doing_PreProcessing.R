#function for diversity calculation
Div <- function(d) {
  x <- (entropy(d)-(length(d)-1)/(2*sum(d)))/log(length(d))
  return(x)
}

#function for creating Functional and NonFunctional set of sequences
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


desc <- read.table("sample_description_Dean.txt", sep="\t", header = T)
source("PreProcessing_human.R")

#list of 587 samples - Dean et al
lf <- list.files(pattern =".tsv")

## outliers
ktoreNA <- which(is.na(desc$SEX) | is.na(desc$AGE))
lf <- lf[-ktoreNA]
desc <- desc[-ktoreNA]

#based on Brouffaerts criterion:
outliers <- c("HIP03111.tsv", "HIP03678.tsv", "HIP14118.tsv", "HIP14138.tsv", "HIP16867.tsv");
ktore <- which(lf %in% outliers)
lf=lf[-ktore]
desc <- desc[-ktoreNA]

summary <- matrix(NA, ncol=4, nrow=length(lf), 
                dimnames=list(lf, c("prod.unique", "prod.total", "non.unique", "non_total")))
status <- matrix(NA, ncol=1, nrow=length(lf), 
                  dimnames=list(lf, c("diversity")))
sequence.prod <- matrix(NA, ncol=1, nrow=length(lf), 
                 dimnames=list(lf, c("diversity")))
sequence.non <- matrix(NA, ncol=1, nrow=length(lf), 
                        dimnames=list(lf, c("diversity")))

for(i in lf) {
  data = read.table(i, sep="\t", header = T)
  wyniki <- PreProcessing(data)
  
  data <- wyniki[[1]]
  
  #saving preprocessed datasets
  write.table(data, file=paste0('preprocessed/xx',i), quote=FALSE, row.names=F, col.names=TRUE, sep="\t", na="")
  
  
  FunIn <- Curate(data)[[1]]
  FunOut <- Curate(data)[[2]]
  
  #status diversity
  counts <- c(sum(as.numeric(FunIn[,3])), sum(as.numeric(FunOut[,3])))
  status[i] <- Div(counts)
  
  
  # sequence diversity
  sequence.prod[i] <- Div(as.numeric(FunIn[,3]))
  sequence.non[i] <- Div(as.numeric(FunOut[,3]))
  
  
  summary[i,1] = nrow(FunIn)
  summary[i,2] = sum(FunIn[,3])
  summary[i,3] = nrow(FunOut)
  summary[i,4] = sum(FunOut[,3])
  write.table(summary, "preprocessed/summaryprodvsnon.txt", sep = "\t", row.names = T, col.names = T)
  write.table(status, "preprocessed/statusDiv.txt", sep = "\t", row.names = T, col.names = T)
  write.table(sequence.prod, "preprocessed/prodSeqDiv.txt", sep = "\t", row.names = T, col.names = T)
  write.table(sequence.non, "preprocessed/nonSeqDiv.txt", sep = "\t", row.names = T, col.names = T)
  write.table(desc, "preprocessed/sample_description_Dean.txt", sep = "\t", row.names = T, col.names = T)
  
}
