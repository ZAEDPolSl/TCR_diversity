
PreProcessing <- function(A) {

    sumup <- matrix(NA, ncol=16, nrow=1, 
                    dimnames=list(c("count"), c("before.u","before.t", "functionality.u", "functionality.t", 
                                    "Vorphon.u", "Vorphon.t",
                                    "V-unres.u","V-unres.t", "J-unres.u","J-unres.t", 
                                    "Vind>N.u","Vind>N.t","V3.u","V3.t", "after.u", "after.t")))
    sumup[1] = nrow(A)
    sumup[2] = sum(as.numeric(A[,3]))
    prod=which(A[,39]=="In")
    
   
    #human Vs
    badV = c("TCRBV01-01", "TCRBV03-02", "TCRBV05-02", "TCRBV05-03", "TCRBV05-07", "TCRBV06-07",
             "TCRBV07-01", "TCRBV07-05", "TCRBV08-01", "TCRBV08-02", "TCRBV12-01", "TCRBV12-02", 
             "TCRBV17-01", "TCRBV21-01", "TCRBV22-01", "TCRBV23-01", "TCRBV26-01", 
             "TCRBVA", "TCRBVB", "TCRBVC") ## tutaj wszystkie mozliwe allele są niefunkcjonalne!
    badV_al2 = c("TCRBV07-03", "TCRBV07-04", "TCRBV16-01")
    badV_al3 = c("TCRBV07-03", "TCRBV10-01", "TCRBV30-01")
    
    badJ1 = "TCRBJ02-02P"
    badJ= "TCRBJ02-07" # ORF only in the second allele
 
    
    #productiveness of genes
    y1=which((A[prod,'vGeneName'] %in% badV) == TRUE)
    if(length(y1)>0) A[prod[y1], 'sequenceStatus']="NonFun"
    y2=which((A[prod,'vGeneName'] %in% badV_al2) == TRUE & A[prod, 'vGeneAllele']==2)
    if(length(y2)>0) A[prod[y2], 'sequenceStatus']="NonFun"
    y3=which((A[prod,'vGeneName'] %in% badV_al3) == TRUE & A[prod, 'vGeneAllele']==3)
    if(length(y3)>0) A[prod[y3], 'sequenceStatus']="NonFun"
    
    j1=which((A[prod,'jGeneName'] %in% badJ1) == TRUE)
    if(length(j1)>0) A[prod[j1], 'sequenceStatus']="NonFun"
    j2=which((A[prod,'jGeneName'] %in% badJ) == TRUE & A[prod, 'jGeneAllele']==2)
    if(length(j2)>0) A[prod[j2], 'sequenceStatus']="NonFun"

    
    nf=which(A[,'sequenceStatus']=="NonFun")
    sumup[3] = length(nf)
    sumup[4] = sum(as.numeric(A[nf,3]))
    sumup[15] = sum(j1,j2)
    
    #orphon V genes
    ktore <- grep("or", A[,'vGeneName'])
    sumup[5] = length(ktore)
    sumup[6] = sum(as.numeric(A[ktore,3]))
    if(length(ktore)>0) A <- A[-ktore,]
    
    
    #V UNRESOLVED
    x<-which(A[,'vGeneName']=="unresolved")
    sumup[7] = length(x)
    sumup[8] = sum(as.numeric(A[x,3]))
    if(length(x)!=0) A <- A[-x,]

    

    #J-unresolved
    x<-which(A[,22]=="unresolved")
    sumup[9] = length(x)
    sumup[10] = sum(as.numeric(A[x,3]))
    if(length(x)!=0) A<-A[-x,]

    
    #V ind >= N2 ind
    x<-which(as.numeric(A[,33])>=as.numeric(A[,34])&as.numeric(A[,34])!=-1)
    sumup[11] = length(x)
    sumup[12] = sum(as.numeric(A[x,3]))
    if(length(x)!=0) A <- A[-x,]
    
    ## TCRBV03-01 gene
    x <- which(A[,8]=="TCRBV03-01")
    sumup[13] = length(x)
    sumup[14] = sum(as.numeric(A[x,3]))
    if(length(x)!=0) A<-A[-x,]
    

    
    sumup[15] = nrow(A)
    sumup[16] = sum(as.numeric(A[,3]))
    
    return(list(A, sumup))

    write.table(A, file=nazwa, quote=FALSE, row.names=F, col.names=TRUE, sep="\t", na="")
    
}