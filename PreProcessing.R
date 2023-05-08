
PreProcessing <- function(A) {

    sumup <- matrix(NA, ncol=2, nrow=8, 
                    dimnames=list(c("before","functionality", "V-unres", "J-unres", "orphon", "Vind>D","TCRBV03-01", "after"), 
                                                      c("unique", "total")))
    sumup[1, 1] = nrow(A)
    sumup[1, 2] = sum(as.numeric(A[,3]))
    prod=which(A[,39]=="In")
    
    A[,39]=factor(A[,39], levels=c(levels(A[,39]), "NonFun"))

    # #mouse_V
    # badV=c("TCRBV06-01","TCRBV07-01","TCRBV08-01","TCRBV09-01","TCRBV10-01","TCRBV11-01","TCRBV12-03","TCRBV18-01","TCRBV22-01","TCRBV24-01","TCRBV25-01","TCRBV27-01","TCRBV28-01")
    # 
    #human Vs
    badV = c("TCRBV01-01", "TCRBV03-02", "TCRBV05-02", "TCRBV05-03", "TCRBV05-07", "TCRBV06-07",
             "TCRBV07-01", "TCRBV07-05", "TCRBV08-02", "TCRBV12-01", "TCRBV12-02", "TCRBV17-01",
             "TCRBV21-01", "TCRBV22-01", "TCRBV23-01", "TCRBV26-01")
    badV_al2 = c("TCRBV06-02", "TCRBV07-03", "TCRBV07-04")
    badV_al3 = c("TCRBV06-02", "TCRBV07-03", "TCRBV10-01", "TCRBV30-01")
    
    
    for (i in 1:length(badV)) {
        ktore=which(A[prod,'vGeneName']==badV[i])
        if (length(ktore)>0) A[prod[ktore],'sequenceStatus']="NonFun"
    }
    for (i in 1:length(badV_al2)) {
        ktore=which(A[prod,'vGeneName']==badV_al2[i] & A[prod,'vGeneAllele']==2)
        if (length(ktore)>0) A[prod[ktore],'sequenceStatus']="NonFun"
    }
    for (i in 1:length(badV_al3)) {
        ktore=which(A[prod,'vGeneName']==badV_al3[i] & A[prod,'vGeneAllele']==3)
        if (length(ktore)>0) A[prod[ktore],'sequenceStatus']="NonFun"
    }
    
    x<-which(A[,'sequenceStatus']=="NonFun")

    # # mouse Js
    # badJ=c("TCRBJ01-07", "TCRBJ02-06")
    
    # human Js
    badJ = c("TCRBJ02-07")
    
    for (i in 1:length(badJ)) {
        ktore=which(A[prod,'jGeneName']==badJ[i] & A[prod,'jGeneAllele']==2)
        if (length(ktore)>0) A[prod[ktore],'sequenceStatus']="NonFun"
    }
    
    x<-which(A[,'sequenceStatus']=="NonFun")
  
    nf=which(A[,'sequenceStatus']=="NonFun")
    sumup[2,1] = length(nf)
    sumup[2,2] = sum(as.numeric(A[nf,3]))
    

    
    #V unresolved genes
    x<-which(A[,8]=="unresolved")
    sumup[3,1] = length(x)
    sumup[3,2] = sum(as.numeric(A[x,3]))
    if(length(x)!=0) A <- A[-x,]

    

    # J unresolved genes
    x<-which(A[,22]=="unresolved")
    sumup[4,1] = length(x)
    sumup[4,2] = sum(as.numeric(A[x,3]))
    if(length(x)!=0) A<-A[-x,]
    
    
    ## orphon genes
    x <- grep("or", A[,8])
    sumup[5,1] = length(x)
    sumup[5,2] = sum(as.numeric(A[x,3]))
    if(length(x)!=0) A<-A[-x,]
    
    
    ## sequences with V ind >= N2 ind
    x<-which(as.numeric(A[,33])>=as.numeric(A[,34])&as.numeric(A[,34])!=-1)
    sumup[6,1] = length(x)
    sumup[6,2] = sum(as.numeric(A[x,3]))
    if(length(x)!=0) A <- A[-x,]
    
    
    ## TCRBV03-01 gene
    x <- which(A[,8]=="TCRBV03-01")
    sumup[7,1] = length(x)
    sumup[7,2] = sum(as.numeric(A[x,3]))
    if(length(x)!=0) A<-A[-x,]  ## odkomentować, żeby usunąc ten gen!!

    
    sumup[8, 1] = nrow(A)
    sumup[8, 2] = sum(as.numeric(A[,3]))
    
    return(list(A, sumup))

    write.table(A, file=nazwa, quote=FALSE, row.names=F, col.names=TRUE, sep="\t", na="")
    
}