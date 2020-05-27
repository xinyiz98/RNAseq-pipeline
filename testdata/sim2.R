library("splatter")

#set parameters
#fixed
ngenes<-1132
nGroups<-4
dropoutpresent<-FALSE
params<-newSplatParams(nGenes=ngenes,dropout.present=dropoutpresent)
#last element of a vector is the default
groupCells<-c(3,9,6)
mean.shape<-c(0.4,0.8,0.6)
mean.rate<-c(0.1,0.5,0.3)
lib.loc<- c(6,16,11)
lib.scale<-c(0.15,0.5,0.35)
out.prob<-c(0.01,0.1,0.05)
de.prob<-c(0.05,0.3,0.1)
de.downProb<-c(0.25,0.75,0.5)
de.facLoc<-c(0.05,0.3,0.1)
de.facScale<-c(0.2,0.6,0.4)
bcv.common<-c(0.05,0.2,0.1)
#construct data frame for variable parameters
varparams<-data.frame(groupCells,mean.shape,mean.rate,lib.loc,lib.scale,out.prob,de.prob,de.downProb,de.facLoc,de.facScale,bcv.common)
rownames(varparams)<-c("l","h","default")

#extract defaults to list
defaults<-as.list(varparams["default",])

#simulation
nsim<-ncol(varparams)*(nrow(varparams)-1)+1 #number of simulations
for(i in 1:nsim){
  print(i)
  parami<-setParams(params,update = defaults) #set back to defaults
  if(i==nsim){
    group<-getParam(parami,"groupCells") #set different treatment groups
    parami<-setParam(parami,"groupCells",rep(group,each=nGroups))
    simi<-splatSimulate(parami,verbose = FALSE,seed=as.integer(runif(1,max=1000000)),method = "groups")
    simi<-addGeneLengths(simi) #simulating gene lengths
    tpm(simi)<-calculateTPM(simi,fData(simi)$Length) #calculate TPM
    fpkm(simi)<-calculateFPKM(simi,fData(simi)$Length) #calculate FPKM
    saveRDS(simi,file = "simresults//sim_default.rds")
    break
  }
  coli<-ceiling(i/2) #update parameter
  if((i%%2)==1){
    rowi<-"l"
  }else{
    rowi<-"h"
  }
  parname<-colnames(varparams)[coli]
  parupd<-varparams[rowi,coli]
  parami<-setParam(parami,parname,parupd)
  group<-getParam(parami,"groupCells") 
  parami<-setParam(parami,"groupCells",rep(group,each=nGroups))
  simi<-splatSimulate(parami,verbose = FALSE,seed=as.integer(runif(1,max=1000000)),method="groups")
  simi<-addGeneLengths(simi)
  tpm(simi)<-calculateTPM(simi,fData(simi)$Length)
  fpkm(simi)<-calculateFPKM(simi,fData(simi)$Length)
  filename<- paste("simresults//","sim_",parname,"_",rowi,".rds",sep = "")
  saveRDS(simi,file = filename)
}

