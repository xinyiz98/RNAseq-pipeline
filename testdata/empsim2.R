library("splatter")

#load data; Bottomly
botcounts <- read.table("emp//bottomly_count_table.txt", header = TRUE, as.is = TRUE)
rownames<-botcounts[,1]
botcounts <- as.matrix(botcounts[, -1])
rownames(botcounts)<-rownames

#estimate parameters
botparams<-splatEstimate(botcounts)
botparams<-setParam(botparams,"dropout.present",FALSE)
nGroups<-4

#set parameters;last element of a vector is the default
groupCells<-c(3,9,6)
de.prob<-c(0.05,0.3,0.1)
de.downProb<-c(0.25,0.75,0.5)
de.facLoc<-c(0.05,0.3,0.1)
de.facScale<-c(0.2,0.6,0.4)
#construct data frame for variable parameters
varparams<-data.frame(groupCells,de.prob,de.downProb,de.facLoc,de.facScale)
rownames(varparams)<-c("l","h","default")

#extract defaults to list
defaults<-as.list(varparams["default",])
botparams<-setParams(botparams,update = defaults)

#simulation
nsim<-ncol(varparams)*(nrow(varparams)-1)+1 #number of simulations
for(i in 1:nsim){
  print(i)
  parami<-botparams #set back to defaults
  if(i==nsim){
    group<-getParam(parami,"groupCells") #set different treatment groups
    parami<-setParam(parami,"groupCells",rep(group,each=nGroups))
    empsimi<-splatSimulate(parami,verbose = FALSE,seed=as.integer(runif(1,max=1000000)),method = "groups")
    saveRDS(empsimi,file = "empsimresults//empsim_default.rds")
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
  empsimi<-splatSimulate(parami,verbose = FALSE,seed=as.integer(runif(1,max=1000000)),method = "groups")
  filename<- paste("empsimresults//","empsim_",parname,"_",rowi,".rds",sep = "")
  saveRDS(empsimi,file = filename)
}
