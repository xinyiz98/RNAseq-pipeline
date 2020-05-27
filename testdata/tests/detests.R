library("TCC")
library("edgeR")
library("DESeq2")
library("IHW")
library("limma")
library("samr")

#load data
setwd("C://UCBERKELEY//ASlab//rna-seq//testdata")
simresults<-list.files("simresults",full.names = TRUE)
simresultsnames<-list.files("simresults")
empsimresults<-list.files("empsimresults",full.names = TRUE)
empsimresultsnames<-list.files("empsimresults")
sims<-c(simresults,empsimresults)
simnames<-c(simresultsnames,empsimresultsnames)
ndata<-length(sims)

#create directory to store test results
if(!dir.exists("DEtestresults")){
  dir.create("DEtestresults")
}

# for(i in 1:ndata){
for(i in 1:ndata){
  setwd("C://UCBERKELEY//ASlab//rna-seq//testdata")
  datai<-readRDS(sims[i])
  counti<-counts(datai) #extract count matrix
  groupi<-pData(datai)$Group #extract treatment group info
  groupi<-as.factor(groupi)
  dge<-DGEList(counts = counti,group = groupi,remove.zeros = TRUE)
  
  #filter
  keep <- rowSums(cpm.DGEList(dge)>1) >= 6
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  counti<-dge$counts
  
  setwd("DEtestresults")
  #create subdirectory for test i
  testname<-strsplit(simnames[i],".rds")[[1]]
  if(!dir.exists(testname)){
    dir.create(testname)
  }
  setwd(testname)
  
  #set parameters
  #construct bayseq parameter, groups; models based on start and end times of DE
  NDE<-c(1,1,1,1)
  DE14<-c(1,2,3,4)
  DE13<-c(1,2,3,3)
  DE12<-c(1,2,2,2)
  DE24<-c(1,1,2,3)
  DE23<-c(1,1,2,2)
  DE34<-c(1,1,1,2)
  groups<-list(NDE=NDE,DE14=DE14,DE13=DE13,DE12=DE12,DE24=DE24,DE23=DE23,DE34=DE34)
  groups<-lapply(groups,rep,each=(length(groupi)/4))
  
  #Normalization w/TMM-edgeR-TMM pipeline in TCC
  design<-model.matrix(~groupi) #for tcc normalization pipeline 
  tcc<-new("TCC",dge$counts,groupi)
  tcc<-calcNormFactors(tcc,norm.method="tmm",test.method="edger",iteration=1,design=design,FDR=0.1)
  tccnorm<-getNormalizedData(tcc)
  dgetcc<-DGEList(tccnorm,group = groupi)
  tccnorm<-apply(tccnorm,c(1,2),as.integer)
  # #tests on different DE methods
  # DEmethods<-c("edger","deseq2","bayseq","voom")
  # for(j in 1:length(DEmethods)){
  #   samplesize<-100
  #   if(DEmethods[j]=="bayseq"){
  #     #bayseq parameter samplesize
  #     if(i>23){
  #       samplesize<-500
  #     }else{
  #       samplesize<-10000
  #     }
  #   }
  #   if(DEmethods[j]=="deseq2"){
  #     design<-~group
  #   }else{
  #     design<-model.matrix(~groupi)
  #   }
  #   deresult<-estimateDE(tcc,test.method = DEmethods[j], FDR = 0.1,coef = tmmcoef,design = design,replicates=groupi,groups=groups,
  #                        samplesize = samplesize,resp.type="Multiclass")
  #   deresult<-getResult(deresult,sort = TRUE)
  #   filename<-paste0("tcc",DEmethods[j],".rds")
  #   saveRDS(deresult,file = filename)
  # }
  
  #Normalization w/default methods
  dge<-calcNormFactors(dge,method="TMM")
  
  #DE analysis
  #edgeR
  dgeedger<-estimateDisp(dge,design)
  fit<-glmFit(dgeedger,design)
  saveRDS(fit,file = "defedger.rds")
  print("edgerdef")
  tccedger<-estimateDisp(dgetcc,design)
  tccfit<-glmFit(tccedger,design)
  saveRDS(tccfit,file = "tccedger.rds")
  print("edgertcc")
  #DESeq2
  coldata<-data.frame(groupi,row.names = colnames(counti))
  colnames(coldata)<-"group"
  dds<-DESeqDataSetFromMatrix(counti,colData = coldata,design = ~group)
  dds <- dds[ rowSums(counts(dds)) > 0, ] #remove rows with zero counts
  dds<-DESeq(dds)
  saveRDS(dds,file = "defdeseq2.rds")
  # deseqres<-results(dds,contrast = as.list(resultsNames(dds)[2]))
  # saveRDS(deseqres,file = "defdeseq2.rds")
  # #DESeq2 + lfcShrink
  # deseqcoef<-length(resultsNames(dds))
  # deseqlfc<-lfcShrink(dds,coef = deseqcoef,res = deseqres)
  # saveRDS(deseqlfc,file = "defdeseq2lfc.rds")
  # #DESeq2 + IHW
  # deseqihw<-results(dds,filterFun = ihw)
  # saveRDS(deseqihw,file = "defdeseq2ihw.rds")
  print("deseq2def")
  tccdds<-DESeqDataSetFromMatrix(tccnorm,colData = coldata,design = ~group)
  tccdds<-tccdds[rowSums(counts(tccdds))>0,]
  tccdds<-DESeq(tccdds)
  saveRDS(tccdds,file = "tccdeseq2.rds")
  print("deseq2tcc")
  #baySeq
  if(i>23){
    samplesize<-500
  }else{
    samplesize<-10000
  }
  cd<-new("countData",data=counti,replicates=groupi,groups=groups)
  libsizes(cd)<-getLibsizes(cd)
  cd<-getPriors.NB(cd,samplesize = samplesize,cl=NULL)
  cd<-getLikelihoods(cd,verbose = FALSE)
  saveRDS(cd,file = "defbayseq.rds")
  print("bayseqdef")
  tcccd<-new("countData",data=tccnorm,replicates=groupi,groups=groups)
  libsizes(tcccd)<-colSums(tccnorm)
  tcccd<-getPriors.NB(tcccd,samplesize = samplesize,cl=NULL)
  tcccd<-getLikelihoods(tcccd,verbose = FALSE)
  saveRDS(tcccd,file = "tccbayseq.rds")
  print("bayseqtcc")
  #voom+limma
  v<-voom(dge,design = design,plot = FALSE)
  vfit<-lmFit(v,design)
  vfit<-eBayes(vfit)
  saveRDS(vfit,file = "deflmvoom.rds")
  print("voomdef")
  tccv<-voom(dgetcc,design = design,plot = FALSE)
  tccvfit<-lmFit(tccv,design)
  tccvfit<-eBayes(tccvfit)
  saveRDS(tccvfit,file = "tcclmvoom.rds")
  print("voomtcc")
  #SAMseq
  samseq<-SAMseq(counti,groupi,resp.type = "Multiclass",fdr.output = 0.1)
  saveRDS(samseq,file = "defsamseq.rds")
  print("samseq")
  
  print(i)
}