library("scater")
library("RUVSeq")
library("TCC")
library("edgeR")
library("splatter")

#load data
setwd("C://UCBERKELEY//ASlab//rna-seq//testdata")
simresults<-list.files("simresults",full.names = TRUE)
simresultsnames<-list.files("simresults")
empsimresults<-list.files("empsimresults",full.names = TRUE)
empsimresultsnames<-list.files("empsimresults")
sims<-c(simresults,empsimresults)
simnames<-c(simresultsnames,empsimresultsnames)
# methodnames<-c("TPM","FPKM","TMM","RLE","UQ","RUVg+TMM","RUVg+RLE","RUVg+UQ","DEGES+TMM","DEGES+RLE","DEGES+TbT")
ndata<-length(sims)

#create directory to store test results
if(!dir.exists("normtestresults")){
  dir.create("normtestresults")
}

#normalization and DE; loop over all data sets
for(i in 1:ndata){
  setwd("C://UCBERKELEY//ASlab//rna-seq//testdata")
  datai<-readRDS(sims[i])
  counti<-counts(datai) #extract count matrix
  groupi<-pData(datai)$Group #extract treatment group info
  groupi<-as.factor(groupi)
  design<-model.matrix(~groupi)
  
  setwd("normtestresults")
  #create subdirectory for test i
  testname<-strsplit(simnames[i],".rds")[[1]]
  if(!dir.exists(testname)){
    dir.create(testname)
  }
  setwd(testname)
  
  #TPM
  if(i<24){
    tpmi<-tpm(datai)
    tpmdge<-DGEList(counts = tpmi,group = groupi,remove.zeros = TRUE)
    tpmkeep <- rowSums(cpm.DGEList(tpmdge)>1) >= 6 #filtering
    tpmdge <- tpmdge[tpmkeep, , keep.lib.sizes=FALSE]
    tpmdge<-estimateDisp(tpmdge,design) #common and tagwise dispersion
    tpmfit<-glmFit(tpmdge,design)
    # tpmLrt<-glmLRT(tpmfit,coef = coef)
    saveRDS(tpmfit,file="tpm.rds")
  }
  
  #FPKM
  if(i<24){
    fpkmi<-fpkm.SCESet(datai)
    fpkmdge<-DGEList(counts = fpkmi,group = groupi,remove.zeros = TRUE)
    fpkmkeep <- rowSums(cpm.DGEList(fpkmdge)>1) >= 6
    fpkmdge <- fpkmdge[fpkmkeep, , keep.lib.sizes=FALSE]
    fpkmdge<-estimateDisp(fpkmdge,design)
    fpkmfit<-glmFit(fpkmdge,design)
    # fpkmLrt<-glmLRT(fpkmfit,coef = coef)
    saveRDS(fpkmfit,file = "fpkm.rds")
  }
  
  #dgelist from raw counts
  dge<-DGEList(counts = counti,group = groupi,remove.zeros = TRUE)
  #filtering
  keep <- rowSums(cpm.DGEList(dge)>1) >= 6
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  
  #TMM
  tmmdge<-calcNormFactors(dge,method="TMM")
  tmmdge<-estimateDisp(tmmdge,design)
  tmmfit<-glmFit(tmmdge,design)
  # tmmLrt<-glmLRT(tmmfit,coef = coef)
  saveRDS(tmmfit,file = "tmm.rds")
  
  #RLE
  rledge<-calcNormFactors(dge,method="RLE")
  rledge<-estimateDisp(rledge,design)
  rlefit<-glmFit(rledge,design)
  # rleLrt<-glmLRT(rlefit,coef = coef)
  saveRDS(rlefit,file = "rle.rds")

  #UQ
  uqdge<-calcNormFactors(dge,method="upperquartile")
  uqdge<-estimateDisp(uqdge,design)
  uqfit<-glmFit(uqdge,design)
  # uqLrt<-glmLRT(uqfit,coef = coef)
  saveRDS(uqfit,file = "uq.rds")

  #RUVg; k=1
  ruvset<-newSeqExpressionSet(counti,phenoData = data.frame(groupi,row.names = colnames(counti)))
  #add empirical control genes - bottom 75% by DE; normalized by TMM
  ngenes<-nrow(counti)
  threshold<-as.integer(ngenes*0.25)
  tagrank<- topTags(tmmLrt,n=ngenes)$table
  cgenes<-rownames(ruvset)[threshold:ngenes]
  #calculate
  ruvgset<-RUVg(ruvset,cgenes,k=1)
  ruvdesign<-model.matrix(~groupi+W_1,data = pData(ruvgset))
  #+TMM
  ruvtmm<-calcNormFactors(dge,method="TMM")
  ruvtmm<-estimateDisp(ruvtmm,ruvdesign)
  ruvtmmfit<-glmFit(ruvtmm,ruvdesign)
  # ruvtmmLrt<-glmLRT(ruvtmmfit,coef = coef)
  saveRDS(ruvtmmfit,file = "ruvtmm.rds")
  #+RLE
  ruvrle<-calcNormFactors(dge,method="RLE")
  ruvrle<-estimateDisp(ruvrle,ruvdesign)
  ruvrlefit<-glmFit(ruvrle,ruvdesign)
  # ruvrleLrt<-glmLRT(ruvrlefit,coef = coef)
  saveRDS(ruvrlefit,file = "ruvrle.rds")
  #+UQ
  ruvuq<-calcNormFactors(dge,method="upperquartile")
  ruvuq<-estimateDisp(ruvuq,ruvdesign)
  ruvuqfit<-glmFit(ruvuq,ruvdesign)
  # ruvuqLrt<-glmLRT(ruvuqfit,coef = coef)
  saveRDS(ruvuqfit,file = "ruvuq.rds")
  
  #DEGES
  tcc<-new("TCC",dge$counts,groupi)
  #+tmm
  degestmm<-calcNormFactors(tcc,norm.method="tmm",test.method="edger",iteration=1,design=design)
  degestmm<-getNormalizedData(degestmm)
  degestmm<-DGEList(degestmm,group = groupi)
  degestmm<-estimateDisp(degestmm,design)
  degestmmresult<-glmFit(degestmm,design)
  saveRDS(degestmmresult,file = "degestmm.rds")
  #+rle(deseq)
  degesrle<-calcNormFactors(tcc,norm.method = "deseq", test.method = "deseq",iteration = 1)
  degesrle<-getNormalizedData(degesrle)
  degesrle<-DGEList(degesrle,group = groupi)
  degesrle<-estimateDisp(degesrle,design)
  degesrleresult<-glmFit(degesrle,design)
  saveRDS(degesrleresult,file = "degesrle.rds")
  #+TbT
  samplesize<-10000
  if(i>23){
    samplesize<-500
  }
  degestbt<- calcNormFactors(tcc, norm.method = "tmm", test.method = "bayseq",iteration = 1, samplesize = samplesize)
  degestbt<-getNormalizedData(degestbt)
  degestbt<-DGEList(degestbt,group = groupi)
  degestbt<-estimateDisp(degestbt,design)
  degestbtresult<-glmFit(degestbt,design)
  saveRDS(degestbtresult,file = "degestbt.rds")
  
  print(i)
}
