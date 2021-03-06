library("tximport")
library("biomaRt")
library("limma")
library("VennDiagram")
library('edgeR')

tx2gene<-function(alignmethod,reftxpt,batch,condition,allconditions){
  #return gene count matrix for one condition in one batch
  #col: sample in specified condition; row: gene
  #batch: name of batch directory (string); condition: one condition in the batch (string)
  
  # setwd("/Users/Xinyi")
  setwd(sourcedir)
  method_id<-lapply(alignmethod,paste0,reftxpt)
  method_id<-lapply(method_id,paste0,batch)
  sample_id<-lapply(method_id,list.dirs,recursive=FALSE)
  print(sample_id)
  #select condition from batch
  pattern<-vector("list",length(sample_id[[1]]))
  patternlist<-list(condition=pattern)
  patternlevel<-list(condition=allconditions)
  for(j in 1:length(patternlist)){
    pj<-patternlevel[[j]]
    for(i in 1:length(pj)){
      pposition<-grepl(pj[[i]],sample_id[[1]],fixed = TRUE)
      patternlist[[j]][pposition]<-pj[[i]]
    }
  }
  conditionpattern<-as.character(patternlist$condition)
  alignpattern<-patternlist$align
  sample_id[[1]]<-sample_id[[1]][conditionpattern==condition]
  # print(sample_id)
  
  sample_name<-lapply(sample_id,strsplit,"/")
  for(i in 1:length(sample_name)){
    sample_name[[i]]<-lapply(sample_name[[i]],"[[",2)
  }

  
  #import transcript abundance
  sample_path<-lapply(sample_id,file.path,"quant.sf")
  print(sample_path)
  txi<-tximport(sample_path[[1]],type = alignmethod,txOut = FALSE,ignoreTxVersion = TRUE,tx2gene = ttg,
                countsFromAbundance = "lengthScaledTPM")
  txptcounts<-txi$counts
  colnames(txptcounts)<-sample_name[[1]]
  
  return(txptcounts)
}

countdf<-function(alignmethod,reftxpt,batches,conditions,allconditions=condition){
  #return list of count data frames
  #batches: vector of batch names
  #conditions: vector of conditions; the first entry is the reference
  #allconditions: all conditions in batches
  onebatch<-function(alignmethod,reftxpt,batch,conditions,allconditions){
    counts<-vector('list',length = length(conditions))
    for(c in 1:length(conditions)){
      counts[[c]]<-tx2gene(alignmethod,reftxpt,batch,conditions[[c]],allconditions)
    }
    names(counts)<-conditions
    return(counts)
  }
  counts<-vector('list',length = length(batches))
  for(b in 1:length(batches)){
    counts[[b]]<-onebatch(alignmethod,reftxpt,batches[[b]],conditions,allconditions)
  }
  names(counts)<-batches
  return(counts)
}

de<-function(counts,contrast, demethod,adjust_method,nDEG='all',FDR=0.05){
  #DE results for one or multiple batches
  #counts: list of count matrices; each entry is a count matrix for a batch
  #contrast: a condition to be compared against reference
  #demethod: method for DE analysis; choose from limma, edgeR, DEseq2
  #nDEG: a positive integer or 'all'
  #return topTable result, sorted by adjusted p-value
  
  #construct group
  batches<-names(counts)
  conditionpattern<-vector()
  batchpattern<-vector()
  for(b in 1:length(batches)){
    conditionnames<-names(counts[[b]])
    nreplicates<-lapply(counts[[b]],ncol)
    print(nreplicates)
    conditions<-rep(conditionnames,nreplicates)
    conditionpattern<-c(conditionpattern,conditions)
    batchpattern<-c(batchpattern,rep(batches[[b]],length(conditions)))
  }
  
  print(conditionpattern)
  print(batchpattern)
  
  counts<-as.matrix(as.data.frame(counts))
  
  # counts<-normalizeLib(counts)
  if(nDEG=='all'){
    nDEG=nrow(counts)
  }
  
  conditionpattern<-factor(conditionpattern,levels = unique(conditionpattern))
  batchpattern<-factor(batchpattern,levels = unique(batches))
  #design matrix; includes batch pattern only when multiple batches exist
  if(length(batches)>1){
    design<-model.matrix(~batchpattern+conditionpattern)
  }else{
    design<-model.matrix(~conditionpattern)
  }
  # design<-model.matrix(~batchpattern)
  print(design)
  
  #DE
  if(demethod=='voom'){
    dge<-DGEList(counts,remove.zeros = TRUE)
    dge<-calcNormFactors(dge)
    v<-voom(dge,design=design)
    vfit<-lmFit(v,design)
    vfit<-eBayes(vfit)
    print(colnames(vfit))
    alldeg<-topTable(vfit,p.value = FDR,number = Inf,coef=paste0('conditionpattern',contrast),adjust.method = adjust_method)
    alldeg<-alldeg[order(alldeg$adj.P.Val),]
    alldeg<-alldeg[1:min(nDEG,nrow(alldeg)),]
  }else if (demethod=='DEseq2'){
    mode(counts)<-'integer'
    coldata<-data.frame(conditionpattern=conditionpattern,row.names = colnames(counts), batchpattern=batchpattern)
    dds<-DESeqDataSetFromMatrix(counts,colData = coldata,design =~batchpattern+conditionpattern)
    dds <- dds[ rowSums(counts(dds)) > 0, ] #remove rows with zero counts
    dds<-DESeq(dds)
    alldeg<-results(dds,contrast = c('conditionpattern',refcondition,contrast),pAdjustMethod =adjust_method)
    alldeg<-alldeg[which(alldeg[['padj']]<FDR),]
    alldeg<-alldeg[order(alldeg$padj),]
    alldeg<-alldeg[1:min(nDEG,nrow(alldeg)),]
  }else if (demethod=='edgeR'){
    dge<-DGEList(counts,remove.zeros = TRUE)
    dge<-calcNormFactors(dge)
    dge<-estimateDisp(dge,design)
    fit<-glmFit(dge,design)
    res<-glmLRT(fit,coef=paste0('conditionpattern',contrast))
    alldeg<-res$table
    alldeg$FDR<- p.adjust(alldeg$PValue, method=adjust_method)
    alldeg<-alldeg[which(alldeg[['FDR']]<FDR),]
    alldeg<-alldeg[order(alldeg$FDR),]
    alldeg<-alldeg[1:min(nDEG,nrow(alldeg)),]
  }else{
    print(paste('DE method',demethod,'not supported. Choose from edgeR, voom, DEseq2'))
  }
    
  return(alldeg)
}

drawvenn<-function(deresults,savepath,reference,contrast){
  # save venn diagrams of batches to savepath;
  # deresults: list of topTable results; results should be named
  
  temp<-vector('list',length = length(deresults))
  names(temp)<-names(deresults)
  venninputs<-list(allDEG=temp,down=temp,up=temp)
  for(r in 1:length(deresults)){
    venninputs$allDEG[[r]]<-rownames(deresults[[r]])
    venninputs$down[[r]]<-rownames(deresults[[r]][deresults[[r]]$logFC<0,])
    venninputs$up[[r]]<-rownames(deresults[[r]][deresults[[r]]$logFC>0,])
  }
  
  #create labels
  fn<-paste('REF',reference,'CONTRAST',contrast,sep = '_')
  
  #draw venn
  setwd(savepath)
  for(v in 1:length(venninputs)){
    venn.diagram(venninputs[[v]],filename = paste0(fn,'_',names(venninputs)[[v]],".tiff"),width = 4000,
                 main=paste0(fn,'_',names(venninputs)[[v]]))
  }
  
}

genelist<-function(deresult,savepath,filename, id2symb=NA, convID=NA){
  # gene id's with multiple symbols are annotated by the first symbol; other symbols are saved separately
  # add gene symbols to corresponding id's; store in the last column of deresult matrix
  # deresult: matrix; row names are gene id's
  if(is.na(id2symb)&is.na(convID)){
    setwd(savepath)
    write.csv(deresult,paste0(filename,'.csv'))
  }
  if(is.na(id2symb)){
    id2symb<-biomaRt::getBM(attributes = convID,mart = ensembl)
    # id2symb<- dplyr::rename(id2symb, id = ensembl_gene_id, symbol=mgi_symbol)  
  }
  duplicatedID<-id2symb[duplicated(id2symb[,convID[['id']]]),] #save other symbols
  id2symb<-id2symb[!duplicated(id2symb[,convID[['id']]]),] #remove duplicated id's
  rownames(id2symb)<-id2symb$id

  # add gene symbols to deresults
  symb<-id2symb[rownames(deresult),]
  deresult$mgi_symbol<-symb[,convID[['symbol']]]

  #save results
  setwd(savepath)
  write.csv(deresult,paste0(filename,'.csv'))
  write.csv(duplicatedID,paste0(filename,'_altsymbols.csv'))
}

filterCPM<-function(countsdf,threshCPM_byAvgColSum){
  countsdfCPM<-vector('list',length = length(countsdf))
  names(countsdfCPM)<-names(countsdf)
  
  passedGenes<-rep(TRUE,nrow(countsdf[[1]][[1]]))
  for(b in 1:length(countsdf)){
    countsdf_b<-countsdf[[b]]
    counts_b<-as.matrix(as.data.frame(countsdf_b))
    colsum<-colSums(counts_b)
    avgColsum<-mean(colsum)
    threshCPM<-avgColsum*threshCPM_byAvgColSum
    
    for(c in 1:length(countsdf_b)){
      passed_c<-(rowMeans(countsdf_b[[c]])>threshCPM)
      passedGenes<-(passed_c & passedGenes)
    }
    
  }
  for(b in 1:length(countsdf)){
    countsdfCPM_b<-vector('list',length = length(countsdf_b))
    names(countsdfCPM_b)<-names(countsdf_b)
    
    for(c in 1:length(countsdfCPM_b)){
      countsdfCPM_b[[c]]<-countsdf_b[[c]][passedGenes,]
    }
    
    countsdfCPM[[b]]<-countsdfCPM_b
  }
  return(countsdfCPM)
}

generateSimulatedCounts<-function(normalizationMethod,deMethod,groupCells,mean.shape,mean.rate,lib.loc,lib.scale,out.prob,de.prob,de.downProb,de.facLoc,bcv.common){}
