library("tximport")
library("sleuth")
library("biomaRt")
library("limma")
library("DESeq2")
library("edgeR")
library("VennDiagram")

setwd('adipocytes_early')

alignmethod<-c("salmon","kallisto")
reftxpt<-c("ucsc","ensembl")
batch<-c("B1","B2")
condition<-c("Day.1","Day.5")
method_id0<-lapply(alignmethod,paste0,reftxpt)
method_id0<-c(method_id0[[1]],method_id0[[2]])
method_id<-lapply(method_id0,paste0,batch)
# method_id<-c(lapply(method_id,"[[",1),lapply(method_id,"[[",2))
sample_id<-lapply(method_id,list.dirs,recursive=FALSE)
for(i in 1:length(sample_id)){
  sample_id_i<-vector("list",length = 0)
  for(j in 1:length(condition)){
    conditionj<-condition[[j]]
    idx<-grepl(conditionj,sample_id[[i]],fixed=TRUE)
    sample_id_j<-sample_id[[i]][idx]
    sample_id_i<-c(sample_id_i,sample_id_j)
  }
  sample_id[[i]]<-sample_id_i
}

#transcripts to genes
ensembl=useMart("ensembl")
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
ttgucsc <- biomaRt::getBM(attributes = c("refseq_mrna", "ensembl_gene_id"),mart = ensembl)
ttgucsc <- dplyr::rename(ttgucsc, target_id = refseq_mrna,ens_gene = ensembl_gene_id)
ttgensembl<-biomaRt::getBM(attributes = c("ensembl_transcript_id","ensembl_gene_id"),mart = ensembl)
ttgensembl <- dplyr::rename(ttgensembl, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id)

#import transcripts and DE
de_res<-vector("list",length = length(sample_id))
de_method<-c("voom","deseq2")
de_res<-list(de_res,de_res)
names(de_res)<-de_method
for(i in 1:length(sample_id)){
  #import transcripts
  #select transcriptome
  print(method_id[[i]])
  if(grepl("ensembl",sample_id[[i]][[1]],fixed = TRUE)){
    ttg<-ttgensembl
    print("ensembl")
  }else{
    ttg<-ttgucsc
    print("ucsc")
  }
  ttg<-ttg[!duplicated(ttg[,1]),]
  #select alignment method
  if((i-1)%%4<2){
    aligni<-"salmon"
    sample_path<-lapply(sample_id[[i]],file.path,"quant.sf")
  }else{
    aligni<-"kallisto"
    sample_path<-lapply(sample_id[[i]],file.path,"abundance.tsv")
  }
  print(aligni)
    
  #DE
  #construct design matrix
  dmtx<-data.frame(condition=rep(condition,c(2,3)),batch=c('B1','B2','B1','B1','B2'))
  #loop over de methods
  for(j in 1:length(de_method)){
    if(de_method[[j]]=="voom"){
      #import
      txi<-tximport(as.character(sample_path),type = aligni,ignoreTxVersion = TRUE,tx2gene = ttg,countsFromAbundance = "lengthScaledTPM")
      dge<-DGEList(txi$counts,group = dmtx$condition,remove.zeros = TRUE)
      dge<-calcNormFactors(dge)
      design=model.matrix(~condition+batch,data = dmtx)
      v<-voom(dge,design=design)
      vfit<-lmFit(v,design)
      vfit<-eBayes(vfit)
      de_res$voom[[i]]<-vfit
    }
    if(de_method[[j]]=="deseq2"){
      #import
      txi<-tximport(as.character(sample_path),type = aligni,ignoreTxVersion = TRUE,tx2gene = ttg)
      dds<-DESeqDataSetFromTximport(txi,colData = dmtx,design = ~condition+batch)
      dds <- dds[ rowSums(counts(dds)) > 0, ] #remove rows with zero counts
      dds<-DESeq(dds)
      # ddsres<-results(dds,contrast = c('condition','Day.1','Day.5'))
      # ddspvalue<-as.data.frame(ddsres[,"padj"])
      # rownames(ddspvalue)<-rownames(ddsres)
      de_res$deseq2[[i]]<-dds
    }
  }
}

setwd('..')
setwd("simAlign/results/3T3L1/")

#DE agreement excluding genes not shared
for(i in 1:length(de_res)){
  names(de_res[[i]])<-method_id0
  dei<-de_res[[i]]
  deiname<-names(de_res)[[i]]
  ncomp<-length(dei)/length(batch) #number of components in the diagram
  comparison<-c(pvalue05=0.05,pvalue1=0.1,lNDE=200)
  comparison_lfc<-c(lfc0=0,lfc1=1)
  if(deiname=="deseq2"){
    comparison<-c(comparison,outlier="outlier")
  }

  vennlist<-dei
  commonGene<-rownames(vennlist[[1]])
  for(k in 2:length(vennlist)){
    commonGene<-commonGene[which(commonGene %in% rownames(vennlist[[k]]))]
  }
  for(c in 1:length(comparison)){
    vennlist_input<-vector("list",length = length(vennlist))
    comparisonc<-comparison[[c]]
    comparisonname<-names(comparison)[[c]]
    if(grepl("pvalue",comparisonname,fixed = TRUE)){
      for(k in 1:length(vennlist)){
        vk<-vennlist[[k]]
        notna<-!is.na(vk)
        vk_input<-as.data.frame(vk[notna,]) #filter pvalue = NA
        rownames(vk_input)<-rownames(vk)[notna]
        bpv<-vk_input<comparisonc # get DE genes
        bpvrownames<-rownames(vk_input)[bpv]
        vk_input<-as.data.frame(vk_input[bpv,])
        rownames(vk_input)<-bpvrownames
        vennlist_input[[k]]<-rownames(vk_input)
      }
    }else if(grepl("lNDE",comparisonname,fixed = TRUE)){
      for(k in 1:length(vennlist)){
        vk<-vennlist[[k]]
        notna<-!is.na(vk)
        vk_input<-as.data.frame(vk[notna,]) #filter pvalue = NA
        rownames(vk_input)<-rownames(vk)[notna]
        pvoder<-order(-vk_input[,1]) #descending order wrt p-values
        pvrownames<-rownames(vk_input)[pvoder]
        vk_input<-data.frame(pvalue=vk_input[pvoder,],row.names = pvrownames)
        vennlist_input[[k]]<-rownames(vk_input)[1:comparisonc]
      }
    }else if(comparisonc=="outlier"){
      for(k in 1:length(vennlist)){
        vk<-vennlist[[k]]
        na<-is.na(vk)
        vk_input<-as.data.frame(vk[na,]) #filter pvalue = NA
        rownames(vk_input)<-rownames(vk)[na]
        vennlist_input[[k]]<-rownames(vk_input)
      }
    }
    
    #draw venn diagram
    names(vennlist_input)<-c(names(dei))
    venn.diagram(vennlist_input,filename = paste0('deagreement/',comparisonname,deiname,".jpeg"),width = 4000,
                 fill = c("cornflowerblue", "green", "yellow", "darkorchid1"))
    
    ##subset by commonGene
    for(k in 1:length(vennlist_input)){
      vk<-vennlist_input[[k]]
      vennlist_input[[k]]<-vk[which(vk %in% commonGene)]
    }
    #draw venn diagram after subset
    names(vennlist_input)<-c(names(dei))
    venn.diagram(vennlist_input,filename = paste0('deagreement_commonGene/',comparisonname,deiname,".jpeg"),width = 4000,
                 fill = c("cornflowerblue", "green", "yellow", "darkorchid1"))
    
  }
}
