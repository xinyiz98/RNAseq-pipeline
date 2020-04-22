library("tximport")
library("sleuth")
library("biomaRt")
library("limma")
library("DESeq2")
library("edgeR")
library("VennDiagram")

# setwd("/Users/Xinyi")

alignmethod<-c("salmon","kallisto")
reftxpt<-c("ucsc","ensembl")
batch<-c("Batch1","Batch3")
condition<-list(c("T2","T7"),c("T7","T24"))
method_id0<-lapply(alignmethod,paste0,reftxpt)
method_id0<-c(method_id0[[1]],method_id0[[2]])
method_id<-lapply(method_id0,paste0,batch)
method_id<-c(lapply(method_id,"[[",1),lapply(method_id,"[[",2))
sample_id<-lapply(method_id,list.dirs,recursive=FALSE)
for(i in 1:length(sample_id)){
  conditioni<-condition[[floor((i-1)/4)+1]]
  sample_id_i<-vector("list",length = 0)
  for(j in 1:length(conditioni)){
    conditionij<-conditioni[[j]]
    if(conditionij=="T2"){
      wt2<-grepl("T2",sample_id[[i]],fixed=TRUE)
      t24<-grepl("T24",sample_id[[i]],fixed=TRUE)
      idx<-which(sample_id[[i]][wt2] %in% sample_id[[i]][(!t24)])
      sample_id_j<-sample_id[[i]][wt2][idx]
      sample_id_i<-c(sample_id_i,sample_id_j[!is.na(sample_id_j)])
    }else{
      idx<-grepl(conditionij,sample_id[[i]],fixed=TRUE)
      sample_id_j<-sample_id[[i]][idx]
      sample_id_i<-c(sample_id_i,sample_id_j)
    }
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
    
  #DE
  #construct design matrix
  dmtx<-data.frame(condition=rep(condition[[floor((i-1)/4)+1]],each=length(sample_id[[i]])/2))
  #loop over de methods
  for(j in 1:length(de_method)){
    if(de_method[[j]]=="voom"){
      #import
      txi<-tximport(as.character(sample_path),type = aligni,ignoreTxVersion = TRUE,tx2gene = ttg,countsFromAbundance = "lengthScaledTPM")
      dge<-DGEList(txi$counts,group = dmtx$condition,remove.zeros = TRUE)
      dge<-calcNormFactors(dge)
      design=model.matrix(~condition,data = dmtx)
      v<-voom(dge,design=design)
      vfit<-lmFit(v,design)
      vfit<-eBayes(vfit)
      de_res$voom[[i]]<-as.data.frame(vfit$p.value[,2])
    }
    if(de_method[[j]]=="deseq2"){
      #import
      txi<-tximport(as.character(sample_path),type = aligni,ignoreTxVersion = TRUE,tx2gene = ttg)
      dds<-DESeqDataSetFromTximport(txi,colData = dmtx,design = ~condition)
      dds <- dds[ rowSums(counts(dds)) > 0, ] #remove rows with zero counts
      dds<-DESeq(dds)
      ddsres<-results(dds)
      ddspvalue<-as.data.frame(ddsres[,"padj"])
      rownames(ddspvalue)<-rownames(ddsres)
      de_res$deseq2[[i]]<-ddspvalue
    }
  }
}

setwd("simAlign/results/deagreement")
#DE agreement
for(i in 1:length(de_res)){
  names(de_res[[i]])<-method_id
  dei<-de_res[[i]]
  deiname<-names(de_res)[[i]]
  ncomp<-length(dei)/length(batch) #number of components in the diagram
  comparison<-c(pvalue05=0.05,pvalue1=0.1,lNDE=200)
  if(deiname=="deseq2"){
    comparison<-c(comparison,outlier="outlier")
  }
  for(j in 1:length(batch)){
    vennlist<-dei[((j-1)*ncomp+1):(j*ncomp)]
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
      names(vennlist_input)<-names(dei)[((j-1)*ncomp+1):(j*ncomp)]
      venn.diagram(vennlist_input,filename = paste0(comparisonname,batch[[j]],deiname,".tiff"),width = 4000,
                   fill = c("cornflowerblue", "green", "yellow", "darkorchid1"))
    }
  }
}
