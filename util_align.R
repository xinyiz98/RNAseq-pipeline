suppressPackageStartupMessages(library("polyester"))
suppressPackageStartupMessages(library("seqinr"))
suppressPackageStartupMessages(library("Biostrings"))

simulateLib_polyester<-function(refPath,ngroup,ndata,deprob,savepath_fcmtx,savepath_lib){
  
  # read in transcriptome
  ucscfasta<-readDNAStringSet(refPath)
  
  
  ntxpt<-length(ucscfasta)
  if(!dir.exists("simfcmtx")){
    dir.create("simfcmtx")
  }
  for(i in 1:ndata){
    fc<-rep(1,ntxpt*ngroup) # create fold change matrix
    for(j in 1:length(fc)){
      urn<-runif(1,min = 0,max = 1)
      if(urn<deprob){
        fcj<-rpois(1,demean)
        fcj<-ceiling(fcj)
        if(fcj<1){
          fcj<-1
        }
        fc[j]<-fcj
      }
    }
    fc<-matrix(fc,nrow = ntxpt)
    saveRDS(fc,file = paste0(savepath_fcmtx,"/",toString(i),".rds")) #save fc mtx as truth
    # simulate_experiment(fasta = "reftx/reftxmrna_ucsc_refgene.fasta.gz",
    # outdir = "simfasta",num_reps = rep(1,ngroup),
    # meanmodel=TRUE,fold_changes = fc)
    simulate_experiment(fasta = refPath,
                        outdir = paste0(savepath_lib,toString(i)),num_reps = rep(1,ngroup),
                        reads_per_transcript = round(20 * width(ucscfasta) / 100),
                        fold_changes = fc,error_model="illumina5")
  }
  
  
  
}

getLibPairs<-function(libID_len,dataDir_fq){}