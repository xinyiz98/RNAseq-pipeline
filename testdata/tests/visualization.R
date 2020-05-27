library("scater")
library("RUVSeq")
library("TCC")
library("edgeR")
library("splatter")
library("ggplot2")
library("reshape")

#load test data
setwd("C://UCBERKELEY//ASlab//rna-seq//testdata")
normtests<-dir("normtestresults",full.names = TRUE,no..=TRUE)
normtestsnames<-dir("normtestresults",full.names = FALSE,no..=TRUE)
ntests<-length(normtests)

#load simulated data
simfiles<-list.files("simresults",full.names = TRUE)
simfilenames<-list.files("simresults")
empsimfiles<-list.files("empsimresults",full.names = TRUE)
empsimfilenames<-list.files("empsimresults")
sims<-c(empsimfiles,simfiles)
simnames<-c(empsimfilenames,simfilenames)

#extract defaults
#test data
d<-c("sim_default","empsim_default")
didx<-which(normtestsnames %in% d)
defaultnames<-normtestsnames[didx]
defaults<-normtests[didx]
tidx<-which(!(normtestsnames %in% d))
normtests<-normtests[tidx]
normtestsnames<-normtestsnames[tidx]
#sim data
simdefaults<-sims[didx]
sims<-sims[tidx]
simnames<-simnames[tidx]

#create directory to store test results
if(!dir.exists("visualizedresults")){
  dir.create("visualizedresults")
}

i<-1
while(i<=(ntests-2)){
# while(i<=2){
  setwd("C://UCBERKELEY//ASlab//rna-seq//testdata")
  testi<-c(normtests[i],normtests[i+1])
  testnamei<-c(normtestsnames[i],normtestsnames[i+1])
  coeflist<-2:4#edger coefs
  # contrastlist<-list(c("group","Group1","Group2"),c("group","Group1","Group3"),c("group","Group1","Group4")) #deseq2 contrasts
  
  #get corresponding default
  if(strsplit(testnamei[1],"_")[[1]][[1]]=="empsim"){
      defaultipath<-defaults[1]
      simdefaultipath<-simdefaults[1]
  }else{
    defaultipath<-defaults[2]
    simdefaultipath<-simdefaults[2]
  }
  defaulti<-list.files(defaultipath,full.names = TRUE)
  defaultnamesi<-list.files(defaultipath)
  
  #get sim data
  simi<-readRDS(sims[i])
  simi1<-readRDS(sims[i+1])
  simidef<-readRDS(simdefaultipath)
  simnamei<-c(simnames[i],simnames[i+1])
  
  #get truths
  truthi<-simi@featureData@data[,c("DEFacGroup2","DEFacGroup3","DEFacGroup4")]
  refi<-simi@featureData@data[,"DEFacGroup1"]
  for(g in 1:ncol(truthi)){
    nde<- which(truthi[,g]==refi)
    de<- which(truthi[,g]!=refi)
    truthi[,g][nde]<-rep(0,length(nde)) #non DE genes are marked as 0, DE genes as 1
    truthi[,g][de]<-rep(1,length(de))
  }
  truthi1<-simi1@featureData@data[,c("DEFacGroup2","DEFacGroup3","DEFacGroup4")]
  refi1<-simi1@featureData@data[,"DEFacGroup1"]
  for(g in 1:ncol(truthi1)){
    nde<- which(truthi1[,g]==refi1)
    de<- which(truthi1[,g]!=refi1)
    truthi1[,g][nde]<-rep(0,length(nde))
    truthi1[,g][de]<-rep(1,length(de))
  }
  truthidef<-simidef@featureData@data[,c("DEFacGroup2","DEFacGroup3","DEFacGroup4")]
  refidef<-simidef@featureData@data[,"DEFacGroup1"]
  for(g in 1:ncol(truthidef)){
    nde<- which(truthidef[,g]==refidef)
    de<- which(truthidef[,g]!=refidef)
    truthidef[,g][nde]<-rep(0,length(nde))
    truthidef[,g][de]<-rep(1,length(de))
  }
  truthlist<-list(truthi,truthidef,truthi1)
  
  methodnames<-list.files(testi[1])
  methodi<-list.files(testi[1],full.names = TRUE)
  methodi1<-list.files(testi[2],full.names = TRUE)
  methodsi<-list(methodi,methodi1)
  
  #expected fdr levels for testing
  interval<-0.001
  efdrlevels<-seq(from=0,to=0.1,by=interval)
  
  #construct list of data frames to store values vs efdr
  x<-list(rep(0,length(efdrlevels)))
  TFDR<-data.frame(rep(x,length(methodnames)))
  colnames(TFDR)<-methodnames
  TFDR<-list("h"=TFDR,"def"=TFDR,"l"=TFDR)
  TPR<-FPR<-ACC<-TFDR
  
  ##construct lists to store ppv, tfdr at fdr threshold
  # ppvlist0<-rep(0,length(methodnames))
  # ppvlist<-list(ppvlist0,ppvlist0,ppvlist0)
  # tfdrlist<-ppvlist
  
  predlist<-list(methodi,defaulti,methodi1)
  for(j in 1:length(predlist)){
    predj<-predlist[[j]]
    truthj<-truthlist[[j]]
    for(k in 1:length(methodnames)){
      predjk<-predj[k]
      predjk<-readRDS(predjk) #done glmFit
      methodnamek<-methodnames[k]
      
      
      #omit filtered genes in truthi
      truthidx<-which(rownames(truthj) %in% rownames(predjk$counts))
      truthk<-truthj[truthidx,]
      
      #construct prediction dataframe
      pred_pval<-data.frame(matrix(ncol = length(coeflist), nrow = nrow(truthk)))
      #edger lrt
      for(c in 1:length(coeflist)){
        predjklrt<-glmLRT(predjk,coeflist[c])
        pred_pval[,c]<-predjklrt$table[,"PValue"]
      }
      
      #calculate TPR, FPR, ACC, TFDR
      for(f in 1:length(efdrlevels)){
        p<-truthk==1 #condition positive
        np<-sum(p)
        predp<-pred_pval<efdrlevels[f] #prediction positive
        tp<-predp&p #true positive
        ntp<-sum(tp)
        tpr<-ntp/np #true positive rate
        TPR[[j]][f,k]<-tpr
        n<-truthk==0 #condition negative
        nn<-sum(n)
        fp<-predp&n #false positive
        nfp<-sum(fp)
        fpr<-nfp/nn #false positive rate
        FPR[[j]][f,k]<-fpr
        predn<-pred_pval>=efdrlevels[f] #prediction negative
        tn<-predn&n #true negative
        ntn<-sum(tn)
        acc<-(ntp+ntn)/(nn+np) #accuracy
        ACC[[j]][f,k]<-acc
        if((ntp+nfp)==0){
          tfdr<-0
        }else{
          tfdr<-nfp/(ntp+nfp) #true false discovery rate; tfdr=0 if denominator is 0
        }
        TFDR[[j]][f,k]<-tfdr
      }
      
      # # calculate positive predictive value at fdr 0.1; ppv=pp/p
      # fdrthreshold<-0.1
      # pp<-(pred_pval<0.1)&(truthk==1)
      # npp<-sum(pp==TRUE)
      # pp<-pred_pval<0.1
      # npp<-sum(pp==TRUE)
      # ppv<-npp/np
      # ppvlist[[j]][k]<-ppv
      # 
      # # calculate true fdr at fdr=0.1; tfdr=fp/p
      # fp<-(pred_pval<0.1)&(truthk==0)
      # nfp<-sum(fp==TRUE)
      # tfdr<-nfp/npp
      # tfdrlist[[j]][k]<-tfdr
      
    }
  }
  
  #plot roc curve
  setwd("C://UCBERKELEY//ASlab//rna-seq//testdata//visualizedresults")
  if(!dir.exists("roc")){
    dir.create("roc")
  }
  setwd("roc")
  for(pl in 1:length(TPR)){
    tprplot<-TPR[[pl]]
    fprplot<-FPR[[pl]]
    fprlong<-melt(fprplot)
    tprlong<-melt(tprplot)
    rocvar<-data.frame(method=fprlong[,"variable"],fprvalue=fprlong[,"value"],tprvalue=tprlong[,"value"])
    rocplot<-ggplot(data=rocvar,aes(x=fprvalue, y=tprvalue, colour=method)) +geom_line()
    rocplot<-rocplot+ggtitle(paste0(strsplit(testnamei[[1]],"_h")[[1]],"_",names(TPR)[pl]))
    ggsave(paste0(rocplot$labels$title,".jpeg"),plot = rocplot)
  }
  
  #plot comparison against efdr
  plotlist<-list(accuracy=ACC,trueFDR=TFDR,tprvalue=TPR)
  plotnames<-names(plotlist)
  for(y in 1:length(plotlist)){
    plotnamey<-plotnames[[y]]
    plotdata<-plotlist[[y]]
    setwd("C://UCBERKELEY//ASlab//rna-seq//testdata//visualizedresults")
    if(!dir.exists(plotnamey)){
      dir.create(plotnamey)
    }
    setwd(plotnamey)
    for(pl in 1:length(plotdata)){
      value<-plotdata[[pl]]
      valuelong<-melt(value)
      var<-data.frame(method=valuelong[,"variable"],efdrvalue=efdrlevels,value=valuelong[,"value"])
      valueplot<-ggplot(data=var,aes(x=efdrvalue, y=value, colour=method)) +geom_line()
      valueplot<-valueplot+ggtitle(paste0(strsplit(testnamei[[1]],"_h")[[1]],"_",names(plotdata)[pl])) 
      ggsave(paste0(valueplot$labels$title,".jpeg"),plot = valueplot)
    }
  }
 
  # # plot ppv
  # setwd("C://UCBERKELEY//ASlab//rna-seq//testdata//visualizedresults")
  # if(!dir.exists("ppv0.1")){
  #   dir.create("ppv0.1")
  # }
  # setwd("ppv0.1")
  # ppvlist<-as.data.frame(ppvlist)
  # colnames(ppvlist)<-c("h","def","l")
  # value<-list(ppvlist)
  # ppvlist<-data.frame(Methods=methodnames,value=value)
  # ppv.m<- melt(ppvlist, id.vars='Methods')
  # ppvplot<-ggplot(ppv.m, aes(Methods, value)) +geom_bar(aes(fill = variable), position = "dodge", stat="identity")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # ppvplot<-ppvplot+ggtitle(strsplit(testnamei[[1]],"_h")[[1]])
  # ggsave(paste0(ppvplot$labels$title,".jpeg"),plot = ppvplot)
  # 
  # #plot tfdr
  # setwd("C://UCBERKELEY//ASlab//rna-seq//testdata//visualizedresults")
  # if(!dir.exists("tfdr0.1")){
  #   dir.create("tfdr0.1")
  # }
  # setwd("tfdr0.1")
  # tfdrlist<-as.data.frame(tfdrlist)
  # colnames(tfdrlist)<-c("h","def","l")
  # value<-list(tfdrlist)
  # tfdrlist<-data.frame(Methods=methodnames,value=value)
  # tfdr.m<- melt(tfdrlist, id.vars='Methods')
  # tfdrplot<-ggplot(tfdr.m, aes(Methods, value)) +geom_bar(aes(fill = variable), position = "dodge", stat="identity")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # tfdrplot<-tfdrplot+ggtitle(strsplit(testnamei[[1]],"_h")[[1]])
  # ggsave(paste0(tfdrplot$labels$title,".jpeg"),plot = tfdrplot)
  
  i<-i+2
}
