
analysis<-function(sp1, sp2) {
  
  u_sp1_ind<-unique(get(paste0(sp1, "_ind")))
  u_sp2_ind<-unique(get(paste0(sp2, "_ind")))
  commonind<-u_sp2_ind[which(u_sp2_ind %fin% u_sp1_ind)]
  
  sp1_common_aa<-c()
  sp1_common_refseq<-c()
  sp2_common_aa<-c()
  sp2_common_refseq<-c()
  pb = txtProgressBar(min = 0, max = length(commonind), initial = 0)
  for (i in commonind){
    
    # Mouse - Celegans
    anumber<- Reduce(intersect, list(get(paste0(sp1, "_sqsqsq"))[which(get(paste0(sp1, "_ind")) %fin% i)], get(paste0(sp2, "_sqsqsq"))[which(get(paste0(sp2, "_ind")) %fin% i)]))
    
    xsp1_common_aa<-get(paste0(sp1, "_aapos"))[which(get(paste0(sp1, "_sqsqsq")) %fin% anumber)]
    xsp1_common_refseq<-get(paste0(sp1, "_rfsq"))[which(get(paste0(sp1, "_sqsqsq")) %fin% anumber)]
    
    xsp2_common_aa<-get(paste0(sp2, "_aapos"))[which(get(paste0(sp2, "_sqsqsq")) %fin% anumber)]
    xsp2_common_refseq<-get(paste0(sp2, "_rfsq"))[which(get(paste0(sp2, "_sqsqsq")) %fin% anumber)]
    
    sp1_common_aa<-c(sp1_common_aa, xsp1_common_aa)
    sp1_common_refseq<-c(sp1_common_refseq, xsp1_common_refseq)
    sp2_common_aa<-c(sp2_common_aa, xsp2_common_aa)
    sp2_common_refseq<-c(sp2_common_refseq, xsp2_common_refseq)
    setTxtProgressBar(pb,i)
  }
  
  lst<-list(sp1_common_aa, sp1_common_refseq, sp2_common_aa, sp2_common_refseq)
  return(lst)
}
