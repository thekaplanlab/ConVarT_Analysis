
# !NOT RUN! #

analysis2<-function(analysis, sp1, sp2) {
  
  mouse_aap<-analysis[[1]]
  mouse_rfsq<-analysis[[2]]
  mouse_sqsqsq<-analysis[[3]]
  mouse_ind<-analysis[[4]]
  
  mutagen_aap<-analysis[[5]]
  mutagen_rfsq<-analysis[[6]]
  mutagen_sqsqsq<-analysis[[7]]
  mutagen_ind<-analysis[[8]]
  
  clinvar_aap<-analysis[[9]]
  clinvar_rfsq<-analysis[[10]]
  clinvar_sqsqsq<-analysis[[11]]
  clinvar_ind<-analysis[[12]]
  
  dbsnp_aap<-analysis[[13]]
  dbsnp_rfsq<-analysis[[14]]
  dbsnp_sqsqsq<-analysis[[15]]
  dbsnp_ind<-analysis[[16]]
  
  cosmic_aap<-analysis[[17]]
  cosmic_rfsq<-analysis[[18]]
  cosmic_sqsqsq<-analysis[[19]]
  cosmic_ind<-analysis[[20]]
  
  gnomad_aap<-analysis[[21]]
  gnomad_rfsq<-analysis[[22]]
  gnomad_sqsqsq<-analysis[[23]]
  gnomad_ind<-analysis[[24]]
  
  celegans_aap<-analysis[[25]]
  celegans_rfsq<-analysis[[26]]
  celegans_sqsqsq<-analysis[[27]]
  celegans_ind<-analysis[[28]]
  
  
  u_sp1_ind<-unique(get(paste0(sp1, "_ind")))
  u_sp2_ind<-unique(get(paste0(sp2, "_ind")))
  commonind<-u_sp2_ind[which(u_sp2_ind %fin% u_sp1_ind)]
  
  sp1_common_aa<-c()
  sp1_common_refseq<-c()
  sp2_common_aa<-c()
  sp2_common_refseq<-c()
  pb = txtProgressBar(min = 0, max = max(commonind), initial = 0)
  for (i in commonind){
    
    # Mouse - Celegans
    anumber<- Reduce(intersect, list(get(paste0(sp1, "_sqsqsq"))[which(get(paste0(sp1, "_ind")) %fin% i)], get(paste0(sp2, "_sqsqsq"))[which(get(paste0(sp2, "_ind")) %fin% i)]))
    
    if (length(anumber > 0)){
      xsp1_common_aa<-get(paste0(sp1, "_aap"))[which((get(paste0(sp1, "_ind")) %fin% i) & (get(paste0(sp1, "_sqsqsq")) %fin% anumber))]
      xsp1_common_refseq<-get(paste0(sp1, "_rfsq"))[which((get(paste0(sp1, "_ind")) %fin% i) & (get(paste0(sp1, "_sqsqsq")) %fin% anumber))]
      
      xsp2_common_aa<-get(paste0(sp2, "_aap"))[which((get(paste0(sp2, "_ind")) %fin% i) & (get(paste0(sp2, "_sqsqsq")) %fin% anumber))]
      xsp2_common_refseq<-get(paste0(sp2, "_rfsq"))[which((get(paste0(sp2, "_ind")) %fin% i) & (get(paste0(sp2, "_sqsqsq")) %fin% anumber))]
      
      sp1_common_aa<-c(sp1_common_aa, xsp1_common_aa)
      sp1_common_refseq<-c(sp1_common_refseq, xsp1_common_refseq)
      sp2_common_aa<-c(sp2_common_aa, xsp2_common_aa)
      sp2_common_refseq<-c(sp2_common_refseq, xsp2_common_refseq)
    }
    
    progress(i, max(commonind), init = TRUE)
    if (i == (max(commonind)+1)) message("Done!")
    #setTxtProgressBar(pb,i)
  }
  
  lst<-list(sp1_common_aa, sp1_common_refseq, sp2_common_aa, sp2_common_refseq)
  return(lst)
}
