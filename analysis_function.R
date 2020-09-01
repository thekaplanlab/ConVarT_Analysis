
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
  
  mouse_aap<-mouse_aap[!is.na(mouse_sqsqsq)]
  mouse_rfsq<-mouse_rfsq[!is.na(mouse_sqsqsq)]
  mouse_sqsqsq<-mouse_sqsqsq[!is.na(mouse_sqsqsq)]
  mouse_ind<-mouse_ind[!is.na(mouse_sqsqsq)]
  
  celegans_aap<-celegans_aap[!is.na(celegans_sqsqsq)]
  celegans_rfsq<-celegans_rfsq[!is.na(celegans_sqsqsq)]
  celegans_sqsqsq<-celegans_sqsqsq[!is.na(celegans_sqsqsq)]
  celegans_ind<-celegans_ind[!is.na(celegans_sqsqsq)]
  
  u_sp1_ind<-unique(get(paste0(sp1, "_ind")))
  u_sp2_ind<-unique(get(paste0(sp2, "_ind")))
  commonind<-u_sp2_ind[which(u_sp2_ind %fin% u_sp1_ind)]
  
  u1<-get(paste0(sp1, "_aap"))
  u2<-get(paste0(sp1, "_rfsq"))
  u3<-get(paste0(sp2, "_aap"))
  u4<-get(paste0(sp2, "_rfsq"))
  
  sp1_common_aa<-numeric(length(u1))
  sp1_common_refseq<-character(length(u2))
  sp2_common_aa<-numeric(length(u3))
  sp2_common_refseq<-character(length(u4))
  
  sp1_common_aa_length<-0
  sp1_common_refseq_length<-0
  sp2_common_aa_length<-0
  sp2_common_refseq_length<-0
  
  for (i in commonind){
    
    # Mouse - Celegans
    anumber<- Reduce(intersect, list(get(paste0(sp1, "_sqsqsq"))[which(get(paste0(sp1, "_ind")) %fin% i)], get(paste0(sp2, "_sqsqsq"))[which(get(paste0(sp2, "_ind")) %fin% i)]))
    
    if (length(anumber > 0)){
      xsp1_common_aa<-get(paste0(sp1, "_aap"))[which((get(paste0(sp1, "_ind")) %fin% i) & (get(paste0(sp1, "_sqsqsq")) %fin% anumber))]
      xsp1_common_refseq<-get(paste0(sp1, "_rfsq"))[which((get(paste0(sp1, "_ind")) %fin% i) & (get(paste0(sp1, "_sqsqsq")) %fin% anumber))]
      
      xsp2_common_aa<-get(paste0(sp2, "_aap"))[which((get(paste0(sp2, "_ind")) %fin% i) & (get(paste0(sp2, "_sqsqsq")) %fin% anumber))]
      xsp2_common_refseq<-get(paste0(sp2, "_rfsq"))[which((get(paste0(sp2, "_ind")) %fin% i) & (get(paste0(sp2, "_sqsqsq")) %fin% anumber))]
      
      sp1_common_aa[(sp1_common_aa_length+1) : (sp1_common_aa_length+length(xsp1_common_aa))]<-xsp1_common_aa
      sp1_common_refseq[(sp1_common_refseq_length+1) : (sp1_common_refseq_length+length(xsp1_common_refseq))]<-xsp1_common_refseq
      sp2_common_aa[(sp2_common_aa_length+1) : (sp2_common_aa_length+length(xsp2_common_aa))]<-xsp2_common_aa
      sp2_common_refseq[(sp2_common_refseq_length+1) : (sp2_common_refseq_length+length(xsp2_common_refseq))]<-xsp2_common_refseq
      
      sp1_common_aa_length<-sp1_common_aa_length + length(xsp1_common_aa)
      sp1_common_refseq_length<-sp1_common_refseq_length + length(xsp1_common_refseq)
      sp2_common_aa_length<-sp2_common_aa_length + length(xsp2_common_aa)
      sp2_common_refseq_length<-sp2_common_refseq_length + length(xsp2_common_refseq)
    }
    
    progress(i, max(commonind), init = TRUE)
    if (i == max(commonind)) message("Done!")
    #setTxtProgressBar(pb,i)
  }
  sp1_common_aa<-sp1_common_aa[sp1_common_aa != 0]
  sp1_common_refseq<-sp1_common_refseq[sp1_common_refseq != ""]
  sp2_common_aa<-sp2_common_aa[sp2_common_aa != 0]
  sp2_common_refseq<-sp2_common_refseq[sp2_common_refseq != ""]
  
  lst<-list(sp1_common_aa, sp1_common_refseq, sp2_common_aa, sp2_common_refseq)
  return(lst)
}
