
# !NOT RUN! #

library(dplyr)
library(fastmatch)
source("convart_helper.R")


orthologous_analysis<-function(analysis, sp1, sp2) {
  
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
  mouse_ind<-mouse_ind[!is.na(mouse_sqsqsq)]
  mouse_sqsqsq<-mouse_sqsqsq[!is.na(mouse_sqsqsq)]

  mutagen_aap<-mutagen_aap[!is.na(mutagen_sqsqsq)]
  mutagen_rfsq<-mutagen_rfsq[!is.na(mutagen_sqsqsq)]
  mutagen_ind<-mutagen_ind[!is.na(mutagen_sqsqsq)]
  mutagen_sqsqsq<-mutagen_sqsqsq[!is.na(mutagen_sqsqsq)]


  clinvar_aap<-clinvar_aap[!is.na(clinvar_sqsqsq)]
  clinvar_rfsq<-clinvar_rfsq[!is.na(clinvar_sqsqsq)]
  clinvar_ind<-clinvar_ind[!is.na(clinvar_sqsqsq)]
  clinvar_sqsqsq<-clinvar_sqsqsq[!is.na(clinvar_sqsqsq)]

  dbsnp_aap<-dbsnp_aap[!is.na(dbsnp_sqsqsq)]
  dbsnp_rfsq<-dbsnp_rfsq[!is.na(dbsnp_sqsqsq)]
  dbsnp_ind<-dbsnp_ind[!is.na(dbsnp_sqsqsq)]
  dbsnp_sqsqsq<-dbsnp_sqsqsq[!is.na(dbsnp_sqsqsq)]

  cosmic_aap<-cosmic_aap[!is.na(cosmic_sqsqsq)]
  cosmic_rfsq<-cosmic_rfsq[!is.na(cosmic_sqsqsq)]
  cosmic_ind<-cosmic_ind[!is.na(cosmic_sqsqsq)]
  cosmic_sqsqsq<-cosmic_sqsqsq[!is.na(cosmic_sqsqsq)]

  gnomad_aap<-gnomad_aap[!is.na(gnomad_sqsqsq)]
  gnomad_rfsq<-gnomad_rfsq[!is.na(gnomad_sqsqsq)]
  gnomad_ind<-gnomad_ind[!is.na(gnomad_sqsqsq)]
  gnomad_sqsqsq<-gnomad_sqsqsq[!is.na(gnomad_sqsqsq)]

  celegans_aap<-celegans_aap[!is.na(celegans_sqsqsq)]
  celegans_rfsq<-celegans_rfsq[!is.na(celegans_sqsqsq)]
  celegans_ind<-celegans_ind[!is.na(celegans_sqsqsq)]
  celegans_sqsqsq<-celegans_sqsqsq[!is.na(celegans_sqsqsq)]

  
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
  common_ind<-numeric(length(u4))
  
  sp1_common_aa_length<-0
  sp1_common_refseq_length<-0
  sp2_common_aa_length<-0
  sp2_common_refseq_length<-0
  common_ind_length<-0
  
  pb <- txtProgressBar(min = 0, max = max(commonind), style = 3)
  for (i in commonind){
    
    # Mouse - Celegans
    anumber<- Reduce(intersect, list(get(paste0(sp1, "_sqsqsq"))[which(get(paste0(sp1, "_ind")) %fin% i)], get(paste0(sp2, "_sqsqsq"))[which(get(paste0(sp2, "_ind")) %fin% i)]))
    
    if (length(anumber > 0)){
      
      a1<-get(paste0(sp1, "_sqsqsq"))[get(paste0(sp1, "_ind")) == i]
      b1<-get(paste0(sp1, "_aap"))[get(paste0(sp1, "_ind")) == i]
      c1<-get(paste0(sp1, "_rfsq"))[get(paste0(sp1, "_ind")) == i]
      
      a2<-get(paste0(sp2, "_sqsqsq"))[get(paste0(sp2, "_ind")) == i]
      b2<-get(paste0(sp2, "_aap"))[get(paste0(sp2, "_ind")) == i]
      c2<-get(paste0(sp2, "_rfsq"))[get(paste0(sp2, "_ind")) == i]
      
      
      xsp1_common_aa<-b1[fmatch(anumber, a1)]
      xsp1_common_refseq<-c1[fmatch(anumber, a1)]
      
      xsp2_common_aa<-b2[fmatch(anumber, a2)]
      xsp2_common_refseq<-c2[fmatch(anumber, a2)]
      
      
      xsp1_common_aa<-xsp1_common_aa[!is.na(xsp1_common_aa)]
      xsp1_common_refseq<-xsp1_common_refseq[!is.na(xsp1_common_refseq)]
      
      xsp2_common_aa<-xsp2_common_aa[!is.na(xsp2_common_aa)]
      xsp2_common_refseq<-xsp2_common_refseq[!is.na(xsp2_common_refseq)]
      
      sp1_common_aa[(sp1_common_aa_length+1) : (sp1_common_aa_length+length(xsp1_common_aa))]<-xsp1_common_aa
      sp1_common_refseq[(sp1_common_refseq_length+1) : (sp1_common_refseq_length+length(xsp1_common_refseq))]<-xsp1_common_refseq
      sp2_common_aa[(sp2_common_aa_length+1) : (sp2_common_aa_length+length(xsp2_common_aa))]<-xsp2_common_aa
      sp2_common_refseq[(sp2_common_refseq_length+1) : (sp2_common_refseq_length+length(xsp2_common_refseq))]<-xsp2_common_refseq
      common_ind[(sp2_common_refseq_length+1) : (sp2_common_refseq_length+length(xsp2_common_refseq))]<-rep(i,length(xsp2_common_refseq))
      
      sp1_common_aa_length<-sp1_common_aa_length + length(xsp1_common_aa)
      sp1_common_refseq_length<-sp1_common_refseq_length + length(xsp1_common_refseq)
      sp2_common_aa_length<-sp2_common_aa_length + length(xsp2_common_aa)
      sp2_common_refseq_length<-sp2_common_refseq_length + length(xsp2_common_refseq)
      
    }
    
    setTxtProgressBar(pb, i)
  }
  sp1_common_aa<-sp1_common_aa[sp1_common_aa != 0]
  sp1_common_refseq<-sp1_common_refseq[sp1_common_refseq != ""]
  sp2_common_aa<-sp2_common_aa[sp2_common_aa != 0]
  sp2_common_refseq<-sp2_common_refseq[sp2_common_refseq != ""]
  common_ind<-common_ind[common_ind != 0]
  
  lst<-list(sp1_common_aa, sp1_common_refseq, sp2_common_aa, sp2_common_refseq, common_ind)
  return(lst)
}

