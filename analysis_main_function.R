
# !NOT RUN! #

# aapos = aa position of variations of one protein as given in the table
# sqsqsq = index of aa positions in msa
# rfsq = NP numbers of aminoacids of given variations
# aap = aa positions of variations
# ind = index of protein in msa

analysis_main<-function(...){
  
  mouse_aap<-numeric(300000)
  mouse_rfsq<-character(300000)
  mouse_sqsqsq<-numeric(300000)
  mouse_ind<-numeric(300000)
  mutagen_aap<-numeric(800000)
  mutagen_rfsq<-character(800000)
  mutagen_sqsqsq<-numeric(800000)
  mutagen_ind<-numeric(800000)
  clinvar_aap<-numeric(8000000)
  clinvar_rfsq<-character(8000000)
  clinvar_sqsqsq<-numeric(8000000)
  clinvar_ind<-numeric(8000000)
  dbsnp_aap<-numeric(8000000)
  dbsnp_rfsq<-character(8000000)
  dbsnp_sqsqsq<-numeric(8000000)
  dbsnp_ind<-numeric(8000000)
  cosmic_aap<-numeric(66000000)
  cosmic_rfsq<-character(66000000)
  cosmic_sqsqsq<-numeric(66000000)
  cosmic_ind<-numeric(66000000)
  gnomad_aap<-numeric(80000000)
  gnomad_rfsq<-character(80000000)
  gnomad_sqsqsq<-numeric(80000000)
  gnomad_ind<-numeric(80000000)
  celegans_aap<-numeric(11000000)
  celegans_rfsq<-character(11000000)
  celegans_sqsqsq<-numeric(11000000)
  celegans_ind<-numeric(11000000)
  
  mouse_aap_length<-0
  mouse_rfsq_length<-0
  mouse_sqsqsq_length<-0
  mouse_ind_length<-0
  mutagen_aap_length<-0
  mutagen_rfsq_length<-0
  mutagen_sqsqsq_length<-0
  mutagen_ind_length<-0
  clinvar_aap_length<-0
  clinvar_rfsq_length<-0
  clinvar_sqsqsq_length<-0
  clinvar_ind_length<-0
  dbsnp_aap_length<-0
  dbsnp_rfsq_length<-0
  dbsnp_sqsqsq_length<-0
  dbsnp_ind_length<-0
  cosmic_aap_length<-0
  cosmic_rfsq_length<-0
  cosmic_sqsqsq_length<-0
  cosmic_ind_length<-0
  gnomad_aap_length<-0
  gnomad_rfsq_length<-0
  gnomad_sqsqsq_length<-0
  gnomad_ind_length<-0
  celegans_aap_length<-0
  celegans_rfsq_length<-0
  celegans_sqsqsq_length<-0
  celegans_ind_length<-0
  
  for (i in 1:length(table23)){
    mouse_an<-0
    mouse_anumber<-c()
    mouse_indx<-c()
    mouse_a<-c(unlist(table23[[i]][2]))
    mouse_a<-c(unlist(strsplit(mouse_a, "]")))
    mouse_ag<-c(unlist(strsplit(mouse_a, " \\["))[[1]][1])
    mouse_ag<-c(unlist(strsplit(mouse_ag, "\\."))[[1]][1])
    
    mouse_a1<-c(unlist(strsplit(mouse_a[[2]][1], "")))         #split msa into individual aa
    
    mouse_a12<-mouse_a1
    mouse_a12[mouse_a12 != "-"]<-1
    mouse_a12[mouse_a12 == "-"]<-0
    mouse_a12<-as.numeric(mouse_a12)
    mouse_an<-cumsum(c(mouse_a12 != 0))
    mouse_anumber<-mouse_an*mouse_a12
    
    mouse_aapos<-c(mouse$aa_position[which(mouse$Refseq_ID == mouse_ag)])
    mouse_sssqq<-fmatch(mouse_aapos, mouse_anumber)

    if(length(mouse_aapos)>0){
    mouse_sqsqsq[(mouse_sqsqsq_length+1) : (mouse_sqsqsq_length+length(mouse_sssqq))]<-mouse_sssqq
    mouse_sqsqsq_length<-mouse_sqsqsq_length + length(mouse_sssqq)
    
    mouse_rrffssq<-mouse$Refseq_ID[which(mouse$Refseq_ID == mouse_ag)]
    mouse_rfsq[(mouse_rfsq_length+1) : (mouse_rfsq_length+length(mouse_rrffssq))]<-mouse_rrffssq
    mouse_rfsq_length<-mouse_rfsq_length + length(mouse_rrffssq)

    mouse_aap[(mouse_aap_length+1) : (mouse_aap_length+length(mouse_aapos))]<-mouse_aapos
    mouse_aap_length<-mouse_aap_length + length(mouse_aapos)
    }
    
    if (length(mouse_aapos) != 0){
      mouse_indx <- rep(i, length(mouse_aapos))
      mouse_ind[(mouse_ind_length+1) : (mouse_ind_length+length(mouse_indx))]<-mouse_indx
      mouse_ind_length<-mouse_ind_length + length(mouse_indx)
    }
    
    
    mutagen_aapos<-c(mutagen$aa_position[which(mutagen$Refseq_ID == mouse_ag)])
    mutagen_sssqq<-fmatch(mutagen_aapos, mouse_anumber)
    
    if(length(mutagen_aapos>0)){
    mutagen_sqsqsq[(mutagen_sqsqsq_length+1) : (mutagen_sqsqsq_length+length(mutagen_sssqq))]<-mutagen_sssqq
    mutagen_sqsqsq_length<-mutagen_sqsqsq_length + length(mutagen_sssqq)
    mutagen_rrffssq<-mutagen$Refseq_ID[which(mutagen$Refseq_ID == mouse_ag)]

    mutagen_rfsq[(mutagen_rfsq_length+1) : (mutagen_rfsq_length+length(mutagen_rrffssq))]<-mutagen_rrffssq
    mutagen_rfsq_length<-mutagen_rfsq_length + length(mutagen_rrffssq)

    mutagen_aap[(mutagen_aap_length+1) : (mutagen_aap_length+length(mutagen_aapos))]<-mutagen_aapos
    mutagen_aap_length<-mutagen_aap_length + length(mutagen_aapos)
    }
    
    if (length(mutagen_aapos) != 0){
      mutagen_indx <- rep(i, length(mutagen_aapos))
      mutagen_ind[(mutagen_ind_length+1) : (mutagen_ind_length+length(mutagen_indx))]<-mutagen_indx
      mutagen_ind_length<-mutagen_ind_length + length(mutagen_indx)
    }
    
    
    # Human # 
    
    human_an<-0
    human_anumber<-c()
    clinvar_indx<-c()
    dbsnp_indx<-c()
    cosmic_indx<-c()
    gnomad_indx<-c()
    human_a<-c(unlist(table23[[i]][1]))
    human_a<-c(unlist(strsplit(human_a, "]")))
    human_ag<-c(unlist(strsplit(human_a, " \\["))[[1]][1])
    human_ag<-c(unlist(strsplit(human_ag, "\\."))[[1]][1])
    human_ag<-gsub(">","",human_ag)
    
    human_a1<-c(unlist(strsplit(human_a[[2]][1], "")))         #split msa into individual aa
    
    
    human_a12<-human_a1
    human_a12[human_a12 != "-"]<-1
    human_a12[human_a12 == "-"]<-0
    human_an<-cumsum(c(as.numeric(human_a12 != 0)))
    human_a12<-as.numeric(human_a12)
    human_anumber<-human_an*human_a12
    
    
    clinvar_aapos<-c(clinvar$aapos[which(clinvar$Refseq_ID == human_ag)])
    clinvar_sssqq<-fmatch(clinvar_aapos, human_anumber)

    if(length(clinvar_aapos>0)){
  
    clinvar_sqsqsq[(clinvar_sqsqsq_length+1) : (clinvar_sqsqsq_length+length(clinvar_sssqq))]<-clinvar_sssqq
    clinvar_sqsqsq_length<-clinvar_sqsqsq_length + length(clinvar_sssqq)
    
    clinvar_rrffssq<-clinvar$Refseq_ID[which(clinvar$Refseq_ID == human_ag)]
    clinvar_rfsq[(clinvar_rfsq_length+1) : (clinvar_rfsq_length+length(clinvar_rrffssq))]<-clinvar_rrffssq
    clinvar_rfsq_length<-clinvar_rfsq_length + length(clinvar_rrffssq)

    clinvar_aap[(clinvar_aap_length+1) : (clinvar_aap_length+length(clinvar_aapos))]<-clinvar_aapos
    clinvar_aap_length<-clinvar_aap_length + length(clinvar_aapos)
    }
    
    dbsnp_aapos<-c(dbsnp$aapos[which(dbsnp$Refseq_ID == human_ag)])
    dbsnp_sssqq<-fmatch(dbsnp_aapos, human_anumber)

    if(length(dbsnp_aapos>0)){
    dbsnp_sqsqsq[(dbsnp_sqsqsq_length+1) : (dbsnp_sqsqsq_length+length(dbsnp_sssqq))]<-dbsnp_sssqq
    dbsnp_sqsqsq_length<-dbsnp_sqsqsq_length + length(dbsnp_sssqq)
    
    dbsnp_rrffssq<-dbsnp$Refseq_ID[which(dbsnp$Refseq_ID == human_ag)]
    dbsnp_rfsq[(dbsnp_rfsq_length+1) : (dbsnp_rfsq_length+length(dbsnp_rrffssq))]<-dbsnp_rrffssq
    dbsnp_rfsq_length<-dbsnp_rfsq_length + length(dbsnp_rrffssq)
    
    dbsnp_aap[(dbsnp_aap_length+1) : (dbsnp_aap_length+length(dbsnp_aapos))]<-dbsnp_aapos
    dbsnp_aap_length<-dbsnp_aap_length + length(dbsnp_aapos)
    }
    
    cosmic_aapos<-c(cosmic$aapos[which(cosmic$Refseq_ID == human_ag)])
    cosmic_sssqq<-fmatch(cosmic_aapos, human_anumber)
  
    if(length(cosmic_aapos>0)){
    cosmic_sqsqsq[(cosmic_sqsqsq_length+1) : (cosmic_sqsqsq_length+length(cosmic_sssqq))]<-cosmic_sssqq
    cosmic_sqsqsq_length<-cosmic_sqsqsq_length + length(cosmic_sssqq)
    
    cosmic_rrffssq<-cosmic$Refseq_ID[which(cosmic$Refseq_ID == human_ag)]
    cosmic_rfsq[(cosmic_rfsq_length+1) : (cosmic_rfsq_length+length(cosmic_rrffssq))]<-cosmic_rrffssq
    cosmic_rfsq_length<-cosmic_rfsq_length + length(cosmic_rrffssq)

    cosmic_aap[(cosmic_aap_length+1) : (cosmic_aap_length+length(cosmic_aapos))]<-cosmic_aapos
    cosmic_aap_length<-cosmic_aap_length + length(cosmic_aapos)
    }
    gnomad_aapos<-c(gnomad$aapos[which(gnomad$Refseq_ID == human_ag)])
    gnomad_sssqq<-fmatch(gnomad_aapos, human_anumber)

    if(length(gnomad_aapos>0)){
    gnomad_sqsqsq[(gnomad_sqsqsq_length+1) : (gnomad_sqsqsq_length+length(gnomad_sssqq))]<-gnomad_sssqq
    gnomad_sqsqsq_length<-gnomad_sqsqsq_length + length(gnomad_sssqq)
    
    gnomad_rrffssq<-gnomad$Refseq_ID[which(gnomad$Refseq_ID == human_ag)]
    gnomad_rfsq[(gnomad_rfsq_length+1) : (gnomad_rfsq_length+length(gnomad_rrffssq))]<-gnomad_rrffssq
    gnomad_rfsq_length<-gnomad_rfsq_length + length(gnomad_rrffssq)

    gnomad_aap[(gnomad_aap_length+1) : (gnomad_aap_length+length(gnomad_aapos))]<-gnomad_aapos
    gnomad_aap_length<-gnomad_aap_length + length(gnomad_aapos)
    }
    
    if (length(clinvar_aapos) != 0){
      clinvar_indx <- rep(i, length(clinvar_aapos))
      clinvar_ind[(clinvar_ind_length+1) : (clinvar_ind_length+length(clinvar_indx))]<-clinvar_indx
      clinvar_ind_length<-clinvar_ind_length + length(clinvar_indx)
    }
    
    if (length(dbsnp_aapos) != 0){
      dbsnp_indx <- rep(i, length(dbsnp_aapos))
      dbsnp_ind[(dbsnp_ind_length+1) : (dbsnp_ind_length+length(dbsnp_indx))]<-dbsnp_indx
      dbsnp_ind_length<-dbsnp_ind_length + length(dbsnp_indx)
    }
    
    if (length(cosmic_aapos) != 0){
      cosmic_indx <- rep(i, length(cosmic_aapos))
      cosmic_ind[(cosmic_ind_length+1) : (cosmic_ind_length+length(cosmic_indx))]<-cosmic_indx
      cosmic_ind_length<-cosmic_ind_length + length(cosmic_indx)
    }
    
    
    if (length(gnomad_aapos) != 0){
      gnomad_indx <- rep(i, length(gnomad_aapos))
      gnomad_ind[(gnomad_ind_length+1) : (gnomad_ind_length+length(gnomad_indx))]<-gnomad_indx
      gnomad_ind_length<-gnomad_ind_length + length(gnomad_indx)
    }
    
    
    # Celegans #
    
    celegans_an<-0
    celegans_anumber<-c()
    celegans_indx<-c()
    celegans_a<-c(unlist(table23[[i]][3]))
    celegans_a<-c(unlist(strsplit(celegans_a, "]")))
    celegans_ag<-c(unlist(strsplit(celegans_a, " \\["))[[1]][1])
    celegans_ag<-c(unlist(strsplit(celegans_ag, "\\."))[[1]][1])
    celegans_a1<-c(unlist(strsplit(celegans_a[[2]][1], "")))         #split msa into individual aa
    
    
    celegans_a12<-celegans_a1
    celegans_a12[celegans_a12 != "-"]<-1
    celegans_a12[celegans_a12 == "-"]<-0
    celegans_an<-cumsum(c(as.numeric(celegans_a12 != 0)))
    celegans_a12<-as.numeric(celegans_a12)
    celegans_anumber<-celegans_an*celegans_a12
    
    
    celegans_aapos<-c(celegans$V6[which(celegans$V9 == celegans_ag)])
    celegans_sssqq<-fmatch(celegans_aapos, celegans_anumber)

    if(length(celegans_aapos>0)){
    celegans_sqsqsq[(celegans_sqsqsq_length+1) : (celegans_sqsqsq_length+length(celegans_sssqq))]<-celegans_sssqq
    celegans_sqsqsq_length<-celegans_sqsqsq_length + length(celegans_sssqq)
    
    celegans_rrffssq<-celegans$V9[which(celegans$V9 == celegans_ag)]
    celegans_rfsq[(celegans_rfsq_length+1) : (celegans_rfsq_length+length(celegans_rrffssq))]<-celegans_rrffssq
    celegans_rfsq_length<-celegans_rfsq_length + length(celegans_rrffssq)

    celegans_aap[(celegans_aap_length+1) : (celegans_aap_length+length(celegans_aapos))]<-celegans_aapos
    celegans_aap_length<-celegans_aap_length + length(celegans_aapos)
    }
    
    if (length(celegans_aapos) != 0){
      celegans_indx <- rep(i, length(celegans_aapos))
      celegans_ind[(celegans_ind_length+1) : (celegans_ind_length+length(celegans_indx))]<-celegans_indx
      celegans_ind_length<-celegans_ind_length + length(celegans_indx)
    }
    
    progress(i, 411642, init = (value == 1))
    if (i == 411642) cat("Done!\n")
  }
  
  mouse_aap<-mouse_aap[1:mouse_aap_length]
  mouse_rfsq<-mouse_rfsq[1:mouse_rfsq_length]
  mouse_sqsqsq<-mouse_sqsqsq[1:mouse_sqsqsq_length]
  mouse_ind<-mouse_ind[1:mouse_ind_length]
  mutagen_aap<-mutagen_aap[1:mutagen_aap_length]
  mutagen_rfsq<-mutagen_rfsq[1:mutagen_rfsq_length]
  mutagen_sqsqsq<-mutagen_sqsqsq[1:mutagen_sqsqsq_length]
  mutagen_ind<-mutagen_ind[1:mutagen_ind_length]
  clinvar_aap<-clinvar_aap[1:clinvar_aap_length]
  clinvar_rfsq<-clinvar_rfsq[1:clinvar_rfsq_length]
  clinvar_sqsqsq<-clinvar_sqsqsq[1:clinvar_sqsqsq_length]
  clinvar_ind<-clinvar_ind[1:clinvar_ind_length]
  dbsnp_aap<-dbsnp_aap[1:dbsnp_aap_length]
  dbsnp_rfsq<-dbsnp_rfsq[1:dbsnp_rfsq_length]
  dbsnp_sqsqsq<-dbsnp_sqsqsq[1:dbsnp_sqsqsq_length]
  dbsnp_ind<-dbsnp_ind[1:dbsnp_ind_length]
  cosmic_aap<-cosmic_aap[1:cosmic_aap_length]
  cosmic_rfsq<-cosmic_rfsq[1:cosmic_rfsq_length]
  cosmic_sqsqsq<-cosmic_sqsqsq[1:cosmic_sqsqsq_length]
  cosmic_ind<-cosmic_ind[1:cosmic_ind_length]
  gnomad_aap<-gnomad_aap[1:gnomad_aap_length]
  gnomad_rfsq<-gnomad_rfsq[1:gnomad_rfsq_length]
  gnomad_sqsqsq<-gnomad_sqsqsq[1:gnomad_sqsqsq_length]
  gnomad_ind<-gnomad_ind[1:gnomad_ind_length]
  celegans_aap<-celegans_aap[1:celegans_aap_length]
  celegans_rfsq<-celegans_rfsq[1:celegans_rfsq_length]
  celegans_sqsqsq<-celegans_sqsqsq[1:celegans_sqsqsq_length]
  celegans_ind<-celegans_ind[1:celegans_ind_length]
  
  lstt<-list(mouse_aap,
             mouse_rfsq,
             mouse_sqsqsq,
             mouse_ind,
             mutagen_aap,
             mutagen_rfsq,
             mutagen_sqsqsq,
             mutagen_ind,
             clinvar_aap,
             clinvar_rfsq,
             clinvar_sqsqsq,
             clinvar_ind,
             dbsnp_aap,
             dbsnp_rfsq,
             dbsnp_sqsqsq,
             dbsnp_ind,
             cosmic_aap,
             cosmic_rfsq,
             cosmic_sqsqsq,
             cosmic_ind,
             gnomad_aap,
             gnomad_rfsq,
             gnomad_sqsqsq,
             gnomad_ind,
             celegans_aap,
             celegans_rfsq,
             celegans_sqsqsq,
             celegans_ind,
             mouse_aap_length,
             mouse_rfsq_length,
             mouse_sqsqsq_length,
             mouse_ind_length,
             mutagen_aap_length,
             mutagen_rfsq_length,
             mutagen_sqsqsq_length,
             mutagen_ind_length,
             clinvar_aap_length,
             clinvar_rfsq_length,
             clinvar_sqsqsq_length,
             clinvar_ind_length,
             dbsnp_aap_length,
             dbsnp_rfsq_length,
             dbsnp_sqsqsq_length,
             dbsnp_ind_length,
             cosmic_aap_length,
             cosmic_rfsq_length,
             cosmic_sqsqsq_length,
             cosmic_ind_length,
             gnomad_aap_length,
             gnomad_rfsq_length,
             gnomad_sqsqsq_length,
             gnomad_ind_length,
             celegans_aap_length,
             celegans_rfsq_length,
             celegans_sqsqsq_length,
             celegans_ind_length)
  return(lstt)
}
