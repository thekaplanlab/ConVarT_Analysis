library(data.table)
library(tidyr)
library(splitstackshape)
library(pbapply)
library(stringr)
library(svMisc)

source("convart_helper.R")
source("data.R")
source("analysis_function.R")


Sys.time()

mouse_aap<-c()
mouse_rfsq<-c()
mouse_sqsqsq<-c()
mouse_ind<-c()
mutagen_aap<-c()
mutagen_rfsq<-c()
mutagen_sqsqsq<-c()
mutagen_ind<-c()
clinvar_aap<-c()
clinvar_rfsq<-c()
clinvar_sqsqsq<-c()
clinvar_ind<-c()
dbsnp_aap<-c()
dbsnp_rfsq<-c()
dbsnp_sqsqsq<-c()
dbsnp_ind<-c()
cosmic_aap<-c()
cosmic_rfsq<-c()
cosmic_sqsqsq<-c()
cosmic_ind<-c()
gnomad_aap<-c()
gnomad_rfsq<-c()
gnomad_sqsqsq<-c()
gnomad_ind<-c()
celegans_aap<-c()
celegans_rfsq<-c()
celegans_sqsqsq<-c()
celegans_ind<-c()
pb = txtProgressBar(min = 0, max = length(table23), initial = 0)
for (i in 1:length(table23)){
    mouse_an<-0
    mouse_anumber<-c()
    mouse_indx<-c()
    mouse_a<-c(unlist(table23[[i]][2]))
    mouse_a<-c(unlist(strsplit(mouse_a, "]")))
    mouse_ag<-c(unlist(strsplit(mouse_a, " \\["))[[1]][1])
    mouse_ag<-c(unlist(strsplit(mouse_ag, "\\."))[[1]][1])
    
    mouse_a1<-c(unlist(strsplit(mouse_a[[2]][1], "")))         #split msa into individual aa
    for (j in 1:length(mouse_a1)){
      if (mouse_a1[j] != "-"){
        mouse_an<-mouse_an+1
        mouse_anumber[j]<-mouse_an
      }
      else {mouse_anumber[j]<-NA}                        # aa number of individual aa
    }
    
    mouse_aapos<-c(mouse$aa_position[which(mouse$Refseq_ID == mouse_ag)])
    mouse_sssqq<-which(mouse_anumber %fin% mouse_aapos)
    mouse_sqsqsq<-c(mouse_sqsqsq, mouse_sssqq)
    mouse_rrffssq<-mouse$Refseq_ID[which(mouse$Refseq_ID == mouse_ag)]
    mouse_rfsq<-c(mouse_rfsq, mouse_rrffssq)
    mouse_aap<-c(mouse_aap, mouse_aapos)
    
    if (length(mouse_aapos) != 0){
      mouse_indx <- rep(i, length(mouse_aapos))
      mouse_ind <- c(mouse_ind, mouse_indx)
    }

    
    mutagen_aapos<-c(mutagen$aa_position[which(mutagen$Refseq_ID == mouse_ag)])
    mutagen_sssqq<-which(mouse_anumber %fin% mutagen_aapos)
    mutagen_sqsqsq<-c(mutagen_sqsqsq, mutagen_sssqq)
    mutagen_rrffssq<-mutagen$Refseq_ID[which(mutagen$Refseq_ID == mouse_ag)]
    mutagen_rfsq<-c(mutagen_rfsq, mutagen_rrffssq)
    mutagen_aap<-c(mutagen_aap, mutagen_aapos)
    
    if (length(mutagen_aapos) != 0){
      mutagen_indx <- rep(i, length(mutagen_aapos))
      mutagen_ind <- c(mutagen_ind, mutagen_indx)
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
     for (j in 1:length(human_a1)){
       if (human_a1[j] != "-"){
         human_an<-human_an+1
         human_anumber[j]<-human_an
       }
       else {human_anumber[j]<-NA}                        # aa number of individual aa
     }
     
     clinvar_aapos<-c(clinvar$aapos[which(clinvar$Refseq_ID == human_ag)])
     clinvar_sssqq<-which(human_anumber %fin% clinvar_aapos)
     clinvar_sqsqsq<-c(clinvar_sqsqsq, clinvar_sssqq)
     clinvar_rrffssq<-clinvar$Refseq_ID[which(clinvar$Refseq_ID == human_ag)]
     clinvar_rfsq<-c(clinvar_rfsq, clinvar_rrffssq)
     clinvar_aap<-c(clinvar_aap, clinvar_aapos)
     
     dbsnp_aapos<-c(dbsnp$aapos[which(dbsnp$Refseq_ID == human_ag)])
     dbsnp_sssqq<-which(human_anumber %fin% dbsnp_aapos)
     dbsnp_sqsqsq<-c(dbsnp_sqsqsq, dbsnp_sssqq)
     dbsnp_rrffssq<-dbsnp$Refseq_ID[which(dbsnp$Refseq_ID == human_ag)]
     dbsnp_rfsq<-c(dbsnp_rfsq, dbsnp_rrffssq)
     dbsnp_aap<-c(dbsnp_aap, dbsnp_aapos)
     
     cosmic_aapos<-c(cosmic$aapos[which(cosmic$Refseq_ID == human_ag)])
     cosmic_sssqq<-which(human_anumber %fin% cosmic_aapos)
     cosmic_sqsqsq<-c(cosmic_sqsqsq, cosmic_sssqq)
     cosmic_rrffssq<-cosmic$Refseq_ID[which(cosmic$Refseq_ID == human_ag)]
     cosmic_rfsq<-c(cosmic_rfsq, cosmic_rrffssq)
     cosmic_aap<-c(cosmic_aap, cosmic_aapos)
     
     gnomad_aapos<-c(gnomad$aapos[which(gnomad$Refseq_ID == human_ag)])
     gnomad_sssqq<-which(human_anumber %fin% gnomad_aapos)
     gnomad_sqsqsq<-c(gnomad_sqsqsq, gnomad_sssqq)
     gnomad_rrffssq<-gnomad$Refseq_ID[which(gnomad$Refseq_ID == human_ag)]
     gnomad_rfsq<-c(gnomad_rfsq, gnomad_rrffssq)
     gnomad_aap<-c(gnomad_aap, gnomad_aapos)
    
     if (length(clinvar_aapos) != 0){
       clinvar_indx <- rep(i, length(clinvar_aapos))
       clinvar_ind <- c(clinvar_ind, clinvar_indx)
     }
     
     
     if (length(dbsnp_aapos) != 0){
       dbsnp_indx <- rep(i, length(dbsnp_aapos))
       dbsnp_ind <- c(dbsnp_ind, dbsnp_indx)
     }
     
     
     if (length(cosmic_aapos) != 0){
       cosmic_indx <- rep(i, length(cosmic_aapos))
       cosmic_ind <- c(cosmic_ind, cosmic_indx)
     }
     
     
     if (length(gnomad_aapos) != 0){
       gnomad_indx <- rep(i, length(gnomad_aapos))
       gnomad_ind <- c(gnomad_ind, gnomad_indx)
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
    for (j in 1:length(celegans_a1)){
      if (celegans_a1[j] != "-"){
        celegans_an<-celegans_an+1
        celegans_anumber[j]<-celegans_an
      }
      else {celegans_anumber[j]<-NA}                        # aa number of individual aa
    }
    
    celegans_aapos<-c(celegans$V6[which(celegans$V9 == celegans_ag)])
    celegans_sssqq<-which(celegans_anumber %fin% celegans_aapos)
    celegans_sqsqsq<-c(celegans_sqsqsq, celegans_sssqq)
    celegans_rrffssq<-celegans$V9[which(celegans$V9 == celegans_ag)]
    celegans_rfsq<-c(celegans_rfsq, celegans_rrffssq)
    celegans_aap<-c(celegans_aap, celegans_aapos)
    
    if (length(celegans_aapos) != 0){
      celegans_indx <- rep(i, length(celegans_aapos))
      celegans_ind <- c(celegans_ind, celegans_indx)
    }
    
    progress(i, 411668, init = (value == 1))
    if (i == 411669) cat("Done!\n")
    #setTxtProgressBar(pb, i)
}
Sys.time()

# aapos = aa position of variations of one protein as given in the table
# sqsqsq = index of aa positions in msa
# rfsq = NP numbers of aminoacids of given variations
# aap = aa positions of variations
# ind = index of protein in msa


# Analysis

mousecelegans<-analysis("mouse", "celegans")
mutagencelegans<-analysis("mutagen", "celegans")

