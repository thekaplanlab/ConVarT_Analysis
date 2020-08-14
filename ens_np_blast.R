
library(seqinr)
library(dplyr)
library(data.table)
library(stringr)

  seq<-read.fasta(file = "Homo_sapiens.GRCh38.pep.all.fa", seqtype = "AA", as.string = TRUE, set.attributes = TRUE) %>%
    unlist() %>%
    as.data.frame()
  seq$ens<-rownames(seq)
  
  
  ncbiseq<-read.fasta(file = "GCF_000001405.39_GRCh38.p13_protein.faa", seqtype = "AA", as.string = TRUE, set.attributes = TRUE) %>%
    unlist() %>%
    as.data.frame()
  ncbiseq$refseq<-rownames(ncbiseq)
  
  
  all<-merge(ncbiseq, seq, by = ".", all = TRUE)
  all1<-na.omit(all)
  all2<-unique(setDT(all1), by= "ens")
  
  tempall<-str_split_fixed(all1$ens, "\\.", 2)
  all1$ens<-tempall[,1]
  tempall1<-str_split_fixed(all1$refseq, "\\.", 2)
  all1$refseq<-tempall1[,1]
