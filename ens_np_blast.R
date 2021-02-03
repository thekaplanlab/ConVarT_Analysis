
# !NOT RUN! #
# This script generates the tables required for protein id conversions #

library(seqinr)
library(dplyr)
library(data.table)
library(stringr)

enstogname<-fread("ens_to_gname.txt")

nptogname<-fread("np_to_gname_human.txt")
colnames(nptogname)<-c("Gene_name","Refseq_ID")

ensttoensp<-fread("enst-ensp.txt")
ensttoensp<-ensttoensp[!ensttoensp$`Protein stable ID` == "",]

seq1<-read.fasta(file = "CCDS_protein.current.faa", seqtype = "AA", as.string = TRUE, set.attributes = TRUE) %>%
  unlist() %>%
  as.data.frame()
seq1$ens<-rownames(seq1)
ss<-str_split(seq1$ens, "\\|", simplify = TRUE)
ss<-data.frame(ss, stringsAsFactors = FALSE)
ss1<-ss[,1]
ss1<-str_split(ss1, "\\.", simplify = TRUE)
seq1$ens<-ss1[,1]


ccdsids<-fread("CCDS.current.txt")
ccdsids<-ccdsids[ccdsids$match_type == "Identical"]
ccdsids1<-str_split(ccdsids$ccds_id, "\\.", simplify = TRUE)
ccdsids$ccds_id<-ccdsids1[,1]


mouseseq<-read.fasta(file = "GCF_000001635.27_GRCm39_protein.faa", seqtype = "AA", as.string = TRUE, set.attributes = TRUE) %>%
  unlist() %>%
  as.data.frame()
mouseseq$ens<-rownames(mouseseq)
tempmouseseq<-str_split_fixed(mouseseq$ens, "\\.", 2)
mouseseq$ens<-tempmouseseq[,1]

mouseall<-merge(mouseseq, seq1, by = ".")
mouseall<-mouseall[,2:3]
colnames(mouseall)<-c("Refseq_ID","CCDS_ID")


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

all1_temp<-all1
colnames(all1_temp)[2:3]<-c("RefSeq Protein ID","Protein stable ID")
all3<-merge(all1_temp[,2:3], ensttoensp[,2:3], by = "Protein stable ID", all = TRUE)
all3<-all3[!is.na(all3$`RefSeq Protein ID`),] %>% unique()

write.table(all3, "est_esp_np_table.txt", row.names = FALSE, quote = FALSE, sep = "\t")
