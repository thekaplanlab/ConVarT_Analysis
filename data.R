

# !NOT RUN! #

library(data.table)
library(stringr)
library(dplyr)
library(seqinr)
library(biomaRt)

source("ens_np_blast.R")


# msa

table<-fread("msa_table.txt")
table33<-table[grep("musculus", table$fasta, ignore.case = TRUE),]
table33<-table[grep("sapiens", table33$fasta, ignore.case = TRUE),]
table33<-table[grep("elegans", table33$fasta, ignore.case = TRUE),]

table2<-strsplit(as.character(table33$fasta), "\n>")

table2<-pblapply(X = table2, FUN = function(t) gsub(pattern = "\n", replacement = "", x = t))

table233<-pblapply(table2, function(x) length(x) == 3)
table23<-table2[unlist(table233)]
table23<-table23[1:length(table23)]


# mouse

mouse<-fread("mouse_variations.txt") %>% separate("Refseq_ID", "Refseq_ID", sep = "\\.")
mouse1<-str_split_fixed(mouse$`Amino Acid Change`, "->", 2)
mouse$from<-mouse1[,1]
mouse$to<-mouse1[,2]
colnames(mouse)[26:27]<-c("from", "to")


# celegans

celegans<-fread("c_elegans_variations.csv")
celegans1<-str_split_fixed(celegans$V7, " to ", 2)
colnames(celegans1)<-c("from", "to")
celegans<-cbind(celegans, celegans1)


# clinvar

clinvar<-fread("clinvar.csv")
clinvar$V16<-gsub("\\.", "", clinvar$V16)
clinvar$V16<-gsub("\\(p", "", clinvar$V16)
clinvar$V16<-gsub("\\)", "", clinvar$V16)

clinvar1<-str_split_fixed(clinvar$V16, "[0-9]", 4)
clinvar1<-data.frame(clinvar1, stringsAsFactors = FALSE)

clinvar1$to<-apply(clinvar1, 1, function(x) paste0(x[2],x[3],x[4]))
clinvar$from<-clinvar1$X1
clinvar$to<-clinvar1$to

clinvar$from<-a(clinvar$from)
clinvar$to<-a(clinvar$to)
colnames(clinvar)[17]<-"aapos"
colnames(clinvar)[18]<-"Refseq_ID"


# dbsnp

dbsnp<-fread("dbsnp.csv")
dbsnp1<-str_split(dbsnp$V12, "\\.", simplify = TRUE)
dbsnp$Refseq_ID<-dbsnp1[,1]

dbsnp2<-str_split(dbsnp1[,3], "[0-9]", simplify = TRUE)
dbsnp2<-data.frame(dbsnp2, stringsAsFactors = FALSE)
dbsnp$to<-apply(dbsnp2, 1, function(x) paste0(x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11]))
dbsnp$from<-dbsnp2$X1
colnames(dbsnp)[11]<-"aapos"

dbsnp$to<-a(dbsnp$to)
dbsnp$from<-a(dbsnp$from)
dbsnp<-dplyr::select(dbsnp, V2, aapos, Refseq_ID, from, to)

dbsnp$refseq<-all1$refseq[match(dbsnp$Refseq_ID, all1$ens)]
dbsnp<-dbsnp[!is.na(dbsnp$refseq),]
colnames(dbsnp)[c(3,6)]<-c("ens","Refseq_ID")


# cosmic

cosmic<-fread("CosmicMutantExport.csv", select = c(3, 20, 37)) # ENST, var, aapos

cosmic$V20<-gsub("p\\.", "", cosmic$V20)
cosmic1<-str_split(cosmic$V20, "[0-9999]", simplify = TRUE)
cosmic1<-data.frame(cosmic1, stringsAsFactors = FALSE)
cosmic$from<-cosmic1[,1]
cosmic$to<-apply(cosmic1, 1, function(x) paste0(x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13]))
colnames(cosmic)[c(1,3)]<-c("Refseq_ID", "aapos")

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
glist<-getBM(filters = "ensembl_transcript_id", 
              attributes = c("ensembl_transcript_id","refseq_peptide"),
              values = unique(cosmic$Refseq_ID), mart = mart)
cosmic$refseq<-glist$refseq_peptide[match(cosmic$Refseq_ID, glist$ensembl_transcript_id)]
cosmic$refseq[which(cosmic$refseq == "")]<-NA
cosmic<-cosmic[!is.na(cosmic$refseq),]
colnames(cosmic)[c(1,6)]<-c("ens","Refseq_ID")


# gnomad

gnomad<-fread("gnomad.csv", select = c(4,5,7,9)) # ENST, var, RS_numb, aapos

gnomad1<-gnomad[which(gnomad$V9 != 0)]
gnomad1$V5<-gsub("p\\.", "", gnomad1$V5)
gnomad11<-str_split(gnomad1$V5, "[0-9]", simplify = TRUE)
gnomad11<-data.frame(gnomad11, stringsAsFactors = FALSE)
gnomad1$from<-gnomad11[,1]
gnomad1$to<-apply(gnomad11, 1, function(x) paste0(x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11]))
gnomad1$from<-a(gnomad1$from)
gnomad1$to<-a(gnomad1$to)
colnames(gnomad1)[c(1,3,4)]<-c("Refseq_ID", "RS_ID", "aapos")

glist<-getBM(filters = "ensembl_transcript_id", 
             attributes = c("ensembl_transcript_id","refseq_peptide"),
             values = unique(gnomad1$Refseq_ID), mart = mart)

gnomad1$refseq<-glist$refseq_peptide[match(gnomad1$Refseq_ID, glist$ensembl_transcript_id)]
gnomad1$refseq[which(gnomad1$refseq == "")]<-NA
gnomad<-gnomad1[!is.na(gnomad1$refseq),]
colnames(gnomad)[c(1,7)]<-c("ens","Refseq_ID")


# mutagen

mutagen<-fread("mutagenetix_incidental_mutations_20190904.txt") %>%
  dplyr::select(aa_from, aa_to, aa_pos, fasta_header, sequence)
ncbimouseseq<-read.fasta(file = "GCF_000001635.26_GRCm38.p6_protein.faa", seqtype = "AA", as.string = TRUE, set.attributes = TRUE) %>%
  unlist() %>%
  as.data.frame()
ncbimouseseq$refseq<-rownames(ncbimouseseq)
mutagen$Refseq_ID<-ncbimouseseq$refseq[match(mutagen$sequence, ncbimouseseq$.)]
mutagen<-separate(mutagen, "Refseq_ID", "Refseq_ID", sep = "\\.")
colnames(mutagen)[1:3]<-c("from","to","aa_position")
mutagen<-mutagen[-which(mutagen$from == ""),]
mutagen<-mutagen[!is.na(mutagen$Refseq_ID),]
mutagen$type<-"incidental"

mutagenphe<-fread("mutagenetix_phenotypic_mutations_20190909.txt")%>%
  dplyr::select(aa_from, aa_to, aa_pos, fasta_header, sequence)
mutagenphe$Refseq_ID<-ncbimouseseq$refseq[match(mutagenphe$sequence, ncbimouseseq$.)]
mutagenphe<-separate(mutagenphe, "Refseq_ID", "Refseq_ID", sep = "\\.")
colnames(mutagenphe)[1:3]<-c("from","to","aa_position")
mutagenphe<-mutagenphe[-which(mutagenphe$from == ""),]
mutagenphe<-mutagenphe[!is.na(mutagenphe$Refseq_ID),]
mutagenphe$type<-"phenotypic"

mutagen<-rbind(mutagen, mutagenphe)

tempmutagen<-str_split_fixed(mutagen$Refseq_ID , "\\.", 2)
mutagen$Refseq_ID<-tempmutagen[,1]

