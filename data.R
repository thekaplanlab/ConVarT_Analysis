

library(data.table)
library(stringr)
library(dplyr)
library(seqinr)
library(biomaRt)


table<-fread("msa_table.txt")

table33<-table[grep("musculus", table$fasta, ignore.case = TRUE),]
table33<-table[grep("sapiens", table33$fasta, ignore.case = TRUE),]
table33<-table[grep("elegans", table33$fasta, ignore.case = TRUE),]

table2<-strsplit(as.character(table33$fasta), "\n>")
table2<-pblapply(X = table2, FUN = function(t) gsub(pattern = "\n", replacement = "", x = t))

table233<-pblapply(table2, function(x) length(x) == 3)
table23<-table2[unlist(table233)]
table23<-table23[1:length(table23)]



mouse<-fread("mouse_variations.txt") %>% separate("Refseq_ID", "Refseq_ID", sep = "\\.")
mouse1<-str_split_fixed(mouse$`Amino Acid Change`, "->", 2)
mouse$from<-mouse1[,1]
mouse$to<-mouse1[,2]
colnames(mouse)[26:27]<-c("from", "to")

celegans<-fread("c_elegans_variations.csv")
celegans1<-str_split_fixed(celegans$V7, " to ", 2)
colnames(celegans1)<-c("from", "to")
celegans<-cbind(celegans, celegans1)



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

# clinvar

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

dbsnpw<-dbsnp
dbsnpw$refseq<-all1$refseq[match(dbsnpw$Refseq_ID, all1$ens)]
dbsnpwNoNA<-dbsnpw[!is.na(dbsnpw$refseq)]
dbsnpwNA<-dbsnpw[is.na(dbsnpw$refseq)]
enswo<-unique(dbsnpwNA$Refseq_ID)
ensw<-unique(dbsnpwNoNA$Refseq_ID)
totalens<-unique(dbsnp$Refseq_ID)

glist21<-glist2[which(glist2$refseq_peptide != ""),]
dbsnpw2<-merge(dbsnpw, glist21, by= "Refseq_ID", all = TRUE)
dbsnpw$refseq2<-glist21$refseq_peptide[match(dbsnpw$Refseq_ID, glist21$Refseq_ID)]
dbsnpw<- dbsnpw %>% mutate(refseq1= coalesce(refseq,refseq2))

gene2ensembl<-fread("gene2ensembl", select = c(6,7))
tempens<-str_split_fixed(gene2ensembl$Ensembl_protein_identifier, "\\.", 2)
gene2ensembl$Ensembl_protein_identifier<-tempens[,1]
gene2ensembl<-gene2ensembl[which(gene2ensembl$Ensembl_protein_identifier %in% totalens)]


lsts<-gconvert(c(totalens), organism = "hsapiens", target = "REFSEQ_PEPTIDE")

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
listAttributes(ensembl)
glist3<-getBM(filters = "ensembl_transcript_id", 
      attributes = c("ensembl_transcript_id","refseq_peptide"),
      values = unique(cosmic$Refseq_ID) , mart = mart)
colnames(glist2)[1]<-"Refseq_ID"
dbsnp_f<-merge(dbsnp, glist2, by = "Refseq_ID", allow.cartesian = TRUE, all = TRUE)
dbsnp_f<-dplyr::select(dbsnp_f, Refseq_ID, V2, aapos, V13, to, from, refseq_peptide)
dbsnp_fu<-unique(dbsnp_f)
dbsnp_funa<-dbsnp_fu[-which(dbsnp_fu$refseq_peptide == "")]

# dbsnp_f


cosmic<-fread("CosmicMutantExport.csv", select = c(3, 20, 37)) # ENST, var, aapos

cosmic$V20<-gsub("p\\.", "", cosmic$V20)
cosmic1<-str_split(cosmic$V20, "[0-9999]", simplify = TRUE)
cosmic1<-data.frame(cosmic1, stringsAsFactors = FALSE)
cosmic$from<-cosmic1[,1]
cosmic$to<-apply(cosmic1, 1, function(x) paste0(x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13]))
colnames(cosmic)[c(1,3)]<-c("Refseq_ID", "aapos")

# cosmic 

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

glist1<-getBM(filters = "ensembl_transcript_id", 
              attributes = c("ensembl_transcript_id","refseq_peptide"),
              values = unique(gnomad1$Refseq_ID), mart = mart)
colnames(glist_gnomad)[1]<-"Refseq_ID"
gnomad_f<-merge(gnomad1, glist1_gnomad, by = "Refseq_ID", allow.cartesian = TRUE, all = TRUE)



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
