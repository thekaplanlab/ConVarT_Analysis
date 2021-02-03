

library(data.table)
library(stringr)
library(pbapply)
library(dplyr)
library(seqinr)
library(tidyr)
library(geneName)
library(svMisc)
library(readxl)
library(gdata)
library(stringi)
library(fastmatch)


source("convart_helper.R")
source("ens_np_blast.R")
ensttoensp<-fread("enst-ensp.txt")
ensttoensp<-ensttoensp[!ensttoensp$`Protein stable ID` == "",]
enstogname<-fread("ens_to_gname.txt")
nptogname<-fread("np_to_gname_human.txt")
colnames(nptogname)<-c("Gene_name","Refseq_ID")

# Triple msa


# Mouse double msa

table<-fread("msa_table.txt") %>%
  filter(grepl('musculus', fasta)) %>%
  filter(grepl('sapiens', fasta))

table2<-strsplit(as.character(table$fasta), "\n>")
table2<-pblapply(X = table2, FUN = function(t) gsub(pattern = "\n", replacement = "", x = t))

table233<-pblapply(table2, function(x) length(x) == 2)
table23_mouse<-table2[unlist(table233)]


# C. elegans double msa

table<-fread("msa_table.txt") %>%
  filter(grepl('elegans', fasta)) %>%
  filter(grepl('sapiens', fasta))

table2<-strsplit(as.character(table$fasta), "\n>")
table2<-pblapply(X = table2, FUN = function(t) gsub(pattern = "\n", replacement = "", x = t))

table233<-pblapply(table2, function(x) length(x) == 2)
table23_celegans<-table2[unlist(table233)]


# celegans

celegans<-fread("c_elegans_variations.csv")
celegans1<-str_split_fixed(celegans$V7, " to ", 2)
colnames(celegans1)<-c("from", "to")
celegans<-cbind(celegans, celegans1) %>%
  dplyr::select(V3, V6, V9, V10, V11, from, to) %>%
  unique()
colnames(celegans)[1:5]<-c("Id","aapos","Refseq_ID","Gene_ID","Gene_name")
celegans_phenotype<-fread("c_elegans_phenotypic_variations.csv")
celegans$type<-"Unknown"
celegans$type[which(celegans$Id %in% celegans_phenotype$V2)]<-"Phenotypic"

celegans$to[celegans$to == "opal stop"]<-"*"
celegans$to[celegans$to == "amber stop"]<-"*"
celegans$to[celegans$to == "ochre stop"]<-"*"

celegans$from<-aaa(celegans$from)
celegans$to<-aaa(celegans$to)

celegans<-celegans[!is.na(celegans$from)]
celegans<-celegans[!is.na(celegans$to)]

celegans$from<-a(celegans$from)
celegans$to<-a(celegans$to)

celegans$Phenotype<-celegans_phenotype$V3[match(celegans$Id, celegans_phenotype$V2)]
celegans<-celegans[,-4]



# mouse

mouse<-fread("mouse_variations.txt")
mouse1<-str_split_fixed(mouse$`Amino Acid Change`, "->", 2)
mouse$from<-mouse1[,1]
mouse$to<-mouse1[,2]
colnames(mouse)[26:27]<-c("from", "to")
mouse$to[mouse$to == "Stop"]<-"*"
#mouse<-unique(setDT(mouse), by= c("aa_position","Refseq_ID"))
mouse<-mouse[,c(1,12,15,23,24,26,27)]
colnames(mouse)<-c("Id","Gene_name","Phenotype","aapos","CCDS_ID","from","to")
mouse<-mouse %>% 
  mutate(CCDS_ID = strsplit(as.character(CCDS_ID), ",")) %>%
  unnest(CCDS_ID)
tempmouse<-str_split_fixed(mouse$CCDS_ID, "\\.", 2)
mouse$CCDS_ID<-tempmouse[,1]

mouse<-merge(mouse, mouseall, by = "CCDS_ID", all = TRUE)
mouse<-mouse[!is.na(mouse$aapos),]
mouse<-unique(mouse)
mouse$aapos<-as.integer(mouse$aapos)
mouse<-mouse[mouse$to != "",]
mouse<-mouse[mouse$from != "Disrupted splicing",]
mouse<-mouse[mouse$from != "Stop",]

mouse$type<-"phenotypic"
mouse<-mouse[mouse$to != "",]
mouse$Variation<-paste0(mouse$from, mouse$aapos, mouse$to)
mouse<-mouse[!is.na(mouse$Refseq_ID),]


# mutagen

mutagen<-fread("mutagenetix_incidental_mutations_20190904.txt") %>%
  dplyr::select(gene_symbol, aa_from, aa_to, aa_pos, fasta_header, sequence, predicted_effect, mutation_type, coordinate)
ncbimouseseq<-read.fasta(file = "GCF_000001635.26_GRCm38.p6_protein.faa", seqtype = "AA", as.string = TRUE, set.attributes = TRUE) %>%
  unlist() %>%
  as.data.frame()
ncbimouseseq$refseq<-rownames(ncbimouseseq)
mutagen$Refseq_ID<-ncbimouseseq$refseq[match(mutagen$sequence, ncbimouseseq$.)]
mutagen<-separate(mutagen, "Refseq_ID", "Refseq_ID", sep = "\\.")
colnames(mutagen)[2:4]<-c("from","to","aapos")
mutagen<-mutagen[!mutagen$from == "",]
mutagen<-mutagen[!is.na(mutagen$Refseq_ID),]
mutagen$type<-"incidental"
mutagen<-mutagen[mutagen$mutation_type != "makesense"]

mutagenphe<-fread("mutagenetix_phenotypic_mutations_20190909.txt")%>%
  dplyr::select(gene_symbol, aa_from, aa_to, aa_pos, fasta_header, sequence, predicted_effect, coordinate, mutation_type)
mutagenphe$Refseq_ID<-ncbimouseseq$refseq[match(mutagenphe$sequence, ncbimouseseq$.)]
mutagenphe<-separate(mutagenphe, "Refseq_ID", "Refseq_ID", sep = "\\.")
colnames(mutagenphe)[2:4]<-c("from","to","aapos")
mutagenphe<-mutagenphe[!mutagenphe$from == "",]
mutagenphe<-mutagenphe[!is.na(mutagenphe$Refseq_ID),]
mutagenphe<-mutagenphe[!is.na(mutagenphe$aapos),]
mutagenphe$type<-"phenotypic"
mutagenphe<-mutagenphe[mutagenphe$mutation_type !="makesense"]

mutagenphe2<-fread("mutagenetix_phenotypic_mutations_20201024.txt")
mutagenphe$Phenotype<-mutagenphe2$phenotypes[match(mutagenphe$coordinate, mutagenphe2$coordinate)]
mutagen$Phenotype<-NA


mutagen<-rbind(mutagen, mutagenphe)
mutagen<-dplyr::select(mutagen, gene_symbol, Refseq_ID, aapos, from, to, type, predicted_effect, Phenotype, coordinate)
mutagen<-unique(mutagen)
colnames(mutagen)[c(1,3,9)]<-c("Gene_name","aapos","Id")
mutagen<-mutagen[,-7]
mutagen$Variation<-paste0(mutagen$from, mutagen$aapos, mutagen$to)


# gnomad

gnomad<-fread("gnomad.csv", select = c(2,3,4,5,7,9,10,11)) # ENST, var, RS_numb, aapos
gnomad$Allele_frequency<-gnomad$V10/gnomad$V11
gnomad1<-gnomad[which(gnomad$V9 != 0)]
gnomad1$V5<-gsub("p\\.", "", gnomad1$V5)
gnomad11<-str_split(gnomad1$V5, "[0-9]", simplify = TRUE)
gnomad11<-data.frame(gnomad11, stringsAsFactors = FALSE)
gnomad1$from<-gnomad11[,1]
gnomad1$to<-apply(gnomad11, 1, function(x) paste0(x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11]))
gnomad1$from[gnomad1$from == "Ter"]<-"Stp"
gnomad1$to[gnomad1$to == "Ter"]<-"Stp"
gnomad1$from<-a(gnomad1$from)
gnomad1$to<-a(gnomad1$to)
colnames(gnomad1)[c(3,5,6)]<-c("Refseq_ID", "Id", "aapos")
gnomad1<-gnomad1[!is.na(gnomad1$from)]
gnomad1<-gnomad1[!is.na(gnomad1$to)]
gnomad1<-gnomad1[gnomad1$from !="*"]


gnomad1$refseq<-ensttoensp$`Protein stable ID`[match(gnomad1$Refseq_ID, ensttoensp$`Transcript stable ID`)]
gnomad<-gnomad1[!is.na(gnomad1$refseq),]
gnomad$refseq2<-all1$refseq[match(gnomad$refseq, all1$ens)]
gnomad<-gnomad[!is.na(gnomad$refseq2),]
gnomad<-dplyr::select(gnomad, -Refseq_ID, -refseq)
colnames(gnomad)[11]<-"Refseq_ID"
gnomad<-unique(gnomad)
gnomad<-gnomad[!gnomad$aapos == 0,]
gnomad<-gnomad[!gnomad$from == "",]
gnomad$aapos<-as.integer(gnomad$aapos)
gnomad<-gnomad[!is.na(gnomad$from),]
gnomad<-dplyr::select(gnomad, -V5)
gnomad$Id[gnomad$Id == ""]<-paste0(seq(1,length(gnomad$Id[gnomad$Id == ""]),1), "grs-?")

gnomad<-unique(gnomad)
colnames(gnomad)[2]<-"Gene_ID"

gnomad$Gene_name<-enstogname$`Gene name`[match(gnomad$Gene_ID, enstogname$`Gene stable ID`)]
gnomad<-dplyr::select(gnomad, -V10, -V11)
gnomad<-gnomad[,3:9]
gnomad$type<-"unknown"
gnomad$Phenotype<-"Unknown"
gnomad$Source<-"gnomAD"



# cosmic

cosmic<-fread("CosmicMutantExport.csv", select = c(2, 3, 18, 20, 21, 29, 37)) # ENST, var, aapos

cosmic$V20<-gsub("p\\.", "", cosmic$V20)
cosmic1<-str_split(cosmic$V20, "[0-9999]", simplify = TRUE)
cosmic1<-data.frame(cosmic1, stringsAsFactors = FALSE)
cosmic$from<-cosmic1[,1]
cosmic$to<-apply(cosmic1, 1, function(x) paste0(x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13]))
cosmic<-dplyr::select(cosmic, -V20)
colnames(cosmic)[1:6]<-c("Gene_name", "Refseq_ID", "Id", "mut", "type", "aapos")


cosmic$refseq<-ensttoensp$`Protein stable ID`[match(cosmic$Refseq_ID, ensttoensp$`Transcript stable ID`)]
cosmic<-cosmic[!is.na(cosmic$refseq),]
cosmic$refseq2<-all1$refseq[match(cosmic$refseq, all1$ens)]
cosmic<-cosmic[!is.na(cosmic$refseq2),]
cosmic<-dplyr::select(cosmic, -Refseq_ID, -refseq)
colnames(cosmic)[8]<-"Refseq_ID"
cosmic<-cosmic[cosmic$mut != "Deletion - In frame",]
cosmic$from<-aaa(cosmic$from)
cosmic$to<-aaa(cosmic$to)
cosmic<-cosmic[!is.na(cosmic$from)]
cosmic<-cosmic[!is.na(cosmic$to)]
cosmic$from<-a(cosmic$from)
cosmic$to<-a(cosmic$to)
cosmic<-cosmic[cosmic$from != "*"]

cosmic<-cosmic[!cosmic$aapos == 0,]
cosmic<-cosmic[!cosmic$from == "",]
cosmic$aapos<-as.integer(cosmic$aapos)
cosmic_e1<-str_split(cosmic$Gene_name, "_", simplify = TRUE)
cosmic_e<-cosmic
cosmic_e$Gene_name<-cosmic_e1[,1]
cosmic<-cosmic_e
cosmic<-unique(cosmic)
cosmic$Phenotype<-"Unknown"
cosmic$Source<-"COSMIC"
cosmic$Allele_frequency<-NA


# clinvar

clinvar<-fread("clinvar.csv")#, select = c(4,5,10,16,17,18))
clinvar2<-fread("clinvar_2020_11_21.csv")
colnames(clinvar)[7]<-"variation_id"
clinvar2<-anti_join(clinvar2, clinvar, by = "variation_id")
colnames(clinvar2)<-colnames(clinvar)
clinvar<-rbind(clinvar, clinvar2)

clinvar<-clinvar[clinvar$V8 == "single nucleotide variant"]
clinvar$V16<-gsub("\\.", "", clinvar$V16)
clinvar$V16<-gsub("\\(p", "", clinvar$V16)
clinvar$V16<-gsub("\\)", "", clinvar$V16)

clinvar1<-str_split_fixed(clinvar$V16, "[0-9]", 7)
clinvar1<-data.frame(clinvar1, stringsAsFactors = FALSE)

clinvar1$to<-apply(clinvar1, 1, function(x) paste0(x[2],x[3],x[4],x[5],x[6],x[7]))
clinvar$from<-clinvar1$X1
clinvar$to<-clinvar1$to
clinvar$from<-gsub("\\(", "", clinvar$from)
clinvar_stop<-clinvar[clinvar$to == "*"]
clinvar_synonymous<-clinvar[clinvar$to == "="]
clinvar<-anti_join(clinvar, clinvar_stop, by = "to")
clinvar<-anti_join(clinvar, clinvar_synonymous, by = "to")

clinvar$to<-a(clinvar$to)
clinvar<-clinvar[!is.na(clinvar$to),]
clinvar$from<-a(clinvar$from)
clinvar<-clinvar[!is.na(clinvar$from),]

clinvar_stop$from[clinvar_stop$from == "Glu"]<-"E"
clinvar_stop$from[clinvar_stop$from == "Gln"]<-"Q"
clinvar_synonymous$from<-a(clinvar_synonymous$from)
clinvar_synonymous<-clinvar_synonymous[!is.na(clinvar_synonymous$from),]
clinvar_synonymous$to<-clinvar_synonymous$from
clinvar<-rbind(clinvar, clinvar_stop, clinvar_synonymous)

colnames(clinvar)[17]<-"aapos"
colnames(clinvar)[18]<-"Refseq_ID"
clinvar<-dplyr::select(clinvar, V4, V5, variation_id, aapos, Refseq_ID, from, to, V10, V14)
clinvar<-unique(clinvar)
clinvar$aapos<-as.integer(clinvar$aapos)

clinvar$V5[clinvar$V5 == "rs-1"]<-paste0(seq(1,length(clinvar$V5[clinvar$V5 == "rs-1"]),1), "rs-?")
colnames(clinvar)[1:2]<-c("Gene_name","Id")
colnames(clinvar)[8:9]<-c("type","Phenotype")
colnames(all_clinvar_frequency)[1]<-"variation_id"
clinvarxy<-merge(clinvar, all_clinvar_frequency, by = "variation_id", all = TRUE)
clinvarxy<-clinvarxy[!is.na(clinvarxy$aapos),]
colnames(clinvarxy)[10]<-"Allele_frequency"
clinvarxy$Source<-"ClinVar"



# dbsnp

dbsnp<-fread("dbsnp.csv")
dbsnp1<-str_split(dbsnp$V12, "\\.", simplify = TRUE)
dbsnp$Refseq_ID<-dbsnp1[,1]

dbsnp2<-str_split(dbsnp1[,3], "[0-9]", simplify = TRUE)
dbsnp2<-data.frame(dbsnp2, stringsAsFactors = FALSE)
dbsnp$to<-apply(dbsnp2, 1, function(x) paste0(x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11]))
dbsnp$from<-dbsnp2$X1
colnames(dbsnp)[11]<-"aapos"

dbsnp$from[dbsnp$from == "Ter"]<-"Stp"
dbsnp$to[dbsnp$to == "Ter"]<-"Stp"
dbsnp$to[dbsnp$to == "%D"]<-dbsnp$from[dbsnp$to == "%D"]

dbsnp$to<-a(dbsnp$to)
dbsnp$from<-a(dbsnp$from)
dbsnp<-dbsnp[!is.na(dbsnp$from),]
dbsnp<-dbsnp[!is.na(dbsnp$to),]
dbsnp<-dbsnp[dbsnp$from != "*"]

dbsnp<-dplyr::select(dbsnp, V2, V5, aapos, Refseq_ID, from, to, V13)

dbsnp$refseq<-all1$refseq[match(dbsnp$Refseq_ID, all1$ens)]
dbsnp<-dbsnp[!is.na(dbsnp$refseq),]
colnames(dbsnp)[c(4,8)]<-c("ens","Refseq_ID")
dbsnp<-unique(dbsnp)
dbsnp$aapos<-as.integer(dbsnp$aapos)
dbsnp<-dbsnp[!is.na(dbsnp$aapos),]
colnames(dbsnp)[c(1,2,7)]<-c("Id","Gene_ID","type")

dbsnp$Gene_name<-enstogname$`Gene name`[match(dbsnp$Gene_ID, enstogname$`Gene stable ID`)]
dbsnp<-dbsnp[,c(1,3,5,6,7,8,9)]
dbsnp$Phenotype<-"Unknown"
dbsnp$Allele_frequency<-NA
dbsnp$Source<-"dbSNP"


