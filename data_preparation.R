
somatic<-read.vcfR("homo_sapiens_somatic_incl_consequences1.vcf")
ch1<-read.vcfR("homo_sapiens-chr1x.vcf")
ch1<-fread("ch1_30.txt")
smtc<-fread("homo_sapiens-chr1x.vcf")
write.vcf(somatic, "chr1_1m.vcf")
colnames(output)<-as.character(output[1,])
output2<-output[output$Consequence == "coding_sequence_variant",]
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(readxl)
library(geneName)
library(pbapply)
library(stringi)

source("msa_index.R")
source("convart_helper.R")
source("ens_np_blast.R")

# Read ClinVar frequency data
clinvar_frequency<-fread("clinvar_frequency.tsv")

# Read amino acid conservation data
aaconservation<-read_xlsx("aa_conservation.xlsx")
aaconservation$to[is.na(aaconservation$to)]<-"B"

# Get all necessary frequency data from ClinVar frequency data
clinvar_frequency_topmed<-clinvar_frequency[grep("topmed", clinvar_frequency$frequ, ignore.case = TRUE),]
clinvar_frequency_topmed<-clinvar_frequency_topmed[!grep("gnomad", clinvar_frequency_topmed$frequ, ignore.case = TRUE),]
clinvar_frequency_topmed_temp<-str_split(clinvar_frequency_topmed$frequ,"TOPMed\\)",simplify = TRUE)
clinvar_frequency_topmed_temp<-str_split(clinvar_frequency_topmed_temp[,2], ",", simplify = TRUE)
clinvar_frequency_topmed_temp[,1]<-gsub(" ","",clinvar_frequency_topmed_temp[,1])
clinvar_frequency_topmed$frequ1<-clinvar_frequency_topmed_temp[,1]
clinvar_frequency_topmed<-select(clinvar_frequency_topmed,id,frequ1)

clinvar_frequency_1000<-clinvar_frequency[grep("1000 Genomes Project", clinvar_frequency$frequ, ignore.case = TRUE),]
clinvar_frequency_1000<-clinvar_frequency_1000[!grep("gnomad", clinvar_frequency_1000$frequ, ignore.case = TRUE),]
clinvar_frequency_1000<-clinvar_frequency_1000[!grep("topmed", clinvar_frequency_1000$frequ, ignore.case = TRUE),]

clinvar_frequency_1000_temp<-str_split(clinvar_frequency_1000$frequ,"1000 Genomes Project ",simplify = TRUE)
clinvar_frequency_1000_temp<-str_split(clinvar_frequency_1000_temp[,2], ",", simplify = TRUE)
clinvar_frequency_1000$frequ1<-clinvar_frequency_1000_temp[,1]
clinvar_frequency_1000<-select(clinvar_frequency_1000,id,frequ1)

clinvar_frequency_rest<-clinvar_frequency[!grep("gnomad", clinvar_frequency$frequ, ignore.case = TRUE),]
clinvar_frequency_rest<-clinvar_frequency_rest[!grep("topmed", clinvar_frequency_rest$frequ, ignore.case = TRUE),]
clinvar_frequency_rest<-clinvar_frequency_rest[!grep("1000 Genomes Project", clinvar_frequency_rest$frequ, ignore.case = TRUE),]
clinvar_frequency_exac<-clinvar_frequency_rest[grep("exac", clinvar_frequency_rest$frequ, ignore.case = TRUE),]

clinvar_frequency_exac_temp<-str_split(clinvar_frequency_exac$frequ,"ExAC\\)",simplify = TRUE)
clinvar_frequency_exac_temp<-str_split(clinvar_frequency_exac_temp[,2], ",", simplify = TRUE)
clinvar_frequency_exac_temp[,1]<-gsub(" ","",clinvar_frequency_exac_temp[,1])
clinvar_frequency_exac$frequ1<-clinvar_frequency_exac_temp[,1]
clinvar_frequency_exac<-select(clinvar_frequency_exac,id,frequ1)

clinvar_frequency_rest<-clinvar_frequency_rest[!grep("exac", clinvar_frequency_rest$frequ, ignore.case = TRUE),]
clinvar_frequency_evs<-clinvar_frequency_rest[grep("Exome Variant Server", clinvar_frequency_rest$frequ, ignore.case = TRUE),]

clinvar_frequency_evs_temp<-str_split(clinvar_frequency_evs$frequ,"Exome Variant Server",simplify = TRUE)
clinvar_frequency_evs_temp<-str_split(clinvar_frequency_evs_temp[,2], ",", simplify = TRUE)
clinvar_frequency_evs_temp[,1]<-gsub(" ","",clinvar_frequency_evs_temp[,1])
clinvar_frequency_evs$frequ1<-clinvar_frequency_evs_temp[,1]
clinvar_frequency_evs<-select(clinvar_frequency_evs,id,frequ1)

clinvar_frequency_rest<-clinvar_frequency_rest[!grep("Exome Variant Server", clinvar_frequency_rest$frequ, ignore.case = TRUE),]
clinvar_frequency_T<-clinvar_frequency_rest[grep("\\(T\\)", clinvar_frequency_rest$frequ, ignore.case = TRUE),]
clinvar_frequency_T_temp<-str_match(clinvar_frequency_T$frequ, "(.*?)\\(T\\)")
clinvar_frequency_T$frequ1<-clinvar_frequency_T_temp[,2]
clinvar_frequency_T_temp<-str_split(clinvar_frequency_T$frequ1, "0\\.", simplify = TRUE)
clinvar_frequency_T$frequ1<-paste0("0.",clinvar_frequency_T_temp[,2])
clinvar_frequency_T<-select(clinvar_frequency_T,id,frequ1)

clinvar_frequency_A<-clinvar_frequency_rest[grep("\\(A\\)", clinvar_frequency_rest$frequ, ignore.case = TRUE),]
clinvar_frequency_A_temp<-str_match(clinvar_frequency_A$frequ, "(.*?)\\(A\\)")
clinvar_frequency_A_temp<-gsub(",","",clinvar_frequency_A_temp[,2])
clinvar_frequency_A$frequ1<-clinvar_frequency_A_temp
clinvar_frequency_A<-select(clinvar_frequency_A,id,frequ1)

clinvar_frequency_gnomad<-clinvar_frequency[grep("gnomad",clinvar_frequency$frequ, ignore.case = TRUE),]
clinvar_frequency_gnomad1<-separate(clinvar_frequency_gnomad, frequ, c("freq1","freq2","freq3","freq4","freq5","freq6","freq7"), sep = ",")
clinvar_frequency_gnomad1t<-clinvar_frequency_gnomad1[grep("exomes",clinvar_frequency_gnomad1[,2:8]),]
clinvar_frequency_gnomad1t<-with(clinvar_frequency_gnomad1, clinvar_frequency_gnomad1[grepl('exomes',freq1) | grepl('exomes',freq2) |
                                                                                        grepl('exomes',freq3) | grepl('exomes',freq4) |
                                                                                        grepl('exomes',freq5) | grepl('exomes',freq6), ])
clinvar_frequency_gnomad1b<-anti_join(clinvar_frequency_gnomad1, clinvar_frequency_gnomad1t, by = "id")
rclinvar_frequency_gnomad1b<-clinvar_frequency_gnomad1b[,1:2]
rclinvar_frequency_gnomad1b<-with(clinvar_frequency_gnomad1b, clinvar_frequency_gnomad1b[grepl('gnomAD',freq1) | grepl('gnomAD',freq2) |
                                                                                           grepl('gnomAD',freq3) | grepl('gnomAD',freq4) |
                                                                                           grepl('gnomAD',freq5) | grepl('gnomAD',freq6), ])
rclinvar_frequency_gnomad1b1<-clinvar_frequency_gnomad1b[grep("gnomAD", clinvar_frequency_gnomad1b$freq1, value = FALSE),]
rclinvar_frequency_gnomad1bf1<-str_split(rclinvar_frequency_gnomad1b1$freq1, "\\)", simplify = TRUE)
rclinvar_frequency_gnomad1b1$frequ1<-rclinvar_frequency_gnomad1bf1[,2]
rclinvar_frequency_gnomad1b1$frequ1<-gsub(" ","",rclinvar_frequency_gnomad1b1$frequ1)

rclinvar_frequency_gnomad1b2<-clinvar_frequency_gnomad1b[grep("gnomAD", clinvar_frequency_gnomad1b$freq2, value = FALSE),]
rclinvar_frequency_gnomad1bf2<-str_split(rclinvar_frequency_gnomad1b2$freq2, "\\)", simplify = TRUE)
rclinvar_frequency_gnomad1b2$frequ2<-rclinvar_frequency_gnomad1bf2[,2]
rclinvar_frequency_gnomad1b2$frequ2<-gsub(" ","",rclinvar_frequency_gnomad1b2$frequ2)

rclinvar_frequency_gnomad1b3<-clinvar_frequency_gnomad1b[grep("gnomAD", clinvar_frequency_gnomad1b$freq3, value = FALSE),]
rclinvar_frequency_gnomad1bf3<-str_split(rclinvar_frequency_gnomad1b3$freq3, "\\)", simplify = TRUE)
rclinvar_frequency_gnomad1b3$frequ3<-rclinvar_frequency_gnomad1bf3[,2]
rclinvar_frequency_gnomad1b3$frequ3<-gsub(" ","",rclinvar_frequency_gnomad1b3$frequ3)

rclinvar_frequency_gnomad1b4<-clinvar_frequency_gnomad1b[grep("gnomAD", clinvar_frequency_gnomad1b$freq4, value = FALSE),]
rclinvar_frequency_gnomad1bf4<-str_split(rclinvar_frequency_gnomad1b4$freq4, "\\)", simplify = TRUE)
rclinvar_frequency_gnomad1b4$frequ4<-rclinvar_frequency_gnomad1bf4[,2]
rclinvar_frequency_gnomad1b4$frequ4<-gsub(" ","",rclinvar_frequency_gnomad1b4$frequ4)

rclinvar_frequency_gnomad1b5<-clinvar_frequency_gnomad1b[grep("gnomAD", clinvar_frequency_gnomad1b$freq5, value = FALSE),]
rclinvar_frequency_gnomad1bf5<-str_split(rclinvar_frequency_gnomad1b5$freq5, "\\)", simplify = TRUE)
rclinvar_frequency_gnomad1b5$frequ5<-rclinvar_frequency_gnomad1bf5[,2]
rclinvar_frequency_gnomad1b5$frequ5<-gsub(" ","",rclinvar_frequency_gnomad1b5$frequ5)

rclinvar_frequency_gnomad1b1<-select(rclinvar_frequency_gnomad1b1, id, frequ1)
rclinvar_frequency_gnomad1b2<-select(rclinvar_frequency_gnomad1b2, id, frequ2)
rclinvar_frequency_gnomad1b3<-select(rclinvar_frequency_gnomad1b3, id, frequ3)
rclinvar_frequency_gnomad1b4<-select(rclinvar_frequency_gnomad1b4, id, frequ4)
rclinvar_frequency_gnomad1b5<-select(rclinvar_frequency_gnomad1b5, id, frequ5)

colnames(rclinvar_frequency_gnomad1b2)[2]<-"frequ1"
colnames(rclinvar_frequency_gnomad1b3)[2]<-"frequ1"
colnames(rclinvar_frequency_gnomad1b4)[2]<-"frequ1"
colnames(rclinvar_frequency_gnomad1b5)[2]<-"frequ1"

onlygnomad_clinvar_frequency<-rbind(rclinvar_frequency_gnomad1b1, rclinvar_frequency_gnomad1b2, rclinvar_frequency_gnomad1b3,
                                    rclinvar_frequency_gnomad1b4, rclinvar_frequency_gnomad1b5)
onlygnomad_clinvar_frequency$source<-"only_gnomAD"

exomes2<-clinvar_frequency_gnomad1t[grep("exomes", clinvar_frequency_gnomad1t$freq2),c(1,3)]
exomes3<-clinvar_frequency_gnomad1t[grep("exomes", clinvar_frequency_gnomad1t$freq3),c(1,4)]
exomes4<-clinvar_frequency_gnomad1t[grep("exomes", clinvar_frequency_gnomad1t$freq4),c(1,5)]
exomes5<-clinvar_frequency_gnomad1t[grep("exomes", clinvar_frequency_gnomad1t$freq5),c(1,6)]
exomes6<-clinvar_frequency_gnomad1t[grep("exomes", clinvar_frequency_gnomad1t$freq6),c(1,7)]

exomes2$freq2<-gsub(" exomes ","", exomes2$freq2)
exomes3$freq3<-gsub(" exomes ","", exomes3$freq3)
exomes4$freq4<-gsub(" exomes ","", exomes4$freq4)
exomes5$freq5<-gsub(" exomes ","", exomes5$freq5)
exomes6$freq6<-gsub(" exomes ","", exomes6$freq6)

colnames(exomes2)[2]<-"frequ1"
colnames(exomes3)[2]<-"frequ1"
colnames(exomes4)[2]<-"frequ1"
colnames(exomes5)[2]<-"frequ1"
colnames(exomes6)[2]<-"frequ1"

gnomadexome_clinvar_frequency<-rbind(exomes2, exomes3, exomes4, exomes5, exomes6)

all_clinvar_frequency<-rbind(gnomadexome_clinvar_frequency, clinvar_frequency_topmed, clinvar_frequency_1000,
                             clinvar_frequency_exac, clinvar_frequency_evs, clinvar_frequency_T, clinvar_frequency_A)

all_clinvar_frequency$frequ1<-as.numeric(all_clinvar_frequency$frequ1)


# Mouse and C. elegans Homology file

celegans_homology_blast<-fread("C_elegans_Human_Homologs_20.01.2019.csv") %>%
  dplyr::select(2,4,5,6,7) %>%
  unique()
celegans_homology_blast$Celegans_Symbol[celegans_homology_blast$Celegans_Symbol == ""]<-celegans_homology_blast$Celegans_Tag[celegans_homology_blast$Celegans_Symbol == ""]
celegans_homology_blast<-celegans_homology_blast[,1:4]
colnames(celegans_homology_blast)[c(1,4)]<-c("Gene_stable_ID", "human_gene_name")
celegans_homology_blast<-celegans_homology_blast[!celegans_homology_blast$human_gene_name == "-"]
celegans_homology_blast<-gnameConverter(celegans_homology_blast, "human_gene_name")

colnames(celegans_homology_blast)[2]<-"Gene_ID"
celegans_homology_blast<-celegans_homology_blast[,1:3]
celegans_homology_blast<-unique(celegans_homology_blast)

mouse_homology<-fread("mouse_human_homology_all.csv")
colnames(mouse_homology)[2]<-"human_gene_name"
mouse_homology<-gnameConverter(mouse_homology, "human_gene_name")
mouse_homology<-select(mouse_homology, c(1,4)) %>% 
  unique()

triplehomology<-merge(celegans_homology_blast, mouse_homology, by = "human_gene_name", all = TRUE)
triplehomology<-na.omit(triplehomology)
colnames(triplehomology)[c(3,4)]<-c("celegans_gene_name", "mouse_gene_name")


# Mouse Variants

mouse<-fread("mouse_variations.txt") %>% separate("Refseq_ID", "Refseq_ID", sep = "\\.")
mouse1<-str_split_fixed(mouse$`Amino Acid Change`, "->", 2)
mouse$from<-mouse1[,1]
mouse$to<-mouse1[,2]
colnames(mouse)[26:27]<-c("from", "to")
mouse$to[mouse$to == "Stop"]<-"*"
mouse<-mouse[,c(1,12,15,23,24,26,27)]
colnames(mouse)[1:7]<-c("Id","Gene_name","Phenotype","aapos","CCDS_ID","from","to")
mouse<-mouse %>% 
  mutate(CCDS_ID = strsplit(as.character(CCDS_ID), ",")) %>%
  unnest(CCDS_ID)
tempmouse<-str_split_fixed(mouse$CCDS_ID, "\\.", 2)
mouse$CCDS_ID<-tempmouse[,1]

mouse<-merge(mouse, mouseall, by = "CCDS_ID", all = TRUE)
mouse<-mouse[!is.na(mouse$aapos),]
mouse<-unique(mouse)
mouse$type<-"phenotypic"
mouse<-mouse[mouse$to != "",]
mouse$Variation<-paste0(mouse$from, mouse$aapos, mouse$to)

mouse<-mouse[,-1]
mouse<-mouse[!is.na(mouse$Refseq_ID),]

mouse<-mouse[mouse$from != "Stop",]

mouse_wconservation<-mouse
mouse_wconservation<-merge(mouse_wconservation, aaconservation, by = "from")
mouse_wconservation$conservation<-"radical"
mouse_wconservation$conservation[which(stri_detect_fixed(mouse_wconservation$to.y, mouse_wconservation$to.x))]<-"conservative"
mouse_wconservation<-mouse_wconservation[,-10]
colnames(mouse_wconservation)[6]<-"to"
mouse_wconservation$conservation[which(mouse_wconservation$from == mouse_wconservation$to)]<-"conservative"
mouse_wconservation$Source<-"APF"

# Mutagen Variants

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

mutagen2<-fread("mutagenetix_incidental_mutations_20210109.txt")
mutagen2<-mutagen2[order(quality_score)]
mutagen2<-unique(setDT(mutagen2), by = c("coordinate","ref_amino_acid","var_amino_acid"))
mutagen2<-dplyr::select(mutagen2, mutagenetix_id, coordinate, ref_amino_acid, var_amino_acid)
colnames(mutagen2)[3:4]<-c("from","to")

mutagen<-merge(mutagen, mutagen2, by = c("coordinate","from","to"), all = TRUE)
mutagen<-mutagen[!is.na(mutagen$Refseq_ID)]
colnames(mutagen)[12]<-"Variant_ID"
mutagen$Phenotype<-NA

mutagenphe<-fread("mutagenetix_phenotypic_mutations_20190909.txt") %>%
  dplyr::select(gene_symbol, aa_from, aa_to, aa_pos, fasta_header, sequence, predicted_effect, coordinate, mutation_type)
mutagenphe$Refseq_ID<-ncbimouseseq$refseq[match(mutagenphe$sequence, ncbimouseseq$.)]
mutagenphe<-separate(mutagenphe, "Refseq_ID", "Refseq_ID", sep = "\\.")
colnames(mutagenphe)[2:4]<-c("from","to","aapos")
mutagenphe<-mutagenphe[!mutagenphe$from == "",]
mutagenphe<-mutagenphe[!is.na(mutagenphe$Refseq_ID),]
mutagenphe<-mutagenphe[!is.na(mutagenphe$aapos),]
mutagenphe$type<-"phenotypic"
mutagenphe<-mutagenphe[mutagenphe$mutation_type !="makesense"]

mutagenphe2<-fread("mutagenetix_phenotypic_mutations_20201024.txt") %>% 
  dplyr::select(phenotypic_id, coordinate, ref_amino_acid, var_amino_acid, phenotypes)
colnames(mutagenphe2)[3:4]<-c("from","to")
mutagenphe<-merge(mutagenphe, mutagenphe2, by = c("coordinate","from","to"))
mutagenphe<-mutagenphe[!is.na(mutagenphe$Refseq_ID)]
colnames(mutagenphe)[12:13]<-c("Variant_ID","Phenotype")


mutagen<-rbind(mutagen, mutagenphe)
mutagen<-dplyr::select(mutagen, gene_symbol, Refseq_ID, aapos, from, to, type, predicted_effect, Phenotype, Variant_ID)
mutagen<-unique(mutagen)
colnames(mutagen)[c(1,3,9)]<-c("Gene_name","aapos","Id")
mutagen<-mutagen[,-7]
mutagen$Variation<-paste0(mutagen$from, mutagen$aapos, mutagen$to)

mutagen_wconservation<-mutagen
mutagen_wconservation<-merge(mutagen_wconservation, aaconservation, by = "from")
mutagen_wconservation$conservation<-"radical"
mutagen_wconservation$conservation[which(stri_detect_fixed(mutagen_wconservation$to.y, mutagen_wconservation$to.x))]<-"conservative"
mutagen_wconservation<-mutagen_wconservation[,-10]
colnames(mutagen_wconservation)[5]<-"to"
mutagen_wconservation$conservation[which(mutagen_wconservation$from == mutagen_wconservation$to)]<-"conservative"
mutagen_wconservation$Source<-"Mutagenetix"

# length(allmouse_wcons$from[allmouse_wcons$type == "phenotypic"])
# length(allmouse_wcons$from[allmouse_wcons$type == "incidental"])


# gnomAD Variants

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

gnomad_wconservation<-gnomad
gnomad_wconservation<-merge(gnomad_wconservation, aaconservation, by = "from")
gnomad_wconservation$conservation<-"radical"
gnomad_wconservation$conservation[which(stri_detect_fixed(gnomad_wconservation$to.y, gnomad_wconservation$to.x))]<-"conservative"
gnomad_wconservation<-gnomad_wconservation[,-11]
colnames(gnomad_wconservation)[5]<-"to"
gnomad_wconservation$conservation[which(gnomad_wconservation$from == gnomad_wconservation$to)]<-"conservative"


# COSMIC Variants

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

cosmic_wconservation<-cosmic
cosmic_wconservation<-merge(cosmic_wconservation, aaconservation, by = "from")
cosmic_wconservation$conservation<-"radical"
cosmic_wconservation$conservation[which(stri_detect_fixed(cosmic_wconservation$to.y, cosmic_wconservation$to.x))]<-"conservative"
cosmic_wconservation<-cosmic_wconservation[,c(-12,-4)]
colnames(cosmic_wconservation)[6]<-"to"
cosmic_wconservation$conservation[which(cosmic_wconservation$from == cosmic_wconservation$to)]<-"conservative"


# ClinVar Variants

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

clinvarxy_wconservation<-clinvarxy
clinvarxy_wconservation<-merge(clinvarxy_wconservation, aaconservation, by = "from")
clinvarxy_wconservation$conservation<-"radical"
clinvarxy_wconservation$conservation[which(stri_detect_fixed(clinvarxy_wconservation$to.y, clinvarxy_wconservation$to.x))]<-"conservative"
clinvarxy_wconservation<-clinvarxy_wconservation[,-12]
colnames(clinvarxy_wconservation)[7]<-"to"
clinvarxy_wconservation$conservation[which(clinvarxy_wconservation$from == clinvarxy_wconservation$to)]<-"conservative"
clinvarxy_wconservation<-clinvarxy_wconservation[,-2]


# dbSNP Variants

dbsnpx<-fread("dbsnp.csv")
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

dbsnp_wconservation<-dbsnp
dbsnp_wconservation<-merge(dbsnp_wconservation, aaconservation, by = "from")
dbsnp_wconservation$conservation<-"radical"
dbsnp_wconservation$conservation[which(stri_detect_fixed(dbsnp_wconservation$to.y, dbsnp_wconservation$to.x))]<-"conservative"
dbsnp_wconservation<-dbsnp_wconservation[,-11]
colnames(dbsnp_wconservation)[4]<-"to"
dbsnp_wconservation$conservation[which(dbsnp_wconservation$from == dbsnp_wconservation$to)]<-"conservative"


# C. elegans Variants

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

celegans_wconservation<-celegans
celegans_wconservation<-merge(celegans_wconservation, aaconservation, by = "from")
celegans_wconservation$conservation<-"radical"
celegans_wconservation$conservation[which(stri_detect_fixed(celegans_wconservation$to.y, celegans_wconservation$to.x))]<-"conservative"
celegans_wconservation<-celegans_wconservation[,-9]
colnames(celegans_wconservation)[6]<-"to"
celegans_wconservation$conservation[which(celegans_wconservation$from == celegans_wconservation$to)]<-"conservative"
celegans_wconservation$Source<-"Wormbase"


celegans_manual<-fread("manual_variant_list.txt", select = c(1:7), header = TRUE)

celegans_manualx<-str_split(celegans_manual[[2]], "[0-9]", simplify = TRUE)
celegans_manual$from<-celegans_manualx[,1]
celegans_manual$to<-apply(celegans_manualx, 1, function(x) paste0(x[3],x[4],x[5]))
celegans_manualx1<-str_split(celegans_manual[[2]], "[A-Z]", simplify = TRUE)
celegans_manualx1<-str_split(celegans_manualx1[,2], "\\*", simplify = TRUE)
celegans_manual$aapos<-celegans_manualx1[,1]
celegans_manual$aapos<-as.integer(celegans_manual$aapos)
colnames(celegans_manual)[c(1,3,4,5)]<-c("Refseq_ID","type","Id","Gene_name")
celegans_manual<-celegans_manual[,-2]
celegans_manual_wconservation<-merge(celegans_manual, aaconservation, by = "from")
celegans_manual_wconservation$conservation<-"radical"
celegans_manual_wconservation$conservation[which(stri_detect_fixed(celegans_manual_wconservation$to.y, celegans_manual_wconservation$to.x))]<-"conservative"
celegans_manual_wconservation<-celegans_manual_wconservation[,-10]
colnames(celegans_manual_wconservation)[8]<-"to"

celegans$Source<-"Wormbase"
#celegans_wm<-anti_join(celegans_manual, celegans, by = c("Refseq_ID","from","to","aapos"))
celegans<-rbind(celegans, celegans_manual)

# celegans_manual_wconservation<-fread("wan_id_edited.csv")
# celegans_manual_wconservation$type<-"Phenotypic"
# 
# celegans_manual_wconservation<-merge(celegans_manual_wconservation, aaconservation, by = "from")
# celegans_manual_wconservation$conservation<-"radical"
# celegans_manual_wconservation$to.x[10]<-"V"
# celegans_manual_wconservation$conservation[which(stri_detect_fixed(celegans_manual_wconservation$to.y, celegans_manual_wconservation$to.x))]<-"conservative"
# celegans_manual_wconservation<-celegans_manual_wconservation[,c(-3,-5,-11)]
# colnames(celegans_manual_wconservation)[5]<-"to"
# celegans_manual_wconservation$conservation[which(celegans_manual_wconservation$from == celegans_manual_wconservation$to)]<-"conservative"
# colnames(celegans_manual_wconservation)[7]<-"Id"

# Combine variants by organism

human_wconservation<-rbind(clinvarxy_wconservation, dbsnp_wconservation, cosmic_wconservation, gnomad_wconservation)

allmouse_wcons<-rbind(mouse_wconservation, mutagen_wconservation) %>% unique()
# allmouse_wcons<-allmouse_wcons[,c(-1,-5,-6)]
# colnames(allmouse_wcons)[c(1,4,5,6,7)]<-c("Variant_ID","Protein_ID","Phenotypic_Significance","Variant","Conservation")
# write.table(allmouse_wcons,"mouse_variants_id.txt", quote = FALSE, row.names = FALSE, sep = "\t")
# 
# allmouse_wcons<-unique(setDT(allmouse_wcons), by = c(7,8,9))
colnames(allmouse_wcons)[3]<-"Mouse_Gene_name"

celegans_wconservation<-rbind(celegans_wconservation, celegans_manual_wconservation) %>% unique()
celegans_wconservation$aapos<-as.numeric(celegans_wconservation$aapos)



# Supp. Tables




