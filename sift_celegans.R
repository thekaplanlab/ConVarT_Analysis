

nohuman_celegans<-fread("no_human_celegans_phenotypic_variations_.tsv")
nohuman_celegans_sift<-data.frame(var = paste0(nohuman_celegans$C_elegans_Protein_ID, ":p.", nohuman_celegans$C_elegans_Variation),
                                  stringsAsFactors = FALSE) %>% unique()
write.table(nohuman_celegans_sift, "nohuman_celegans_sift.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


celegans_sift<-fread("celegans_nohuman_sift_results.txt")
celegans_ref_sift<-fread("nohuman_celegans_ens_ref.txt")

celegans_ref_sift<-celegans_ref_sift[,c(1,2,7,9)]
celegans_ref_sift1<-str_split(celegans_ref_sift$Feature, "\\.", simplify = TRUE)

celegans_ref_sift$transcript<-apply(celegans_ref_sift1[,1:2] , 1 , paste , collapse = ".")
celegans_ref_sift2<-str_split(celegans_ref_sift$`#Uploaded_variation`, ":", simplify = TRUE)
celegans_ref_sift$`#Uploaded_variation`<-celegans_ref_sift2[,1]
celegans_ref_sift3<-str_split(celegans_ref_sift2[,2], "\\.", simplify = TRUE)
celegans_ref_sift$Variation<-celegans_ref_sift3[,2]

celegans_ref_sift4<-str_split(celegans_ref_sift$Location, "-", simplify = TRUE)
celegans_ref_sift$Location<-celegans_ref_sift4[,2]

celegans_sift1<-celegans_sift[,c(2,5,9:17)]
colnames(celegans_sift1)[c(1,2)]<-c("Location","transcript")
celegans_sift1$Location<-as.integer(celegans_sift1$Location)
celegans_ref_sift$Location<-as.integer(celegans_ref_sift$Location)
celegans_sift_ref<-merge(celegans_ref_sift, celegans_sift1, by = c("Location","transcript"), all = TRUE, allow.cartesian = TRUE) %>%
  unique()
celegans_sift_ref<-celegans_sift_ref[!is.na(celegans_sift_ref$AMINO_POS)]
celegans_sift_ref11<-celegans_sift_ref[!is.na(celegans_sift_ref$`#Uploaded_variation`)]


# ********************************* #







celegans_sift<-fread("celegans_nohuman_sift_results.txt")
celegans_ref_sift<-fread("nohuman_celegans_ens_ref.txt")

celegans_ref_sift<-celegans_ref_sift[,c(1,2,7,9)]
celegans_ref_sift$transcript<-celegans_ref_sift$Feature
celegans_ref_sift2<-str_split(celegans_ref_sift$`#Uploaded_variation`, ":", simplify = TRUE)
celegans_ref_sift$`#Uploaded_variation`<-celegans_ref_sift2[,1]
celegans_ref_sift3<-str_split(celegans_ref_sift2[,2], "\\.", simplify = TRUE)
celegans_ref_sift$Variation<-celegans_ref_sift3[,2]

celegans_ref_sift4<-str_split(celegans_ref_sift$Location, "-", simplify = TRUE)
celegans_ref_sift$Location<-celegans_ref_sift4[,2]

celegans_sift1<-celegans_sift[,c(2,5,9:17)]
colnames(celegans_sift1)[c(1,2)]<-c("Location","transcript")
celegans_sift1$Location<-as.integer(celegans_sift1$Location)
celegans_ref_sift$Location<-as.integer(celegans_ref_sift$Location)
celegans_sift_ref2<-merge(celegans_ref_sift, celegans_sift1, by = c("Location","transcript"), all = TRUE, allow.cartesian = TRUE) %>%
  unique()
celegans_sift_ref2<-celegans_sift_ref2[!is.na(celegans_sift_ref2$AMINO_POS)]

celegans_sift_ref21<-celegans_sift_ref2[!is.na(celegans_sift_ref2$`#Uploaded_variation`)]


# ****************************************************** #


celegans_sift_ref3<-rbind(celegans_sift_ref11, celegans_sift_ref21) %>% unique()
celegans_sift_scores<-celegans_sift_ref3[,c(3,4,6,11:15)]
colnames(celegans_sift_scores)[c(1,3)]<-c("C_elegans_Protein_ID","C_elegans_Variation")
celegans_sift_scores<-merge(nohuman_celegans, celegans_sift_scores, by = c("C_elegans_Protein_ID","C_elegans_Variation"), 
                            all = TRUE) %>% unique()


celegans_homology_blast_temp<-celegans_homology_blast[,1:2]
colnames(celegans_homology_blast_temp)[2]<-"Gene"
celegans_sift_scores$Gene[is.na(celegans_sift_scores$Gene)]<-celegans_sift_scores$C_elegans_Gene_name[is.na(celegans_sift_scores$Gene)]
celegans_sift_scores_hom<-merge(celegans_sift_scores, celegans_homology_blast_temp, by = "Gene", all = TRUE)
celegans_sift_scores_hom<-celegans_sift_scores_hom[!is.na(celegans_sift_scores_hom$C_elegans_Protein_ID)]

celegans_homology_blast_temp2<-celegans_homology_blast[,c(1,3)]
colnames(celegans_homology_blast_temp2)[2]<-"Gene"
celegans_sift_scores_hom2<-merge(celegans_sift_scores, celegans_homology_blast_temp2, by = "Gene", all = TRUE)
celegans_sift_scores_hom2<-celegans_sift_scores_hom2[!is.na(celegans_sift_scores_hom2$C_elegans_Protein_ID)]

celegans_sift_scores_hom<-celegans_sift_scores_hom[!is.na(celegans_sift_scores_hom$human_gene_name)]
celegans_sift_scores_hom2<-celegans_sift_scores_hom2[!is.na(celegans_sift_scores_hom2$human_gene_name)]
celegans_sift_scores_hom3<-rbind(celegans_sift_scores_hom, celegans_sift_scores_hom2)



# ******************************************************** #


diseases_alliance_human<-fread("DISEASE-ALLIANCE_HUMAN_37.tsv") %>% select(c(5,8))
colnames(diseases_alliance_human)<-c("human_gene_name","Disease")
diseases_alliance_human1<-aggregate(Disease ~ ., diseases_alliance_human, toString)

celegans_sift_scores_hom3_disease<-merge(celegans_sift_scores_hom3, diseases_alliance_human1, by = "human_gene_name", all = TRUE)
celegans_sift_scores_hom3_disease<-celegans_sift_scores_hom3_disease[!is.na(celegans_sift_scores_hom3_disease$C_elegans_Protein_ID)] %>%
  unique()


celegans_sift_scores_hom3_disease_na<-celegans_sift_scores_hom3_disease[is.na(celegans_sift_scores_hom3_disease$Disease),-10]

diseases_alliance_celegans<-fread("DISEASE-ALLIANCE_WB_37.tsv") %>% select(c(5,8))
colnames(diseases_alliance_celegans)<-c("C_elegans_Gene_name","Disease")
diseases_alliance_celegans1<-aggregate(Disease ~ ., diseases_alliance_celegans, toString)

celegans_sift_scores_hom3_disease1<-merge(celegans_sift_scores_hom3_disease_na, diseases_alliance_celegans1, 
                                          by = "C_elegans_Gene_name", all = TRUE)
celegans_sift_scores_hom3_disease1<-celegans_sift_scores_hom3_disease1[!is.na(celegans_sift_scores_hom3_disease1$C_elegans_Protein_ID)] %>%
  unique()


celegans_sift_scores_hom3_disease2<-rbind(celegans_sift_scores_hom3_disease, celegans_sift_scores_hom3_disease1) %>% unique()

celegans_sift_scores_hom3_disease2<-celegans_sift_scores_hom3_disease2[,c(-2,-9,-13)]
colnames(celegans_sift_scores_hom3_disease2)[c(1,12)]<-c("Human_Gene_name","Human_Disease")
celegans_sift_scores_hom3_disease2<-setcolorder(celegans_sift_scores_hom3_disease2, c(2:7,1,12,8:11))
celegans_sift_scores_hom3_disease2[is.na(celegans_sift_scores_hom3_disease2)]<-""

write.table(celegans_sift_scores_hom3_disease2, "celegans_nohuman_wSIFT_disease.txt", row.names = FALSE, quote = FALSE, sep = "\t")


xxyyzz<-semi_join(celegans_sift_scores_hom3_disease2, ort_humancelegans, by = c("C_elegans_Protein_ID","C_elegans_Variation"))
xxyyzz<-semi_join(ort_humancelegans, celegans_sift_scores_hom3_disease2, by = c("C_elegans_Protein_ID","C_elegans_Variation"))


# C. elegans Table

celegans_table<-fread("celegans_table.txt")



# Mouse 

nohuman_mouse<-fread("mouse_phenotypic_variations_no_human.txt")

nohuman_mouse<-allmouse_nohuman_phenotypic


mousex<-fread("mouse_variations.txt") %>% separate("Refseq_ID", "Refseq_ID", sep = "\\.")
mousex1<-str_split_fixed(mousex$`Amino Acid Change`, "->", 2)
mousex$from<-mousex1[,1]
mousex$to<-mousex1[,2]
colnames(mousex)[26:27]<-c("from", "to")
mousex$to[mousex$to == "Stop"]<-"*"
mousex<-mousex[,c(1,12,15,23,24,26,27,7,8,9,10)]
colnames(mousex)[1:7]<-c("Id","Gene_name","Phenotype","aapos","CCDS_ID","from","to")
mousex<-mousex %>% 
  mutate(CCDS_ID = strsplit(as.character(CCDS_ID), ",")) %>%
  unnest(CCDS_ID)
tempmousex<-str_split_fixed(mousex$CCDS_ID, "\\.", 2)
mousex$CCDS_ID<-tempmousex[,1]

mousex<-merge(mousex, mouseall, by = "CCDS_ID", all = TRUE)
mousex<-mousex[!is.na(mousex$aapos),]
mousex<-unique(mousex)
mousex$type<-"phenotypic"
mousex<-mousex[mousex$to != "",]
mousex$Variation<-paste0(mousex$from, mousex$aapos, mousex$to)


mousex<-mousex[,2:14]
mousex<-mousex[!is.na(mousex$Refseq_ID),]
mousex<-mousex[mousex$from != "Stop",]

mousex_wconservation<-mousex
mousex_wconservation<-merge(mousex_wconservation, aaconservation, by = "from")
mousex_wconservation$conservation<-"radical"
mousex_wconservation$conservation[which(stri_detect_fixed(mousex_wconservation$to.y, mousex_wconservation$to.x))]<-"conservative"
mousex_wconservation<-mousex_wconservation[,-14]
colnames(mousex_wconservation)[6]<-"to"
mousex_wconservation$conservation[which(mousex_wconservation$from == mousex_wconservation$to)]<-"conservative"





# mutagen

mutagenx<-fread("mutagenetix_incidental_mutations_20190904.txt") %>%
  dplyr::select(gene_symbol, aa_from, aa_to, aa_pos, fasta_header, sequence, predicted_effect, mutation_type, coordinate, pph_score)
ncbimouseseq<-read.fasta(file = "GCF_000001635.26_GRCm38.p6_protein.faa", seqtype = "AA", as.string = TRUE, set.attributes = TRUE) %>%
  unlist() %>%
  as.data.frame()
ncbimouseseq$refseq<-rownames(ncbimouseseq)
mutagenx$Refseq_ID<-ncbimouseseq$refseq[match(mutagenx$sequence, ncbimouseseq$.)]
mutagenx<-separate(mutagenx, "Refseq_ID", "Refseq_ID", sep = "\\.")
colnames(mutagenx)[2:4]<-c("from","to","aapos")
mutagenx<-mutagenx[!mutagenx$from == "",]
mutagenx<-mutagenx[!is.na(mutagenx$Refseq_ID),]
mutagenx$type<-"incidental"
mutagenx<-mutagenx[mutagenx$mutation_type != "makesense"]

mutagenxphe<-fread("mutagenetix_phenotypic_mutations_20190909.txt")%>%
  dplyr::select(gene_symbol, aa_from, aa_to, aa_pos, fasta_header, sequence, predicted_effect, coordinate, mutation_type, pph_score)
mutagenxphe$Refseq_ID<-ncbimouseseq$refseq[match(mutagenxphe$sequence, ncbimouseseq$.)]
mutagenxphe<-separate(mutagenxphe, "Refseq_ID", "Refseq_ID", sep = "\\.")
colnames(mutagenxphe)[2:4]<-c("from","to","aapos")
mutagenxphe<-mutagenxphe[!mutagenxphe$from == "",]
mutagenxphe<-mutagenxphe[!is.na(mutagenxphe$Refseq_ID),]
mutagenxphe<-mutagenxphe[!is.na(mutagenxphe$aapos),]
mutagenxphe$type<-"phenotypic"
mutagenxphe<-mutagenxphe[mutagenxphe$mutation_type !="makesense"]

mutagenxphe2<-fread("mutagenetix_phenotypic_mutations_20201024.txt")
mutagenxphe$Phenotype<-mutagenxphe2$phenotypes[match(mutagenxphe$coordinate, mutagenxphe2$coordinate)]
mutagenx$Phenotype<-NA


mutagenx<-rbind(mutagenx, mutagenxphe)
mutagenx<-dplyr::select(mutagenx, gene_symbol, Refseq_ID, aapos, from, to, type, predicted_effect, Phenotype, coordinate, pph_score)
mutagenx<-unique(mutagenx)
colnames(mutagenx)[c(1,3,9)]<-c("Gene_name","aapos","Id")
mutagenx$Variation<-paste0(mutagenx$from, mutagenx$aapos, mutagenx$to)

mutagenx_wconservation<-mutagenx
mutagenx_wconservation<-merge(mutagenx_wconservation, aaconservation, by = "from")
mutagenx_wconservation$conservation<-"radical"
mutagenx_wconservation$conservation[which(stri_detect_fixed(mutagenx_wconservation$to.y, mutagenx_wconservation$to.x))]<-"conservative"
mutagenx_wconservation<-mutagenx_wconservation[,-12]
colnames(mutagenx_wconservation)[5]<-"to"
mutagenx_wconservation$conservation[which(mutagenx_wconservation$from == mutagenx_wconservation$to)]<-"conservative"

colnames(mutagenx_wconservation)[c(7,10)]<-colnames(mousex_wconservation)[c(8,7)]
btl<-paste0(mutagenx_wconservation$Refseq_ID, ":p.", mutagenx_wconservation$Variation)

btl2<-data.frame(btl, stringsAsFactors = FALSE)
write.table(btl2, "mutagenforsift.txt", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")



sift_mouse1<-fread("mutagen_sift_1.txt", select = c(1,14,17,26,30))
sift_mouse2<-fread("mutagen_sift_2.txt", select = c(1,14,17,27,31))
sift_mouse3<-fread("mutagen_sift_3.txt", select = c(1,14,17,26,30))
sift_mouse4<-fread("mutagen_sift_4.txt", select = c(1,14,17,26,30))
sift_mouse5<-fread("mutagen_sift_5.txt", select = c(1,14,17,26,30))

sift_mouse<-rbind(sift_mouse1,sift_mouse2,sift_mouse3,sift_mouse4,sift_mouse5)

sift_mousex<-str_split(sift_mouse$`#Uploaded_variation`, ":p.", simplify = TRUE)
sift_mouse$Mouse_Protein_ID<-sift_mousex[,1]
sift_mouse$Mouse_Variation<-sift_mousex[,2]

sift_mousey<-str_split(sift_mouse$ENSP, "\\.", simplify = TRUE)
sift_mouse$ENSP<-sift_mousey[,1]

sift_mousez<-str_split(sift_mouse$SIFT, "\\(", simplify = TRUE)
sift_mouse$`SIFT Prediction`<-sift_mousez[,1]
sift_mousez2<-str_split(sift_mousez[,2], "\\)", simplify = TRUE)
sift_mouse$`SIFT Score`<-sift_mousez2[,1]

sift_mouse_o<-sift_mouse[sift_mouse$ENSP == sift_mouse$Mouse_Protein_ID]
sift_mouse_o2<-sift_mouse_o[,6:9]

write.table(sift_mouse_o2, "sift_mutagenetix_processed.txt", row.names = FALSE, quote = FALSE, sep = "\t")

mutagenx_wconservation_1<-mutagenx_wconservation[mutagenx_wconservation$type == "phenotypic",c(3,7,10,11)]

mutagenx_wconservation_1$`Polyphen Prediction`[mutagenx_wconservation_1$`Polyphen Prediction` == "1"]<-"Unknown"
mutagenx_wconservation_1$`Polyphen Prediction`[mutagenx_wconservation_1$`Polyphen Prediction` == "2"]<-"Possibly damaging"
mutagenx_wconservation_1$`Polyphen Prediction`[mutagenx_wconservation_1$`Polyphen Prediction` == "3"]<-"Probably damaging"
mutagenx_wconservation_1$`Polyphen Prediction`[mutagenx_wconservation_1$`Polyphen Prediction` == "4"]<-"Probably null"
mutagenx_wconservation_1$`Polyphen Prediction`[mutagenx_wconservation_1$`Polyphen Prediction` == "5"]<-"Probably benign"
mutagenx_wconservation_1$`Polyphen Prediction`[mutagenx_wconservation_1$`Polyphen Prediction` == "9"]<-"Silent"

colnames(mutagenx_wconservation_1)[c(1,4)]<-c("Mouse_Protein_ID", "Mouse_Variation")
mutagenx_wconservation_1<-merge(mutagenx_wconservation_1, sift_mouse_o2, by = c("Mouse_Protein_ID", "Mouse_Variation"), all = TRUE)
mutagenx_wconservation_1<-mutagenx_wconservation_1[!is.na(mutagenx_wconservation_1$`Polyphen Prediction`)]

mousex_wconservation_1<-mousex_wconservation[,c(7,8,9,10,11,13)]
colnames(mousex_wconservation_1)[c(5,6)]<-c("Mouse_Protein_ID", "Mouse_Variation")

allmousex_wconservation_1<-rbind(mousex_wconservation_1, mutagenx_wconservation_1)

nohuman_mouse_sift<-merge(nohuman_mouse, allmousex_wconservation_1, by = c("Mouse_Protein_ID", "Mouse_Variation"), all = TRUE)
nohuman_mouse_sift<-nohuman_mouse_sift[!is.na(nohuman_mouse_sift$Mouse_Gene_name)]

nohuman_mouse_sift[is.na(nohuman_mouse_sift)]<-""
nohuman_mouse_sift<-unique(nohuman_mouse_sift)


mouse_homology<-fread("mouse_human_homology_all.csv")
colnames(mouse_homology)[c(2,4)]<-c("Human_Gene_name", "Mouse_Gene_name")
mouse_homology<-gnameConverter(mouse_homology, "Human_Gene_name")
mouse_homology<-select(mouse_homology, c(1,4)) %>% 
  unique()

nohuman_mouse_sift_whom<-merge(nohuman_mouse_sift, mouse_homology, by = "Mouse_Gene_name")

diseases_alliance_mouse<-fread("DISEASE-ALLIANCE_MGI_37.tsv")
diseases_alliance_mouse<-diseases_alliance_mouse[diseases_alliance_mouse$DBobjectType == "gene"]
diseases_alliance_mouse<-diseases_alliance_mouse[,c(5,8)] %>% unique()
colnames(diseases_alliance_mouse)<-c("Mouse_Gene_name","Disease")
diseases_alliance_mouse1<-aggregate(Disease ~ ., diseases_alliance_mouse, toString)

nohuman_mouse_sift_whom_dis<-merge(nohuman_mouse_sift_whom, diseases_alliance_mouse1, by = "Mouse_Gene_name", all = TRUE)
nohuman_mouse_sift_whom_dis<-nohuman_mouse_sift_whom_dis[!is.na(nohuman_mouse_sift_whom_dis$Mouse_Variation)]

na_nohuman_mouse_sift_whom_dis<-nohuman_mouse_sift_whom_dis[is.na(nohuman_mouse_sift_whom_dis$Disease),-13]
colnames(diseases_alliance_human1)[1]<-"Human_Gene_name"
na_nohuman_mouse_sift_whom_dis<-merge(na_nohuman_mouse_sift_whom_dis, diseases_alliance_human1, by = "Human_Gene_name", all = TRUE)
na_nohuman_mouse_sift_whom_dis<-na_nohuman_mouse_sift_whom_dis[!is.na(na_nohuman_mouse_sift_whom_dis$Mouse_Variation)]

nohuman_mouse_sift_whom_dis<-rbind(nohuman_mouse_sift_whom_dis, na_nohuman_mouse_sift_whom_dis)
nohuman_mouse_sift_whom_dis<-unique(nohuman_mouse_sift_whom_dis)
nohuman_mouse_sift_whom_dis<-setcolorder(nohuman_mouse_sift_whom_dis, c(1:6,12,13))
colnames(nohuman_mouse_sift_whom_dis)[8]<-"Human_Disease"
nohuman_mouse_sift_whom_dis<-nohuman_mouse_sift_whom_dis[,-9]


write.table(nohuman_mouse_sift_whom_dis, "mouse_nohuman_wSIFT_PP_disease.txt", sep = "\t", row.names = FALSE, quote = FALSE)



