

mouse_celegans<-readRDS("mouse_celegans_05.11.20.RDS")
mutagen_celegans<-readRDS("mutagen_celegans_03.11.20.RDS")
clinvar_celegans<-readRDS("clinvar_celegans_02.11.20.RDS")
clinvar2_celegans<-readRDS("clinvar2_celegans_22.11.20.RDS")
dbsnp_celegans<-readRDS("dbsnp_celegans_03.11.20.RDS")
cosmic_celegans<-readRDS("cosmic_celegans_30.10.20.RDS")
gnomad_celegans<-readRDS("gnomad_celegans_02.11.20.RDS")
mouse_clinvar<-readRDS("mouse_clinvar_05.11.20.RDS")
mouse_clinvar2<-readRDS("mouse_clinvar2_23.11.20.RDS")
mouse_dbsnp<-readRDS("mouse_dbsnp_05.11.20.RDS")
mouse_cosmic<-readRDS("mouse_cosmic_04.11.20.RDS")
mouse_gnomad<-readRDS("mouse_gnomad_04.11.20.RDS")

mutagen_clinvar<-readRDS("mutagen_clinvar_30.10.20.RDS")
mutagen_clinvar2<-readRDS("mutagen_clinvar2_23.11.20.RDS")
mutagen_dbsnp<-readRDS("mutagen_dbsnp_03.11.20.RDS")
mutagen_cosmic<-readRDS("mutagen_cosmic_28.10.20.RDS")
mutagen_gnomad<-readRDS("mutagen_gnomad_28.10.20.RDS")




mousecelegans<-data.frame(mouse_Refseq = mouse_celegans[[2]], mouse_aa = mouse_celegans[[1]], 
                          celegans_Refseq = mouse_celegans[[4]], celegans_aa = mouse_celegans[[3]], 
                          msa_index = mouse_celegans[[5]], stringsAsFactors = FALSE) %>% unique()
mutagencelegans<-data.frame(mouse_Refseq = mutagen_celegans[[2]], mouse_aa = mutagen_celegans[[1]], 
                            celegans_Refseq = mutagen_celegans[[4]], celegans_aa = mutagen_celegans[[3]], 
                            msa_index = mutagen_celegans[[5]], stringsAsFactors = FALSE) %>% unique()

clinvarcelegans<-data.frame(human_Refseq = clinvar_celegans[[2]], human_aa = as.numeric(clinvar_celegans[[1]]), 
                            celegans_Refseq = clinvar_celegans[[4]], celegans_aa = clinvar_celegans[[3]], 
                            msa_index = clinvar_celegans[[5]], stringsAsFactors = FALSE) %>% unique()
clinvar2celegans<-data.frame(human_Refseq = clinvar2_celegans[[2]], human_aa = as.numeric(clinvar2_celegans[[1]]), 
                             celegans_Refseq = clinvar2_celegans[[4]], celegans_aa = clinvar2_celegans[[3]], 
                             msa_index = clinvar2_celegans[[5]], stringsAsFactors = FALSE) %>% unique()
clinvarcelegans<-rbind(clinvarcelegans, clinvar2celegans) %>% unique()

dbsnpcelegans<-data.frame(human_Refseq = dbsnp_celegans[[2]], human_aa = dbsnp_celegans[[1]], 
                          celegans_Refseq = dbsnp_celegans[[4]], celegans_aa = dbsnp_celegans[[3]], 
                          msa_index = dbsnp_celegans[[5]], stringsAsFactors = FALSE) %>% unique()
dbsnpcelegans$human_aa<-as.integer(dbsnpcelegans$human_aa)

cosmiccelegans<-data.frame(human_Refseq = cosmic_celegans[[2]], human_aa = cosmic_celegans[[1]], 
                           celegans_Refseq = cosmic_celegans[[4]], celegans_aa = cosmic_celegans[[3]], 
                           msa_index = cosmic_celegans[[5]], stringsAsFactors = FALSE) %>% unique()

gnomadcelegans<-data.frame(human_Refseq = gnomad_celegans[[2]], human_aa = gnomad_celegans[[1]], 
                           celegans_Refseq = gnomad_celegans[[4]], celegans_aa = gnomad_celegans[[3]], 
                           msa_index = gnomad_celegans[[5]], stringsAsFactors = FALSE) %>% unique()

mouseclinvar<-data.frame(mouse_Refseq = mouse_clinvar[[2]], mouse_aa = mouse_clinvar[[1]], 
                         human_Refseq = mouse_clinvar[[4]], human_aa = mouse_clinvar[[3]], 
                         msa_index = mouse_clinvar[[5]], stringsAsFactors = FALSE) %>% unique()
mouseclinvar$human_aa<-as.integer(mouseclinvar$human_aa)

mouseclinvar2<-data.frame(mouse_Refseq = mouse_clinvar2[[2]], mouse_aa = mouse_clinvar2[[1]], 
                          human_Refseq = mouse_clinvar2[[4]], human_aa = mouse_clinvar2[[3]], 
                          msa_index = mouse_clinvar2[[5]], stringsAsFactors = FALSE) %>% unique()
mouseclinvar2$human_aa<-as.integer(mouseclinvar2$human_aa)
mouseclinvar<-rbind(mouseclinvar, mouseclinvar2)

mousedbsnp<-data.frame(mouse_Refseq = mouse_dbsnp[[2]], mouse_aa = mouse_dbsnp[[1]], 
                       human_Refseq = mouse_dbsnp[[4]], human_aa = mouse_dbsnp[[3]], 
                       msa_index = mouse_dbsnp[[5]], stringsAsFactors = FALSE) %>% unique()
mousedbsnp$human_aa<-as.integer(mousedbsnp$human_aa)

mousecosmic<-data.frame(mouse_Refseq = mouse_cosmic[[2]], mouse_aa = mouse_cosmic[[1]], 
                        human_Refseq = mouse_cosmic[[4]], human_aa = mouse_cosmic[[3]], 
                        msa_index = mouse_cosmic[[5]], stringsAsFactors = FALSE) %>% unique()

mousegnomad<-data.frame(mouse_Refseq = mouse_gnomad[[2]], mouse_aa = mouse_gnomad[[1]], 
                        human_Refseq = mouse_gnomad[[4]], human_aa = mouse_gnomad[[3]], 
                        msa_index = mouse_gnomad[[5]], stringsAsFactors = FALSE) %>% unique()

mutagenclinvar<-data.frame(mouse_Refseq = mutagen_clinvar[[2]], mouse_aa = mutagen_clinvar[[1]], 
                           human_Refseq = mutagen_clinvar[[4]], human_aa = mutagen_clinvar[[3]], 
                           msa_index = mutagen_clinvar[[5]], stringsAsFactors = FALSE) %>% unique()
mutagenclinvar$human_aa<-as.integer(mutagenclinvar$human_aa)

mutagenclinvar2<-data.frame(mouse_Refseq = mutagen_clinvar2[[2]], mouse_aa = mutagen_clinvar2[[1]], 
                            human_Refseq = mutagen_clinvar2[[4]], human_aa = mutagen_clinvar2[[3]], 
                            msa_index = mutagen_clinvar2[[5]], stringsAsFactors = FALSE) %>% unique()
mutagenclinvar2$human_aa<-as.integer(mutagenclinvar2$human_aa)
mutagenclinvar<-rbind(mutagenclinvar, mutagenclinvar2)


mutagendbsnp<-data.frame(mouse_Refseq = mutagen_dbsnp[[2]], mouse_aa = mutagen_dbsnp[[1]], 
                         human_Refseq = mutagen_dbsnp[[4]], human_aa = mutagen_dbsnp[[3]], 
                         msa_index = mutagen_dbsnp[[5]], stringsAsFactors = FALSE) %>% unique()
mutagendbsnp$human_aa<-as.integer(mutagendbsnp$human_aa)

mutagencosmic<-data.frame(mouse_Refseq = mutagen_cosmic[[2]], mouse_aa = mutagen_cosmic[[1]], 
                          human_Refseq = mutagen_cosmic[[4]], human_aa = mutagen_cosmic[[3]], 
                          msa_index = mutagen_cosmic[[5]], stringsAsFactors = FALSE) %>% unique()

mutagengnomad<-data.frame(mouse_Refseq = mutagen_gnomad[[2]], mouse_aa = mutagen_gnomad[[1]], 
                          human_Refseq = mutagen_gnomad[[4]], human_aa = mutagen_gnomad[[3]], 
                          msa_index = mutagen_gnomad[[5]], stringsAsFactors = FALSE) %>% unique()




mousecelegans1<-mousecelegans
mousecelegans1$human_Refseq<-indexed$human_refseq[match(mousecelegans1$msa_index, indexed$index)]
colnames(mousecelegans1)[3:4]<-c("Refseq_ID","aapos")
mousecelegans1_m<-merge(mousecelegans1, celegans_wconservation, by = c("Refseq_ID","aapos"))
# mousecelegans1_m<-mousecelegans1_m[!is.na(mousecelegans1_m$mouse_aa),]
colnames(mousecelegans1_m)[1:4]<-c("celegans_Refseq","celegans_aa","Refseq_ID","aapos")
mousecelegans2_m<-merge(mousecelegans1_m, mouse_wconservation, by = c("Refseq_ID","aapos"))

ort_mousecelegans<-mousecelegans2_m[mousecelegans2_m$from.x == mousecelegans2_m$from.y,]
ort_mousecelegans<-ort_mousecelegans[ort_mousecelegans$conservation.x == ort_mousecelegans$conservation.y,]

ort_mousecelegans<-ort_mousecelegans[,c(-5,-13,-20)]
colnames(ort_mousecelegans)[c(1,3,5,13,7,17,10)]<-c("Mouse_Protein_ID","C_elegans_Protein_ID","Human_Protein_ID","Mouse_Variant_ID",
                                                    "C_elegans_Variant_ID","Mouse_Clinical_Significance","C_elegans_Clinical_Significance")
ort_mousecelegans<-unique(ort_mousecelegans)



mouseclinvar1<-mouseclinvar
mouseclinvar1$celegans_Refseq<-indexed$celegans_refseq[match(mouseclinvar1$msa_index, indexed$index)]
colnames(mouseclinvar1)[3:4]<-c("Refseq_ID","aapos")
mouseclinvar1_m<-merge(mouseclinvar1, clinvarxy_wconservation, by = c("Refseq_ID","aapos"))

colnames(mouseclinvar1_m)[1:4]<-c("human_Refseq","human_aa","Refseq_ID","aapos")
mouseclinvar2_m<-merge(mouseclinvar1_m, mouse_wconservation, by = c("Refseq_ID","aapos"))

ort_mouseclinvar<-mouseclinvar2_m[mouseclinvar2_m$from.x == mouseclinvar2_m$from.y,]
ort_mouseclinvar<-ort_mouseclinvar[ort_mouseclinvar$conservation.x == ort_mouseclinvar$conservation.y,]
ort_mouseclinvar<-ort_mouseclinvar[,c(-15,-23)]
ort_mouseclinvar<-unique(ort_mouseclinvar)



mousedbsnp1<-mousedbsnp
mousedbsnp1$celegans_Refseq<-indexed$celegans_refseq[match(mousedbsnp1$msa_index, indexed$index)]
colnames(mousedbsnp1)[3:4]<-c("Refseq_ID","aapos")
mousedbsnp1_m<-merge(mousedbsnp1, dbsnp_wconservation, by = c("Refseq_ID","aapos"))
colnames(mousedbsnp1_m)[1:4]<-c("human_Refseq","human_aa","Refseq_ID","aapos")
mousedbsnp2_m<-merge(mousedbsnp1_m, mouse_wconservation, by = c("Refseq_ID","aapos"))

ort_mousedbsnp<-mousedbsnp2_m[mousedbsnp2_m$from.x == mousedbsnp2_m$from.y,]
ort_mousedbsnp<-ort_mousedbsnp[ort_mousedbsnp$conservation.x == ort_mousedbsnp$conservation.y,]
ort_mousedbsnp<-ort_mousedbsnp[,c(-15,-23)]



mousecosmic1<-mousecosmic
mousecosmic1$celegans_Refseq<-indexed$celegans_refseq[match(mousecosmic1$msa_index, indexed$index)]
colnames(mousecosmic1)[3:4]<-c("Refseq_ID","aapos")
mousecosmic1_m<-merge(mousecosmic1, cosmic_wconservation, by = c("Refseq_ID","aapos"))
colnames(mousecosmic1_m)[1:4]<-c("human_Refseq","human_aa","Refseq_ID","aapos")
mousecosmic2_m<-merge(mousecosmic1_m, mouse_wconservation, by = c("Refseq_ID","aapos"))

ort_mousecosmic<-mousecosmic2_m[mousecosmic2_m$from.x == mousecosmic2_m$from.y,]
ort_mousecosmicx<-ort_mousecosmic[ort_mousecosmic$conservation.x == ort_mousecosmic$conservation.y,]
ort_mousecosmic<-ort_mousecosmic[,c(-15,-23)]




mousegnomad1<-mousegnomad
mousegnomad1$celegans_Refseq<-indexed$celegans_refseq[match(mousegnomad1$msa_index, indexed$index)]
colnames(mousegnomad1)[3:4]<-c("Refseq_ID","aapos")
mousegnomad1_m<-merge(mousegnomad1, gnomad_wconservation, by = c("Refseq_ID","aapos"))
colnames(mousegnomad1_m)[1:4]<-c("human_Refseq","human_aa","Refseq_ID","aapos")
mousegnomad2_m<-merge(mousegnomad1_m, mouse_wconservation, by = c("Refseq_ID","aapos"))

ort_mousegnomad<-mousegnomad2_m[mousegnomad2_m$from.x == mousegnomad2_m$from.y,]
ort_mousegnomad<-ort_mousegnomad[ort_mousegnomad$conservation.x == ort_mousegnomad$conservation.y,]
ort_mousegnomad<-ort_mousegnomad[,c(-15,-23)]
ort_mousegnomad<-unique(ort_mousegnomad)




ort_mousehuman<-rbind(ort_mouseclinvar, ort_mousedbsnp, ort_mousecosmic, ort_mousegnomad)
ort_mousehuman$Human_Variation<-paste0(ort_mousehuman$from.x, ort_mousehuman$human_aa, ort_mousehuman$to.x)
ort_mousehuman$Mouse_Variation<-paste0(ort_mousehuman$from.y, ort_mousehuman$aapos, ort_mousehuman$to.y)
ort_mousehuman<-ort_mousehuman[,c(-2,-4,-7,-10,-15,-19,-21)]
colnames(ort_mousehuman)[c(1,2,4,5,6,7,8,9,10,11,12,13,14,15)]<-c("Mouse_Protein_ID","Human_Protein_ID","C_elegans_Protein_ID","Human_Gene_name",
                                                            "Human_Variant_ID","Human_Clinical_Significance","Human_Disease","Allele_Frequency",
                                                            "Human_Variant_Source","Mouse_Variant_ID","Mouse_Gene_Name","Mouse_Phenotype",
                                                            "Mouse_Significance","Mouse_Variant_Source")
ort_mousehuman<-select(ort_mousehuman, Human_Protein_ID, Mouse_Protein_ID, C_elegans_Protein_ID, Human_Variation, Mouse_Variation,
                       Human_Clinical_Significance, Mouse_Significance, Human_Variant_ID, Mouse_Variant_ID, Human_Gene_name,
                       Mouse_Gene_Name, Human_Disease, Mouse_Phenotype, Allele_Frequency, Human_Variant_Source, 
                       Mouse_Variant_Source, msa_index)



# Mutagen #

mutagencosmic1<-mutagencosmic
mutagencosmic1$celegans_Refseq<-indexed$celegans_refseq[match(mutagencosmic1$msa_index, indexed$index)]
colnames(mutagencosmic1)[3:4]<-c("Refseq_ID","aapos")
mutagencosmic1_m<-merge(mutagencosmic1, cosmic_wconservation, by = c("Refseq_ID","aapos"))
colnames(mutagencosmic1_m)[1:4]<-c("human_Refseq","human_aa","Refseq_ID","aapos")
mutagencosmic2_m<-merge(mutagencosmic1_m, mutagen_wconservation, by = c("Refseq_ID","aapos"))

ort_mutagencosmic<-mutagencosmic2_m[mutagencosmic2_m$from.x == mutagencosmic2_m$from.y,]
ort_mutagencosmic<-ort_mutagencosmic[ort_mutagencosmic$conservation.x == ort_mutagencosmic$conservation.y,]
ort_mutagencosmic<-ort_mutagencosmic[,c(-15,-23)]
ort_mutagencosmic<-unique(ort_mutagencosmic)



mutagengnomad1<-mutagengnomad
mutagengnomad1$celegans_Refseq<-indexed$celegans_refseq[match(mutagengnomad1$msa_index, indexed$index)]
colnames(mutagengnomad1)[3:4]<-c("Refseq_ID","aapos")
mutagengnomad1_m<-merge(mutagengnomad1, gnomad_wconservation, by = c("Refseq_ID","aapos"))
colnames(mutagengnomad1_m)[1:4]<-c("human_Refseq","human_aa","Refseq_ID","aapos")
mutagengnomad2_m<-merge(mutagengnomad1_m, mutagen_wconservation, by = c("Refseq_ID","aapos"))

ort_mutagengnomad<-mutagengnomad2_m[mutagengnomad2_m$from.x == mutagengnomad2_m$from.y,]
ort_mutagengnomad<-ort_mutagengnomad[ort_mutagengnomad$conservation.x == ort_mutagengnomad$conservation.y,]
ort_mutagengnomad<-ort_mutagengnomad[,c(-15,-23)]
ort_mutagengnomad<-unique(ort_mutagengnomad)



mutagendbsnp1<-mutagendbsnp
mutagendbsnp1$celegans_Refseq<-indexed$celegans_refseq[match(mutagendbsnp1$msa_index, indexed$index)]
colnames(mutagendbsnp1)[3:4]<-c("Refseq_ID","aapos")
mutagendbsnp1_m<-merge(mutagendbsnp1, dbsnp_wconservation, by = c("Refseq_ID","aapos"))
colnames(mutagendbsnp1_m)[1:4]<-c("human_Refseq","human_aa","Refseq_ID","aapos")
mutagendbsnp2_m<-merge(mutagendbsnp1_m, mutagen_wconservation, by = c("Refseq_ID","aapos"))

ort_mutagendbsnp<-mutagendbsnp2_m[mutagendbsnp2_m$from.x == mutagendbsnp2_m$from.y,]
ort_mutagendbsnp<-ort_mutagendbsnp[ort_mutagendbsnp$conservation.x == ort_mutagendbsnp$conservation.y,]
ort_mutagendbsnp<-ort_mutagendbsnp[,c(-15,-23)]
ort_mutagendbsnp<-unique(ort_mutagendbsnp)



mutagenclinvar1<-mutagenclinvar
mutagenclinvar1$celegans_Refseq<-indexed$celegans_refseq[match(mutagenclinvar1$msa_index, indexed$index)]
colnames(mutagenclinvar1)[3:4]<-c("Refseq_ID","aapos")
mutagenclinvar1_m<-merge(mutagenclinvar1, clinvarxy_wconservation, by = c("Refseq_ID","aapos"))
colnames(mutagenclinvar1_m)[1:4]<-c("human_Refseq","human_aa","Refseq_ID","aapos")
mutagenclinvar2_m<-merge(mutagenclinvar1_m, mutagen_wconservation, by = c("Refseq_ID","aapos"))

ort_mutagenclinvar<-mutagenclinvar2_m[mutagenclinvar2_m$from.x == mutagenclinvar2_m$from.y,]
ort_mutagenclinvar<-ort_mutagenclinvar[ort_mutagenclinvar$conservation.x == ort_mutagenclinvar$conservation.y,]
ort_mutagenclinvar<-ort_mutagenclinvar[,c(-15,-23)]
ort_mutagenclinvar<-unique(ort_mutagenclinvar)




mutagencelegans1<-mutagencelegans
mutagencelegans1$human_Refseq<-indexed$human_refseq[match(mutagencelegans1$msa_index, indexed$index)]
colnames(mutagencelegans1)[1:2]<-c("Refseq_ID","aapos")
mutagencelegans1_m<-merge(mutagencelegans1, mutagen, by = c("Refseq_ID","aapos"))
colnames(mutagencelegans1_m)[1:4]<-c("mouse_Refseq","mouse_aa","Refseq_ID","aapos")
mutagencelegans2_m<-merge(mutagencelegans1_m, celegans, by = c("Refseq_ID","aapos"))

ort_mutagencelegans<-mutagencelegans2_m[mutagencelegans2_m$from.x == mutagencelegans2_m$from.y,]



ort_mutagenhuman<-rbind(ort_mutagenclinvar, ort_mutagendbsnp, ort_mutagencosmic, ort_mutagengnomad)
ort_mutagenhuman$Human_Variation<-paste0(ort_mutagenhuman$from.x, ort_mutagenhuman$human_aa, ort_mutagenhuman$to.x)
ort_mutagenhuman$Mouse_Variation<-paste0(ort_mutagenhuman$from.y, ort_mutagenhuman$aapos, ort_mutagenhuman$to.y)
ort_mutagenhuman<-ort_mutagenhuman[,c(-2,-4,-7,-10,-15,-17,-21)]
colnames(ort_mutagenhuman)[c(1,2,4,5,6,7,8,9,10,11,12,13,14,15)]<-c("Mouse_Protein_ID","Human_Protein_ID","C_elegans_Protein_ID","Human_Gene_name",
                                                              "Human_Variant_ID","Human_Clinical_Significance","Human_Disease","Allele_Frequency",
                                                              "Human_Variant_Source","Mouse_Gene_Name","Mouse_Significance","Mouse_Phenotype",
                                                              "Mouse_Variant_ID","Mouse_Variant_Source")
ort_mutagenhuman<-select(ort_mutagenhuman, Human_Protein_ID, Mouse_Protein_ID, C_elegans_Protein_ID, Human_Variation, Mouse_Variation,
                         Human_Clinical_Significance, Mouse_Significance, Human_Variant_ID, Mouse_Variant_ID, Human_Gene_name,
                         Mouse_Gene_Name, Human_Disease, Mouse_Phenotype, Allele_Frequency, Human_Variant_Source, 
                         Mouse_Variant_Source, msa_index)


ort_mouseallhuman<-rbind(ort_mousehuman, ort_mutagenhuman)
ort_mouseallhuman1<-ort_mouseallhuman
ort_mouseallhuman1$Human_Clinical_Significance[ort_mouseallhuman1$Human_Clinical_Significance == "PATHOGENIC"]<-"Unknown"

ort_mouseallhuman2<-gnameConverter(ort_mouseallhuman1, "Human_Gene_name")
ort_mouseallhuman2<-unique(ort_mouseallhuman2)
ort_mouseallhuman3<-data.table(ort_mouseallhuman2)
ort_mouseallhuman3<-setcolorder(ort_mouseallhuman3, c(2:10))


# Human - C. elegans Orthologous Variants
clinvarcelegans1<-clinvarcelegans
clinvarcelegans1$mouse_Refseq<-indexed$mouse_refseq[match(clinvarcelegans1$msa_index, indexed$index)]
colnames(clinvarcelegans1)[3:4]<-c("Refseq_ID","aapos")
clinvarcelegans1_m<-merge(clinvarcelegans1, celegans_wconservation, by = c("Refseq_ID","aapos"))
colnames(clinvarcelegans1_m)[1:4]<-c("celegans_Refseq","celegans_aa","Refseq_ID","aapos")
clinvarcelegans2_m<-merge(clinvarcelegans1_m, clinvarxy_wconservation, by = c("Refseq_ID","aapos"))

ort_clinvarcelegans<-clinvarcelegans2_m[clinvarcelegans2_m$from.x == clinvarcelegans2_m$from.y,]
ort_clinvarcelegans<-ort_clinvarcelegans[ort_clinvarcelegans$conservation.x == ort_clinvarcelegans$conservation.y,]
ort_clinvarcelegans<-ort_clinvarcelegans[,c(-13,-23)]
ort_clinvarcelegans<-unique(ort_clinvarcelegans)


dbsnpcelegans1<-dbsnpcelegans
dbsnpcelegans1$mouse_Refseq<-indexed$mouse_refseq[match(dbsnpcelegans1$msa_index, indexed$index)]
colnames(dbsnpcelegans1)[3:4]<-c("Refseq_ID","aapos")
dbsnpcelegans1_m<-merge(dbsnpcelegans1, celegans_wconservation, by = c("Refseq_ID","aapos"))
colnames(dbsnpcelegans1_m)[1:4]<-c("celegans_Refseq","celegans_aa","Refseq_ID","aapos")
dbsnpcelegans2_m<-merge(dbsnpcelegans1_m, dbsnp_wconservation, by = c("Refseq_ID","aapos"))

ort_dbsnpcelegans<-dbsnpcelegans2_m[dbsnpcelegans2_m$from.x == dbsnpcelegans2_m$from.y,]
ort_dbsnpcelegans<-ort_dbsnpcelegans[ort_dbsnpcelegans$conservation.x == ort_dbsnpcelegans$conservation.y,]
ort_dbsnpcelegans<-ort_dbsnpcelegans[,c(-13,-23)]
ort_dbsnpcelegans<-unique(ort_dbsnpcelegans)



cosmiccelegans1<-cosmiccelegans
cosmiccelegans1$mouse_Refseq<-indexed$mouse_refseq[match(cosmiccelegans1$msa_index, indexed$index)]
colnames(cosmiccelegans1)[3:4]<-c("Refseq_ID","aapos")
cosmiccelegans1_m<-merge(cosmiccelegans1, celegans_wconservation, by = c("Refseq_ID","aapos"))
colnames(cosmiccelegans1_m)[1:4]<-c("celegans_Refseq","celegans_aa","Refseq_ID","aapos")
cosmiccelegans2_m<-merge(cosmiccelegans1_m, cosmic_wconservation, by = c("Refseq_ID","aapos"))

ort_cosmiccelegans<-cosmiccelegans2_m[cosmiccelegans2_m$from.x == cosmiccelegans2_m$from.y,]
ort_cosmiccelegans<-ort_cosmiccelegans[ort_cosmiccelegans$conservation.x == ort_cosmiccelegans$conservation.y,]
ort_cosmiccelegans<-ort_cosmiccelegans[,c(-13,-23)]
ort_cosmiccelegans<-unique(ort_cosmiccelegans)



gnomadcelegans1<-gnomadcelegans
gnomadcelegans1$mouse_Refseq<-indexed$mouse_refseq[match(gnomadcelegans1$msa_index, indexed$index)]
colnames(gnomadcelegans1)[3:4]<-c("Refseq_ID","aapos")
gnomadcelegans1_m<-merge(gnomadcelegans1, celegans_wconservation, by = c("Refseq_ID","aapos"))
colnames(gnomadcelegans1_m)[1:4]<-c("celegans_Refseq","celegans_aa","Refseq_ID","aapos")
gnomadcelegans2_m<-merge(gnomadcelegans1_m, gnomad_wconservation, by = c("Refseq_ID","aapos"))

ort_gnomadcelegans<-gnomadcelegans2_m[gnomadcelegans2_m$from.x == gnomadcelegans2_m$from.y,]
ort_gnomadcelegans<-ort_gnomadcelegans[ort_gnomadcelegans$conservation.x == ort_gnomadcelegans$conservation.y,]
ort_gnomadcelegans<-ort_gnomadcelegans[,c(-13,-23)]
ort_gnomadcelegans<-unique(ort_gnomadcelegans)



ort_humancelegans<-rbind(ort_clinvarcelegans, ort_dbsnpcelegans, ort_cosmiccelegans, ort_gnomadcelegans)
ort_humancelegans<-unique(ort_humancelegans)

ort_humancelegans$Human_Variation<-paste0(ort_humancelegans$from.y, ort_humancelegans$aapos, ort_humancelegans$to.y)
ort_humancelegans$C_elegans_Variation<-paste0(ort_humancelegans$from.x, ort_humancelegans$celegans_aa, ort_humancelegans$to.x)
ort_humancelegans<-ort_humancelegans[,c(-2,-4,-7,-10,-14,-17)]
colnames(ort_humancelegans)[c(1,2,4,5,6,7,8,9,10,11,12,13,14,15)]<-c("Human_Protein_ID","C_elegans_Protein_ID","Mouse_Protein_ID","C_elegans_Variant_ID",
                                                                     "C_elegans_Gene_Name","C_elegans_Significance","C_elegans_Phenotype","C_elegans_Variant_Source","Human_Gene_Name",
                                                                     "Human_Variant_ID","Human_Clinical_Significance","Human_Disease","Allele_Frequency","Human_Variant_Source")
ort_humancelegans<-select(ort_humancelegans, Human_Protein_ID, C_elegans_Protein_ID, Mouse_Protein_ID, Human_Variation, C_elegans_Variation,
                          Human_Clinical_Significance, C_elegans_Significance, Human_Variant_ID, C_elegans_Variant_ID, Human_Gene_Name,
                          C_elegans_Gene_Name, Human_Disease, C_elegans_Phenotype, Allele_Frequency, Human_Variant_Source, C_elegans_Variant_Source, msa_index)

ort_humancelegans1<-ort_humancelegans
ort_humancelegans1$Human_Clinical_Significance[ort_humancelegans1$Human_Clinical_Significance == "PATHOGENIC"]<-"Unknown"
ort_humancelegans2<-ort_humancelegans1
ort_humancelegans2$C_elegans_Phenotype<-gsub(pattern = '""', "", ort_humancelegans1$C_elegans_Phenotype)
ort_humancelegans2<-gnameConverter(ort_humancelegans2, "Human_Gene_Name")
ort_humancelegans3<-data.table(ort_humancelegans2)
ort_humancelegans3<-setcolorder(ort_humancelegans3, c(2:10))



# Human - Mouse Double orthologous Variants

d_mutagengnomad<-readRDS("mutagen_gnomad_double_09.12.20.RDS")
d_mutagencosmic<-readRDS("mutagen_cosmic_double_09.12.20.RDS")
d_mutagenclinvar<-readRDS("mutagen_clinvar_double_09.12.20.RDS")
d_mutagendbsnp<-readRDS("mutagen_dbsnp_double_09.12.20.RDS")

d_mousegnomad<-readRDS("mouse_gnomad_double_09.12.20.RDS")
d_mousecosmic<-readRDS("mouse_cosmic_double_09.12.20.RDS")
d_mouseclinvar<-readRDS("mouse_clinvar_double_09.12.20.RDS")
d_mousedbsnp<-readRDS("mouse_dbsnp_double_09.12.20.RDS")

d_gnomadcelegans<-readRDS("gnomad_celegans_double_10.12.20.RDS")
d_clinvarcelegans<-readRDS("clinvar_celegans_double_10.12.20.RDS")
d_dbsnpcelegans<-readRDS("dbsnp_celegans_double_10.12.20.RDS")
d_cosmiccelegans<-readRDS("cosmic_celegans_double_10.12.20.RDS")


mousedbsnp<-data.frame(mouse_Refseq = d_mousedbsnp[[2]], mouse_aa = d_mousedbsnp[[1]], 
                       human_Refseq = d_mousedbsnp[[4]], human_aa = d_mousedbsnp[[3]], 
                       msa_index = d_mousedbsnp[[5]], stringsAsFactors = FALSE) %>% unique()

mouseclinvar<-data.frame(mouse_Refseq = d_mouseclinvar[[2]], mouse_aa = d_mouseclinvar[[1]], 
                       human_Refseq = d_mouseclinvar[[4]], human_aa = d_mouseclinvar[[3]], 
                       msa_index = d_mouseclinvar[[5]], stringsAsFactors = FALSE) %>% unique()

mousegnomad<-data.frame(mouse_Refseq = d_mousegnomad[[2]], mouse_aa = d_mousegnomad[[1]], 
                       human_Refseq = d_mousegnomad[[4]], human_aa = d_mousegnomad[[3]], 
                       msa_index = d_mousegnomad[[5]], stringsAsFactors = FALSE) %>% unique()

mousecosmic<-data.frame(mouse_Refseq = d_mousecosmic[[2]], mouse_aa = d_mousecosmic[[1]], 
                       human_Refseq = d_mousecosmic[[4]], human_aa = d_mousecosmic[[3]], 
                       msa_index = d_mousecosmic[[5]], stringsAsFactors = FALSE) %>% unique()

mutagendbsnp<-data.frame(mouse_Refseq = d_mutagendbsnp[[2]], mouse_aa = d_mutagendbsnp[[1]], 
                       human_Refseq = d_mutagendbsnp[[4]], human_aa = d_mutagendbsnp[[3]], 
                       msa_index = d_mutagendbsnp[[5]], stringsAsFactors = FALSE) %>% unique()

mutagenclinvar<-data.frame(mouse_Refseq = d_mutagenclinvar[[2]], mouse_aa = d_mutagenclinvar[[1]], 
                         human_Refseq = d_mutagenclinvar[[4]], human_aa = d_mutagenclinvar[[3]], 
                         msa_index = d_mutagenclinvar[[5]], stringsAsFactors = FALSE) %>% unique()

mutagengnomad<-data.frame(mouse_Refseq = d_mutagengnomad[[2]], mouse_aa = d_mutagengnomad[[1]], 
                         human_Refseq = d_mutagengnomad[[4]], human_aa = d_mutagengnomad[[3]], 
                         msa_index = d_mutagengnomad[[5]], stringsAsFactors = FALSE) %>% unique()

mutagencosmic<-data.frame(mouse_Refseq = d_mutagencosmic[[2]], mouse_aa = d_mutagencosmic[[1]], 
                         human_Refseq = d_mutagencosmic[[4]], human_aa = d_mutagencosmic[[3]], 
                         msa_index = d_mutagencosmic[[5]], stringsAsFactors = FALSE) %>% unique()


clinvarcelegans<-data.frame(human_Refseq = d_clinvarcelegans[[2]], human_aa = as.numeric(d_clinvarcelegans[[1]]), 
                            celegans_Refseq = d_clinvarcelegans[[4]], celegans_aa = d_clinvarcelegans[[3]], 
                            msa_index = d_clinvarcelegans[[5]], stringsAsFactors = FALSE) %>% unique()

dbsnpcelegans<-data.frame(human_Refseq = d_dbsnpcelegans[[2]], human_aa = as.numeric(d_dbsnpcelegans[[1]]), 
                            celegans_Refseq = d_dbsnpcelegans[[4]], celegans_aa = d_dbsnpcelegans[[3]], 
                            msa_index = d_dbsnpcelegans[[5]], stringsAsFactors = FALSE) %>% unique()

cosmiccelegans<-data.frame(human_Refseq = d_cosmiccelegans[[2]], human_aa = as.numeric(d_cosmiccelegans[[1]]), 
                            celegans_Refseq = d_cosmiccelegans[[4]], celegans_aa = d_cosmiccelegans[[3]], 
                            msa_index = d_cosmiccelegans[[5]], stringsAsFactors = FALSE) %>% unique()

gnomadcelegans<-data.frame(human_Refseq = d_gnomadcelegans[[2]], human_aa = as.numeric(d_gnomadcelegans[[1]]), 
                            celegans_Refseq = d_gnomadcelegans[[4]], celegans_aa = d_gnomadcelegans[[3]], 
                            msa_index = d_gnomadcelegans[[5]], stringsAsFactors = FALSE) %>% unique()



mouseclinvar1<-mouseclinvar
colnames(mouseclinvar1)[3:4]<-c("Refseq_ID","aapos")
mouseclinvar1_m<-merge(mouseclinvar1, clinvarxy_wconservation, by = c("Refseq_ID","aapos"))

colnames(mouseclinvar1_m)[1:4]<-c("human_Refseq","human_aa","Refseq_ID","aapos")
mouseclinvar2_m<-merge(mouseclinvar1_m, mouse_wconservation, by = c("Refseq_ID","aapos"))

ort_mouseclinvar<-mouseclinvar2_m[mouseclinvar2_m$from.x == mouseclinvar2_m$from.y,]
ort_mouseclinvar<-ort_mouseclinvar[ort_mouseclinvar$conservation.x == ort_mouseclinvar$conservation.y,]
ort_mouseclinvar<-ort_mouseclinvar[,c(-14,-22)]
ort_mouseclinvar<-unique(ort_mouseclinvar)


mousedbsnp1<-mousedbsnp
colnames(mousedbsnp1)[3:4]<-c("Refseq_ID","aapos")
mousedbsnp1_m<-merge(mousedbsnp1, dbsnp_wconservation, by = c("Refseq_ID","aapos"))
colnames(mousedbsnp1_m)[1:4]<-c("human_Refseq","human_aa","Refseq_ID","aapos")
mousedbsnp2_m<-merge(mousedbsnp1_m, mouse_wconservation, by = c("Refseq_ID","aapos"))

ort_mousedbsnp<-mousedbsnp2_m[mousedbsnp2_m$from.x == mousedbsnp2_m$from.y,]
ort_mousedbsnp<-ort_mousedbsnp[ort_mousedbsnp$conservation.x == ort_mousedbsnp$conservation.y,]
ort_mousedbsnp<-ort_mousedbsnp[,c(-14,-22)]


mousecosmic1<-mousecosmic
colnames(mousecosmic1)[3:4]<-c("Refseq_ID","aapos")
mousecosmic1_m<-merge(mousecosmic1, cosmic_wconservation, by = c("Refseq_ID","aapos"))
colnames(mousecosmic1_m)[1:4]<-c("human_Refseq","human_aa","Refseq_ID","aapos")
mousecosmic2_m<-merge(mousecosmic1_m, mouse_wconservation, by = c("Refseq_ID","aapos"))

ort_mousecosmic<-mousecosmic2_m[mousecosmic2_m$from.x == mousecosmic2_m$from.y,]
ort_mousecosmic<-ort_mousecosmic[ort_mousecosmic$conservation.x == ort_mousecosmic$conservation.y,]
ort_mousecosmic<-ort_mousecosmic[,c(-14,-22)]


mousegnomad1<-mousegnomad
colnames(mousegnomad1)[3:4]<-c("Refseq_ID","aapos")
mousegnomad1_m<-merge(mousegnomad1, gnomad_wconservation, by = c("Refseq_ID","aapos"))
colnames(mousegnomad1_m)[1:4]<-c("human_Refseq","human_aa","Refseq_ID","aapos")
mousegnomad2_m<-merge(mousegnomad1_m, mouse_wconservation, by = c("Refseq_ID","aapos"))

ort_mousegnomad<-mousegnomad2_m[mousegnomad2_m$from.x == mousegnomad2_m$from.y,]
ort_mousegnomad<-ort_mousegnomad[ort_mousegnomad$conservation.x == ort_mousegnomad$conservation.y,]
ort_mousegnomad<-ort_mousegnomad[,c(-14,-22)]
ort_mousegnomad<-unique(ort_mousegnomad)



ort_mousehuman<-rbind(ort_mouseclinvar, ort_mousedbsnp, ort_mousecosmic, ort_mousegnomad)
ort_mousehuman$Human_Variation<-paste0(ort_mousehuman$from.x, ort_mousehuman$human_aa, ort_mousehuman$to.x)
ort_mousehuman$Mouse_Variation<-paste0(ort_mousehuman$from.y, ort_mousehuman$aapos, ort_mousehuman$to.y)
ort_mousehuman<-ort_mousehuman[,c(-2,-4,-6,-9,-14,-18,-20)]
colnames(ort_mousehuman)[c(1,2,4,5,6,7,8,9,10,11,12,13,14)]<-c("Mouse_Protein_ID","Human_Protein_ID","Human_Gene_name",
                                                             "Human_Variant_ID","Human_Clinical_Significance","Human_Disease","Allele_Frequency",
                                                             "Human_Variant_Source","Mouse_Variant_ID","Mouse_Gene_Name","Mouse_Phenotype",
                                                             "Mouse_Significance","Mouse_Variant_Source")
ort_mousehuman<-select(ort_mousehuman, Human_Protein_ID, Mouse_Protein_ID, Human_Variation, Mouse_Variation,
                       Human_Clinical_Significance, Mouse_Significance, Human_Variant_ID, Mouse_Variant_ID, Human_Gene_name,
                       Mouse_Gene_Name, Human_Disease, Mouse_Phenotype, Allele_Frequency, 
                       Human_Variant_Source, Mouse_Variant_Source, msa_index)



mutagencosmic1<-mutagencosmic
colnames(mutagencosmic1)[3:4]<-c("Refseq_ID","aapos")
mutagencosmic1_m<-merge(mutagencosmic1, cosmic_wconservation, by = c("Refseq_ID","aapos"))
colnames(mutagencosmic1_m)[1:4]<-c("human_Refseq","human_aa","Refseq_ID","aapos")
mutagencosmic2_m<-merge(mutagencosmic1_m, mutagen_wconservation, by = c("Refseq_ID","aapos"))

ort_mutagencosmic<-mutagencosmic2_m[mutagencosmic2_m$from.x == mutagencosmic2_m$from.y,]
ort_mutagencosmic<-ort_mutagencosmic[ort_mutagencosmic$conservation.x == ort_mutagencosmic$conservation.y,]
ort_mutagencosmic<-ort_mutagencosmic[,c(-14,-22)]
ort_mutagencosmic<-unique(ort_mutagencosmic)


mutagengnomad1<-mutagengnomad
colnames(mutagengnomad1)[3:4]<-c("Refseq_ID","aapos")
mutagengnomad1_m<-merge(mutagengnomad1, gnomad_wconservation, by = c("Refseq_ID","aapos"))
colnames(mutagengnomad1_m)[1:4]<-c("human_Refseq","human_aa","Refseq_ID","aapos")
mutagengnomad2_m<-merge(mutagengnomad1_m, mutagen_wconservation, by = c("Refseq_ID","aapos"))

ort_mutagengnomad<-mutagengnomad2_m[mutagengnomad2_m$from.x == mutagengnomad2_m$from.y,]
ort_mutagengnomad<-ort_mutagengnomad[ort_mutagengnomad$conservation.x == ort_mutagengnomad$conservation.y,]
ort_mutagengnomad<-ort_mutagengnomad[,c(-14,-22)]
ort_mutagengnomad<-unique(ort_mutagengnomad)


mutagendbsnp1<-mutagendbsnp
colnames(mutagendbsnp1)[3:4]<-c("Refseq_ID","aapos")
mutagendbsnp1_m<-merge(mutagendbsnp1, dbsnp_wconservation, by = c("Refseq_ID","aapos"))
colnames(mutagendbsnp1_m)[1:4]<-c("human_Refseq","human_aa","Refseq_ID","aapos")
mutagendbsnp2_m<-merge(mutagendbsnp1_m, mutagen_wconservation, by = c("Refseq_ID","aapos"))

ort_mutagendbsnp<-mutagendbsnp2_m[mutagendbsnp2_m$from.x == mutagendbsnp2_m$from.y,]
ort_mutagendbsnp<-ort_mutagendbsnp[ort_mutagendbsnp$conservation.x == ort_mutagendbsnp$conservation.y,]
ort_mutagendbsnp<-ort_mutagendbsnp[,c(-14,-22)]
ort_mutagendbsnp<-unique(ort_mutagendbsnp)


mutagenclinvar1<-mutagenclinvar
colnames(mutagenclinvar1)[3:4]<-c("Refseq_ID","aapos")
mutagenclinvar1_m<-merge(mutagenclinvar1, clinvarxy_wconservation, by = c("Refseq_ID","aapos"))
colnames(mutagenclinvar1_m)[1:4]<-c("human_Refseq","human_aa","Refseq_ID","aapos")
mutagenclinvar2_m<-merge(mutagenclinvar1_m, mutagen_wconservation, by = c("Refseq_ID","aapos"))

ort_mutagenclinvar<-mutagenclinvar2_m[mutagenclinvar2_m$from.x == mutagenclinvar2_m$from.y,]
ort_mutagenclinvar<-ort_mutagenclinvar[ort_mutagenclinvar$conservation.x == ort_mutagenclinvar$conservation.y,]
ort_mutagenclinvar<-ort_mutagenclinvar[,c(-14,-22)]
ort_mutagenclinvar<-unique(ort_mutagenclinvar)


ort_mutagenhuman<-rbind(ort_mutagenclinvar, ort_mutagendbsnp, ort_mutagencosmic, ort_mutagengnomad)
ort_mutagenhuman$Human_Variation<-paste0(ort_mutagenhuman$from.x, ort_mutagenhuman$human_aa, ort_mutagenhuman$to.x)
ort_mutagenhuman$Mouse_Variation<-paste0(ort_mutagenhuman$from.y, ort_mutagenhuman$aapos, ort_mutagenhuman$to.y)
ort_mutagenhuman<-ort_mutagenhuman[,c(-2,-4,-6,-9,-14,-16,-20)]
colnames(ort_mutagenhuman)[c(1,2,4,5,6,7,8,9,10,11,12,13,14)]<-c("Mouse_Protein_ID","Human_Protein_ID","Human_Gene_name",
                                                            "Human_Variant_ID","Human_Clinical_Significance","Human_Disease","Allele_Frequency",
                                                            "Human_Variant_Source","Mouse_Gene_Name","Mouse_Significance","Mouse_Phenotype",
                                                            "Mouse_Variant_ID","Mouse_Variant_Source")
ort_mutagenhuman<-select(ort_mutagenhuman, Human_Protein_ID, Mouse_Protein_ID, Human_Variation, Mouse_Variation,
                         Human_Clinical_Significance, Mouse_Significance, Human_Variant_ID, Mouse_Variant_ID, Human_Gene_name,
                         Mouse_Gene_Name, Human_Disease, Mouse_Phenotype, Allele_Frequency, 
                         Human_Variant_Source, Mouse_Variant_Source, msa_index)


ort_mouseallhuman<-rbind(ort_mousehuman, ort_mutagenhuman)
ort_mouseallhuman<-unique(ort_mouseallhuman)
ort_mouseallhuman1<-ort_mouseallhuman
ort_mouseallhuman1$Human_Clinical_Significance[ort_mouseallhuman1$Human_Clinical_Significance == "PATHOGENIC"]<-"Unknown"

ort_mouseallhuman2<-gnameConverter(ort_mouseallhuman1, "Human_Gene_name")
ort_mouseallhuman4<-data.table(ort_mouseallhuman2)
ort_mouseallhuman4<-setcolorder(ort_mouseallhuman4, c(2:10))

ort_mouseallhuman4$C_elegans_Protein_ID<-""
ort_mouseallhuman4<-setcolorder(ort_mouseallhuman4, c(1,2,16))

ort_mouseallhuman34<-rbind(ort_mouseallhuman3, ort_mouseallhuman4) %>% unique()

#ort_mouseallhuman34<-unique(setDT(ort_mouseallhuman34), by = c(1,2,4,5,6,7))

# m_ort_mouseallhuman34<-unique(setDT(ort_mouseallhuman34), by = c(2,5,7))
# h_ort_mouseallhuman34<-unique(setDT(ort_mouseallhuman34), by = c(1,4,6))
# 
# length(m_ort_mouseallhuman34$Human_Protein_ID[m_ort_mouseallhuman34$Mouse_Significance == "phenotypic"])
# length(m_ort_mouseallhuman34$Human_Protein_ID[m_ort_mouseallhuman34$Mouse_Significance == "incidental"])

#write.table(ort_mouseallhuman4, "human_mouse_double_orthologous_variations.txt", quote = FALSE, row.names = FALSE, sep = "\t")





# C. elegans

clinvarcelegans1<-clinvarcelegans
colnames(clinvarcelegans1)[3:4]<-c("Refseq_ID","aapos")
clinvarcelegans1_m<-merge(clinvarcelegans1, celegans_wconservation, by = c("Refseq_ID","aapos"))
colnames(clinvarcelegans1_m)[1:4]<-c("celegans_Refseq","celegans_aa","Refseq_ID","aapos")
clinvarcelegans2_m<-merge(clinvarcelegans1_m, clinvarxy_wconservation, by = c("Refseq_ID","aapos"))

ort_clinvarcelegans<-clinvarcelegans2_m[clinvarcelegans2_m$from.x == clinvarcelegans2_m$from.y,]
ort_clinvarcelegans<-ort_clinvarcelegans[ort_clinvarcelegans$conservation.x == ort_clinvarcelegans$conservation.y,]
ort_clinvarcelegans<-ort_clinvarcelegans[,c(-12,-22)]
ort_clinvarcelegans<-unique(ort_clinvarcelegans)



dbsnpcelegans1<-dbsnpcelegans
colnames(dbsnpcelegans1)[3:4]<-c("Refseq_ID","aapos")
dbsnpcelegans1_m<-merge(dbsnpcelegans1, celegans_wconservation, by = c("Refseq_ID","aapos"))
colnames(dbsnpcelegans1_m)[1:4]<-c("celegans_Refseq","celegans_aa","Refseq_ID","aapos")
dbsnpcelegans2_m<-merge(dbsnpcelegans1_m, dbsnp_wconservation, by = c("Refseq_ID","aapos"))

ort_dbsnpcelegans<-dbsnpcelegans2_m[dbsnpcelegans2_m$from.x == dbsnpcelegans2_m$from.y,]
ort_dbsnpcelegans<-ort_dbsnpcelegans[ort_dbsnpcelegans$conservation.x == ort_dbsnpcelegans$conservation.y,]
ort_dbsnpcelegans<-ort_dbsnpcelegans[,c(-12,-22)]
ort_dbsnpcelegans<-unique(ort_dbsnpcelegans)



cosmiccelegans1<-cosmiccelegans
colnames(cosmiccelegans1)[3:4]<-c("Refseq_ID","aapos")
cosmiccelegans1_m<-merge(cosmiccelegans1, celegans_wconservation, by = c("Refseq_ID","aapos"))
colnames(cosmiccelegans1_m)[1:4]<-c("celegans_Refseq","celegans_aa","Refseq_ID","aapos")
cosmiccelegans2_m<-merge(cosmiccelegans1_m, cosmic_wconservation, by = c("Refseq_ID","aapos"))

ort_cosmiccelegans<-cosmiccelegans2_m[cosmiccelegans2_m$from.x == cosmiccelegans2_m$from.y,]
ort_cosmiccelegans<-ort_cosmiccelegans[ort_cosmiccelegans$conservation.x == ort_cosmiccelegans$conservation.y,]
ort_cosmiccelegans<-ort_cosmiccelegans[,c(-12,-22)]
ort_cosmiccelegans<-unique(ort_cosmiccelegans)



gnomadcelegans1<-gnomadcelegans
colnames(gnomadcelegans1)[3:4]<-c("Refseq_ID","aapos")
gnomadcelegans1_m<-merge(gnomadcelegans1, celegans_wconservation, by = c("Refseq_ID","aapos"))
colnames(gnomadcelegans1_m)[1:4]<-c("celegans_Refseq","celegans_aa","Refseq_ID","aapos")
gnomadcelegans2_m<-merge(gnomadcelegans1_m, gnomad_wconservation, by = c("Refseq_ID","aapos"))

ort_gnomadcelegans<-gnomadcelegans2_m[gnomadcelegans2_m$from.x == gnomadcelegans2_m$from.y,]
ort_gnomadcelegans<-ort_gnomadcelegans[ort_gnomadcelegans$conservation.x == ort_gnomadcelegans$conservation.y,]
ort_gnomadcelegans<-ort_gnomadcelegans[,c(-12,-22)]
ort_gnomadcelegans<-unique(ort_gnomadcelegans)




ort_humancelegans<-rbind(ort_clinvarcelegans, ort_dbsnpcelegans, ort_cosmiccelegans, ort_gnomadcelegans)
ort_humancelegans<-unique(ort_humancelegans)

ort_humancelegans$Human_Variation<-paste0(ort_humancelegans$from.y, ort_humancelegans$aapos, ort_humancelegans$to.y)
ort_humancelegans$C_elegans_Variation<-paste0(ort_humancelegans$from.x, ort_humancelegans$celegans_aa, ort_humancelegans$to.x)
ort_humancelegans<-ort_humancelegans[,c(-2,-4,-6,-9,-13,-16)]
colnames(ort_humancelegans)[c(1,2,4,5,6,7,8,9,10,11,12,13,14)]<-c("Human_Protein_ID","C_elegans_Protein_ID","C_elegans_Variant_ID",
                                                               "C_elegans_Gene_Name","C_elegans_Significance","C_elegans_Phenotype","C_elegans_Variant_Source",
                                                            "Human_Gene_Name","Human_Variant_ID","Human_Clinical_Significance","Human_Disease","Allele_Frequency",
                                                            "Human_Variant_Source")
ort_humancelegans<-select(ort_humancelegans, Human_Protein_ID, C_elegans_Protein_ID, Human_Variation, C_elegans_Variation,
                          Human_Clinical_Significance, C_elegans_Significance, Human_Variant_ID, C_elegans_Variant_ID, Human_Gene_Name,
                          C_elegans_Gene_Name, Human_Disease, C_elegans_Phenotype, Allele_Frequency, 
                          Human_Variant_Source, C_elegans_Variant_Source, msa_index)

ort_humancelegans12<-ort_humancelegans
ort_humancelegans12$Human_Clinical_Significance[ort_humancelegans12$Human_Clinical_Significance == "PATHOGENIC"]<-"Unknown"
ort_humancelegans22<-ort_humancelegans12
ort_humancelegans22$C_elegans_Phenotype<-gsub(pattern = '""', "", ort_humancelegans12$C_elegans_Phenotype)
ort_humancelegans22<-gnameConverter(ort_humancelegans22, "Human_Gene_Name")
ort_humancelegans32<-data.table(ort_humancelegans22)
ort_humancelegans32<-setcolorder(ort_humancelegans32, c(2:9))
ort_humancelegans32$Mouse_Protein_ID<-""
ort_humancelegans32<-setcolorder(ort_humancelegans32, c(1,2,17))

ort_humancelegansCombined<-rbind(ort_humancelegans3, ort_humancelegans32) %>% unique()
ort_humancelegansCombined$C_elegans_Phenotype<-gsub("\n", "", ort_humancelegansCombined$C_elegans_Phenotype)




# c_ort_humancelegansCombined<-unique(setDT(ort_humancelegansCombined), by = c(2,5,7))
# length(c_ort_humancelegansCombined$Human_Protein_ID[c_ort_humancelegansCombined$C_elegans_Significance == "Phenotypic"])
# length(c_ort_humancelegansCombined$Human_Protein_ID[c_ort_humancelegansCombined$C_elegans_Significance == "Unknown"])
# length(c_ort_humancelegansCombined$Human_Protein_ID)
# 
# 
# 
# write.table(ort_humancelegans3, "human_celegans_double_orthologous_variations.txt", quote = FALSE, row.names = FALSE, sep = "\t")
# 
# 
# 
# phe_ort_mouseallhuman3<-ort_mouseallhuman3[ort_mouseallhuman3$Mouse_Significance == "phenotypic" | 
#                                              ort_mouseallhuman3$Mouse_Significance == "Phenotypic"]
# uphe_ort_mouseallhuman3<-unique(setDT(phe_ort_mouseallhuman3), by = c("Mouse_Protein_ID","Mouse_Variation","Mouse_Phenotype"))
# 
# 
# phe_ort_humancelegans3<-ort_humancelegans3[ort_humancelegans3$C_elegans_Significance == "phenotypic" | 
#                                              ort_humancelegans3$C_elegans_Significance == "Phenotypic"]
# uphe_ort_humancelegans3<-unique(setDT(phe_ort_humancelegans3), by = c("C_elegans_Protein_ID","C_elegans_Variation","C_elegans_Phenotype"))
