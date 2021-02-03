

library(VennDiagram)
library(tidyverse)

dflist<-list(clinvar, dbsnp, cosmic, gnomad)
x<-map(dflist, ~.$Gene_name) %>% #creates a list just of the column of interest
  reduce(intersect)  #applies the intersect function cumulatively to the list

map(aaconservation, ~.$to)

int_clinvar<-data.frame(Refseq_ID = clinvarxy_wconservation$Refseq_ID, 
                        Variant = paste0(clinvarxy_wconservation$from, clinvarxy_wconservation$aapos, clinvarxy_wconservation$to), 
                        stringsAsFactors = FALSE) %>% unique()

int_dbsnp<-data.frame(Refseq_ID = dbsnp_wconservation$Refseq_ID, 
                        Variant = paste0(dbsnp_wconservation$from, dbsnp_wconservation$aapos, dbsnp_wconservation$to), 
                        stringsAsFactors = FALSE) %>% unique()

int_cosmic<-data.frame(Refseq_ID = cosmic_wconservation$Refseq_ID, 
                      Variant = paste0(cosmic_wconservation$from, cosmic_wconservation$aapos, cosmic_wconservation$to), 
                      stringsAsFactors = FALSE) %>% unique()

int_gnomad<-data.frame(Refseq_ID = gnomad_wconservation$Refseq_ID, 
                      Variant = paste0(gnomad_wconservation$from, gnomad_wconservation$aapos, gnomad_wconservation$to), 
                      stringsAsFactors = FALSE) %>% unique()



# dflist<-list(int_clinvar, int_dbsnp, int_cosmic, int_gnomad)
# 
# df<-setNames(combn(1:3, 2, FUN = function(i) sum(df[i[1],] - df[i[2],]==0)), 
#          combn(rownames(df), 2, toString))


#co<-fintersect(setDT(dflist[xa[1,i]]),setDT(dflist[xa[2,i]]))
  
co1<-dplyr::intersect(int_clinvar,int_dbsnp)
co2<-dplyr::intersect(int_clinvar,int_cosmic)
co3<-dplyr::intersect(int_clinvar,int_gnomad)
co4<-dplyr::intersect(int_cosmic,int_dbsnp)
co5<-dplyr::intersect(int_gnomad,int_dbsnp)
co6<-dplyr::intersect(int_cosmic,int_gnomad)

cot1<-Reduce(dplyr::intersect, list(int_clinvar, int_dbsnp, int_cosmic))
cot2<-Reduce(dplyr::intersect, list(int_clinvar, int_dbsnp, int_gnomad))
cot3<-Reduce(dplyr::intersect, list(int_clinvar, int_cosmic, int_gnomad))
cot4<-Reduce(dplyr::intersect, list(int_dbsnp, int_cosmic, int_gnomad))

coq<-Reduce(dplyr::intersect, list(int_clinvar, int_dbsnp, int_cosmic, int_gnomad))

onlyclinvar<-anti_join(int_clinvar, int_dbsnp, by = c("Refseq_ID","Variant"))
onlyclinvar<-anti_join(onlyclinvar, int_cosmic, by = c("Refseq_ID","Variant"))
onlyclinvar<-anti_join(onlyclinvar, int_gnomad, by = c("Refseq_ID","Variant"))

onlydbsnp<-anti_join(int_dbsnp, int_clinvar, by = c("Refseq_ID","Variant"))
onlydbsnp<-anti_join(onlydbsnp, int_cosmic, by = c("Refseq_ID","Variant"))
onlydbsnp<-anti_join(onlydbsnp, int_gnomad, by = c("Refseq_ID","Variant"))

onlycosmic<-anti_join(int_cosmic, int_dbsnp, by = c("Refseq_ID","Variant"))
onlycosmic<-anti_join(onlycosmic, int_clinvar, by = c("Refseq_ID","Variant"))
onlycosmic<-anti_join(onlycosmic, int_gnomad, by = c("Refseq_ID","Variant"))

onlygnomad<-anti_join(int_gnomad, int_dbsnp, by = c("Refseq_ID","Variant"))
onlygnomad<-anti_join(onlygnomad, int_cosmic, by = c("Refseq_ID","Variant"))
onlygnomad<-anti_join(onlygnomad, int_clinvar, by = c("Refseq_ID","Variant"))



draw.quad.venn(494609,600476,3485814,7435907,52256,56584,131935,92979,197327,578192,13982,23250,26186,43521,6715,
               category = c("ClinVar","dbSNP","COSMIC","gnomAD"), fill = c("darkseagreen2","deepskyblue2","firebrick2","bisque2"))

draw.quad.venn(length(int_clinvar[[1]]),length(int_dbsnp[[1]]),length(int_cosmic[[1]]),length(int_gnomad[[1]]),
               length(co1[[1]]),length(co2[[1]]),length(co3[[1]]),length(co4[[1]]),length(co5[[1]]),length(co6[[1]]),
               length(cot1[[1]]),length(cot2[[1]]),length(cot3[[1]]),length(cot4[[1]]),
               length(coq[[1]]),
               category = c("ClinVar","dbSNP","COSMIC","gnomAD"), lty = "blank", 
               fill = c("skyblue", "pink1", "mediumorchid", "orange"))

