
clinvar_newest<-fread("variant_summary.txt")
clinvar_newest_snp<-clinvar_newest[clinvar_newest$Type == "single nucleotide variant"]
clinvar_newest_snp<-unique(setDT(clinvar_newest_snp), by = "VariationID")

cpatho<-clinvar_newest_snp[grep("^(?=.*pathogenic)(?!.*Conflicting)", clinvar_newest_snp$ClinicalSignificance, perl = TRUE, ignore.case = TRUE)]
cbenig<-clinvar_newest_snp[grep("^(?=.*benign)(?!.*Conflicting)", clinvar_newest_snp$ClinicalSignificance, perl = TRUE, ignore.case = TRUE)]
cconfl<-clinvar_newest_snp[grep("^(?=.*conflicting)", clinvar_newest_snp$ClinicalSignificance, perl = TRUE, ignore.case = TRUE)]
cvus<-clinvar_newest_snp[grep("^(?=.*uncertain)|(?=.*no interpretation)", clinvar_newest_snp$ClinicalSignificance, perl = TRUE, ignore.case = TRUE)]

cother<-anti_join(clinvar_newest_snp, cpatho, by = "VariationID")
cother<-anti_join(cother, cbenig, by = "VariationID")
cother<-anti_join(cother, cconfl, by = "VariationID")
cother<-anti_join(cother, cvus, by = "VariationID")

snp_types<-data.frame(types = c("Pathogenic","Benign","Conflicting","VUS","Others","Total"), 
                      numbers = c(length(cpatho[[1]]),length(cbenig[[1]]),length(cconfl[[1]]),length(cvus[[1]]),length(cother[[1]]),
                      length(clinvar_newest_snp[[1]])))
snp_types$Percentage<-snp_types$numbers/(snp_types$numbers[6])
snp_types$types<-factor(snp_types$types, levels =c("Others","VUS","Conflicting","Benign","Pathogenic","Total"))

snp_types$Variant_type<-"Single Nucleotide Variants"
ggplot(data = snp_types[-6,], aes(x = types, y = Percentage, fill = types)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_brewer(palette="Set3") +
  theme_minimal() + coord_flip() + theme(legend.position = "none", plot.title = element_text(hjust = 0.4, size = 11),
                                         axis.title.x = element_text(size = 10),
                                         axis.title.y = element_text(size = 10)) +
  ggtitle("ClinVar") +
  xlab("Variant Type") + ylab("Ratio in All SNPs")



clinvar_types<-unique(clinvar3$V8)
clinvar_type_numbers<-c()
for (i in clinvar_types){
  cn<-length(clinvar3$V8[clinvar3$V8 == i])
  clinvar_type_numbers<-c(clinvar_type_numbers, cn)
}

clinvar_type_num1<-clinvar_type_num
clinvar_type_num1<-clinvar_type_num1[c(1:5,11:14,18),]
clinvar_type_num$numbers<-as.numeric(clinvar_type_num$numbers)
clin_numb<-sum(clinvar_type_num$numbers[c(6:10,15:17)])
clinvar_type_num1[(length(clinvar_type_num1$types)+1),]<-c("Others",clin_numb)

clinvar_type_num1$Pathogenic<-0
clinvar_type_num1$Benign<-0
clinvar_type_num1$Conflicting<-0
clinvar_type_num1$Vus<-0
clinvar_type_num1$Others<-0

clinvar_type_num1$Pathogenic[clinvar_type_num1$types == "single nucleotide variant"]<-snp_types$numbers[snp_types$types == "Pathogenic"]
clinvar_type_num1$Benign[clinvar_type_num1$types == "single nucleotide variant"]<-snp_types$numbers[snp_types$types == "Benign"]
clinvar_type_num1$Conflicting[clinvar_type_num1$types == "single nucleotide variant"]<-snp_types$numbers[snp_types$types == "Conflicting"]
clinvar_type_num1$Vus[clinvar_type_num1$types == "single nucleotide variant"]<-snp_types$numbers[snp_types$types == "VUS"]
clinvar_type_num1$Others[clinvar_type_num1$types == "single nucleotide variant"]<-snp_types$numbers[snp_types$types == "Others"]

clinvar_type_num2<-pivot_longer(clinvar_type_num1, c(Pathogenic, Benign, Conflicting, Vus, Others), names_to = "SNP_type", values_to = "SNP_numbers")

cpatho<-clinvar3[grep("^(?=.*pathogenic)(?!.*Conflicting)", clinvar3$V10, perl = TRUE, ignore.case = TRUE)]
cbenig<-clinvar3[grep("^(?=.*benign)(?!.*Conflicting)", clinvar3$V10, perl = TRUE, ignore.case = TRUE)]
cconfl<-clinvar3[grep("^(?=.*conflicting)", clinvar3$V10, perl = TRUE, ignore.case = TRUE)]
cvus<-clinvar3[grep("^(?=.*uncertain)|(?=.*no interpretation)", clinvar3$V10, perl = TRUE, ignore.case = TRUE)]

cother<-anti_join(clinvar3, cpatho, by = "variation_id")
cother<-anti_join(cother, cbenig, by = "variation_id")
cother<-anti_join(cother, cconfl, by = "variation_id")
cother<-anti_join(cother, cvus, by = "variation_id")



length(clinvar3$V10[grep("^(?=.*pathogenic)(?!.*Conflicting)", clinvar3$V10, perl = TRUE, ignore.case = TRUE)])
length(clinvar3$V10[grep("^(?=.*benign)(?!.*Conflicting)", clinvar3$V10, perl = TRUE, ignore.case = TRUE)])
length(clinvar3$V10[grep("^(?=.*conflicting)", clinvar3$V10, perl = TRUE, ignore.case = TRUE)])
length(clinvar3$V10[grep("^(?=.*uncertain)|(?=.*no interpretation)", clinvar3$V10, perl = TRUE, ignore.case = TRUE)])

snp_types<-data.frame(types = c("Pathogenic","Benign","Conflicting","VUS","Others","Total"), 
                      numbers = c(107339,172792,21982,270584,11439,584136), 
                      stringsAsFactors = FALSE)
snp_types$Variant_type<-"Single Nucleotide Variants"


clinvar_type_num<-data.frame(types = clinvar_types, numbers = clinvar_type_numbers, stringsAsFactors = FALSE)
clinvar_type_num[18,]<-c("Total",sum(clinvar_type_num$numbers))



var_numbers<-fread("convart_clinvar_variant_numbers.txt")
var_numbers$Variant_type[8]<-"Single Nucleotide Variants"



var_numbers$ConVart_Percent<-var_numbers$ConVarT/(var_numbers$ConVarT + var_numbers$ClinVar)
var_numbers$ClinVar_Percent<-var_numbers$ClinVar/(var_numbers$ConVarT + var_numbers$ClinVar)
var_numbers$Percent<-var_numbers$ClinVar/var_numbers$ClinVar[var_numbers$Variant_type == "TOTAL"]
var_numbersx<-var_numbers
var_numbersx$Percent<-var_numbersx$ConVarT/var_numbersx$ConVarT[var_numbersx$Variant_type == "TOTAL"]

var_numbers1<-pivot_longer(var_numbers[,c(1,4,5)], !Variant_type, names_to = "Source", values_to = "Number")
var_numbers2<-pivot_longer(var_numbers[,c(1,2)], !Variant_type, names_to = "Source", values_to = "Number")
var_numbers$types<-var_numbers$Variant_type



ggplot(data = var_numbers1, aes(x = Variant_type, y = Number, fill = Source)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Number of Variants") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#999999", "#E69F00")) + coord_flip()



ggplot(data = var_numbers2[c(-8,-9),], aes(x = Variant_type, y = Number, fill = Variant_type)) +
  geom_bar(stat = "identity", color = "black") +
  #geom_label(aes(label = Number), size = 1, position = position_stack(vjust = 0.5)) +
  #scale_y_continuous("Number of Variants") +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9","#999999", "#E69F00", "#56B4E9","#999999", "#E69F00"))+
  theme_minimal()
  
# theme(panel.background = element_rect(fill = "white"),
  #       axis.line = element_line(colour = "black"), text = element_text(size=12),
  #       plot.title = element_text(hjust = 0.5)) + coord_flip()

ggplot(data = var_numbersx[-9,], aes(x = Variant_type, y = Percent, fill = Variant_type)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_brewer(palette="Set3") +
  theme_minimal() + coord_flip() + theme(legend.position = "none", plot.title = element_text(hjust = 0.2, size = 11),
                                         axis.title.x = element_text(size = 10),
                                         axis.title.y = element_text(size = 10)) +
  ggtitle("ConVarT") +
  xlab("Variant Type") + ylab("Ratio in All Variants")

ggplot(data = var_numbers[-9,], aes(x = Variant_type, y = Percent, fill = Variant_type)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_brewer(palette="Set3") +
  theme_minimal() + coord_flip() + theme(legend.position = "none", plot.title = element_text(hjust = 0.2, size = 11),
                                         axis.title.x = element_text(size = 10),
                                         axis.title.y = element_text(size = 10)) +
  ggtitle("ClinVar") +
  xlab("Variant Type") + ylab("Ratio in All Variants")



clinvar_type_num2$SNP_numbers[c(1:10,16:55)]<-clinvar_type_num2$numbers[c(1:10,16:55)]
clinvar_type_num2$SNP_type[c(1:10,16:55)]<-clinvar_type_num2$types[c(1:10,16:55)]
clinvar_type_num3<-unique(clinvar_type_num2)
mycolors = c(brewer.pal(name="Set3", n = 5), rep("#cccc00",9))

ggplot(data = clinvar_type_num3, aes(x = types, y = SNP_numbers, fill = SNP_type)) +
  geom_bar(stat = "identity", color = "black") +
  scale_color_manual("mycolors") +
  theme_minimal() + coord_flip() + theme(legend.position = "right", plot.title = element_text(hjust = 0.2, size = 11),
                                         axis.title.x = element_text(size = 10),
                                         axis.title.y = element_text(size = 10)) +
  ggtitle("ClinVar") +
  xlab("Variant Type") + ylab("Ratio in All Variants")
