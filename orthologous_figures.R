
library(ggplot2)
library(patchwork)


# Clinvar figures

# Mouse

ort_mouse_pat<-ort_mouseallhuman34[grep("^(?=.*pathogenic)(?!.*Conflicting)", ort_mouseallhuman34$Human_Clinical_Significance, perl=TRUE, ignore.case = TRUE),]
ort_mouse_ben<-ort_mouseallhuman34[grep("^(?=.*benign)(?!.*Conflicting)", ort_mouseallhuman34$Human_Clinical_Significance, perl=TRUE, ignore.case = TRUE),]
ort_mouse_conf<-ort_mouseallhuman34[grep("^(?=.*conflicting)", ort_mouseallhuman34$Human_Clinical_Significance, perl=TRUE, ignore.case = TRUE),]
ort_mouse_vus<-ort_mouseallhuman34[grep("^(?=.*uncertain)|(?=.*not provided)|(?=.*no interpretation)", ort_mouseallhuman34$Human_Clinical_Significance, perl=TRUE, ignore.case = TRUE),]

ort_pat_phe_m<-ort_mouse_pat[ort_mouse_pat$Mouse_Significance == "phenotypic"]
ort_ben_phe_m<-ort_mouse_ben[ort_mouse_ben$Mouse_Significance == "phenotypic"]
ort_conf_phe_m<-ort_mouse_conf[ort_mouse_conf$Mouse_Significance == "phenotypic"]
ort_vus_phe_m<-ort_mouse_vus[ort_mouse_vus$Mouse_Significance == "phenotypic"]

ort_pat_inc_m<-ort_mouse_pat[ort_mouse_pat$Mouse_Significance == "incidental"]
ort_ben_inc_m<-ort_mouse_ben[ort_mouse_ben$Mouse_Significance == "incidental"]
ort_conf_inc_m<-ort_mouse_conf[ort_mouse_conf$Mouse_Significance == "incidental"]
ort_vus_inc_m<-ort_mouse_vus[ort_mouse_vus$Mouse_Significance == "incidental"]

ort_pat_phe_m<-unique(setDT(ort_pat_phe_m), by = c(1,4,6))
ort_ben_phe_m<-unique(setDT(ort_ben_phe_m), by = c(1,4,6))
ort_conf_phe_m<-unique(setDT(ort_conf_phe_m), by = c(1,4,6))
ort_vus_phe_m<-unique(setDT(ort_vus_phe_m), by = c(1,4,6))

ort_pat_inc_m<-unique(setDT(ort_pat_inc_m), by = c(1,4,6))
ort_ben_inc_m<-unique(setDT(ort_ben_inc_m), by = c(1,4,6))
ort_conf_inc_m<-unique(setDT(ort_conf_inc_m), by = c(1,4,6))
ort_vus_inc_m<-unique(setDT(ort_vus_inc_m), by = c(1,4,6))

pp<-length(ort_pat_phe_m[[1]])
pi<-length(ort_pat_inc_m[[1]])
bp<-length(ort_ben_phe_m[[1]])
bi<-length(ort_ben_inc_m[[1]])

cp<-length(ort_conf_phe_m[[1]])
ci<-length(ort_conf_inc_m[[1]])
vp<-length(ort_vus_phe_m[[1]])
vi<-length(ort_vus_inc_m[[1]])


ort_mouse_cospat<-ort_mouseallhuman34[grep("^(?=.*Unknown)", ort_mouseallhuman34$Human_Clinical_Significance, perl=TRUE, ignore.case = FALSE),]
ort_mouse_cosunk<-ort_mouseallhuman34[ort_mouseallhuman34$Human_Clinical_Significance == "NEUTRAL" | ort_mouseallhuman34$Human_Clinical_Significance == ""]
ort_mouse_gnodb<-ort_mouseallhuman34[ort_mouseallhuman34$Source == "gnomAD" | ort_mouseallhuman34$Source == "dbSNP"]


ort_cospat_phe_m<-ort_mouse_cospat[ort_mouse_cospat$Mouse_Significance == "phenotypic"]
ort_cospat_inc_m<-ort_mouse_cospat[ort_mouse_cospat$Mouse_Significance == "incidental"]

ort_cosunk_phe_m<-ort_mouse_cosunk[ort_mouse_cosunk$Mouse_Significance == "phenotypic"]
ort_cosunk_inc_m<-ort_mouse_cosunk[ort_mouse_cosunk$Mouse_Significance == "incidental"]

ort_gnodb_phe_m<-ort_mouse_gnodb[ort_mouse_gnodb$Mouse_Significance == "phenotypic"]
ort_gnodb_inc_m<-ort_mouse_gnodb[ort_mouse_gnodb$Mouse_Significance == "incidental"]

ort_cospat_phe_m<-unique(setDT(ort_cospat_phe_m), by = c(1,4,6))
ort_cospat_inc_m<-unique(setDT(ort_cospat_inc_m), by = c(1,4,6))
ort_cosunk_phe_m<-unique(setDT(ort_cosunk_phe_m), by = c(1,4,6))
ort_cosunk_inc_m<-unique(setDT(ort_cosunk_inc_m), by = c(1,4,6))
ort_gnodb_phe_m<-unique(setDT(ort_gnodb_phe_m), by = c(1,4,6))
ort_gnodb_inc_m<-unique(setDT(ort_gnodb_inc_m), by = c(1,4,6))



cospatp<-length(ort_cospat_phe_m[[1]])
cospati<-length(ort_cospat_inc_m[[1]])
cosunkp<-length(ort_cosunk_phe_m[[1]])
cosunki<-length(ort_cosunk_inc_m[[1]])
gnop<-length(ort_gnodb_phe_m[[1]])
gnoi<-length(ort_gnodb_inc_m[[1]])




df<-data.frame(Source = c("ClinVar(Pathogenic)","ClinVar(Benign)","ClinVar(Conflicting)","ClinVar(VUS)","COSMIC(PATHOGENIC)","COSMIC(UNKNOWN)","gnomAD,dbSNP", 
                          "ClinVar(Pathogenic)","ClinVar(Benign)","ClinVar(Conflicting)","ClinVar(VUS)","COSMIC(PATHOGENIC)","COSMIC(UNKNOWN)","gnomAD,dbSNP"), 
               Number = c(59, 138, 25, 333, 1465, 492, 4061, 353, 1047, 140, 2055, 13732, 6014, 48487),
               Type = c(rep("Phenotypic",7), rep("Unknown", 7)), stringsAsFactors = FALSE)

df1<-data.frame(Source = c("ClinVar(Pathogenic)","ClinVar(Benign)","ClinVar(Conflicting)","ClinVar(VUS)","COSMIC(PATHOGENIC)","COSMIC(UNKNOWN)","gnomAD,dbSNP", 
                           "ClinVar(Pathogenic)","ClinVar(Benign)","ClinVar(Conflicting)","ClinVar(VUS)","COSMIC(PATHOGENIC)","COSMIC(UNKNOWN)","gnomAD,dbSNP"), 
                Number = c(59, 138, 25, 333, 1465, 492, 4061, 353, 1047, 140, 2055, 13732, 6014, 48487),
                Number_Percent = c(df$Number[1:7]/(df$Number[1:7]+df$Number[8:14]), df$Number[8:14]/(df$Number[1:7]+df$Number[8:14])),
                Type = c(rep("Phenotypic",7), rep("Unknown", 7)), stringsAsFactors = FALSE)



df<-data.frame(Source = c("ClinVar(Pathogenic)","ClinVar(Benign)","ClinVar(Conflicting)","ClinVar(VUS)","COSMIC(PATHOGENIC)","COSMIC(UNKNOWN)","gnomAD,dbSNP", 
                          "ClinVar(Pathogenic)","ClinVar(Benign)","ClinVar(Conflicting)","ClinVar(VUS)","COSMIC(PATHOGENIC)","COSMIC(UNKNOWN)","gnomAD,dbSNP"), 
               Number = c(pp, bp, cp, vp, cospatp, cosunkp, gnop, pi, bi, ci, vi, cospati, cosunki, gnoi),
               Type = c(rep("Phenotypic",7), rep("Unknown", 7)), stringsAsFactors = FALSE)

df1<-data.frame(Source = c("ClinVar(Pathogenic)","ClinVar(Benign)","ClinVar(Conflicting)","ClinVar(VUS)","COSMIC(PATHOGENIC)","COSMIC(UNKNOWN)","gnomAD,dbSNP", 
                           "ClinVar(Pathogenic)","ClinVar(Benign)","ClinVar(Conflicting)","ClinVar(VUS)","COSMIC(PATHOGENIC)","COSMIC(UNKNOWN)","gnomAD,dbSNP"), 
                Number = c(pp, bp, cp, vp, cospatp, cosunkp, gnop, pi, bi, ci, vi, cospati, cosunki, gnoi),
                Number_Percent = c(df$Number[1:7]/(df$Number[1:7]+df$Number[8:14]), df$Number[8:14]/(df$Number[1:7]+df$Number[8:14])),
                Type = c(rep("Phenotypic",7), rep("Unknown", 7)), stringsAsFactors = FALSE)



1130/637

# Percentages

a<-ggplot(data = df1[c(1,8),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")



b<-ggplot(data = df1[c(2,9),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


c<-ggplot(data = df1[c(3,10),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


d<-ggplot(data = df1[c(4,11),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


e<-ggplot(data = df1[c(5,12),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


f<-ggplot(data = df1[c(6,13),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


g<-ggplot(data = df1[c(7,14),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.4, size = 10)) +
  labs(title = "Mouse-dbSNP & gnomAD Orthologous Variants") +
  scale_fill_brewer(palette="Blues")


# Numbers

a<-ggplot(data = df1[c(1,8),], aes(x = Source, y = Number, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")



b<-ggplot(data = df1[c(2,9),], aes(x = Source, y = Number, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


c<-ggplot(data = df1[c(3,10),], aes(x = Source, y = Number, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


d<-ggplot(data = df1[c(4,11),], aes(x = Source, y = Number, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


e<-ggplot(data = df1[c(5,12),], aes(x = Source, y = Number, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


f<-ggplot(data = df1[c(6,13),], aes(x = Source, y = Number, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


g<-ggplot(data = df1[c(7,14),], aes(x = Source, y = Number, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.4, size = 10)) +
  labs(title = "Mouse-dbSNP & gnomAD Orthologous Variants") +
  scale_fill_brewer(palette="Blues")


# Combine plots
combined1 <-a + b + c + d & theme(legend.position = "right", plot.title = element_text(hjust = 4, size = 10))
combined1 + plot_layout(ncol = 4, nrow = 1, guides = "collect") +
  labs(title = "Mouse-ClinVar Orthologous Variants")


combined2 <-e + f  & theme(legend.position = "right", plot.title = element_text(hjust = 1.3, size = 10))
combined2 + plot_layout(ncol = 2, nrow = 1, guides = "collect") +
  labs(title = "Mouse-COSMIC Orthologous Variants")




# C. elegans #
ort_humancelegansCombined$Human_Clinical_Significance[ort_humancelegansCombined$Human_Clinical_Significance == "PATHOGENIC"]<-"Unknown"

ort_celegans_pat<-ort_humancelegansCombined[grep("^(?=.*pathogenic)(?!.*Conflicting)", ort_humancelegansCombined$Human_Clinical_Significance, perl=TRUE, ignore.case = TRUE),]
ort_celegans_ben<-ort_humancelegansCombined[grep("^(?=.*benign)(?!.*Conflicting)", ort_humancelegansCombined$Human_Clinical_Significance, perl=TRUE, ignore.case = TRUE),]
ort_celegans_conf<-ort_humancelegansCombined[grep("^(?=.*conflicting)", ort_humancelegansCombined$Human_Clinical_Significance, perl=TRUE, ignore.case = TRUE),]
ort_celegans_vus<-ort_humancelegansCombined[grep("^(?=.*uncertain)|(?=.*not provided)|(?=.*no interpretation)", ort_humancelegansCombined$Human_Clinical_Significance, perl=TRUE, ignore.case = TRUE),]

ort_pat_phe_c<-ort_celegans_pat[ort_celegans_pat$C_elegans_Significance == "Phenotypic"]
ort_ben_phe_c<-ort_celegans_ben[ort_celegans_ben$C_elegans_Significance == "Phenotypic"]
ort_conf_phe_c<-ort_celegans_conf[ort_celegans_conf$C_elegans_Significance == "Phenotypic"]
ort_vus_phe_c<-ort_celegans_vus[ort_celegans_vus$C_elegans_Significance == "Phenotypic"]

ort_pat_inc_c<-ort_celegans_pat[ort_celegans_pat$C_elegans_Significance == "Unknown"]
ort_ben_inc_c<-ort_celegans_ben[ort_celegans_ben$C_elegans_Significance == "Unknown"]
ort_conf_inc_c<-ort_celegans_conf[ort_celegans_conf$C_elegans_Significance == "Unknown"]
ort_vus_inc_c<-ort_celegans_vus[ort_celegans_vus$C_elegans_Significance == "Unknown"]

ort_pat_phe_c<-unique(setDT(ort_pat_phe_c), by = c(1,4,6))
ort_ben_phe_c<-unique(setDT(ort_ben_phe_c), by = c(1,4,6))
ort_conf_phe_c<-unique(setDT(ort_conf_phe_c), by = c(1,4,6))
ort_vus_phe_c<-unique(setDT(ort_vus_phe_c), by = c(1,4,6))

ort_pat_inc_c<-unique(setDT(ort_pat_inc_c), by = c(1,4,6))
ort_ben_inc_c<-unique(setDT(ort_ben_inc_c), by = c(1,4,6))
ort_conf_inc_c<-unique(setDT(ort_conf_inc_c), by = c(1,4,6))
ort_vus_inc_c<-unique(setDT(ort_vus_inc_c), by = c(1,4,6))

c_pp<-length(ort_pat_phe_c[[1]])
c_pi<-length(ort_pat_inc_c[[1]])
c_bp<-length(ort_ben_phe_c[[1]])
c_bi<-length(ort_ben_inc_c[[1]])

c_cp<-length(ort_conf_phe_c[[1]])
c_ci<-length(ort_conf_inc_c[[1]])
c_vp<-length(ort_vus_phe_c[[1]])
c_vi<-length(ort_vus_inc_c[[1]])


ort_celegans_cospat<-ort_humancelegansCombined[grep("^(?=.*Unknown)", ort_humancelegansCombined$Human_Clinical_Significance, perl=TRUE, ignore.case = FALSE),]
ort_celegans_cosunk<-ort_humancelegansCombined[ort_humancelegansCombined$Human_Clinical_Significance == "NEUTRAL" | ort_humancelegansCombined$Human_Clinical_Significance == ""]
ort_celegans_gnodb<-ort_humancelegansCombined[ort_humancelegansCombined$Human_Source == "gnomAD" | ort_humancelegansCombined$Human_Source == "dbSNP"]


ort_cospat_phe_c<-ort_celegans_cospat[ort_celegans_cospat$C_elegans_Significance == "Phenotypic"]
ort_cospat_inc_c<-ort_celegans_cospat[ort_celegans_cospat$C_elegans_Significance == "Unknown"]

ort_cosunk_phe_c<-ort_celegans_cosunk[ort_celegans_cosunk$C_elegans_Significance == "Phenotypic"]
ort_cosunk_inc_c<-ort_celegans_cosunk[ort_celegans_cosunk$C_elegans_Significance == "Unknown"]

ort_gnodb_phe_c<-ort_celegans_gnodb[ort_celegans_gnodb$C_elegans_Significance == "Phenotypic"]
ort_gnodb_inc_c<-ort_celegans_gnodb[ort_celegans_gnodb$C_elegans_Significance == "Unknown"]

ort_cospat_phe_c<-unique(setDT(ort_cospat_phe_c), by = c(1,4,6))
ort_cospat_inc_c<-unique(setDT(ort_cospat_inc_c), by = c(1,4,6))
ort_cosunk_phe_c<-unique(setDT(ort_cosunk_phe_c), by = c(1,4,6))
ort_cosunk_inc_c<-unique(setDT(ort_cosunk_inc_c), by = c(1,4,6))
ort_gnodb_phe_c<-unique(setDT(ort_gnodb_phe_c), by = c(1,4,6))
ort_gnodb_inc_c<-unique(setDT(ort_gnodb_inc_c), by = c(1,4,6))



c_cospatp<-length(ort_cospat_phe_c[[1]])
c_cospati<-length(ort_cospat_inc_c[[1]])
c_cosunkp<-length(ort_cosunk_phe_c[[1]])
c_cosunki<-length(ort_cosunk_inc_c[[1]])
c_gnop<-length(ort_gnodb_phe_c[[1]])
c_gnoi<-length(ort_gnodb_inc_c[[1]])


df<-data.frame(Source = c("ClinVar(Pathogenic)","ClinVar(Benign)","ClinVar(Conflicting)","ClinVar(VUS)","COSMIC(PATHOGENIC)","COSMIC(UNKNOWN)","gnomAD,dbSNP", 
                          "ClinVar(Pathogenic)","ClinVar(Benign)","ClinVar(Conflicting)","ClinVar(VUS)","COSMIC(PATHOGENIC)","COSMIC(UNKNOWN)","gnomAD,dbSNP"), 
               Number = c(71, 11, 6, 80, 479, 55, 509, 378, 584, 111, 1290, 11372, 2788, 26894),
               Type = c(rep("Phenotypic",7), rep("Unknown", 7)), stringsAsFactors = FALSE)

df1<-data.frame(Source = c("ClinVar(Pathogenic)","ClinVar(Benign)","ClinVar(Conflicting)","ClinVar(VUS)","COSMIC(PATHOGENIC)","COSMIC(UNKNOWN)","gnomAD,dbSNP", 
                           "ClinVar(Pathogenic)","ClinVar(Benign)","ClinVar(Conflicting)","ClinVar(VUS)","COSMIC(PATHOGENIC)","COSMIC(UNKNOWN)","gnomAD,dbSNP"), 
                Number = c(71, 11, 6, 80, 479, 55, 509, 378, 584, 111, 1290, 11372, 2788, 26894),
                Number_Percent = c(df$Number[1:7]/(df$Number[1:7]+df$Number[8:14]), df$Number[8:14]/(df$Number[1:7]+df$Number[8:14])),
                Type = c(rep("Phenotypic",7), rep("Unknown", 7)), stringsAsFactors = FALSE)



df<-data.frame(Source = c("ClinVar(Pathogenic)","ClinVar(Benign)","ClinVar(Conflicting)","ClinVar(VUS)","COSMIC(PATHOGENIC)","COSMIC(UNKNOWN)","gnomAD,dbSNP", 
                          "ClinVar(Pathogenic)","ClinVar(Benign)","ClinVar(Conflicting)","ClinVar(VUS)","COSMIC(PATHOGENIC)","COSMIC(UNKNOWN)","gnomAD,dbSNP"), 
               Number = c(c_pp, c_bp, c_cp, c_vp, c_cospatp, c_cosunkp, c_gnop, c_pi, c_bi, c_ci, c_vi, c_cospati, c_cosunki, c_gnoi),
               Type = c(rep("Phenotypic",7), rep("Unknown", 7)), stringsAsFactors = FALSE)

df1<-data.frame(Source = c("ClinVar(Pathogenic)","ClinVar(Benign)","ClinVar(Conflicting)","ClinVar(VUS)","COSMIC(PATHOGENIC)","COSMIC(UNKNOWN)","gnomAD,dbSNP", 
                           "ClinVar(Pathogenic)","ClinVar(Benign)","ClinVar(Conflicting)","ClinVar(VUS)","COSMIC(PATHOGENIC)","COSMIC(UNKNOWN)","gnomAD,dbSNP"), 
                Number = c(c_pp, c_bp, c_cp, c_vp, c_cospatp, c_cosunkp, c_gnop, c_pi, c_bi, c_ci, c_vi, c_cospati, c_cosunki, c_gnoi),
                Number_Percent = c(df$Number[1:7]/(df$Number[1:7]+df$Number[8:14]), df$Number[8:14]/(df$Number[1:7]+df$Number[8:14])),
                Type = c(rep("Phenotypic",7), rep("Unknown", 7)), stringsAsFactors = FALSE)




1130/637

# Percent
a<-ggplot(data = df1[c(1,8),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_brewer(palette="Blues")
  


b<-ggplot(data = df1[c(2,9),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


c<-ggplot(data = df1[c(3,10),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


d<-ggplot(data = df1[c(4,11),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


e<-ggplot(data = df1[c(5,12),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


f<-ggplot(data = df1[c(6,13),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


g<-ggplot(data = df1[c(7,14),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.35, size = 10)) +
  labs(title = "C. elegans-dbSNP & gnomAD Orthologous Variants")+
  scale_fill_brewer(palette="Blues")



# Numbers

a<-ggplot(data = df1[c(1,8),], aes(x = Source, y = Number, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_brewer(palette="Blues")



b<-ggplot(data = df1[c(2,9),], aes(x = Source, y = Number, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


c<-ggplot(data = df1[c(3,10),], aes(x = Source, y = Number, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


d<-ggplot(data = df1[c(4,11),], aes(x = Source, y = Number, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


e<-ggplot(data = df1[c(5,12),], aes(x = Source, y = Number, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


f<-ggplot(data = df1[c(6,13),], aes(x = Source, y = Number, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


g<-ggplot(data = df1[c(7,14),], aes(x = Source, y = Number, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.35, size = 10)) +
  labs(title = "C. elegans-dbSNP & gnomAD Orthologous Variants")+
  scale_fill_brewer(palette="Blues")

# Combine plots
combined2 <- a + b + c + d & theme(legend.position = "right", plot.title = element_text(hjust = 3.5, size = 10))
combined2 + plot_layout(ncol = 4, nrow = 1, guides = "collect") +
  labs(title = "C. elegans-ClinVar Orthologous Variants")

combined2 <- e + f & theme(legend.position = "right", plot.title = element_text(hjust = 1.5, size = 10))
combined2 + plot_layout(ncol = 2, nrow = 1, guides = "collect") +
  labs(title = "C. elegans-COSMIC Orthologous Variants")





# Some non-sense


# Mouse

df<-data.frame(Source = c("ClinVar(Orthologs)","ClinVar(Orthologs)","Mouse(Pathogenic)","Mouse(Pathogenic)"), 
               Number = c(406, 1230, 86, 320),
               Type = c("Pathogenic","Benign","Mouse(Phenotypic)","Mouse(Unknown)"))
df$Type<-factor(df$Type, levels = c("ClinVar(Pathogenic)","ClinVar(Benign)","Mouse(Phenotypic)","Mouse(Unknown)"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")

df<-data.frame(Source = c("ClinVar(Orthologs)","ClinVar(Orthologs)","Mouse(Pathogenic)","Mouse(Pathogenic)"), 
               Number = c(406, 1230, 86, 320),
               Number_Percent = c(df$Number[1]/(df$Number[1]+df$Number[2]), df$Number[2]/(df$Number[1]+df$Number[2]), 
                                  df$Number[3]/(df$Number[3]+df$Number[4]), df$Number[4]/(df$Number[3]+df$Number[4])),
               Type = c("Pathogenic","Benign","Phenotypic","Unknown"))
df$Type<-factor(df$Type, levels = c("Pathogenic","Benign","Phenotypic","Unknown"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")

df1<-data.frame(Source = c("ClinVar(Orthologs)","ClinVar(Orthologs)","Mouse(Benign)","Mouse(Benign)"), 
               Number = c(406, 1230, 203, 1027),
               Type = c("Pathogenic","Benign","Phenotypic","Unknown"))
df1$Type<-factor(df1$Type, levels = c("Pathogenic","Benign","Phenotypic","Unknown"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")

df1<-data.frame(Source = c("ClinVar(Orthologs)","ClinVar(Orthologs)","Mouse(Benign)","Mouse(Benign)"), 
               Number = c(406, 1230, 203, 1027),
               Number_Percent = c(df1$Number[1]/(df1$Number[1]+df1$Number[2]), df1$Number[2]/(df1$Number[1]+df1$Number[2]), 
                                  df1$Number[3]/(df1$Number[3]+df1$Number[4]), df1$Number[4]/(df1$Number[3]+df1$Number[4])),
               Type = c("Pathogenic","Benign","Phenotypic","Unknown"))
df1$Type<-factor(df1$Type, levels = c("Pathogenic","Benign","Phenotypic","Unknown"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")



a<-ggplot(data = df[c(1,2),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")



b<-ggplot(data = df[c(3,4),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Greens")

c<-ggplot(data = df1[c(3,4),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Greens")



combined2 <- a + b + c & theme(legend.position = "right", plot.title = element_text(size = 12, hjust = 5))
combined2 + plot_layout(ncol = 3, nrow = 1, guides = "collect") +
  labs(title = "ClinVar - Mouse orthologs")









df<-data.frame(Source = c("ClinVar(Orthologs)","ClinVar(Orthologs)","C. elegans(Pathogenic)","C. elegans(Pathogenic)"), 
               Number = c(850, 1265, 148, 702),
               Number_Percent = c(df$Number[1]/(df$Number[1]+df$Number[2]), df$Number[2]/(df$Number[1]+df$Number[2]), 
                                  df$Number[3]/(df$Number[3]+df$Number[4]), df$Number[4]/(df$Number[3]+df$Number[4])),
               Type = c("Pathogenic","Benign","Phenotypic","Unknown"))
df$Type<-factor(df$Type, levels = c("Pathogenic","Benign","Phenotypic","Unknown"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")



a<-ggplot(data = df[c(1,2),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")



b<-ggplot(data = df[c(3,4),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


combined2 <- a + b & theme(legend.position = "right", plot.title = element_text(hjust = 3.5, size = 10))
combined2 + plot_layout(ncol = 2, nrow = 1, guides = "collect") +
  labs(title = "Mouse - ClinVar")





df<-data.frame(Source = c("ClinVar(Orthologs)","ClinVar(Orthologs)","C. elegans(Pathogenic)","C. elegans(Pathogenic)"), 
               Number = c(850, 1265, 148, 702),
               Type = c("Pathogenic","Benign","Phenotypic","Unknown"))
df$Type<-factor(df$Type, levels = c("Pathogenic","Benign","Phenotypic","Unknown"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")

df<-data.frame(Source = c("ClinVar(Orthologs)","ClinVar(Orthologs)","C. elegans(Pathogenic)","C. elegans(Pathogenic)"), 
               Number = c(850, 1265, 148, 702),
               Number_Percent = c(df$Number[1]/(df$Number[1]+df$Number[2]), df$Number[2]/(df$Number[1]+df$Number[2]), 
                                  df$Number[3]/(df$Number[3]+df$Number[4]), df$Number[4]/(df$Number[3]+df$Number[4])),
               Type = c("Pathogenic","Benign","Phenotypic","Unknown"))
df$Type<-factor(df$Type, levels = c("Pathogenic","Benign","Phenotypic","Unknown"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")

df1<-data.frame(Source = c("ClinVar(Orthologs)","ClinVar(Orthologs)","C. elegans(Benign)","C. elegans(Benign)"), 
                Number = c(850, 1265, 22, 1243),
                Type = c("Pathogenic","Benign","Phenotypic","Unknown"))
df1$Type<-factor(df1$Type, levels = c("Pathogenic","Benign","Phenotypic","Unknown"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")

df1<-data.frame(Source = c("ClinVar(Orthologs)","ClinVar(Orthologs)","C. elegans(Benign)","C. elegans(Benign)"), 
               Number = c(850, 1265, 22, 1243),
               Number_Percent = c(df1$Number[1]/(df1$Number[1]+df1$Number[2]), df1$Number[2]/(df1$Number[1]+df1$Number[2]), 
                                  df1$Number[3]/(df1$Number[3]+df1$Number[4]), df1$Number[4]/(df1$Number[3]+df1$Number[4])),
               Type = c("Pathogenic","Benign","Phenotypic","Unknown"))
df1$Type<-factor(df1$Type, levels = c("Pathogenic","Benign","Phenotypic","Unknown"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")



a<-ggplot(data = df1[c(1,2),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")



b<-ggplot(data = df[c(3,4),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Greens")


c<-ggplot(data = df1[c(3,4),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Greens")


# combined2 <- a + b & theme(legend.position = "right", plot.title = element_text(hjust = 3.5, size = 10))
# combined2 + plot_layout(ncol = 2, nrow = 1, guides = "collect") +
#   labs(title = "Mouse - ClinVar")

combined2 <- a + b + c & theme(legend.position = "right", plot.title = element_text(hjust = 5, size = 12))
combined2 + plot_layout(ncol = 3, nrow = 1, guides = "collect") +
  labs(title = "ClinVar-Mouse orthologs")




# All clinvar


df<-data.frame(Source = c("ClinVar(All)","ClinVar(All)","Mouse(All)","Mouse(All)"), 
               Number = c(17548, 120346, 29521, 179261),
               Type = c("ClinVar(All)","ClinVar(All)","Mouse(All)","Mouse(All)"))
df$Type<-factor(df$Type, levels = c("ClinVar(Pathogenic)","ClinVar(Benign)","Phenotypic","Unknown"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")

df<-data.frame(Source = c("ClinVar(All/Mouse homologs)","ClinVar(All/Mouse homologs)","Mouse(All)","Mouse(All)"), 
               Number = c(17548, 120346, 29521, 179261),
               Number_Percent = c(df$Number[1]/(df$Number[1]+df$Number[2]), df$Number[2]/(df$Number[1]+df$Number[2]), 
                                  df$Number[3]/(df$Number[3]+df$Number[4]), df$Number[4]/(df$Number[3]+df$Number[4])),
               Type = c("ClinVar(Pathogenic)","ClinVar(Benign)","Phenotypic","Unknown"))
df$Type<-factor(df$Type, levels = c("ClinVar(Pathogenic)","ClinVar(Benign)","Phenotypic","Unknown"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")

df1<-data.frame(Source = c("ClinVar(All/C. elegans homologs)","ClinVar(All/C. elegans homologs)","C. elegans(All)","C. elegans(All)"), 
                Number = c(10732, 73199, 4421, 245358),
                Type = c("ClinVar(All)","ClinVar(All)","C. elegans(All)","C. elegans(All)"))
df1$Type<-factor(df1$Type, levels = c("ClinVar(Pathogenic)","ClinVar(Benign)","Phenotypic","Unknown"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")

df1<-data.frame(Source = c("ClinVar(All/C. elegans homologs)","ClinVar(All/C. elegans homologs)","C. elegans(All)","C. elegans(All)"), 
                Number = c(10732, 73199, 4421, 245358),
                Number_Percent = c(df1$Number[1]/(df1$Number[1]+df1$Number[2]), df1$Number[2]/(df1$Number[1]+df1$Number[2]), 
                                   df1$Number[3]/(df1$Number[3]+df1$Number[4]), df1$Number[4]/(df1$Number[3]+df1$Number[4])),
                Type = c("ClinVar(Pathogenic)","ClinVar(Benign)","Phenotypic","Unknown"))
df1$Type<-factor(df1$Type, levels = c("ClinVar(Pathogenic)","ClinVar(Benign)","Phenotypic","Unknown"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")



a<-ggplot(data = df[c(1,2),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")



b<-ggplot(data = df[c(3,4),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Greens")

c<-ggplot(data = df1[c(3,4),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Greens")

d<-ggplot(data = df1[c(1,2),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


combined2 <- a + b + d + c & theme(legend.position = "right", plot.title = element_text(hjust = 3, size = 12, vjust = 80))
combined2 + plot_layout(ncol = 2, nrow = 2, guides = "collect") +
  labs(title = "ClinVar-Mouse-C.elegans all variations")







# Some sense

ort_phe_m<-ort_mouseallhuman34[ort_mouseallhuman34$Mouse_Significance == "phenotypic"]
ort_inc_m<-ort_mouseallhuman34[ort_mouseallhuman34$Mouse_Significance == "incidental"]

ort_phe_pat_m<-ort_phe_m[grep("^(?=.*pathogenic)(?!.*Conflicting)", ort_phe_m$Human_Clinical_Significance, perl=TRUE, ignore.case = TRUE),]
ort_phe_ben_m<-ort_phe_m[grep("^(?=.*benign)(?!.*Conflicting)", ort_phe_m$Human_Clinical_Significance, perl=TRUE, ignore.case = TRUE),]
ort_inc_pat_m<-ort_inc_m[grep("^(?=.*pathogenic)(?!.*Conflicting)", ort_inc_m$Human_Clinical_Significance, perl=TRUE, ignore.case = TRUE),]
ort_inc_ben_m<-ort_inc_m[grep("^(?=.*benign)(?!.*Conflicting)", ort_inc_m$Human_Clinical_Significance, perl=TRUE, ignore.case = TRUE),]

ort_phe_pat_m<-unique(setDT(ort_phe_pat_m), by = c(2,5,7))
ort_phe_ben_m<-unique(setDT(ort_phe_ben_m), by = c(2,5,7))
ort_inc_pat_m<-unique(setDT(ort_inc_pat_m), by = c(2,5,7))
ort_inc_ben_m<-unique(setDT(ort_inc_ben_m), by = c(2,5,7))


df<-data.frame(Source = c("ClinVar(All)","ClinVar(All)","Mouse(All)","Mouse(All)"), 
               Number = c(86, 203, 323, 1027),
               Type = c("ClinVar(All)","ClinVar(All)","Mouse(All)","Mouse(All)"))
df$Type<-factor(df$Type, levels = c("ClinVar(Pathogenic)","ClinVar(Benign)","Phenotypic","Unknown"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")

df<-data.frame(Source = c("Mouse(Phenotypic)","Mouse(Phenotypic)","Mouse(Unknown)","Mouse(Unknown)"), 
               Number = c(86, 203, 323, 1027),
               Number_Percent = c(df$Number[1]/(df$Number[1]+df$Number[2]), df$Number[2]/(df$Number[1]+df$Number[2]), 
                                  df$Number[3]/(df$Number[3]+df$Number[4]), df$Number[4]/(df$Number[3]+df$Number[4])),
               Type = c("Pathogenic","Benign","Pathogenic","Benign"))
df$Type<-factor(df$Type, levels = c("Pathogenic","Benign"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")

chisq_df<-data.frame(Pathogenic = c(86,323), Benign = c(203,1027))
row.names(chisq_df)<-c("Phenotypic","Unknown")
test<-chisq.test(chisq_df)



a<-ggplot(data = df[c(1,2),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")



b<-ggplot(data = df[c(3,4),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


combined2 <- a + b & theme(legend.position = "right", plot.title = element_text(hjust = 4.5, size = 10))
combined2 + plot_layout(ncol = 2, nrow = 1, guides = "collect") +
  labs(title = "Mouse Orthologous Variants")


a<-ggplot(data = df[c(1,2),], aes(x = Source, y = Number, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Number of Variants") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")



b<-ggplot(data = df[c(3,4),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Number of Variants") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


combined2 <- a + b & theme(legend.position = "right", plot.title = element_text(hjust = 4.5, size = 10))
combined2 + plot_layout(ncol = 2, nrow = 1, guides = "collect") +
  labs(title = "Mouse Orthologous Variants")




# C. elegans

ort_phe_m<-ort_humancelegansCombined[ort_humancelegansCombined$C_elegans_Significance == "Phenotypic"]
ort_inc_m<-ort_humancelegansCombined[ort_humancelegansCombined$C_elegans_Significance == "Unknown" | ort_humancelegansCombined$C_elegans_Significance == "Non-phenotypic"]

ort_phe_pat_m<-ort_phe_m[grep("^(?=.*pathogenic)(?!.*Conflicting)", ort_phe_m$Human_Clinical_Significance, perl=TRUE, ignore.case = TRUE),]
ort_phe_ben_m<-ort_phe_m[grep("^(?=.*benign)(?!.*Conflicting)", ort_phe_m$Human_Clinical_Significance, perl=TRUE, ignore.case = TRUE),]
ort_inc_pat_m<-ort_inc_m[grep("^(?=.*pathogenic)(?!.*Conflicting)", ort_inc_m$Human_Clinical_Significance, perl=TRUE, ignore.case = TRUE),]
ort_inc_ben_m<-ort_inc_m[grep("^(?=.*benign)(?!.*Conflicting)", ort_inc_m$Human_Clinical_Significance, perl=TRUE, ignore.case = TRUE),]

ort_phe_pat_m<-unique(setDT(ort_phe_pat_m), by = c(2,5,7))
ort_phe_ben_m<-unique(setDT(ort_phe_ben_m), by = c(2,5,7))
ort_inc_pat_m<-unique(setDT(ort_inc_pat_m), by = c(2,5,7))
ort_inc_ben_m<-unique(setDT(ort_inc_ben_m), by = c(2,5,7))

pp<-length(ort_phe_pat_m[[1]])
pb<-length(ort_phe_ben_m[[1]])
ip<-length(ort_inc_pat_m[[1]])
ib<-length(ort_inc_ben_m[[1]])

df<-data.frame(Source = c("ClinVar(All)","ClinVar(All)","Mouse(All)","Mouse(All)"), 
               Number = c(154, 22, 702, 1243),
               Type = c("ClinVar(All)","ClinVar(All)","Mouse(All)","Mouse(All)"))
df$Type<-factor(df$Type, levels = c("ClinVar(Pathogenic)","ClinVar(Benign)","Phenotypic","Unknown"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")

df<-data.frame(Source = c("C. elegans(Phenotypic)","C. elegans(Phenotypic)","C. elegans(Unknown)","C. elegans(Unknown)"), 
               Number = c(154, 22, 702, 1243),
               Number_Percent = c(df$Number[1]/(df$Number[1]+df$Number[2]), df$Number[2]/(df$Number[1]+df$Number[2]), 
                                  df$Number[3]/(df$Number[3]+df$Number[4]), df$Number[4]/(df$Number[3]+df$Number[4])),
               Type = c("Pathogenic","Benign","Pathogenic","Benign"))
df$Type<-factor(df$Type, levels = c("Pathogenic","Benign"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")


df<-data.frame(Source = c("ClinVar(All)","ClinVar(All)","Mouse(All)","Mouse(All)"), 
               Number = c(pp, pb, ip, ib),
               Type = c("ClinVar(All)","ClinVar(All)","Mouse(All)","Mouse(All)"))
df$Type<-factor(df$Type, levels = c("ClinVar(Pathogenic)","ClinVar(Benign)","Phenotypic","Unknown"))
#df$Type<-relevel(df$Type, "ClinVar(Pathogenic)")

df<-data.frame(Source = c("C. elegans(Phenotypic)","C. elegans(Phenotypic)","C. elegans(Unknown)","C. elegans(Unknown)"), 
               Number = c(pp, pb, ip, ib),
               Number_Percent = c(df$Number[1]/(df$Number[1]+df$Number[2]), df$Number[2]/(df$Number[1]+df$Number[2]), 
                                  df$Number[3]/(df$Number[3]+df$Number[4]), df$Number[4]/(df$Number[3]+df$Number[4])),
               Type = c("Pathogenic","Benign","Pathogenic","Benign"))
df$Type<-factor(df$Type, levels = c("Pathogenic","Benign"))


chisq_df<-data.frame(Pathogenic = c(pp,ip), Benign = c(pb,ib))
row.names(chisq_df)<-c("Phenotypic","Unknown")
chisq.test(chisq_df)

a<-ggplot(data = df[c(1,2),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")



b<-ggplot(data = df[c(3,4),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = paste0(round(100*Number_Percent, digits = 1),"%")), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Variants (%)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


combined2 <- a + b & theme(legend.position = "right", plot.title = element_text(hjust = 2.5, size = 10))
combined2 + plot_layout(ncol = 2, nrow = 1, guides = "collect") +
  labs(title = "C. elegans Orthologous Variants")


a<-ggplot(data = df[c(1,2),], aes(x = Source, y = Number, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Number of Variants") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")



b<-ggplot(data = df[c(3,4),], aes(x = Source, y = Number_Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_label(aes(label = Number), size = 5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous("Number of Variants") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), text = element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues")


combined2 <- a + b & theme(legend.position = "right", plot.title = element_text(hjust = 2.5, size = 10))
combined2 + plot_layout(ncol = 2, nrow = 1, guides = "collect") +
  labs(title = "C. elegans Orthologous Variants")

