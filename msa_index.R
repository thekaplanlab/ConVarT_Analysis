
library(pbapply)
library(data.table)
library(pbapply)

# msa file #

table<-fread("msa_table.txt") %>%
  filter(grepl('musculus', fasta)) %>%
  filter(grepl('sapiens', fasta)) %>%
  filter(grepl('elegans', fasta))

table<-table[grep("sapiens", table$fasta, ignore.case = TRUE),]
table<-table[grep("elegans", table$fasta, ignore.case = TRUE),]

table2<-strsplit(as.character(table$fasta), "\n>")

table2<-pblapply(X = table2, FUN = function(t) gsub(pattern = "\n", replacement = "", x = t))

table233<-pblapply(table2, function(x) length(x) == 3)
table23<-table2[unlist(table233)]

h_ag<-c()
human_i<-c()
m_ag<-c()
mouse_i<-c()
c_ag<-c()
celegans_i<-c()

for (i in 1:length(table23)){
  
  human_a<-c(unlist(table23[[i]][1]))
  human_a<-c(unlist(strsplit(human_a, "]")))
  human_ag<-c(unlist(strsplit(human_a, " \\["))[[1]][1])
  human_ag<-c(unlist(strsplit(human_ag, "\\."))[[1]][1])
  human_ag<-gsub(">","",human_ag)
  h_ag<-c(h_ag, human_ag)
  human_i<-c(human_i, i)
  
  mouse_a<-c(unlist(table23[[i]][2]))
  mouse_a<-c(unlist(strsplit(mouse_a, "]")))
  mouse_ag<-c(unlist(strsplit(mouse_a, " \\["))[[1]][1])
  mouse_ag<-c(unlist(strsplit(mouse_ag, "\\."))[[1]][1])
  mouse_ag<-gsub(">","",mouse_ag)
  m_ag<-c(m_ag, mouse_ag)
  mouse_i<-c(mouse_i, i)
  
  celegans_a<-c(unlist(table23[[i]][3]))
  celegans_a<-c(unlist(strsplit(celegans_a, "]")))
  celegans_ag<-c(unlist(strsplit(celegans_a, " \\["))[[1]][1])
  celegans_ag<-c(unlist(strsplit(celegans_ag, "\\."))[[1]][1])
  celegans_ag<-gsub(">","",celegans_ag)
  c_ag<-c(c_ag, celegans_ag)
  celegans_i<-c(celegans_i, i)
}

indexed<-data.table(human_refseq = h_ag, mouse_refseq = m_ag, celegans_refseq = c_ag, index = human_i)


