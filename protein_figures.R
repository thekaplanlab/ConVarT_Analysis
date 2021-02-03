
library(trackViewer)

SNP <- c(69,131,145,167,185,200,257,359,386,404)
gene1 <- GRanges("chr1", IRanges(SNP, width=1, names=c("E69G", "Q131L", "S145F", "F167Y", "A185P", "F200Y", "M257I", "R359H", "T386I", "D404N")))
gene1$border <- sample(c("gray30"), length(SNP), replace=TRUE)
x <- c(0,0,0,50,0,0,50,50,0,0)
y <- c(100,100,100,50,100,100,50,0,0,100)
z <- c(0,0,0,0,0,0,0,50,100,0)
gene1$value1 <- x
gene1$value2 <- y
gene1$value3 <- z
gene1$color <- rep(list(c("#AF58BA","#FFC61E","#F28522")), length(SNP))
gene1.rot <- gene1
gene1.rot$label.parameter.rot <- 45
gene1.rot$SNPsideID<-"top"
features1 <- GRanges("chr1", IRanges(c(1, 3, 261, 444), 
                                     width=c(0, 209, 122, 0),
                                     names=c("", "Tubulin","Tubulin_C", "")))
features1$fill <- c("#FF1F5B", "#FF1F5B","#009ADE", "#009ADE")
features1$height <- c(0, 0.2, 0.2, 0)

# lolliplot(gene1.rot, 
#           features1,
#           yaxis=FALSE, ylab = "ben-1\nNP_497728", xaxis = xaxis, legend = legend, type = "pie")
# 


SNP2<- c(41,156,167,185,251,257,277,359,386,391,63,189,345,395)
gene2 <- GRanges("chr1", IRanges(SNP2, width=1, names=c("D41NA","R156Q", "F167L", "A185T", "R251W", "M257I", "G277A", "R359Q", "T386A",
                                                        "R391H","A63A","V189V","I345I","L395L")))
gene2$border <- sample(c("gray30"), length(SNP2), replace=TRUE)
x1 <- c(0,0,50,0,0,50,0,50,0,0,0,0,0,0)
y1 <- c(0,0,0,0,0,0,0,0,0,100,0,0,0,0)
z1 <- c(100,100,50,100,100,50,100,50,100,0,0,0,0,0)
t1 <- c(0,0,0,0,0,0,0,0,0,0,100,100,100,100)

gene2$value1 <- x1
gene2$value2 <- y1
gene2$value3 <- z1
gene2$value4 <- t1

## the length of the color should be no less than that of value1 or value2
gene2$color <- rep(list(c("#AF58BA","#FFC61E","#F28522","#7CCBA2")), length(SNP2))
gene2.rot <- gene2
gene2.rot$label.parameter.rot <- 45
gene2.rot$SNPsideID<-"top"
features2 <- GRanges("chr1", IRanges(c(1, 3, 261, 445), 
                                    width=c(0, 209, 122, 0),
                                    names=c("", "Tubulin","Tubulin_C", "")))
features2$fill <- c("#FF1F5B", "#FF1F5B","#009ADE","#009ADE")
features2$height <- c(0, 0.2, 0.2, 0)

lolliplot(gene2.rot, 
          features2,
          yaxis=FALSE, ylab = "TUBB4B\nNP_006079", xaxis = xaxis, type = "pie")



legend <- c("#458CB0","#FF3333","#330000","#E0E0E0") ## legend fill color
names(legend) <- c("Phenotypic","Orthologous","Pathogenic","Unknown")
xaxis <- c(1, 100, 200, 300, 400)




SNP3<- c(41,77)
gene3 <- GRanges("chr1", IRanges(SNP3, width=1, names=c("D41G","R77H")))
gene3$color <- sample.int(6, length(SNP3), replace=TRUE)
gene3$border <- sample(c("gray30"), length(SNP3), replace=TRUE)
#gene3$score <- sample.int(5, length(gene3), replace = TRUE)
g<-c(50,0)
gene3$value1 <- g
gene3$value2 <- 100-g

gene3$color <- rep(list(c("#FF3333", "#E0E0E0")), length(SNP3))
gene3.rot <- gene3
gene3.rot$label.parameter.rot <- 45
gene3.rot$SNPsideID<-"top"
features3 <- GRanges("chr1", IRanges(c(1, 3, 261, 445), 
                                     width=c(0, 209, 122, 0),
                                     names=c("", "Tubulin","Tubulin_C", "")))
features3$fill <- c("#FF8833", "#51C6E6","#FF8833", "#51C6E6")
features3$height <- c(0, 0.2, 0.2, 0)
legend <- c("#AF58BA","#FFC61E","#F28522","#7CCBA2") ## legend fill color
names(legend) <- c("Orthologous","Phenotypic/Pathogenic","Unknown","Benign/Likely benign")


lolliplot(list(B=gene2.rot, C=gene1.rot), 
          list(y=features2, z=features1),
          yaxis=FALSE, ylab = c("TUBB4B\nNP_006079","ben-1\nNP_497728"), xaxis = xaxis, legend = legend, type = "pie")

# #FF1F5B  #00CD6C  #009ADE  #AF58BA  #FFC61E  #F28522







# DYF-18 #


SNP <- c(26,114,173,187,197,207,210,245)
gene1 <- GRanges("chr1", IRanges(SNP, width=1, names=c("R26K","G114R","P173L","D187N","E197K","G207R","E210K","P245L")))
gene1$border <- sample(c("gray30"), length(SNP), replace=TRUE)
x <- c(50,0,0,0,0,50,0,50)
y <- c(0,100,100,100,100,50,100,0)
z <- c(50,0,0,0,0,0,0,50)
gene1$value1 <- x
gene1$value2 <- y
gene1$value3 <- z

## the length of the color should be no less than that of value1 or value2
gene1$color <- rep(list(c("#AF58BA","#FFC61E","#F28522")), length(SNP))
gene1.rot <- gene1
gene1.rot$label.parameter.rot <- 45
gene1.rot$SNPsideID<-"top"
features1 <- GRanges("chr1", IRanges(c(1, 8, 323), 
                                     width=c(0, 282, 0),
                                     names=c("", "Pkinase", "")))
features1$fill <- c("#FF1F5B", "#FF1F5B","#009ADE")
features1$height <- c(0, 0.09, 0)


legend <- c("#AF58BA","#FFC61E","#F28522") ## legend fill color
names(legend) <- c("Orthologous","Phenotypic/Pathogenic","Unknown")

xaxis1 <- c(1, 100, 200, 300, 325)

lolliplot(gene1.rot, 
          features1,
          yaxis=FALSE, ylab = "dyf-18\nNP_502232\n", xaxis = xaxis1, legend = legend, type = "pie")



SNP2<- c(30,92,180,215,253,283,342)
gene2 <- GRanges("chr1", IRanges(SNP2, width=1, names=c("R30K","D92N", "A180V", "G215E", "P253L", "R283Q", "K342N")))
#gene2$color <- sample.int(1, length(SNP2), replace=TRUE)
gene2$border <- sample(c("gray30"), length(SNP2), replace=TRUE)
#gene1$score <- sample.int(5, length(gene1), replace = TRUE)
x1 <- c(50,0,0,50,50,0,0)
y1 <- c(0,0,0,0,0,0,0)
z1 <- c(50,100,100,50,50,100,100)

gene2$value1 <- x1
gene2$value2 <- y1
gene2$value3 <- z1

## the length of the color should be no less than that of value1 or value2
gene2$color <- rep(list(c("#AF58BA","#FFC61E","#F28522")), length(SNP2))
gene2.rot <- gene2
gene2.rot$label.parameter.rot <- 45
gene2.rot$SNPsideID<-"top"
features2 <- GRanges("chr1", IRanges(c(1, 12, 346), 
                                     width=c(0, 283, 0),
                                     names=c("", "Pkinase", "")))
features2$fill <- c("#FF1F5B", "#FF1F5B","#009ADE")
features2$height <- c(0, 0.08, 0)
xaxis2 <- c(1, 100, 200, 300, 350)


lolliplot(gene2.rot, 
          features2,
          yaxis=FALSE, ylab = "CDK7\nNP_001790", xaxis = xaxis2, type = "pie")


lolliplot(list(B=gene2.rot, C=gene1.rot), 
          list(y=features2, z=features1),
          yaxis=FALSE, ylab = c("CDK7\nNP_001790","dyf-18\nNP_502232"), xaxis = xaxis2, legend = legend, type = "pie")




# OSM-3

SNP <- c(89,161,208,329,444,464,572,649,
         98,316,346,579,439)
gene1 <- GRanges("chr1", IRanges(SNP, width=1, names=c("T89I","P161L","S208L","N329H","G444E","A464T","W572*","S649N",
                                                       "Q98P","S316F","Q346*","W579*","Q439*")))
gene1$border <- sample(c("gray30"), length(SNP), replace=TRUE)
x <- c(50,50,50,50,50,50,50,50,0,0,0,0,0)
y <- c(0,0,0,50,50,0,0,0,100,100,100,100,100)
z <- c(50,50,50,0,0,50,50,50,0,0,0,0,0)
gene1$value1 <- x
gene1$value2 <- y
gene1$value3 <- z

## the length of the color should be no less than that of value1 or value2
gene1$color <- rep(list(c("#AF58BA","#FFC61E","#F28522")), length(SNP))
gene1.rot <- gene1
gene1.rot$label.parameter.rot <- 45
gene1.rot$SNPsideID<-"top"
gene1.rot$cex<-1.3
features1 <- GRanges("chr1", IRanges(c(1, 10, 699), 
                                     width=c(0, 317, 0),
                                     names=c("", "Kinesin", "")))
features1$fill <- c("#FF1F5B", "#FF1F5B","#009ADE")
features1$height <- c(0, 0.12, 0)


legend <- c("#AF58BA","#FFC61E","#F28522") ## legend fill color
names(legend) <- c("Orthologous","Phenotypic/Pathogenic","Unknown")

xaxis1 <- c(1,100,200,300,400,500,600,700)

lolliplot(gene1.rot, 
          features1,
          yaxis=FALSE, ylab = "osm-3\nNP_001023308\n\n\n", xaxis = xaxis1, legend = legend, type = "pie", cex = 1)



SNP2<- c(93,165,212,337,754,774,886,975,
         586,511,748,402,32)
gene2 <- GRanges("chr1", IRanges(SNP2, width=1, names=c("T93K","P165S", "S212C", "N337K", "G754V", "A774T", "W886*","S975S",
                                                        "A586P","D511Y","L748M","V402M","C32G")))
#gene2$color <- sample.int(1, length(SNP2), replace=TRUE)
gene2$border <- sample(c("gray30"), length(SNP2), replace=TRUE)
#gene1$score <- sample.int(5, length(gene1), replace = TRUE)
x1 <- c(50,50,50,50,50,50,50,50,0,0,0,0,0)
y1 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0)
z1 <- c(50,50,50,50,50,50,50,50,100,100,100,100,100)

gene2$value1 <- x1
gene2$value2 <- y1
gene2$value3 <- z1

## the length of the color should be no less than that of value1 or value2
gene2$color <- rep(list(c("#AF58BA","#FFC61E","#F28522")), length(SNP2))
gene2.rot <- gene2
gene2.rot$label.parameter.rot <- 45
gene2.rot$SNPsideID<-"top"
features2 <- GRanges("chr1", IRanges(c(1, 11, 1029), 
                                     width=c(0, 324, 0),
                                     names=c("", "Kinesin", "")))
features2$fill <- c("#FF1F5B", "#FF1F5B","#009ADE")
features2$height <- c(0, 0.08, 0)
xaxis2 <- c(1,100,200,300,400,500,600,700,800,900,1000,1030)


lolliplot(gene2.rot, 
          features2,
          yaxis=FALSE, ylab = "KIF17\nNP_065867", xaxis = xaxis2, type = "pie")


lolliplot(list(B=gene2.rot, C=gene1.rot), 
          list(y=features2, z=features1),
          yaxis=FALSE, ylab = c("CDK7\nNP_001790","osm-3\nNP_001023308"), xaxis = xaxis2, legend = legend, type = "pie", rescale = TRUE)




# MOUSE #

# ROS1



SNP <- c(79,220,238,394,490,531,669,763,766,907,1008,1227,1277,1402,1470,1500,1877,2050,2187,47,1158,2111)
gene1 <- GRanges("chr1", IRanges(SNP, width=1, names=c("D79A","T220A","Y238F","I394F","T490M","D531G","T669P","Y763H","Q766*",
                                                       "S907T","S1008T","F1227C","Y1277N","I1402F","T1470K","M1500T","E1877D","V2050A","M2187K",
                                                       "D47G","N1158Y","R2111G")))
gene1$border <- sample(c("gray30"), length(SNP), replace=TRUE)
x <- c(0,50,0,50,50,50,0,0,50,0,0,0,50,0,0,0,0,50,0,0,0,0)
y <- c(100,0,100,50,50,50,100,100,50,100,100,100,50,100,100,100,100,0,100,0,0,0)
z <- c(0,50,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,50,0,100,100,100)
gene1$value1 <- x
gene1$value2 <- y
gene1$value3 <- z

## the length of the color should be no less than that of value1 or value2
gene1$color <- rep(list(c("#AF58BA","#FFC61E","#F28522")), length(SNP))
gene1.rot <- gene1
gene1.rot$label.parameter.rot <- 45
gene1.rot$SNPsideID<-"top"
features1 <- GRanges("chr1", IRanges(c(1, 111, 206, 1650, 1938, 2340), 
                                     width=c(0, 75, 79, 84, 270, 0),
                                     names=c("", "fn3", "fn3", "fn3", "PK_Tyr_Ser-Thr", "")))
features1$fill <- c("#FF1F5B", "#FF1F5B","#FF1F5B","#FF1F5B","#009ADE","#FF1F5B")
features1$height <- c(0, 0.2, 0.2, 0.2, 0.2, 0)

legend <- c("#AF58BA","#FFC61E","#F28522") ## legend fill color
names(legend) <- c("Orthologous","Phenotypic/Pathogenic","Unknown")
xaxis1 <- c(1, 100, 200, 300, 325)





SNP2<- c(210,384,480,521,770,1281,2057,365,2032)
gene2 <- GRanges("chr1", IRanges(SNP2, width=1, names=c("T210A","I384F","T480I","D521V","Q770P","Y1281D","V2057I","G365A","G2032R")))
#gene2$color <- sample.int(1, length(SNP2), replace=TRUE)
gene2$border <- sample(c("gray30"), length(SNP2), replace=TRUE)
#gene1$score <- sample.int(5, length(gene1), replace = TRUE)
x1 <- c(50,50,50,50,50,50,50,0,0)
y1 <- c(0,0,0,0,0,0,0,100,100)
z1 <- c(50,50,50,50,50,50,50,0,0)

gene2$value1 <- x1
gene2$value2 <- y1
gene2$value3 <- z1

## the length of the color should be no less than that of value1 or value2
gene2$color <- rep(list(c("#AF58BA","#FFC61E","#F28522")), length(SNP2))
#gene2$color <- rep(list(c("#009ADE","#FFC61E","#FF1F5B")), length(SNP2))
gene2.rot <- gene2
gene2.rot$label.parameter.rot <- 45
gene2.rot$SNPsideID<-"top"
features2 <- GRanges("chr1", IRanges(c(1, 100, 196, 1659, 1945, 2347), 
                                     width=c(0, 76, 79, 82, 270, 0),
                                     names=c("", "fn3", "fn3", "fn3", "PK_Tyr_Ser-Thr", "")))
features2$fill <- c("#FF1F5B", "#FF1F5B","#FF1F5B","#FF1F5B","#009ADE","#FF1F5B")
#features2$fill <- c("#F28522", "#F28522","#F28522","#F28522","#AF58BA","#FF1F5B")
features2$height <- c(0, 0.2, 0.2, 0.2, 0.2, 0)
xaxis2 <- c(1, 500, 1000, 1500, 2000, 2350)


lolliplot(list(B=gene2.rot, C=gene1.rot), 
          list(y=features2, z=features1),
          yaxis=FALSE, ylab = c("ROS1\nNP_002935","Ros1\nNP_035412"), xaxis = xaxis2, legend = legend, type = "pie")



# Tg



SNP <- c(207,1342,2154,2212,2249,38,53,233,785,1256,1507,1640,1856,2110,164,1555,1993)
gene1 <- GRanges("chr1", IRanges(SNP, width=1, names=c("D207G","V1342I","F2154L","S2212P","N2249D","Q38L","C53Y","S233T","Y785N",
                                                       "S1256G","T1507A","T1640I","V1856L","V2110A","S164G","T1555S","C1993R")))
gene1$border <- sample(c("gray30"), length(SNP), replace=TRUE)
x <- c(50,50,50,50,50,0,0,0,0,0,0,0,0,0,50,50,50)
y <- c(50,50,50,50,50,100,100,100,100,100,100,100,100,100,0,0,0)
z <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,50,50,50)
gene1$value1 <- x
gene1$value2 <- y
gene1$value3 <- z

## the length of the color should be no less than that of value1 or value2
gene1$color <- rep(list(c("#AF58BA","#FFC61E","#F28522")), length(SNP))
gene1.rot <- gene1
gene1.rot$label.parameter.rot <- 45
gene1.rot$SNPsideID<-"top"
features1 <- GRanges("chr1", IRanges(c(1, 34, 97, 169, 301, 598, 662, 743, 879, 1007, 1150, 1464, 2196, 2766), 
                                     width=c(0, 59, 83, 84, 58, 60, 64, 55, 43, 67, 61, 45, 521, 0),
                                     names=c("", "Thyroglobulin_1", "Thyroglobulin_1", "Thyroglobulin_1",
                                             "Thyroglobulin_1", "Thyroglobulin_1", "Thyroglobulin_1", "Thyroglobulin_1",
                                             "Thyroglobulin_1", "Thyroglobulin_1", "Thyroglobulin_1", "Ephrin_rec_like", 
                                             "COesterase", "")))
features1$fill <- c("#FF1F5B","#FF1F5B","#FF1F5B","#FF1F5B","#FF1F5B","#FF1F5B","#FF1F5B","#FF1F5B","#FF1F5B",
                    "#FF1F5B","#FF1F5B","#009ADE","#A0B1BA","#FF1F5B")
features1$height <- c(0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0)

legend <- c("#AF58BA","#FFC61E","#F28522","#7CCBA2","#9EC9E2") ## legend fill color
names(legend) <- c("Orthologous","Phenotypic/Pathogenic","Unknown","Benign/Likely benign","VUS")

xaxis1 <- c(1, 100, 200, 300, 325)




SNP2<- c(37,206,232,1255,1342,2156,2213,2250,160,1077,1245,1897,1997,2242,2264,2336,2341,2375,163,1556,1996)
gene2 <- GRanges("chr1", IRanges(SNP2, width=1, names=c("Q37H","D206Y","S232C","S1255R","V1342I","F2156L","S2213*","N2250H",
                                                        "C160S","C1077R","C1245R","C1897Y","C1997S","R2242H","C2264Y",
                                                        "R2336Q","G2341S","G2375R","S163N","T1556T","C1996S")))
#gene2$color <- sample.int(1, length(SNP2), replace=TRUE)
gene2$border <- sample(c("gray30"), length(SNP2), replace=TRUE)
#gene1$score <- sample.int(5, length(gene1), replace = TRUE)
x1 <- c(0,50,0,0,50,50,50,50,0,0,0,0,0,0,0,0,0,0,50,50,50)
y1 <- c(0,0,0,0,0,0,0,0,100,100,100,100,100,100,100,100,100,100,0,0,50)
z1 <- c(100,50,100,100,50,50,50,50,0,0,0,0,0,0,0,0,0,0,0,0,0)
t1 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,50,0)
k1 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,50,0,0)

gene2$value1 <- x1
gene2$value2 <- y1
gene2$value3 <- z1
gene2$value4 <- t1
gene2$value5 <- k1

## the length of the color should be no less than that of value1 or value2
gene2$color <- rep(list(c("#AF58BA","#FFC61E","#F28522","#7CCBA2","#9EC9E2")), length(SNP2))
#gene2$color <- rep(list(c("#009ADE","#FFC61E","#FF1F5B")), length(SNP2))
gene2.rot <- gene2
gene2.rot$label.parameter.rot <- 45
gene2.rot$SNPsideID<-"top"
features2 <- GRanges("chr1", IRanges(c(1, 34, 96, 169, 301, 597, 662, 730, 1006, 1077, 1149, 1465, 2197, 2768), 
                                     width=c(0, 58, 83, 83, 57, 61, 64, 71, 67, 72, 61, 45, 521, 0),
                                     names=c("", "Thyroglobulin_1", "Thyroglobulin_1", "Thyroglobulin_1",
                                             "Thyroglobulin_1", "Thyroglobulin_1", "Thyroglobulin_1", "Thyroglobulin_1",
                                             "Thyroglobulin_1", "Thyroglobulin_1", "Thyroglobulin_1", "Ephrin_rec_like", 
                                             "COesterase", "")))
features2$fill <- c("#FF1F5B","#FF1F5B","#FF1F5B","#FF1F5B","#FF1F5B","#FF1F5B","#FF1F5B","#FF1F5B","#FF1F5B",
                    "#FF1F5B","#FF1F5B","#009ADE","#A0B1BA","#FF1F5B")
#features2$fill <- c("#F28522", "#F28522","#F28522","#F28522","#AF58BA","#FF1F5B")
features2$height <- c(0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0)
xaxis2 <- c(1, 500, 1000, 1500, 2000, 2350)



lolliplot(list(B=gene2.rot, C=gene1.rot), 
          list(y=features2, z=features1),
          yaxis=FALSE, ylab = c("TG\nNP_003226","Tg\nNP_033401"), xaxis = xaxis2, legend = legend, type = "pie")




