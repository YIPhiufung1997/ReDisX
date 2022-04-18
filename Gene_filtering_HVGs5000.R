
library(readr)
library(edgeR)
library(limma)
library(Glimma)
library(org.Hs.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)
library(pheatmap)
library(clusterProfiler)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(dplyr)
library(grid)
setwd("/home/yip/Documents/29_12_2021_ReDisX_analysis")

phenotype <- read_csv("mergepheonotypeGSE59867and93272.csv")

edata=read.csv("mergedataGSE59867and93272.csv",header=TRUE)

gene_SYMBOL=read_csv('genename.csv')

phenotype_CAD=phenotype[1:436,]
phenotype_RA=phenotype[437:711,]


edata_CAD=edata[,2:437]
row.names(edata_CAD)=gene_SYMBOL$SYMBOL
edata_RA=edata[,438:712]
row.names(edata_RA)=gene_SYMBOL$SYMBOL

### filter top 500 HVGs CAD
myvars_CAD <- apply(edata_CAD,1, var,na.rm=TRUE)
myvars_CAD <- sort(myvars_CAD,decreasing=TRUE)
myvars_CAD <- myvars_CAD[1:5000]
edata_CAD_top5000 <- edata_CAD[names(myvars_CAD),]
dim(edata_CAD_top5000)

reorder_data_CAD=edata_CAD_top5000[ order(row.names(edata_CAD_top5000)), ]



### filter top 500 HVGs RA
myvars_RA <- apply(edata_RA,1, var,na.rm=TRUE)
myvars_RA <- sort(myvars_RA,decreasing=TRUE)
myvars_RA <- myvars_RA[1:5000]
edata_RA_top5000 <- edata_RA[names(myvars_RA),]
dim(edata_RA_top5000)

reorder_data_RA=edata_RA_top5000[ order(row.names(edata_RA_top5000)), ]



write_csv(edata_RA_top5000,file = "Top5000RA_data.csv")
write_csv(edata_CAD_top5000,file = "Top5000CAD_data.csv")

genename_RA_top5000=data.frame(row.names(edata_RA_top5000))
genename_CAD_top5000=data.frame(row.names(edata_CAD_top5000))
write_csv(genename_RA_top5000,file = "Top5000RA_genename.csv")
write_csv(genename_CAD_top5000,file = "Top5000CAD_genename.csv")


write_csv(phenotype_RA,file = "phenotype_RA_data.csv")
write_csv(phenotype_CAD,file = "phenotype_CAD_data.csv")

common_CAD<-edata_CAD_top5000[match(row.names(edata_RA_top5000),row.names(edata_CAD_top5000), nomatch=0),]
common_RA<-edata_RA_top5000[match(row.names(edata_CAD_top5000),row.names(edata_RA_top5000), nomatch=0),]





