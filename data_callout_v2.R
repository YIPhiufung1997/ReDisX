library(dplyr)
library(readr)
library(GEOquery)
library(clusterProfiler)
library(org.Hs.eg.db)
library(sva)

library(pamr)
library(limma)

setwd("/home/yip/Documents/Reclassification_customcdf")

all_no_CAD1 <- read_csv("GSE59867_RAW/all_no.na.csv")

all_no_RA1<-read.csv('GSE93272_RAW/all_no.na.csv')


geneSYM_in_CAD1<-data.frame(all_no_CAD1[,4])
geneSYM_in_RA1<-data.frame(all_no_RA1[,4])

colnames(geneSYM_in_RA1)<-"SYMBOL"
gene_in_common<-data.frame(SYMBOL=intersect(geneSYM_in_CAD1$SYMBOL,geneSYM_in_RA1$SYMBOL))


common_CAD<-all_no_CAD1[match(gene_in_common$SYMBOL,all_no_CAD1$SYMBOL, nomatch=0),]
common_RA<-all_no_RA1[match(gene_in_common$SYMBOL,all_no_RA1$SYMBOL, nomatch=0),]

gse_RA_dataset1 <- getGEO(filename = "GSE93272_customcdf_prerequire/GSE93272_family.soft.gz")
sample_id_RA<-colnames(all_no_RA1)[6:ncol(all_no_RA1)]
sample_id_RA<-gsub("_.*", "", sample_id_RA) 


diseasetype_GSE93272<-c()
for(id in sample_id_RA)
  diseasetype_GSE93272<-append(diseasetype_GSE93272,gse_RA_dataset1@gsms[[toString(id)]]@header[["source_name_ch1"]])


diseasetype_GSE93272<-gsub('Whole blood from healthy control', 'Control', diseasetype_GSE93272)
diseasetype_GSE93272<-gsub('Whole blood from rheumatoid arthritis', 'RA', diseasetype_GSE93272)
phenotype_GSE93272<-data.frame(sample_id=sample_id_RA,GSEseries=rep(1,length(sample_id_RA)),disease_type=diseasetype_GSE93272)

#########################################
######################GSE59867

gse_CAD_dataset1 <- getGEO(filename = "GSE59867_customcdf_prerequire/GSE59867_family.soft.gz")
sample_id_CAD<-colnames(all_no_CAD1)[6:ncol(all_no_CAD1)]
sample_id_CAD<-gsub("_.*", "", sample_id_CAD) 


diseasetype_GSE59867<-c()
for(id in sample_id_CAD)
  diseasetype_GSE59867<-append(diseasetype_GSE59867,gse_CAD_dataset1@gsms[[toString(id)]]@header[["title"]])

diseasetype_GSE59867<-gsub('[0-9]+', '', diseasetype_GSE59867)
diseasetype_GSE59867<-gsub('Patient Z,sampling', 'CAD', diseasetype_GSE59867)
diseasetype_GSE59867<-gsub('patient', '', diseasetype_GSE59867)

phenotype_GSE59867<-data.frame(sample_id=sample_id_CAD,GSEseries=rep(2,length(sample_id_CAD)),disease_type=diseasetype_GSE59867)


### merge two dataframe
sub_common_CAD<-common_CAD[,6:ncol(all_no_CAD1)]
total_data <- cbind(common_CAD[,6:ncol(all_no_CAD1)], common_RA[,6:ncol(all_no_RA1)]) 
total_pheonotype <- rbind(phenotype_GSE59867, phenotype_GSE93272)
colnames(total_data)<-total_pheonotype$sample_id

## remove missing

missing<-is.na(gene_in_common) 
missing_row<-which(missing == "TRUE")

gene_in_common_withoutmissing<-gene_in_common %>% dplyr::slice(-c(missing_row))
total_data_withoutmissing<-total_data %>% dplyr::slice(-c(missing_row))
row.names(total_data_withoutmissing)<-gene_in_common_withoutmissing$SYMBOL


write.csv(total_data_withoutmissing,file="mergedataGSE59867and93272.csv")
write.csv(total_pheonotype,file="mergepheonotypeGSE59867and93272.csv")
### sva batch effect remove

pheno = total_pheonotype
edata=total_data_withoutmissing
mod = model.matrix(~as.factor(disease_type), data=pheno)

mod0 = model.matrix(~as.factor(GSEseries),data=pheno)


svobj = sva(edata,mod,mod0)
