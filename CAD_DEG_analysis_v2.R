library(umap)
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

clustering_result <- read_csv("matlab_GSE59867_ReDisX/CAD_clustering_result_pluscontrol.csv")

edata=read.csv("matlab_GSE59867_ReDisX/CAD_expression_final_pluscontrol.csv",header=FALSE)

gene_SYMBOL=read_csv('Top5000CAD_genename.csv')

colnames(gene_SYMBOL)="SYMBOL"
colnames(clustering_result)[1]="sampleid"
colnames(clustering_result)[2]="final_clus"



clustering_result$final_clus=as.factor(clustering_result$final_clus)


clustering_result[,5]=clustering_result[,2]
colnames(clustering_result)[5]="cluster_for_limma"

clustering_result$cluster_for_limma<-gsub('1', 'ReDisX_CAD_clus1', clustering_result$cluster_for_limma)
clustering_result$cluster_for_limma<-gsub('2', 'ReDisX_CAD_clus2', clustering_result$cluster_for_limma)
clustering_result$cluster_for_limma<-gsub('3', 'ReDisX_CAD_clus3', clustering_result$cluster_for_limma)
clustering_result$cluster_for_limma<-gsub('4', 'ReDisX_CAD_clus4', clustering_result$cluster_for_limma)
clustering_result$cluster_for_limma<-gsub('5', 'ReDisX_CAD_clus5', clustering_result$cluster_for_limma)
clustering_result$cluster_for_limma<-gsub('6', 'Control', clustering_result$cluster_for_limma)


a=factor(clustering_result$cluster_for_limma)
disease_type_num=as.numeric(a)

row.names(edata)=gene_SYMBOL$SYMBOL
colnames(edata)=clustering_result$sampleid
edata[is.na(edata)] <- 0
y <- DGEList(edata)
y$genes <- gene_SYMBOL


#barplot(y$samples$lib.size[1:10],names=colnames(y)[1:10],las=2)
# Add a title to the plot
#title("Barplot of library sizes")




#pdf(file="result/tradis_MDS.pdf")
### MDS plot
#a=plotMDS(edata)
#dev.off()



#### traditional disease diagnosis
disease_target <- factor(clustering_result$cluster_for_limma, levels=c("ReDisX_CAD_clus1","ReDisX_CAD_clus2","ReDisX_CAD_clus3","ReDisX_CAD_clus4","ReDisX_CAD_clus5",'Control'))
design <- model.matrix(~ 0 + disease_target)





pdf(file="result_v2/ReDisX_GSE59867voom:Mean_variance_trend.pdf")
par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE,save.plot = TRUE)
dev.off()

write.csv(v[["E"]], 'CAD_matrix_after_voom.csv')

fit <- lmFit(v)
names(fit)

cont.matrix <- makeContrasts(ReDisX_CAD_clus1_Control_DEG= disease_targetReDisX_CAD_clus1-disease_targetControl, ReDisX_CAD_clus2_Control_DEG= disease_targetReDisX_CAD_clus2-disease_targetControl,
                             ReDisX_CAD_clus3_Control_DEG= disease_targetReDisX_CAD_clus3-disease_targetControl,
                             ReDisX_CAD_clus4_Control_DEG= disease_targetReDisX_CAD_clus4-disease_targetControl,
                             ReDisX_CAD_clus5_Control_DEG= disease_targetReDisX_CAD_clus5-disease_targetControl,
                                levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)


toplist_CAD=topTable(fit.cont,adjust.method="BH",p.value=0.05,number=20000)


summa.fit <-decideTests(fit.cont,method="separate",adjust.method="BH",p.value=0.05,lfc=0.05)
summary=data.frame(summa.fit)
summary(summa.fit)
write.csv(summary, 'ReDisX_DEG_result_CAD.csv')
write.csv(toplist_CAD, 'ReDisX_DEG_toplist_CAD.csv')

##### ReDisX diagnosis

#disease_target_ReDisX <- factor(clustering_result$final_clus, levels=c("1","2",'3','4'))
#design_ReDisX <- model.matrix(~ 0 + disease_target_ReDisX)


#pdf(file="result/ReDisX_voom_plot.pdf")
#par(mfrow=c(1,1))
#v_ReDisX <- voom(y,design_ReDisX,plot = TRUE,save.plot = TRUE)
#dev.off()
#pdf(file="result/ReDisX_boxplot.pdf")
#boxplot(v_ReDisX$E[1:10,1:10], xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
#abline(h=median(v_ReDisX$E),col="blue")
#dev.off()



#fit_ReDisX <- lmFit(v_ReDisX)
#names(fit_ReDisX)

#cont.matrix_ReDisX <- makeContrasts(ReDisXclus2_healthy_DEG= disease_target_ReDisX2-disease_target_ReDisX4,levels=design_ReDisX)
#fit.cont_ReDisX <- contrasts.fit(fit_ReDisX, cont.matrix_ReDisX)
#fit.cont_ReDisX <- eBayes(fit.cont_ReDisX)
#dim(fit.cont_ReDisX)


#toplist_tradisease_ReDisX=topTable(fit.cont_ReDisX,adjust.method="none",p.value=1,number=200000)


#summa.fit_ReDisX <-decideTests(fit.cont_ReDisX,method="global",adjust.method="none",p.value=0.05)
#summary_ReDisX=data.frame(summa.fit_ReDisX)
#summary(summa.fit_ReDisX)
#write.csv(summary_ReDisX, 'ReDisX_DEG_result.csv')




####  DEG for RA, CAD, ReDisXclus2

degs_index_ReDisX_clus1= which(summary[,1] %in% c(1,-1))
degs_index_ReDisX_clus2= which(summary[,2] %in% c(1,-1))
degs_index_ReDisX_clus3= which(summary[,3] %in% c(1,-1))
degs_index_ReDisX_clus4= which(summary[,4] %in% c(1,-1))
degs_index_ReDisX_clus5= which(summary[,5] %in% c(1,-1))
#degs_index_ReDisX2= which(summary_ReDisX[,1] %in% c(1,-1))
#degs_index_ReDisX1= which(summary_ReDisX[,2] %in% c(1,-1))
#degs_index_ReDisX3= which(summary_ReDisX[,3] %in% c(1,-1))

### distinct deg clus1
distinct_ReDisX_clus1DEG=setdiff(degs_index_ReDisX_clus1,degs_index_ReDisX_clus2)
distinct_ReDisX_clus1DEG=setdiff(distinct_ReDisX_clus1DEG,degs_index_ReDisX_clus3)
distinct_ReDisX_clus1DEG=setdiff(distinct_ReDisX_clus1DEG,degs_index_ReDisX_clus4)
distinct_ReDisX_clus1DEG=setdiff(distinct_ReDisX_clus1DEG,degs_index_ReDisX_clus5)
### distinct deg clus2
distinct_ReDisX_clus2DEG=setdiff(degs_index_ReDisX_clus2,degs_index_ReDisX_clus1)
distinct_ReDisX_clus2DEG=setdiff(distinct_ReDisX_clus2DEG,degs_index_ReDisX_clus3)
distinct_ReDisX_clus2DEG=setdiff(distinct_ReDisX_clus2DEG,degs_index_ReDisX_clus4)
distinct_ReDisX_clus2DEG=setdiff(distinct_ReDisX_clus2DEG,degs_index_ReDisX_clus5)
### distinct deg clus3
distinct_ReDisX_clus3DEG=setdiff(degs_index_ReDisX_clus3,degs_index_ReDisX_clus1)
distinct_ReDisX_clus3DEG=setdiff(distinct_ReDisX_clus3DEG,degs_index_ReDisX_clus2)
distinct_ReDisX_clus3DEG=setdiff(distinct_ReDisX_clus3DEG,degs_index_ReDisX_clus4)
distinct_ReDisX_clus3DEG=setdiff(distinct_ReDisX_clus3DEG,degs_index_ReDisX_clus5)
### distinct deg clus4
distinct_ReDisX_clus4DEG=setdiff(degs_index_ReDisX_clus4,degs_index_ReDisX_clus1)
distinct_ReDisX_clus4DEG=setdiff(distinct_ReDisX_clus4DEG,degs_index_ReDisX_clus2)
distinct_ReDisX_clus4DEG=setdiff(distinct_ReDisX_clus4DEG,degs_index_ReDisX_clus3)
distinct_ReDisX_clus4DEG=setdiff(distinct_ReDisX_clus4DEG,degs_index_ReDisX_clus5)
### distinct deg clus5
distinct_ReDisX_clus5DEG=setdiff(degs_index_ReDisX_clus5,degs_index_ReDisX_clus1)
distinct_ReDisX_clus5DEG=setdiff(distinct_ReDisX_clus5DEG,degs_index_ReDisX_clus2)
distinct_ReDisX_clus5DEG=setdiff(distinct_ReDisX_clus5DEG,degs_index_ReDisX_clus3)
distinct_ReDisX_clus5DEG=setdiff(distinct_ReDisX_clus5DEG,degs_index_ReDisX_clus4)



### save distinct DEG
write.csv(gene_SYMBOL[distinct_ReDisX_clus1DEG,],file = 'result_v2/GSE59867_result/distict_gene_DisRedXclus1_GSE59867.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus2DEG,],file = 'result_v2/GSE59867_result/distict_gene_DisRedXclus2_GSE59867.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus3DEG,],file = 'result_v2/GSE59867_result/distict_gene_DisRedXclus3_GSE59867.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus4DEG,],file = 'result_v2/GSE59867_result/distict_gene_DisRedXclus4_GSE59867.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus5DEG,],file = 'result_v2/GSE59867_result/distict_gene_DisRedXclus5_GSE59867.csv')


### MD plot after DE


pdf(file="result_v2/ReDisXclus1_GSE59867_MD.pdf")


plotMD(fit.cont,coef=1,status=summa.fit[,"ReDisX_CAD_clus1_Control_DEG"], values = c(-1, 1), hl.col=c("blue","red"))
dev.off()
pdf(file="result_v2/ReDisXclus2_GSE59867_MD.pdf")
plotMD(fit.cont,coef=2,status=summa.fit[,"ReDisX_CAD_clus2_Control_DEG"], values = c(-1, 1), hl.col=c("blue","red"))
dev.off()
pdf(file="result_v2/ReDisXclus3_GSE59867_MD.pdf")
plotMD(fit.cont,coef=3,status=summa.fit[,"ReDisX_CAD_clus3_Control_DEG"], values = c(-1, 1), hl.col=c("blue","red"))
dev.off()
pdf(file="result_v2/ReDisXclus4_GSE59867_MD.pdf")
plotMD(fit.cont,coef=4,status=summa.fit[,"ReDisX_CAD_clus4_Control_DEG"], values = c(-1, 1), hl.col=c("blue","red"))
dev.off()
pdf(file="result_v2/ReDisXclus5_GSE59867_MD.pdf")
plotMD(fit.cont,coef=5,status=summa.fit[,"ReDisX_CAD_clus5_Control_DEG"], values = c(-1, 1), hl.col=c("blue","red"))
dev.off()
# For the volcano plot we have to specify how many of the top genes to highlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
pdf(file="result_v2/ReDisXclus1_GSE59867_volcano.pdf")
volcanoplot(fit.cont,coef=1,highlight=50,names=fit.cont$genes$SYMBOL, main="ReDisX_clus1_Control")
dev.off()
pdf(file="result_v2/ReDisXclus2_GSE59867_volcano.pdf")
volcanoplot(fit.cont,coef=2,highlight=50,names=fit.cont$genes$SYMBOL, main="ReDisX_clus2_Control")
dev.off()
pdf(file="result_v2/ReDisXclus3_GSE59867_volcano.pdf")
volcanoplot(fit.cont,coef=3,highlight=50,names=fit.cont$genes$SYMBOL, main="ReDisX_clus3_Control")
dev.off()
pdf(file="result_v2/ReDisXclus4_GSE59867_volcano.pdf")
volcanoplot(fit.cont,coef=4,highlight=50,names=fit.cont$genes$SYMBOL, main="ReDisX_clus4_Control")
dev.off()
pdf(file="result_v2/ReDisXclus5_GSE59867_volcano.pdf")
volcanoplot(fit.cont,coef=5,highlight=50,names=fit.cont$genes$SYMBOL, main="ReDisX_clus5_Control")
dev.off()

png(file="result_v2/ReDisX_clus2_MD+volcano_GSE59867.png",
    width=600, height=350)
#pdf(file="result/tradisRA_MD+volcano.pdf")
# We want to highlight the significant genes. We can get this from decideTests.
par(mfrow=c(1,2))
plotMD(fit.cont,coef=1,status=summa.fit[,"ReDisX_clus2_Control_DEG"], values = c(-1, 1), hl.col=c("blue","red"))

# For the volcano plot we have to specify how many of the top genes to highlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
volcanoplot(fit.cont,coef=1,highlight=100,names=fit.cont$genes$SYMBOL, main="ReDisX_clus2_Control_type")

dev.off()





### heatmap  ReDisX_clus1



degs_clus1=gene_SYMBOL[distinct_ReDisX_clus1DEG,]


vst_clus1 = v[["E"]]

data_for_clus1 = vst_clus1[degs_clus1$SYMBOL,]
data_for_clus1_label1=data_for_clus1[,disease_target=="ReDisX_CAD_clus1"]
data_for_clus1_label2=data_for_clus1[,disease_target=="ReDisX_CAD_clus2"]
data_for_clus1_label3=data_for_clus1[,disease_target=="ReDisX_CAD_clus3"]
data_for_clus1_label4=data_for_clus1[,disease_target=="ReDisX_CAD_clus4"]
data_for_clus1_label5=data_for_clus1[,disease_target=="ReDisX_CAD_clus5"]
data_for_clus1_label6=data_for_clus1[,disease_target=="Control"]
data_for_clus1_all =cbind(data_for_clus1_label1,data_for_clus1_label2,data_for_clus1_label3,data_for_clus1_label4,data_for_clus1_label5,data_for_clus1_label6)





my_breaks = c(seq(-2, -0.5, length.out=25),
              0,
              seq(0.5, 2, length.out=25))
my_colors = c("green", "black", "red")
my_colors = colorRampPalette(my_colors)(50)


annotation_col_clus1 = data.frame(
  Cluster = disease_target)
rownames(annotation_col_clus1)=clustering_result$sampleid


pdf(file="result_v2/ReDisX_clus1_GSE59867.pdf",width=30, height=30)
pheatmap(data_for_clus1_all, fontsize_row=4,fontsize_col = 2, scale='row',breaks = my_breaks, annotation_col=annotation_col_clus1, cluster_cols = F,color = my_colors,col_split = annotation_col_clus1$Cluster,
         gaps_col = cumsum(c(30,70,12,97,181))
         )
dev.off()







### heatmap  ReDisX_clus2



degs_clus2=gene_SYMBOL[distinct_ReDisX_clus2DEG,]


vst_clus2 = v[["E"]]

data_for_clus2 = vst_clus2[degs_clus2$SYMBOL,]
data_for_clus2_label1=data_for_clus2[,disease_target=="ReDisX_CAD_clus1"]
data_for_clus2_label2=data_for_clus2[,disease_target=="ReDisX_CAD_clus2"]
data_for_clus2_label3=data_for_clus2[,disease_target=="ReDisX_CAD_clus3"]
data_for_clus2_label4=data_for_clus2[,disease_target=="ReDisX_CAD_clus4"]
data_for_clus2_label5=data_for_clus2[,disease_target=="ReDisX_CAD_clus5"]
data_for_clus2_label6=data_for_clus2[,disease_target=="Control"]
data_for_clus2_all =cbind(data_for_clus2_label1,data_for_clus2_label2,data_for_clus2_label3,data_for_clus2_label4,data_for_clus2_label5,data_for_clus2_label6)





my_breaks = c(seq(-2, -0.5, length.out=25),
              0,
              seq(0.5, 2, length.out=25))
my_colors = c("green", "black", "red")
my_colors = colorRampPalette(my_colors)(50)


annotation_col_clus2 = data.frame(
  Cluster = disease_target)
rownames(annotation_col_clus2)=clustering_result$sampleid


pdf(file="result_v2/ReDisX_clus2_GSE59867.pdf",width=30, height=30)
pheatmap(data_for_clus2_all, fontsize_row=4,fontsize_col = 2, scale='row',breaks = my_breaks, annotation_col=annotation_col_clus2, cluster_cols = F,color = my_colors,col_split = annotation_col_clus2$Cluster,
         gaps_col = cumsum(c(30,70,12,97,181))
)
dev.off()
### heatmap  ReDisX_clus3



degs_clus3=gene_SYMBOL[distinct_ReDisX_clus3DEG,]


vst_clus3 = v[["E"]]

data_for_clus3 = vst_clus3[degs_clus3$SYMBOL,]
data_for_clus3_label1=data_for_clus3[,disease_target=="ReDisX_CAD_clus1"]
data_for_clus3_label2=data_for_clus3[,disease_target=="ReDisX_CAD_clus2"]
data_for_clus3_label3=data_for_clus3[,disease_target=="ReDisX_CAD_clus3"]
data_for_clus3_label4=data_for_clus3[,disease_target=="ReDisX_CAD_clus4"]
data_for_clus3_label5=data_for_clus3[,disease_target=="ReDisX_CAD_clus5"]
data_for_clus3_label6=data_for_clus3[,disease_target=="Control"]
data_for_clus3_all =cbind(data_for_clus3_label1,data_for_clus3_label2,data_for_clus3_label3,data_for_clus3_label4,data_for_clus3_label5,data_for_clus3_label6)





my_breaks = c(seq(-2, -0.5, length.out=25),
              0,
              seq(0.5, 2, length.out=25))
my_colors = c("green", "black", "red")
my_colors = colorRampPalette(my_colors)(50)


annotation_col_clus3 = data.frame(
  Cluster = disease_target)
rownames(annotation_col_clus3)=clustering_result$sampleid


pdf(file="result_v2/ReDisX_clus3_GSE59867.pdf",width=30, height=30)
pheatmap(data_for_clus3_all, fontsize_row=4,fontsize_col = 2, scale='row',breaks = my_breaks, annotation_col=annotation_col_clus3, cluster_cols = F,color = my_colors,col_split = annotation_col_clus3$Cluster,
         gaps_col = cumsum(c(30,70,12,97,181))
)
dev.off()


### heatmap  ReDisX_clus4



degs_clus4=gene_SYMBOL[distinct_ReDisX_clus4DEG,]


vst_clus4 = v[["E"]]

data_for_clus4 = vst_clus4[degs_clus4$SYMBOL,]
data_for_clus4_label1=data_for_clus4[,disease_target=="ReDisX_CAD_clus1"]
data_for_clus4_label2=data_for_clus4[,disease_target=="ReDisX_CAD_clus2"]
data_for_clus4_label3=data_for_clus4[,disease_target=="ReDisX_CAD_clus3"]
data_for_clus4_label4=data_for_clus4[,disease_target=="ReDisX_CAD_clus4"]
data_for_clus4_label5=data_for_clus4[,disease_target=="ReDisX_CAD_clus5"]
data_for_clus4_label6=data_for_clus4[,disease_target=="Control"]
data_for_clus4_all =cbind(data_for_clus4_label1,data_for_clus4_label2,data_for_clus4_label3,data_for_clus4_label4,data_for_clus4_label5,data_for_clus4_label6)





my_breaks = c(seq(-2, -0.5, length.out=25),
              0,
              seq(0.5, 2, length.out=25))
my_colors = c("green", "black", "red")
my_colors = colorRampPalette(my_colors)(50)


annotation_col_clus4 = data.frame(
  Cluster = disease_target)
rownames(annotation_col_clus4)=clustering_result$sampleid


pdf(file="result_v2/ReDisX_clus4_GSE59867.pdf",width=30, height=30)
pheatmap(data_for_clus4_all, fontsize_row=4,fontsize_col = 2, scale='row',breaks = my_breaks, annotation_col=annotation_col_clus4, cluster_cols = F,color = my_colors,col_split = annotation_col_clus4$Cluster,
         gaps_col = cumsum(c(30,70,12,97,181))
)
dev.off()

### heatmap  ReDisX_clus5



degs_clus5=gene_SYMBOL[distinct_ReDisX_clus5DEG,]


vst_clus5 = v[["E"]]

data_for_clus5 = vst_clus5[degs_clus5$SYMBOL,]
data_for_clus5_label1=data_for_clus5[,disease_target=="ReDisX_CAD_clus1"]
data_for_clus5_label2=data_for_clus5[,disease_target=="ReDisX_CAD_clus2"]
data_for_clus5_label3=data_for_clus5[,disease_target=="ReDisX_CAD_clus3"]
data_for_clus5_label4=data_for_clus5[,disease_target=="ReDisX_CAD_clus4"]
data_for_clus5_label5=data_for_clus5[,disease_target=="ReDisX_CAD_clus5"]
data_for_clus5_label6=data_for_clus5[,disease_target=="Control"]
data_for_clus5_all =cbind(data_for_clus5_label1,data_for_clus5_label2,data_for_clus5_label3,data_for_clus5_label4,data_for_clus5_label5,data_for_clus5_label6)





my_breaks = c(seq(-2, -0.5, length.out=25),
              0,
              seq(0.5, 2, length.out=25))
my_colors = c("green", "black", "red")
my_colors = colorRampPalette(my_colors)(50)


annotation_col_clus5 = data.frame(
  Cluster = disease_target)
rownames(annotation_col_clus5)=clustering_result$sampleid


pdf(file="result_v2/ReDisX_clus5_GSE59867.pdf",width=30, height=30)
pheatmap(data_for_clus5_all, fontsize_row=4,fontsize_col = 2, scale='row',breaks = my_breaks, annotation_col=annotation_col_clus5, cluster_cols = F,color = my_colors,col_split = annotation_col_clus5$Cluster,
         gaps_col = cumsum(c(30,70,12,97,181))
)
dev.off()


#### hubgenes in query heatmap


summa_p0_005.fit <-decideTests(fit.cont,method="separate",adjust.method="BH",p.value=0.005)
summary_p0_005=data.frame(summa_p0_005.fit)
summary(summa_p0_005.fit)
write.csv(summary_p0_005, 'result_v2/GSE59867_result/p_005_distinct_result/ReDisX_DEG_result_CAD_p0_005.csv')



degs_index_ReDisX_clus1_p0_005= which(summary_p0_005[,1] %in% c(1))
degs_index_ReDisX_clus2_p0_005= which(summary_p0_005[,2] %in% c(1))
degs_index_ReDisX_clus3_p0_005= which(summary_p0_005[,3] %in% c(1))
degs_index_ReDisX_clus4_p0_005= which(summary_p0_005[,4] %in% c(1))
degs_index_ReDisX_clus5_p0_005= which(summary_p0_005[,5] %in% c(1))


### distinct deg clus1
distinct_ReDisX_clus1DEG_p0_005=setdiff(degs_index_ReDisX_clus1_p0_005,degs_index_ReDisX_clus2_p0_005)
distinct_ReDisX_clus1DEG_p0_005=setdiff(distinct_ReDisX_clus1DEG_p0_005,degs_index_ReDisX_clus3_p0_005)
distinct_ReDisX_clus1DEG_p0_005=setdiff(distinct_ReDisX_clus1DEG_p0_005,degs_index_ReDisX_clus4_p0_005)
distinct_ReDisX_clus1DEG_p0_005=setdiff(distinct_ReDisX_clus1DEG_p0_005,degs_index_ReDisX_clus5_p0_005)
### distinct deg clus2
distinct_ReDisX_clus2DEG_p0_005=setdiff(degs_index_ReDisX_clus2_p0_005,degs_index_ReDisX_clus1_p0_005)
distinct_ReDisX_clus2DEG_p0_005=setdiff(distinct_ReDisX_clus2DEG_p0_005,degs_index_ReDisX_clus3_p0_005)
distinct_ReDisX_clus2DEG_p0_005=setdiff(distinct_ReDisX_clus2DEG_p0_005,degs_index_ReDisX_clus4_p0_005)
distinct_ReDisX_clus2DEG_p0_005=setdiff(distinct_ReDisX_clus2DEG_p0_005,degs_index_ReDisX_clus5_p0_005)
### distinct deg clus3
distinct_ReDisX_clus3DEG_p0_005=setdiff(degs_index_ReDisX_clus3_p0_005,degs_index_ReDisX_clus1_p0_005)
distinct_ReDisX_clus3DEG_p0_005=setdiff(distinct_ReDisX_clus3DEG_p0_005,degs_index_ReDisX_clus2_p0_005)
distinct_ReDisX_clus3DEG_p0_005=setdiff(distinct_ReDisX_clus3DEG_p0_005,degs_index_ReDisX_clus4_p0_005)
distinct_ReDisX_clus3DEG_p0_005=setdiff(distinct_ReDisX_clus3DEG_p0_005,degs_index_ReDisX_clus5_p0_005)
### distinct deg clus4
distinct_ReDisX_clus4DEG_p0_005=setdiff(degs_index_ReDisX_clus4_p0_005,degs_index_ReDisX_clus1_p0_005)
distinct_ReDisX_clus4DEG_p0_005=setdiff(distinct_ReDisX_clus4DEG_p0_005,degs_index_ReDisX_clus2_p0_005)
distinct_ReDisX_clus4DEG_p0_005=setdiff(distinct_ReDisX_clus4DEG_p0_005,degs_index_ReDisX_clus3_p0_005)
distinct_ReDisX_clus4DEG_p0_005=setdiff(distinct_ReDisX_clus4DEG_p0_005,degs_index_ReDisX_clus5_p0_005)
### distinct deg clus5
distinct_ReDisX_clus5DEG_p0_005=setdiff(degs_index_ReDisX_clus5_p0_005,degs_index_ReDisX_clus1_p0_005)
distinct_ReDisX_clus5DEG_p0_005=setdiff(distinct_ReDisX_clus5DEG_p0_005,degs_index_ReDisX_clus2_p0_005)
distinct_ReDisX_clus5DEG_p0_005=setdiff(distinct_ReDisX_clus5DEG_p0_005,degs_index_ReDisX_clus3_p0_005)
distinct_ReDisX_clus5DEG_p0_005=setdiff(distinct_ReDisX_clus5DEG_p0_005,degs_index_ReDisX_clus4_p0_005)



### save distinct DEG
write.csv(gene_SYMBOL[distinct_ReDisX_clus1DEG_p0_005,],file = 'result_v2/GSE59867_result/p_005_distinct_result/distict_gene_DisRedXclus1_GSE59867.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus2DEG_p0_005,],file = 'result_v2/GSE59867_result/p_005_distinct_result/distict_gene_DisRedXclus2_GSE59867.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus3DEG_p0_005,],file = 'result_v2/GSE59867_result/p_005_distinct_result/distict_gene_DisRedXclus3_GSE59867.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus4DEG_p0_005,],file = 'result_v2/GSE59867_result/p_005_distinct_result/distict_gene_DisRedXclus4_GSE59867.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus5DEG_p0_005,],file = 'result_v2/GSE59867_result/p_005_distinct_result/distict_gene_DisRedXclus5_GSE59867.csv')


############## distinct DEGs down regulated
degs_index_ReDisX_clus1_p0_05_down= which(summary_p0_005[,1] %in% c(-1))
degs_index_ReDisX_clus2_p0_05_down= which(summary_p0_005[,2] %in% c(-1))
degs_index_ReDisX_clus3_p0_05_down= which(summary_p0_005[,3] %in% c(-1))
degs_index_ReDisX_clus4_p0_05_down= which(summary_p0_005[,4] %in% c(-1))
degs_index_ReDisX_clus5_p0_05_down= which(summary_p0_005[,5] %in% c(-1))


### distinct deg clus1
distinct_ReDisX_clus1DEG_p0_05_down=setdiff(degs_index_ReDisX_clus1_p0_05_down,degs_index_ReDisX_clus2_p0_05_down)
distinct_ReDisX_clus1DEG_p0_05_down=setdiff(distinct_ReDisX_clus1DEG_p0_05_down,degs_index_ReDisX_clus3_p0_05_down)
distinct_ReDisX_clus1DEG_p0_05_down=setdiff(distinct_ReDisX_clus1DEG_p0_05_down,degs_index_ReDisX_clus4_p0_05_down)
distinct_ReDisX_clus1DEG_p0_05_down=setdiff(distinct_ReDisX_clus1DEG_p0_05_down,degs_index_ReDisX_clus5_p0_05_down)
### distinct deg clus2
distinct_ReDisX_clus2DEG_p0_05_down=setdiff(degs_index_ReDisX_clus2_p0_05_down,degs_index_ReDisX_clus1_p0_05_down)
distinct_ReDisX_clus2DEG_p0_05_down=setdiff(distinct_ReDisX_clus2DEG_p0_05_down,degs_index_ReDisX_clus3_p0_05_down)
distinct_ReDisX_clus2DEG_p0_05_down=setdiff(distinct_ReDisX_clus2DEG_p0_05_down,degs_index_ReDisX_clus4_p0_05_down)
distinct_ReDisX_clus2DEG_p0_05_down=setdiff(distinct_ReDisX_clus2DEG_p0_05_down,degs_index_ReDisX_clus5_p0_05_down)
### distinct deg clus3
distinct_ReDisX_clus3DEG_p0_05_down=setdiff(degs_index_ReDisX_clus3_p0_05_down,degs_index_ReDisX_clus1_p0_05_down)
distinct_ReDisX_clus3DEG_p0_05_down=setdiff(distinct_ReDisX_clus3DEG_p0_05_down,degs_index_ReDisX_clus2_p0_05_down)
distinct_ReDisX_clus3DEG_p0_05_down=setdiff(distinct_ReDisX_clus3DEG_p0_05_down,degs_index_ReDisX_clus4_p0_05_down)
distinct_ReDisX_clus3DEG_p0_05_down=setdiff(distinct_ReDisX_clus3DEG_p0_05_down,degs_index_ReDisX_clus5_p0_05_down)
### distinct deg clus4
distinct_ReDisX_clus4DEG_p0_05_down=setdiff(degs_index_ReDisX_clus4_p0_05_down,degs_index_ReDisX_clus1_p0_05_down)
distinct_ReDisX_clus4DEG_p0_05_down=setdiff(distinct_ReDisX_clus4DEG_p0_05_down,degs_index_ReDisX_clus2_p0_05_down)
distinct_ReDisX_clus4DEG_p0_05_down=setdiff(distinct_ReDisX_clus4DEG_p0_05_down,degs_index_ReDisX_clus3_p0_05_down)
distinct_ReDisX_clus4DEG_p0_05_down=setdiff(distinct_ReDisX_clus4DEG_p0_05_down,degs_index_ReDisX_clus5_p0_05_down)
### distinct deg clus5
distinct_ReDisX_clus5DEG_p0_05_down=setdiff(degs_index_ReDisX_clus5_p0_05_down,degs_index_ReDisX_clus1_p0_05_down)
distinct_ReDisX_clus5DEG_p0_05_down=setdiff(distinct_ReDisX_clus5DEG_p0_05_down,degs_index_ReDisX_clus2_p0_05_down)
distinct_ReDisX_clus5DEG_p0_05_down=setdiff(distinct_ReDisX_clus5DEG_p0_05_down,degs_index_ReDisX_clus3_p0_05_down)
distinct_ReDisX_clus5DEG_p0_05_down=setdiff(distinct_ReDisX_clus5DEG_p0_05_down,degs_index_ReDisX_clus4_p0_05_down)



### save distinct DEG
write.csv(gene_SYMBOL[distinct_ReDisX_clus1DEG_p0_05_down,],file = 'result_v2/GSE59867_result/p_005_distinct_result/distict_gene_down_DisRedXclus1_GSE59867.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus2DEG_p0_05_down,],file = 'result_v2/GSE59867_result/p_005_distinct_result/distict_gene_down_DisRedXclus2_GSE59867.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus3DEG_p0_05_down,],file = 'result_v2/GSE59867_result/p_005_distinct_result/distict_gene_down_DisRedXclus3_GSE59867.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus4DEG_p0_05_down,],file = 'result_v2/GSE59867_result/p_005_distinct_result/distict_gene_down_DisRedXclus4_GSE59867.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus5DEG_p0_05_down,],file = 'result_v2/GSE59867_result/p_005_distinct_result/distict_gene_down_DisRedXclus5_GSE59867.csv')







distinct_ReDisX_all_p0_005=c(distinct_ReDisX_clus1DEG_p0_005[1:10],distinct_ReDisX_clus1DEG_p0_05_down[1:10],distinct_ReDisX_clus2DEG_p0_005[21:30],distinct_ReDisX_clus2DEG_p0_05_down[21:30],distinct_ReDisX_clus3DEG_p0_005[1:10],distinct_ReDisX_clus3DEG_p0_05_down[11:20],
                             distinct_ReDisX_clus4DEG_p0_005[21:30],distinct_ReDisX_clus4DEG_p0_05_down[21:30],distinct_ReDisX_clus5DEG_p0_005[1:10],distinct_ReDisX_clus5DEG_p0_05_down[1:10] )


degs_hub=gene_SYMBOL[distinct_ReDisX_all_p0_005,]


vst_hub = v[["E"]]

data_for_hub = vst_hub[degs_hub$SYMBOL,]
data_for_hub_label1=data_for_hub[,disease_target=="ReDisX_CAD_clus1"]
data_for_hub_label2=data_for_hub[,disease_target=="ReDisX_CAD_clus2"]
data_for_hub_label3=data_for_hub[,disease_target=="ReDisX_CAD_clus3"]
data_for_hub_label4=data_for_hub[,disease_target=="ReDisX_CAD_clus4"]
data_for_hub_label5=data_for_hub[,disease_target=="ReDisX_CAD_clus5"]
data_for_hub_label6=data_for_hub[,disease_target=="Control"]

data_for_hub_all =cbind(data_for_hub_label1,data_for_hub_label2,data_for_hub_label3,data_for_hub_label4,
                        data_for_hub_label5,data_for_hub_label6)

write.csv(data_for_hub_all, 'result_v2/GSE59867_result/p_005_distinct_result/data_for_hub_all_matrix.csv')


df_hub <- data.frame(x = data_for_hub_all.umap$layout[,1],
                 y = data_for_hub_all.umap$layout[,2],
                 Species = disease_target)

pdf(file="result_v2/GSE59867_result/p_005_distinct_result/UMAP_PLotReDisX_hub_GSE59867.pdf",width=5, height=5)
ggplot(df_hub, aes(x, y, colour = Species)) +
  geom_point()

dev.off()

my_breaks = c(seq(-2, -0.5, length.out=25),
              0,
              seq(0.5, 2, length.out=25))
my_colors = c("green", "black", "red")
my_colors = colorRampPalette(my_colors)(50)


annotation_col_hub = data.frame(
   Cluster = disease_target,disease_state=clustering_result$all_disease_type)
rownames(annotation_col_hub)=clustering_result$sampleid

annotation_row_hub = data.frame(DEGs_anno=rep(c("ReDisX_clus1_DEGs","ReDisX_clus2_DEGs","ReDisX_clus3_DEGs","ReDisX_clus4_DEGs","ReDisX_clus5_DEGs"),each=20))
rownames(annotation_row_hub)=degs_hub$SYMBOL

pdf(file="result_v2/GSE59867_result/p_005_distinct_result/ReDisX_hub_GSE59867.pdf")
pheatmap(data_for_hub_all, fontsize_row=2,fontsize_col = 2, scale='row',breaks = my_breaks, annotation_col=annotation_col_hub,annotation_row = annotation_row_hub, cluster_cols = F,cluster_rows = F,color = my_colors,col_split = annotation_col_hub$Cluster,
         gaps_col = cumsum(c(30,70,12,97,181))
)
dev.off()


############################ heatmap_all_in_one

distinct_ReDisX_all_p0_005_allineone=c(distinct_ReDisX_clus1DEG_p0_005[1:50],distinct_ReDisX_clus1DEG_p0_05_down[1:50],distinct_ReDisX_clus2DEG_p0_005[1:50],distinct_ReDisX_clus2DEG_p0_05_down[1:50],distinct_ReDisX_clus3DEG_p0_005[1:50],distinct_ReDisX_clus3DEG_p0_05_down[1:50],
                             distinct_ReDisX_clus4DEG_p0_005[1:50],distinct_ReDisX_clus4DEG_p0_05_down[1:50],distinct_ReDisX_clus5DEG_p0_005[1:11],distinct_ReDisX_clus5DEG_p0_05_down[1:25] )


degs_hub_allineone=gene_SYMBOL[distinct_ReDisX_all_p0_005_allineone,]


vst_hub_allineone = v[["E"]]

data_for_hub_allineone = vst_hub[degs_hub_allineone$SYMBOL,]
data_for_hub_label1_allineone=data_for_hub_allineone[,disease_target=="ReDisX_CAD_clus1"]
data_for_hub_label2_allineone=data_for_hub_allineone[,disease_target=="ReDisX_CAD_clus2"]
data_for_hub_label3_allineone=data_for_hub_allineone[,disease_target=="ReDisX_CAD_clus3"]
data_for_hub_label4_allineone=data_for_hub_allineone[,disease_target=="ReDisX_CAD_clus4"]
data_for_hub_label5_allineone=data_for_hub_allineone[,disease_target=="ReDisX_CAD_clus5"]
data_for_hub_label6_allineone=data_for_hub_allineone[,disease_target=="Control"]

data_for_hub_all_allineone =cbind(data_for_hub_label1_allineone,data_for_hub_label2_allineone,data_for_hub_label3_allineone,data_for_hub_label4_allineone,
                        data_for_hub_label5_allineone,data_for_hub_label6_allineone)



my_breaks = c(seq(-2, -0.5, length.out=25),
              0,
              seq(0.5, 2, length.out=25))
my_colors = c("green", "black", "red")
my_colors = colorRampPalette(my_colors)(50)


annotation_col_hub_allineone = data.frame(
  Cluster = disease_target,disease_state=clustering_result$all_disease_type)
rownames(annotation_col_hub_allineone)=clustering_result$sampleid

pdf(file="result_v2/GSE59867_result/p_005_distinct_result/ReDisX_hub_GSE59867_allineone.pdf")
pheatmap(data_for_hub_all_allineone, fontsize_row=2,fontsize_col = 2, scale='row',breaks = my_breaks, annotation_col=annotation_col_hub_allineone, cluster_cols = F,color = my_colors,col_split = annotation_col_hub_allineone$Cluster,
         gaps_col = cumsum(c(30,70,12,97,181))
)
dev.off()


