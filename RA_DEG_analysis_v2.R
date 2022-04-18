
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

clustering_result <- read_csv("matlab_GSE93272_ReDisX/RA_clustering_result_pluscontrol.csv")

edata=read.csv("matlab_GSE93272_ReDisX/RA_expression_final_pluscontrol.csv",header=FALSE)

gene_SYMBOL=read_csv('Top5000RA_genename.csv')

colnames(gene_SYMBOL)="SYMBOL"
colnames(clustering_result)[1]="sampleid"
colnames(clustering_result)[2]="final_clus"



clustering_result$final_clus=as.factor(clustering_result$final_clus)


clustering_result[,5]=clustering_result[,2]
colnames(clustering_result)[5]="cluster_for_limma"

clustering_result$cluster_for_limma<-gsub('1', 'ReDisX_clus1', clustering_result$cluster_for_limma)
clustering_result$cluster_for_limma<-gsub('2', 'ReDisX_clus2', clustering_result$cluster_for_limma)
clustering_result$cluster_for_limma<-gsub('3', 'ReDisX_clus3', clustering_result$cluster_for_limma)
clustering_result$cluster_for_limma<-gsub('4', 'Control', clustering_result$cluster_for_limma)


a=factor(clustering_result$cluster_for_limma)
disease_type_num=as.numeric(a)

row.names(edata)=gene_SYMBOL$SYMBOL
colnames(edata)=clustering_result$sampleid

y <- DGEList(edata)
y$genes <- gene_SYMBOL




#pdf(file="result/tradis_MDS.pdf")
### MDS plot
#a=plotMDS(edata)
#dev.off()




disease_target <- factor(clustering_result$cluster_for_limma, levels=c("ReDisX_clus1","ReDisX_clus2","ReDisX_clus3",'Control'))
design <- model.matrix(~ 0 + disease_target)





pdf(file="result_v2/GSE93272voom:Mean_variance_trend.pdf")
par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE,save.plot = TRUE)
dev.off()
#pdf(file="result_v2/tradis_boxplot.pdf")
##boxplot(v$E[1:10,1:10], xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
#abline(h=median(v$E),col="blue")
#dev.off()

write.csv(v[["E"]], 'RA_matrix_after_voom.csv')

fit <- lmFit(v)
names(fit)

cont.matrix <- makeContrasts(ReDisX_clus1_Control_DEG= disease_targetReDisX_clus1-disease_targetControl, ReDisX_clus2_Control_DEG=disease_targetReDisX_clus2-disease_targetControl
                             , ReDisX_clus3_Control_DEG=disease_targetReDisX_clus3-disease_targetControl,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)


toplist_RA=topTable(fit.cont,adjust.method="BH",p.value=0.05,number=20000)


summa.fit <-decideTests(fit.cont,method="separate",adjust.method="BH",p.value=0.05,lfc=0.05)
summary=data.frame(summa.fit)
summary(summa.fit)
write.csv(summary, 'ReDisX_DEG_result_RA.csv')
write.csv(toplist_RA, 'ReDisX_DEG_toplist_RA.csv')



####  DEG for RA, CAD, ReDisXclus2

degs_index_ReDisX_clus1= which(summary[,1] %in% c(1,-1))
degs_index_ReDisX_clus2= which(summary[,2] %in% c(1,-1))
degs_index_ReDisX_clus3= which(summary[,3] %in% c(1,-1))
#degs_index_ReDisX2= which(summary_ReDisX[,1] %in% c(1,-1))
#degs_index_ReDisX1= which(summary_ReDisX[,2] %in% c(1,-1))
#degs_index_ReDisX3= which(summary_ReDisX[,3] %in% c(1,-1))


distinct_ReDisX_clus1=setdiff(degs_index_ReDisX_clus1,degs_index_ReDisX_clus2)
distinct_ReDisX_clus1=setdiff(degs_index_ReDisX_clus1,degs_index_ReDisX_clus3)
#distinct_RADEG_CAD_DisRedX2=setdiff(distinct_RADEG_CAD,degs_index_ReDisX2)

distinct_ReDisX_clus2=setdiff(degs_index_ReDisX_clus2,degs_index_ReDisX_clus1)
distinct_ReDisX_clus2=setdiff(distinct_ReDisX_clus2,degs_index_ReDisX_clus3)
#distinct_CADDEG_RA_DisRedX2=setdiff(distinct_CADDEG_RA,degs_index_ReDisX2)
distinct_ReDisX_clus3=setdiff(degs_index_ReDisX_clus3,degs_index_ReDisX_clus2)
distinct_ReDisX_clus3=setdiff(distinct_ReDisX_clus3,degs_index_ReDisX_clus1)
#distinct_DisRedX2_RA=setdiff(degs_index_ReDisX2,degs_index_RA)
#distinct_DisRedX2_RA_CAD=setdiff(distinct_DisRedX2_RA,degs_index_CAD)
#distinct_DisRedX2_RA_CAD_clus1=setdiff(distinct_DisRedX2_RA_CAD,degs_index_ReDisX1)
#distinct_DisRedX2_RA_CAD_clus1_clus3=setdiff(distinct_DisRedX2_RA_CAD_clus1,degs_index_ReDisX3)

write.csv(gene_SYMBOL[distinct_ReDisX_clus1,],file = 'result_v2/GSE93272_result/distict_gene_DisRedXclus1_GSE93272.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus2,],file = 'result_v2/GSE93272_result/distict_gene_DisRedXclus2_GSE93272.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus3,],file = 'result_v2/GSE93272_result/distict_gene_DisRedXclus3_GSE93272.csv')


### set Venn color
myCol <- brewer.pal(3, "Pastel2")

### VENN DIAGRAM for three set
venn.diagram(
  x = list(degs_index_ReDisX_clus1, degs_index_ReDisX_clus2,  degs_index_ReDisX_clus3),
  category.names = c("ReDisX_clus1" , "ReDisX_clus2" , "ReDisX_clus3"),
  filename = 'result_v2/14_venn_diagramm_Three_Set_repeat.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 900 , 
  width = 900 , 
  resolution = 400,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)








### MD plot after DE

pdf(file="result_v2/ReDisX_GSE93272clus1_MD.pdf")


plotMD(fit.cont,coef=1,status=summa.fit[,"ReDisX_clus1_Control_DEG"], values = c(-1, 1), hl.col=c("blue","red"))
dev.off()
pdf(file="result_v2/ReDisX_GSE93272clus2_MD.pdf")
plotMD(fit.cont,coef=2,status=summa.fit[,"ReDisX_clus2_Control_DEG"], values = c(-1, 1), hl.col=c("blue","red"))
dev.off()
pdf(file="result_v2/ReDisX_GSE93272clus3_MD.pdf")
plotMD(fit.cont,coef=3,status=summa.fit[,"ReDisX_clus3_Control_DEG"], values = c(-1, 1), hl.col=c("blue","red"))
dev.off()
pdf(file="result_v2/ReDisX_GSE93272clus1_volcanoplot.pdf")
volcanoplot(fit.cont,coef=1,highlight=50,names=fit.cont$genes$SYMBOL, main="ReDisX_clus1_Control")
dev.off()
pdf(file="result_v2/ReDisX_GSE93272clus2_volcanoplot.pdf")
volcanoplot(fit.cont,coef=2,highlight=50,names=fit.cont$genes$SYMBOL, main="ReDisX_clus2_Control")
dev.off()
pdf(file="result_v2/ReDisX_GSE93272clus3_volcanoplot.pdf")
volcanoplot(fit.cont,coef=3,highlight=50,names=fit.cont$genes$SYMBOL, main="ReDisX_clus3_Control")

dev.off()

### heatmap  ReDisX_clus1



degs_clus1=gene_SYMBOL[distinct_ReDisX_clus1,]


vst_clus1 = v[["E"]]

data_for_clus1 = vst_clus1[degs_clus1$SYMBOL,]
data_for_clus1_label1=data_for_clus1[,disease_target=="ReDisX_clus1"]
data_for_clus1_label2=data_for_clus1[,disease_target=="ReDisX_clus2"]
data_for_clus1_label3=data_for_clus1[,disease_target=="ReDisX_clus3"]
data_for_clus1_label4=data_for_clus1[,disease_target=="Control"]
data_for_clus1_all =cbind(data_for_clus1_label1,data_for_clus1_label2,data_for_clus1_label3,data_for_clus1_label4)





my_breaks = c(seq(-2, -0.5, length.out=25),
              0,
              seq(0.5, 2, length.out=25))
my_colors = c("green", "black", "red")
my_colors = colorRampPalette(my_colors)(50)


annotation_col_clus1 = data.frame(
  Cluster = disease_target)
rownames(annotation_col_clus1)=clustering_result$sampleid


pdf(file="result_v2/ReDisX_clus1_GSE93272.pdf",width=30, height=30)
pheatmap(data_for_clus1_all, fontsize_row=4,fontsize_col = 2, scale='row',breaks = my_breaks, annotation_col=annotation_col_clus1, cluster_cols = F,color = my_colors,col_split = annotation_col_clus1$Cluster,
         gaps_col = cumsum(c(26,49,157))
         )
dev.off()





### heatmap  ReDisX_clus2



degs_clus2=gene_SYMBOL[distinct_ReDisX_clus2,]


vst_clus2 = v[["E"]]

data_for_clus2 = vst_clus2[degs_clus2$SYMBOL,]
data_for_clus2_label1=data_for_clus2[,disease_target=="ReDisX_clus1"]
data_for_clus2_label2=data_for_clus2[,disease_target=="ReDisX_clus2"]
data_for_clus2_label3=data_for_clus2[,disease_target=="ReDisX_clus3"]
data_for_clus2_label4=data_for_clus2[,disease_target=="Control"]
data_for_clus2_all =cbind(data_for_clus2_label1,data_for_clus2_label2,data_for_clus2_label3,data_for_clus2_label4)





my_breaks = c(seq(-2, -0.5, length.out=25),
              0,
              seq(0.5, 2, length.out=25))
my_colors = c("green", "black", "red")
my_colors = colorRampPalette(my_colors)(50)


annotation_col_clus2 = data.frame(
  Cluster = disease_target)
rownames(annotation_col_clus2)=clustering_result$sampleid


pdf(file="result_v2/ReDisX_clus2_GSE93272.pdf",width=30, height=30)
pheatmap(data_for_clus2_all, fontsize_row=4,fontsize_col = 2, scale='row',breaks = my_breaks, annotation_col=annotation_col_clus2, cluster_cols = F,color = my_colors,col_split = annotation_col_clus2$Cluster,
         gaps_col = cumsum(c(26,49,157))
)
dev.off()


### heatmap  ReDisX_clus3



degs_clus3=gene_SYMBOL[distinct_ReDisX_clus3,]


vst_clus3 = v[["E"]]

data_for_clus3 = vst_clus3[degs_clus3$SYMBOL,]
data_for_clus3_label1=data_for_clus3[,disease_target=="ReDisX_clus1"]
data_for_clus3_label2=data_for_clus3[,disease_target=="ReDisX_clus2"]
data_for_clus3_label3=data_for_clus3[,disease_target=="ReDisX_clus3"]
data_for_clus3_label4=data_for_clus3[,disease_target=="Control"]
data_for_clus3_all =cbind(data_for_clus3_label1,data_for_clus3_label2,data_for_clus3_label3,data_for_clus3_label4)





my_breaks = c(seq(-2, -0.5, length.out=25),
              0,
              seq(0.5, 2, length.out=25))
my_colors = c("green", "black", "red")
my_colors = colorRampPalette(my_colors)(50)


annotation_col_clus3 = data.frame(
  Cluster = disease_target)
rownames(annotation_col_clus3)=clustering_result$sampleid


pdf(file="result_v2/ReDisX_clus3_GSE93272.pdf",width=30, height=30)
pheatmap(data_for_clus3_all, fontsize_row=4,fontsize_col = 2, scale='row',breaks = my_breaks, annotation_col=annotation_col_clus3, cluster_cols = F,color = my_colors,col_split = annotation_col_clus3$Cluster,
         gaps_col = cumsum(c(26,49,157))
)
dev.off()



#### hubgenes in query heatmap


summa_p0_005.fit <-decideTests(fit.cont,method="separate",adjust.method="BH",p.value=0.005,lfc=0.07)
summary_p0_005=data.frame(summa_p0_005.fit)
summary(summa_p0_005.fit)
write.csv(summary_p0_005, 'result_v2/GSE93272_result/p_005_distinct_result/ReDisX_DEG_result_RA_p0_005.csv')




####  DEG for RA, CAD, ReDisXclus2

degs_index_ReDisX_clus1_p0_005= which(summary_p0_005[,1] %in% c(1))
degs_index_ReDisX_clus2_p0_005= which(summary_p0_005[,2] %in% c(1))
degs_index_ReDisX_clus3_p0_005= which(summary_p0_005[,3] %in% c(1))
#degs_index_ReDisX2= which(summary_ReDisX[,1] %in% c(1,-1))
#degs_index_ReDisX1= which(summary_ReDisX[,2] %in% c(1,-1))
#degs_index_ReDisX3= which(summary_ReDisX[,3] %in% c(1,-1))


distinct_ReDisX_clus1_p0_005=setdiff(degs_index_ReDisX_clus1_p0_005,degs_index_ReDisX_clus2_p0_005)
distinct_ReDisX_clus1_p0_005=setdiff(degs_index_ReDisX_clus1_p0_005,degs_index_ReDisX_clus3_p0_005)
#distinct_RADEG_CAD_DisRedX2=setdiff(distinct_RADEG_CAD,degs_index_ReDisX2)

distinct_ReDisX_clus2_p0_005=setdiff(degs_index_ReDisX_clus2_p0_005,degs_index_ReDisX_clus1_p0_005)
distinct_ReDisX_clus2_p0_005=setdiff(distinct_ReDisX_clus2_p0_005,degs_index_ReDisX_clus3_p0_005)
#distinct_CADDEG_RA_DisRedX2=setdiff(distinct_CADDEG_RA,degs_index_ReDisX2)
distinct_ReDisX_clus3_p0_005=setdiff(degs_index_ReDisX_clus3_p0_005,degs_index_ReDisX_clus2_p0_005)
distinct_ReDisX_clus3_p0_005=setdiff(distinct_ReDisX_clus3_p0_005,degs_index_ReDisX_clus1_p0_005)
#distinct_DisRedX2_RA=setdiff(degs_index_ReDisX2,degs_index_RA)
#distinct_DisRedX2_RA_CAD=setdiff(distinct_DisRedX2_RA,degs_index_CAD)
#distinct_DisRedX2_RA_CAD_clus1=setdiff(distinct_DisRedX2_RA_CAD,degs_index_ReDisX1)
#distinct_DisRedX2_RA_CAD_clus1_clus3=setdiff(distinct_DisRedX2_RA_CAD_clus1,degs_index_ReDisX3)

write.csv(gene_SYMBOL[distinct_ReDisX_clus1_p0_005,],file = 'result_v2/GSE93272_result/p_005_distinct_result/distict_gene_DisRedXclus1_GSE93272_p0_005.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus2_p0_005,],file = 'result_v2/GSE93272_result/p_005_distinct_result/distict_gene_DisRedXclus2_GSE93272_p0_005.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus3_p0_005,],file = 'result_v2/GSE93272_result/p_005_distinct_result/distict_gene_DisRedXclus3_GSE93272_p0_005.csv')

####  DEG for RA down regulate

degs_index_ReDisX_clus1_p0_005_down= which(summary_p0_005[,1] %in% c(-1))
degs_index_ReDisX_clus2_p0_005_down= which(summary_p0_005[,2] %in% c(-1))
degs_index_ReDisX_clus3_p0_005_down= which(summary_p0_005[,3] %in% c(-1))
#degs_index_ReDisX2= which(summary_ReDisX[,1] %in% c(1,-1))
#degs_index_ReDisX1= which(summary_ReDisX[,2] %in% c(1,-1))
#degs_index_ReDisX3= which(summary_ReDisX[,3] %in% c(1,-1))


distinct_ReDisX_clus1_p0_005_down=setdiff(degs_index_ReDisX_clus1_p0_005_down,degs_index_ReDisX_clus2_p0_005_down)
distinct_ReDisX_clus1_p0_005_down=setdiff(degs_index_ReDisX_clus1_p0_005_down,degs_index_ReDisX_clus3_p0_005_down)
#distinct_RADEG_CAD_DisRedX2=setdiff(distinct_RADEG_CAD,degs_index_ReDisX2)

distinct_ReDisX_clus2_p0_005_down=setdiff(degs_index_ReDisX_clus2_p0_005_down,degs_index_ReDisX_clus1_p0_005_down)
distinct_ReDisX_clus2_p0_005_down=setdiff(distinct_ReDisX_clus2_p0_005_down,degs_index_ReDisX_clus3_p0_005_down)
#distinct_CADDEG_RA_DisRedX2=setdiff(distinct_CADDEG_RA,degs_index_ReDisX2)
distinct_ReDisX_clus3_p0_005_down=setdiff(degs_index_ReDisX_clus3_p0_005_down,degs_index_ReDisX_clus2_p0_005_down)
distinct_ReDisX_clus3_p0_005_down=setdiff(distinct_ReDisX_clus3_p0_005_down,degs_index_ReDisX_clus1_p0_005_down)
#distinct_DisRedX2_RA=setdiff(degs_index_ReDisX2,degs_index_RA)
#distinct_DisRedX2_RA_CAD=setdiff(distinct_DisRedX2_RA,degs_index_CAD)
#distinct_DisRedX2_RA_CAD_clus1=setdiff(distinct_DisRedX2_RA_CAD,degs_index_ReDisX1)
#distinct_DisRedX2_RA_CAD_clus1_clus3=setdiff(distinct_DisRedX2_RA_CAD_clus1,degs_index_ReDisX3)

write.csv(gene_SYMBOL[distinct_ReDisX_clus1_p0_005_down,],file = 'result_v2/GSE93272_result/p_005_distinct_result/distict_gene_DisRedXclus1_GSE93272_p0_005_down.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus2_p0_005_down,],file = 'result_v2/GSE93272_result/p_005_distinct_result/distict_gene_DisRedXclus2_GSE93272_p0_005_down.csv')
write.csv(gene_SYMBOL[distinct_ReDisX_clus3_p0_005_down,],file = 'result_v2/GSE93272_result/p_005_distinct_result/distict_gene_DisRedXclus3_GSE93272_p0_005_down.csv')



distinct_ReDisX_all_p0_005=c(distinct_ReDisX_clus1_p0_005[1:10],distinct_ReDisX_clus1_p0_005_down[1:10],distinct_ReDisX_clus2_p0_005[11:20],distinct_ReDisX_clus2_p0_005_down[1:10],
                             distinct_ReDisX_clus3_p0_005[11:20],distinct_ReDisX_clus3_p0_005_down[21:30])


degs_hub=gene_SYMBOL[distinct_ReDisX_all_p0_005,]


vst_hub = v[["E"]]

data_for_hub = vst_hub[degs_hub$SYMBOL,]
data_for_hub_label1=data_for_hub[,disease_target=="ReDisX_clus1"]
data_for_hub_label2=data_for_hub[,disease_target=="ReDisX_clus2"]
data_for_hub_label3=data_for_hub[,disease_target=="ReDisX_clus3"]
data_for_hub_label4=data_for_hub[,disease_target=="Control"]
data_for_hub_all =cbind(data_for_hub_label1,data_for_hub_label2,data_for_hub_label3,data_for_hub_label4)





my_breaks = c(seq(-2, -0.5, length.out=25),
              0,
              seq(0.5, 2, length.out=25))
my_colors = c("green", "black", "red")
my_colors = colorRampPalette(my_colors)(50)


annotation_col_hub = data.frame(
  Cluster = disease_target,disease_state=clustering_result$all_disease_type)
rownames(annotation_col_hub)=clustering_result$sampleid

annotation_row_hub = data.frame(DEGs_partial=rep(c("ReDisX_clus1_DEGs","ReDisX_clus2_DEGs","ReDisX_clus3_DEGs"),each=20))
rownames(annotation_row_hub)=degs_hub$SYMBOL


pdf(file="result_v2/GSE93272_result/p_005_distinct_result/ReDisX_hub_GSE93272.pdf")
pheatmap(data_for_hub_all, fontsize_row=4,fontsize_col = 2, scale='row',breaks = my_breaks, annotation_col=annotation_col_hub,annotation_row=annotation_row_hub, cluster_cols = F, cluster_rows = F,color = my_colors,col_split = annotation_col_hub$Cluster,
         gaps_col = cumsum(c(26,49,157))
)
dev.off()

#################### whole big heatmap

distinct_ReDisX_allinone=c(distinct_ReDisX_clus1_p0_005,distinct_ReDisX_clus1_p0_005_down,distinct_ReDisX_clus2_p0_005,distinct_ReDisX_clus2_p0_005_down,
                             distinct_ReDisX_clus3_p0_005,distinct_ReDisX_clus3_p0_005_down)


degs_hub_allinone=gene_SYMBOL[distinct_ReDisX_allinone,]


vst_hub_allinone = v[["E"]]

data_for_hub_allinone = vst_hub[degs_hub_allinone$SYMBOL,]
data_for_hub_label1_allinone=data_for_hub_allinone[,disease_target=="ReDisX_clus1"]
data_for_hub_label2_allinone=data_for_hub_allinone[,disease_target=="ReDisX_clus2"]
data_for_hub_label3_allinone=data_for_hub_allinone[,disease_target=="ReDisX_clus3"]
data_for_hub_label4_allinone=data_for_hub_allinone[,disease_target=="Control"]
data_for_hub_all_allinone =cbind(data_for_hub_label1_allinone,data_for_hub_label2_allinone,data_for_hub_label3_allinone,data_for_hub_label4_allinone)





my_breaks = c(seq(-2, -0.5, length.out=25),
              0,
              seq(0.5, 2, length.out=25))
my_colors = c("green", "black", "red")
my_colors = colorRampPalette(my_colors)(50)


annotation_col_hub_allinone = data.frame(
  Cluster = disease_target,disease_state=clustering_result$all_disease_type)
rownames(annotation_col_hub_allinone)=clustering_result$sampleid




pdf(file="result_v2/GSE93272_result/p_005_distinct_result/ReDisX_hub_GSE93272_allinone.pdf")
pheatmap(data_for_hub_all_allinone, fontsize_row=2,fontsize_col = 2, scale='row',breaks = my_breaks, annotation_col=annotation_col_hub, cluster_cols = F,color = my_colors,col_split = annotation_col_hub$Cluster,
         gaps_col = cumsum(c(26,49,157))
)
dev.off()
