###
#8 writers (CBLL1, KIAA1429, METTL14,
#METTL3, RBM15, RBM15B, WTAP, and ZC3H13)
#erasers (ALKBH5 and FTO)
#13 readers (ELAVL1,FMR1, HNRNPA2B1, HNRNPC, IGF2BP1, IGF2BP2,IGF2BP3, LRPPRC, YTHDC1, YTHDC2, YTHDF1,YTHDF2, and YTHDF3).
#seruat.R for single cell of total cells
###for total cluster CMS cluster
library(data.table)
library(GEOquery)
library(dplyr)
library(Seurat)
library(ggplot2)
library(umap)
library(hypeR)
library(msigdbr)
library(tidyverse)
library(pheatmap)
library(Seurat)
library(clustree)
library(NMF)##
m6agenes<-c("CBLL1","KIAA1429","METTL14","METTL3", "RBM15","RBM15B", "WTAP",
            "ZC3H13","ELAVL1","FMR1", "HNRNPA2B1", "HNRNPC", "IGF2BP1", 
            "IGF2BP2", "IGF2BP3", "LRPPRC", "YTHDC1", "YTHDC2", "YTHDF1", 
            "YTHDF2","YTHDF3","ALKBH5","FTO")
##
mycol <- c("#223D6C","#D20A13","#FFD121","#088247",
           "#58CDD9","#5D90BA","#431A3D","#11AA4D","#91612D","#6E568C","#7A142C",
           "#E0367A","#D8D155","#64495D","#7CC767")
#################################
GSE132465_single<-fread("GSE132465_GEO_processed_CRC_10X_natural_log_TPM_matrix.txt")
GSE132465_single[1:5,1:5]
GSE132465_single<-data.frame(GSE132465_single)
rownames(GSE132465_single)<-GSE132465_single[,1]
GSE132465_single<-GSE132465_single[,-1]
# input phe data for single cell data
GSE132465_phe<-read.delim2("GSE132465_GEO_processed_CRC_10X_cell_annotation.txt")
table(GSE132465_phe$Cell_subtype)
table(GSE132465_phe$Cell_type)
GSE132465_phe$Index<-gsub("-",".",GSE132465_phe$Index)
rownames(GSE132465_phe)<-GSE132465_phe$Index
####  construct the seruat object
##
####################
#for GSE132465
#single_dat<-GSE132465_single
#single_dat[1:5,1:5]
#single_dat_sce <- CreateSeuratObject(counts =single_dat, project = "GSE144735_CAF_colon")
#dim(single_dat_sce) #results are 61 cells and 14494 features
#
###
#saveRDS(single_dat_sce, "single_dat_sce_GSE132465.rds")
#saveRDS(GSE132465_phe,"GSE132465_phe.rds")
#new_markers <- readRDS("../../scFT-paper_rds/20190213new_markers.rds")
## 过滤细胞 和基因  
#rownames(single_dat_sce)[grepl('^MT-',rownames(single_dat_sce))]
#rownames(single_dat_sce)[grepl('^RP[SL]',rownames(single_dat_sce))]
#single_dat_sce[["percent.mt"]] <- PercentageFeatureSet(single_dat_sce, pattern = "^MT-")
#fivenum(single_dat_sce[["percent.mt"]][,1])
#rb.genes <- rownames(single_dat_sce)[grep("^RP[SL]",rownames(single_dat_sce))]
#C<-GetAssayData(object = single_dat_sce, slot = "counts")
#percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
#single_dat_sce <- AddMetaData(single_dat_sce, percent.ribo, col.name = "percent.ribo")
#plot1 <- FeatureScatter(single_dat_sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(single_dat_sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1, plot2))
#VlnPlot(single_dat_sce, features = c("percent.ribo", "percent.mt"), ncol = 2)
#VlnPlot(single_dat_sce, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#VlnPlot(single_dat_sce, features = c("percent.ribo", "nCount_RNA"), ncol = 2)
#dim(single_dat_sce)
#single_dat_sce <- subset(single_dat_sce, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 10 & percent.ribo<20)

single_dat_sce <- readRDS('single_dat_sce_GSE132465.rds')
GSE132465_phe <- readRDS('GSE132465_phe.rds')
dim(single_dat_sce)
single_dat_sce <- FindVariableFeatures(single_dat_sce, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(single_dat_sce), 10)
top10
plot1 <- VariableFeaturePlot(single_dat_sce)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot2))
all.genes <- rownames(single_dat_sce)
single_dat_sce <- ScaleData(single_dat_sce, features = all.genes)
single_dat_sce <- RunPCA(single_dat_sce, features = VariableFeatures(object = single_dat_sce))
DimHeatmap(single_dat_sce, dims = 1:10, cells = 500, balanced = TRUE)
DimHeatmap(single_dat_sce, dims = 11:20, cells = 500, balanced = TRUE)
#DimHeatmap(single_dat_sce, dims = 21:30, cells = 500, balanced = TRUE)
#DimHeatmap(single_dat_sce, dims = 31:40, cells = 500, balanced = TRUE)
#DimHeatmap(single_dat_sce, dims = 41:50, cells = 500, balanced = TRUE)
#single_dat_sce <- JackStraw(single_dat_sce, num.replicate = 20)
#single_dat_sce <- ScoreJackStraw(single_dat_sce, dims = 1:20)
#JackStrawPlot(single_dat_sce, dims = 1:20)
ElbowPlot(single_dat_sce, ndims = 50, reduction = "pca")
single_dat_sce <- FindNeighbors(single_dat_sce, dims = 1:20)
##
single_dat_sce  <- FindClusters(
  object = single_dat_sce,
  resolution = c(seq(.1,1.6,.2))
)
clustree(single_dat_sce@meta.data, prefix = "RNA_snn_res.")
colnames(single_dat_sce@meta.data)

##########
single_dat_sce<- FindClusters(single_dat_sce, resolution = 0.1)
##################
single_dat_sce<- RunUMAP(single_dat_sce, dims = 1:20)
DimPlot(single_dat_sce, reduction = "umap", cols=mycol)

single_dat_sce<- RunTSNE(object = single_dat_sce, dims = 1:20, do.fast = TRUE)
DimPlot(single_dat_sce,reduction = "tsne",label=T)
###
phe<-GSE132465_phe
single_dat_sce=AddMetaData(object=single_dat_sce, metadata=phe, col.name=colnames(phe))
library(readxl)
newphe<-data.frame(read_excel("clinicalphe.xlsx",1, col_names= TRUE, col_types= NULL, na=" "))##chip data
colnames(phe)
colnames(newphe)
newphe<-merge(phe,newphe,by="Sample",all=T)
rownames(newphe)<-newphe$Index
single_dat_sce<-AddMetaData(single_dat_sce,metadata=newphe, col.name=colnames(newphe))
#saveRDS(single_dat_sce, "single_dat_sce_GSE132465_20210612.rds")
#saveRDS(GSE132465_phe,"GSE132465_phe.rds")

#
#BiocManager::install("dittoSeq")#
library(scater)
library(dittoSeq)
dittoDimPlot(single_dat_sce,reduction = "tsne", "Patient")
ggsave("Patient_GSE132465.pdf",height=4.66,width=6.83)
dittoDimPlot(single_dat_sce,reduction = "tsne","Cell_type",color.panel=mycol)
ggsave("Cell_type_TSNE_GSE132465.pdf",height=4.66,width=6.2)
dittoDimPlot(single_dat_sce,reduction = "tsne", "Class",color.panel=mycol)
ggsave("dittoDimPlot_Class_GSE132465.pdf",height=4.6,width=5.9)
dittoDimPlot(single_dat_sce, reduction = "tsne", "Sample")
ggsave("dittoDimPlot_sample_GSE132465.pdf",height=4.6,width=7.1)
dittoDimPlot(single_dat_sce, reduction = "tsne", "Cell_type",color.panel=mycol)
ggsave("dittoDimPlot_Cell_type_GSE132465.pdf",height=4.6,width=6.1)
dittoDimPlot(single_dat_sce, reduction = "tsne", "MSI",color.panel=mycol)
ggsave("dittoDimPlot_MSI_GSE132465.pdf",height=4.6,width=6.1)
dittoDimPlot(single_dat_sce, reduction = "tsne", "Gender",color.panel=mycol)
ggsave("dittoDimPlot_Gender_GSE132465.pdf",height=4.6,width=6.1)
dittoDimPlot(single_dat_sce, reduction = "tsne", "Stage",color.panel=mycol)
ggsave("dittoDimPlot_Stage_GSE132465.pdf",height=4.6,width=5.8)

single_dat_sce@meta.data$Agegroup<-ifelse(single_dat_sce@meta.data$Age>60,"Older","Young")
dittoDimPlot(single_dat_sce, reduction = "tsne", "Agegroup",color.panel=mycol)
ggsave("dittoDimPlot_Agegroup_GSE132465.pdf",height=4.6,width=5.8)

dittoBarPlot(single_dat_sce, "Cell_type", group.by = "Sample",color.panel=mycol)+coord_flip()
####
colnames(single_dat_sce@meta.data)
dittoDimPlot(single_dat_sce,reduction = "tsne", "Cell_type",split.by = "Class",color.panel=ann_cols)
ggsave("Cell_type_Class_dimplot_GSE132465.pdf",height=3.8,width=7.52)
###for immune vs no-immune
single_dat_sce@meta.data$cell_immune<-ifelse(single_dat_sce@meta.data$Cell_type %in% c("Epithelial cells","Stromal cells"), "Non-Immune cells","Immune cells")
dittoDimPlot(single_dat_sce,reduction = "tsne", "cell_immune",color.panel=ann_cols)
ggsave("Cell_type_Cell_immune_dimplot_GSE132465.pdf",height=4.6,width=6.2)
table(single_dat_sce@meta.data$Cell_type)
###
dittoDimPlot(single_dat_sce,reduction = "tsne", "ImmuneScore")
dittoDimPlot(single_dat_sce,reduction = "tsne", "TumorPurity",split.by = "Class")
dittoDimPlot(single_dat_sce,reduction = "tsne", "StromalScore")
#####
##
library(dittoSeq)
colnames(phe)
ann_cols<-c("#223D6C" ,"#D20A13", "#FFD121" ,"#088247", "#58CDD9","#5D90BA","#223D6C","#D20A13")
pdf("Heatmap_m6A_GSE132465.pdf",height=3.35, width=6.3)
#,annot.colors=ann_cols,
dittoHeatmap(subset(single_dat_sce, downsample = 10), m6agenes, use_raster = TRUE,
             annot.by = c("Sample","Cell_type","Class"),color.by=c("red","blue")) 
dev.off()
ggsave("Heatmap_m6A_GSE132465.pdf",height=3.35,width=6.3)
##
##seurat 包 画图
DimPlot(single_dat_sce, reduction = "umap",group.by = "Cell_type", cols=mycol)& NoLegend() & NoAxes()##总体无标签和坐标
DimPlot(single_dat_sce, reduction = "umap",group.by = "Cell_type",split.by="Class",label=T, cols=mycol,seed=123)
DimPlot(single_dat_sce,reduction = "tsne",group.by = "Cell_type",label=T, cols=mycol,seed=123)+theme(position="none")
ggsave("Tsne2_GSE132465.pdf",height=5.31,width=7.78)
DimPlot(single_dat_sce,reduction = "tsne",group.by = "Patient",label=T, seed=123)
DimPlot(single_dat_sce,reduction = "tsne",group.by = "Class",label=T, cols=mycol,seed=123)
ggsave("Tsne_tumor_normal_GSE132465.pdf",height=4.51,width=7.78)
DimPlot(single_dat_sce,reduction = "tsne",group.by = "Cell_type",split.by="Class",label=T, cols=mycol,seed=123)
ggsave("Tsne_tumor_normal_cell_GSE132465.pdf",height=4.31,width=8.78)
DimPlot(single_dat_sce,reduction = "pca",group.by = "Cell_type",label=T, cols=mycol)
ggsave("pca__GSE132465.pdf",height=5.31,width=7.78)
###################
GSE132465_tsne<-data.frame(single_dat_sce@reductions$tsne@cell.embeddings)
GSE132465_umap<-data.frame(single_dat_sce@reductions$umap@cell.embeddings)

###saving 保存数据
saveRDS(GSE132465_tsne, "single_dat_sce_GSE132465_tsne.rds")
saveRDS(GSE132465_umap, "single_dat_sce_GSE132465_umap.rds")
########
########
#??Dotplot
#set.seed(11000)
#reducedDim(GSE132465_colon_Epis_resolution_0.5, "force") <- igraph::layout_with_fr(g)
#plotReducedDim(single_dat_sce, colour_by="label", dimred="force")
##
#saveRDS(single_dat_sce, "single_dat_sce_GSE132465.rds")
#saveRDS(GSE132465_phe,"GSE132465_phe.rds")
###
####seurat 包运行
features<-m6agenes
FeaturePlot(single_dat_sce,reduction = "tsne", features = c("ALKBH5", "FTO"), ncol = 2)
writers<-c("CBLL1"    , "KIAA1429" , "METTL14", "METTL3", "RBM15" ,"RBM15B"  ,"WTAP"  , "ZC3H13")
FeaturePlot(single_dat_sce,reduction = "tsne", features = writers, ncol = 3)
ggsave("FeaturePlot__GSE132465.pdf",height=8.15,width=9)
readers<-m6agenes[c(9:21)]
FeaturePlot(single_dat_sce,reduction = "tsne", features = readers, ncol = 4)
ggsave("FeaturePlot_readers__GSE132465.pdf",height=8.15,width=12)
##散点图seurat 包运行
FeatureScatter(object = single_dat_sce, group.by = "Cell_type", 
               feature1 = 'NNMT', feature2 = 'CD3E', plot.cor = TRUE,
               cols=mycol)
######
colnames(single_dat_sce@meta.data)
FeaturePlot(single_dat_sce,reduction = "tsne", features = c("StromalScore", "ImmuneScore","ESTIMATEScore"), ncol = 3)

###saveRDS(single_cell_estimateScores,file="single_cell_estimateScores.rds")
single_cell_estimateScores<-readRDS("single_cell_estimateScores.rds")
single_cell_estimateScores$StromalScore<-as.numeric(single_cell_estimateScores$StromalScore)
single_cell_estimateScores$ImmuneScore<-as.numeric(single_cell_estimateScores$ImmuneScore)
single_cell_estimateScores$ESTIMATEScore<-as.numeric(single_cell_estimateScores$ESTIMATEScore)
single_cell_estimateScores$TumorPurity<-as.numeric(single_cell_estimateScores$TumorPurity)
rownames(single_cell_estimateScores)<-rownames(single_dat_sce@meta.data)
single_dat_sce=AddMetaData(object=single_dat_sce, metadata=single_cell_estimateScores, col.name=colnames(single_cell_estimateScores))
colnames(single_cell_estimateScores)
m6agenes
gene<-m6agenes[1]
for( gene in m6agenes){
dittoScatterPlot(
  object = single_dat_sce,
  x.var = gene, y.var = 'ImmuneScore',color.panel = mycol,
  color.var = "Cell_type") & NoLegend()
#ggsave(paste0("ImmuneScore",gene ,"dittoScatterplot.pdf"),width = 6.41, height=4.34)
ggsave(paste0("ImmuneScore",gene ,"dittoScatterplot.png"),width = 6.41, height=4.34)
dittoScatterPlot(
  object = single_dat_sce,
  x.var = gene, y.var = 'StromalScore',color.panel = mycol,
  color.var = "Cell_type")
#ggsave(paste0("StromalScore",gene ,"dittoScatterplot.pdf"),width = 6.41, height=4.34)
ggsave(paste0("StromalScore",gene ,"dittoScatterplot.png"),width = 6.41, height=4.34)
dittoScatterPlot(
  object = single_dat_sce,
  x.var = gene, y.var = 'TumorPurity',color.panel = mycol,
  color.var = "Cell_type")
#ggsave(paste0("TumorPurity",gene ,"dittoScatterplot.pdf"),width = 6.41,height=4.34)
ggsave(paste0("TumorPurity",gene ,"dittoScatterplot.png"),width = 6.41,height=4.34)
}
##
for( gene in m6agenes){
  dittoScatterPlot(
    object = single_dat_sce,
    x.var = gene, y.var = 'ImmuneScore',color.panel = mycol,
    color.var = "Cell_type") & NoLegend()
  #ggsave(paste0("ImmuneScore",gene ,"dittoScatterplot.pdf"),width = 6.41, height=4.34)
  ggsave(paste0("ImmuneScore",gene ,"dittoScatterplot.png"),width = 4.8, height=4.34)
  dittoScatterPlot(
    object = single_dat_sce,
    x.var = gene, y.var = 'StromalScore',color.panel = mycol,
    color.var = "Cell_type") & NoLegend()
  #ggsave(paste0("StromalScore",gene ,"dittoScatterplot.pdf"),width = 6.41, height=4.34)
  ggsave(paste0("StromalScore",gene ,"dittoScatterplot.png"),width = 4.8, height=4.34)
  dittoScatterPlot(
    object = single_dat_sce,
    x.var = gene, y.var = 'TumorPurity',color.panel = mycol,
    color.var = "Cell_type") & NoLegend()
  #ggsave(paste0("TumorPurity",gene ,"dittoScatterplot.pdf"),width = 6.41,height=4.34)
  ggsave(paste0("TumorPurity",gene ,"dittoScatterplot.png"),width = 4.8, height=4.34)
}

###
###计算分数合并后 
#metascore_writers=colSums(single_dat_sce@assays$RNA@scale.data[rownames(single_dat_sce@assays$RNA@scale.data) %in% writers,])
#metascore_erasers=colSums(single_dat_sce@assays$RNA@scale.data[rownames(single_dat_sce@assays$RNA@scale.data) %in% erasers,])
#metascore_readers=colSums(single_dat_sce@assays$RNA@scale.data[rownames(single_dat_sce@assays$RNA@scale.data) %in% readers,])
#methy_metascore<-data.frame(metascore_writers,metascore_erasers,metascore_readers)
str(methy_metascore)
writers<-list(c("CBLL1"    , "KIAA1429" , "METTL14"  , 
           "METTL3", "RBM15" ,    "RBM15B"  ,  
           "WTAP"  ,    "ZC3H13"))
erasers<-list(c("ALKBH5","FTO"))
readers<-list(m6agenes[9:21])
totalm6agenes<-list(m6agenes)
#single_dat_sce=AddMetaData(object=single_dat_sce, metadata=methy_metascore, col.name=colnames(methy_metascore))
single_dat_sce=AddModuleScore(object=single_dat_sce, features=writers,name="Writers")
single_dat_sce=AddModuleScore(object=single_dat_sce, features=erasers,name="Erasers")
single_dat_sce=AddModuleScore(object=single_dat_sce, features=readers,name="Readers")
single_dat_sce=AddModuleScore(object=single_dat_sce, features=totalm6agenes,name="Total_m6A")
#single_dat_sce@meta.data[,c(16:38)]<-c()
colnames(single_dat_sce@meta.data)
colnames(single_dat_sce@meta.data)[c(16:19)]<-c("WritersScore", "ErasersScore","ReadersScore","Total_m6AScore")

methy<-colnames(methy_metascore)
methy<-c("WritersScore", "ErasersScore","ReadersScore","Total_m6AScore")
gene<-methy[4]
for( gene in methy){
  dittoScatterPlot(
    object = single_dat_sce,
    x.var = gene, y.var = 'ImmuneScore',color.panel = mycol,
    color.var = "Cell_type")
  ggsave(paste0("ImmuneScore",gene ,"dittoScatterplot.pdf"),width = 6.41, height=4.34)
  ggsave(paste0("ImmuneScore",gene ,"dittoScatterplot.png"),width = 6.41, height=4.34)
  dittoScatterPlot(
    object = single_dat_sce,
    x.var = gene, y.var = 'StromalScore',color.panel = mycol,
    color.var = "Cell_type")
  ggsave(paste0("StromalScore",gene ,"dittoScatterplot.pdf"),width = 6.41, height=4.34)
  ggsave(paste0("StromalScore",gene ,"dittoScatterplot.png"),width = 6.41, height=4.34)
  dittoScatterPlot(
    object = single_dat_sce,
    x.var = gene, y.var = 'TumorPurity',color.panel = mycol,
    color.var = "Cell_type")
  ggsave(paste0("TumorPurity",gene ,"dittoScatterplot.pdf"),width = 6.41,height=4.34)
  ggsave(paste0("TumorPurity",gene ,"dittoScatterplot.png"),width = 6.41,height=4.34)
}
#### for immune checkpoints 
## ,CTLA-4 CD274,PDCD1,PDCD1LG2
dittoScatterPlot(
  object = single_dat_sce,
  x.var = "WTAP", y.var = 'CD274',color.panel = mycol,
  color.var = "Class")
ggsave(paste0("TumorPurity",gene ,"dittoScatterplot.pdf"),width = 6.41,height=4.34)
ggsave(paste0("TumorPurity",gene ,"dittoScatterplot.png"),width = 6.41,height=4.34)

?VlnPlot
features<-m6agenes
RidgePlot(single_dat_sce, features = c("TP53", "NNMT"), ncol = 3)

for( gene in features){
  VlnPlot(object = single_dat_sce, features = gene,split.by= "Class",group.by = "Cell_type", ncol = 1,pt.size = 0,cols = mycol) & xlab("") & NoLegend()
  ggsave(paste0("vlnplot",gene,".pdf"),height=3.48,width = 3.80)
  ggsave(paste0("vlnplot",gene,".png"),height=3.48,width = 3.80)
}
VlnPlot(object = single_dat_sce, features = c("WTAP"),split.by= "Class",group.by = "Cell_type", ncol = 1,pt.size = 0,cols = mycol) & xlab("")
VlnPlot(single_dat_sce, features = features[1:8],split.by= "Class",group.by = "Cell_type", ncol = 4,pt.size = 0)
VlnPlot(single_dat_sce, features = features[9:16],split.by= "Class",group.by = "Cell_type", ncol = 4,pt.size = 0) & xlab("")
VlnPlot(single_dat_sce, features = features[17:23],split.by= "Class",group.by = "Cell_type", ncol = 4,pt.size = 0) & xlab("")

dittoDimHex(single_dat_sce,reduction = "tsne")
dittoPlot(single_dat_sce, "WTAP", group.by = "Cell_type",
          plots = c("jitter", "vlnplot", "boxplot"), # <- order matters
          # change the color and size of jitter points
          jitter.color = "blue", jitter.size = 0.5,
          # change the outline color and width, and remove the fill of boxplots
          boxplot.color = "white", boxplot.width = 0.1,
          boxplot.fill = FALSE
          # change how the violin plot widths are normalized across groups
         # vlnplot.scaling = "count"
)

#plotExpression(single_dat_sce, features = c("TP53", "CD79A"), x = "source", ncol = 2, colour_by = "Cell_type")
dittoFeaturePlot(single_dat_sce, features = c("NNMT", "CD79A"), ncol = 3)

dittoDimPlot(single_dat_sce, "NNMT")
dittoHeatmap(subset(single_dat_sce,downsample=50) ,m6agenes, scaled.to.max = TRUE,
             complex = TRUE,
             use_raster = TRUE)


DotPlot(single_dat_sce, features = unique(features),group.by =  "Class") & RotatedAxis() & xlab("") 
DotPlot(single_dat_sce, features = unique(features),group.by =  "Cell_type") + coord_flip()
ggsave("DotPlot_m6Agenes.pdf",height=2.8,width = 7.66)
ggsave("DotPlot_m6Agenes_Cell_type2.pdf",height=5.2,width = 9.56)
DotPlot(single_dat_sce, features = unique(features),group.by="Patient") + coord_flip()


DotPlot(single_dat_sce, features = unique(features),size = 3, cols = c("lightgrey", "red"),group.by = "Cell_type") + RotatedAxis()+coord_flip()
ggsave("dotplotofm6Agenes.pdf",height=4.1,width=9.5)
DoHeatmap(subset(single_dat_sce, downsample = 50), group.by = "Cell_type",features = features, size = 3)
ggsave("DoHeatmapofm6Agenes.pdf",height=3.36,width=8.05)
################ averageexpression
library(scales)
library(ggplot2)
library(ggpubr)
library(ggplotify)
library(pheatmap)
AveExpression <- AverageExpression(single_dat_sce, assays = "RNA",group.by = "Cell_type", split.by="Class",features = features,verbose = TRUE) %>% .$RNA
AveExpression<-data.frame(AveExpression)
AveExpression$Gene <- rownames(AveExpression)
Ave_df <- reshape2::melt(AveExpression,id.vars= "Gene")
colnames(Ave_df) <- c("Gene", "Cluster", "Expression")
###
pB1 <- ggplot(Ave_df,aes(x=Gene, y=Expression))+
  geom_hline(yintercept = seq(0, 10, 2.5),linetype = 2, color = "lightgray",size=1)+
  geom_line()+
  geom_segment(aes(x=Gene,xend=Gene,y=0,yend=Expression),color="lightgray",size = 1.5)+
  geom_point(size=3,aes(color=Gene))+
  #scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07", "#41ab5d")) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Gene exp.")
pB1 <- facet(pB1, facet.by = "Cluster",ncol = length(unique(Ave_df$Cluster)),panel.labs.font = list(size = 12),panel.labs.background = list(fill = "#a6cee3"))
pB1 <- pB1 + scale_y_continuous(position = "left")+  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_blank())+ ## 删去所有刻度标签# theme(axis.text.y = element_blank())   ## 设置 axis.text.y 则只删去 Y 轴的刻度标签，X 轴同理。
  guides(fill=guide_legend(ncol=4))+
  #ylim(0,8)+
  theme(legend.position = "top")
       # panel.border = element_blank(),## 去掉最外层的正方形边框 
       # axis.ticks.x = element_line(color =  NA))
pB1
ggsave("m6Agene_averageexpression.pdf",height=6,width=10.51)
#####paitents with cluster
table(single_dat_sce@meta.data$Cell_type)
single_dat_sce@meta.data$Patients<- str_split(colnames(single_dat_sce), "_") %>% sapply(., "[[", 1) 
pB2_df <- table(single_dat_sce@meta.data$Cell_type,single_dat_sce@meta.data$Sample) %>% reshape2::melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")
clusters<-c( 'Epithelial cells', 'Stromal cells', 'T cells','B cells','Mast cells','Myeloids')
pB2_df$Cluster <- factor(pB2_df$Cluster,levels = clusters)

#sample_color <- c("#d95f02","#66a61e","#1b9e77","#e7298a","#386cb0")
pB2 <- ggplot(data = pB2_df, aes(x =Sample , y = Number, fill =Cluster )) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cluster_cell_Ratio")+
  scale_y_continuous(position = "right",labels = percent)+  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12,  hjust = 0,angle=0, colour = "black"))+ coord_flip()##翻转过来
pB2
ggsave("patients with cell barplot.pdf")## 5.96 x 7.26
####
saveRDS(estimateScores,file="estimateScores.rds")
saveRDS(single_cell_estimateScores,file="single_cell_estimateScores.rds")
##
single_dat_sce.markers <- FindAllMarkers(single_dat_sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
single_dat_sce.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 2, wt = avg_log2FC)
single_dat_sce.markers_top20 <- single_dat_sce.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)

DoHeatmap(single_dat_sce, features = single_dat_sce.markers_top20$gene, label = F) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white")
feature <- single_dat_sce.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 2, wt = avg_log2FC)
feature<-data.frame(feature)
feature<-feature$gene
colon2_dot_plot <- DotPlot(single_dat_sce, features = feature, cols = c("lightgrey", "red")) + RotatedAxis()
colon2_dot_plot

message("top contributing genes (by percentage) contributing to signature")
##NMF
#new_markers <- readRDS("../../scFT-paper_rds/20190213new_markers.rds")
#ftsc.integrated$KRT17_score <- colSums(ftsc.integrated@assays$RNA@data[new_markers$C4,])
#saveRDS(ftsc.integrated, "rds/ftsc.integrated_20191013.rds")
#
#

library(NMF)
table(NMFgroup)
#########
single_dat_sce=AddMetaData(object=single_dat_sce, metadata=NMFgroup, col.name="NMFgroup")
DimPlot(single_dat_sce, reduction = "tsne", group.by="NMFgroup",cols=mycol)
DimPlot(single_dat_sce, reduction = "tsne",cols=mycol)
table(single_dat_sce$NMFgroup,single_dat_sce$seurat_clusters)
###################################################

####dot plot of each features
feature <- m1agenes
dot_plot <- DotPlot(single_dat_sce, features = feature, group.by='NMFgroup', cols = c("lightgrey", "red")) + RotatedAxis()
dot_plot 
###
##
