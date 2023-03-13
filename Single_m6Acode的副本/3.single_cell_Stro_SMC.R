#
library(stringr)
library(magrittr)
library(readxl)
library(forcats)
library(NMF)
library(data.table)
library(ggsci)
library(pheatmap)
library(AUCell)
library(GEOquery)
library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(corrplot)#
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
###for stromal cell m6A regulators
#####previous step  准备细胞轨 和 NMF差异基因
#gmt=read.gmt(misfile) # 读入 gmt 文件   input from the MsiDB data.    
m6agenes<-c("CBLL1","KIAA1429","METTL14","METTL3", "RBM15","RBM15B", "WTAP",
            "ZC3H13","ELAVL1","FMR1", "HNRNPA2B1", "HNRNPC", "IGF2BP1", 
            "IGF2BP2", "IGF2BP3", "LRPPRC", "YTHDC1", "YTHDC2", "YTHDF1", 
            "YTHDF2","YTHDF3","ALKBH5","FTO")
mycol <- c("#223D6C","#D20A13","#FFD121","#088247",
           "#58CDD9","#5D90BA","#431A3D","#11AA4D","#91612D","#6E568C","#7A142C",
           "#E0367A","#D8D155","#64495D","#7CC767")
#load("Pan_epi_cells_CRC.Rdata")
###1
single_dat<-readRDS("Stromal cellssingle_dat_GSE132465.rds")
rowGenenames<-readRDS("rowGenenames.rds")
rownames(single_dat)<-rowGenenames
single_dat[1:5,1:5]
GSE132465_phe <- readRDS('GSE132465_phe.rds')#####
table(GSE132465_phe$Cell_subtype)
###################################
patients<-str_split(colnames(single_dat), "_") %>% sapply(., "[[", 1) %>% unique()
patients
single_dat_Tumor<-single_dat[,str_split(colnames(single_dat), "_") %>% sapply(., "[[", 1) %in% patients[1:23]]
#single_dat_Tumor
########################
sce <- CreateSeuratObject(counts =single_dat, project = "GSE132465_Epi_colon")
dim(sce)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(sce), 10)
top10
plot1 <- VariableFeaturePlot(sce)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot2))
all.genes <- rownames(sce)
sce <- ScaleData(sce, features = all.genes)
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
DimHeatmap(sce, dims = 1:10, cells = 500, balanced = TRUE)
DimHeatmap(sce, dims = 11:20, cells = 500, balanced = TRUE)
#DimHeatmap(sce, dims = 21:30, cells = 500, balanced = TRUE)
#DimHeatmap(sce, dims = 31:40, cells = 500, balanced = TRUE)
#DimHeatmap(sce, dims = 41:50, cells = 500, balanced = TRUE)
#sce <- JackStraw(sce, num.replicate = 20)
#sce <- ScoreJackStraw(sce, dims = 1:20)
#JackStrawPlot(sce, dims = 1:20)
ElbowPlot(sce, ndims = 50, reduction = "pca")
sce <- FindNeighbors(sce, dims = 1:20)
##
sce  <- FindClusters(
  object = sce,
  resolution = c(seq(.1,1.6,.2))
)
clustree(sce@meta.data, prefix = "RNA_snn_res.")
colnames(sce@meta.data)
##########
sce<- FindClusters(sce, resolution = 0.1)
##################RunUMAP
sce<- RunUMAP(sce, dims = 1:20)
#RunTSNE
sce<- RunTSNE(object = sce, dims = 1:20, do.fast = TRUE)
DimPlot(sce, reduction = "umap",label=T, cols=mycol)
DimPlot(sce,reduction = "tsne",label=T,cols=mycol)
DimPlot(sce,reduction = "tsne",group.by="Sample")
DimPlot(sce,reduction = "tsne",group.by="Cell_subtype",label=T,cols=mycol)
DimPlot(sce,reduction = "tsne",cols=mycol)
##
phe<-GSE132465_phe
sce=AddMetaData(object=sce, metadata=phe, col.name=colnames(phe))
newphe<-data.frame(read_excel("clinicalphe.xlsx",1, col_names= TRUE, col_types= NULL, na=" "))##chip data
newphe<-merge(phe,newphe,by="Sample",all=T)
rownames(newphe)<-newphe$Index
sce<-AddMetaData(sce,metadata=newphe, col.name=colnames(newphe))###
####
#
library(dittoSeq)
sce@meta.data$Fibroblast<-ifelse(sce@meta.data$Cell_subtype %in% c("Stromal 1",
                                                                   "Stromal 2","Stromal 3","Myofibroblasts"),
                                                                   "Fibroblast","Non-Fibroblast"  )
sce@meta.data$Fibroblast2<-ifelse(sce@meta.data$Cell_subtype %in% c("Stromal 1","Stromal 2","Stromal 3"),
                                  "LipoFibroblast",ifelse(sce@meta.data$Cell_subtype=="Myofibroblasts",
                                                          "Myofibroblasts","Non-Fibroblast"))
#dittoDimPlot(sce,reduction = "tsne","Cell_subtype",color.panel = mycol)##
#ggsave("Cell_subtype_strimal_GSE132465.pdf",height=4.6,width=6.3)###
#dittoDimPlot(sce,reduction = "tsne","Fibroblast",color.panel = mycol)##
#ggsave("FibroblastvsNonFibroblast_GSE132465.pdf",height=4.6,width=6.1)###
dittoDimPlot(sce,reduction = "tsne","Fibroblast",split.by="Class",color.panel = mycol)##
ggsave("FibroblastvsNonFibroblast_Class_GSE132465.pdf",height=2.6,width=5.86)###
##############################################
#DotPlot(sce, features =m6agenes,group.by="Fibroblast2", cols = c("green", "red")) + RotatedAxis()+ coord_flip()
#ggsave("Dotplot_m6A_Fibroblasts_GSE132465.pdf",height=4.9,width=4.62)###
#DotPlot(sce, features =m6agenes,group.by="Fibroblast", cols = c("green", "red")) + RotatedAxis()+ coord_flip()
#ggsave("Dotplot_m6A_FibroblastvsNonFibroblast_GSE132465.pdf",height=4.9,width=4.22)###
##check thenames of meta.data##########
colnames(sce@meta.data)
##############################################
DotPlot(sce, features =m6agenes,group.by="Sample", cols = c("green", "red")) + RotatedAxis()
DotPlot(sce, features =m6agenes,group.by="Cell_subtype", cols = c("green", "red")) + RotatedAxis()+ coord_flip()
ggsave("dotplot_stromal_SMC.pdf")
#注释细胞类型###################################
#cg=BlueprintEncodeData()
#cg=DatabaseImmuneCellExpressionData()
#cg=NovershternHematopoieticData()
#cg=MonacoImmuneData()
#cg=ImmGenData()
#cg=MouseRNAseqData()
#cg=HumanPrimaryCellAtlasData()
##
#counts<-sce@assays$RNA@counts
#furhter analysis for stromal cells CAF firstly. ###
#############################################
######
#######################
single_dat_stromal_names<-c("Stromal 1","Stromal 2","Stromal 3","Myofibroblasts")
single_dat_stromal<-single_dat[,colnames(single_dat) %in% phe$Index[phe$Cell_subtype %in% single_dat_stromal_names]]
#single_dat_stromal<-single_dat_stromal[,colnames(single_dat_stromal) %in% phe$Index[phe$Class=="Tumor"]]
single_dat_stromal[1:5,1:5]
#################
sce <- CreateSeuratObject(counts =single_dat_stromal, project = "GSE132465_stromal_colon")
dim(sce)
## 过滤细胞 和基因  
#rownames(sce)[grepl('^MT-',rownames(sce))]
#rownames(sce)[grepl('^RP[SL]',rownames(sce))]
#sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
#fivenum(sce[["percent.mt"]][,1])
#rb.genes <- rownames(sce)[grep("^RP[SL]",rownames(sce))]
#C<-GetAssayData(object = sce, slot = "counts")
#percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
#sce <- AddMetaData(sce, percent.ribo, col.name = "percent.ribo")
#plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1, plot2))
##VlnPlot(sce, features = c("percent.ribo", "percent.mt"), ncol = 2)
##VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#VlnPlot(sce, features = c("percent.ribo", "nCount_RNA"), ncol = 2)
#dim(sce)
#sce <- subset(sce, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 20 & percent.ribo< 20)
###########
patients<-str_split(colnames(single_dat), "_") %>% sapply(., "[[", 1) %>% unique()
patients
sce <- CreateSeuratObject(counts =single_dat_stromal, project = "GSE132465_stromal_colon")
dim(sce)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
#top10 <- head(VariableFeatures(sce), 10)
#top10
#plot1 <- VariableFeaturePlot(sce)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot2))
all.genes <- rownames(sce)
sce <- ScaleData(sce, features = all.genes)
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
DimHeatmap(sce, dims = 1:10, cells = 500, balanced = TRUE)
DimHeatmap(sce, dims = 11:20, cells = 500, balanced = TRUE)
ElbowPlot(sce, ndims = 50, reduction = "pca")
sce <- FindNeighbors(sce, dims = 1:20)
##########
sce<- FindClusters(sce, resolution = 0.1)
#RunUMAP
sce<- RunUMAP(sce, dims = 1:20)
DimPlot(sce, reduction = "umap",label=T, cols=mycol)
#RunTSNE
sce<- RunTSNE(object = sce, dims = 1:20, do.fast = TRUE)
DimPlot(sce,reduction = "tsne",label=T,cols=mycol)

####### for clinical data
## re annotation for stromal cells 
newphe<-data.frame(read_excel("clinicalphe.xlsx",1, col_names= TRUE, col_types= NULL, na=" "))##chip data
newphe<-merge(phe,newphe,by="Sample",all=T)
rownames(newphe)<-newphe$Index
sce<-AddMetaData(sce,metadata=newphe, col.name=colnames(newphe))###
##########
table(sce@meta.data$Cell_subtype)
sce@meta.data$Fibroblast2<-ifelse(sce@meta.data$Cell_subtype %in% c("Stromal 1","Stromal 2","Stromal 3"),"LipoFibroblast", 
                                             ifelse(sce@meta.data$Cell_subtype=="Myofibroblasts","Myofibroblasts","Non-Fibroblast"))
#VlnPlot(sce,m6agenes,group.by="Cell_subtype",pt.size = 0,ncol=1)
DotPlot(sce, features =m6agenes,group.by="Sample", cols = c("green", "red")) + RotatedAxis()
DotPlot(sce, features =m6agenes,group.by="Fibroblast2", cols = c("green", "red")) + RotatedAxis()+ coord_flip()
###
######################find markers in the TUmor or Normal Fibroblasts
sce@meta.data$Class_type<-ifelse(sce@meta.data$Class =="Tumor",1,0)
DotPlot(sce, features =m6agenes,group.by="Class_type", cols = c("green", "red")) + RotatedAxis()+ coord_flip()    
ggsave("dotpolot_m6A_class_CAF.pdf")
####
sce$Class_type<-as.factor(sce$Class_type)
Idents(sce)<-sce$Class_type
sce.markers <- FindAllMarkers(sce, group.by="Class_type",
                                         only.pos = TRUE, 
                                         min.pct = 0.25, logfc.threshold = 0.25)
write.csv(sce.markers,"sce.markers_SMC_CAF_tumor_normal.csv")###
sce.markers  %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
sce.markers_top5 <- sce.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)
###
row_metadata <- as.data.frame(sce.markers_top5[,c("cluster","gene")])
row_metadata$celltype <- "celltype"
row_metadata$gene <- factor(row_metadata$gene,levels = row_metadata$gene)
gene <- as.character(row_metadata$gene)
label <- c("Normal Fibroblas",NA,gene[1:4],
           "Primary Fibroblast ",NA,gene[5:8]
)
p1<- ggplot(row_metadata,mapping = aes(x=rev(gene),y=celltype,fill=cluster))+
  geom_tile()+
  scale_fill_manual(values =mycol)+###ggsci::pal_d3()(9)
  theme(legend.position = "none",
        axis.title = element_blank(),axis.text = element_blank(),
        axis.ticks = element_blank(),panel.background = element_blank())+
  coord_flip()+
  annotate("text",
           x=rev(rep(seq(1.3,8,4/3),each=2)),
           y=rep(c(0.8,1.3),times=6),
           label=label,
           size=rep(c(3.5,NA,3,3,3,3),2),
           fontface=rep(c("bold",NA,"italic","italic","italic","italic"),2))
p1##
ggsave("annotate_DoHeatmap_Class_type_top10.pdf",height =3.17,width=3.11)
#p %>%aplot::insert_left(p1,0.3)
#logFCfilter<-0
#adjPvalFilter<-0.01
#sig.markers=sce.markers[(abs(as.numeric(as.vector(sce.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(sce.markers$p_val_adj))<adjPvalFilter),]
p<-DoHeatmap(sce,features = sce.markers_top5$gene) & NoLegend()
ggsave("DoHeatmap_Class_type_top10.pdf",height =3.17,width=9.31)
#################
#先鉴定 CAF  tumor associated Fibroblasts CAFs then cluster for CAF
single_dat_Tumor<-single_dat[,str_split(colnames(single_dat), "_") %>% sapply(., "[[", 1) %in% patients[1:23]]
single_dat_stromal_tumor<-single_dat_stromal[,colnames(single_dat_stromal) %in% phe$Index[phe$Class=="Tumor"]]
#NMF for single_dat_stromal
#nmf.input<-single_dat_myoCAF[m6agenes,]

nmf.input<-single_dat_stromal[m6agenes,]
nmf.input<-nmf.input[rowSums(nmf.input!= 0)>=1,colSums(nmf.input != 0) >= 1]
############
ranks <- 2:5
estim <- lapply(ranks, function(r){
  fit <- nmf(nmf.input, r, nrun = 5, seed = 4, method = "lee") #nrun
  list(fit = fit, consensus = consensus(fit), .opt = "vp",coph = cophcor(fit))
})
names(estim) <- paste('rank', ranks)

#pheneticϵ
#pdf("Cophenetic coefficient for seleting optimal nmf rank.pdf")
par(cex.axis=1.5)
plot(ranks, sapply(estim, '[[', 'coph'), xlab="", ylab="", type="b", col="red", lwd=4,xaxt="n")
axis(side = 1, at=1:5)
title(xlab="number of clusters", ylab="Cophenetic coefficient", cex.lab=1.5)
#invisible(dev.off())

######
rank <- 3
seed <- 2019620
#rownames(nmf.input) <- gsub("Signature","Sig",rownames(nmf.input)) #
mut.nmf <- nmf(nmf.input, 
               rank = rank, 
               seed = seed,  
               #method="nsNMF",
               method = "lee") 
group_CAF<- predict(mut.nmf) # 
save(group_CAF,file="group_CAF_SMC_NMF.Rdata")
save(group_CAF,file="group_CAF_SMC_NMFallcaf.Rdata")
load("group_CAF_SMC_NMF.Rdata")
load("group_CAF_SMC_NMFallcaf.Rdata")
table(group_CAF)
library(NMF)
#jco <- c("#2874C5","#EABF00","#C6524A","#868686")
jco<-c("#223D6C" ,"#D20A13", "#FFD121", "#088247")
mycol
pdf(file = "consensusmap.pdf",width = 4,height = 4)
consensusmap(mut.nmf,
             labRow = NA,
             labCol = NA,
             annCol = data.frame("cluster"=group[colnames(nmf.input)]),
             annColors = list(cluster=c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4])))
invisible(dev.off())
###########
pdf(file = "basismap_CAF_NMF.pdf",width = 4,height = 4,onefile = FALSE)
basismap(mut.nmf,
         cexCol = 1,
         cexRow = 1,
         annColors=list(c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4])))
invisible(dev.off())
##################################################################################
load("group_CAF_SMC_NMF.Rdata")###
#####

load("group_CAF_SMC_NMFallcaf.Rdata")
###########creat a seruat for Cancer associated fibroblasts cells
###########single_dat_stromal_tumor
sce<-sce[,sce$Class=="Tumor"]
sce<-sce[,sce$Cell_subtype %in% c("Stromal 1","Stromal 2","Stromal 3","Myofibroblasts")]
sce<-AddMetaData(sce,metadata=group_CAF, col.name="NMFcluster")
sce@meta.data$NMFcluster<-ifelse(is.na(sce@meta.data$NMFcluster),"4",sce@meta.data$NMFcluster)
DimPlot(sce,reduction = "tsne",group.by="NMFcluster",cols=mycol)
##########
DimPlot(sce,reduction = "umap",group.by="NMFcluster",cols=mycol)
DotPlot(sce,features =m6agenes,group.by="NMFcluster", cols = c("green", "red")) + RotatedAxis()+ coord_flip()
ggsave("CAF_NMF_dotplot_tumor.pdf")
####
sce$NMFcluster<-as.factor(sce$NMFcluster)
Idents(sce)<-sce$NMFcluster
DimPlot(sce,reduction = "tsne",cols=mycol)
###
sce.markers <- FindAllMarkers(sce, 
                              only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
sce.markers<-sce.markers[sce.markers$p_val_adj<0.05,]
VlnPlot(sce,features = sce.markers$gene[sce.markers$gene %in% m6agenes],group.by="NMFcluster",pt.size=0)+ RotatedAxis()
ggsave("CAF_m6Agenes.dotplot.pdf")
DotPlot(sce,features = sce.markers$gene[sce.markers$gene %in% m6agenes],group.by="NMFcluster")+ RotatedAxis()
sce.markers  %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
sce.markers_top5 <- sce.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)
########
newphe<-data.frame(read_excel("clinicalphe.xlsx",1, col_names= TRUE, col_types= NULL, na=" "))##chip data
newphe<-merge(phe,newphe,by="Sample",all=T)
newphe$Index<-gsub("-",".",newphe$Index)
rownames(newphe)<-newphe$Index
sce<-AddMetaData(sce,metadata=newphe, col.name=colnames(newphe))###
colnames(sce@meta.data)
library(dittoSeq)
dittoBarPlot(sce, "NMFcluster", group.by = "MSI",color.panel=mycol)+coord_flip()
ggsave("NMFcluster In_Class.pdf",width=4.36,height =1.8 )

###
nmfscore<-split(sce.markers$gene,sce.markers$cluster)
FeaturePlot(sce,reduction="tsne","nmfscore21") & viridis::scale_color_viridis(option="H")
ggsave("nmfscore21_Tumor_CAF_Class.pdf" )
#nmfscore<-nmfscore[2]
sce<- AddModuleScore(object = sce, features = nmfscore, name = "nmfscore") #AddModuleScore计算得分
colnames(sce@meta.data)
cordata<-data.frame(sce@meta.data[,39:42],sce@meta.data[,33:37])
r<-cor(cordata)
colnames(r)
pdf("CAFCCRType_CAF_NMF_tumor_correlation.pdf",height = 4.2,width=4)
pheatmap(r[c(5:9),c(1:4)],
         cluster_rows =T, 
         cluster_cols = F,  
         show_rownames = T, show_colnames = T,  
         fontsize_row =12,
         color=colorRampPalette(c("blue","white","red"))(50))
dev.off()
##
####single   density
tt_df<-data.frame(sce$NMFcluster,sce$Cluster1)
rownames(tt_df)<-sce$Index
colnames(tt_df)[1]<-c("clusterID")
data=tt_df
data$clusterID<-as.factor(data$clusterID)
data$sce.Cluster1<-as.numeric(data$sce.Cluster1)
#pdf(file = "density_cluster.pdf", width = 16, height = 6)
#par(mfrow=c(2,4))
#dev.off()
table(data$clusterID)
#
col <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02') 
par(mfrow=c(2,4))
plot(density(data[which(data$clusterID == 1),2]), 
     xlim = c(-1,2),ylim = c(0, 1.5), lwd=5,
     las=2, main="", xlab="", col=col[1])
den <- density(data[which(data$clusterID ==2),2])
lines(den$x, den$y, col=col[2], lwd=5)
den <- density(data[which(data$clusterID == 3),2])
lines(den$x, den$y, col=col[3], lwd=5)
den <- density(data[which(data$clusterID == 4),2])
lines(den$x, den$y, col=col[4], lwd=5)
den <- density(data[which(data$clusterID == 5),2])
lines(den$x, den$y, col=col[5], lwd=5)

########### ##################1. comparison with CCR caf type in the data
CAFtypeScore<-data.frame(read_excel("CAFtypeScore.xlsx",1, col_names= TRUE, col_types= NULL, na=" "))##chip data
CAFtypeScore<-split(CAFtypeScore$Genes,CAFtypeScore$Cell_type)
sce<- AddModuleScore(object = sce, features = CAFtypeScore, name = names(CAFtypeScore)) #AddModuleScore计算得分
colnames(sce@meta.data)
features_caf<-names(CAFtypeScore)
##################
#
###############################
#features_caf_for_heatmap <- AverageExpression(sce, group.by = "NMFcluster", features = features_caf,verbose = TRUE)
features_caf_for_heatmap<-aggregate(x=sce@meta.data[,c("pan.dCAF1","pan.iCAF2","pan.iCAF.23","pan.myCAF4","pan.pCAF5")],     
          # Specify group indicator
          by = list(sce@meta.data$NMFcluster),      
          # Specify function (i.e. mean)
          FUN = mean)
features_caf_for_heatmap<-data.frame(features_caf_for_heatmap)
features_caf_for_heatmap<-features_caf_for_heatmap[,-1]
rownames(features_caf_for_heatmap)<-c("C1", "C2" ,"C3","C4")
colnames(features_caf_for_heatmap)<-names(CAFtypeScore)
#features_caf_for_heatmap <- colMeans(sce@meta.data$`pan-dCAF`, "NMFcluster")
#pan_caf_file_proteins_for_heatmap <- pan_caf_file[match(proteins, row.names(pan_caf_file)),]
features_caf_for_heatmap <- na.omit(features_caf_for_heatmap)
pdf("CAFCCRType_CAF_NMF_tumor.pdf",height = 4.2,width=4)
pheatmap(t(features_caf_for_heatmap), 
         scale = "row", 
         cluster_rows =T, 
         cluster_cols = T,  
         show_rownames = T, show_colnames = T,  
         fontsize_row =12,
         color=colorRampPalette(c("blue","white","red"))(50))
dev.off()
##################################2.calcuated the cell number for class CAFs
library(janitor)
meta.datas<-sce@meta.data
#a<-list()
#a[[1]]<-tabyl(meta.datas, NMFcluster)
#a[[2]]<-tabyl(meta.datas, Class)
#a[[3]]<-tabyl(meta.datas, MSI)
#as<-rbindlist(a,use.names=FALSE, fill=FALSE, idcol=NULL)
#as<-cbind(as,splitgroups)
as<-tabyl(meta.datas, NMFcluster)
colnames(as)[1]<-"Subgroup"
library(ComplexHeatmap)
library(circlize)
AveExpression <- AverageExpression(sce, assays = "RNA",group.by = "NMFcluster", features =m6agenes,verbose = TRUE) %>% .$RNA
AveExpressions<-t(AveExpression)
AveExpressions<-scale(AveExpressions)
rownames(AveExpressions)
col_fun = colorRamp2(c(0, 0.3, 0.8), c("blue", "white", "red"))
rowAnn<-HeatmapAnnotation(Percent = as$percent, No.cells = anno_barplot(as$n),col=list(Percent = col_fun))
##################################
pdf("CAF_cell_m6aGenes_Heatmap_m6aGroup.pdf",width = 5.42,height=4.83)
Heatmap(t(AveExpressions),name = "AveE.",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        #row_names_size=6,
        #row_split = as$splitgroups ,
        row_names_side =  "left",
        top_annotation = rowAnn
)
dev.off()
###################

#########
###224 marker Figure YA
celltype <- c("Zero","High_M","Median","Low_M")
colourCount = length(unique(sce@meta.data$NMFcluster))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
celltype_colors <- getPalette(colourCount)
heatmap_gene<-m6agenes
heatmap_gene <- c("CD3E", "CD8A","CD14", "LYZ","FCGR3A", "MS4A7","MS4A1","GNLY", "NKG7","FCER1A", "CST3", "PPBP")
heatmap_AveE <- AverageExpression(sce, assays = "RNA",group.by = "NMFcluster", features = heatmap_gene,verbose = TRUE) %>% .$RNA
gene_num <- c(8,6,5,4)###number of genes
gaps_row <- cumsum(gene_num)
cluster_num <- c(1,1,1,1)##numer of each clusters for cluster
gaps_col <- cumsum(cluster_num)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
annotation_row <- data.frame(row.names = rownames(heatmap_AveE),
                             `CellType` = rep(factor(celltype,levels = celltype),gene_num))
annotation_col <- data.frame(row.names = colnames(heatmap_AveE),
                             `CellType` = rep(factor(celltype,levels = celltype),cluster_num))
annotation_colors = list(`CellType` = celltype_colors)
names(annotation_colors$`CellType`) = celltype
ph3 <- pheatmap(heatmap_AveE,
                cluster_cols =T,
                cluster_rows =T,
                show_colnames=F,
                show_rownames=T,
                border=F,#border_color = "white",
                color = c(colorRampPalette(colors = c("#2166ac","#f7fbff"))(length(bk)/2),
                          colorRampPalette(colors = c("#f7fbff","#b2182b"))(length(bk)/2)),
                breaks=bk,
                scale="row",
                legend_breaks=seq(-2,2,2),
                gaps_row = gaps_row,
                gaps_col = gaps_col,
                annotation_row = annotation_row,annotation_col = annotation_col,
                annotation_colors = annotation_colors,
                annotation_names_row = F,annotation_names_col = T)
require(ggplotify)
pB3 = as.ggplot(ph3) ### 将pheatmap对象转为ggplot对象，便于后续拼图
pB_1_2 <- pB1 + pB2 + plot_layout(ncol = 1, heights = c(1, 2))
#pB_all <- pB1 + pB2 + pB3 + plot_layout(ncol = 1, heights = c(1, 2,3))
###########################
############
####
sce$NMFcluster<-as.factor(sce$NMFcluster)
Idents(sce)<-sce$NMFcluster
sce.markers <- FindAllMarkers(sce,
                                         only.pos = TRUE, 
                                         min.pct = 0.15, logfc.threshold = 0.15)

write.csv(sce.markers,"sce.markers_SMC_CAF_NMF_Tumor.csv")###
sce.markers  %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
sce.markers_top5 <- sce.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)
#######
###
sigmethy<-split(sce.markers_top20$gene,sce.markers_top20$cluster)
sce<- AddModuleScore(object = sce, features = sigmethy, name = "Cluster") #AddModuleScore计算得分
colnames(sce@meta.data)
data.frame(table(sce$Sample))
table(sce$Sample=="NA")
####
sce@meta.data$NMFcluster<-as.factor(sce@meta.data$NMFcluster)
C1<-dittoPlot(sce, "Cluster1", group.by = "Sample",
              #plots = c("jitter", "vlnplot", "boxplot"), # <- order matters
              plots = c( "vlnplot", "boxplot"), # <- order matters
              vlnplot.lineweight = 0.6,
              # change the color and size of jitter points
              #jitter.color = "blue", jitter.size = 0.0,
              # change the outline color and width, and remove the fill of boxplots
              boxplot.color = "white", boxplot.width = 0.2,
              boxplot.fill = FALSE
              # change how the violin plot widths are normalized across groups
              # vlnplot.scaling = "count"
)& NoLegend()
C2<-dittoPlot(sce, "Cluster2", group.by = "Sample",
              #plots = c("jitter", "vlnplot", "boxplot"), # <- order matters
              plots = c( "vlnplot", "boxplot"), # <- order matters
              vlnplot.lineweight = 0.6,
              # change the color and size of jitter points
              #jitter.color = "blue", jitter.size = 0.0,
              # change the outline color and width, and remove the fill of boxplots
              boxplot.color = "white", boxplot.width = 0.2,
              boxplot.fill = FALSE
              # change how the violin plot widths are normalized across groups
              # vlnplot.scaling = "count"
)& NoLegend()
C3<-dittoPlot(sce, "Cluster3", group.by = "Sample",
              #plots = c("jitter", "vlnplot", "boxplot"), # <- order matters
              plots = c( "vlnplot", "boxplot"), # <- order matters
              vlnplot.lineweight = 0.6,
              # change the color and size of jitter points
              #jitter.color = "blue", jitter.size = 0.0,
              # change the outline color and width, and remove the fill of boxplots
              boxplot.color = "white", boxplot.width = 0.2,
              boxplot.fill = FALSE
              # change how the violin plot widths are normalized across groups
              # vlnplot.scaling = "count"
)& NoLegend()
C4<-dittoPlot(sce, "Cluster4", group.by = "Sample",
              #plots = c("jitter", "vlnplot", "boxplot"), # <- order matters
              plots = c( "vlnplot", "boxplot"), # <- order matters
              vlnplot.lineweight = 0.6,
              # change the color and size of jitter points
              #jitter.color = "blue", jitter.size = 0.0,
              # change the outline color and width, and remove the fill of boxplots
              boxplot.color = "white", boxplot.width = 0.2,
              boxplot.fill = FALSE
              # change how the violin plot widths are normalized across groups
              # vlnplot.scaling = "count"
)& NoLegend()
C1+C2+C3+C4
ggsave("NMF_group_score_CAF_plot.pdf")
FeaturePlot(sce,reduction="tsne",c("Cluster1","Cluster2", "Cluster3","Cluster4")) & viridis::scale_color_viridis(option="H")
ggsave("Featureplot_for_CAF_ClusterScore.pdf")
RidgePlot(sce,c("Cluster1","Cluster2", "Cluster3","Cluster4"),ncol=2) & viridis::scale_color_viridis(option="H")
####
#logFCfilter<-0
#adjPvalFilter<-0.01
#sig.markers=sce.markers[(abs(as.numeric(as.vector(sce.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(sce.markers$p_val_adj))<adjPvalFilter),]
DoHeatmap(sce,features = sce.markers_top5$gene,group.by="NMFcluster")
ggsave("DoHeatmap_CAF_NMF.pdf",height =2.8,width=6.6)
DotPlot(sce,features = sce.markers_top5$gene,group.by="NMFcluster")+ RotatedAxis()
ggsave("DoDotPlot_CAF_NMF_tumor_5.pdf",height =2.91,width=6.6)
###
#a<-dittoPlot(sce,"HNRNPA2B1",group.by="Sample") & NoLegend()
#b<-dittoPlot(sce,"WTAP",group.by="Sample") & NoLegend()
#c<-dittoPlot(sce,"HNRNPC",group.by="Sample") & NoLegend()
#d<-dittoPlot(sce,"NDUFA4L2",group.by="Sample") & NoLegend()
#a+b+c+d
#dittoHeatmap(subset(sce, downsample = 4), m6agenes, use_raster = TRUE,
         #    annot.by = "NMFcluster",color.by=c("red","blue"))

##pathway analysis
Myenrich <- function(genes, category = c("kegg", "gobp"), 
                     geneid = c("SYMBOL", "ENTREZID", "ENSEMBL", "UNIPROT")){
  library(clusterProfiler)
  category <- match.arg(category)
  geneid <- match.arg(geneid)
  if (geneid != "ENTREZID"){
    genes <- bitr(genes, fromType = geneid,toType ="ENTREZID",OrgDb="org.Hs.eg.db") %>% .[, 2] %>% as.character()
  }else{genes <- genes}
  if (category == "kegg"){
    enrich <- enrichKEGG(genes, organism = "hsa",keyType = "kegg",
                         pAdjustMethod = "BH",pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05, maxGSSize = 5000)
  }
  if (category == "gobp"){
    enrich <- enrichGO(genes, OrgDb="org.Hs.eg.db",ont= "BP",pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,qvalueCutoff  = 0.05, readable = TRUE)
  }
  return(enrich)
}
##output the gene set for each clusters of NMF in CAFs 
ProSig <- split(sce.markers$gene, sce.markers$cluster) 
ProSig<-ProSig[2]


Sigbp <- lapply(ProSig, Myenrich, category = "gobp", geneid = "SYMBOL")
#
library(ggplot2)
library(ggpubr)
#Sigbpplot <- lapply(names(Sigbp), function(z) barplot(Sigbp[[z]], showCategory = 10) +scale_x_discrete(labels=function(x) str_wrap(x, width=40)) +
 #                     ggtitle(z)  + theme(axis.text.x = element_text(angle = 90)))
#Sigbpplot
###

names(Sigbp)<-c("C1","C2","C3","C4")
GOs<-c()
for (i in names(Sigbp)){
  GO<-data.frame(Sigbp[[i]])
  GO<-GO[order(GO$Count,decreasing = T),][c(1:5),]
  GO<-GO[,c("Description","p.adjust")]
  GO<-cbind( GO,set=i)
  GOs<-rbind(GOs,GO)
}
#######
orderDescription<-GOs$Description
GOs<-GOs[GOs$p.adjust<0.01,]
#dev.off()
library(forcats)
GOs %>% 
  mutate(Description = fct_reorder(Description, set)) %>%
  #mutate(set = fct_relevel(set, "C4","C3","C2","C1")) %>%
  ggplot(aes(set, Description, fill=-log10(p.adjust))) + 
  geom_tile(colour = "white") + 
  scale_y_discrete(labels=function(y) str_wrap(y, width=50))+
  #facet_grid(year~monthf) + 
  #theme_bw()+# 不要背景
  theme_classic()+
  #theme_minimal()+
  #theme_light()+
  #theme_dark()+
  #theme_void()+
  scale_fill_gradient(low="yellow", high="red") +
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(angle = 35, hjust = 1),# 调整x轴文字
        axis.text.y = element_text(size = 10),#调整y轴文字
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  coord_trans(y ="reverse")+
  labs(x="Week of Month",
       y="",
       title = "Funcation Heatmap", 
       #subtitle="Yahoo Closing Price", 
       fill="-log(adjustP)")
ggsave("CAF_GOmap_NMF.pdf",height=3.8,width = 5.87)
  #theme(aspect.ratio=1)
#colnames(GO)
#GOlists<-rbindlist(GOlist,use.names=FALSE, fill=FALSE, idcol=NULL)
#as<-cbind(as,splitgroups)

Sigkegg <- lapply(ProSig, Myenrich, category = "kegg", geneid = "SYMBOL")
names(Sigkegg)<-c("C1","C2","C3","C4")
Keggs<-c()
for (i in names(Sigkegg)){
  Kegg<-data.frame(Sigkegg[[i]])
  Kegg<-Kegg[order(Kegg$Count,decreasing = T),][c(1:5),]
  Kegg<-Kegg[,c("Description","p.adjust")]
  Kegg<-cbind(Kegg,set=i)
  Keggs<-rbind(Keggs,Kegg)
}
#dev.off()
#Keggs<-data.frame(Keggs)
Keggs<-na.omit(Keggs)
Keggs<-Keggs[-10,]
#Keggs$p.adjust<-as.numeric(Keggs$p.adjust)
Keggs<-Keggs[Keggs$p.adjust<0.01,]
#Keggs$Description<-as.factor(Keggs$Description)
library(forcats)
Keggs %>% 
  mutate(Description = fct_reorder(Description, set)) %>%
  #mutate(set = fct_relevel(set, "C1","C2","C3","C4")) %>%
  ggplot(aes(x=set, y=Description, fill=-log10(p.adjust))) + 
  geom_tile(colour = "white") + 
  scale_y_discrete(labels=function(y) str_wrap(y, width=40))+
  #facet_grid(year~monthf) + 
  #theme_bw()+# 不要背景
  theme_classic()+
  #theme_minimal()+
  #theme_light()+
  #theme_dark()+
  #theme_void()+
  scale_fill_gradient(low="yellow", high="red") +
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(angle = 35, hjust = 1),# 调整x轴文字
        axis.text.y = element_text(size = 10),#调整y轴文字
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  coord_trans(y ="reverse")+
  labs(x="",
       y="",
       title = "KEGG Heatmap", 
       #subtitle="Yahoo Closing Price", 
       fill="-log(adjustP)")
ggsave("CAF_KEGGmap_NMF.pdf",height=3.8,width = 5)
########################################################


#Sigkeggplot <- lapply(names(Sigkegg), function(z)barplot(Sigkegg[[z]], showCategory = 10) +
                        #ggtitle(z) + scale_x_discrete(labels=function(x) str_wrap(x, width=30)))
#####################################
pan_caf_file=as.matrix(sce@assays$RNA@data)
##Qian et al functional gene sets
#pan_caf_file <- read.table(file.choose(), header = T, row.names = 1, sep = '\t')
proteins <- c('COL1A1', 'COL1A2', 'COL3A1', 'COL4A1', 'COL4A2', 'COL5A1', 'COL5A2', 'COL6A1', 
              'COL7A1', 'COL8A1', 'COL10A1', 'COL11A1', 'COL12A1', 'COL13A1', 'COL14A1',
              'COL15A1', 'COL16A1', 'COL18A1', 'BGN', 'DCN', 'LUM', 'TAGLN','ELN', 'FN1',
              'MMP1', 'MMP2', 'MMP3', 'MMP9', 'MMP10', 'MMP11', 'MMP14', 'MMP19', 'SERPINE1', 
              'CTHRC1', 'THBS2', 'SULF1', 'TGFBI', 'COMP', 'INHBA', 'EGFL6', 'ANGPT2', 'PDGFA',
              'PDGFC', 'VEGFA', 'ACTA2', 'MYL6', 'MYH9', 'MYH11', 'PLN', 'TPM1', 'TMP2', 'SORBS2', 
              'RRAS', 'RASL12', 'RASGRP2', 'CFD', 'CFI', 'C3', 'C7', 'CCL21', 'CXCL14', 'CXCL12',
              'IL33', 'CXCL3', 'CXCL2', 'CXCL1', 'CCL2', 'CCXL26',  'IL6', 'IL7')
##
pan_caf_file_proteins_for_heatmap <- AverageExpression(sce, assays = "RNA",group.by = "NMFcluster", features = proteins,verbose = TRUE) %>% .$RNA
#pan_caf_file_proteins_for_heatmap <- pan_caf_file[match(proteins, row.names(pan_caf_file)),]
pan_caf_file_proteins_for_heatmap <- na.omit(pan_caf_file_proteins_for_heatmap)

pdf("pathway_CAF_NMF_tumor.pdf",height = 8.15,width=2.55)
pheatmap(pan_caf_file_proteins_for_heatmap, 
         scale = "row", 
         cluster_rows =F, 
         cluster_cols = T,  
         show_rownames = T, show_colnames = T,  
         fontsize_row = 7,
         color=colorRampPalette(c("blue","white","red"))(50))
dev.off()
#Cell surface
cell_surface_markers <- c('PARM1', 'APOC3', 'CSPG4', 'CD36', 'SUSD2', 'SDC1', 'TREM1', 
                          'THY1', 'BAMBI', 'P2RY6', 'CD34', 'PI16', 'GPC3', 'RAMP2', 'LRRN3', 
                          'CLDN1', 'ICAM1', 'MUSK', 'CXCR4', 'NCAM1', 'EFNB1', 'EBP')
#pan_caf_file <- read.table(file.choose(), header = T, row.names = 1, sep = '\t')
#colnames_to_declare <- c('Pan-myCAF', 'Pan-dCAF', 'Pan-iCAF', 'Pan-iCAF-2', 'Pan-nCAF', 'Pan-nCAF-2', 'Pan-pCAF')

#colnames(pan_caf_file) <- colnames_to_declare
##pan_caf_file_receptors <- pan_caf_file[match(cell_surface_markers, row.names(pan_caf_file)),]
#pan_caf_file_receptors <- na.omit(pan_caf_file_receptors)
#pan_caf_file_receptors <- pan_caf_file_receptors[, -6]
##
pan_caf_epic_file_receptors_heatmap <- AverageExpression(sce, assays = "RNA",
                                                         group.by = "NMFcluster",
                                                         features = cell_surface_markers ,verbose = TRUE) %>% .$RNA
pan_caf_epic_file_receptors_heatmap <- na.omit(pan_caf_epic_file_receptors_heatmap)
###
pdf("cell_surface_CAF_NMF.pdf",height = 3.94,width=3.02)
pheatmap(pan_caf_epic_file_receptors_heatmap[-2,],
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = T,  
         show_rownames = T, angle_col = c('45'),
         color=colorRampPalette(c("blue","white","red"))(20))
dev.off()
#ggsave("cell_surface_CAF_NMF.pdf",height = 4.94,width=3.02)
#TF analysis
#TFgenes<-read_xlsx("TFgenes.xlsx",1)
TFgenes<-data.frame(read_excel("TFgenes.xlsx",1, col_names= TRUE, col_types= NULL, na=" "))##TFgenes
#####
caf_tfgenes_file_receptors_heatmap <- AverageExpression(sce, assays = "RNA",
                                                         group.by = "NMFcluster",
                                                         features = sce.markers$gene ,verbose = TRUE) %>% .$RNA
caf_tfgenes_file_receptors_heatmap<-caf_tfgenes_file_receptors_heatmap[rownames(caf_tfgenes_file_receptors_heatmap) %in% TFgenes$Name,]
##
caf_tfgenes_file_receptors_heatmap <- na.omit(caf_tfgenes_file_receptors_heatmap )
###
pdf("TFgenes_CAF_NMF.pdf",height =6,width=3.02)
pheatmap(caf_tfgenes_file_receptors_heatmap,
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = T,  
         show_rownames = T, angle_col = c('45'),
         color=colorRampPalette(c("blue","white","red"))(20))
dev.off()

############
colnames(sce@meta.data)
###Featureplot

#############
library(scater)
library(dittoSeq)
dittoDimPlot(sce,reduction = "tsne", "NMFcluster")
ggsave("NMFgroup_CAF_tumor_GSE132465.pdf",height=4.66,width=6.03)
################
save(ProSig,file="ProSig_CAF_NNMF.Rdata")
