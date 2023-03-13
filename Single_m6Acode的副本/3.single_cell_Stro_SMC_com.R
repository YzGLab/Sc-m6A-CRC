##
#extract the fbi cells for research
##
library(data.table)
#################################
load("GSE132465_single_stromal.Rdata")
GSE132465_single<-fread("GSE132465_GEO_processed_CRC_10X_natural_log_TPM_matrix.txt")
phe<-read.delim2("GSE132465_GEO_processed_CRC_10X_cell_annotation.txt")
table(phe$Cell_subtype)
table(phe$Cell_type)
#'Stromal 1','Stromal 2','Stromal 3' ,'Myofibroblasts'
phe<-phe[phe$Cell_subtype %in% c('Stromal 1','Stromal 2','Stromal 3' ,'Myofibroblasts'),]
phe$Index<-gsub("-",".",phe$Index)
rownames(phe)<-phe$Index
GSE132465_single_cafs<-GSE132465_single_stromal[,colnames(GSE132465_single_stromal) %in% phe$Index ]
GSE132465_single_cafs[1:5,1:5]
#GSE144735
################################
GSE144735_single <- fread("GSE144735_processed_KUL3_CRC_10X_natural_log_TPM_matrix.txt")##
GSE144735_single[1:5,1:5]
##
GSE144735_phe<-read.delim2("GSE144735_processed_KUL3_CRC_10X_annotation.txt")
GSE144735_phe$Index<-gsub("-",".",GSE144735_phe$Index)
rownames(GSE144735_phe)<-GSE144735_phe$Index
colnames(GSE144735_phe)
table(GSE144735_phe$Cell_type)
table(GSE144735_phe$Cell_subtype)
#GSE144735_phe<-GSE144735_phe[GSE144735_phe$Sample %like% "T" & GSE144735_phe$Cell_type=="Myofibroblasts",]
GSE144735_phe<-GSE144735_phe[GSE144735_phe$Cell_subtype %in% c('Stromal 1','Stromal 2','Stromal 3' ,'Myofibroblasts'),]
GSE144735_single<-data.frame(GSE144735_single)
GSE144735_single_cafs<-GSE144735_single[,colnames(GSE144735_single) %in% GSE144735_phe$Index]
rownames(GSE144735_single_cafs)<-GSE144735_single$Index
####GSE146771_single_cafs,
save(GSE132465_single_cafs,GSE144735_single_cafs,file="Pan_single_cafs_CRC.Rdata")###############
################################
###################进行 seruat
load("Pan_single_cafs_CRC.Rdata")
GSE144735_phe<-read.delim2("GSE144735_processed_KUL3_CRC_10X_annotation.txt")
GSE132465_phe<-read.delim2("GSE132465_GEO_processed_CRC_10X_cell_annotation.txt")############


#######
GSE132465_colon_caf <- CreateSeuratObject(counts = GSE132465_single_cafs, project = "GSE132465_CAF_colon")
dim(GSE132465_colon_caf) #results are 61 cells and 14494 features
##
GSE132465_colon_caf <- FindVariableFeatures(GSE132465_colon_caf, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(GSE132465_colon_caf), 10)
top10
plot1 <- VariableFeaturePlot(GSE132465_colon_caf)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot2))
all.genes <- rownames(GSE132465_colon_caf)
GSE132465_colon_caf <- ScaleData(GSE132465_colon_caf, features = all.genes)
GSE132465_colon_caf <- RunPCA(GSE132465_colon_caf, features = VariableFeatures(object = GSE132465_colon_caf))
DimHeatmap(GSE132465_colon_caf, dims = 1:10, cells = 500, balanced = TRUE)
DimHeatmap(GSE132465_colon_caf, dims = 11:20, cells = 500, balanced = TRUE)
DimHeatmap(GSE132465_colon_caf, dims = 21:30, cells = 500, balanced = TRUE)
DimHeatmap(GSE132465_colon_caf, dims = 31:40, cells = 500, balanced = TRUE)
DimHeatmap(GSE132465_colon_caf, dims = 41:50, cells = 500, balanced = TRUE)
GSE132465_colon_caf <- JackStraw(GSE132465_colon_caf, num.replicate = 10)
GSE132465_colon_caf <- ScoreJackStraw(GSE132465_colon_caf, dims = 1:10)
#JackStrawPlot(GSE132465_colon_caf, dims = 1:10)
ElbowPlot(GSE132465_colon_caf, ndims = 50, reduction = "pca")
GSE132465_colon_caf <- FindNeighbors(GSE132465_colon_caf, dims = 1:20)
##
GSE132465_colon_caf  <- FindClusters(
  object = GSE132465_colon_caf,
  resolution = c(seq(.1,1.6,.2))
)
clustree(GSE132465_colon_caf@meta.data, prefix = "RNA_snn_res.")
colnames(GSE132465_colon_caf@meta.data)
table(GSE132465_colon_caf$seurat_clusters)
##########
GSE132465_colon_caf_resolution_0.5 <- FindClusters(GSE132465_colon_caf, resolution = 0.2)
##################
GSE132465_colon_caf_resolution_0.5 <- RunUMAP(GSE132465_colon_caf_resolution_0.5, dims = 1:10)
DimPlot(GSE132465_colon_caf_resolution_0.5, reduction = "umap", cols=mycol)

GSE132465_colon_caf_resolution_0.5 <- RunTSNE(object = GSE132465_colon_caf_resolution_0.5, dims = 1:10, do.fast = TRUE)
DimPlot(GSE132465_colon_caf_resolution_0.5,reduction = "tsne",label=T)
#set.seed(11000)
#reducedDim(SE132465_colon_cafs_resolution_0.5, "force") <- igraph::layout_with_fr(g)
#plotReducedDim(sce.caf.combined.resolution.02, colour_by="label", dimred="force")
GSE132465_colon_caf_resolution_0.5.markers <- FindAllMarkers(GSE132465_colon_caf_resolution_0.5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
GSE132465_colon_caf_resolution_0.5.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 2, wt = avg_log2FC)
GSE132465_colon_caf_resolution_0.5.markers_top20 <- GSE132465_colon_caf_resolution_0.5.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)
DoHeatmap(GSE132465_colon_caf_resolution_0.5, features = GSE132465_colon_caf_resolution_0.5.markers_top20$gene, label = F) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white")
feature <- GSE132465_colon_caf_resolution_0.5.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 2, wt = avg_log2FC)
feature<-data.frame(feature)
feature<-feature$gene
colon1_dot_plot <- DotPlot(GSE132465_colon_caf_resolution_0.5, features = feature, cols = c("lightgrey", "red")) + RotatedAxis()
colon1_dot_plot
DotPlot(GSE132465_colon_caf_resolution_0.5, features = unique(feature)) + RotatedAxis()
DoHeatmap(subset(GSE132465_colon_caf_resolution_0.5, downsample = 100), features = feature, size = 3)
##################################
###GSE144735
#######
GSE144735_colon_caf <- CreateSeuratObject(counts = GSE144735_single_cafs, project = "GSE144735_CAF_colon")
dim(GSE144735_colon_caf) #results are 61 cells and 14494 features
GSE144735_colon_caf <- FindVariableFeatures(GSE144735_colon_caf, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(GSE144735_colon_caf), 10)
top10
plot1 <- VariableFeaturePlot(GSE144735_colon_caf)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot2))
all.genes <- rownames(GSE144735_colon_caf)
GSE144735_colon_caf <- ScaleData(GSE144735_colon_caf, features = all.genes)
GSE144735_colon_caf <- RunPCA(GSE144735_colon_caf, features = VariableFeatures(object = GSE144735_colon_caf))
DimHeatmap(GSE144735_colon_caf, dims = 1:10, cells = 500, balanced = TRUE)
DimHeatmap(GSE144735_colon_caf, dims = 11:20, cells = 500, balanced = TRUE)
DimHeatmap(GSE144735_colon_caf, dims = 21:30, cells = 500, balanced = TRUE)
DimHeatmap(GSE144735_colon_caf, dims = 31:40, cells = 500, balanced = TRUE)
DimHeatmap(GSE144735_colon_caf, dims = 41:50, cells = 500, balanced = TRUE)
GSE144735_colon_caf <- JackStraw(GSE144735_colon_caf, num.replicate = 10)
GSE144735_colon_caf <- ScoreJackStraw(GSE144735_colon_caf, dims = 1:8)
JackStrawPlot(GSE144735_colon_caf, dims = 1:8)
ElbowPlot(GSE144735_colon_caf, ndims = 50, reduction = "pca")
GSE144735_colon_caf <- FindNeighbors(GSE144735_colon_caf, dims = 1:20)
##########
GSE144735_colon_caf_resolution_0.5 <- FindClusters(GSE144735_colon_caf, resolution = 0.3)
##################
GSE144735_colon_caf_resolution_0.5 <- RunUMAP(GSE144735_colon_caf_resolution_0.5, dims = 1:8)
DimPlot(GSE144735_colon_caf_resolution_0.5, reduction = "umap", cols=mycol)
GSE144735_colon_caf_resolution_0.5 <- RunTSNE(object = GSE144735_colon_caf_resolution_0.5, dims = 1:10, do.fast = TRUE)
DimPlot(GSE144735_colon_caf_resolution_0.5,reduction = "tsne",label=T, cols=mycol)
#
GSE144735_colon_caf_resolution_0.5.markers <- FindAllMarkers(GSE144735_colon_caf_resolution_0.5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
GSE144735_colon_caf_resolution_0.5.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 2, wt = avg_log2FC)
GSE144735_colon_caf_resolution_0.5.markers_top20 <- GSE144735_colon_caf_resolution_0.5.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)
DoHeatmap(GSE144735_colon_caf_resolution_0.5, features = GSE144735_colon_caf_resolution_0.5.markers_top20$gene, label = F) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white")
feature <- GSE144735_colon_caf_resolution_0.5.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 2, wt = avg_log2FC)
feature<-data.frame(feature)
feature<-feature$gene
colon2_dot_plot <- DotPlot(GSE144735_colon_caf_resolution_0.5, features = feature, cols = c("lightgrey", "red")) + RotatedAxis()
colon2_dot_plot
####
DotPlot(GSE144735_colon_caf_resolution_0.5, features = unique(feature)) + RotatedAxis()
DoHeatmap(subset(GSE144735_colon_caf_resolution_0.5, downsample = 100), features = feature, size = 3)
##############
######################################################

#merge
###
DimPlot(GSE132465_colon_caf_resolution_0.5, reduction = "umap", cols=mycol)
DimPlot(GSE144735_colon_caf_resolution_0.5, reduction = "umap", cols=mycol)
##CAF Integration
#hnsc_human_caf <- readRDS(file.choose())
colon1_human_caf<-GSE132465_colon_caf_resolution_0.5
#colon1_cafs_ids <- c("colon1 CAFs 1", "colon1 CAFs 2", "colon1 CAFs 3")
#names(colon1_cafs_ids) <- levels(colon1_human_caf)
#colon1_human_caf <- RenameIdents(colon1_human_caf, colon1_cafs_ids)
DimPlot(colon1_human_caf, cols=mycol)
#caf2_human_caf <- readRDS(file.choose())
colon2_human_caf<-GSE144735_colon_caf_resolution_0.5
#colon2_cafs_ids <- c("colon2 CAFs 1", "colon2 CAFs 2","colon2 CAFs 3")
#names(colon2_cafs_ids) <- levels(colon2_human_caf)
#colon2_human_caf<-RenameIdents(colon2_human_caf,colon2_cafs_ids)
DimPlot(colon2_human_caf, cols=mycol)

###############
save(colon1_human_caf,colon2_human_caf,file="CRC_caf.combined.Rdata")
####merge 合并
caf.anchors <- FindIntegrationAnchors(object.list = list(colon1_human_caf,colon2_human_caf),dims = 1:10, k.filter = 150)
caf.combined <- IntegrateData(anchorset =caf.anchors, dims = 1:30,k.weight=30)
#features <- SelectIntegrationFeatures(object.list =  list(colon1_human_caf,colon2_human_caf),                                      
#                                     nfeatures = 2000)
#seob_list <- PrepSCTIntegration(object.list =  list(colon1_human_caf,colon2_human_caf),                                 
#                              anchor.features = features)
## 找 anchors，15到30分钟
#anchors <- FindIntegrationAnchors(object.list = seob_list,                                  
# reference = 3 # 当有多个样本时，制定一个作为参考可加快速度                                 
# normalization.method = "SCT", # 选择方法“SCT”                                  
#anchor.features = features)
#caf.combined<- IntegrateData(anchorset = anchors,                       
#normalization.method = "SCT")
#caf.combined <- IntegrateData(anchorset =anchors, dims = 1:50,k.weight=40)
DefaultAssay(caf.combined) <- "integrated"
caf.combined <- ScaleData(caf.combined, verbose = FALSE)
caf.combined <- RunPCA(caf.combined, features = VariableFeatures(object = caf.combined))
DimPlot(caf.combined, reduction = "pca")
DimHeatmap(caf.combined, dims = 1:10, cells = 500, balanced = TRUE)
DimHeatmap(caf.combined, dims = 11:20, cells = 500, balanced = TRUE)
DimHeatmap(caf.combined, dims = 21:30, cells = 500, balanced = TRUE)
DimHeatmap(caf.combined, dims = 31:40, cells = 500, balanced = TRUE)
DimHeatmap(caf.combined, dims = 41:50, cells = 500, balanced = TRUE)
caf.combined <- JackStraw(caf.combined, num.replicate = 20)
caf.combined <- ScoreJackStraw(caf.combined, dims = 1:20)
JackStrawPlot(caf.combined, dims = 1:20)
ElbowPlot(caf.combined, ndims = 50)
caf.combined <- FindNeighbors(caf.combined, reduction = "pca", dims = 1:20)
caf.combined <- RunTSNE(object = caf.combined, dims = 1:20, do.fast = TRUE)
caf.combined <- RunUMAP(caf.combined, reduction = "pca", dims = 1:20)
##
#Resolution 0.2 PCs 1 through 30
caf.combined.resolution.02 <- FindClusters(caf.combined, resolution = 0.1)
DimPlot(caf.combined.resolution.02, reduction = "umap", label = TRUE, cols=mycol)
DimPlot(caf.combined.resolution.02, reduction = "tsne", label = TRUE, cols=mycol)
##
phe<-rbind(GSE132465_phe,GSE144735_phe)
phe$Index<-gsub("-",".",phe$Index)
rownames(phe)<-phe$Index
colnames(caf.combined.resolution.02@meta.data)
caf.combined.resolution.02 <- AddMetaData(caf.combined.resolution.02, phe, col.name = colnames(phe))
DimPlot(caf.combined.resolution.02, reduction = "tsne", group.by= "Cell_type", label = TRUE, cols=mycol)
DimPlot(caf.combined.resolution.02,reduction = "tsne",group.by="Cell_subtype",label=T, cols=mycol)
DimPlot(caf.combined.resolution.02,reduction = "tsne",label=T, cols=mycol)
#############################################################################################
##library 
library(readxl)
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
############
single_dat_sce<-caf.combined.resolution.02
##
m6agenes<-c("CBLL1","KIAA1429","METTL14","METTL3", "RBM15","RBM15B", "WTAP",
            "ZC3H13","ELAVL1","FMR1", "HNRNPA2B1", "HNRNPC", "IGF2BP1", 
            "IGF2BP2", "IGF2BP3", "LRPPRC", "YTHDC1", "YTHDC2", "YTHDF1", 
            "YTHDF2","YTHDF3","ALKBH5","FTO")
mycol <- c("#223D6C","#D20A13","#FFD121","#088247",
           "#58CDD9","#5D90BA","#431A3D","#11AA4D","#91612D","#6E568C","#7A142C",
           "#E0367A","#D8D155","#64495D","#7CC767")
#NMF for single_dat_stromal
pan_caf_file=as.matrix(single_dat_sce@assays$RNA@data)
##使用NMF对不同类型的caf细胞聚类
m6agenes<-c("CBLL1","KIAA1429","METTL14","METTL3", "RBM15","RBM15B", "WTAP",
            "ZC3H13","ELAVL1","FMR1", "HNRNPA2B1", "HNRNPC", "IGF2BP1", 
            "IGF2BP2", "IGF2BP3", "LRPPRC", "YTHDC1", "YTHDC2", "YTHDF1", 
            "YTHDF2","YTHDF3","ALKBH5","FTO")

single_dat<-pan_caf_file[m6agenes,]
####
phe[phe$Cell_subtype== patients,]
phe$Cell_subtype2<-ifelse(phe$Cell_subtype=="Myofibroblasts",'Stromal 4',phe$Cell_subtype)
topn <- "m6A"
ranks <- 4##
patients<-c('Stromal 1','Stromal 2','Stromal 3' ,'Stromal 4')
i<-patients[1]
# NMF
for (i in patients){
  if (!dir.exists("nmfSingle_colon_caf")){
    dir.create("./nmfSingle_colon_caf")
  }
  if (!dir.exists(paste0("nmfSingle_colon_caf/", i))){
    dir.create(paste0("nmfSingle_colon_caf/", i))
  }
  print(i)
  #nmfdat <- single_dat[, str_detect(colnames(single_dat), i)]
  nmfdat <- single_dat[, phe[phe$Cell_subtype2 %in% i,]$Index]
  ## Step 02: chose top 6000 genes with highest sds
  nmfdat[nmfdat < 0] <- 0 
  #genesd <- apply(nmfdat, 1, sd)
  #topNgene <- genesd[order(genesd, decreasing = T)][1:topn]
  #nmfdat <- nmfdat[rownames(nmfdat) %in% m6agenes, ]
  nmfdat<-nmfdat[rowSums(nmfdat != 0)>=1,colSums(nmfdat != 0) >= 1]
  # ?"brunet"
  print(system.time(res_4 <- nmf(nmfdat, rank = ranks, method="lee", nrun=1, seed=123456))) 
  signature <- NMF::basis(res_4)
  #consensusdat <- res_4@consensus  #consensus matrix
  si <- silhouette(res_4, what = 'samples')
  sicluster <- si[1:ncol(nmfdat), ]
  #tiff(paste0("./nmfSingle_colon_caf/", i, "/consensus_", topn, ".tiff"), 
  #  width = 6*480, height = 4*480, res = 100)
  #consensusmap(res_4)
  #dev.off()
  #tiff(paste0("./nmfSingle_colon_caf/", i, "/basicmap_", topn, ".tiff"), 
  #  width = 8 * 480, height = 4*480, res = 100)
  #basismap(res_4)
  #dev.off()
  #tiff(paste0("./nmfSingle_colon_caf/", i, "/coefmap_", topn, ".tiff"), 
  # width = 8 * 480, height = 4*480, res = 100)
  # coefmap(res_4)
  # dev.off()
  colnames(signature) <- paste0(i, "_", 1:ranks)
  signature <- as.data.frame(signature)
  # 
  write.table(signature, paste0("./nmfSingle_colon_caf/", i, "/singature", topn, ".txt"),
              sep = "\t")
  write.table(sicluster, paste0("./nmfSingle_colon_caf/", i, "/sicluster", topn, ".txt"),
              sep = "\t")
  #saveRDS(res_4, file = paste0("./nmfSingle_colon_caf/", i, "/res", ranks, "_", topn, ".rds"))
  print(paste0("NMF for the", i, "is Done!"))
}
#####
#NMF program聚类，画图
# 提取program
topRank <- 10
programG <- list()
for (i in 1:length(patients)){
  filedir <- paste0("./nmfSingle_colon_caf/", patients[i], "/singature", topn, ".txt")
  geneloading <- read.table(filedir, header = T, sep = "\t")
  geneloading$maxC <- apply(geneloading, 1, which.max) %>% 
    paste0(patients[i], "_", .)
  topgenelist <- rownames_to_column(geneloading, var = "gene") %>%
    pivot_longer(., cols = starts_with("S" ), 
                 names_to = "program", values_to = "loading")
  #topgenelist <- dplyr::filter(topgenelist, maxC == program) %>% 
    #group_by(maxC) %>% top_n(n = topRank, wt = loading)
  topgenelist <- split(topgenelist$gene, topgenelist$maxC)
  programG <- c(programG, topgenelist)
}
#对fb细胞进行打分，这里我们选择AUCell
single_dat<-single_dat[rownames(single_dat) %in% c(m6agenes),]
cells_rankings <- AUCell_buildRankings(as.matrix(single_dat), nCores = 10, plotStats=TRUE)
#cells_AUC <- AUCell_calcAUC(geneSets = programG, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
cells_AUC <- AUCell_calcAUC(geneSets = programG, cells_rankings, aucMaxRank=5)
programAUC <- getAUC(cells_AUC)
programAUC[1:5,1:5]
###
#save(single_dat,programAUC,cells_AUC,cells_rankings,clus,file="programAUC.Rdata")

clusterK = 4
M <- cor(t(programAUC), method = "pearson")

#用corrplot画图，便于标出各cluster的黑色方框
pdf("scNMF_CAF_combined.pdf")
corrplot(M, 
         method = "color", #用颜色展示相关系数，还可以改为"circle" (default), "square", "ellipse", "number", "pie", "shade"
         order = "hclust", 
         hclust.method = "ward.D2", 
         addrect = clusterK, #画黑色方框
         tl.pos = "n", #"n" means don't add textlabel
         col = rev(brewer.pal(n = 8, name = "RdBu"))) # 运行?brewer.pal查看更多配色方案

dev.off()
##提取program
cororder <- corrMatOrder(M, order = "hclust", hclust.method = "ward.D2")
M.hc <- M[cororder, cororder]
tree <- hclust(as.dist(1 - M.hc), method = "ward.D2")
clus <- cutree(tree, clusterK)
table(clus)
head(clus)
clus
#save(programAUC,clus,file="programAUC.Rdata")##建立seruat>R

# ??ȡsignature
ProSig <- split(names(clus), clus) 
names(ProSig) <- paste0("cafSig", names(ProSig))
ProSig <- lapply(ProSig, function(z){
  programG[which(names(programG) %in% z)] %>% unlist() %>% as.character() %>% 
    unique()
})
sapply(ProSig, length)
###

###
metalist <- split(names(clus), clus) 
patientLoading <- lapply(patients, function(z){
  filedir <- paste0("./nmfSingle_colon_caf/", z, "/singature", topn, ".txt")
  geneloading <- read.table(filedir, header = T, sep = "\t")
  data.frame(Gene = rownames(geneloading), geneloading)
})
AllLoading <- Reduce(function(x, y)merge(x = x, y = y, by = "Gene", all = T), patientLoading)
head(AllLoading)

##############
library(ggplot2)
library(ggpubr)
require(ggplotify)
library(scales)
pB2_df<-data.frame(clus)
patients <- str_split(rownames(pB2_df), "_") %>% sapply(., "[[", 1) 
Number<-str_split(rownames(pB2_df), "_") %>% sapply(., "[[", 2) 
pB2_df<-cbind(pB2_df,Sample=patients,Number)
head(pB2_df)
pB2_df$Number<-as.numeric(pB2_df$Number)
pB2_df$clus<-as.factor(pB2_df$clus)

###extract the data of cluster for all cells in NMF
cell_cluster <- lapply(patients, function(z){
  filedir <- paste0("./nmfSingle_colon_caf/", z, "/sicluster", topn, ".txt")
  cell_clustering <- read.table(filedir, header = T, sep = "\t")
  data.frame(cell_name = rownames(cell_clustering), cell_clustering,z)
})

Allcell_clusters<-c()
for(i in c(1:length(cell_cluster))){
  mycell_cluster<-cell_cluster[[i]]
  mycell_cluster<-data.frame(mycell_cluster)
  rownames(mycell_cluster)<-c()
  Allcell_clusters<- rbind(Allcell_clusters,mycell_cluster)
}
head(Allcell_clusters)
dim(Allcell_clusters)
Allcell_clusters$clusters<-paste0(Allcell_clusters$z,"_",Allcell_clusters$cluster)
head(Allcell_clusters)
##
pB2_df<-data.frame(clus)
patients <- str_split(rownames(pB2_df), "_") %>% sapply(., "[[", 1) 
Number<-str_split(rownames(pB2_df), "_") %>% sapply(., "[[", 2) 
pB2_df<-cbind(pB2_df,Sample=patients,Number)
head(pB2_df)
pB2_df$clusters<-rownames(pB2_df)
Allcell_clusterss<-merge(Allcell_clusters,pB2_df,by="clusters",all=T)
Allcell_clusterss<-Allcell_clusterss[!duplicated(Allcell_clusterss),]

############################
group_CAF<-Allcell_clusterss[,c("cell_name","clus")]
rownames(group_CAF)<-group_CAF$cell_name
group_CAF<-group_CAF[,2]
#save(group_CAF,file="single_dat_CAF_NMF_m6aGroup.Rdata")
#load("single_dat_CAF_NMF_m6aGroup.Rdata")
dim(single_dat_sce)
single_dat_sce<-AddMetaData(single_dat_sce,metadata=group_CAF, col.name=c("cell","CAF_m6aGroup"))
DimPlot(single_dat_sce,reduction = "tsne",group.by="CAF_m6aGroup",cols=mycol)
#####################
##########
single_dat_sce@meta.data$CAF_m6aGroup<-ifelse(is.na(single_dat_sce@meta.data$CAF_m6aGroup),"5",single_dat_sce@meta.data$CAF_m6aGroup)
DimPlot(single_dat_sce,reduction = "tsne",group.by="CAF_m6aGroup",cols=mycol)
DotPlot(single_dat_sce,features =m6agenes,group.by="CAF_m6aGroup", cols = c("green", "red")) + RotatedAxis()+ coord_flip()
####
single_dat_sce.markers <- FindAllMarkers(single_dat_sce, group.by="CAF_m6aGroup",
                                         only.pos = TRUE, 
                                         min.pct = 0.25, logfc.threshold = 0.25)

write.csv(single_dat_sce.markers,"single_dat_sce.markers_SMC_CAF_NMF_combined.csv")###
single_dat_sce.markers  %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
single_dat_sce.markers_top5 <- single_dat_sce.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

logFCfilter<-0
adjPvalFilter<-0.01
sig.markers=single_dat_sce.markers[(abs(as.numeric(as.vector(single_dat_sce.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(single_dat_sce.markers$p_val_adj))<adjPvalFilter),]

DoHeatmap(single_dat_sce,features = single_dat_sce.markers_top5$gene,group.by="CAF_m6aGroup")
ggsave("DoHeatmap_CAF_NMF.pdf",height =2.8,width=6.6)

DotPlot(single_dat_sce,features = single_dat_sce.markers_top5$gene,group.by="CAF_m6aGroup")+ RotatedAxis()
ggsave("DoDotPlot_CAF_NMF.pdf",height =2.91,width=6.6)
###"CBLL1","KIAA1429","METTL14","METTL3", "RBM15","RBM15B", "WTAP", "IGF2BP1", 
#"IGF2BP2", "IGF2BP3", "LRPPRC", "YTHDC1", "YTHDC2", "YTHDF1", 
#"YTHDF2","YTHDF3","ALKBH5","FTO"
#"ZC3H13","ELAVL1","FMR1", 
Featuresgenes<-c("HNRNPA2B1","WTAP", "HNRNPC")
VlnPlot(single_dat_sce,Featuresgenes,group.by="CAF_m6aGroup",pt.size = 0,ncol=3)
ggsave("VlnPlot_Featuresgenes_CAF_NMF.pdf",height =2.51,width=5.68)

###
###########################
library(dittoSeq)
a<-dittoPlot(single_dat_sce,"HNRNPA2B1",group.by="Sample") & NoLegend()
b<-dittoPlot(single_dat_sce,"WTAP",group.by="Sample") & NoLegend()
c<-dittoPlot(single_dat_sce,"HNRNPC",group.by="Sample") & NoLegend()
d<-dittoPlot(single_dat_sce,"cluster33",group.by="Sample") & NoLegend()
a+b+c
dittoHeatmap(subset(single_dat_sce, downsample = 4), m6agenes, use_raster = TRUE,
             annot.by = "CAF_m6aGroup",color.by=c("red","blue"))

dittoHeatmap(single_dat_sce, single_dat_sce.markers_top5[6:16,]$gene, 
             #use_raster = TRUE,
             scaled.to.max = TRUE,
             #annot.colors=ann_cols,
             #scale = "row",
             fontsize_number = 0.2,
             fontsize = 10,
             rownames="left",
             annot.by = c("CAF_m6aGroup"),color.by=c("red","blue"))
#####################################
pan_caf_file=as.matrix(single_dat_sce@assays$RNA@data)
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
pan_caf_file_proteins_for_heatmap <- AverageExpression(single_dat_sce, assays = "RNA",group.by = "CAF_m6aGroup", features = proteins,verbose = TRUE) %>% .$RNA
#pan_caf_file_proteins_for_heatmap <- pan_caf_file[match(proteins, row.names(pan_caf_file)),]
pan_caf_file_proteins_for_heatmap <- na.omit(pan_caf_file_proteins_for_heatmap)

pdf("pathway_CAF_NMF.pdf",height = 8.15,width=2.55)
pheatmap(pan_caf_file_proteins_for_heatmap, 
         scale = "row", 
         cluster_rows = T, 
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
pan_caf_epic_file_receptors_heatmap <- AverageExpression(single_dat_sce, assays = "RNA",
                                                         group.by = "CAF_m6aGroup",
                                                         features = cell_surface_markers ,verbose = TRUE) %>% .$RNA
pan_caf_epic_file_receptors_heatmap <- na.omit(pan_caf_epic_file_receptors_heatmap)
###
pdf("cell_surface_CAF_NMF.pdf",height = 4.94,width=3.02)
pheatmap(pan_caf_epic_file_receptors_heatmap[-2,],
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = T,  
         show_rownames = T, angle_col = c('45'),
         color=colorRampPalette(c("blue","white","red"))(20))
dev.off()
#ggsave("cell_surface_CAF_NMF.pdf",height = 4.94,width=3.02)
#TF analysis
TFgenes<-read_xlsx("TFgenes.xlsx",1)
TFgenes<-data.frame(read_excel("TFgenes.xlsx",1, col_names= TRUE, col_types= NULL, na=" "))##TFgenes
##


##########################################
caf.combined.resolution.02.markers <- FindAllMarkers(caf.combined.resolution.02, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.6)
caf.combined.resolution.02.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

caf.combined.resolution.02.markers.top20 <-   caf.combined.resolution.02.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(caf.combined.resolution.02, features = caf.combined.resolution.02.markers.top20$gene, label = F) + NoLegend()
#integrated_cafs_ids <- c("C1", "C2", "C3", "C4")
#names(integrated_cafs_ids) <- levels(caf.combined.resolution.02)
#caf.combined.resolution.02 <-RenameIdents(caf.combined.resolution.02, integrated_cafs_ids)
feature <- caf.combined.resolution.02.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
feature<-data.frame(feature)
feature<-feature$gene
colon4_dot_plot <- DotPlot(caf.combined.resolution.02, features = feature, cols = c("lightgrey", "red")) + RotatedAxis()
colon4_dot_plot
#################################
