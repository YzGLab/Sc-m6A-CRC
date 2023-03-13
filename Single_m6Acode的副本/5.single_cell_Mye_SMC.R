####myeloid for SMC not com
#input data
###
############
#library()
#devtools::install_github("cran/randomSurvivalForest")
#library("randomSurvivalForest")
library(dittoSeq)
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
#gmt=read.gmt(misfile) #   input from the MsiDB data.    
m6agenes<-c("CBLL1","KIAA1429","METTL14","METTL3", "RBM15","RBM15B", "WTAP",
            "ZC3H13","ELAVL1","FMR1", "HNRNPA2B1", "HNRNPC", "IGF2BP1", 
            "IGF2BP2", "IGF2BP3", "LRPPRC", "YTHDC1", "YTHDC2", "YTHDF1", 
            "YTHDF2","YTHDF3","ALKBH5","FTO")
mycol <- c("#223D6C","#D20A13","#FFD121","#088247",
           "#58CDD9","#5D90BA","#431A3D","#11AA4D","#91612D","#6E568C","#7A142C",
           "#E0367A","#D8D155","#64495D","#7CC767")
#load("Pan_myeloid_cells_CRC.Rdata")
#single_dat<-readRDS("O_data/Myeloidssingle_dat_GSE132465.rds")
single_dat<-readRDS("Myeloidssingle_dat_GSE132465.rds")
rowGenenames<-readRDS("rowGenenames.rds")
rownames(single_dat)<-rowGenenames
single_dat[1:5,1:5]
GSE132465_phe <- readRDS('GSE132465_phe.rds')#####

##
single_dat<-single_dat[rownames(single_dat) %in% c(m6agenes),]
topn <- "m6A"
ranks <- 4##
i<-Cellsubtype[1]
phe<-GSE132465_phe
Celltype<-unique(phe$Cell_type)[3]
Cellsubtype<-unique(phe$Cell_subtype[phe$Cell_type=="Myeloids"])
Cellphe<-phe[phe$Cell_type=="Myeloids",]
# NMF
for (i in Cellsubtype){
  if (!dir.exists("nmfSingle_colon_myeloid")){
    dir.create("./nmfSingle_colon_myeloid")
  }
  if (!dir.exists(paste0("nmfSingle_colon_myeloid/", i))){
    dir.create(paste0("nmfSingle_colon_myeloid/", i))
  }
  print(i)
  #nmfdat <- single_dat[, str_detect(colnames(single_dat), i)]
  nmfdat <- single_dat[, colnames(single_dat) %in% Cellphe$Index[Cellphe$Cell_subtype==i]]
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
  #tiff(paste0("./nmfSingle_colon_myeloid/", i, "/consensus_", topn, ".tiff"), 
  #  width = 6*480, height = 4*480, res = 100)
  #consensusmap(res_4)
  #dev.off()
  #tiff(paste0("./nmfSingle_colon_myeloid/", i, "/basicmap_", topn, ".tiff"), 
  #  width = 8 * 480, height = 4*480, res = 100)
  #basismap(res_4)
  #dev.off()
  #tiff(paste0("./nmfSingle_colon_myeloid/", i, "/coefmap_", topn, ".tiff"), 
  # width = 8 * 480, height = 4*480, res = 100)
  # coefmap(res_4)
  # dev.off()
  colnames(signature) <- paste0(i, "_", 1:ranks)
  signature <- as.data.frame(signature)
  # 
  write.table(signature, paste0("./nmfSingle_colon_myeloid/", i, "/singature", topn, ".txt"),
              sep = "\t")
  write.table(sicluster, paste0("./nmfSingle_colon_myeloid/", i, "/sicluster", topn, ".txt"),
              sep = "\t")
  #saveRDS(res_4, file = paste0("./nmfSingle_colon_myeloid/", i, "/res", ranks, "_", topn, ".rds"))
  print(paste0("NMF for the", i, "is Done!"))
}
#####
#NMF program???�???ͼ
library(dplyr)
# ??ȡprogram
topRank <- 5
programG <- list()
for (i in 1:length(Cellsubtype)){
  filedir <- paste0("./nmfSingle_colon_myeloid/", Cellsubtype[i], "/singature", topn, ".txt")
  geneloading <- read.delim2(filedir, header = T, sep = "\t",check.names = F)
  geneloading$maxC <- apply(geneloading, 1, which.max) %>% 
    paste0(Cellsubtype[i], "_", .)
  topgenelist <- rownames_to_column(geneloading, var = "gene")  %>% 
    pivot_longer(., cols = names(geneloading)[1:4], 
                 names_to = "program", values_to = "loading")
  topgenelist <- dplyr::filter(topgenelist, maxC == program) %>% 
    group_by(maxC) %>% top_n(n = topRank, wt = loading)
  topgenelist <- split(topgenelist$gene, topgenelist$maxC)
  programG <- c(programG, topgenelist)
}
#??fbϸ?????д??֣?????????ѡ??AUCell
#single_dat<-single_dat[rownames(single_dat) %in% c(m6agenes),]

single_dat<-single_dat[rowSums(single_dat!= 0)>=1,colSums(single_dat != 0) >= 1]
cells_rankings <- AUCell_buildRankings(as.matrix(single_dat), nCores = 10, plotStats=TRUE)
#cells_AUC <- AUCell_calcAUC(geneSets = programG, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
cells_AUC <- AUCell_calcAUC(geneSets = programG, cells_rankings, aucMaxRank=3)
programAUC <- getAUC(cells_AUC)
programAUC[1:5,1:5]
###
#save(single_dat,programAUC,cells_AUC,cells_rankings,clus,file="programAUC.Rdata")
#Cellphe_NMF<-Cellphe[Cellphe$Index %in% colnames(df),]
#head(Cellphe_NMF)
#Cellphe_NMF<-Cellphe_NMF[match(colnames(df),row.names(Cellphe_NMF)),]
#head(Cellphe_NMF)
#pan_caf_file_receptors <- pan_caf_file[match(cell_surface_markers, row.names(pan_caf_file)),]

clusterK = 4
M <- cor(t(programAUC), method = "pearson")

#??corrplot??ͼ?????ڱ?????cluster?ĺ?ɫ????
pdf("scNMF_myeloid.pdf")

corrplot(M, 
         method = "color", #????ɫչʾ????ϵ?????????Ը?Ϊ"circle" (default), "square", "ellipse", "number", "pie", "shade"
         order = "hclust", 
         hclust.method = "ward.D2", 
         addrect = clusterK, #????ɫ????
         tl.pos = "n", #"n" means don't add textlabel
         col = rev(brewer.pal(n = 8, name = "RdBu"))) # ?????brewer.pal?�???????ɫ????

dev.off()

library(ggplot2)
library(ggpubr)
require(ggplotify)
library(scales)
##??ȡprogram
cororder <- corrMatOrder(M, order = "hclust", hclust.method = "ward.D2")
M.hc <- M[cororder, cororder]
tree <- hclust(as.dist(1 - M.hc), method = "ward.D2")
clus <- cutree(tree, clusterK)
table(clus)
head(clus)
clus
clus<-data.frame(clus)
cellsubtypes2 <- str_split(rownames(clus), "_") %>% sapply(., "[[", 1) 
Number<-str_split(rownames(clus), "_") %>% sapply(., "[[", 2) 
clus<-cbind(clus,cellsubtype=cellsubtypes2,Number)
head(clus)
clus$Number<-as.numeric(clus$Number)
clus$clus<-as.factor(clus$clus)
###extract the data of cluster for all cells in NMF
cell_cluster <- lapply(Cellsubtype, function(z){
  filedir <- paste0("./nmfSingle_colon_myeloid/", z, "/sicluster", topn, ".txt")
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
head(clus)

######
clus$clusters<-rownames(clus)
clus$clusters<-gsub("\\."," ",clus$clusters)
Allcell_clusterss<-merge(Allcell_clusters,clus,by="clusters",all=T)
Allcell_clusterss<-Allcell_clusterss[!duplicated(Allcell_clusterss),]
##########################################
head(Allcell_clusterss)
Allcell_clusterss$NMFcluster<-Allcell_clusterss$clus
#####################################################
save(Allcell_clusterss,file="myeloids_NNFclusters.Rdata")
################################### ??ѡ ????????��Դ??ϸ??Ⱥ??
patients <- str_split(colnames(single_dat), "_") %>% sapply(., "[[", 1) %>% unique()
#single_dat_Tumor<-single_dat[,str_split(colnames(single_dat), "_") %>% sapply(., "[[", 1) %in% patients[1:23]]
########################????seurat object
single_dat<-readRDS("Myeloidssingle_dat_GSE132465.rds")
rowGenenames<-readRDS("rowGenenames.rds")
rownames(single_dat)<-rowGenenames
single_dat[1:5,1:5]
sce <- CreateSeuratObject(counts =single_dat, project = "GSE132465_myeloid_colon")
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
ElbowPlot(sce, ndims = 50, reduction = "pca")
sce <- FindNeighbors(sce, dims = 1:20)
##########
sce<- FindClusters(sce, resolution = 0.2)
##################RunUMAP
sce<- RunUMAP(sce, dims = 1:20)
DimPlot(sce, reduction = "umap",label=T,cols=mycol)
#RunTSNE
sce<- RunTSNE(object = sce, dims = 1:20, do.fast = TRUE)
DimPlot(sce,reduction = "tsne",label=T,cols=mycol)
DimPlot(sce,reduction = "tsne",label=T,cols=mycol,group.by = "Cell_subtype")
##
phe<-GSE132465_phe
sce=AddMetaData(object=sce, metadata=phe, col.name=colnames(phe))##
#
Idents(sce)<-sce$Cell_subtype

#FeaturePlot(sce,reduction = "tsne",features = m6agenes[1:18],ncol = 6) & viridis::scale_color_viridis(option="H")
##########
FeaturePlot(sce,reduction = "tsne",features = m6agenes[1:6]) & viridis::scale_color_viridis(option="H")
ggsave("m6Agenes_featureplot1_6.pdf",height =5.73 ,width =6.02 )
FeaturePlot(sce,reduction = "tsne",features = m6agenes[7:12]) & viridis::scale_color_viridis(option="H")
ggsave("m6Agenes_featureplot7_12.pdf",height =5.73 ,width =6.02 )
FeaturePlot(sce,reduction = "tsne",features = m6agenes[13:18]) & viridis::scale_color_viridis(option="H")
ggsave("m6Agenes_featureplot13_18.pdf",height =5.73 ,width =6.02 )
FeaturePlot(sce,reduction = "tsne",features = m6agenes[19:23]) & viridis::scale_color_viridis(option="H")
ggsave("m6Agenes_featureplot19_23.pdf",height =5.73 ,width =6.02 )
####################
table(sce$Cell_subtype)
sce<-sce[,sce$Cell_subtype %in% c('Pro-inflammatory','Proliferating' ,'SPP1+')]
sce<-sce[,sce$Class=="Tumor"]
##
DimPlot(sce,reduction = "tsne",cols=mycol)
DimPlot(sce,reduction = "tsne",group.by="Class",label=T,cols=mycol)
ggsave("myeloid_Class.SMC.pdf")
DimPlot(sce,reduction = "tsne",group.by="Sample")
ggsave("myeloid_Sample.SMC.pdf")
DimPlot(sce,reduction = "tsne",group.by="Cell_subtype",label=T,cols=mycol)
ggsave("myeloidType.SMC.pdf")
DimPlot(sce,reduction = "tsne",cols=mycol)
##add the phe data
newphe<-data.frame(read_excel("clinicalphe.xlsx",1, col_names= TRUE, col_types= NULL, na=" "))##chip data
newphe<-merge(phe,newphe,by="Sample",all=T)
rownames(newphe)<-newphe$Index
sce<-AddMetaData(sce,metadata=newphe, col.name=colnames(newphe))######
colnames(sce@meta.data)
#################
DimPlot(sce,reduction = "tsne",cols=mycol,group.by="MSI")
load("myeloids_NNFclusters.Rdata")
rownames(Allcell_clusterss)<-Allcell_clusterss$cell_name
sce<-AddMetaData(sce,metadata=Allcell_clusterss, col.name=colnames(Allcell_clusterss))######
sce@meta.data$clusters<-ifelse(is.na(sce@meta.data$clusters),0,sce@meta.data$clusters)
sce@meta.data$NMFcluster<-ifelse(is.na(sce@meta.data$NMFcluster),5,sce@meta.data$NMFcluster)
sce@meta.data$NMFcluster<-as.factor(sce@meta.data$NMFcluster)
table(sce@meta.data$NMFcluster,sce@meta.data$Cell_subtype)
#####
####################################
colnames(sce@meta.data)
dittoBarPlot(sce, "NMFcluster",group.by="Class",color.panel=mycol)+coord_flip()
dittoBarPlot(sce, "NMFcluster",group.by="MSI",color.panel=mycol)+coord_flip()
####

#####################
inflammatoryscore<-c("IL1B","IL6","S100A8","S100A9")
anti_inflammatory<-c('CD163','SEPP1','SELENOP', 'APOE','MAF')
inflammatorylist<-list()
inflammatorylist[["inflammatory"]]<-inflammatoryscore
inflammatorylist[["anti_inflammatory"]]<-anti_inflammatory
sce<-AddModuleScore(sce,features=inflammatorylist,col.names="inflam")
colnames(sce@meta.data)
FeatureScatter(sce,feature1="Cluster1",feature2 ="Cluster2",group.by = "NMFcluster" )
VlnPlot(sce,features=c("Cluster2","Cluster1"),group.by="NMFcluster",pt.size = 0)
RidgePlot(sce,features=c("Cluster2","Cluster1"),group.by="NMFcluster")

#####################
inflammatoryscores<-cbind(sce$NMFcluster,sce$Cluster1,sce$Cluster2)
inflammatoryscoresmeans<-aggregate(inflammatoryscores[,c(2:3)],list(inflammatoryscores[,1]),mean)######
inflammatoryscoresmeans<-t(inflammatoryscoresmeans)
colnames(inflammatoryscoresmeans)<-inflammatoryscoresmeans[1,]
inflammatoryscoresmeans<-inflammatoryscoresmeans[-1,,drop=F]
pheatmap(inflammatoryscoresmeans,
         scale="row")

##############
library(ComplexHeatmap)
library(circlize)
library(janitor)
meta.datas<-sce@meta.data
as<-tabyl(meta.datas, NMFcluster)
as<-as[-5,]
colnames(as)[1]<-"Subgroup"
AveExpression <- AverageExpression(sce, assays = "RNA",group.by = "NMFcluster", features =m6agenes,verbose = TRUE) %>% .$RNA
AveExpressions<-t(AveExpression)
AveExpressions<-scale(AveExpressions)
rownames(AveExpressions)
col_fun = colorRamp2(c(0, 0.3, 0.8), c("blue", "white", "red"))
rowAnn<-HeatmapAnnotation(Percent = as$percent, No.cells = anno_barplot(as$n),col=list(Percent = col_fun))
##################################
AveExpressions<-AveExpressions[,-13]
pdf("TAM_cell_m6aGenes_Heatmap_m6aGroup.pdf",width = 5.42,height=4.83)
Heatmap(t(AveExpressions[-5,]),name = "AveE.",
        cluster_rows = TRUE,
        cluster_columns = F,
        #row_names_size=6,
        #row_split = as$splitgroups ,
        row_names_side =  "left",
        top_annotation = rowAnn
)
dev.off()###

###################
data_plotC <- table(sce@meta.data$Cell_subtype, sce@meta.data$NMFcluster) %>% melt()
colnames(data_plotC) <- c("Cell_subtype", "NMFcluster","Number")
data_plotC$NMFcluster<-as.factor(data_plotC$NMFcluster)
pC1 <- ggplot(data = data_plotC, aes(x = Cell_subtype, y = Number, fill =NMFcluster )) +
  geom_bar(stat = "identity", width=0.8,aes(group=NMFcluster),position="stack")+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Average number")+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  #theme(axis.title.y=element_blank(),
  #  axis.ticks.y=element_blank(),
  # axis.text.y=element_blank()  )+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6)) + 
  coord_flip()
 # scale_x_discrete(position = "top") +
 # scale_y_reverse()+theme(legend.position="none")
#scale_y_discrete(position = "left") 
#
pC1
pC2 <- ggplot(data = data_plotC, aes(x = Cell_subtype, y = Number, fill =NMFcluster )) +
  geom_bar(stat = "identity", width=0.8,aes(group=NMFcluster),position="fill")+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  scale_y_continuous(labels = percent)+  ####用来将y轴移动位�?
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6)) + 
  coord_flip()
 #+theme(legend.position="top")
# + scale_x_discrete(position = "top") 
#让横轴上的标签倾斜45�?
pC2
ggsave("TAM_cell_subtype_NMF2.pdf")
library(ggpubr)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(scales)
pC <- pC1 + pC2 + plot_layout(ncol = 1, widths = c(1,1),guides = 'collect')

pC
ggsave("TAM_cell_subtype_NMF.pdf")
###############
library(dittoSeq)
#levels(Idents(sce))
Idents(sce)<-sce$NMFcluster
#sce = sce[, Idents(sce) %in% 
#   c( "FCGR3A+ Mono", "CD14+ Mono"  )] # CD16 
dittoDimPlot(sce,reduction = "tsne","Cell_subtype",color.panel = mycol)
ggsave("myeloid_Cell_subtype.SMC.pdf")
dittoDimPlot(sce,reduction = "tsne","NMFcluster")
ggsave("myeloid_NMFcluster.SMC.pdf")
dittoDimPlot(sce,reduction = "tsne","Class",color.panel = mycol)
ggsave("myeloid_Class.SMC.pdf")
dittoDimPlot(sce,reduction = "tsne","Cell_subtype",color.panel = mycol)
ggsave("myeloidType.SMC2.pdf")
dittoDimPlot(sce,reduction = "tsne","Sample",split.by="Cell_subtype")
ggsave("sample.myeloid_subtyp_smc2.pdf")
FeaturePlot(sce,reduction = "tsne",features = "WTAP")+ viridis::scale_color_viridis(option="H")
########################
table(sce@meta.data$Cell_subtype)##
## do table for 
####
########################
##
saveRDS(sce,file="myeloid_sce.rds")
sce<-readRDS("myeloid_sce.rds")
colnames(sce@meta.data)
#################################################################
sce$NMFcluster<-as.factor(sce$NMFcluster)
Idents(sce)<-sce$NMFcluster
sce<-sce[,sce$Class=="Tumor"]
sce<-sce[,sce$Cell_subtype %in% c('Pro-inflammatory','Proliferating' ,'SPP1+')]

sce.markers <- FindAllMarkers(sce, 
                              only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
sce.markers<-sce.markers[sce.markers$p_val_adj<0.05,]
##
save(sce.markers,file="sce.markers_mye_SMC_TAM.Rdata")
##############################################################
sce.markers_top20 <- sce.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)
sigmethy<-split(sce.markers_top20$gene,sce.markers_top20$cluster)
sigmethy<-sigmethy[c(1,3,4)]
sce<- AddModuleScore(object = sce, features = sigmethy, name = "Cluster") #AddModuleScore计算得分
colnames(sce@meta.data)
####
sce@meta.data$NMFcluster<-as.factor(sce@meta.data$NMFcluster)
FeaturePlot(sce,features ="Cluster1",reduction = "tsne")+ viridis::scale_color_viridis(option="H")
ggsave("FeaturePlotNMF_C1_TAM_plot.pdf")
FeaturePlot(sce,features = "Cluster2",reduction = "tsne")+ viridis::scale_color_viridis(option="H")
#FeaturePlot(sce,features = "Cluster3",reduction = "tsne")+ viridis::scale_color_viridis(option="H")
ggsave("FeaturePlotNMF_C3_TAM_plot.pdf")
###
sce.markers_top10 <- sce.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
sce.markers_top5<- sce.markers_top10[sce.markers_top10$cluster %in% c("1","3"),]
DotPlot(sce,features = sce.markers_top5$gene,group.by ="NMFcluster") + RotatedAxis()
ggsave("FeaturePlotNMF_C3_TAM_plot.pdf")
sce.markers_mye_list<-split(sce.markers$gene,sce.markers$cluster)
####################################################################################
####
write.csv(sce.markers,"sce.markers_SMC_TAM_tumor_NMF.csv")###
DoHeatmap(subset(sce, downsample = 50),features = sce.markers_top5[sce.markers_top5$gene %in% m6agenes,]$gene,
          group.by ="NMFcluster") 
ggsave("TAM_Heatmap_m6Agenes.pdf")
DoHeatmap(subset(sce, downsample = 100),features = sce.markers_top5$gene,group.by ="NMFcluster") 
DotPlot(sce,features = sce.markers_top5$gene,group.by ="NMFcluster") + RotatedAxis()
ggsave("dotplot_TAM_SMC.pdf")
#FeaturePlot(sce,reduction = "tsne",features = "WTAP")+ viridis::scale_color_viridis(option="H")
FeaturePlot(sce,features = sce.markers_top5$gene[2],reduction = "tsne")+ viridis::scale_color_viridis(option="H")
VlnPlot(sce,features=sce.markers_top5$gene[2],group.by = "NMFcluster",pt.size = 0)
## WTAP HNRNPC HNRNPA2B1 YTHDC1 YTHDF3 FTH1
VlnPlot(sce,features="HNRNPC",group.by = "NMFcluster",pt.size = 0)

########################################################################
##构建 score list
names(sigmethy)
sce<-sce[,sce$Cell_subtype %in% c('Pro-inflammatory','Proliferating' ,'SPP1+')]
sce$cellsubtype<-as.numeric(factor(sce$Cell_subtype))
table(sce$cellsubtype,sce$Cell_subtype)
sce.markers_naturegenetic <- FindAllMarkers(sce, group.by="cellsubtype",
                              only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
sce.markers_naturegenetic<-sce.markers_naturegenetic[sce.markers_naturegenetic$p_val_adj<0.05,]
sce.markers_naturegenetic_top5 <- sce.markers_naturegenetic %>% 
  group_by(cluster) %>% 
  top_n(n =20, wt = avg_log2FC)
sce.markers_naturegenetic_top5<-split(sce.markers_naturegenetic_top5$gene,sce.markers_naturegenetic_top5$cluster)
sce.markers_naturegenetic_top5<-sce.markers_naturegenetic_top5[c(1:3)]
names(sce.markers_naturegenetic_top5)<-c("Pro-inflammatory","Proliferating","SPP1+")
#####
names(sigmethy)<-c("C1","C3","C4")
###
c1agenes<-c("CD74","H2-AA","H2-AB1","H2-EB1","PF4","CD63","C1QB","C1QC","C1QC","C1QA","PLTP","CX3CR1","CD81")
c1agenes<-list(c1agenes)
names(c1agenes)<-"c1agenes"
####
#################read m1 m2
m2m2genes<-read.csv("immune genes markers list_Macrophage.csv",header = T)
m2m2genes<-m2m2genes[m2m2genes$type %in% rownames(sce) &m2m2genes$iord !="Tregs genes" ,]
m2m2genes<-split(m2m2genes$type,m2m2genes$iord)
names(m2m2genes)<-c("m1","m2")####
######
scoretype<-list()
for(set in names(sigmethy)){
  scoretype[[set]]<-sigmethy[[set]]
}
for(set in names(m2m2genes)){
  scoretype[[set]]<-m2m2genes[[set]]
}
for(set in names(c1agenes)){
  scoretype[[set]]<-c1agenes[[set]]
}
#
for(set in names(sce.markers_naturegenetic_top5)){
  scoretype[[set]]<-sce.markers_naturegenetic_top5[[set]]
}
###################
names(scoretype)
scoretype[["C4"]]<-c()
################################
names(scoretype)
save(scoretype,file="8typesofscoretype.mac.Rdata")
###########
colnames(sce@meta.data)
sce<- AddModuleScore(object = sce, features = scoretype, name = names(scoretype)) #AddModuleScore计算得分
#colnames(sce@meta.data)[c(25:32)]<- names(scoretype)
colnames(sce@meta.data)[c(22:29)]<- names(scoretype)
########################################################
Idents(sce)<-sce$NMFcluster
table(Idents(sce) )
####
cg=sce@assays$RNA@var.features

##########
RidgePlot(sce,features=names(scoretype),group.by="NMFcluster")
FeaturePlot(sce,features=names(scoretype),reduction = "tsne") & viridis::scale_color_viridis(option="H")
ggsave("scorelistformac_featureplot.pdf")
###
library(ggcorrplot)
library(ggthemes)
#??corrplot??ͼ?????ڱ?????cluster?ĺ?ɫ????
results<-cor(sce@meta.data[c(22:29)])
pdf("cor_c1,c3_myeloid.pdf")
ggcorrplot(results, method = "circle",
           hc.order = TRUE, hc.method = "ward.D",
           outline.col = "white", ggtheme = theme_bw())
dev.off()
??ggcorrplot
##################################################################
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

##how to do next step? 
###

#sce = sce[, Idents(sce) %in% 
# c( "FCGR3A+ Mono", "CD14+ Mono"  )] # CD16 
##
colnames(sce@meta.data)
DoHeatmap(subset(sce, downsample = 100),features = sce.markers_top5$gene,group.by ="NMFcluster")
####
data_plotC <- table(sce@meta.data$NMFcluster, sce@meta.data$Cell_subtype) %>% melt()
colnames(data_plotC) <- c("Cell_subtype", "NMFcluster","Number")
data_plotC$NMFcluster<-as.factor(data_plotC$NMFcluster)
pC1 <- ggplot(data = data_plotC, aes(x = Cell_subtype, y = Number, fill = NMFcluster)) +
  geom_bar(stat = "identity", width=0.8,aes(group=NMFcluster),position="stack")+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Average number")+
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6)) 
pC1
#ggsave("cellnumberofNMF_TAM.pdf")
pC2 <- ggplot(data = data_plotC, aes(x = Cell_subtype, y = Number, fill = NMFcluster)) +
  geom_bar(stat = "identity", width=0.8,aes(group=NMFcluster),position="fill")+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  scale_y_continuous(labels = percent)+  ####??��??y???ƶ?λ??
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6))     

library(ggpubr)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(scales)
pC <- pC1 + pC2 + plot_layout(ncol = 2, widths = c(1,1),guides = 'collect')

pC
ggsave("myeloid_subtype_NMF.pdf")
####metabolism for mac
metabolism<-read.gmt("Immunemetabolism_geneset.txt")
metabolism<-metabolism[metabolism$gene!="",]
metabolism<-split(metabolism$gene,metabolism$term)#########
rt<-t(as.matrix(sce@assays$RNA@data))
rt<-as.matrix(t(rt))
rt[1:5,1:5]
library(GSVA)
gsva_es <- gsva(rt,metabolism,method="ssgsea",abs.ranking=F,kcdf="Gaussian",parallel.sz=40)###Array data,ssgsea.norm=TRUE
#gsva_es <- gsva(norm.expMat,gmt.list,method="ssgsea",abs.ranking=F,kcdf="Poisson",parallel.sz=10)#RNA-Seq 
save(gsva_es,file="TAM_metabolism_GSVA.Rdata")
load("TAM_metabolism_GSVA.Rdata")
gsva_es<-t(gsva_es)
gsva_es<-data.frame(gsva_es)
sce <- AddMetaData(sce, gsva_es, col.name =colnames(gsva_es))
####
library(ggpubr)
#########################
gsva_es<-cbind(gsva_es,group= as.character(Idents(sce)))
colnames(gsva_es)<-gsub("\\."," ",colnames(gsva_es))
colnames(gsva_es)
#####
head(wilcox)
group= as.character(Idents(sce))
data<-gsva_es
outTab<-c()
for (i in colnames(gsva_es)){
  geneName=i
  data[,i]<-as.numeric(data[,i])
  rt=rbind(expression=data[,i],group=group)
  rt=as.matrix(t(rt))
  rt<-data.frame(rt)
  #Means=rowMeans(data[,i])
  #if(is.numeric(Means) & Means>0){
    kTest<-kruskal.test(expression ~ group, data=rt)
    pvalue=kTest$p.value
    outTab=rbind(outTab,cbind(gene=i,pValue=pvalue))
    message(i,"is done now!!!")
  #}
}#
outTab<-data.frame(outTab)
outTab<-outTab[outTab$pValue<0.0001,]

#aucs<-data.frame(aucs)
#group_by(aucs, group) %>% summarize_each(funs(mean), colnames(aucs)[1:81])
gsva_ess<-apply(gsva_es,2,as.numeric)
rownames(gsva_ess)<-rownames(gsva_es)
gsva_esmenas<-aggregate(gsva_ess[,c(1:113)],list(gsva_ess[,114]),mean)#########根据最后确定的TF数量进�?�构�? heatmap
gsva_esmenas<-t(gsva_esmenas)
colnames(gsva_esmenas)<-gsva_esmenas[1,]
gsva_esmenas<-gsva_esmenas[-1,]
gsva_esmenas<-gsva_esmenas[rownames(gsva_esmenas) %in% outTab$gene,]
######
pdf("Metabolism_NMF_TAM.pdf",height = 6.74,width=6.04)
pheatmap(gsva_esmenas,
         scale = "row",
         show_colnames =T,
         show_rownames = T,
         cluster_cols = F,
         angle_col = c('0'),
         color=colorRampPalette(c("blue","white","red"))(20)
)
dev.off()

#####checkpoints
checkpoints<-c("CD274", "CTLA4","LAG3", "TIM3", "TNFRSF9","TIGIT","CD226","CD7",
               "CD4","CD8A","CD8B","FOXP3","IL2",
               "CXCL8","PDCD1","HAVCR2","GZMB","PRF1","NLG1","IFNG","TNFRSF18","TNFRSF4")
checkpoints_T_heatmap <- AverageExpression(sce, assays = "RNA",
                                           group.by = "NMFcluster",
                                           features = checkpoints ,verbose = TRUE) %>% .$RNA
checkpoints_T_heatmap <- na.omit(checkpoints_T_heatmap)
colnames(sce@meta.data)
#annotations<-data.frame(sce@meta.data$clusters,sce@meta.data$NMFcluster)
#annotations<-annotations[!duplicated(annotations$sce.meta.data.clusters),]
#rownames(annotations)<-annotations$sce.meta.data.clusters
#annotations<-annotations[order(annotations$sce.meta.data.NMFcluster),]
#checkpoints_T_heatmap<-checkpoints_T_heatmap[,annotations$sce.meta.data.clusters]
pdf("cell_checkpoints_NMF_macrophage.pdf",height = 4.44,width=3.02)
pheatmap(checkpoints_T_heatmap,
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = F,  
         show_rownames = T, 
         angle_col = c('45'),
         #annotation_col = annotations,annotation_legend = FALSE,
         color=colorRampPalette(c("blue","white","red"))(20))
dev.off()
###c1Q macrophage
c1agenes<-c("CD74","H2-AA","H2-AB1","H2-EB1","PF4","CD63","C1QB","C1QC","C1QC","C1QA","PLTP","CX3CR1","CD81")
c1agenes<-list(c1agenes)
sce<-AddModuleScore(sce,features =c1agenes,names="c1agenes")
###
VlnPlot(sce,features = "Cluster1",group.by = "NMFcluster",pt.size = 0)
FeaturePlot(sce,features = "Cluster1",reduction = "tsne") & viridis::scale_color_viridis(option="H")
ggsave("c1qliksTAM_methylation.pdf")
data_plot_mean <- cbind(cell=sce@meta.data$Index,NMFcluster=sce@meta.data$NMFcluster,c1agenes=sce$Cluster1) %>% data.frame() 
data_plot_mean <- melt(data_plot_mean,
                       id.vars = c("cell","NMFcluster"),
                       variable.name ="Signature",
                       value.name = "Signature_score")
data_plot_mean$Signature_score<-as.numeric(data_plot_mean$Signature_score)
#data_plot_mean <- data_frame(sce@meta.data$NMFcluster, sce@meta.data$Cluster2)
#colnames(data_plot_mean) <- c( "NMFcluster","gene")
ggplot(data = data_plot_mean, aes(x = NMFcluster, y = Signature_score, fill = NMFcluster)) +
  #geom_bar(stat = "identity", width=0.8,aes(group=NMFcluster),position="fill")+
  geom_boxplot()+
  #facet_wrap(data_plot_mean$Signature)+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Expression")+
  #scale_y_continuous(labels = percent)+  ####
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6)) +
  stat_compare_means(method = 'wilcox.test',comparisons = list(c("1", "2"),
                                                                c("2","3"),
                                                                c("3","4"),
                                                                c("4","5")))  + NoLegend()  #?, kruskal.test
ggsave("c1qliksTAM_methylationwith NMF.pdf")
#col_fun = colorRamp2(c(0, 0.3, 0.8), c("blue", "white", "red"))
#rowAnn<-rowAnnotation(Percent = as$percent, No.cells = anno_barplot(as$n),col=list(Percent = col_fun))
#####
c1agenes<-c("CD74","H2-AA","H2-AB1","H2-EB1","PF4","CD63","C1QB","C1QC","C1QC","C1QA","PLTP","CX3CR1","CD81")
c1agenes_heatmap <- AverageExpression(sce, assays = "RNA",
                                           group.by = "NMFcluster",
                                           features = c1agenes ,verbose = TRUE) %>% .$RNA
c1agenes_heatmap <- na.omit(c1agenes_heatmap)
colnames(sce@meta.data)
#annotations<-data.frame(sce@meta.data$clusters,sce@meta.data$NMFcluster)
#annotations<-annotations[!duplicated(annotations$sce.meta.data.clusters),]
#rownames(annotations)<-annotations$sce.meta.data.clusters
#annotations<-annotations[order(annotations$sce.meta.data.NMFcluster),]
#checkpoints_T_heatmap<-checkpoints_T_heatmap[,annotations$sce.meta.data.clusters]
pdf("cell_c1q_CancerCell_NMF_macrophage.pdf",height = 3.11,width=3.04)
pheatmap(c1agenes_heatmap,
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = F,  
         show_rownames = T, angle_col = c('45'),
         #annotation_col = annotations,annotation_legend = FALSE,
         color=colorRampPalette(c("blue","white","red"))(20))
dev.off()

####AUC 计算 通路激活分�?
load("hallmark.gsva.symbol.RData")
#####
cells_rankings <- AUCell_buildRankings(as.matrix(single_dat), nCores = 10, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(geneSets = gs, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
#cells_AUC <- AUCell_calcAUC(geneSets = gs, cells_rankings, aucMaxRank=5)
programAUC <- getAUC(cells_AUC)
programAUC[1:5,1:5]
####################
library(ggpubr)
library(pheatmap)
aucs<-getAUC(cells_AUC)
aucs<-t(aucs)
aucs<-aucs[rownames(aucs) %in% colnames(sce),]
aucs<-cbind(aucs,group= as.character(sce$NMFcluster))
#aucs<-data.frame(aucs)
dim(aucs)
#group_by(aucs, group) %>% summarize_each(funs(mean), colnames(aucs)[1:50])
aucss<-apply(aucs,2,as.numeric)
rownames(aucss)<-rownames(aucs)
aucsmenas<-aggregate(aucss[,c(1:50)],list(aucss[,51]),mean)#########根据最后确定的TF数量进�?�构�? heatmap
aucsmenas<-t(aucsmenas)
colnames(aucsmenas)<-c("C1","C2","C3","C4","C5")
aucsmenas<-aucsmenas[-1,]
#aucsmenas<-aucsmenas[!rownames(aucsmenas) %like% "extended",]
library(ComplexHeatmap)
Heatmap(aucsmenas)
rownames(aucsmenas)<-gsub("HALLMARK_","",rownames(aucsmenas))
pdf("HALLMARK_NMF_TAM.pdf",height = 6.74,width=4.54)
pheatmap(aucsmenas,
         scale = "row",
         show_colnames =T,
         show_rownames = T,
         cluster_cols = F,
         fontsize = 8,
         angle_col = c('45'),
         color=colorRampPalette(c("blue","white","red"))(20)
)
dev.off()
#################read m1 m2
m2m2genes<-read.csv("immune genes markers list_Macrophage.csv",header = T)
m2m2genes<-m2m2genes[m2m2genes$type %in% rownames(sce) &m2m2genes$iord !="Tregs genes" ,]
m2m2genes_heatmap <- AverageExpression(sce, assays = "RNA",
                                      group.by = "NMFcluster",
                                      features = m2m2genes$type ,verbose = TRUE) %>% .$RNA
m2m2genes_heatmap<- na.omit(m2m2genes_heatmap)
#m2m2genes_heatmap<-m2m2genes_heatmap[m2m2genes_heatmap!=0]
colnames(sce@meta.data)
pdf("cell_m2m2genes_NMF_macrophage.pdf",height = 3.11,width=3.04)
pheatmap(m2m2genes_heatmap,
         scale = "row", 
         cluster_rows =F, 
         cluster_cols = F,  
         show_rownames = T, angle_col = c('45'),
         #annotation_col = annotations,annotation_legend = FALSE,
         color=colorRampPalette(c("blue","white","red"))(20))
dev.off()
##score
m2m2genes<-split(m2m2genes$type,m2m2genes$iord)
sce<-AddModuleScore(sce,features = m2m2genes,colnames=c("M1","M2"))
colnames(sce@meta.data)
FeatureScatter(sce,feature1="Cluster1",feature2 ="Cluster2",group.by = "NMFcluster" )
VlnPlot(sce,features=c("Cluster2","Cluster1"),pt.size = 0)

data_plot_mean <- cbind(cell=rownames(rt),NMFcluster=sce@meta.data$NMFcluster,M1=sce$Cluster1,M2=sce$Cluster2) %>% data.frame() 
data_plot_mean <- melt(data_plot_mean,
                       id.vars = c("cell","NMFcluster"),
                       variable.name ="Signature",
                       value.name = "Signature_score")
data_plot_mean$Signature_score<-as.numeric(data_plot_mean$Signature_score)
#data_plot_mean <- data_frame(sce@meta.data$NMFcluster, sce@meta.data$Cluster2)
#colnames(data_plot_mean) <- c( "NMFcluster","gene")
ggplot(data = data_plot_mean, aes(x = NMFcluster, y = Signature_score, fill = Signature)) +
  #geom_bar(stat = "identity", width=0.8,aes(group=NMFcluster),position="fill")+
  geom_boxplot()+
  facet_wrap(data_plot_mean$Signature)+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Expression")+
  #scale_y_continuous(labels = percent)+  ####??��??y???ƶ?λ??
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6))      #?ú????ϵı?ǩ??б45??
#col_fun = colorRamp2(c(0, 0.3, 0.8), c("blue", "white", "red"))
#rowAnn<-rowAnnotation(Percent = as$percent, No.cells = anno_barplot(as$n),col=list(Percent = col_fun))
##################################
library(ComplexHeatmap)
library(circlize)
m2m2genes<-m2m2genes[m2m2genes$type %in% rownames(m2m2genes_heatmap),]
m2m2genes_heatmap<-t(m2m2genes_heatmap)
m2m2genes_heatmap<-scale(m2m2genes_heatmap)
m2m2genes_heatmap<-t(m2m2genes_heatmap)
pdf("cell_m2m2genes_NMF_macrophage.pdf",height = 5.11,width=3.04)
Heatmap(m2m2genes_heatmap,
        name = "Zscore",
        cluster_rows = TRUE,
        cluster_columns = FALSE, 
        split = m2m2genes$iord ,
        row_names_side =  "left"
        #right_annotation = rowAnn
)
dev.off()
#####################
#HNRNPA2B1  C1QB
#####################
inflammatoryscore<-c("IL1B","IL6","S100A8","S100A9")
anti_inflammatory<-c('CD163','SEPP1','SELENOP', 'APOE','MAF')
inflammatorylist<-list()
inflammatorylist[["inflammatory"]]<-inflammatoryscore
inflammatorylist[["anti_inflammatory"]]<-anti_inflammatory
sce<-AddModuleScore(sce,features=inflammatorylist,names=names(inflammatorylist))
colnames(sce@meta.data)
FeatureScatter(sce,feature1="Cluster1",feature2 ="Cluster2",group.by = "NMFcluster" )
VlnPlot(sce,features=c("Cluster2","Cluster1"),pt.size = 0)
#####################
inflammatoryscores<-cbind(sce$NMFcluster,sce$Cluster1,sce$Cluster2)
inflammatoryscoresmeans<-aggregate(inflammatoryscores[,c(2:3)],list(inflammatoryscores[,1]),mean)#########根据最后确定的TF数量进�?�构�? heatmap
inflammatoryscoresmeans<-t(inflammatoryscoresmeans)
colnames(inflammatoryscoresmeans)<-inflammatoryscoresmeans[1,]
inflammatoryscoresmeans<-inflammatoryscoresmeans[-1,]
rownames(inflammatoryscoresmeans)<-c("inflammatory","anti-inflammatory")
pheatmap(t(inflammatoryscoresmeans),
         #scale = "column", 
         cluster_rows =F, 
         cluster_cols = F,  
         show_rownames = T, angle_col = c('45'),
         #annotation_col = annotations,annotation_legend = FALSE,
         color=colorRampPalette(c("blue","white","red"))(20))
#####################
rt<-t(as.matrix(sce@assays$RNA@data))
markergenes<-sce.markers_top5[sce.markers_top5$gene %in% m6agenes,]$gene

data_plot_mean <- cbind(rownames(rt), sce@meta.data$NMFcluster,rt[,colnames(rt)%in% markergenes], FTH1=rt[,"FTH1"]) %>% data.frame() 
data_plot_mean <- melt(data_plot_mean,
                   id.vars = c("V1","V2"),
                   variable.name ="Signature",
                   value.name = "Signature_score")
data_plot_mean$Signature_score<-as.numeric(data_plot_mean$Signature_score)
#data_plot_mean <- data_frame(sce@meta.data$NMFcluster, sce@meta.data$Cluster2)
#colnames(data_plot_mean) <- c( "NMFcluster","gene")
ggplot(data = data_plot_mean, aes(x = V2, y = Signature_score, fill = Signature)) +
  #geom_bar(stat = "identity", width=0.8,aes(group=NMFcluster),position="fill")+
  geom_boxplot()+
  facet_wrap(data_plot_mean$Signature)+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Expression")+
  #scale_y_continuous(labels = percent)+  ####??��??y???ƶ?λ??
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6))      #?ú????ϵı?ǩ??б45??

ggsave("TOP_GENES_TAM.pdf")

