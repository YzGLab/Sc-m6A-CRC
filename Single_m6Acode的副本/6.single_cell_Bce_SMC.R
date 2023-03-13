###for B cell
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
#load("Pan_bcell_cells_CRC.Rdata")
dir("O_data")
#single_dat<-readRDS("O_data/bcellssingle_dat_GSE132465.rds")
single_dat<-readRDS("B cellssingle_dat_GSE132465.rds")
rowGenenames<-readRDS("rowGenenames.rds")
rownames(single_dat)<-rowGenenames
single_dat[1:5,1:5]
GSE132465_phe <- readRDS('GSE132465_phe.rds')#####
table(GSE132465_phe$Cell_type,GSE132465_phe$Cell_subtype)
#gmt=read.gmt(misfile) #   input from the MsiDB data.    
m6agenes<-c("CBLL1","KIAA1429","METTL14","METTL3", "RBM15","RBM15B", "WTAP",
            "ZC3H13","ELAVL1","FMR1", "HNRNPA2B1", "HNRNPC", "IGF2BP1", 
            "IGF2BP2", "IGF2BP3", "LRPPRC", "YTHDC1", "YTHDC2", "YTHDF1", 
            "YTHDF2","YTHDF3","ALKBH5","FTO")
mycol <- c("#223D6C","#D20A13","#FFD121","#088247",
           "#58CDD9","#5D90BA","#431A3D","#11AA4D","#91612D","#6E568C","#7A142C",
           "#E0367A","#D8D155","#64495D","#7CC767")
##
sce <- CreateSeuratObject(counts =single_dat, project = "GSE132465_bcell_colon")
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
sce <- FindNeighbors(sce, dims = 1:15)
##########
sce<- FindClusters(sce, resolution = 0.1)
##################RunUMAP
sce<- RunUMAP(sce, dims = 1:20)
DimPlot(sce, reduction = "umap",label=T,cols=mycol)
#RunTSNE
sce<- RunTSNE(object = sce, dims = 1:15, do.fast = TRUE)
DimPlot(sce,reduction = "tsne",label=T,cols=mycol)
##
phe<-GSE132465_phe
sce=AddMetaData(object=sce, metadata=phe, col.name=colnames(phe))##
####################
table(sce$Cell_subtype)
sce<-sce[,sce$Cell_subtype %in% c('CD19+CD20+ B','IgA+ Plasma' ,'IgG+ Plasma')]
#sce<-sce[,sce$Class=="Tumor"]
######
colnames(sce@meta.data)
##########
DimPlot(sce,reduction = "tsne",group.by="Cell_subtype",label=T,cols=mycol)
DimPlot(sce,reduction = "tsne",group.by="Class",label=T,cols=mycol)
ggsave("bcell_Class.SMC.pdf")
DimPlot(sce,reduction = "tsne",group.by="Sample")
ggsave("bcell_Sample.SMC.pdf")
DimPlot(sce,reduction = "tsne",group.by="Cell_subtype",label=T,cols=mycol)
ggsave("bcellType.SMC.pdf")
DimPlot(sce,reduction = "tsne",cols=mycol)
DimPlot(sce,reduction = "umap",cols=mycol)
##add the phe data
newphe<-data.frame(read_excel("clinicalphe.xlsx",1, col_names= TRUE, col_types= NULL, na=" "))##chip data
newphe<-merge(phe,newphe,by="Sample",all=T)
rownames(newphe)<-newphe$Index
sce<-AddMetaData(sce,metadata=newphe, col.name=colnames(newphe))######
colnames(sce@meta.data)
####
##############
#nmf.input<-single_dat_myoBcell[m6agenes,]
nmf.input<-single_dat[m6agenes,colnames(sce)]
nmf.input<-nmf.input[rowSums(nmf.input!= 0)>=1,colSums(nmf.input != 0) >= 1]
rank <- 3
seed <- 2019620
#rownames(nmf.input) <- gsub("Signature","Sig",rownames(nmf.input)) #
mut.nmf <- nmf(nmf.input, 
               rank = rank, 
               seed = seed,  
               #method="nsNMF",
               method = "lee") 
group_Bcell<- predict(mut.nmf) # 
save(group_Bcell,file="group_Bcell_SMC_NMF.Rdata")
##
load("group_Bcell_SMC_NMF.Rdata")
sce<-AddMetaData(sce,metadata=group_Bcell, col.name="NMFcluster")
sce@meta.data$NMFcluster<-ifelse(is.na(sce@meta.data$NMFcluster),"4",sce@meta.data$NMFcluster)
DimPlot(sce,reduction = "tsne",group.by="NMFcluster",cols=mycol)

###
data_plotC <- table(sce@meta.data$Cell_subtype, sce@meta.data$NMFcluster) %>% melt()
colnames(data_plotC) <- c("Cell_subtype", "NMFcluster","Number")
data_plotC$NMFcluster<-as.factor(data_plotC$NMFcluster)
pC1 <- ggplot(data = data_plotC, aes(x = Cell_subtype , y = Number, fill =NMFcluster)) +
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
  #scale_x_discrete(position = "top") +
  #scale_y_reverse()+theme(legend.position="none")
#scale_y_discrete(position = "left") 
#
#

pC1
library(dplyr)
library(plyr)
library(scales)
pC2 <- ggplot(data = data_plotC, aes(x =Cell_subtype , y = Number, fill =NMFcluster )) +
  geom_bar(stat = "identity", width=0.8,aes(group=NMFcluster),position="fill")+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  scale_y_continuous(labels = percent)+  ####用来将y轴移动位置
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6)) + coord_flip()
# + scale_x_discrete(position = "top") 
#让横轴上的标签倾斜45度
pC2
library(ggpubr)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(scales)
pC <- pC1 + pC2 + plot_layout(ncol = 1, widths = c(1,1),guides = 'collect')
pC
######################
ggsave("B_cell_subtype_NMF.pdf")
##
library(ComplexHeatmap)
library(circlize)
library(janitor)
meta.datas<-sce@meta.data
as<-tabyl(meta.datas, NMFcluster)
colnames(as)[1]<-"Subgroup"
AveExpression <- AverageExpression(sce, assays = "RNA",group.by = "NMFcluster", features =m6agenes,verbose = TRUE) %>% .$RNA
AveExpressions<-t(AveExpression)
AveExpressions<-scale(AveExpressions)
rownames(AveExpressions)
col_fun = colorRamp2(c(0, 0.3, 0.8), c("blue", "white", "red"))
rowAnn<-HeatmapAnnotation(Percent = as$percent, No.cells = anno_barplot(as$n),col=list(Percent = col_fun))
##################################
AveExpressions<-AveExpressions[,-13]
pdf("T_cell_m6aGenes_Heatmap_m6aGroup.pdf",width = 5.42,height=4.83)
Heatmap(t(AveExpressions),name = "AveE.",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        #row_names_size=6,
        #row_split = as$splitgroups ,
        row_names_side ="left",
        top_annotation = rowAnn
)
dev.off()
###################################
sce$NMFcluster<-as.factor(sce$NMFcluster)
Idents(sce)<-sce$NMFcluster
sce.markers <- FindAllMarkers(sce,
                              only.pos = TRUE, 
                              min.pct = 0.15, logfc.threshold = 0.15)
sce.markers<-sce.markers[sce.markers$p_val_adj<0.05,]
sce.markers<-split(sce.markers$gene,sce.markers$cluster)
names(sce.markers)<-c("B_C1","B_C2","B_C3","B_C4")
cell_subtype_signatures[["B Cells"]]<-sce.markers
sce.markers  %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
sce.markers_top5 <- sce.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 3, wt = avg_log2FC)
DoHeatmap(sce,features = sce.markers_top5$gene) & NoLegend()

sce.markers<-sce.markers[sce.markers$p_val_adj<0.05,]
write.csv(sce.markers,file="Bcell_sce.markers.csv")

features<-c('HNRNPA2B1',
            'WTAP',
            'TPT1',
            'PFDN5',
            'HNRNPC',
            'SERTAD1',
            'TUBB4B',
            'TCF25',
            'GLA',
            'C12orf57')

DoHeatmap(sce,features = features) & NoLegend()
DotPlot(sce,features = features) & NoLegend() + RotatedAxis()
ggsave("bcell_dotplot_sig.marker.pdf")
###
clustersscore<-split(sce.markers$gene,sce.markers$cluster)
sce<- AddModuleScore(object = sce, features = clustersscore, name = "C") #AddModuleScore计算得分
colnames(sce@meta.data)
FeaturePlot(sce,reduction = "tsne",features = c("C1","C2","C3","C4")) & viridis::scale_color_viridis(option="H")
###############
library(dittoSeq)
#levels(Idents(sce))
#sce = sce[, Idents(sce) %in% 
#   c( "FCGR3A+ Mono", "CD14+ Mono"  )] # CD16 
dittoDimPlot(sce,reduction = "tsne","NMFcluster",split.by="Class")
###\
###########准备细胞轨
library(monocle)
library(scRNAseq)
library(celldex)
library(SingleR)
library(monocle)
#准备细胞轨迹分析需要的文件
####
## do table for 
dim(sce)
selcet_sce<-sce
Idents(selcet_sce)<-selcet_sce$NMFcluster
########################
monocle.matrix=as.matrix(selcet_sce@assays$RNA@data)
monocle.sample=selcet_sce@meta.data
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), 
                           row.names = row.names(monocle.matrix))
#monocle.clusterAnn=clusterAnn
#selcet_sce <- FindNeighbors(selcet_sce, dims = 1:20)
#selcet_sce<- FindClusters(selcet_sce, resolution = 0.1)
#selcet_sce.markers <- FindAllMarkers(selcet_sce,only.pos = TRUE,min.pct = 0.1, logfc.threshold = 0.1)
#monocle.markers=selcet_sce.markers
#将Seurat结果转换为monocle需要的细胞矩阵，细胞注释表格和基因注释表格
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])

#添加细胞聚类数据
#clusterAnn=as.character(monocle.clusterAnn[,2])
#names(clusterAnn)=paste0("cluster",monocle.clusterAnn[,1])
#pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)

#细胞轨迹分析流程
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
##
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
#plot_ordering_genes(cds)
#cds <- setOrderingFilter(cds,as.vector(monocle.markers$gene))
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree')
cds <- orderCells(cds)
#保存树枝的细胞轨迹图
#pdf(file="05.trajectory.State.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "State",backbone_color="red",cell_size = 0.8)+  scale_color_manual(values = mycol)
ggsave("monocle_majorCluster_State_TcellSMC_Bcell.pdf", width = 5, height = 5)

#pdf(file="05.trajectory.Pseudotime.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "Pseudotime")
ggsave("monocle_majorCluster_Pseudotime_Bcell.pdf", width = 5, height = 5)

#dev.off()
##有时候（大多数时候），拟时序的方向或是根节点弄错了，还需要手动更改
cds=orderCells(cds,root_state = 4) 
plot_cell_trajectory(cds,color_by = "Pseudotime")
#保存时间的细胞轨迹图
#dev.off()

#保存聚类的细胞轨迹图
#pdf(file="05.trajectory.cluster.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "Cluster")
#dev.off()
ggsave("monocle_majorCluster_Bcell.pdf", width = 5, height = 5)
###
plot_cell_trajectory(cds,color_by = "NMFcluster",backbone_color="red",cell_stroke = I(cell_size/2),cell_size = 0.8)+  scale_color_manual(values = mycol)
ggsave("monocle_majorCluster_NMFclusterBcellSMC.pdf", width = 5, height = 5)
plot_cell_trajectory(cds, color_by = "Cell_subtype",cell_size = 0.8)+  scale_color_manual(values = mycol)
ggsave("monocle_majorCluster_Cell_subtypeBcellSMC.pdf", width = 5, height = 5)
plot_cell_trajectory(cds,markers="WTAP")
plot_cell_trajectory(cds)
##简单的一个函数就可以绘制热图：
##
BEAM_res=BEAM(cds,branch_point = 1,cores = 1)
#会返回每个基因的显著性，显著的基因就是那些随不同branch变化的基因
#这一步很慢
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]
saveRDS(BEAM_res, file = "Tcell_BEAM_res.rds")
cds$NMFcluster
??plot_pseudotime_heatmap
#
pdf("Bcell_pseudotime_Heatmap_m6a.pdf",width = 4.47,height=5.53)
plot_pseudotime_heatmap(cds[m6agenes,],
                        #branch_point = 1,
                        #add_annotation_row=cds$NMFcluster,
                        #cluster_rows=cds$NMFcluster,
                        #branch_labels = cds$NMFcluster,
                        num_clusters =5,
                        cores = 1,
                        #hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                        use_gene_short_name = T,
                        return_heatmap = T,
                        show_rownames = T)
dev.off()
####
###
plot_genes_branched_heatmap(cds[m6agenes,],
                            #branch_point = 1,
                            num_clusters = 4, #这些基因被分成几个group
                            cores = 3,
                            branch_labels = c("NMF1", "NMF2"),
                            #hmcols = NULL, #默认值
                            hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                            #branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                            use_gene_short_name = T,
                            show_rownames = T,
                            return_heatmap = F)
#是否返回一些重要信息
plot_genes_branched_pseudotime(cds[sce.markers_top5$gene[3:5],],
                               branch_point = 1,
                               color_by = "NMFcluster",
                               cell_size=2,
                               ncol = 2)

##提取热图基因
p=plot_pseudotime_heatmap(cds[ordering_genes,],
                          num_clusters = 3,
                          cores = 1,return_heatmap=T,
                          show_rownames = T)
p$tree_row
clusters <- cutree(p$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)####
###
###write.csv(t(as.matrix(sce@assays$RNA@counts)),file = "Bcell.csv")
