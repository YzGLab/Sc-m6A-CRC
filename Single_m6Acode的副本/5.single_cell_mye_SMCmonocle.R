#####
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
ggsave("monocle_majorCluster_State_MacSMC.pdf", width = 5, height = 5)

#pdf(file="05.trajectory.Pseudotime.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "Pseudotime")
ggsave("monocle_majorCluster_PseudotimeMAC.pdf", width = 5, height = 5)

#dev.off()
##有时候（大多数时候），拟时序的方向或是根节点弄错了，还需要手动更改
cds=orderCells(cds,root_state = 4) 
plot_cell_trajectory(cds,color_by = "Pseudotime")
#保存时间的细胞轨迹图

#dev.off()
#保存细胞名称的细胞轨迹图
#pdf(file="05.trajectory.cellType.pdf",width=6.5,height=6)
#plot_cell_trajectory(cds,color_by = "cell_type2")
#dev.off()
#保存聚类的细胞轨迹图
#pdf(file="05.trajectory.cluster.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "Cluster")
#dev.off()
ggsave("monocle_majorCluster_Mac.pdf", width = 5, height = 5)
###
plot_cell_trajectory(cds,color_by = "NMFcluster",backbone_color="red",cell_stroke = I(cell_size/2),cell_size = 0.8)+  scale_color_manual(values = mycol)
ggsave("monocle_majorCluster_NMFclusterMacSMC.pdf", width = 5, height = 5)
plot_cell_trajectory(cds, color_by = "Cell_subtype",cell_size = 0.8)+  scale_color_manual(values = mycol)
ggsave("monocle_majorCluster_Cell_subtypeMacSMC.pdf", width = 5, height = 5)
plot_cell_trajectory(cds,markers="WTAP")
plot_cell_trajectory(cds)
##简单的一个函数就可以绘制热图：
##
BEAM_res=BEAM(cds,branch_point = 1,cores = 1)
#会返回每个基因的显著性，显著的基因就是那些随不同branch变化的基因
#这一步很慢
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]
saveRDS(BEAM_res, file = "Mac_BEAM_res.rds")

??plot_pseudotime_heatmap
#
pdf("Mac_pseudotime_Heatmap_m6a.pdf",width = 4.47,height=5.53)
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
                            branch_point = 2,
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
plot_genes_branched_pseudotime(cds[m6agenes[c(1:2)],],
                               branch_point = 1,
                               color_by = "State",
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