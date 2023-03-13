table(sce$Cell_subtype)
colnames(sce@meta.data)
sce_CD8 = sce[, sce$Cell_subtype %in% 
                c( "CD8+ T cells" )] # CD16
sce_CD4 = sce[, sce$Cell_subtype %in% 
                c( "CD4+ T cells" )] # CD16
sce_NK = sce[, sce$Cell_subtype %in% 
               c( "NK cells" )] # CD16
sce_gamma = sce[, sce$Cell_subtype %in% 
                  c( "gamma delta T cells" )] # CD16
sce_Regulatory = sce[, sce$Cell_subtype %in% 
                       c( "Regulatory T cells" )] # CD16
sce_Tfoll = sce[, sce$Cell_subtype %in% 
                  c( "T follicular helper cells" )] # CD16
sce_Thelper = sce[, sce$Cell_subtype %in% 
                    c( "T helper 17 cells" )] # CD16
###############
checkpoints<-c("CD274", "CTLA4","LAG3", "TIM3", "TNFRSF9","TIGIT","CD226","CD7",
               "CD4","CD8A","CD8B","FOXP3","IL2",
               "CXCL8","PDCD1","HAVCR2","GZMB","PRF1","NLG1","IFNG","TNFRSF18","TNFRSF4")
checkpoints_T_heatmap <- AverageExpression(sce_CD8, assays = "RNA",
                                           group.by = "NMFcluster",
                                           features = checkpoints ,verbose = TRUE) %>% .$RNA
checkpoints_T_heatmap2 <- AverageExpression(sce_CD4, assays = "RNA",
                                            group.by = "NMFcluster",
                                            features = checkpoints ,verbose = TRUE) %>% .$RNA
checkpoints_T_heatmap3 <- AverageExpression(sce_NK, assays = "RNA",
                                            group.by = "NMFcluster",
                                            features = checkpoints ,verbose = TRUE) %>% .$RNA
checkpoints_T_heatmap4 <- AverageExpression(sce_gamma, assays = "RNA",
                                            group.by = "NMFcluster",
                                            features = checkpoints ,verbose = TRUE) %>% .$RNA
checkpoints_T_heatmap5 <- AverageExpression(sce_Regulatory, assays = "RNA",
                                            group.by = "NMFcluster",
                                            features = checkpoints ,verbose = TRUE) %>% .$RNA
checkpoints_T_heatmap6 <- AverageExpression(sce_Tfoll, assays = "RNA",
                                            group.by = "NMFcluster",
                                            features = checkpoints ,verbose = TRUE) %>% .$RNA
checkpoints_T_heatmap7 <- AverageExpression(sce_Thelper, assays = "RNA",
                                            group.by = "NMFcluster",
                                            features = checkpoints ,verbose = TRUE) %>% .$RNA
#################################################
checkpoints_T_heatmap <- na.omit(checkpoints_T_heatmap)
colnames(sce@meta.data)
#annotations<-data.frame(sce@meta.data$clusters,sce@meta.data$NMFcluster)
#annotations<-annotations[!duplicated(annotations$sce.meta.data.clusters),]
#rownames(annotations)<-annotations$sce.meta.data.clusters
#annotations<-annotations[order(annotations$sce.meta.data.NMFcluster),]
#checkpoints_T_heatmap<-checkpoints_T_heatmap[,annotations$sce.meta.data.clusters]
pdf("cell_checkpoints_NMF_CD8.pdf",height = 4.44,width=3.02)
pheatmap(checkpoints_T_heatmap,
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = F,  
         show_rownames = T, angle_col = c('45'),
         #annotation_col = annotations,annotation_legend = FALSE,
         color=colorRampPalette(c("blue","white","red"))(20))
dev.off()
pdf("cell_checkpoints_NMF_CD4.pdf",height = 4.44,width=3.02)
pheatmap(checkpoints_T_heatmap2,
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = F,  
         show_rownames = T, angle_col = c('45'),
         #annotation_col = annotations,annotation_legend = FALSE,
         color=colorRampPalette(c("blue","white","red"))(20))
dev.off()
pdf("cell_checkpoints_NMF_NK.pdf",height = 4.44,width=3.02)
pheatmap(checkpoints_T_heatmap3,
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = F,  
         show_rownames = T, angle_col = c('45'),
         #annotation_col = annotations,annotation_legend = FALSE,
         color=colorRampPalette(c("blue","white","red"))(20))
dev.off()
pdf("cell_checkpoints_NMF_Regulatory.pdf",height = 4.44,width=3.02)
pheatmap(checkpoints_T_heatmap5,
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = F,  
         show_rownames = T, angle_col = c('45'),
         #annotation_col = annotations,annotation_legend = FALSE,
         color=colorRampPalette(c("blue","white","red"))(20))
dev.off()
pdf("cell_checkpoints_NMF_Tfoll.pdf",height = 4.44,width=3.02)
pheatmap(checkpoints_T_heatmap6,
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = F,  
         show_rownames = T, angle_col = c('45'),
         #annotation_col = annotations,annotation_legend = FALSE,
         color=colorRampPalette(c("blue","white","red"))(20))
dev.off()
pdf("cell_checkpoints_NMF_Thelper.pdf",height = 4.44,width=3.02)
pheatmap(checkpoints_T_heatmap7,
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = F,  
         show_rownames = T, angle_col = c('45'),
         #annotation_col = annotations,annotation_legend = FALSE,
         color=colorRampPalette(c("blue","white","red"))(20))
dev.off()
###############
library(ComplexHeatmap)
comparsiongenes<-read.delim2("ImmunecomparisonGene2.txt",header = T)
#checkpoints<-c("CD274", "CTLA4","LAG3", "TIM3", "TNFRSF9","TIGIT","CD226","CD7",
             #  "CD4","CD8A","CD8B","FOXP3","IL2","CXCL8","PDCD1","HAVCR2","GZMB",
            # "PRF1","NLG1","IFNG","TNFRSF18","TNFRSF4")
sce$newNMFcluster<-paste0(sce$Cell_subtype,"_",sce$NMFcluster)
cellnames<-c("CD8+ T cells","NK cells", "Regulatory T cells", "CD4+ T cells")
T_heatmap <- AverageExpression(sce[,sce$Cell_subtype %in% cellnames], assays = "RNA",
                                        group.by = "newNMFcluster",
                                        features = comparsiongenes$Genes ,verbose = TRUE) %>% .$RNA
#T_heatmap<-T_heatmap[rownames(T_heatmap) %in% sce.markers$gene,]
T_heatmap<-T_heatmap[rowMeans(T_heatmap)!=0,]
#T_heatmap <- na.omit(T_heatmap)
####
colunmsplit<-table(sce$NMFcluster,sce$Cell_subtype,sce$newNMFcluster) %>% data.frame()
colunmsplit<-colunmsplit[colunmsplit$Freq!=0,]
table(colunmsplit$Var2)
colunmsplit<-colunmsplit[colunmsplit$Var2 %in% cellnames,]
comparsiongenes<- comparsiongenes[ comparsiongenes$Genes %in% rownames(T_heatmap),]
T_heatmap<-t(scale(t(T_heatmap)))
#col_fun = colorRamp2(c(0, 0.3, 0.8), c("blue", "white", "red"))
#rowAnn<-HeatmapAnnotation(Percent = as$percent, No.cells = anno_barplot(as$n),col=list(Percent = col_fun))
##################################
pdf("T_cell_Heatmap_m6aGroup.pdf",width = 6.12,height=7.13)
Heatmap(T_heatmap,
        name = "AveE.",
        cluster_rows = T,
        cluster_columns = F,
        border = T,
        #fontsize=6,
        row_gap = unit(2, "mm"),
        row_split = comparsiongenes$types ,
        row_title_gp = gpar(col = c("red", "blue","orange"), font = 1:3),
        row_names_gp = gpar(col = c("red",  "blue","orange"), fontsize =6),
        row_names_side =  "right",
        column_split = colunmsplit$Var2,
        column_title_gp = gpar(fill = mycol, font = 1:4,fontsize =10,angle="60"),
        #column_names_gp = gpar(col = mycol, fontsize =8))
        #top_annotation = rowAnn
)
dev.off()
####

###
#####dittoBarPlot(sce, "Cell_subtype",group.by="NMFcluster")
#aheatmap(checkpoints_T_heatmap, Rowv=Rowv, Colv=NA, annCol=tmp, annColors=annColors, color=heatCol, revC=TRUE, breaks=0, fontsize=fontsize, labCol = NA, labRow = NA)
############################################################################AUCell
####
#immune chekpoint
Teffectscore<-c('CD8A', 'CXCL10', 'CXCL9', 'GZMA', 'GZMB', 'IFNG', 'PRF1', 'TBX21')
Exhaustionscore<-c('CTLA4','HAVCR2','LAG3', 'PDCD1','TIGIT')
Cytotoxicscore<-c('CST7','GZMA', 'GZMB', 'IFNG', 'NKG7','PRF1')
#T_function<-c("CD3E","CD4","CD8B","FOXP3","GZMB","PRF1","TBX21")
ievgenes<-read.delim2("immune_ecasion_genes.txt")
Immune.sig<-list()
Immune.sig[["Exhaustionscore"]]<-c('CTLA4','HAVCR2','LAG3', 'PDCD1','TIGIT')
Immune.sig[["Cytotoxicscore"]]<-c('CST7','GZMA', 'GZMB', 'IFNG', 'NKG7','PRF1')
Immune.sig[["Teffectscore"]]<-c('CD8A', 'CXCL10', 'CXCL9', 'GZMA', 'GZMB', 'IFNG', 'PRF1', 'TBX21')
Immune.sig[["Tevasionscore"]]<-ievgenes$Genes
####using AUCell to calculate the score for each signature
single_dat<-readRDS("T cellssingle_dat_GSE132465.rds")
rowGenenames<-readRDS("rowGenenames.rds")
rownames(single_dat)<-rowGenenames
single_dat[1:5,1:5]
####
expSet<-as.matrix(single_dat)
library(GSVA)
gsva_es <- gsva(expSet,Immune.sig,method="ssgsea",abs.ranking=F,kcdf="Poisson",parallel.sz=40)###RNA-Seq  data,ssgsea.norm=TRUE
#gsva_es <- gsva(norm.expMat,gmt.list,method="ssgsea",abs.ranking=F,kcdf="Gaussian","Poisson",parallel.sz=10)# Array
gsva_es[1:3,1:4]
load("Tcell_SMA_GSVA.Rdata")
gsva_es<-t(gsva_es)
sce<-AddMetaData(sce, metadata=gsva_es, col.name=colnames(gsva_es))
save(gsva_es,file="Tcell_SMA_GSVA.Rdata")
load("Tcell_SMA_GSVA.Rdata")
##
cells_rankings <- AUCell_buildRankings(as.matrix(single_dat), nCores = 20, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(geneSets = Immune.sig, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
#cells_AUC <- AUCell_calcAUC(geneSets = programG, cells_rankings, aucMaxRank=4)
programAUC <- getAUC(cells_AUC)
programAUC[1:4,1:3]
programAUC[1,1:3]
colnames(sce)
sce<-AddMetaData(sce, metadata=programAUC[1,], col.name=c("Exhaustionscore"))
sce<-AddMetaData(sce, metadata=programAUC[2,], col.name=c("Cytotoxicscore"))
sce<-AddMetaData(sce, metadata=programAUC[3,], col.name=c("Teffectscore"))
sce<-AddMetaData(sce, metadata=programAUC[4,], col.name=c("Tevasionscore"))
###ggplot plot 
colnames(sce@meta.data)
df<- data.frame(sce@meta.data, sce@reductions$umap@cell.embeddings)
head(df)
df$Exhaustionscore<-as.numeric(df$Exhaustionscore)
###
class_avg <- df %>%
  group_by(Cell_subtype) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
ggplot(df, aes(UMAP_1, UMAP_2))  +
  geom_point(aes(colour  = Exhaustionscore)) + 
  viridis::scale_color_viridis(option="D") +
  ggrepel::geom_label_repel(aes(label = Cell_subtype),
                            data = class_avg,
                            size = 6,
                            label.size = 0,
                            segment.color = NA
  )+  
  theme(legend.position = "none") + 
  theme_bw()
###
cellnames<-c("CD8+ T cells","NK cells", "Regulatory T cells", "CD4+ T cells")
selcet_sce<-sce[,sce$Cell_subtype %in% "NK cells"]
FeaturePlot(selcet_sce,feature=c("Exhaustionscore","Cytotoxicscore","Teffectscore","Tevasionscore")) & viridis::scale_color_viridis(option="H")
VlnPlot(selcet_sce,group.by="NMFcluster",features=c("Exhaustionscore","Cytotoxicscore","Teffectscore","Tevasionscore"),pt.size=0)#####
####
###
data_plot_mean <- data.frame(sce@meta.data, sce@reductions$tsne@cell.embeddings)
colnames(data_plot_mean)
data_plot_mean<-data_plot_mean[,c("Index",'NMFcluster',"Exhaustionscore","Cytotoxicscore","Teffectscore","Tevasionscore")]
data_plot_mean <- melt(data_plot_mean,
                       id.vars = c("Index","NMFcluster"),
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
#####

cellnames<-c("CD8+ T cells","NK cells", "Regulatory T cells", "CD4+ T cells")
selcet_sce<-sce[,sce$Cell_subtype %in% "NK cells"]
FeaturePlot(selcet_sce,feature=c("Exhaustionscore","Cytotoxicscore","Teffectscore","Tevasionscore")) & viridis::scale_color_viridis(option="H")
VlnPlot(selcet_sce,group.by="NMFcluster",features=c("Exhaustionscore","Cytotoxicscore","Teffectscore","Tevasionscore"),pt.size=0)#####

scemetadata$newNMFcluster<-paste0(scemetadata$Cell_subtype,"_",scemetadata$NMFcluster)
scemetadata_cells<-scemetadata[scemetadata$Cell_subtype %in% cellnames, ]

colnames(scemetadata_cells)
scemetadata_cells<-scemetadata_cells[,c(46,42:45)]

#aucss<-apply(aucs,2,as.numeric)
#rownames(aucss)<-rownames(aucs)
aucsmenas<-aggregate(scemetadata_cells[,c(2:5)],list(scemetadata_cells[,1]),mean)####
rownames(aucsmenas)<-aucsmenas$Group.1

#aucsmenas<-aucsmenas[!rownames(aucsmenas) %like% "extended",]
pdf("CD4_effector_NMF.pdf",height = 2.3,width=4.32)
pheatmap(t(aucsmenas[1:5,c(2:5)]),
         scale = "row",
         show_colnames =F,
         show_rownames = T,
         cluster_cols = F,
         cluster_rows = F,
         angle_col = c('45'),
         color=colorRampPalette(c("blue","white","red"))(20)
)
dev.off()
pdf("CD8_effector_NMF.pdf",height = 2.3,width=4.32)
pheatmap(t(aucsmenas[6:10,c(2:5)]),
         scale = "row",
         show_colnames =F,
         show_rownames = T,
         cluster_cols = F,
         cluster_rows = F,
         angle_col = c('45'),
         color=colorRampPalette(c("blue","white","red"))(20)
)
dev.off()
pdf("NK_effector_NMF.pdf",height = 2.3,width=4.32)
pheatmap(t(aucsmenas[11:15,c(2:5)]),
         scale = "row",
         show_colnames =F,
         show_rownames = T,
         cluster_cols = F,
         cluster_rows = F,
         angle_col = c('45'),
         color=colorRampPalette(c("blue","white","red"))(20)
)
dev.off()
pdf("Treg_effector_NMF.pdf",height = 2.3,width=4.32)
pheatmap(t(aucsmenas[16:20,c(2:5)]),
         scale = "row",
         show_colnames =F,
         show_rownames = T,
         cluster_cols = F,
         cluster_rows = F,
         angle_col = c('45'),
         color=colorRampPalette(c("blue","white","red"))(20)
)
dev.off()

### tf genes analysis for each type genes 
load("SMC_sce.meta_Tcell.Rdata")
cellnames<-c("CD8+ T cells","NK cells", "Regulatory T cells", "CD4+ T cells")
TFgenes<-read.csv("Tcell_TF_idente.csv",check.names = F,row.names = 1)
#
colnames(TFgenes)
##CD4+ T cells1:5, CD8+ T cells 6-10 gamma delta T cells11-15 NK cells16-20 Regulatory T cells21-25 T follicular helper cells26-30 T helper 17 cells31-35 Unknown36-40
TFgenes<-cbind(group=paste0(scemetadata$Cell_subtype,"_",scemetadata$NMFcluster),TFgenes)
TFgenes<-TFgenes[scemetadata[scemetadata$Cell_subtype %in% cellnames,]$Index,]
aucs<-TFgenes
aucs$group<-as.numeric(factor(aucs$group))
table(aucs$group,TFgenes$group)
#aucs<-data.frame(aucs)
#group_by(aucs, group) %>% summarize_each(funs(mean), colnames(aucs)[1:81])
aucss<-apply(aucs,2,as.numeric)
rownames(aucss)<-rownames(aucs)
aucss[,1]<-aucs[,1]
aucsmenas<-aggregate(aucss[,c(2:227)],list(aucss[,1]),mean)#########根据最后确定的TF数量进行构建 heatmap
aucsmenas<-t(aucsmenas)
colnames(aucsmenas)<-aucsmenas[1,]
aucsmenas<-aucsmenas[-c(1),]
aucsmenas[1:5,1:4]
#aucsmenas<-rbind(unique(TFgenes$group),aucsmenas)
#aucsmenas<-aucsmenas[!rownames(aucsmenas) %like% "extended",]
write.csv(aucsmenas,file="aucsmenas_TFgenes_AUC.csv")
sce.markers_top5 <- sce.markers %>% 
########################CD8+ T cells"
data<-TFgenes[scemetadata[scemetadata$Cell_subtype %in% cellnames[1],]$Index,]
outTab<-c()
for (i in colnames(aucs)[-1]){
  geneName=i
  data[,i]<-as.numeric(data[,i])
  rt=rbind(expression=data[,i],group=data[,"group"])
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
#gene.cox.i[,11] <- p.adjust(gene.cox.i[,9], method = "fdr")
outTab$fdr<- p.adjust(outTab$pValue, method = "fdr")
outTab<-outTab[outTab$fdr<0.0001,]
TF_heatmap<-aucsmenas[outTab$gene,6:10]
pdf("TF_scenic_NMF_CD8SMC.pdf",height = 7.83,width=2.89)
pheatmap(TF_heatmap,
         scale = "row",
         show_colnames =T,
         show_rownames = T,
         cluster_cols = F,
         fontsize = 6,
         angle_col = c('0'),
         color=colorRampPalette(c("blue","white","red"))(20)
)
dev.off()
###CD8 tfgenes
cd8tfgenes<-outTab$gene

###################"NK cells"
data<-TFgenes[scemetadata[scemetadata$Cell_subtype %in% cellnames[2],]$Index,]
outTab<-c()
for (i in colnames(aucs)[-1]){
  geneName=i
  data[,i]<-as.numeric(data[,i])
  rt=rbind(expression=data[,i],group=data[,"group"])
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
#gene.cox.i[,11] <- p.adjust(gene.cox.i[,9], method = "fdr")
outTab$fdr<- p.adjust(outTab$pValue, method = "fdr")
outTab<-outTab[outTab$fdr<0.0001,]
TF_heatmap<-aucsmenas[outTab$gene,11:15]
pdf("TF_scenic_NMF_NKSMC.pdf",height = 3.84,width=3.52)
pheatmap(TF_heatmap,
         scale = "row",
         show_colnames =T,
         show_rownames = T,
         cluster_cols = F,
         fontsize = 6,
         angle_col = c('0'),
         color=colorRampPalette(c("blue","white","red"))(20)
)
dev.off()
nktfgenes<-outTab$gene
#####################################Regulatory T cells
data<-TFgenes[scemetadata[scemetadata$Cell_subtype %in% cellnames[3],]$Index,]
outTab<-c()
for (i in colnames(aucs)[-1]){
  geneName=i
  data[,i]<-as.numeric(data[,i])
  rt=rbind(expression=data[,i],group=data[,"group"])
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
#gene.cox.i[,11] <- p.adjust(gene.cox.i[,9], method = "fdr")
outTab$fdr<- p.adjust(outTab$pValue, method = "fdr")
outTab<-outTab[outTab$fdr<0.0001,]
TF_heatmap<-aucsmenas[outTab$gene,16:20]
pdf("TF_scenic_NMF_Treg.pdf",height = 3.924,width=2.89)
pheatmap(TF_heatmap,
         scale = "row",
         show_colnames =T,
         show_rownames = T,
         cluster_cols = F,
         fontsize = 6,
         angle_col = c('0'),
         color=colorRampPalette(c("blue","white","red"))(20)
)
dev.off()

regttfgenes<-outTab$gene
##########################################"CD4+ T cells"
data<-TFgenes[scemetadata[scemetadata$Cell_subtype %in% cellnames[4],]$Index,]
outTab<-c()
for (i in colnames(aucs)[-1]){
  geneName=i
  data[,i]<-as.numeric(data[,i])
  rt=rbind(expression=data[,i],group=data[,"group"])
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
#gene.cox.i[,11] <- p.adjust(gene.cox.i[,9], method = "fdr")
outTab$fdr<- p.adjust(outTab$pValue, method = "fdr")
outTab<-outTab[outTab$fdr<0.0001,]
TF_heatmap<-aucsmenas[outTab$gene,1:5]

pdf("TF_scenic_NMF_CD4.pdf",height = 7.73,width=2.89)
pheatmap(TF_heatmap,
         scale = "row",
         show_colnames =T,
         show_rownames = T,
         cluster_cols = F,
         fontsize = 6,
         angle_col = c('0'),
         color=colorRampPalette(c("blue","white","red"))(20)
)
dev.off()
cd4tfgenes<-outTab$gene###
###
##
library(ComplexHeatmap)
library(circlize)
library(janitor)
###
tfgeneslist<-list()
tfgeneslist[["cd8"]]<-cd8tfgenes
tfgeneslist[["cd4"]]<-cd4tfgenes
tfgeneslist[["nk"]]<-nktfgenes
tfgeneslist[["treg"]]<-regttfgenes
###
names(tfgeneslist)
tfgeneslist2<-tfgeneslist[c(1,3,4)]
a<-Reduce(intersect,tfgeneslist2)
###
TF_heatmap<-aucsmenas[a,]
TF_heatmap<-t(scale(t(TF_heatmap)))

#as<-tabyl(meta.datas, NMFcluster)
cellnames
#combined<-scemetadata[scemetadata$Cell_subtype %in% cellnames,]
splits<-cbind(rep("cd4",5),rep("cd8",5),rep("nk",5),rep("treg",5)) %>% melt()
splits<-data.frame(splits)
splits<-cbind(TFgenes$group,splits)
##################################
pdf("TF_heatmap_Tcell_m6aGroup.pdf",width = 6,height=6)
ht1<-Heatmap(TF_heatmap,
        name ="TF activity",
        cluster_rows = T,
        cluster_columns = F,
        border = T,
        #fontsize=6,
        row_gap = unit(2, "mm"),
        heatmap_legend_param = list(
          legend_direction = "vertical",
          legend_width = unit(2.5, "cm")
        ),
        #row_split = comparsiongenes$types ,
        #row_title_gp = gpar(col = c("red", "blue","orange"), font = 1:3),
        #row_names_gp = gpar(col = c("red",  "blue","orange"), fontsize =6),
        row_names_side =  "right",
        column_split = splits$value,
        #column_title_gp = gpar(fill = mycol, font = 1:4,fontsize =10,angle="60"),
        #column_names_gp = gpar(col = mycol, fontsize =8))
        #top_annotation = rowAnn
)
draw(ht1, heatmap_legend_side = "top")
dev.off()

