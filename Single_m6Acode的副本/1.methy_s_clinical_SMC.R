##clincal data difference of m6Agenes
library(dittoSeq)
colnames(phe)
ann_cols<-c("#223D6C" ,"#D20A13", "#FFD121" ,"#088247", "#58CDD9","#5D90BA","#223D6C","#D20A13")
pdf("Heatmap_m6A_GSE132465.pdf",height=2.2, width=6.8)
#,annot.colors=ann_cols,
dittoHeatmap(subset(single_dat_sce, downsample = 10), m6agenes, use_raster = TRUE, scaled.to.max = TRUE,annot.colors=ann_cols,scale = "row",
             fontsize_number = 0.2,
             fontsize = 5,
             rownames="left",
             annot.by = c("Cell_type"),color.by=c("red","blue"))
dev.off()

pdf("Heatmap_m6A_Class_GSE132465.pdf",height=2.2, width=6.8)
dittoHeatmap(subset(single_dat_sce, downsample = 10),
             m6agenes, 
             use_raster = TRUE, 
             scaled.to.max = TRUE,
             annot.colors=ann_cols,scale = "row",
             fontsize_number = 0.2,
             fontsize = 5,
             rownames="left",
             annot.by ="Class",color.by=c("red","blue"))
dev.off()
##
pdf("Heatmap_m6A_Class_GSE132465.pdf",height=2.2, width=6.8)
dittoHeatmap(subset(single_dat_sce, downsample = 10),
             m6agenes, 
             use_raster = TRUE, 
             scaled.to.max = TRUE,
             annot.colors=ann_cols,scale = "row",
             fontsize_number = 0.2,
             fontsize = 5,
             rownames="left",
             annot.by ="Gender",color.by=c("red","blue"))
dev.off()
dittoHeatmap(subset(single_dat_sce,downsample=50) ,m6agenes, scaled.to.max = TRUE,
             complex = TRUE,
             use_raster = TRUE)
?dittoHeatmap
###
AveExpression <- AverageExpression(single_dat_sce, assays = "RNA",group.by = "Sample",features = features,verbose = TRUE) %>% .$RNA
AveExpression<-data.frame(AveExpression)
AveExpression$Gene <- rownames(AveExpression)
Ave_df <- reshape2::melt(AveExpression,id.vars= "Gene")
colnames(Ave_df) <- c("Gene", "Cluster", "Expression")


features<-m6agenes
AveExpression<-list()
colnames(single_dat_sce@meta.data)
newphe$Agegroup<-ifelse(newphe$Age>60,"Older","Young")
single_dat_sce@meta.data$Agegroup<-ifelse(single_dat_sce@meta.data$Age>60,"Older","Young")
####extract the Average of expression in each group.
AveExpression[[1]] <- AverageExpression(single_dat_sce, assays = "RNA",group.by = "Cell_type",features = features,verbose = TRUE) %>% .$RNA
AveExpression[[2]] <- AverageExpression(single_dat_sce, assays = "RNA",group.by = "Class",features = features,verbose = TRUE) %>% .$RNA
AveExpression[[3]] <- AverageExpression(single_dat_sce, assays = "RNA",group.by = "MSI",features = features,verbose = TRUE) %>% .$RNA
AveExpression[[4]] <- AverageExpression(single_dat_sce, assays = "RNA",group.by = "Gender",features = features,verbose = TRUE) %>% .$RNA
AveExpression[[5]] <- AverageExpression(single_dat_sce, assays = "RNA",group.by = "Agegroup",features = features,verbose = TRUE) %>% .$RNA
AveExpression[[6]] <- AverageExpression(single_dat_sce, assays = "RNA",group.by = "Stage",features = features,verbose = TRUE) %>% .$RNA
AveExpression[[7]] <- AverageExpression(single_dat_sce, assays = "RNA",group.by = "nearestCMS..RF.",features = features,verbose = TRUE) %>% .$RNA
AveExpressions<-c()
for( i in c(1:7)){
  AveExpressions<-cbind(AveExpressions,AveExpression[[i]])
}
####
table(single_dat_sce@meta.data$Stage)
splitgroups<-c(rep("Cell_type",6),rep("Location",2),rep("MSI",2),rep("Gender",2),rep("Age",2),rep("Stage",6),rep("sc-CMS",4))
splitgroups<-data.frame(splitgroups)
#####
#install.packages("janitor")
library(janitor)
a<-list()
a[[1]]<-tabyl(newphe, Cell_type)
a[[2]]<-tabyl(newphe, Class)
a[[3]]<-tabyl(newphe, MSI)
a[[4]]<-tabyl(newphe, Gender)
a[[5]]<-tabyl(newphe, Agegroup)
a[[6]]<-tabyl(newphe, Stage)
a[[7]]<-tabyl(newphe, nearestCMS..RF.)

as<-rbindlist(a,use.names=FALSE, fill=FALSE, idcol=NULL)
as<-cbind(as,splitgroups)
colnames(as)[1]<-"Subgroup"
#####

library(ComplexHeatmap)
library(circlize)
AveExpressions<-t(AveExpressions)
AveExpressions<-scale(AveExpressions)
rownames(AveExpressions)

col_fun = colorRamp2(c(0, 0.3, 0.8), c("blue", "white", "red"))
rowAnn<-rowAnnotation(Percent = as$percent, No.cells = anno_barplot(as$n),col=list(Percent = col_fun))

pdf("cell_m6aGenes_clinical.pdf",width = 7.81,height=6.69)
Heatmap(AveExpressions,name = "Raltive AveExpression",
        cluster_rows = TRUE,
        cluster_columns = TRUE, 
        row_split = as$splitgroups ,
        row_names_side =  "left",
        right_annotation = rowAnn)
dev.off()

#################



library(readxl)
newphe<-data.frame(read_excel("clinicalphe.xlsx",1, col_names= TRUE, col_types= NULL, na=" "))##chip data

colnames(phe)
colnames(newphe)
newphe<-merge(phe,newphe,by="Sample",all=T)
rownames(newphe)<-newphe$Index
single_dat_sce<-AddMetaData(single_dat_sce,metadata=newphe, col.name=colnames(newphe))
dim(newphe)
colnames(single_dat_sce@meta.data)
d<-single_dat_sce@meta.data
save(d,file="meta.data.rda")
saveRDS(newphe,file="newGSE132465_newphe.rds")

dittoDotPlot(single_dat_sce, vars = m6agenes, group.by = "Sample",scale = TRUE)

DotPlot(single_dat_sce, features = m6agenes, group.by = "Class")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12,  hjust = 0.8,angle=45, colour = "black"))+
  theme(legend.position = "right",
        panel.border = element_blank(),## 去掉最外层的正方形边框 
        axis.ticks.x = element_line(color =  NA))

DotPlot(single_dat_sce, features = m6agenes, group.by = "Gender")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12,  hjust = 0.8,angle=45, colour = "black"))+
  theme(legend.position = "right",
        panel.border = element_blank(),## 去掉最外层的正方形边框 
        axis.ticks.x = element_line(color =  NA))

DotPlot(single_dat_sce, features = m6agenes, group.by = "Stage")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12,  hjust = 0.8,angle=45, colour = "black"))+
  theme(legend.position = "right",
        panel.border = element_blank(),## 去掉最外层的正方形边框 
        axis.ticks.x = element_line(color =  NA))
DotPlot(single_dat_sce, features = m6agenes, group.by = "nearestCMS..RF.")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12,  hjust = 0.8,angle=45, colour = "black"))+
  theme(legend.position = "right",
        panel.border = element_blank(),## 去掉最外层的正方形边框 
        axis.ticks.x = element_line(color =  NA))
DotPlot(single_dat_sce, features = m6agenes, group.by = "MSI")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12,  hjust = 0.8,angle=45, colour = "black"))+
  theme(legend.position = "right",
        panel.border = element_blank(),## 去掉最外层的正方形边框 
        axis.ticks.x = element_line(color =  NA))


DotPlot(single_dat_sce, features = m6agenes, group.by = "Sample")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12,  hjust = 0.8,angle=45, colour = "black"))+
  theme(legend.position = "right",
        panel.border = element_blank(),## 去掉最外层的正方形边框 
        axis.ticks.x = element_line(color =  NA))
###
colnames(single_dat_sce@meta.data)
dittoScatterPlot(
  object = single_dat_sce,
  x.var = "WritersScore", y.var = 'ImmuneScore',color.panel = mycol,
  color.var = "Class")
dittoScatterPlot(
  object = single_dat_sce,
  x.var = "ErasersScore", y.var = 'ImmuneScore',color.panel = mycol,
  color.var = "Class")
dittoScatterPlot(
  object = single_dat_sce,
  x.var = "ReadersScore", y.var = 'ImmuneScore',color.panel = mycol,
  color.var = "Class")
######################### combined sample with m6A expression 
library(scales)
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
  scale_y_continuous(position = "left",labels = percent)+  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(legend.position="top")+
  guides(fill=guide_legend(ncol=6))+
  theme(axis.text.x = element_text(size=12,  hjust = 0.5, vjust=0.5,angle=0, colour = "black"))+ coord_flip()##翻转过来
pB2
ggsave("patients with cell barplot.pdf")## 5.96 x 7.26

library(ggplot2)
AveExpression_sample_or <- AverageExpression(single_dat_sce, assays = "RNA",group.by = "Sample",features = features,verbose = TRUE) %>% .$RNA
AveExpression_sample_or<-t(scale(t(AveExpression_sample_or)))
pheatmap(t(AveExpression_sample_or),cluster_rows = F)
ord <- hclust(dist(AveExpression_sample_or, method = "euclidean"), method = "ward.D" )$order
ord
table(rownames(AveExpression_sample_or))
AveExpression_sample<-t(AveExpression_sample_or)
AveExpression_sample <-scale(AveExpression_sample)

AveExpression_sample<-AveExpression_sample  %>% reshape2::melt()
colnames(AveExpression_sample) <- c("Sample","m6Agenes","expression")
AveExpression_sample$m6Agenes<-factor(AveExpression_sample$m6Agenes,levels=rownames(AveExpression_sample_or)[ord])
table(AveExpression_sample$m6Agenes)
AveExpression_sample<-data.frame(AveExpression_sample)
AveExpression_sample$expression<-as.numeric(AveExpression_sample$expression)
pB1 <- ggplot(data =AveExpression_sample, aes(x =Sample , y =m6Agenes , fill = expression)) +
  #scale_fill_distiller(palette = "RdPu") +
  scale_fill_gradient2(low="darkblue", mid="white",high="red") +
  geom_tile()+
  #, "red"
  #theme(legend.position = "")+
  #guides(fill=guide_legend(ncol=6))+
  #geom_box(stat = "identity", width=0.8,position="fill")+
  #scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="m6AGenes")+
  #scale_y_continuous(position = "right",labels = percent)+  ####用来将y轴移动位置
  theme(axis.text.y =element_blank())+
  #theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=8,  hjust = 1,angle=45, colour = "black"))+ coord_flip()##翻转过来
pB1
pB2+pB1
ggsave("patients with cell barplotpB2pB1.pdf")## 5.96 x 7.26
library(ggpubr)
ggarrange(pB2, pB1, 
          labels = c("A", "B"),
          ncol = 3, nrow = 1)


#p1<-Heatmap(AveExpression_sample ,
     #   name = "Raltive AveExpression",
      #  cluster_rows = FALSE,
       # cluster_columns = TRUE, 
        #row_split = as$splitgroups ,
       # row_names_side =  "left")
        #right_annotation = rowAnn)
