
#############
#1::::###
rm(list = ls()) 
library(Seurat) 
# data from SMC 
#load(file = 'sce-monocyte.Rdata')
sce <-single_dat_sce
table(Idents(sce))
phe=sce@meta.data   
mat=sce@assays$RNA@counts

mat[1:4,1:4]
exprMat =as.matrix(mat) 
write.table(mat,file="SampleExpression.txt")

write.table(mat,file="SampleExpression.txt",sep = "\t")

dim(exprMat)
exprMat[1:4,1:4] 
head(phe)

cellInfo <-  phe[,c('CAF_m6aGroup','nCount_RNA' ,'nFeature_RNA' )]
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')
head(cellInfo)
table(cellInfo$CellType)
cellInfo$CellType=Idents(sce)
table(cellInfo$CellType)

### Initialize settings
# https://github.com/aertslab/SCENIC
# https://pyscenic.readthedocs.io/en/latest/

###
#  https://resources.aertslab.org/cistarget/
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
# mc9nr: Motif collection version 9: 24k motifs
dbFiles
for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
  #  (1041.7 MB)
  # 
}

###
library(SCENIC)
# https://resources.aertslab.org/cistarget/

library(SCENIC)
db='cisTarget_databases/'
list.files(db)
# 保证cisTarget_databases 文件夹下面有下载好 的文件
scenicOptions <- initializeScenic(org="hgnc", 
                                  dbDir=db , nCores=20) 
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
saveRDS(cellInfo, file="int/cellInfo.Rds")

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
length(genesKept)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
dim(exprMat_filtered)
############
runCorrelation(exprMat_filtered, scenicOptions)
#############
exprMat_filtered_log <- log2(exprMat_filtered+1) 
# 最耗费时间的就是这个步骤
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
# 这个步骤也很耗时
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget")) # Toy run settings
library(doParallel)
# 因为莫名其妙的错误，需要把 多线程重新设置成为 1 个线程
scenicOptions <- initializeScenic(org="hgnc", 
                                  dbDir=db , nCores=1) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
# 
# Binary regulon activity: 29 TF regulons x 93 cells.
# (34 regulons including 'extended' versions)
# 29 regulons are active in more than 1% (0.93) cells.

tsneAUC(scenicOptions, aucType="AUC") # choose settings

export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# 运行 



rm(list = ls())  
library(SCENIC)

library(SCopeLoomR)
scenicLoomPath='output/scenic.loom'
loom <- open_loom(scenicLoomPath)
# Read information from loom file:
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
exprMat_log[1:4,1:4] 

scenicOptions=readRDS(file="int/scenicOptions.Rds")
scenicOptions
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
library(tidyverse)
savedSelections <- shiny::runApp(aucellApp)

###
rm(list = ls())  
library(SCENIC)
packageVersion("SCENIC")

scenicOptions=readRDS(file="int/scenicOptions.Rds")
scenicOptions

scenicOptions@inputDatasetInfo

library(SCopeLoomR)
scenicLoomPath='output/scenic.loom'
loom <- open_loom(scenicLoomPath)
# Read information from loom file:
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
exprMat_log[1:4,1:4] 
dim(exprMat_log)

regulons_incidMat <- get_regulons(loom, column.attr.name="MotifRegulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom)
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom) 
#cellClusters <- get_clusterings(loom)
close_loom(loom)

rownames(regulonAUC)
names(regulons)

head(names(regulons))
regulons[[1]]
#load(file = 'sce-monocyte.Rdata')
sce <-single_dat_sce###构建好的identy sce seurat 对象

DimPlot(sce, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()

# 检查任意一个转录因子 
sg_list=lapply(regulons, function(x) x[x%in% rownames(sce)])
sg_list
th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))
DotPlot(object = sce, 
        features=sg_list[[28]] ,  # [28] "HES1"
        assay = "RNA") + th

length(sg_list[[28]])
# 30 , "HES1 (16g)"  
getAUC(regulonAUC[30,])
hist(getAUC(regulonAUC[30,]))
boxplot( as.numeric(getAUC(regulonAUC[30,]) ) ~ 
           as.character( Idents(sce)))


# [7] "JUN_extended"
DotPlot(object = sce, 
        features=sg_list[[7]],  
        assay = "RNA") + th
getAUC(regulonAUC[1,])
hist(getAUC(regulonAUC[1,]))
boxplot( as.numeric(getAUC(regulonAUC[1,]) ) ~ 
           as.character( Idents(sce)))
library(ggpubr)
df=data.frame(value=as.numeric(getAUC(regulonAUC[1,]))  ,
              group= as.character( Idents(sce)))
ggboxplot(df, "group", "value",
          color = "group", 
          add = "jitter", shape = "group")+ stat_compare_means(method = "t.test")


library(pheatmap)
pheatmap(  getAUC(regulonAUC[,] ),show_colnames = F)
n=t(scale(t( getAUC(regulonAUC[,] )))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = T)
ac=data.frame(group= as.character( Idents(sce)))
rownames(ac)=colnames(n)

n2<-n[!rownames(n) %like% 'extended',]
pheatmap(n2,show_colnames =F,show_rownames = T,cluster_cols = T,
         annotation_col=ac)
pheatmap(n,show_colnames =F,show_rownames = T,
         annotation_col=ac,
         filename = 'heatmap_top_regulon.png')

dev.off()

####
library(ggpubr)
aucs<-getAUC(regulonAUC)
aucs<-t(aucs)
aucs<-cbind(aucs,group= as.character( Idents(sce)))
#aucs<-data.frame(aucs)
group_by(aucs, group) %>% summarize_each(funs(mean), colnames(aucs)[1:81])
aucss<-apply(aucs,2,as.numeric)
rownames(aucss)<-rownames(aucs)
aucsmenas<-aggregate(aucss[,c(1:81)],list(aucss[,82]),mean)#########根据最后确定的TF数量进行构建 heatmap
aucsmenas<-t(aucsmenas)
colnames(aucsmenas)<-aucsmenas[1,]
aucsmenas<-aucsmenas[-1,]
aucsmenas<-aucsmenas[!rownames(aucsmenas) %like% "extended",]

pdf("TF_scenic_NMF_Stromal.pdf",height = 5.84,width=3.52)
pheatmap(aucsmenas,scale = "row",
         show_colnames =T,
         show_rownames = T,
         cluster_cols = T,
         color=colorRampPalette(c("blue","white","red"))(20)
         )
dev.off()

