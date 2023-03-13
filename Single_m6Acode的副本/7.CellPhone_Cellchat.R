####input the nmfcluster foreach cell~
m6agenes<-c("CBLL1","KIAA1429","METTL14","METTL3", "RBM15","RBM15B", "WTAP",
            "ZC3H13","ELAVL1","FMR1", "HNRNPA2B1", "HNRNPC", "IGF2BP1", 
            "IGF2BP2", "IGF2BP3", "LRPPRC", "YTHDC1", "YTHDC2", "YTHDF1", 
            "YTHDF2","YTHDF3","ALKBH5","FTO")
mycol <- c("#223D6C","#D20A13","#FFD121","#088247",
           "#58CDD9","#5D90BA","#431A3D","#11AA4D","#91612D","#6E568C","#7A142C",
           "#E0367A","#D8D155","#64495D","#7CC767")
###input the phe data 
GSE132465_phe <- readRDS('GSE132465_phe.rds')#####
####mye mac
load("myeloids_NNFclusters.Rdata")
mac_cluster<-data.frame(index=Allcell_clusterss$cell_name,NMFcluster=paste0("mac_",Allcell_clusterss$NMFcluster))
#####epi cell
load("meta.data_epi_SMC.Rdata")
epi_cluster<-data.frame(index=meta.data$Index,NMFcluster=paste0("epi_",meta.data$NMFcluster))
#####CAF
load("group_CAF_SMC_NMF.Rdata")
group_CAF<-data.frame(group_CAF)
caf_cluster<-data.frame(index=rownames(group_CAF),NMFcluster=paste0("caf_",group_CAF$group_CAF))
####Tcells
load("Tcells_NNFclusters.Rdata")
rownames(Allcell_clusterss)<-Allcell_clusterss$cell_name
t_cluster<-data.frame(index=Allcell_clusterss$cell_name,NMFcluster=paste0("t_",Allcell_clusterss$NMFcluster))
##Bcells
load("group_Bcell_SMC_NMF.Rdata")
group_Bcell<-data.frame(group_Bcell)
b_cluster<-data.frame(index=rownames(group_Bcell),NMFcluster=paste0("b_",group_Bcell$group_Bcell))
###########################################################################################
NMF_cluster<-rbind(epi_cluster,caf_cluster,t_cluster,b_cluster,mac_cluster)
save(NMF_cluster,file = "SMC_NMF_clusterof_each_celltype.Rdata")
#####devtools::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)
options(stringsAsFactors = FALSE)
#####
##

load("meta_methylation_SMC.Rdata")

####################################
single_dat_sce <- readRDS('single_dat_sce_GSE132465.rds')
setwd("/home/data/vip19/myprojects")##cellchat test
###
data.input<-single_dat_sce@assays$RNA@data
data.input<-data.input[,colnames(data.input) %in% meta$Index]
single_dat_sce <- CreateSeuratObject(counts =data.input, project = "GSE144735_CAF_colon")
table(single_dat_sce$orig.ident)
single_dat_sce<-AddMetaData(single_dat_sce,metadata = meta,col.name = colnames(meta))
data.input<-single_dat_sce@assays$RNA@data
data.input<-as.matrix(data.input)
data.input[1:4,1:5]
write.table(data.input, 'cellphonedb_count3.txt', sep='\t', quote=F)
###########
colnames(meta)<-c("Cell","cell_type")
meta<-as.matrix(meta)
write.table(meta, 'cellphonedb_meta2.txt', sep='\t', quote=F, row.names=F)
#######################

sce<-single_dat_sce[,single_dat_sce$orig.ident=="SMC01.T"]

meta<-meta[meta$Index %in% colnames(sce),]
meta<-as.matrix(meta)
write.table(meta, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
