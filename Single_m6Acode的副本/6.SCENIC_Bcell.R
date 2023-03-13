###
library(SCENIC)
# https://resources.aertslab.org/cistarget/

dir()
getwd()
library(SCENIC)
db='cisTarget_databases/'
list.files(db)
# 保证cisTarget_databases 文件夹下面有下载好 的文件
scenicOptions <- initializeScenic(org="hgnc", 
                                  dbDir=db , nCores=4) 
#saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
#saveRDS(cellInfo, file="int/cellInfo.Rds")

##
exprMat<-as.matrix(single_dat) 
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered<-t(exprMat_filtered)
write.csv(exprMat_filtered,file="Bcells_SMC.csv")


rm(list = ls())  
library(SCENIC)
packageVersion("SCENIC")  
library(SCopeLoomR)
scenicLoomPath='sample_SCENIC_Bcell.loom'
loom <- open_loom(scenicLoomPath)
# Read information from loom file: 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")##
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 
###########################
#run Tcell SMC regulation
sce=readRDS("./TcellSMC.rds")
Idents(sce)<-sce$NMFcluster
sce
library(pheatmap) 
n=t(scale(t(getAUC(regulonAUC[,] )))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
dim(n) 
ac=data.frame(group= as.character(Idents(sce)))
n[1:4,1:4]
n=n[,colnames(n) %in% colnames(sce)]
rownames(ac)=colnames(n) ##
lengths(regulons)
library(ggpubr)
#aucs<-getAUC(regulonAUC)
aucs<-getAUC(regulonAUC)
aucs<-t(aucs)
regulonss<-regulons
for(i in c(1:198)){
  names(regulonss)[i]<-gsub("\\+",paste0(lengths(regulons)[i],"g"),names(regulons)[i])
}
names(regulonss)
colnames(aucs)<-names(regulonss)
#aucs<-t(n)
#aucs<-data.frame(aucs)
aucs<-cbind(aucs,group= as.character(Idents(sce)))
dim(aucs)
write.csv(aucs,file="Bcell_TF_idente.csv")
aucs<-read.csv("Bcell_TF_idente.csv",check.names = F)
head(aucs)
#aucs<-data.frame(aucs)
#group_by(aucs, group) %>% summarize_each(funs(mean), colnames(aucs)[1:81])
aucss<-apply(aucs,2,as.numeric)
rownames(aucss)<-rownames(aucs)
aucsmenas<-aggregate(aucss[,c(1:198)],list(aucss[,199]),mean)#########根据最后确定的TF数量进行构建 heatmap
aucsmenas<-t(aucsmenas)
colnames(aucsmenas)<-aucsmenas[1,]
aucsmenas<-aucsmenas[-c(1),]
aucsmenas[1:5,1:4]
######
data<-aucs
outTab<-c()
for (i in colnames(aucs)[c(1:198)]){
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
outTabs<-outTab[outTab$pValue<0.05,]
TF_heatmap<-aucsmenas[outTabs$gene,]
##
pdf("TF_scenic_NMF_Bcell.pdf",height = 3.24,width=3.52)
pheatmap(TF_heatmap,
         scale = "row",
         show_colnames =T,
         show_rownames = T,
         cluster_cols = F,
         fontsize = 10,
         angle_col = c('45'),
         color=colorRampPalette(c("blue","white","red"))(20)
)
dev.off()
##############################################
#library(KernSmooth)
#library(RColorBrewer)
#dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
#image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
#contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)