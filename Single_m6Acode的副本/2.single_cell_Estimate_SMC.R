
###
library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
####################
#################################
GSE132465_single<-fread("GSE132465_GEO_processed_CRC_10X_natural_log_TPM_matrix.txt")
GSE132465_single[1:5,1:5]
GSE132465_single<-data.frame(GSE132465_single)
rownames(GSE132465_single)<-GSE132465_single[,1]
GSE132465_single<-GSE132465_single[,-1]
# input phe data for single cell data
GSE132465_phe<-read.delim2("GSE132465_GEO_processed_CRC_10X_cell_annotation.txt")
table(GSE132465_phe$Cell_subtype)
table(GSE132465_phe$Cell_type)
GSE132465_phe$Index<-gsub("-",".",GSE132465_phe$Index)
rownames(GSE132465_phe)<-GSE132465_phe$Index
#################
my_data_frame<-GSE132465_single
# split for speed
chunk <- 1000
n <- ncol(my_data_frame)
r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
r<-data.frame(r,colnames(my_data_frame))
#d <- split(my_data_frame,r)
sample_split_list<-list()
table(r)
for ( i in unique(r$r)){
  if (!dir.exists("split_sampl_Single")){
    dir.create("./split_sampl_Single")
  }
  print(i)
  mydata<-my_data_frame[,colnames(my_data_frame) %in% r[r$r==i,]$colnames.my_data_frame.]
  sample_split_list[[i]]<-mydata
  write.table(mydata,paste0("./split_sampl_Single/", i, ".txt"),sep="\t",quote=F,col.names=T)
}
####

###############
library(estimate)
#ssgsva<-lapply(sample_split_list, function(x) return(gsva(x,gene.set.H,parallel.sz=1)))
#expMat<-as.matrix(FetchData(epidat,vars.all=epidat@var.genes))
?estimateScore
?filterCommonGenes
#write.table(data,file="CancerExpr.txt",sep="\t",quote=F,col.names=T)
estimateScores<-list()
for ( i in unique(r$r)){
print(i)
in.file<-paste0("./split_sampl_Single/", i, ".txt")
out.file<-paste0("./split_sampl_Single/", i,"outfile.gct")
out.file2<-paste0("./split_sampl_Single/", i,"out.file_score.gct")
filterCommonGenes(input.f= in.file , output.f=out.file, id="GeneSymbol")
estimateScore( out.file, out.file2)
#estimateScore[1:2,1:4]
estimateScore<-read.table(out.file2, skip = 2, header = TRUE, sep = "\t")
estimateScore<-data.frame(t(estimateScore))
colnames(estimateScore)<- c('StromalScore','ImmuneScore', 'ESTIMATEScore',"TumorPurity")
estimateScore<-estimateScore[-c(1,2),]
estimateScores[[i]]<-estimateScore
message(paste0(i,"calculation of Immunescore is done!"))
}


single_cell_estimateScores<-rbindlist(estimateScores,use.names=TRUE, fill=TRUE)
saveRDS(estimateScores,file="estimateScores.rds")
saveRDS(single_cell_estimateScores,file="single_cell_estimateScores.rds")
#write.table(estimateScore,file="TCGA_PAAD_estimate_score.txt",sep="\t",quote=F,col.names=T)
#glycolytic<-data.frame(read_excel("glycolytic.xlsx",1, col_names= TRUE, col_types= NULL, na=" "))##chip data