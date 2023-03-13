##immune score
###prognisis  of pan cancer
library(data.table)
library(GEOquery)
library(dplyr)
library(Seurat)
library(ggplot2)
library(umap)
library(hypeR)
library(msigdbr)
library(tidyverse)
library(pheatmap)
library(Seurat)
library(clustree)
library(survival)
library(survminer)
library(NMF)##
load("CRC_immunecell_signature.Rdata")
load("CRC_cell_subtype_signatures.Rdata")
load("Colon_rectal_2021.Rdata")
load("cell_sutypelist_signature_Scoreforbulksequence.Rdata")###
##################
load("~/Projects/PanNNMT/PanCancerdata/immunecohorts_allRNASeqandPhe20211223.Rdata")###
save(immunecohorts,immunephe,file = "~/Projects/PanNNMT/PanCancerdata/immunecohorts_allRNASeqandPhe20211223.Rdata")
cell_sutypelist<-c() 
for (set in names(cell_subtype_signatures)) {
  print(set)
  celltype<-cell_subtype_signatures[[set]]
  for (subtype in names(celltype)){
    print(subtype)
    mySubcelltype<-celltype[[subtype]]
    if(length(mySubcelltype)>1){
      cell_sutypelist[[subtype]]<-mySubcelltype
    }
  }
}###
##
#######
library(GSVA)
####################
names(immunecohorts)
delet<-c("Chen2016_CTLA4_Melanoma","GSE140901_HCC_2021","GSE93157_Lung_PD1_NIVOLUMAB_2017" )
immune_methylationscore<-list()
names(cell_sutypelist)<-gsub(" ","_",names(cell_sutypelist))
names(cell_sutypelist)<-gsub("\\+","_",names(cell_sutypelist))
###################################################
for(set in names(immunecohorts)){
  myexp<-immunecohorts[[set]]
  #myphe<-phe[[set]]
  print(set)
  norm.expMat<-as.matrix(myexp)
  #gsva_es <- gsva(norm.expMat,cell_sutypelist,method="ssgsea",abs.ranking=F,kcdf="Gaussian",parallel.sz=50)###Array data,ssgsea.norm=TRUE
  gsva_es <- gsva(norm.expMat,cell_sutypelist,method="ssgsea",abs.ranking=F,kcdf="Poisson",parallel.sz=40)#RNA-Seq 
  gsva_es<-t(gsva_es)
  gsva_es<-cbind(rownames(gsva_es),gsva_es)
  immune_methylationscore[[set]]<-gsva_es
  message(set,"GSVA is done!")
}

b<-c()
for (set in names(immunephe)){
  myphe<-immunephe[[set]]
  a<-ifelse("Cat_Response" %in% colnames(myphe),"Cat_Response","NO")
  a<-cbind(a,set)
  b<-rbind(b,a)
}##
b
set<-"Braun_RCC_Nat_medicine2020"

phe<-immunephe[[set]]

rownames(phe)<-phe$patient
immunephe[[set]]<-phe

Univcoxss<-list()
for (set in names(immune_methylationscore)[-c(2,5,10,13,16,18)]){
  print(set)
  rt<-immune_methylationscore[[set]]
  rt<-cbind(sample=rownames(rt),rt)
  phe<-immunephe[[set]]
  phe<-cbind(sample=rownames(phe),phe)
  rt<-merge(phe,rt,by="sample")
  Univcoxs<-c()
  Varname<-names(cell_sutypelist)
  for(index in Varname){
    print(index)
    if(index %in% colnames(rt)){
    rt[,index]<-as.numeric(rt[,index])
    rt$groups<- ifelse(rt[,index]> median(rt[,index]),"High","Low")
    fml<-as.formula(paste0("Cat_Response~","groups"))
    logist<- glm(fml,family=binomial(link='logit'), data =rt)
    p_value<-round(summary(logist)$coefficients[8], 3)
    coeff<-round(summary(logist)$coefficients[2], 3)
    stderr<-round(summary(logist)$coefficients[4], 3)
    or<- round(exp(summary(logist)$coefficients[2]),3)
    up<-round(exp(coeff+1.96*stderr),3)
    low<-round(exp(coeff-1.96*stderr),3) 
    CI95<-paste0(low,'-',up)
    univlogistic<-data.frame('Characteristics'=index,
                             'Odds Ratio'= or,
                             'CI95'= CI95,
                             'P-value'=p_value)
    Univcoxs<-rbind(Univcoxs,univlogistic)
    }
  }
  Univcoxs<-cbind(Univcoxs,set)
  Univcoxss[[set]]<- Univcoxs
 }
Univcoxss2<-do.call(rbind,Univcoxss)
####
write.table(Univcoxss2, file="Univcoxss_pathway.txt",sep = "\t")


OSdata<-read.delim2("Univcoxss_pathway.txt")
library(readxl)
colnames(OSdata)

library(ggplot2)
library(data.table)
library(ggpubr)
library(ggsci)

darkblue <- "#303B7F"
darkred <- "#D51113"
yellow <- "#EECA1F"

mid<-0

my_palette <- colorRampPalette(c(darkblue,yellow,darkred), alpha=TRUE)(n=128)

table(OSdata$Characteristics)

OSdata$OR<-ifelse(OSdata$OR>2,2,OSdata$OR)
OSdata$P.value<-as.numeric(OSdata$P.value)
OSdata$OR<-as.numeric(OSdata$OR)


OSdata %>%
  mutate(Characteristics=fct_relevel(Characteristics,c("Fibroblast",
                                                       "HNRNPA2B1+CAF-C1",
                                                       "WTAP+CAF-C2",
                                                       'HNRNPC+CAF-C3',
                                                       'NoneMethy-CAF-C4',
                                                       "Macrophage",
                                                       'WTAP+mac-C1',
                                                       'HNRNPA2B1+mac-C3',
                                                       'YTHDC1 &YTHDF3+mac_C4',
                                                       "Regulatory T cells",
                                                       'HNRNPA2B1+Reg T cells-C1',
                                                       'YTHDF2+Reg T cells-C4',
                                                       "CD8+ T cells",
                                                       'HNRNPA2B1+CD8+ T cells-C1',
                                                       'ALKBH5+CD8+ T cells-C4',
                                                       "B cells",
                                                       'HNRNPA2B1+B-C1',
                                                       'WTAP+B-C2',
                                                       'NoneMethy+B-C4'))) %>%
  ggplot(aes(x=Characteristics,y=Set)) +
  geom_point(aes(size=-log10(P.value),color=OR)) +
  geom_point(shape=21,aes(size=-log10(P.value)),position =position_dodge(0))+
  scale_color_gradientn('OR',colors=my_palette) + 
  #scale_color_gradient2(midpoint=1,low="#303B7F",mid="#EECA1F",high="#D51113" )+
  theme_bw() +
  #coord_flip()+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8, hjust = 0.95, vjust = 0.2, color = "black"),###对齐  超级棒  卧槽  ！！！！！！！
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
        legend.position = "right",
        plot.margin = unit(c(1,1,1,1), "lines"))
ggsave("prognosis_of_each_type_responseImmunotherapy.pdf",width = 9.38, height=4.3)
