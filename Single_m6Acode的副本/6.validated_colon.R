##validation in GSE or COAD TGCA cohort
##
load("colorectal.Rdata")
load("Fbcell500.Rdata")### 选取前五十个 作为maker 高表达loading 基因~
##########
library(GSVA)
library(data.table)
library(survminer)
library(survival)
names(exp)
names(phe)
####
#for RFS
variables<-c("ACTA2_related", "NNMT_related" ,"ADAM33_related" , "ADAMTS1_related")
##
for(set in names(exp)){
myexp<-exp[[set]]
myphe<-phe[[set]]
print(set)
norm.expMat<-myexp
#gsva_es <- gsva(norm.expMat,caf.sig,method="ssgsea",abs.ranking=F,kcdf="Gaussian",ssgsea.norm=TRUE)###Array data
gsva_es <- gsva(norm.expMat,caf.sig,method="ssgsea",abs.ranking=F,kcdf="Poisson",ssgsea.norm=TRUE)#RNA-Seq 
gsva_es<-t(gsva_es)
gsva_es<-cbind(rownames(gsva_es),gsva_es)
caflist[[set]]<-gsva_es
colnames(gsva_es)[1]<-"geo_accession"
#colnames(PANCAN.survival)
rt<-merge(gsva_es,myphe,by="geo_accession")
rt$OS.time<-as.numeric(as.character(rt$OS.time))
rt$OS<-as.numeric(as.character(rt$OS))
rt<-rt[which(rt$OS!='NA'),]
if(nrow(rt)>10){
  for (sig in variables){
    rt[,sig]<-as.numeric(as.character(rt[,sig]))
    sur.cut<-surv_cutpoint(rt,time = "OS.time", event = "OS", variables = sig)
    rt$score_group=ifelse(rt[,sig]>(sur.cut$cutpoint[,1]),1,0)##bug 负数越大 反方向！！！
    #rt$score_group=rt[,sig]>median(rt[,sig])
    Gcox1<-coxph(Surv(OS.time,OS)~score_group,data=rt)
    GSum<-summary(Gcox1)
    HR<-round(GSum$coefficients[,2],3)
    Pvalue<-round(GSum$coefficients[,5],3)
    CI<-paste0(round(GSum$conf.int[,3:4],3),collapse='-')
    coeff<-round(GSum$coefficients[1],3)
    se<-GSum$coefficients[1,3]
    low<-round(GSum$conf.int[,3],3)
    up<-round(GSum$conf.int[,4],3)
    cox.p<-data.frame('Characteristics'= sig,
                      'Hazard Ratio'= HR,
                      'CI95'=CI,
                      "coeff"=coeff,
                      "se"=se,
                      "low"=low,
                      "up"=up,
                      'P-value'=Pvalue,
                      "set"=set,
                      "numberofpatients"=nrow(rt))
    Pan_cox.sig=rbind(Pan_cox.sig,cox.p)
  }
}
}