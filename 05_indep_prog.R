rm(list= ls())
# dir.create('./05_indep_prog')

library('forestmodel')
library("survival")
library("survminer")
library('MetBrewer')
library('ComplexHeatmap')
library('circlize')
library('ggpubr')
library('TCGAbiolinks')
##############################################################

# query <- GDCquery(
#   project = "TCGA-GBM",
#   data.category = "Clinical",
#   data.type = "Clinical Supplement",
#   data.format = "BCR Biotab"
# )
# 
# GDCdownload(query)
# 
# GDCprepare(query,save = T,save.filename = "TCGA_GBM_Clin.Rdata")
# load('TCGA_GBM_Clin.Rdata')
# clinical.biotab <- data
# names(clinical.biotab)
# 
# drug<-clinical.biotab[["clinical_drug_gbm"]]
# drug<-drug[-(1:2),] #其中有用药信息还挺好玩的，可以放进去
# 这一段后续也没有被用上，不准备放在单因素的图上面了

###############################################################
#下载maf格式的突变数据
# query <- GDCquery(
#   project = "TCGA-GBM",
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
# GDCdownload(query)
# GDCprepare(query, save = T,save.filename = "TCGA_GBM_SNP.Rdata")

# load('TCGA_GBM_SNP.Rdata')

################################################################

load('data_preparation/TCGA_data.RData')

list<-c('submitter_id.samples','gender.demographic',
        'age_at_initial_pathologic_diagnosis','karnofsky_performance_score',
        'race.demographic','GeneExp_Subtype') #

pheno<-as.data.frame(pheno)
multi<-pheno[,colnames(pheno) %in% list]
multi<-multi[,-c(7:8)]

multi$gender.demographic [multi$gender.demographic == ''] <- NA
multi$race.demographic [multi$race.demographic == '' |
                          multi$race.demographic=='not reported'] <- NA
multi$GeneExp_Subtype [multi$GeneExp_Subtype == ''] <- NA

##TCGA给的临床信息好少，我想添加一些突变进去
# mut<-data[,-1]
# mut$sample<-substr(mut$Tumor_Sample_Barcode,1,16)
mut<-mut[mut$Sample_ID %in% colnames(tcgaexp),]
test<-as.data.frame(table(mut$gene))
test<-test[order(test$Freq,decreasing = T),] #前四个吧
mutgene<-as.character(test$Var1[1:4])

cache<-as.data.frame(array(NA,c(nrow(multi),4)))
colnames(cache)<-mutgene
rownames(cache)<-multi$submitter_id.samples

for (i in 1:4){
  inf<-mut[mut$gene==colnames(cache)[i],]
  cache[,i]<-ifelse(rownames(cache) %in% inf$Sample_ID,'Mut','WT')
}

##加入用药信息
# adjname<-substr(multi$submitter_id.samples,1,12)
# drug<-drug[drug$bcr_patient_barcode %in% adjname,]
# 
# freq<-as.data.frame(table(drug$pharmaceutical_therapy_drug_name))
# 
# cache2<-as.data.frame(array(NA,c(nrow(multi),3)))
# colnames(cache2)<-c('Temodar','Avastin','Dexamethasone')
# rownames(cache2)<-multi$submitter_id.samples
# 
# list<-c('Tem','Avas|Beva','Dex')
# 
# for (i in 1:3){
#   inf<-drug[grep(list[i],drug$pharmaceutical_therapy_drug_name),]
#   cache2[,i]<-ifelse(adjname %in% inf$bcr_patient_barcode,'Treat','Un-treat')
# }
# 

##
identical(multi$submitter_id.samples,rownames(cache))
multi<-cbind(multi,cache) #,cache2

risk<-read.csv("./03_Prognosis/risk_score.csv")

multi<-multi[match(risk$X,multi$submitter_id.samples),]
identical(risk$X,multi$submitter_id.samples)

multi$riskScore<-risk$riskScore
multi$OS<-risk$OS
multi$OS.time<-risk$OS.time

colnames(multi)[c(1:5)]<-c('sampleID','Age','K_score','Gender','Race')

#COX，单因素
rownames(multi)<-multi$sampleID
multi<-multi[,-1]

{
  multi$TTN<-factor(multi$TTN,levels = c('WT','Mut'))
  multi$MUC16<-factor(multi$MUC16,levels = c('WT','Mut'))
  multi$TP53<-factor(multi$TP53,levels = c('WT','Mut'))
  multi$EGFR<-factor(multi$EGFR,levels = c('WT','Mut'))
  # multi$Temodar<-factor(multi$Temodar,levels = c('Un-treat','Treat'))
  # multi$Avastin<-factor(multi$Avastin,levels = c('Un-treat','Treat'))
  # multi$Dexamethasone<-factor(multi$Dexamethasone,levels = c('Un-treat','Treat'))
  multi$GeneExp_Subtype<-factor(multi$GeneExp_Subtype,levels = c('Proneural','Neural','Classical','Mesenchymal'))
}

vars_for_table<-colnames(multi)[1:(ncol(multi)-2)]
univ_formulas <- sapply(vars_for_table, function(x) as.formula(paste('Surv(OS.time, OS)~', x))) #单次的cox回归分析循环
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = multi)}) #把循环分析的结果整合在一起

pdf('./05_indep_prog/sig_COX.pdf',width =10,height = 8)
print(forest_model(model_list = univ_models,covariates = vars_for_table,merge_models =T)) 
dev.off()

#多因素

# res.cox<-coxph(Surv(OS.time,OS)~Age+TP53+strata(Temodar)+riskScore, data=multi)
res.cox2<-coxph(Surv(OS.time,OS)~Age+TP53+riskScore, data=multi)

# multiCox = step(res.cox, direction = "both")
# summary(res.cox)
summary(res.cox2)
# summary(multiCox)

test.ph <- cox.zph(res.cox2)
test.ph

pdf('./05_indep_prog/ph_Assumption.pdf',width =8,height = 8)
ggcoxzph(test.ph)
dev.off()

pdf('./05_indep_prog/mul_COX.pdf',width =10,height = 6)
forest_model(res.cox2,limits=c(-1,1.5))
dev.off()
##########################################
#列线图

library(rms)

dd<-datadist(multi)
options(datadist="dd")
options(na.action="na.delete")

coxpbc<-cph(formula = Surv(OS.time,OS) ~  Age+TP53+riskScore,data=multi,
            x=T,y=T,surv = T,na.action=na.delete)

print(coxpbc)
surv<-Survival(coxpbc) 
surv3<-function(x) surv(365,x)
surv4<-function(x) surv(730,x) #如果时间超过数据里最高的生存时间就会报错
surv5<-function(x) surv(1095,x)

x<-nomogram(coxpbc,fun = list(surv3,surv4,surv5),lp=T,
            funlabel = c('1-year survival Probability','2-year survival Probability','3-year survival Probability'),
            maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))

pdf("./05_indep_prog/nomogram.pdf",width = 12,height = 10)
plot(x, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()

#校准曲线

f1<-cph(formula = Surv(OS.time,OS) ~  Age+TP53+riskScore,data=multi,
        x=T,y=T,surv = T,na.action=na.delete,time.inc = 365)
cal1<-calibrate(f1, cmethod="KM", method="randomization",u=365,m=50,B=1000) 


f3<-cph(formula = Surv(OS.time,OS) ~  Age+TP53+riskScore,data=multi,
        x=T,y=T,surv = T,na.action=na.delete,time.inc = 730) 
cal3<-calibrate(f3, cmethod="KM", method="randomization",u=730,m=50,B=1000)


f5<-cph(formula = Surv(OS.time,OS) ~  Age+TP53+riskScore,data=multi,
        x=T,y=T,surv = T,na.action=na.delete,time.inc = 1095) 
cal5<-calibrate(f5, cmethod="KM", method="randomization",u=1095,m=50,B=1000)


#
pdf("./05_indep_prog/calibration_compare.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 1,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal3,lwd = 2,lty = 1,errbar.col = c("#f7aa58"),
     xlim = c(0,1),ylim= c(0,1),col = c("#f7aa58"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#f7aa58"), pch = 16)

plot(cal5,lwd = 2,lty = 1,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","2-year",'3-year'), #图例文字
       col =c("#2166AC",'#f7aa58',"#B2182B"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
dev.off()

#############################################################
#DCA

library(ggDCA)

m1 <- cph(formula = Surv(OS.time,OS) ~  TP53,data=multi,
          x=T,y=T,surv = T,na.action=na.delete)
m2 <- cph(formula = Surv(OS.time,OS) ~  riskScore,data=multi,
          x=T,y=T,surv = T,na.action=na.delete)
m3 <- cph(formula = Surv(OS.time,OS) ~  TP53++riskScore,data=multi,
          x=T,y=T,surv = T,na.action=na.delete)


data  <- dca(m1,m2,m3,
             model.names =c('TP53','riskScore','TP53+riskScore'))


pdf("./05_indep_prog/DCA.pdf",width = 10,height = 8)
ggplot(data,linetype = T,lwd=1.5,
       color = c(MetBrewer::met.brewer("Manet",n=5)))+
  coord_fixed() +theme_classic(base_size = 18)+theme(aspect.ratio=0.85)
dev.off()


##################################################################
#比较风险分数和各个临床特征之间的相关
#先把几个基因提出来
load('./03_Prognosis/model.RData')

tcgaexp<-((2^tcgaexp)-1)
rpkmTOtpm <- function(rpkm){
  exp(log(rpkm) - log(sum(rpkm)) + log(1e6))
}
tcgaexp <- apply(tcgaexp, 2, rpkmTOtpm)
exp<-log2(tcgaexp+1)
exp<-as.data.frame(exp)

gene<-c('UPP1','PLOD3','GALNT11')
pheat<-exp[rownames(exp) %in% gene,]

multi<-multi[order(multi$riskScore),]
pheat<-t(pheat)
pheat<-pheat[match(rownames(multi),rownames(pheat)),]

mat_scaled = t(scale(pheat))

col_fun1 = colorRamp2(c(-1,0,1), c(met.brewer("Hiroshige")[9],'white',met.brewer("Signac")[4]))
col_fun2 = colorRamp2(c(0,40,80), c(met.brewer("Hiroshige")[9],'white',met.brewer("Signac")[4]))
#一堆乱七八糟的注释
ha <- HeatmapAnnotation(
  riskScore = c(multi$riskScore),
  Age = c(multi$Age),
  Gender = c(multi$Gender),
  Race = c(multi$Race),
  TTN = c(multi$TTN),
  MUC16 = c(multi$MUC16),
  TP53 = c(multi$TP53),
  EGFR = c(multi$EGFR),
  'GeneExp subtype'=c(multi$GeneExp_Subtype),
  col = list(
    Age = col_fun2,
    Gender = c("male" = "red", "female" = "blue")),
  na_col = "black"
)

pdf('./05_indep_prog/complexheatmap.pdf',width = 14,height = 6)
draw(
  Heatmap(mat_scaled,
          col = col_fun1,name = "Score",
          column_title = NULL,
          # column_gap = unit(c(0, 0), "mm"),
          cluster_columns = F,
          cluster_rows =T,
          show_row_dend = T,
          show_column_names =F,
          show_heatmap_legend = T,
          row_names_side = "right",
          width = ncol(mat_scaled)*unit(1, "mm"), 
          height = nrow(mat_scaled)*unit(16, "mm"),
          top_annotation = ha))
dev.off()

#############################################################
#这个做不出来可以用类似于卡方检验的方法来做
#画箱线
multi$group<-ifelse(multi$Age>median(multi$Age),'Old','Younger')
multi$group<-factor(multi$group,levels = c('Younger','Old'))

ggplot(multi, aes(x=group, y=riskScore))+
  geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
  labs(title="GSE13159",y="riskScore")+ 
  theme_classic(base_size = 22)+ 
  scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,3)])+
  theme(legend.position = "none")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())+
  xlab('')+
  stat_compare_means(label.x = 1,size=5)+
  theme(aspect.ratio=1.5)

multi$group<-multi$Gender
  
ggplot(multi, aes(x=group, y=Age))+
    geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
    stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
    stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
    stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
    labs(title="GSE13159",y="riskScore")+ 
    theme_classic(base_size = 22)+ 
    scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,3)])+
    theme(legend.position = "none")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank())+
    xlab('')+
    stat_compare_means(label.x = 1,size=5)+
    theme(aspect.ratio=1.5)
  
multi$group<-multi$Race
multi<-multi[!is.na(multi$group),]
  
  ggplot(multi, aes(x=group, y=Age))+
    geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
    stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
    stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
    stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
    labs(title="GSE13159",y="riskScore")+ 
    theme_classic(base_size = 22)+ 
    scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,3,7)])+
    theme(legend.position = "none")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank())+
    xlab('')+
    stat_compare_means(label.x = 1,size=5)+
    theme(aspect.ratio=1.5)

multi$group<-multi$TTN
multi<-multi[!is.na(multi$group),]
  
  ggplot(multi, aes(x=group, y=Age))+
    geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
    stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
    stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
    stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
    labs(title="GSE13159",y="riskScore")+ 
    theme_classic(base_size = 22)+ 
    scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,3)])+
    theme(legend.position = "none")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank())+
    xlab('')+
    stat_compare_means(label.x = 1,size=5)+
    theme(aspect.ratio=1.5)
  
multi$group<-multi$MUC16
multi<-multi[!is.na(multi$group),]
  
ggplot(multi, aes(x=group, y=Age))+
    geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
    stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
    stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
    stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
    labs(title="GSE13159",y="riskScore")+ 
    theme_classic(base_size = 22)+ 
    scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,3)])+
    theme(legend.position = "none")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank())+
    xlab('')+
    stat_compare_means(label.x = 1,size=5)+
    theme(aspect.ratio=1.5)

multi$group<-multi$TP53
multi<-multi[!is.na(multi$group),]
  
ggplot(multi, aes(x=group, y=Age))+
    geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
    stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
    stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
    stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
    labs(title="GSE13159",y="riskScore")+ 
    theme_classic(base_size = 22)+ 
    scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,3)])+
    theme(legend.position = "none")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank())+
    xlab('')+
    stat_compare_means(label.x = 1,size=5)+
    theme(aspect.ratio=1.5)  
  
multi$group<-multi$EGFR
multi<-multi[!is.na(multi$group),]
  
ggplot(multi, aes(x=group, y=Age))+
    geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
    stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
    stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
    stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
    labs(title="GSE13159",y="riskScore")+ 
    theme_classic(base_size = 22)+ 
    scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,3)])+
    theme(legend.position = "none")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank())+
    xlab('')+
    stat_compare_means(label.x = 1,size=5)+
    theme(aspect.ratio=1.5)  
  
#一个相关的特征都没有吗？
library('ggstatsplot')
  
#需要把数据变成两列
multi$group<-ifelse(multi$riskScore>median(multi$riskScore),'risk_High','risk-Low')


ggbarstats(multi,x = group,y = TTN)
ggbarstats(multi,x = group,y = MUC16)
ggbarstats(multi,x = group,y = TP53)
ggbarstats(multi,x = group,y = EGFR)

ggbarstats(multi,x = group,y = Gender)
ggbarstats(multi,x = group,y = Race)

multi$group2<-ifelse(multi$Age>median(multi$Age),'Old','Younger')
ggbarstats(multi,x = group,y = group2)

ggbarstats(multi,x = group,y = Temodar)
ggbarstats(multi,x = group,y = Avastin)
ggbarstats(multi,x = group,y = Dexamethasone)

#这三个药也都没有相关性

plot<-multi[!is.na(multi$K_score),]
ggplot(plot, aes(x=group, y=K_score))+
  geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
  labs(y="riskScore")+ 
  theme_classic(base_size = 22)+ 
  scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,3)])+
  theme(legend.position = "none")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())+
  xlab('')+
  stat_compare_means(label.x = 1,size=5)+
  theme(aspect.ratio=1.5) 

plot<-multi[!is.na(multi$GeneExp_Subtype),]

pdf('./05_indep_prog/score_with_subtype.pdf',width = 8,height = 8)
ggplot(plot, aes(x=GeneExp_Subtype, y=riskScore))+
  geom_boxplot(width=0.1,aes(fill=GeneExp_Subtype),colour='black',alpha = 0.4,linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=GeneExp_Subtype),color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
  labs(y="riskScore")+ 
  theme_classic(base_size = 22)+ 
  scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,5,3,2)])+
  theme(legend.position = "none")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())+
  xlab('')+
  stat_compare_means(label.x = 1,size=5)+
  theme(aspect.ratio=0.8) 
dev.off()

