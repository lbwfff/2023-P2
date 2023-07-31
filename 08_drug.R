setwd('C:/Users/Administrator/Desktop/project/YQDZ-0506-6')
rm(list= ls())

library(oncoPredict)
# dir.create('./08_drug')

############################################################
#最后是两个药物的分析

load('data_preparation/TCGA_data.RData')
risk<-read.csv("./03_Prognosis/risk_score.csv")

save(tcgacount,file = 'exp.RData')

dir='G:/oncopredict/Training Data'
#oncopredict提供的训练集，有GDCS和CTRP两个数据库的数据集

GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))  #什么叫RMA？
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 

GDSC2_Res<-GDSC2_Res

testExpr<- as.matrix(tcgacount)[,1:2]
# testExpr[1:4,1:4]  
# colnames(testExpr)=paste0('test',colnames(testExpr))
dim(testExpr) 

set.seed(23719)
oncoPredict::calcPhenotype(trainingExprData = GDSC2_Expr,
                           trainingPtype = GDSC2_Res,
                           testExprData = testExpr, #目标数据集
                           batchCorrect = 'eb',  #   "eb" for ComBat  
                           powerTransformPhenotype = TRUE,
                           removeLowVaryingGenes = 0.2,
                           minNumSamples = 10, 
                           printOutput = TRUE, 
                           removeLowVaringGenesFrom = 'rawData' )

result<-read.csv('./08_drug/calcPhenotype_Output/DrugPredictions.csv')

result$X<-substr(result$X,5,999)

result<-result[result$X %in% risk$X,]

result<-result[match(risk$X,result$X),]
result$risk<-risk$riskScore

cor<-cor(result[,-1])
cor<-as.data.frame(cor)
cor$drug<-rownames(cor)
cor<-cor[cor$risk>0.3 | cor$risk<(-0.3),]

list<-rownames(cor)

corlist<-data.frame(drug=c(rownames(cor)[-nrow(cor)]),
                    p=c(NA),
                    r=c(NA))

for (i in 1:nrow(corlist)){
  inf<-result[,colnames(result) %in% c(corlist$drug[i],'risk')]
  colnames(inf)<-c('drug','exp')
  corlist$p[i]<-cor.test(inf$drug,inf$exp)$p.value
  corlist$r[i]<-cor.test(inf$drug,inf$exp)$estimate
}

#还挺多相关的，要怎么做展示呢？
corlist$name<-c(NA)

for(i in 1:nrow(corlist)){
  name<-corlist$drug[i]
  name<-unlist(strsplit(name,'_'))
  name<-name[-length(name)]
  name<-paste0(name,collapse = '_')
  corlist$name[i]<-c(name)
}

corlist<-corlist[order(corlist$p),]

write.csv(corlist,file = './08_drug/oncoPredict_cortest.csv')
#
pdf("./08_drug/score_oncoPredict_cor.pdf",width = 8,height = 6)
ggplot(corlist[1:10,],aes(x=reorder(name,r),y=r))+
  geom_point()+
  geom_segment(aes(x=name,xend=name,y=0,yend=r),size=1.5,color="#C9CACA",linetype="solid")+
  geom_point(size=5,aes(color=-log10(p),fill=-log10(p)),shape=21)+
  theme_test(base_size = 22)+theme(panel.grid.major.x=element_blank(),
                                   panel.border=element_blank(),
                                   axis.ticks.y=element_blank())+
  xlab("")+ylab("R value")+scale_y_continuous(breaks = seq(-3, 3, 0.1))+
  coord_flip()+theme(aspect.ratio=1)

dev.off()

#批量画个箱线
library('ggplot2')
library('ggpubr')
library('MetBrewer')

plot<-result
rownames(plot)<-plot$X
plot<-plot[,-1]
plot<-plot[,colnames(plot) %in% c(corlist$drug[1:10],'risk')]
plot<-plot[,match(c(corlist$drug[1:10],'risk'),colnames(plot))]

p=list() 

for (i in 1:(ncol(plot)-1)){
  
  test<-plot[,c(i,ncol(plot))]
  colnames(test)<-c('IC50','Risk')
  test$group<-ifelse(test$Risk>median(test$Risk),'risk_High','risk_Low')
  
  p[[i]]<-
    ggplot(test, aes(x=group, y=IC50))+
    geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
    stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
    stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
    stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
    labs(title=paste0(corlist$name[i]),y="IC50")+ 
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
}

library(patchwork)
pdf("./08_drug/drug_box.pdf",width = 18,height = 12)
wrap_plots(p,nrow=2, guides="collect") 
dev.off()




