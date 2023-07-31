setwd('C:/Users/Administrator/Desktop/project/YQDZ-0506-6')
rm(list= ls())
library('forestmodel')
library('survminer')

# dir.create('03_Prognosis')


################################################################################
#只有20个，感觉很危险啊，试试这20个里面能不能有预后相关的
#先是单因素
genelist<-read.csv('02_genescreen/genelist.csv')
genelist<-genelist$x
load('data_preparation/TCGA_data.RData')

tcgaexp<-((2^tcgaexp)-1)
rpkmTOtpm <- function(rpkm){
  exp(log(rpkm) - log(sum(rpkm)) + log(1e6))
}
tcgaexp <- apply(tcgaexp, 2, rpkmTOtpm)
exp<-log2(tcgaexp+1)

suvexp<-exp[,colnames(exp) %in% suv$sample]
suvexp<-suvexp[rownames(suvexp) %in% genelist,]
suv<-suv[match(colnames(suvexp),suv$sample),]
net<-as.data.frame(suv[,c(2,4)])
rownames(net)<-suv$sample

library(survival)
library(randomForestSRC)

identical(rownames(net),colnames(suvexp))
pbc=cbind(net,t(suvexp))

vars_for_table<-c(genelist)
univ_formulas <- sapply(vars_for_table, function(x) as.formula(paste('Surv(OS.time, OS)~', x))) #单次的cox回归分析循环
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = pbc)}) #把循环分析的结果整合在一起

print(forest_model(model_list = univ_models,covariates = vars_for_table,
                   merge_models =T)) #绘图，forest_model做出来这个还行

#整理为表格
table<-data.frame(gene=c(genelist),
                  p=c(NA),
                  hz=c(NA))
for (i in 1:length(genelist)){
  table$p[i]<-summary(univ_models[[i]])$coefficients[, 5]
  table$hz[i]<-exp(coef(univ_models[[i]]))
}

write.csv(table,file = './03_Prognosis/all_gene_cox.csv')

filter<-table[table$p<0.05,]
filter<-filter[order(filter$hz,decreasing = T),]

vars_for_table<-c(filter$gene)
univ_formulas <- sapply(vars_for_table, function(x) as.formula(paste('Surv(OS.time, OS)~', x))) #单次的cox回归分析循环
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = pbc)}) #把循环分析的结果整合在一起

pdf('./03_Prognosis/filter_forest.pdf',width =8,height = 6)
forest_model(model_list = univ_models,covariates = vars_for_table,
             merge_models =T)
dev.off()

#############################################################
#11个单因数基因后接随机生存森林
pbc2<-pbc[,colnames(pbc) %in% c(colnames(pbc)[1:2],filter$gene)]
save(pbc2,file = 'pbc.RData')
dat_lasso<-pbc2

cox.weights <- function(rfsrc.f, rfsrc.data) {
  event.names <- all.vars(rfsrc.f)[1:2]
  p <- ncol(rfsrc.data) - 2
  event.pt <- match(event.names, names(rfsrc.data))
  xvar.pt <- setdiff(1:ncol(rfsrc.data), event.pt)
  sapply(1:p, function(j) {
    cox.out <- coxph(rfsrc.f, rfsrc.data[, c(event.pt, xvar.pt[j])])
    pvalue <- summary(cox.out)$coef[5]
    if (is.na(pvalue)) 1.0 else 1/(pvalue + 1e-100)
  })
}

data_rf=dat_lasso  #输入数据
rfsrc.f <- as.formula(Surv(OS.time, OS) ~ .)  #输入y，记得改PFS和PFS_T的列名
cox.wts <- cox.weights(rfsrc.f, data_rf)

set.seed(1256)
b.obj <- rfsrc(rfsrc.f, data_rf , nsplit = 10, xvar.wt = cox.wts,
               importance = T,na.action ="na.impute",ntree = 1000)

pdf('03_Prognosis/rfsrc_plot.pdf',width = 12,height = 6)
plot(b.obj)
dev.off()

vimp=as.data.frame(vimp(b.obj)$importance)
vimp$gene<-rownames(vimp)
vimp$vri<-(abs(vimp$`vimp(b.obj)$importance`)/max(vimp$`vimp(b.obj)$importance`))

# jk.obj <- subsample(b.obj)
# plot(jk.obj)

vimp<-vimp[order(vimp$vri,decreasing = T),]

#画一下相对变量强度的图
pdf('03_Prognosis/rfsrc_Importane_plot.pdf',width = 6,height = 4)
ggplot(data = vimp, aes(x = reorder(gene,vri), y = vri)) + 
  geom_col(aes(fill='#d8443c',colour='#d8443c'),show.legend = FALSE,width=0.8) +
  coord_flip() +
  theme_classic(base_size = 13)+ 
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        text = element_text(size=12,face='bold'))+
  labs(y="Variable Relative Importane",x=NULL)+
  scale_fill_manual(values = c('#0073C2FF'))+
  scale_colour_manual(values = c("black"))+
  scale_y_continuous(expand = c(0,0))+
  theme(plot.title = element_text(size=12))
dev.off()

list<-vimp$gene[vimp$vri>0.25]
###########################################################
#对于这三个基因画相关性的热图
corplot<-as.data.frame(t(exp[rownames(exp) %in% list,]))

M<-cor(corplot)
# cor_p <- corr.p(M,n = nrow(corplot),alpha=.05)
# M[lower.tri(cor_p$p)] <- cor_p$p[lower.tri(cor_p$p)]

pdf(file="./03_Prognosis/3gene_corrplot.pdf",width=6,height=6,onefile=FALSE)
corrplot.mixed(M,lower = "number",upper = "circle",upper.col=rev(COL2('RdBu', 200)),
               tl.pos = "lt",lower.col = "black")
dev.off()

###########################################################
library(glmnet)

lassoarray<-t(suvexp[rownames(suvexp) %in% list,])
net<-net[,c(2,1)]
colnames(net)<-c('time','status')

set.seed(23725)

cv<-cv.glmnet(as.matrix(lassoarray),as.matrix(net),nfold=10,family='cox')
fit.train <- cv$glmnet.fit
pdf(file = "./03_Prognosis/lasso1.pdf",width =8,height = 6)
plot(cv,las=1)
dev.off()

fit <- glmnet(as.matrix(lassoarray),as.matrix(net),family = "cox")
pdf(file = "./03_Prognosis/lasso2.pdf",width =8,height = 6)
plot(fit,xvar = "lambda",label = TRUE, las=1)
dev.off()

##

coef.min = coef(cv, s = "lambda.min")  
active.min = which(coef.min@i != 0) ## 找出那些回归系数没有被惩罚为0的
lasso_suv_geneids <- coef.min@Dimnames[[1]][coef.min@i+1] 

lasso_coef<-data.frame(gene=c(lasso_suv_geneids),
                       coef=c(coef.min@x))

###

pbc<-as.data.frame(pbc)
cox.data = pbc[,colnames(pbc) %in% c('OS','OS.time',lasso_coef$gene)]
cox.data.step <- na.omit(cox.data[, c("OS.time", "OS", colnames(cox.data)[!(colnames(cox.data) %in% c('OS','OS.time'))])])

#LASSO版本得分
lasso.prob<-predict(cv,newx = as.matrix(lassoarray),s=c(cv$lambda.min,cv$lambda.1se))
lasso.prob<-lasso.prob[match(rownames(cox.data.step),rownames(lasso.prob)),]
cox.data.plot<-cbind(cox.data.step, lasso.prob[,1])
colnames(cox.data.plot)[ncol(cox.data.plot)]<-c('riskScore')

write.csv(cox.data.plot, file = "./03_Prognosis/risk_score.csv")
save(cv,file = './03_Prognosis/model.RData')

#########################################################
#模型在训练集的评估
#时间依赖ROC

library(survivalROC)
heatmap_train <- cox.data.plot

nobs<-nrow(heatmap_train)

pdf(file="./03_Prognosis/time_ROC.pdf",height = 6,width = 6)
roc=survivalROC(Stime=heatmap_train$OS.time, status=heatmap_train$OS,span = nrow(heatmap_train)^-1,
                marker = heatmap_train$riskScore,predict.time =1*365)
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='#FFC1C1',
     xlab="False positive rate", ylab="True positive rate",
     #main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     main="ROC curve",
     lwd = 3, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
aucText=c()
rocCol <- c('#FFC1C1','#63B8FF','#FA8072','#ADFF2F','#FFFF00')
aucText=c(aucText,paste0("1 years"," (AUC=",sprintf("%.3f",roc$AUC),")"))
j =0
for (i in c(2,3)){
  roc1=survivalROC(Stime=heatmap_train$OS.time, status=heatmap_train$OS,span = nrow(heatmap_train)^-1,
                   marker = heatmap_train$riskScore,predict.time =i*365)
  j=j+1
  aucText=c(aucText,paste0(i," years"," (AUC=",sprintf("%.3f",roc1$AUC),")"))
  lines(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j+1],lwd = 3)
}
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex = 1)
abline(0,1)
dev.off()

#生存
#生存曲线
matt<-cox.data.plot
med.exp<-median(matt$riskScore)
more.med.exp.index<-which(matt$riskScore>=med.exp)
less.med.exp.index<-which(matt$riskScore< med.exp)
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(OS.time/30,OS) ~ status, data = matt)
s.diff<-survdiff(Surv(OS.time/30,OS) ~ status, data = matt)

sdata.plot3<-ggsurvplot(s.fit, data=matt,
                        palette="Pastel1",
                        pval = TRUE,pval.method = TRUE, conf.int = TRUE,
                        xlab = 'Time (Month)',ggtheme = theme_survminer(),
                        surv.median.line = 'hv',size=2,
                        title=paste0("TCGA survival"))

pdf("./03_Prognosis/survplot.pdf",width = 8,height = 6)
sdata.plot3$plot+
  theme_classic(base_size = 22)+ 
  theme(legend.position = "right")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())+
  theme(aspect.ratio=1)
dev.off()
#####################################################
#曲线

plot<-matt[order(matt$riskScore),]
plot$patientid<-1:nrow(plot)

library(ggplot2)
p1=ggplot(plot,aes(x=patientid,y=riskScore))+geom_point(aes(color=status))+
  scale_colour_manual(values = c("#990033","#003366"))+
  theme_bw()+labs(x=" ",y="Risk score")+theme(legend.position = "top")+
  geom_hline(yintercept=median(plot$riskScore),colour="black", linetype="dotted",size=0.8)+
  geom_vline(xintercept=sum(plot$status==names(table(plot$status))[1]),colour="black", linetype="dotted",size=0.8)+
  scale_x_continuous(expand = c(0,0))+  #把X轴和Y轴的间隙去除
  scale_y_continuous(expand = c(0,0))
p1

#散点
plot$OS<-factor(plot$OS,levels = c(1,0))
p2=ggplot(plot,aes(x=patientid,y=OS.time/365))+geom_point(aes(col=OS))+theme_bw()+
  scale_colour_manual(values = c("#990033","#003366"))+theme(legend.position = "top")+
  labs(x="Patient ID(increasing risk score)",y="Survival time(year)")+
  geom_vline(xintercept=sum(plot$status==names(table(plot$status))[1]),colour="black", linetype="dotted",size=0.8)+
  scale_x_continuous(expand = c(0,0))+  #把X轴和Y轴的间隙去除
  scale_y_continuous(expand = c(0,0))
p2

#热图
library(pheatmap)
mycolors <- colorRampPalette(c("white", "#003366", "#990033"), bias = 1.2)(100)

exp_dat<-plot[,colnames(plot) %in% lasso_coef$gene]
tmp=t(scale(exp_dat))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1

lasso_coef<-lasso_coef[order(lasso_coef$coef),]
tmp<-tmp[match(lasso_coef$gene,rownames(tmp)),]

test <- tmp %>% as.data.frame()%>% dplyr::mutate(B=row.names(.)) %>% reshape::melt()

test$B<-factor(test$B,levels = rev(lasso_coef$gene))
test$variable<-factor(test$variable,levels = colnames(tmp))

p3<-ggplot(test,aes(x=variable,y=B,fill=value))+
  geom_raster()+
  scale_fill_gradient2(low="#003366", high="#990033", mid="white")+
  labs(x=NULL, y = NULL)+
  theme(axis.text.x = element_blank(),axis.ticks=element_blank(),
        legend.position = "none")


#拼图

pdf("./03_Prognosis/riskscore.pdf",width = 8,height = 8)
cowplot::plot_grid(p1,p2,p3,ncol=1,align = "hv")
dev.off()



