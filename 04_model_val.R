setwd('C:/Users/Administrator/Desktop/project/YQDZ-0506-6')
rm(list= ls())
library('forestmodel')
library('survminer')

# dir.create('04_model_val')
load('./03_Prognosis/model.RData')
################################################################
#先等方案老师的回答，先进行独立预后
#看看验证集怎么样
#GSE13041

load('./data_preparation/GSE13041_ready.RData')

genes_expr<-limma::normalizeBetweenArrays(exp_1,method="quantile") #前期已经做了分位数中心化
genes_expr<-log2(genes_expr+1)

array<-as.data.frame(
  t(genes_expr[rownames(genes_expr) %in% names(cox$coefficients),]))

#merge array
meta<-meta1[,c(2,63,62)]
colnames(meta)<-c('name','OS','OS.time')
meta<-meta[!is.na(meta$OS),]
meta$OS<-ifelse(meta$OS=='DECEASED',1,0)

array<-array[match(meta$name,rownames(array)),]

identical(rownames(array),meta$name)

val<-cbind(array,meta)

merge<-val
cox.data = merge[,colnames(merge) %in% c('OS','OS.time',names(cox$coefficients))]
cox.data.step <- na.omit(cox.data[, c("OS.time", "OS", colnames(cox.data)[!(colnames(cox.data) %in% c('OS','OS.time'))])])

#LASSO版本得分
lasso.prob<-predict(cv,newx = as.matrix(array),s=c(cv$lambda.min,cv$lambda.1se))
lasso.prob<-lasso.prob[match(rownames(cox.data.step),rownames(lasso.prob)),]
cox.data.plot<-cbind(cox.data.step, lasso.prob[,1])
colnames(cox.data.plot)[ncol(cox.data.plot)]<-c('riskScore')
cox.data.plot$OS.time<-as.numeric(cox.data.plot$OS.time)
#
library(survivalROC)
heatmap_train <- cox.data.plot

# pdf(file="./STEP8_model_val/time_ROC.pdf",height = 6,width = 6)
roc=survivalROC(Stime=heatmap_train$OS.time, status=heatmap_train$OS,span = nrow(heatmap_train)^-1,
                marker = heatmap_train$riskScore,predict.time =1*365)
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='#FA8072',
     xlab="False positive rate", ylab="True positive rate",
     #main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     main="ROC curve",
     lwd = 3, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
aucText=c()
rocCol <- c('#FFC1C1','#63B8FF','#FA8072','#ADFF2F','#FFFF00')
aucText=c(aucText,paste0("1 years"," (AUC=",sprintf("%.3f",roc$AUC),")"))
j =0
for (i in c(3,5)){
  roc1=survivalROC(Stime=heatmap_train$OS.time, status=heatmap_train$OS,span = nrow(heatmap_train)^-1,
                   marker = heatmap_train$riskScore,predict.time =i*365)
  j=j+1
  aucText=c(aucText,paste0(i," years"," (AUC=",sprintf("%.3f",roc1$AUC),")"))
  lines(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j+1],lwd = 3)
}
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex = 1)
abline(0,1)
# dev.off()

###########################################
#完全不行GSE的验证试试CGGA呢？

load('data_preparation/CGGA_data.RData')

rpkmTOtpm <- function(rpkm){
  exp(log(rpkm) - log(sum(rpkm)) + log(1e6))
}
cggaexp <- apply(cggaexp, 2, rpkmTOtpm)

genes_expr<-log2(cggaexp+1)

test<-as.matrix(coef (cv,s = "lambda.min"))

array<-as.data.frame(
  t(genes_expr[rownames(genes_expr) %in% rownames(test),]))

#merge array
meta<-cggapheno[,c(1,8,7)]
colnames(meta)<-c('name','OS','OS.time')
meta<-meta[!is.na(meta$OS)]

array<-array[match(meta$name,rownames(array)),]

identical(rownames(array),meta$name)

val<-cbind(array,meta)
merge<-val
cox.data = merge[,colnames(merge) %in% c('OS','OS.time',rownames(test))]
cox.data.step <- na.omit(cox.data[, c("OS.time", "OS", colnames(cox.data)[!(colnames(cox.data) %in% c('OS','OS.time'))])])

#计算得分

# riskScore = predict(cox, type = "risk", newdata = cox.data.step)
# cox.data.plot <- cbind(cox.data.step, riskScore)
# cox.data.plot$OS <- as.numeric(as.character(cox.data.plot$OS))
# cox.data.plot$riskScore <- as.numeric(cox.data.plot$riskScore)
# cox.data.plot <- cox.data.plot[order(cox.data.plot$riskScore), ]
# cox.data.plot$OS.time<-as.numeric(cox.data.plot$OS.time)
# write.csv(cox.data.plot, file = "./STEP8_model_val/val_risk_score.csv",)

#试试LASSO版本得分
lasso.prob<-predict(cv,newx = as.matrix(array),s=c(cv$lambda.min,cv$lambda.1se))
lasso.prob<-lasso.prob[match(rownames(cox.data.step),rownames(lasso.prob)),]
cox.data.plot<-cbind(cox.data.step, lasso.prob[,1])
colnames(cox.data.plot)[ncol(cox.data.plot)]<-c('riskScore')

#
library(survivalROC)
heatmap_train <- cox.data.plot

pdf(file="./04_model_val/CGGA_time_ROC.pdf",height = 6,width = 6)
roc=survivalROC(Stime=heatmap_train$OS.time, status=heatmap_train$OS,span = nrow(heatmap_train)^-1,
                marker = heatmap_train$riskScore,predict.time =1*365)
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='#FA8072',
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

#
#
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
                        title=paste0("CGGA survival"))

pdf("./04_model_val/CGGA_survplot.pdf",width = 8,height = 6)
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

#
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

exp_dat<-plot[,colnames(plot) %in% rownames(test)]
tmp=t(scale(exp_dat))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1

cox_coef<-as.data.frame(test)
colnames(cox_coef)<-c('coef')
cox_coef$gene<-rownames(cox_coef)
cox_coef<-cox_coef[order(cox_coef$coef),]
tmp<-tmp[match(rownames(cox_coef),rownames(tmp)),]

test <- tmp %>% as.data.frame()%>% dplyr::mutate(B=row.names(.)) %>% reshape::melt()

test$B<-factor(test$B,levels = rev(rownames(cox_coef)))
test$variable<-factor(test$variable,levels = colnames(tmp))

p3<-ggplot(test,aes(x=variable,y=B,fill=value))+
  geom_raster()+
  scale_fill_gradient2(low="#003366", high="#990033", mid="white")+
  labs(x=NULL, y = NULL)+
  theme(axis.text.x = element_blank(),axis.ticks=element_blank(),
        legend.position = "none")


#拼图

pdf("./04_model_val/cgga_riskscore.pdf",width = 8,height = 8)
cowplot::plot_grid(p1,p2,p3,ncol=1,align = "hv")
dev.off()
