setwd('C:/Users/Administrator/Desktop/project/YQDZ-0506-6')
rm(list= ls())

# dir.create('./06_GSEA')

library('ggplot2')
library('MetBrewer')
library('ggpubr')
library('DESeq2')
library('clusterProfiler')
library('tidyr')
library('tibble')
library('corrplot')
library('ggrepel')
library('stringr')
library('patchwork')
library('RColorBrewer')
library('GseaVis')

##############################################################
#高低风险组的GSEA

load('data_preparation/TCGA_data.RData')
risk<-read.csv("./03_Prognosis/risk_score.csv")
risk$group<-ifelse(risk$riskScore>median(risk$riskScore),'risk_high','risk_low')

tcgaexp<-((2^tcgacount)-1)
tcgaexp<-apply(tcgaexp,2,round)

tcgaexp<-tcgaexp[,colnames(tcgaexp) %in% risk$X]
risk<-risk[match(colnames(tcgaexp),risk$X),]

condition <- factor(risk$group)
dds <- DESeqDataSetFromMatrix(tcgaexp, DataFrame(condition), design= ~ condition )
dds

dds$condition<- relevel(dds$condition, ref = "risk_low") # 指定哪一组作为对照组
dds <- DESeq(dds)
dds

res= results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)

write.csv(res,file="./06_GSEA/TCGA_scorehight_VS_scorelow.csv")

deseq<-read.csv('./06_GSEA/TCGA_scorehight_VS_scorelow.csv')
deseq<-deseq[deseq$baseMean>5,]

deseq$group<-factor(ifelse(deseq$padj < 0.05 & abs(deseq$log2FoldChange) >= 1, 
                           ifelse(deseq$log2FoldChange>= 1 ,'Up','Down'),'NoSignifi'),
                    levels=c('Up','Down','NoSignifi'))
################################
#可视化
table<-deseq
table<-table[!is.na(table$padj),]

library('ggplot2')
library('ggrepel')

pdf("./06_GSEA/DEG_volcano.pdf",width = 9,height = 6)
ggplot(table,aes(x=log2FoldChange,y=-log10(padj),color=group))+
  geom_point(size=2)+
  # geom_text_repel(data=table[table$group!='NoSignifi',],min.segment.length = Inf, box.padding = 0.5,
  #                 aes(label=X),size=3,segment.color='black',show.legend=F)+
  theme_classic(base_size = 18)+ 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),panel.border = element_blank())+
  ylab('-log10 (P adj)')+
  xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)+
  scale_color_manual(values=c("#dd5129","#808fe1","#b9b9b8"),
                     breaks=c("Up", "Down", "NoSignifi"),
                     labels=c(paste0('Up',' (',nrow(table[table$group=='Up',]),')'),
                              paste0('Down',' (',nrow(table[table$group=='Down',]),')'),
                              paste0('NoSignifi',' (',nrow(table[table$group=='NoSignifi',]),')')))
dev.off()

#热图
normalized_counts <- counts(dds,normalized=T) 
n<-50
pheatmap<-normalized_counts[rownames(normalized_counts) %in% 
                              deseq$X[1:n],]
genelist<-deseq[1:n,]
genelist<-genelist[order(genelist$log2FoldChange,decreasing = T),]
pheatmap<-pheatmap[match(genelist$X,rownames(pheatmap)),]

library(ComplexHeatmap)
library(circlize)

mat_scaled = t(scale(t(pheatmap)))
col_fun1 = colorRamp2(c(-1,0,1), c(met.brewer("Hiroshige")[9],'white',met.brewer("Signac")[4]))

df <- data.frame(group = c(risk$group))
df$group<-factor(df$group)
df$what<-c('3')

mat_scaled<-mat_scaled[,c(which(df$group=='risk_low'),which(df$group!='risk_low'))]
df<-df[c(which(df$group=='risk_low'),which(df$group!='risk_low')),]

split = c(ifelse(df$group=='risk_low',1,2))
ha <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 3:2),labels = c("risk_low", "risk_high"), 
                                         labels_gp = gpar(col = "white", fontsize = 14)))

# ha <- HeatmapAnnotation(type = c(df$group),
#                         col = list(type = c("normal" = 3, "tumor" = 2), c("white", "red")))

pdf('./06_GSEA/RISK_heatmap.pdf',width = 18,height = 10)
draw(
  Heatmap(mat_scaled,
          col = col_fun1,name = "Score",
          column_split = split,
          column_title = NULL,
          column_gap = unit(c(0.1, 0.1), "mm"),
          cluster_columns = F,
          cluster_rows =F,
          show_row_dend = T,
          show_column_names =F,
          show_heatmap_legend = T,
          row_names_side = "right",
          width = ncol(mat_scaled)*unit(2, "mm"), 
          height = nrow(mat_scaled)*unit(4, "mm"),
          top_annotation = ha))
dev.off()


#为了做imap的预测，把最显著的150个基因提出来
mostup<-deseq[deseq$log2FoldChange>0,]
mostup<-mostup[order(mostup$padj),][1:150,]
write.csv(mostup,file = './06_GSEA/mostup_150_gene.csv')

mostdown<-deseq[deseq$log2FoldChange<0,]
mostdown<-mostdown[order(mostdown$padj),][1:150,]
write.csv(mostdown,file = './06_GSEA/mostdown_150_gene.csv')

#
gsea<-deseq

gsea<-gsea[order(gsea$log2FoldChange,decreasing = T),]

test<-as.numeric(gsea$log2FoldChange)
names(test) = as.character(gsea$X)


mygsea<- function(file,geneset,path){ 
  #选择gmt文件
  gmtfile =file
  
  #GSEA分析
  library(GSEABase) 
  geneset <- clusterProfiler::read.gmt( gmtfile )  
  length(unique(geneset$term))
  egmt <- GSEA(test, TERM2GENE=geneset, 
               minGSSize = 5,eps=0,
               pvalueCutoff = 0.1,
               verbose=FALSE)
  
  head(egmt)
  egmt@result 
  gsea_results_df <- egmt@result 
  rownames(gsea_results_df)
  write.csv(gsea_results_df,file = path)
  return(egmt)
}

p<-list() #这一段写得巨丑无比，找个时间再优化吧

kegg<-mygsea('./data_preparation/c2.cp.kegg.v2023.1.Hs.symbols.gmt',
       test,'./06_GSEA/kegg_gsea.csv')

for (j in 1:length(kegg@result$Description)){
  de<-kegg@result$Description[j]
  kegg@result$Description[j]<-
    paste0(unlist(strsplit(de,'_'))[-1],collapse = ' ')
}

p[[1]]<-
dotplotGsea(data = kegg,topn = 10,order.by = 'NES',add.seg = T)$plot

gobp<-mygsea('./data_preparation/c5.go.bp.v2023.1.Hs.symbols.gmt',
       test,'./06_GSEA/gobp_gsea.csv')

for (j in 1:length(gobp@result$Description)){
  de<-gobp@result$Description[j]
  gobp@result$Description[j]<-
    paste0(unlist(strsplit(de,'_'))[-1],collapse = ' ')
}

p[[2]]<-
  dotplotGsea(data = gobp,topn = 10,order.by = 'NES',add.seg = T)$plot

gocc<-mygsea('./data_preparation/c5.go.cc.v2023.1.Hs.symbols.gmt',
       test,'./06_GSEA/gocc_gsea.csv')

for (j in 1:length(gocc@result$Description)){
  de<-gocc@result$Description[j]
  gocc@result$Description[j]<-
    paste0(unlist(strsplit(de,'_'))[-1],collapse = ' ')
}

p[[3]]<-
  dotplotGsea(data = gocc,topn = 10,order.by = 'NES',add.seg = T)$plot

gomf<-mygsea('./data_preparation/c5.go.mf.v2023.1.Hs.symbols.gmt',
       test,'./06_GSEA/gomf_gsea.csv')

for (j in 1:length(gomf@result$Description)){
  de<-gomf@result$Description[j]
  gomf@result$Description[j]<-
    paste0(unlist(strsplit(de,'_'))[-1],collapse = ' ')
}

p[[4]]<-
  dotplotGsea(data = gomf,topn = 10,order.by = 'NES',add.seg = T)$plot

pdf('06_GSEA/GSEA_2.pdf',width = 22,height = 14)
print(wrap_plots(p,nrow=2)) 
dev.off()

#画个图

p<-list()
list<-c('kegg','gobp','gocc','gomf')
name<-c('KEGG','GO_BP','GO_CC','GO_MF')
color<-c(met.brewer('Hokusai1')[c(3,5,6,7)],met.brewer('Troy')[5])

for(i in 1:4){
  file=paste0('./06_GSEA/',list[i],'_gsea.csv')
  inf<-read.csv(file)
  inf<-inf[order(inf$p.adjust),]
  
  for (j in 1:10){
    de<-inf$Description[j]
    inf$Description[j]<-
      paste0(unlist(strsplit(de,'_'))[-1],collapse = ' ')
  }
  

  p[[i]]<-
    ggplot(data = inf[1:10,], aes(x = reorder(Description,-p.adjust), y = -log10(pvalue))) + #####这里有一个reorder函数，可以对元素进行排序，这里使用y轴值进行了重新排列
    geom_col(aes(fill='#d8443c',colour='#d8443c'),show.legend = FALSE,width=0.8) +
    ggtitle(name[i]) +coord_flip() +
    theme_classic(base_size = 13)+ 
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank(),
          text = element_text(size=12,face='bold'))+
    labs(y="-log10(FDR)",x=NULL)+
    scale_fill_manual(values = c(color[i]))+
    scale_colour_manual(values = c("black"))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
    theme(plot.title = element_text(size=12))
  
}

pdf('06_GSEA/GSEA.pdf',width = 12,height = 8)
print(wrap_plots(p,nrow=2)) 
dev.off()

#不对不对GSEA不能这样表现，这样就没有正负了。



########################################################
#免疫浸润
write.table(tcgaexp,file = 'TCGA_exp.txt',sep = '\t',quote = F)

source('Cibersort_source.R')

result <- CIBERSORT('LM22.txt','TCGA_exp.txt', perm = 100, QN = F)
write.csv(result,file = '06_GSEA/cibersort.csv')
result<-read.csv('06_GSEA/cibersort.csv')

group<-data.frame(name=c(risk$X),
                  exp=c(risk$riskScore),
                  group=c(NA))
group$group<-ifelse(group$exp>median(group$exp),
                    'risk Hight','risk Low')

#合并分组信息
result1 = as.data.frame(result)
result1$name<-rownames(result1)
result1 <- merge(group,result1,by="name")

row.names(result1) <- result1[,1]
result1 <- result1[,-1]
write.csv(result1,"./06_GSEA/CIBERSORT_result_merge.csv")

#根据pf4表达水平组排序
re1 <- result1[order(result1[,1],decreasing = T),]
colnames(re1)[2] <- "Type"

re1 = subset(re1,re1$`P-value`<0.05)
re1 = re1[,-c(1,25:27)]
re2<-re1[,-1]

mypalette <- colorRampPalette(brewer.pal(8,"Set3"))
#提取数据，多行变成多列，要多学习‘tidyr’里面的三个函数
dat_cell <- re2 %>% as.data.frame() %>%rownames_to_column("Sample") %>%gather(key = Cell_type,value = Proportion,-Sample)
#提取数据
dat_group = gather(re1,Cell_type,Proportion,-Type )
#合并分组
dat = cbind(dat_cell,dat_group$Type)
colnames(dat)[4] <- "Type"

##2.2柱状图##############
p1 <- ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) +
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(28))
p1
#画分组bar
Var1 = re1$Type #high risk & low risk的数量
Var2=c(1:length(row.names(re1))) #high risk & low risk长度，从1开始
annotation_col<- data.frame(Var1 = re1$Type,
                            Var2=c(1:length(row.names(re1))),
                            value=1)
p2 <- ggplot(dat, aes(Var2, Var1),size=0.5)+
  geom_raster(aes(Var2,value,fill = Var1), data= annotation_col)+ 
  labs(x = "", y = "")+
  theme_void()+
  theme(legend.position = "top", legend.title = element_blank())+
  scale_x_discrete(name="")+
  scale_fill_manual(labels=c("High risk group(n=15)","Low risk group(n=25)"),
                    values =c("#DC0000FF","#00A087FF"))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))
p2
#拼图
pdf('./06_GSEA/CIBERSORT_radio.pdf',height = 6,width = 10)
cowplot::plot_grid(p2,p1,
                   rel_heights = c(0.08,1),
                   label_x=0,
                   label_y=1,
                   align = 'v',ncol = 1,greedy = F)
dev.off()

##2.3 免疫细胞热图#############
data <- as.data.frame(t(re1[,-1]))

data_1 <- as.data.frame(lapply(data,as.numeric))
row.names(data_1) <- row.names(data)
data=data_1

library(pheatmap)
pdf(file="./06_GSEA/immunecell_pheatmap.pdf",width=10,height=8,onefile=FALSE)
annotation_col=data.frame(Group=rep(c("Hight risk group(n=15)","Low risk group(n=25)"),c(15,25)))
annColors <- list(Group = c("Hight risk group(n=15)" = "#DC0000FF","Low risk group(n=25)" ="#00A087FF"))
rownames(annotation_col)=colnames(data)
print(pheatmap(data,
         annotation_col = annotation_col,
         annotation_colors = annColors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cluster_rows=TRUE,
         show_rownames=TRUE,
         fontsize_row=10,
         show_colnames=FALSE,
         scale="row",
         cluster_cols=FALSE,
         main=""))
dev.off()

##2.4组免疫细胞的差异###################
library(ggpubr)
library(ggplot2)
box=dat

#图片美化
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.line = element_line(colour = "black"),
          #axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.title.x=element_text(face = "bold",size = 14),#X轴标签加粗及字体大小     
          axis.text.y = element_text(face = "bold",size = 12),#y坐标轴刻度标签加粗
          axis.text.x = element_text(face = "bold",size = 10, vjust = 1, hjust = 1, angle = 45),#x坐标轴刻度标签加粗 
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),
          legend.position=c(0.9, 0.8),#图例在绘图区域的位置
          legend.direction = "horizontal",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=8)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
}

e1 <- ggplot(box,aes(x=reorder(Cell_type,-Proportion),y=Proportion),palette = "jco", add = "jitter")+
  geom_boxplot(aes(fill=Type),position=position_dodge(0.5),width=0.6)+
  labs(x = "Cell Type", y = "Estimated Proportion")+
  scale_fill_manual(values = c("#DC0000FF","#00A087FF")) +
  theme_zg()
e1 = e1 + stat_compare_means(aes(group = Type),label = "p.signif")
e1
#保存图片
ggsave('./06_GSEA/box_Immune.pdf', plot = e1,width=15,height = 6)

#############################################################
#ESTIMATE

library(estimate)

estimate <- function(dat,pro){
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")   ## 注意这个platform参数
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}

scores<-estimate(tcgaexp,'GBM')
scores<-as.data.frame(scores)

cache<-risk
cache$X<-gsub("-", ".", cache$X)
cache<-cache[match(rownames(scores),cache$X),]
scores$riskscore<-cache$riskScore
scores$OS<-cache$OS
scores$OS.time<-cache$OS.time

write.csv(scores,file = '06_GSEA/estimate_score.csv')
scores<-read.csv('06_GSEA/estimate_score.csv')
#看三个评分和生存的关系？写的什么玩意？

matt<-scores
p<-list()
#
med.exp<-median(matt$StromalScore)
more.med.exp.index<-which(matt$StromalScore>=med.exp)
less.med.exp.index<-which(matt$StromalScore< med.exp)
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
                        title=paste0("StromalScore"))

p[[1]]<-
sdata.plot3$plot+
  theme_classic(base_size = 22)+ 
  theme(legend.position = "right")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())+
  theme(aspect.ratio=1)

#
med.exp<-median(matt$ImmuneScore)
more.med.exp.index<-which(matt$ImmuneScore>=med.exp)
less.med.exp.index<-which(matt$ImmuneScore< med.exp)
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
                        title=paste0("ImmuneScore"))

p[[2]]<-
  sdata.plot3$plot+
  theme_classic(base_size = 22)+ 
  theme(legend.position = "right")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())+
  theme(aspect.ratio=1)

#
med.exp<-median(matt$ESTIMATEScore)
more.med.exp.index<-which(matt$ESTIMATEScore>=med.exp)
less.med.exp.index<-which(matt$ESTIMATEScore< med.exp)
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
                        title=paste0("ESTIMATEScore"))

p[[3]]<-
  sdata.plot3$plot+
  theme_classic(base_size = 22)+ 
  theme(legend.position = "right")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())+
  theme(aspect.ratio=1)

pdf('06_GSEA/estimate_suv.pdf',width = 14,height = 10)
print(wrap_plots(p,nrow=2)) 
dev.off()

#看一下三个分数和riskscore相关的箱线
#这里绘图的分组非常奇怪，需要调整游戏

p<-list()

scores$group<-ifelse(scores$riskscore>median(scores$riskscore),
                     'risk_High','risk_Low')
scores$group<-factor(scores$group,levels = c('risk_Low','risk_High'))
# scores$group<-factor(scores$group,levels = c('Younger','Old'))

p[[1]]<-
ggplot(scores, aes(x=group, y=StromalScore))+
  geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
  labs(y="StromalScore")+ 
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


p[[2]]<-
  ggplot(scores, aes(x=group, y=ImmuneScore))+
  geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
  labs(y="ImmuneScore")+ 
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

p[[3]]<-
  ggplot(scores, aes(x=group, y=ESTIMATEScore))+
  geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
  labs(y="ESTIMATEScore")+ 
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

pdf('06_GSEA/estimate_corto_score.pdf',width = 14,height = 6)
print(wrap_plots(p,nrow=1)) 
dev.off()

#######################################################################
#几个免疫检查点的相关性
load('./data_preparation/TCGA_data.RData')

genes_expr<-((2^tcgaexp)-1)
rpkmTOtpm <- function(rpkm){
  exp(log(rpkm) - log(sum(rpkm)) + log(1e6))
}
genes_expr <- apply(genes_expr, 2, rpkmTOtpm)

genes_expr<-log2(genes_expr+1)

checklist<-c('CD80','CD86','CD28',
'ICOS','BTLA','CD274','PDCD1','HAVCR2','CTLA4','LAG3')

checkexp<-genes_expr[rownames(genes_expr) %in% checklist,]
checkexp<-as.data.frame(t(checkexp))

rownames(checkexp)<-gsub("-", ".", rownames(checkexp))
checkexp<-checkexp[rownames(checkexp) %in% rownames(scores),]
cache<-scores[match(rownames(checkexp),rownames(scores)),]
cache$sample<-rownames(cache)

dat_cell <- checkexp %>% as.data.frame() %>%
  rownames_to_column("Sample") %>%
  gather(key = checkpoint,value = exp,-Sample)

cache<-cache[match(dat_cell$Sample,cache$sample),]
dat_cell$group<-cache$group

box=dat_cell

#图片美化
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.line = element_line(colour = "black"),
          #axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.title.x=element_text(face = "bold",size = 14),#X轴标签加粗及字体大小     
          axis.text.y = element_text(face = "bold",size = 12),#y坐标轴刻度标签加粗
          axis.text.x = element_text(face = "bold",size = 10, vjust = 1, hjust = 1, angle = 45),#x坐标轴刻度标签加粗 
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),
          legend.position=c(0.9, 0.8),#图例在绘图区域的位置
          legend.direction = "horizontal",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=8)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
}

e1 <- ggplot(box,aes(x=reorder(checkpoint,-exp),y=exp),palette = "jco", add = "jitter")+
  geom_boxplot(aes(fill=group),position=position_dodge(0.5),width=0.6)+
  labs(x = "checkpoint", y = "expression")+
  scale_fill_manual(values = c("#DC0000FF","#00A087FF")) +
  theme_zg()
e1 = e1 + stat_compare_means(aes(group = group),label = "p.signif")
e1
#保存图片
ggsave('./06_GSEA/box_checkpoint.pdf', plot = e1,width=15,height = 6)

#这个差异有点好得过分了，看看相关性呢？

checkexp<-genes_expr[rownames(genes_expr) %in% checklist,]
checkexp<-as.data.frame(t(checkexp))

rownames(checkexp)<-gsub("-", ".", rownames(checkexp))
checkexp<-checkexp[rownames(checkexp) %in% rownames(scores),]
cache<-scores[match(rownames(checkexp),rownames(scores)),]

checkexp$score<-cache$riskscore

corm<-cor(checkexp) #相关性到没有显得很好

