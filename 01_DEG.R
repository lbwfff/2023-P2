setwd('C:/Users/Administrator/Desktop/project/YQDZ-0506-6')
rm(list= ls())
# options(download.file.method="libcurl")
# options(url.method="libcurl")
# Sys.setenv(LANGUAGE = "en")
# dir.create('./01_DEG')

library('ggplot2')
library('MetBrewer')
library('ggpubr')
library('data.table')
library('DESeq2')
library('clusterProfiler')
library('tidyr')
library('tibble')
library('corrplot')
library('ggrepel')

##############################################################
#TCGA数据差异分析后和尿苷相关基因取交集

load('./data_preparation/TCGA_data.RData')

tcgaexp<-((2^tcgacount)-1)
tcgaexp<-apply(tcgaexp,2,round)

pheno<-pheno[match(colnames(tcgaexp),pheno$submitter_id.samples),]
group<-ifelse(pheno$sample_type.samples=='Solid Tissue Normal','normal','tumor')

condition <- factor(group)
dds <- DESeqDataSetFromMatrix(tcgaexp, DataFrame(condition), design= ~ condition )
dds

dds$condition<- relevel(dds$condition, ref = "normal") # 指定哪一组作为对照组
dds <- DESeq(dds)
dds
normalized_counts <- counts(dds,normalized=T) 

res= results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)

write.csv(res,file="./01_DEG/TCGA_tumor_VS_normal.csv")

deseq<-read.csv('./01_DEG/TCGA_tumor_VS_normal.csv')
deseq<-deseq[deseq$baseMean>5,]

deseq$group<-factor(ifelse(deseq$padj < 0.05 & abs(deseq$log2FoldChange) >= 1, 
                           ifelse(deseq$log2FoldChange>= 1 ,'Up','Down'),'NoSignifi'),
                    levels=c('Up','Down','NoSignifi'))
table(deseq$group)
changegene<-deseq[deseq$group!='NoSignifi',]

write.csv(changegene,file="./01_DEG/TCGA_tumor_VS_normal_changegene.csv")

mostup<-deseq[deseq$log2FoldChange>0,]
mostup<-mostup[order(mostup$padj),][1:150,]
write.csv(mostup,file = 'mostup_150_gene.csv')

mostdown<-deseq[deseq$log2FoldChange<0,]
mostdown<-mostdown[order(mostdown$padj),][1:150,]
write.csv(mostdown,file = 'mostdown_150_gene.csv')

################################
#可视化
table<-deseq
table<-table[!is.na(table$padj),]

library(ggplot2)
library('ggrepel')

pdf("./01_DEG/TCGA_volcano.pdf",width = 9,height = 6)
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

df <- data.frame(group = c(group))
df$group<-factor(df$group)
df$what<-c('3')

mat_scaled<-mat_scaled[,c(which(df$group=='normal'),which(df$group!='normal'))]
df<-df[c(which(df$group=='normal'),which(df$group!='normal')),]

split = c(ifelse(df$group=='normal',1,2))
ha <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 3:2),labels = c("normal", "tumor"), 
                                         labels_gp = gpar(col = "white", fontsize = 14)))

ha <- HeatmapAnnotation(type = c(df$group),
                       col = list(type = c("normal" = 3, "tumor" = 2), c("white", "red")))

pdf('./01_DEG/TCGA_heatmap.pdf',width = 18,height = 10)
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


# col_fun1 = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))
# circos.heatmap(mat_scaled, col = col_fun1,dend.side = "inside",rownames.side = "outside")
# lgd = Legend(title = "z-score", col_fun = col_fun1)
# grid.draw(lgd)
# circos.clear()


