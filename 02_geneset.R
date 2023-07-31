setwd('C:/Users/Administrator/Desktop/project/YQDZ-0506-6')
rm(list= ls())
library('forestmodel')
library('org.Hs.eg.db')
library('clusterProfiler')
library('patchwork')
# dir.create('02_genescreen')

################################################################################
urigene<-read.csv('./data_preparation/GeneCards-SearchResults_all.csv')
urigene<-urigene[1:267,]
urigene<-urigene[urigene$Category=='Protein Coding',]

change<-read.csv('./01_DEG/TCGA_tumor_VS_normal_changegene.csv')

library(ggvenn)

venn<-list('GBM-changed genes'=c(change$X),'Uridine-Related genes'=c(urigene$Gene.Symbol))

pdf("./02_genescreen/VENN.pdf",width = 6,height = 4)
ggvenn(venn,fill_color = c("#0073C2FF", "#EFC000FF"),
       stroke_size = 1, set_name_size = 6)
dev.off()


genelist<-change$X[change$X%in% urigene$Gene.Symbol]
write.csv(genelist,file = '02_genescreen/genelist.csv')

#################################################################
#富集分析

enrichlist<-bitr(genelist,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb =org.Hs.eg.db )[2]

p<-list()
enrich<-unique(enrichlist)$ENTREZID
kk.negative <- enrichKEGG(gene  = enrich,
                            organism = "hsa",
                            #universe = gene_all,
                            pvalueCutoff = 0.5,
                            qvalueCutoff =0.5)
head(kk.negative)[,1:6]
kk=kk.negative@result
kk<-kk[order(kk$pvalue,decreasing = F),]
  
library('stringr')
title<-paste0('KEGG')
color<-c(met.brewer('Hokusai1')[c(3,5,6,7)],met.brewer('Troy')[5])
  
p[[1]]<-
    ggplot(data = kk[1:10,], aes(x = reorder(Description,-p.adjust), y = -log10(p.adjust))) + #####这里有一个reorder函数，可以对元素进行排序，这里使用y轴值进行了重新排列
    geom_col(aes(fill='#d8443c',colour='#d8443c'),show.legend = FALSE,width=0.8) +
    ggtitle(title) +coord_flip() +
    theme_classic(base_size = 13)+ 
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank(),
          text = element_text(size=12,face='bold'))+
    labs(y="-log10(FDR)",x=NULL)+
    scale_fill_manual(values = c(color[1]))+
    scale_colour_manual(values = c("black"))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
    theme(plot.title = element_text(size=12))
  
BP <- enrichGO(gene          = enrich,
                 OrgDb         = org.Hs.eg.db,
                 ont           = 'BP' ,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.2,
                 qvalueCutoff  = 0.5,
                 readable      = TRUE)
  
CC <- enrichGO(gene          = enrich,
                 OrgDb         = org.Hs.eg.db,
                 ont           = 'CC' ,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.2,
                 qvalueCutoff  = 0.5,
                 readable      = TRUE)
  
MF <- enrichGO(gene          = enrich,
                 OrgDb         = org.Hs.eg.db,
                 ont           = 'MF' ,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.2,
                 qvalueCutoff  = 0.5,
                 readable      = TRUE)
  
BP=BP@result
BP<-BP[order(BP$pvalue,decreasing = F),]
CC=CC@result
CC<-CC[order(CC$pvalue,decreasing = F),]
MF=MF@result
MF<-MF[order(MF$pvalue,decreasing = F),]
  
title<-paste0('GO_BP')
p[[2]]<-
    ggplot(data = BP[1:10,], aes(x = reorder(Description,-p.adjust), y = -log10(p.adjust))) + #####这里有一个reorder函数，可以对元素进行排序，这里使用y轴值进行了重新排列
    geom_col(aes(fill='#d8443c',colour='#d8443c'),show.legend = FALSE,width=0.8) +
    ggtitle(title) +coord_flip() +
    theme_classic(base_size = 13)+ 
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank(),
          text = element_text(size=12,face='bold'))+
    labs(y="-log10(FDR)",x=NULL)+
    scale_fill_manual(values = c(color[2]))+
    scale_colour_manual(values = c("black"))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
    theme(plot.title = element_text(size=12))
  
title<-paste0('GO_CC')
p[[3]]<-
    ggplot(data = CC[1:10,], aes(x = reorder(Description,-p.adjust), y = -log10(p.adjust))) + #####这里有一个reorder函数，可以对元素进行排序，这里使用y轴值进行了重新排列
    geom_col(aes(fill='#d8443c',colour='#d8443c'),show.legend = FALSE,width=0.8) +
    ggtitle(title) +coord_flip() +
    theme_classic(base_size = 13)+ 
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank(),
          text = element_text(size=12,face='bold'))+
    labs(y="-log10(FDR)",x=NULL)+
    scale_fill_manual(values = c(color[3]))+
    scale_colour_manual(values = c("black"))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
    theme(plot.title = element_text(size=12))
  
title<-paste0('GO_MF')
p[[4]]<-
    ggplot(data = MF[1:10,], aes(x = reorder(Description,-p.adjust), y = -log10(p.adjust))) + #####这里有一个reorder函数，可以对元素进行排序，这里使用y轴值进行了重新排列
    geom_col(aes(fill='#d8443c',colour='#d8443c'),show.legend = FALSE,width=0.8) +
    ggtitle(title) +coord_flip() +
    theme_classic(base_size = 13)+ 
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank(),
          text = element_text(size=12,face='bold'))+
    labs(y="-log10(FDR)",x=NULL)+
    scale_fill_manual(values = c(color[4]))+
    scale_colour_manual(values = c("black"))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
    theme(plot.title = element_text(size=12))
  
pdf('01_DEG/enrich.pdf',width = 12,height = 8)
print(wrap_plots(p,nrow=2))   #这张图也调为fdr，和GSEA的图保持一致
dev.off()
  
library(openxlsx)  #把富集的结果保存为xlsx文件，这个包还挺方便的
sheets = list("KEGG" = kk,"GO_BP" = BP,'GO_CC'=CC,'GO_MF'=MF) #KEGG的geneid没有转换过来
write.xlsx(sheets,'01_DEG/enrich.xlsx')

##################################################################
#PPI
#用cytoscape画好了

