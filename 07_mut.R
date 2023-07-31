setwd('C:/Users/Administrator/Desktop/project/YQDZ-0506-6')
rm(list= ls())

# dir.create('./07_mut')
library('maftools')
library('TCGAbiolinks')
##########################################################################
load('./data_preparation/TCGA_data.RData')
risk<-read.csv("./03_Prognosis/risk_score.csv")

#下载maf格式的突变数据
# query <- GDCquery(
#   project = "TCGA-GBM", 
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
# GDCdownload(query)
# GDCprepare(query, save = T,save.filename = "TCGA_GBM_SNP.Rdata")

load('TCGA_GBM_SNP.Rdata')

tcgamut<-data[,-1]
tcgamut$sample<-substr(tcgamut$Tumor_Sample_Barcode,1,16)

risk$group<-ifelse(risk$riskScore>median(risk$riskScore),'Hight','Low')

muthi<-tcgamut[tcgamut$sample %in% risk$X[risk$group=='Hight'],]
mutlo<-tcgamut[tcgamut$sample %in% risk$X[risk$group=='Low'],]

laml = read.maf(maf = muthi,
                isTCGA = TRUE)

pdf('07_mut/riskscore_hig_summary.pdf',width = 10,height = 8)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del','Missense_Mutation',
  'Nonsense_Mutation','Multi_Hit',
  'Frame_Shift_Ins','In_Frame_Ins',
  'Splice_Site','In_Frame_Del')

pdf('07_mut/riskscore_hig_oncoplot.pdf',width = 10,height = 8)
oncoplot(maf = laml, top = 20,colors = vc_cols)
dev.off()

pdf('07_mut/riskscore_hig_interact.pdf',width = 8,height = 8)
Interact <- somaticInteractions(maf = laml, top = 25, showSum=F,
                                pvalue = c(0.05, 0.1))
dev.off()

laml2 = read.maf(maf = mutlo,
                isTCGA = TRUE)

pdf('07_mut/riskscore_low_summary.pdf',width = 10,height = 8)
plotmafSummary(maf = laml2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

pdf('07_mut/riskscore_low_oncoplot.pdf',width = 10,height = 8)
oncoplot(maf = laml2, top = 20,colors = vc_cols)
dev.off()

pdf('07_mut/riskscore_low_interact.pdf',width = 8,height = 8)
Interact <- somaticInteractions(maf = laml2, top = 25, showSum=F,
                                pvalue = c(0.05, 0.1))
dev.off()

#两组比较
pt.vs.rt <- mafCompare(m1 = laml, m2 = laml2, 
                       m1Name = 'riskScore_Hig', pseudoCount=T,
                       m2Name = 'riskScore_Low', minMut = 5)
print(pt.vs.rt)

pdf('07_mut/mutcompare.pdf',width = 8,height = 6)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, 
           color = c('royalblue', 'maroon'), geneFontSize = 0.8)
dev.off()


