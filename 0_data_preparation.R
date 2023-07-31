# dir.create('/data/project-libw/R/YQDZ-0506-6/data_preparation')
setwd('C:/Users/Administrator/Desktop/project/YQDZ-0506-6/data_preparation')
rm(list= ls())
##################数据下载和前处理####################

##########TCGA##############

library('data.table')

tcgaexp<-fread('TCGA-GBM.htseq_fpkm.tsv.gz')
tcgacount<-fread('TCGA-GBM.htseq_counts.tsv.gz')

pheno<-fread('TCGA-GBM.GDC_phenotype.tsv.gz')
phenosup<-read.csv('tcga_cli_sup.csv')

suv<-fread('TCGA-GBM.survival.tsv')

mut<-fread('TCGA-GBM.muse_snv.tsv.gz')

#ID转化
idtrans<-fread('gencode.v22.annotation.gene.probeMap')
idtrans<-idtrans[match(tcgaexp$Ensembl_ID,idtrans$id),]

idtrans$exp<-apply(tcgaexp[,-1],1,mean) #根据表达量过滤重复基因
idtrans<-idtrans[order(idtrans$exp,decreasing = T),]
idtrans<-idtrans[!duplicated(idtrans$gene),]
idtrans<-idtrans[!is.na(idtrans$id) &idtrans$exp>0,]

tcgaexp<-as.data.frame(tcgaexp)
tcgaexp<-tcgaexp[tcgaexp$Ensembl_ID %in% idtrans$id,]

idtrans<-idtrans[match(tcgaexp$Ensembl_ID,idtrans$id),]
rownames(tcgaexp)<-idtrans$gene
tcgaexp<-tcgaexp[,-1]

idtrans<-fread('gencode.v22.annotation.gene.probeMap')
idtrans<-idtrans[match(tcgacount$Ensembl_ID,idtrans$id),]

idtrans$exp<-apply(tcgacount[,-1],1,mean) #根据表达量过滤重复基因
idtrans<-idtrans[order(idtrans$exp,decreasing = T),]
idtrans<-idtrans[!duplicated(idtrans$gene),]
idtrans<-idtrans[!is.na(idtrans$id) &idtrans$exp>0,]

tcgacount<-as.data.frame(tcgacount)
tcgacount<-tcgacount[tcgacount$Ensembl_ID %in% idtrans$id,]

idtrans<-idtrans[match(tcgacount$Ensembl_ID,idtrans$id),]
rownames(tcgacount)<-idtrans$gene
tcgacount<-tcgacount[,-1]

#表型
pheno<-pheno[pheno$submitter_id.samples %in% colnames(tcgaexp),]
phenosup<-phenosup[phenosup$sampleID %in% substr(pheno$submitter_id.samples,1,15),]
phenosup<-phenosup[match(substr(pheno$submitter_id.samples,1,15),phenosup$sampleID),]
pheno<-cbind(pheno,phenosup)

#生存
suv<-suv[suv$sample %in% colnames(tcgaexp),]

save(tcgaexp,tcgacount,pheno,suv,mut,file = 'TCGA_data.RData')

##########GSE13041###########

library(GEOquery)
library(AnnoProbe)
Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)

GSE13041_1<-getGEO(filename ='./GSE13041-GPL96_series_matrix.txt.gz',getGPL = FALSE)

exp_1<-exprs(GSE13041_1)

(gpl=GSE13041_1@annotation)
checkGPL(gpl)
printGPLInfo(gpl)
probe2gene=idmap(gpl)
head(probe2gene)
exp_1 <- filterEM(exp_1,probe2gene)
head(exp_1)
table(duplicated(rownames(exp_1)))

meta1<-pData(GSE13041_1)

save(exp_1,meta1,file = 'GSE13041_ready.RData')

##########CGGA############

cggaexp<-fread('CGGA.mRNAseq_693.RSEM-genes.20200506.txt')
cggacount<-fread('CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt')
cggapheno<-fread('CGGA.mRNAseq_693_clinical.20200506.txt')

cggaexp<-as.data.frame(cggaexp)
row.names(cggaexp)<-cggaexp$Gene_Name
cggaexp<-cggaexp[,-1]

cggacount<-as.data.frame(cggacount)
row.names(cggacount)<-cggacount$gene_name
cggacount<-cggacount[,-1]

save(cggaexp,cggacount,cggapheno,file = 'CGGA_data.RData')

