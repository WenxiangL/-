rm(list =ls())
getwd()
setwd("E:/生信/组学数据/")
LR <- read.table("celltalk/mouse_lr_pair (3) .txt",sep ="\t",header = T)

#选取GeneID
LR_geneID <- LR[ , c(1,2,3)]

#-------------------------11/21/2022_心肺RNAseq_LR------------------
total_rawdata <- read.csv("natueBulk_RNASeq/分析原始数据/GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv")
#选取肺数据，年轻组在前，老年组在后
library(dplyr)
#肺
rawdata <- total_rawdata %>% select(gene,A20_384Bulk_Plate1_S20.gencode.vM19,K20_384Bulk_Plate1_S260.gencode.vM19,B3_384Bulk_Plate3_S27.gencode.vM19,F4_384Bulk_Plate1_S124.gencode.vM19,
                                    M6_384Bulk_Plate1_S294.gencode.vM19,K9_384Bulk_Plate1_S249.gencode.vM19,L8_384Bulk_Plate2_S272.gencode.vM19,E2_384Bulk_Plate2_S98.gencode.vM19)

#心脏
rawdata <- total_rawdata %>% select(gene,O20_384Bulk_Plate2_S356.gencode.vM19,P11_384Bulk_Plate2_S371.gencode.vM19,K9_384Bulk_Plate2_S249.gencode.vM19,I13_384Bulk_Plate2_S205.gencode.vM19,
                                    F11_384Bulk_Plate1_S131.gencode.vM19,G3_384Bulk_Plate1_S147.gencode.vM19,G8_384Bulk_Plate2_S152.gencode.vM19,C19_384Bulk_Plate2_S67.gencode.vM19)

#更改列名
#K6-O1,E6-O2,C10-O3,G12-O4;A20-Y1,K20-Y2,B3-Y3,F4-Y4
colnames(rawdata) <- c("gene_id","Y1","Y2","Y3","Y4","O1","O2","O3","O4")

#将数据框第一列设置为行名
rownames(rawdata)<-rawdata[,1]
rawdata<-rawdata[,-1] 


#------------------差异表达基因分析--------------
library(DESeq2)
#读取counts数据和meta数据表
mycounts_1 <- rawdata
mymeta <- read.csv("natueBulk_RNASeq/分析原始数据/meta.csv")


dds <- DESeqDataSetFromMatrix(countData = round(mycounts_1),
                              colData = mymeta,
                              design =~dex)
dds <- DESeq(dds)
res <- results(dds)
res_1 <- data.frame(res)
library(dplyr)
res_1 %>%
  mutate(group = case_when(
    log2FoldChange >= 1 & padj <= 0.05 ~ "UP",
    log2FoldChange <= -1 & padj <= 0.05 ~ "DOWN",
    TRUE ~"NOT_CHANGE"
  )) -> res_2
#查看有多少差异表达基因
table(res_2$group)
#输出结果
write.csv(res_2,file = "natueBulk_RNASeq/差异表达分析结果/heart_differented_2.csv",quote = F)

df_gene.result <-read.csv("E:/生信/组学数据/natueBulk_RNASeq/差异表达分析结果",
                          header = T, stringsAsFactors = F )



