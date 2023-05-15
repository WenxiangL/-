rm(list =ls())
getwd()
setwd("E:/生信/组学数据/natueBulk_RNASeq/配体受体分析结果/")
different_gene <- read.csv("heart_recptor.csv",header = T,stringsAsFactors = F)
table(different_gene$group)

head(different_gene)

#2023/5/9
different_gene$log2FoldChange <- different_gene$log2FoldChange*(-1)

different_gene %>%
  mutate(group = case_when(
    log2FoldChange >= 1 & padj <= 0.05 ~ "UP",
    log2FoldChange <= -1 & padj <= 0.05 ~ "DOWN",
    TRUE ~"NOT_CHANGE"
  )) -> different_gene
table(different_gene$group)

#2023/5/9

# 设置p_value和logFC的阈值
cut_off_pvalue = 0.05  #统计显著性
cut_off_logFC = 1           #差异倍数值

# 根据阈值分别为上调基因设置‘up’，下调基因设置‘Down’，无差异设置‘None’，保存到change列
# 这里的change列用来设置火山图点的颜色
different_gene$group = ifelse(different_gene$pvalue < cut_off_pvalue & 
                   abs(different_gene$log2FoldChange) >= cut_off_logFC, 
                 ifelse(different_gene$log2FoldChange> cut_off_logFC ,'UP','DOWN'),'NOT_CHANGE')

"dodgerblue","gray","firebrick"

library(ggplot2)
ggplot(different_gene,aes(x=log2FoldChange,y = -log10(padj)))+#横轴与纵轴
  geom_point(aes(color=group))+
  scale_color_manual(values = c("gray"))+#选择颜色
  xlim(-1,1)+ylim(0,2)+theme_bw()
df1 <- different_gene$padj
