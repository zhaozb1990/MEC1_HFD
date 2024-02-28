
setwd('H:\\lianlab_others\\ZZB\\20230826_mic1') #设置工作目录，以读取数据，此步不需要改 

source('H:\\lianlab_others\\MJL/20220319/myvolcano.r')
source('H:\\lianlab_others\\MJL/20220319/myGO_KEGG.R')

# directory
merge_count = read.table('./merge_count.txt',header = T,sep = '\t',row.names = 1)

colnames(merge_count)[grep('L',colnames(merge_count))]
# 肝脏：CD：（L2 3 6 8??? + HFD：（L9 11 12 13??? + MIC1：（L15 16 17 20???
merged_liver <- merge_count[,grep('L',colnames(merge_count))]
colnames(merged_liver) <- gsub('.bedGraph','',colnames(merged_liver))
merged_liver <- merged_liver[,c('L2','L3','L6','L8','L9','L11','L12','L13','L15','L16','L17','L20')]
head(merged_liver)

# write.table(merged_liver,file = 'merged_liver.txt')
merged_liver <- read.table('./merged_liver.txt',header = T,row.names = 1)
head(merged_liver)

# 1-7 CD 8-14 HFD 15-20 TIG
merge_count <- merged_liver


merge_count <- merge_count[rowSums(merge_count)>0,]
head(merge_count)


colData = data.frame(row.names = colnames(merge_count),
                     group= c(rep('CD',4),rep('HFD',4),rep('MEC1',4)))

colData$group = as.factor(colData$group)
colData

library(tidyverse)
library(DESeq2)
library(ggrepel)
library(homologene)

merge_count <- merge_count[rowSums(merge_count>0) > 4,]

write.table(merge_count,file = 'select_merge_rf.txt',sep = '\t')

dim(merge_count)

#构建dds对象
dds <- DESeqDataSetFromMatrix(merge_count, colData, design= ~ group)
dds <- DESeq(dds)


res = results(dds, contrast=c("group", "MEC1", "HFD"))
res = res[order(res$pvalue),]
head(res)

write.csv(res,file="All_results_MEC1_vs_HFD.csv") # use for gsea anlysis pre-rnk

  
#获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小???-1的差异表达基???
diff_gene_deseq2 <- subset(res, padj < 0.05 & abs(log2FoldChange) > 0.5)
dim(diff_gene_deseq2)

write.csv(diff_gene_deseq2,file= "DEG_MEC1_vs_HFD_logFC_0.5.csv")


count_data = counts(dds, normalized = TRUE)
library(pheatmap)

count_data = log2(count_data+1)

genelist <- read.table('./HALLMARK_ADIPOGENESIS.tsv',header = T,sep = '\t')
genelist$SYMBOL
library(stringr)


diff_gene_deseq2 <- diff_gene_deseq2[which(rownames(diff_gene_deseq2)%in%genelist$SYMBOL),]
count_data = counts(dds, normalized = TRUE)
count_data <- log2(count_data+1)
count_data1 <- count_data[which(rownames(count_data)%in%rownames(diff_gene_deseq2)),]
a <- as.data.frame(res[which(rownames(res)%in%rownames(count_data1)),1:2])
head(a)
a
count_data1 <- as.data.frame(count_data1)
count_data1$gene <- rownames(count_data1)
head(count_data1)
a$gene <- rownames(a)
count_data1 <- merge(count_data1,a,by='gene')

count_data1 <- count_data1[order(count_data1$log2FoldChange,decreasing = T),]
head(count_data1)
rownames(count_data1) <- count_data1$gene
subset(count_data1,log2FoldChange>0.5)
count_data2 <- count_data1[,2:13]


count_data2 = t(scale(t(count_data2)))

library(stringr)

colData = data.frame(row.names = colnames(count_data2),
                     group= c(rep('CD',4),rep('HFD',4),rep('MEC1',4)))

colData$group = as.factor(colData$group)


library(pheatmap)
pheatmap(count_data2,
              annotation_col = colData,
              cluster_rows = F,
              color = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
              cluster_cols = F,
              treeheight_row = 0,
              treeheight_col = 0,
              border_color = 'white',
              show_rownames=T,
              # main = 'MITOPHAGY',
              angle_col = 45,
         cellwidth = 20,cellheight = 20,
              fontsize_row = 8,fontsize_col = 8)
