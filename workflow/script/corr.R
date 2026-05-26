library(tidyverse)
# 1.输入样本信息
sample_info <- read.table("../data/samples.txt") %>% 
  rename(group=V1,sample=V2) %>% 
  column_to_rownames(var = "sample")

# 2.导入表达矩阵
gene_exp <- read.table("../3.Merge_result/genes.TMM.EXPR.matrix",header = T,row.names = 1,check.names = F)

# 3.样本相关分析
sample_cor <- round(cor(gene_exp),digits = 5)

# 相关系数热图
library(pheatmap)
png(file = "heatmap_sample_relationship.png", width = 3000, height = 3000, res = 300)
pheatmap(sample_cor,
         cluster_cols = F,
         cluster_rows = F,
         cellwidth = 50,
         cellheight = 50,
         border_color = "white",
         fontsize = 10,
         angle_col = 45,
         display_numbers = T,
         fontsize_number = 10)
dev.off()
pdf(file = "heatmap_sample_relationship.pdf", width = 10, height = 10)
pheatmap(sample_cor,
         cluster_cols = F,
         cluster_rows = F,
         cellwidth = 50,
         cellheight = 50,
         border_color = "white",
         fontsize = 10,
         angle_col = 45,
         display_numbers = T,
         fontsize_number = 10)
dev.off()



# 计算样本之间的欧式距离
sample_dist <- dist(t(gene_exp))

# 根据距离矩阵进行层次聚类
png(file = "欧式距离层次聚类_样本关系.png", width = 3000, height = 2500, res = 300)
plot(hclust(sample_dist))
dev.off()
pdf(file = "欧式距离层次聚类_样本关系.pdf", width = 10, height = 8)
plot(hclust(sample_dist))
dev.off()





# 主成分分析
library(PCAtools)
pca <- pca(gene_exp,removeVar = 0.3,metadata = sample_info)
# 主成分解释度
png(file = "主成分解释度.png", width = 5000, height = 2500, res = 300)
screeplot(pca)
dev.off()

pdf(file = "主成分解释度.pdf", width = 10, height = 5)
screeplot(pca)
dev.off()

# PC1 PC2 二维图
png(file = "PC1.PC2二维图.png", width = 2700, height = 2700, res = 300)
biplot(pca,
       x="PC1",
       y="PC2",
       colby = "group",
       #legendPosition = "right",
)
dev.off()
pdf(file = "PC1.PC2二维图.pdf", width = 10, height = 10)
biplot(pca,
       x="PC1",
       y="PC2",
       colby = "group",
       #legendPosition = "right",
)
dev.off()

