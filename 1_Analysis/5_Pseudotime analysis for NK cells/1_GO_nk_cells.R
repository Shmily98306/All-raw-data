# 
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# 
deg = read.csv("markers.csv", header = T, row.names = 1)
head(deg)

genes = subset(deg, cluster %in% c(1,3,4))$gene
length(genes)

genes_id = bitr(genes,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID

bp_go = enrichGO(genes_id,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 ont = "BP",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.2)
# write.csv(bp_go,"bp_go.csv")


data = read.csv("bp_go_select.csv",header = T, row.names = 1)
head(data)

ggplot(data, aes(Count,reorder(Description, Count), size = Count, color = pvalue)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low = "#e53313", high = "#bbded7") +
  labs(title = "Biological process of Natural killer cells", x = "Count", y = "Description")
