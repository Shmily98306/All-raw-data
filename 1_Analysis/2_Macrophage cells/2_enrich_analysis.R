library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)

marker = read.csv("markers.csv",header = T,row.names = 1)

# 
marker$group = case_when(
  marker$cluster == 0 ~ "FCAR+ Macrophage cells",
  marker$cluster == 1 ~ "IL32+ Macrophage cells",
  marker$cluster == 2 ~ "ELANE+ Macrophage cells",
  marker$cluster == 3 ~ "FCGR3A+ Macrophage cells",
)

head(marker)
table(marker$group)

# 
facr = subset(marker, group == "FCAR+ Macrophage cells")
il32 = subset(marker, group == "IL32+ Macrophage cells")
elane = subset(marker, group == "ELANE+ Macrophage cells")
fcgr3a = subset(marker, group == "FCGR3A+ Macrophage cells")

# 
facr_genes = facr$gene
il32_genes = il32$gene
elane_genes = elane$gene
fcgr3a_genes = fcgr3a$gene

# 
facr_id = bitr(facr_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID
il32_id = bitr(il32_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID
elane_id = bitr(elane_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID
fcgr3a_id = bitr(fcgr3a_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID

# 
bp_go <- function(gene) {
  enrichGO(gene,
           OrgDb = "org.Hs.eg.db",
           keyType = "ENTREZID",
           ont = "BP",
           pvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           qvalueCutoff = 0.2)
}

bp_facr = bp_go(facr_id)
bp_il32 = bp_go(il32_id)
bp_elane = bp_go(elane_id)
bp_fcgr3a = bp_go(fcgr3a_id)


# 
# write.csv(bp_facr,"bp_facr.csv")
# write.csv(bp_il32,"bp_il32.csv")
# write.csv(bp_elane,"bp_elane.csv")
# write.csv(bp_fcgr3a,"bp_fcgr3a.csv")



################################################################################################################
# 
bp_facr = read.csv("bp_facr_select.csv", header = T, row.names = 1)
bp_il32 = read.csv("bp_il32_select.csv", header = T, row.names = 1)
bp_elane = read.csv("bp_elane_select.csv", header = T, row.names = 1)
bp_fcgr3a = read.csv("bp_fcgr3a_select.csv", header = T, row.names = 1)

p1 = ggplot(bp_facr, aes(Count,reorder(Description,Count))) +
  geom_point(aes(color = pvalue, size = Count)) +
  geom_col(aes(fill = pvalue), width = 0.05) +
  theme_bw() +
  scale_color_gradient(low = "#e53313",high = "#bbded7") +
  scale_fill_gradient(low = "#e53313",high = "#bbded7") +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
  geom_text(aes(x = 0, y = Description,label = Description),vjust = -0.5,hjust = 0) +
  labs(title = "Biological process of FCAR+ Macrophage cells", x = "Count", y = "Description");p1


p2 = ggplot(bp_il32, aes(Count,reorder(Description,Count))) +
  geom_point(aes(color = pvalue, size = Count)) +
  geom_col(aes(fill = pvalue), width = 0.05) +
  theme_bw() +
  scale_color_gradient(low = "#e53313",high = "#bbded7") +
  scale_fill_gradient(low = "#e53313",high = "#bbded7") +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
  geom_text(aes(x = 0, y = Description,label = Description),vjust = -0.5,hjust = 0) +
  labs(title = "Biological process of IL32+ Macrophage cells", x = "Count", y = "Description");p2

p3 = ggplot(bp_elane, aes(Count,reorder(Description,Count))) +
  geom_point(aes(color = pvalue, size = Count)) +
  geom_col(aes(fill = pvalue), width = 0.05) +
  theme_bw() +
  scale_color_gradient(low = "#e53313",high = "#bbded7") +
  scale_fill_gradient(low = "#e53313",high = "#bbded7") +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
  geom_text(aes(x = 0, y = Description,label = Description),vjust = -0.5,hjust = 0) +
  labs(title = "Biological process of ELANE+ Macrophage cells", x = "Count", y = "Description");p3

p4 = ggplot(bp_fcgr3a, aes(Count,reorder(Description,Count))) +
  geom_point(aes(color = pvalue, size = Count)) +
  geom_col(aes(fill = pvalue), width = 0.05) +
  theme_bw() +
  scale_color_gradient(low = "#e53313",high = "#bbded7") +
  scale_fill_gradient(low = "#e53313",high = "#bbded7") +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
  geom_text(aes(x = 0, y = Description,label = Description),vjust = -0.5,hjust = 0) +
  labs(title = "Biological process of FCGR3A+ Macrophage cells", x = "Count", y = "Description");p3
plot_grid(p1,p2,p3,p4, labels = c("E","F","G","H"), ncol = 2)
