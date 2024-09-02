library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)


bp_facr = read.csv("bp_facr_select.csv", header = T, row.names = 1)
bp_il32 = read.csv("bp_il32_select.csv", header = T, row.names = 1)
bp_elane = read.csv("bp_elane_select.csv", header = T, row.names = 1)
bp_fcgr3a = read.csv("bp_fcgr3a_select.csv", header = T, row.names = 1)

bp_facr$cell_type = "FACR+ Macrophage cells"
bp_il32$cell_type = "IL32+ Macrophage cells"
bp_elane$cell_type = "ELANE+ Macrophage cells"
bp_fcgr3a$cell_type = "FCGR3A+ Macrophage cells"

data = rbind(bp_facr,bp_il32,bp_elane,bp_fcgr3a)
dim(data)
head(data)

ggplot(data, aes(cell_type, reorder(Description,Count), size = Count, color = pvalue)) +
  geom_point() +
  theme_bw() +
  # facet_grid(.~cell_type) +
  scale_color_gradient(low = "#e53313",high = "#bbded7") +
  scale_fill_gradient(low = "#e53313",high = "#bbded7") +
  theme(axis.text.x = element_text(angle = 15,hjust = 1))+
  labs(title = "Biological process of Macrophage cells", x = "Count", y = "Description")
