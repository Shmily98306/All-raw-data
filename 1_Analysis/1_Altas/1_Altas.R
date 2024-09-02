#
library(tidyverse)
library(Seurat)
library(patchwork)
library(harmony)
library(cols4all)
library(glmGamPoi) # 
library(cowplot)
library(clustree)

scRNA = readRDS("final_scRNA.rds")

# 
cell_count = scRNA@meta.data %>%
  group_by(group, cell_type) %>%
  count() %>% 
  group_by(group) %>% 
  mutate(Percent=n/sum(n))
head(cell_count)

ggplot(cell_count, aes(reorder(cell_type, -Percent, median), Percent*100, fill = group)) +
  geom_col(position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("#bbded7", "#e53313")) +
  scale_color_manual() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1),legend.position = "top") +
  labs(title = "Cell Ratio", x = "", y = "Ratio (%)")

# write.csv(cell_count, "Table 2.csv", row.names = F)


# 
head(scRNA)
cell_count = scRNA@meta.data %>%
  group_by(patient_ID, cell_type) %>%
  count() %>% 
  group_by(patient_ID) %>% 
  mutate(Percent=n/sum(n))
head(cell_count)
# write.csv(cell_count, "Table 1.csv", row.names = F)
