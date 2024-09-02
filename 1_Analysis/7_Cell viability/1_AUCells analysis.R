# 
library(tidyverse)
library(Seurat)
library(cols4all)
library(glmGamPoi) # 
library(cowplot)
library(clustree)
library(AUCell)
library(clusterProfiler) # 
library(ggsignif)
library(GSVA)
library(ggpubr)


# 
scRNA = readRDS("final_t_cells.rds")
table(scRNA@meta.data$cell_type)

#Build gene expression rankings for each cell
CellRank<- AUCell_buildRankings(as.matrix(scRNA@assays$RNA@data))

# load gene set
go = load("go_data_list.rdata")
head(data_list)

#Calculates the 'AUC' for each gene-set in each cell.
cells_AUC <- AUCell_calcAUC(data_list, CellRank,nCores = 5, aucMaxRank=nrow(CellRank)*0.05)

aucs <- getAUC(cells_AUC)
head(aucs)

# t(aucs) %>%
#   as.data.frame() %>%
#   write.csv("t_cells aucs.csv")

aucs = read.csv("t_cells aucs.csv",header = T)
head(aucs)

head(scRNA)
group = select(scRNA@meta.data,c(group,cell_type))
group$X = rownames(group)
head(group)

aucs = merge(aucs, group, by = "X")
head(aucs)

colors = c("#bbded7","#e53313")
p1 = ggplot(aucs,aes(reorder(cell_type,-natural.killer.cell.activation), natural.killer.cell.activation, fill = group)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = colors) +
  labs(x = "Natural killer cell activation", y = "AUCells Score") +
  theme(legend.position = "top") +
  stat_compare_means(aes(group = group),label = "p.signif",size=3,method = "t.test")


data = filter(aucs, natural.killer.cell.proliferation > 0)
p2 = ggplot(data,aes(reorder(cell_type,-natural.killer.cell.proliferation), natural.killer.cell.proliferation, fill = group)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = colors) +
  labs(x = "Natural killer cell proliferation", y = "AUCells Score") +
  theme(legend.position = "top") +
  stat_compare_means(aes(group = group),label = "p.format",size=3,method = "t.test")









#############macrophage cells AUCell analysis#########################################################
# 
scRNA = readRDS("final_macrophages.rds")
table(scRNA@meta.data$cell_type)

#Build gene expression rankings for each cell
CellRank<- AUCell_buildRankings(as.matrix(scRNA@assays$RNA@data))

# load gene set
go = load("go_data_list.rdata")
head(data_list)

#Calculates the 'AUC' for each gene-set in each cell.
cells_AUC <- AUCell_calcAUC(data_list, CellRank,nCores = 5, aucMaxRank=nrow(CellRank)*0.05)

aucs <- getAUC(cells_AUC)
head(aucs)

t(aucs) %>%
  as.data.frame() %>%
  write.csv("macro_cells aucs.csv")

aucs = read.csv("macro_cells aucs.csv",header = T)
head(aucs)

head(scRNA)
group = select(scRNA@meta.data,c(group,cell_type))
group$X = rownames(group)
head(group)

aucs = merge(aucs, group, by = "X")
head(aucs)

colors = c("#bbded7","#e53313")
p3 = ggplot(aucs,aes(reorder(cell_type,-macrophage.activation), macrophage.activation, fill = group)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = colors) +
  labs(x = "macrophage activation", y = "AUCells Score") +
  theme(legend.position = "top") +
  stat_compare_means(aes(group = group),label = "p.signif",size=3,method = "t.test")


data = filter(aucs, macrophage.migration > 0)
p4 = ggplot(data,aes(reorder(cell_type,-macrophage.migration), macrophage.migration, fill = group)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = colors) +
  labs(x = "macrophage migration", y = "AUCells Score") +
  theme(legend.position = "top") +
  stat_compare_means(aes(group = group),label = "p.signif",size=3,method = "t.test")

plot_grid(p1,p2,p3,p4, ncol = 2, labels = c("A","B","C","D"))
