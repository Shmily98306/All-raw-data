# 
library(tidyverse)
library(Seurat)
library(patchwork)
library(harmony)
library(cols4all)
library(glmGamPoi) # SCT
library(cowplot)
library(clustree)


# 
# scRNA = readRDS("final_scRNA.rds")
# head(scRNA)
# macrophage cells
# macrophages = subset(scRNA, cell_type == "Macrophage cells")
# macrophages@meta.data$cell_type = droplevels(macrophages@meta.data$cell_type)
# table(macrophages@meta.data$cell_type)
# saveRDS(macrophages, "macrophages.rds")

# 
scRNA = readRDS("macrophages.rds")

# 
scRNA = SCTransform(scRNA, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = F)
# saveRDS(scRNA, file = "./sct_macrophages.rds")

#
scRNA = readRDS("sct_macrophages.rds")

# 
scRNA <- RunPCA(scRNA)
colnames(scRNA@meta.data)[1] = c("patient_ID", "nCount_RNA", "nFeature_RNA", "percent.mt") #
# 
scRNA <- RunHarmony(scRNA, group.by.vars="patient_ID", max.iter.harmony=20, lambda=0.5, assay.use = "SCT")
# 
ElbowPlot(scRNA,ndims = 50)

# 
scRNA <- FindNeighbors(scRNA, dims=1:15, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:15, reduction="harmony")

# 
DimPlot(scRNA, reduction="umap", group.by="patient_ID", pt.size=1)+
  theme(legend.position="right", plot.title=element_blank())

obj = FindClusters(scRNA, resolution = seq(0.1, 1,by=0.1))
p4 = clustree(obj)
p4

# resolution = 0.1
scRNA <- FindClusters(scRNA, resolution=0.1)
p5 = UMAPPlot(scRNA, pt.size=1, label=T)+NoLegend()
p5

# 
# scRNA.markers <- FindAllMarkers(scRNA, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
# write.csv(scRNA.markers, "markers.csv")



# 
VlnPlot(scRNA, features=c("FCAR"), pt.size=0,group.by = "seurat_clusters")+
  NoLegend()+
  theme(axis.title.x=element_blank())

# cluster 0: FCAR+ Macrophage cells: FCAR
# cluster 1: IL32+ Macrophage cells: IL32
# cluster 2: ELANE+ Macrophage cells: ELANE
# cluster 3: FCGR3A+ Macrophage cells: FCGR3A

# cluster 
cell_label = c("FCAR+ Macrophage cells","IL32+ Macrophage cells",
               "ELANE+ Macrophage cells","FCGR3A+ Macrophage cells")
## 
names(cell_label) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, cell_label)
scRNA[["cell_type"]] = Idents(scRNA)

# 
colors = c("#00adad","#e5d4c6","#e73981","#acb8b4")
p8 = UMAPPlot(scRNA, pt.size=1, label=T, label.size=5, cols = colors)+
  NoLegend()
p8

anno_markers = c("FCAR","IL32","ELANE","FCGR3A")

p9 = DotPlot(scRNA, features=anno_markers, cols=c("#bbded7", "#e53313"))+coord_flip()+
  theme(
    axis.text.x=element_text(angle=30, hjust=1, size=10, face="bold"), 
    axis.text.y=element_text(face="bold", size=12), 
    axis.title.x=element_blank(), 
    axis.title.y=element_blank()
  )
p9

# 
p10 = VlnPlot(scRNA, features=anno_markers, pt.size=0, ncol = 2, cols = colors)+
  NoLegend()+
  theme(axis.title.x=element_blank())
p10

scRNA = readRDS("final_macrophages.rds")

# 
cell_count = scRNA@meta.data %>%
  group_by(group, cell_type) %>%
  count() %>% 
  group_by(cell_type) %>% 
  mutate(Percent=n/sum(n))
head(cell_count)
write.csv(cell_count, "Table 3.csv", row.names = F)

ggplot(cell_count, aes(reorder(cell_type, -Percent, median), Percent*100, fill = group)) +
  geom_col(position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("#bbded7", "#e53313")) +
  scale_color_manual() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1),legend.position = "top") +
  labs(title = "Cell Ratio", x = "", y = "Ratio (%)")

write.csv(cell_count, "Table 3.csv", row.names = F)

# saveRDS(scRNA,"final_macrophages.rds")





