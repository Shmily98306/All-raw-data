#
library(tidyverse)
library(Seurat)
library(patchwork)
library(harmony)
library(cols4all)
library(glmGamPoi) # 
library(cowplot)
library(clustree)




# 
scRNA = readRDS("sct_scRNA.rds")

# 
scRNA <- RunPCA(scRNA)
colnames(scRNA@meta.data)[1] = c("patient_ID", "nCount_RNA", "nFeature_RNA", "percent.mt") #
# 
scRNA <- RunHarmony(scRNA, group.by.vars="patient_ID", max.iter.harmony=20, lambda=0.5, assay.use = "SCT")
# 
ElbowPlot(scRNA,ndims = 50)

# 
scRNA <- FindNeighbors(scRNA, dims=1:30, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:30, reduction="harmony")

DimPlot(scRNA, reduction="umap", group.by="patient_ID", pt.size=1)+
  theme(legend.position="right", plot.title=element_blank())

# obj = FindClusters(scRNA, resolution = seq(0.1, 1,by=0.1))
# clustree(obj)

# resolution = 0.1
scRNA <- FindClusters(scRNA, resolution=0.1)
UMAPPlot(scRNA, pt.size=1, label=T)+NoLegend()

# 
# scRNA.markers <- FindAllMarkers(scRNA, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
# write.csv(scRNA.markers, "markers.csv")

VlnPlot(scRNA, features=c("AZU1","ELANE","MPO"), pt.size=0,group.by = "seurat_clusters")+
  NoLegend()+
  theme(axis.title.x=element_blank())

# cluster 0: Mast cells: CTSG, HPGD, IGLL1, SMYD3
# cluster 1: Macrophage cells: CD86, MS4A7, CD163
# cluster 2: Mast cells: 
# cluster 3: NK/T cells: NKG7,GZMA,CCL5,GZMB,GZMH,GZMK, CD3D, CD3G, CD3E, CD2
# cluster 4: NK/T cells:
# cluster 5: Red blood cells: ALAS2, GYPA, RHCE
# cluster 6: Red blood cells
# cluster 7: Mast cells:
# cluster 8: B cells: CD79A, MS4A1
# cluster 9: Mast cells:
# cluster 10: Macrophage cells:
# cluster 11: Megakaryocyte: ITGA2B,PF4,TUBB1,GP9
# cluster 12: Mast cells:
# cluster 13: Mast cells:

sct_scRNA = scRNA

mast_cell = c('CTSG', 'HPGD', 'IGLL1', 'SMYD3')
macro_cells = c('CD86', 'MS4A7', 'CD163')
nk_t = c('NKG7','GZMA','CCL5','GZMB','GZMH','GZMK', 'CD3D', 'CD3G', 'CD3E', 'CD2')
rbc = c('ALAS2', 'GYPA', 'RHCE')
b_cells = c('CD79A', 'MS4A1')
meg = c('ITGA2B','PF4','TUBB1','GP9')

anno_markers = c(mast_cell,macro_cells,nk_t,rbc,b_cells,meg)

cell_label = c("Mast cells","Macrophage cells","Mast cells","NK/T cells","NK/T cells",
               "Red blood cells","Red blood cells","Mast cells","B cells","Mast cells","Macrophage cells",
               "Megakaryocyte","Mast cells","Mast cells")

## 
names(cell_label) <- levels(sct_scRNA)
sct_scRNA <- RenameIdents(sct_scRNA, cell_label)
sct_scRNA[["cell_type"]] = Idents(sct_scRNA)

# 
colors = c("#ff85a8","#ffd685","#fff385","#85ffd0","#85c2ff","#ff85f3")
p8 = UMAPPlot(sct_scRNA, pt.size=1, label=T, label.size=5, cols = colors)+
  NoLegend()
p8

# 
p9 = DotPlot(sct_scRNA, features=anno_markers, cols=c("#bbded7", "#e53313"))+coord_flip()+
  theme_bw() +
  theme(
    axis.text.x=element_text(angle=30, hjust=1, size=10, face="bold"), 
    axis.text.y=element_text(face="bold", size=12), 
    axis.title.x=element_blank(), 
    axis.title.y=element_blank()
  )
p9


genes = c("SMYD3","MS4A7","GZMA","ALAS2","CD79A","ITGA2B")
# 
p10 = VlnPlot(sct_scRNA, features=genes, pt.size=0, ncol = 3, cols = colors)+
  NoLegend()+
  theme(axis.title.x=element_blank())
p10


# 
cell_count = sct_scRNA@meta.data %>%
  group_by(patient_ID, cell_type) %>%
  count() %>% 
  group_by(patient_ID) %>% 
  mutate(Percent=n/sum(n) * 100)
head(cell_count)

ggplot(cell_count, aes(patient_ID, Percent*100)) +
  geom_col(aes(fill = cell_type)) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  scale_color_manual() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  labs(title = "Cell Ratio", x = "", y = "Ratio (%)")


# 
cell_count1 = sct_scRNA@meta.data %>%
  group_by(cell_type, group) %>%
  count() %>%
  group_by(cell_type) %>% 
  mutate(Percent=n/sum(n) * 100)
cell_count1

ggplot(cell_count1, aes(reorder(cell_type, -Percent), Percent, fill = group)) +
  geom_col(position = "dodge") +
  theme_bw() + 
  scale_fill_manual(values = c("#bbded7","#e53313")) +
  scale_color_manual() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1), legend.position = "top") +
  labs(title = "Cell Ratio", x = "", y = "Ratio (%)")


# saveRDS(sct_scRNA, "final_scRNA.rds")
