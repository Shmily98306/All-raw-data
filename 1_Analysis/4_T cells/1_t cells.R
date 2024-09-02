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
# scRNA = readRDS("final_scRNA.rds")
# table(scRNA@meta.data$cell_type)
# 
# t_cells = subset(scRNA, cell_type == "NK/T cells")
# t_cells@meta.data$cell_type = droplevels(t_cells@meta.data$cell_type)
# table(t_cells@meta.data$cell_type)
# saveRDS(t_cells, "t_cells.rds")

# 
scRNA = readRDS("t_cells.rds")

# 
# scRNA = SCTransform(scRNA, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T)
# saveRDS(scRNA, file = "./sct_t_cells.rds")


# 
scRNA = readRDS("sct_t_cells.rds")

# 
scRNA <- RunPCA(scRNA)
colnames(scRNA@meta.data)[1] = c("patient_ID", "nCount_RNA", "nFeature_RNA", "percent.mt") #
# 
scRNA <- RunHarmony(scRNA, group.by.vars="patient_ID", max.iter.harmony=20, lambda=0.5, assay.use = "SCT")
# 
ElbowPlot(scRNA,ndims = 50)

# 
scRNA <- FindNeighbors(scRNA, dims=1:16, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:16, reduction="harmony")

# 
DimPlot(scRNA, reduction="umap", group.by="patient_ID", pt.size=1)+
  theme(legend.position="right", plot.title=element_blank())

# resolution = 0.1
scRNA <- FindClusters(scRNA, resolution=0.1)
p5 = UMAPPlot(scRNA, pt.size=1, label=T)+NoLegend()
p5
# 
# scRNA.markers <- FindAllMarkers(scRNA, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
# write.csv(scRNA.markers, "markers.csv")

VlnPlot(scRNA, features=c('MBOAT7','LAT2','SOX4','CAPG','GRN','CD82','HLA-DMA'), pt.size=0,group.by = "seurat_clusters")+
  NoLegend()+
  theme(axis.title.x=element_blank())


# cluster 0: Naive CD8+ T cells: TRABD2A,MAL,FAM102A,SH3YL1
# cluster 1: Natural killer cells: CLIC3,TRDC,KLRG1
# cluster 2: Treg cells: MBOAT7,SOX4,CAPG,GRN,CD82,HLA-DMA
# cluster 3: Natural killer cells
# cluster 4: Natural killer cells

# 
cell_label = c("Naive CD8+ T cells","Natural killer cells","Treg cells",
               "Natural killer cells","Natural killer cells")
## 
names(cell_label) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, cell_label)
scRNA[["cell_type"]] = Idents(scRNA)

colors = c("#dacb79","#b0d794","#f7a894")
p8 = UMAPPlot(scRNA, pt.size=1, label=T, label.size=5, cols = colors)+
  NoLegend()
p8

anno_markers = c('TRABD2A','MAL','FAM102A','SH3YL1',
                 'CLIC3','TRDC','KLRG1',
                 'MBOAT7','SOX4','CAPG','GRN','CD82','HLA-DMA')

p9 = DotPlot(scRNA, features=anno_markers, cols=c("#bbded7", "#e53313"))+coord_flip()+
  theme(
    axis.text.x=element_text(angle=15, hjust=1, size=10, face="bold"), 
    axis.text.y=element_text(face="bold", size=12), 
    axis.title.x=element_blank(), 
    axis.title.y=element_blank()
  )
p9

genes = c('TRABD2A','MAL','CLIC3','MBOAT7')
# 
p10 = VlnPlot(scRNA, features=genes, pt.size=0, ncol = 2, cols = colors)+
  NoLegend()+
  theme(axis.title.x=element_blank())
p10

scRNA = readRDS("final_t_cells.rds")
# 
cell_count = scRNA@meta.data %>%
  group_by(group, cell_type) %>%
  count() %>% 
  group_by(group) %>% 
  mutate(Percent=n/sum(n))
head(cell_count)
# write.csv(cell_count,"Table 4.csv", row.names = F)

ggplot(cell_count, aes(reorder(cell_type, -Percent, median), Percent*100, fill = group)) +
  geom_col(position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("#bbded7", "#e53313")) +
  scale_color_manual() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1),legend.position = "top") +
  labs(title = "Cell Ratio", x = "", y = "Ratio (%)")

# saveRDS(scRNA,"final_t_cells.rds")
