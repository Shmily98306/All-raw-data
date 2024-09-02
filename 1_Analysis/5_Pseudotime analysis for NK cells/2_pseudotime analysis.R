library(Seurat)
library(tidyverse)
library(dplyr)
library(cowplot)
library(monocle) #
library(cols4all)


# 
# scRNA =  read_rds("final_t_cells.rds")
# head(scRNA)
# 
# nk_cells = subset(scRNA, cell_type == "Natural killer cells")
# nk_cells@meta.data$cell_type = droplevels(nk_cells@meta.data$cell_type)
# 
# saveRDS(nk_cells, "nk_cells.rds")

 # 
scRNA = readRDS("nk_cells.rds")

# 
expr_matrix = as(as.matrix(scRNA@assays$RNA@counts), 'sparseMatrix')
# 
p_data = scRNA@meta.data
p_data$cell_type = scRNA@active.ident # 
# 
f_data = data.frame(gene_short_name = row.names(expr_matrix), row.names = row.names(expr_matrix))
## （cell number）
dim(expr_matrix)
dim(p_data)
dim(f_data)

# 
pd = new('AnnotatedDataFrame', data = p_data)
fd = new('AnnotatedDataFrame', data = f_data)
# 
cds = newCellDataSet(expr_matrix,
                     phenoData = pd,
                     featureData = fd,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())

# 
cds = estimateSizeFactors(cds) # 
cds = estimateDispersions(cds)


# 
cds = detectGenes(cds, min_expr = 0.1) #num_cells_expressed
print(head(fData(cds)))
expressed_genes = row.names(subset(fData(cds),
                                   num_cells_expressed >= 10))

# 
# 
# 
head(scRNA)
Idents(scRNA) = scRNA@meta.data$group # 
deg.cluster = FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)
head(deg.cluster)
express_genes = rownames(subset(deg.cluster, p_val_adj < 0.05, abs(avg_log2FC) > 0.25))
head(express_genes)
cds = setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)

# 
cds = reduceDimension(cds, max_components = 2, method = "DDRTrees")
# 
cds = orderCells(cds) # 


# 
plot_cell_trajectory(cds, color_by = "Pseudotime", size = 1, show_backbone = T)

# 
plot_cell_trajectory(cds, color_by = "State", size = 1, show_backbone = T) +
  # scale_color_manual(values = colors) +
  labs(title = "Natural killer cells")

# 
colors = c("#bbded7","#e53313")
plot_cell_trajectory(cds, color_by = "group", size = 1, show_backbone = T) +
  scale_color_manual(values = colors) +
  labs(title = "Natural killer cells")


# 
Time_diff = differentialGeneTest(cds[deg.cluster$gene], cores = 1,
                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")
head(Time_diff)
Time_genes = Time_diff %>% filter(pval < 0.05) %>% pull(gene_short_name) %>%as.character()
Time_genes = unique(Time_genes)
Time_genes
genes = c("HBB","RPS4Y1","HBA2","HBA1","CXCR4","RPS3A","VIM","CMC1","RGCC","CXCL8",
          "S100A9","S100A8","RPS26","IFITM1","SOD2","WDR74","MYOM2","GNLY","IFNG",
          "SNHG9","IER2","ID2","C1orf162","SAT1","FOSB","FCER1G",
          "ACTB","CITED2","JUN","TYROBP","PLEK","CD247","IFITM2","PPP1R15A","FGFBP2")

p3 = plot_pseudotime_heatmap(cds[genes,], 
                             num_clusters = 3, 
                             show_rownames = T, 
                             return_heatmap = T)
p3


DotPlot(scRNA, features=genes, cols=c("#AEC7E8", "#FF9896"),)+coord_flip()+
  theme(
    axis.text.x=element_text(angle=30, hjust=1, size=10, face="bold"), 
    axis.text.y=element_text(face="bold", size=12), 
    axis.title.x=element_blank(), 
    axis.title.y=element_blank()
  )
