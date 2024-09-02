library(Seurat)
library(tidyverse)
library(dplyr)
library(cowplot)
library(monocle) 

# 
# scRNA = readRDS("final_macrophages.rds")
# head(scRNA)
# table(scRNA@meta.data$cell_type)

# 
# sub_scRNA = subset(scRNA, cell_type %in% c("FCAR+ Macrophage cells","FCGR3A+ Macrophage cells"))
# sub_scRNA@meta.data$cell_type = droplevels(sub_scRNA@meta.data$cell_type)
# saveRDS(sub_scRNA,"FCAR+_FCGR3A+ Macrophage cells.rds")
# table(sub_scRNA@meta.data$cell_type)

scRNA = readRDS("FCAR+_FCGR3A+ Macrophage cells.rds")

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
diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~group", cores = 8)
head(diff)
# 
deg = subset(diff, qval < 0.01)
deg = deg[order(deg$qval, decreasing = F),]

# 
# write.csv(deg, "deg.csv")

# 
ordergene = rownames(deg)
cds = setOrderingFilter(cds, ordergene)

#
# pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
# dev.off()

# 
# 
cds = reduceDimension(cds, max_components = 2, method = "DDRTree")

# 
# 
cds = orderCells(cds) # 

## 
plot_cell_trajectory(cds, color_by = "Pseudotime", size = 1, show_backbone = T)

colors = c("#dacb79","#b0d794","#83c1b2","#f0c670","#f4a971","#f7a894","#a97b87")
plot_cell_trajectory(cds, color_by = "State", size = 1, show_backbone = T) +
  scale_color_manual(values = colors) +
  labs(title = "FCAR+ Macrophage cells")

data =pData(cds)
head(data)

data = data %>% 
  group_by(State,group) %>% 
  count() %>% 
  group_by(State) %>% 
  mutate(percent = 100 * n/sum(n))
data


ggplot(data,aes(State, percent, fill = group)) +
  geom_col() +
  theme_bw() +
  scale_fill_manual(values = c("#bbded7", "#e53313"))






# 
colors = c("#bbded7", "#e53313")
plot_cell_trajectory(cds, color_by = "group", size = 1, show_backbone = T) +
  scale_color_manual(values = colors) +
  labs(title = "FCAR+FCGR3A+ Macrophage cells")


# 
Time_diff = differentialGeneTest(cds[deg$gene], cores = 8,
                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")

Time_genes = Time_diff %>% 
  subset(use_for_ordering == "TRUE")%>%
  filter(pval < 0.05) %>%
  pull(gene_short_name) %>%
  as.character()
Time_genes = unique(Time_genes)
length(Time_genes)
p3 = plot_pseudotime_heatmap(cds[Time_genes,], 
                             num_clusters = 3, 
                             show_rownames = F, 
                             return_heatmap = T)

# 

clustering = data.frame(clusters)
head(clustering)
clustering[,1] = as.character(clustering[,1])
colnames(clustering) = "Gene_Clusters"
table(clustering)
# write.csv(clustering, "Time_clustering_all.csv")


# filer)
library(org.Hs.eg.db)

# 
bp_go <- function(gene) {
  enrichGO(
    gene,
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
}

# cluster1
cluster1_gene = rownames(subset(clustering, Gene_Clusters == "1"))
# cluster2
cluster2_gene = rownames(subset(clustering, Gene_Clusters == "2"))
# cluster3
cluster3_gene = rownames(subset(clustering, Gene_Clusters == "3"))

# 
bp_cluster1 = bp_go(cluster1_gene)
bp_cluster2 = bp_go(cluster2_gene)
bp_cluster3 = bp_go(cluster3_gene)

# 
write.table(bp_cluster1,"bp_cluster1.txt", sep = "\t")
write.table(bp_cluster2,"bp_cluster2.txt",sep = "\t")
write.table(bp_cluster3,"bp_cluster3.txt",sep = "\t")


# 
p4 = plot_genes_branched_heatmap(cds[Time_genes,], branch_point = 1,num_clusters = 3, show_rownames = F,return_heatmap = T)
p4$ph_res


# 

clustering = data.frame(clusters)
clustering[,1] = as.character(clustering[,1])
colnames(clustering) = "Gene_Clusters"
table(clustering)
# write.csv(clustering, "Time_clustering_branch.csv")

# cluster1
cluster1_gene = rownames(subset(clustering, Gene_Clusters == "1"))
# cluster2
cluster2_gene = rownames(subset(clustering, Gene_Clusters == "2"))
# cluster3
cluster3_gene = rownames(subset(clustering, Gene_Clusters == "3"))
# 
bp_cluster1 = bp_go(cluster1_gene)
bp_cluster2 = bp_go(cluster2_gene)
bp_cluster3 = bp_go(cluster3_gene)

# 
# write.table(bp_cluster1,"bp_cluster1_branch1.txt", sep = "\t")
# write.table(bp_cluster2,"bp_cluster2_brach2.txt",sep = "\t")
# write.table(bp_cluster3,"bp_cluster3_brach3.txt",sep = "\t")


# 
# PLK2,RHOH,ARHGAP5,LAT,DNAJC27,ARHGAP4,KANK1,IQSEC3,PSD3,MADD,RREB1,RALGDS,ELMO1,MYO9B,TRIM28,CDKN2A,RALGPS2,NUP62,RIPOR1,ABL1,LPAR4,RAB27A,RAB2B,PREX1,SCAI,FBXO8,NET1,LZTR1,OGT,ROCK1,RASAL3,ARHGEF40,DLC1,RASA3
# PLK2,E2F1,ATM,USP10,STK11,ATRX,NBN,CDKN1B,PERP,AKT1,ATAD5,RBBP7,AURKB,TP73,SIRT1,PIDD1,CNOT3,PCNA,KAT6A,EHMT2,BTG2,PRKAG2,BANP,MYBBP1A


# 
# 'ARHGEF40','FCAR','FCGR3A','MYO9B',
# 'PREX1','RAB27A','ROCK1','ATM',
# 'ATRX','BTG2','KAT6A','PRKAG2'

genes = c('ARHGEF40','FCAR','FCGR3A','MYO9B',
          'PREX1','RAB27A','ROCK1','ATM',
          'ATRX','BTG2','KAT6A','PRKAG2')
length(genes) # 12

plot_genes_branched_pseudotime(cds[genes,],
                                    branch_point = 2,
                                    color_by = "group",
                                    branch_labels=c("Cell fate 1", "Cell fate 2"),
                                    ncol = 4) +
  scale_color_manual(values = c("#AEC7E8", "#FF9896"))

