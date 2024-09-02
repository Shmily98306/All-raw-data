# 
library(tidyverse)
library(spatstat.utils)
library(Seurat)
library(patchwork)
library(harmony)
library(cols4all)
library(glmGamPoi) # 
library(cowplot)
library(clustree)
options(future.globals.maxSize= 1000*1024^2)

## 
### 
assays <- dir("./GSE198681_RAW/")
dir <- paste0("./GSE198681_RAW/", assays)
# 
samples_name = assays
samples_name

  
# 
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir=dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts$`Gene Expression`, project=samples_name[i], min.cells=3, min.features=200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id=samples_name[i])   
  if(T){
    scRNAlist[[i]][["percent.mito"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern="^MT-")
  }
}

### 
names(scRNAlist) <- samples_name

# 
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
head(scRNA)

# 
# saveRDS(scRNA, file = "./scRNA.rds")



# 
scRNA@meta.data$sample = substr(rownames(scRNA@meta.data),1,10)
scRNA@meta.data$group = case_when(
  scRNA@meta.data$sample == "GSM5954920" ~ "CR",
  scRNA@meta.data$sample == "GSM5954922" ~ "CR",
  scRNA@meta.data$sample == "GSM5954924" ~ "CR",
  scRNA@meta.data$sample == "GSM5954926" ~ "non-CR",
  scRNA@meta.data$sample == "GSM5954928" ~ "non-CR"
)
head(scRNA@meta.data)

#
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
head(scRNA@meta.data)

VlnPlot(
  scRNA, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
  # cols = colors,
  pt.size = 0.1, 
  ncol = 3
)

# 
scRNA = subset(scRNA, subset=nFeature_RNA>200&nFeature_RNA<7500&percent.mt<50)
VlnPlot(
  scRNA, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  # cols = colors,
  pt.size = 0, 
  ncol = 3
)


# 
scRNA = SCTransform(scRNA, vars.to.regress = "percent.mt", verbose = F)
saveRDS(scRNA, file = "./sct_scRNA.rds")

# 
scRNA <- RunPCA(scRNA)
colnames(scRNA@meta.data)[1] = c("patient_ID", "nCount_RNA", "nFeature_RNA", "percent.mt") #
head(scRNA@meta.data)

# 
scRNA <- RunHarmony(scRNA, group.by.vars="patient_ID", max.iter.harmony=20, lambda=0.5, assay.use = "SCT")

# 
ElbowPlot(scRNA,ndims = 50)

# 
scRNA <- FindNeighbors(scRNA, dims=1:20, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:20, reduction="harmony")










