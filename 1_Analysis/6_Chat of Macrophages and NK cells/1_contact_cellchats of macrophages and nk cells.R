# 
library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
options(stringAsFactors = F)

# 
macro = readRDS("FCAR+_FCGR3A+ Macrophage cells.rds")
nk_cells = readRDS("nk_cells.rds")

# 
scRNA = merge(macro,nk_cells)


# 
cellchat = createCellChat(object = scRNA, group.by = "cell_type")
cellchat
summary(cellchat)
str(cellchat)
levels(cellchat@idents)

# 
groupSize = as.numeric(table(cellchat@idents))
groupSize

# 
CellChatDB = CellChatDB.human
# CellChatDB = CellChatDB.mouse
str(CellChatDB)

# interaction, complex, cofactorgeneInfo dataframe
colnames((CellChatDB$interaction))
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
# dev.new() #
showDatabaseCategory(CellChatDB)
unique(CellChatDB$interaction$annotation) # 
# “Cell-Cell Contact" 
# Secreted Signaling, ECM-Receptor, Cell-Cell Contact
CellChatDB.use = subsetDB(CellChatDB, search = "Cell-Cell Contact")
cellchat@DB = CellChatDB.use

# 
# 
cellchat = subsetData(cellchat)
# future::plan("multiprocess", workers = 4)
# 
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)
# 
cellchat = projectData(cellchat, PPI.human)


# 
unique(cellchat@idents)
# 
cellchat = computeCommunProb(cellchat, raw.use = F, population.size = T) #

cellchat = filterCommunication(cellchat, min.cells = 10)
df.net = subsetCommunication(cellchat)
# write.csv(df.net, "net_lr.csv", row.names = F)

# 
cellchat = computeCommunProbPathway(cellchat)
df.netp = subsetCommunication(cellchat, slot.name = "netP")
# write.csv(df.netp,"net_patway.csv", row.names = F)
# 

# 
# 
# 
cellchat = aggregateNet(cellchat)
#
groupSize = as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd = T)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Interaction weights/strength")


#
mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  # i = 1
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}


mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=T)
for (i in 1:nrow(mat)) {
  # i = 1
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}


# 
levels(cellchat@idents)
# show all the significant interactions (L-R pairs)
#
p1 = netVisual_bubble(cellchat, sources.use = c(1,2),
                      targets.use = c(3), remove.isolate = FALSE)
p1 + theme(axis.text.x = element_text(angle = 15, vjust = 1))

p2 = netVisual_bubble(cellchat, sources.use = c(3),
                      targets.use = c(1,2), remove.isolate = FALSE)
p2 + theme(axis.text.x = element_text(angle = 15, vjust = 1))
