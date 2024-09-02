library(Seurat)
library(tidyverse)
library(dplyr)

# 
scRNA = readRDS("nk_cells.rds")


genes = c("HBB","RPS4Y1","HBA2","HBA1","CXCR4","RPS3A","VIM","CMC1","RGCC","CXCL8",
          "S100A9","S100A8","RPS26","IFITM1","SOD2","WDR74","MYOM2","GNLY","IFNG",
          "SNHG9","IER2","ID2","C1orf162","SAT1","FOSB","FCER1G",
          "ACTB","CITED2","JUN","TYROBP","PLEK","CD247","IFITM2","PPP1R15A","FGFBP2")


DotPlot(scRNA, features=genes, cols=c("#bbded7", "#e53313"),group.by = "group")+coord_flip()+
  theme_bw() +
  theme(
    axis.text.x=element_text(angle=30, hjust=1, size=10, face="bold"), 
    axis.text.y=element_text(face="bold", size=12), 
    axis.title.x=element_blank(), 
    axis.title.y=element_blank()
  )
