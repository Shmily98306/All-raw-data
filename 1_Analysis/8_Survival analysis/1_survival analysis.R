library(TCGAbiolinks) #
library(SummarizedExperiment)
library(tidyverse)
library(limma)
library(SummarizedExperiment)
# BiocManager::install("EDASeq")
library(EDASeq)  #
library(edgeR) #


# 
query <- GDCquery(
  project = c("TCGA-LAML"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
# 
GDCdownload(query)
# 
data = GDCprepare(query)
# 
saveRDS(data, "data.rds")



# 
data = readRDS("data.rds")

data_fpkm <- assay(data, i = 5) %>% 
  as.data.frame()

anno <- rowData(data)
anno_info = anno %>% 
  as.data.frame() %>% 
  subset(gene_type == "protein_coding") %>% 
  select(gene_id, gene_name)
head(anno_info)

data_fpkm[1:5,1:3]

# 
data_fpkm$gene_id = rownames(data_fpkm)

expr = merge(data_fpkm, anno_info, by = "gene_id")
expr = expr %>% 
  select(-gene_id) %>% 
  select(gene_name, everything())
expr[1:3,1:5]
head(expr)

# 
expr_mean=aggregate(.~gene_name,mean,data=expr)

# 
# write.csv(expr_mean, "expr.csv", row.names = F)

# 
expr =  read.csv("expr.csv", header = T, row.names = 1)
length(colnames(expr)) # 

# 
pd <- colData(data)
# 
pd_info  = as.data.frame(pd)
# 
barcode = pd_info$barcode
length(barcode) # 

# TB: Primary Blood Derived Cancer-Peripheral Blood
samplesTB <- TCGAquery_SampleTypes(barcode = barcode,
                                   typesample = "TB")

length(samplesTB) # 151个
# 

expr[1:3,1:5]







# 'ARHGEF40','FCAR','FCGR3A','MYO9B','PREX1','RAB27A','ROCK1','ATM','ATRX','BTG2','KAT6A','PRKAG2'
# "HBB","RPS4Y1","HBA2","HBA1","CXCR4","RPS3A","VIM","CMC1","RGCC","CXCL8","S100A9","S100A8","RPS26","IFITM1","SOD2","WDR74","MYOM2","GNLY","IFNG",
# "SNHG9","IER2","ID2","C1orf162","SAT1","FOSB","FCER1G",
# "ACTB","CITED2","JUN","TYROBP","PLEK","CD247","IFITM2","PPP1R15A","FGFBP2"
genes = c('ARHGEF40','FCAR','FCGR3A','MYO9B',
          'PREX1','RAB27A','ROCK1','ATM',
          'ATRX','BTG2','KAT6A','PRKAG2',
          "HBB","RPS4Y1","HBA2","HBA1","CXCR4","RPS3A","VIM","CMC1","RGCC","CXCL8",
          "S100A9","S100A8","RPS26","IFITM1","SOD2","WDR74","MYOM2","GNLY","IFNG",
          "SNHG9","IER2","ID2","C1orf162","SAT1","FOSB","FCER1G",
          "ACTB","CITED2","JUN","TYROBP","PLEK","CD247","IFITM2","PPP1R15A","FGFBP2")

length(unique(genes)) # 47
expr_select = expr[rownames(expr) %in% genes,]
dim(expr_select) # 46 151

# 
data = expr_select %>% 
  t() %>% 
  as.data.frame()
data$sample = substr(rownames(data), 1,16)
head(data)

# 
sur = read.delim("TCGA-LAML.survival.tsv")
head(sur)
sur = select(sur, c(sample,OS,OS.time))
colnames(sur) = c("sample","status","time")
sur$sample = gsub("-",".",sur$sample)
head(sur)

df = merge(data,sur, by = "sample")
head(df)

# 
# write.csv(df, "data for sur analysis.csv", row.names = F)







###########################################################################
library(survival)
library(survminer)

# 
df = read.csv("data for sur analysis.csv", header = T, row.names = 1)
head(df)

dim(df) # 132   48
for (i in colnames(df)[1:46]){
  name = paste0(i,"_","group")
  # value = data[,i]
  df[, name] = ifelse(df[,i] > median(df[,i]), "high", "low")
}
head(df)


# 
data_value = data.frame()
for (i in colnames(df)[1:46]){
  name = paste0(i,"_","group")
  formula_str <- paste("Surv(time, status) ~", name)
  cox_model1 <- coxph(as.formula(formula_str), data = df)
  data_value[i,"pvalue"] = summary(cox_model1)$coefficients[,5]
}
head(data_value)
dim(data_value) # 
data_value_diff = filter(data_value, pvalue < 0.05)
dim(data_value_diff) #
head(data_value_diff)


# 
gene = c("HBA1","PREX1","S100A8","S100A9")
for (i in gene){
  name = paste0(i,"_","group")
  formula_str <- paste("Surv(time, status) ~", name)
  cox_model2 <- survfit(as.formula(formula_str), data = df)
  title = paste0(i," Survival Curve")
  p = ggsurvplot(cox_model2, # 
                 data = df,  # 
                 conf.int = TRUE, # 
                 pval = TRUE, # 
                 risk.table = TRUE, # 
                 surv.median.line = "hv", # 
                 # add.all = TRUE, # 
                 palette = "hue") +  # 
    ggtitle(title)
  print(p)
} 


# https://ashpublications.org/blood/article/122/21/2610/11805/High-Transcription-Levels-Of-S100A8-and-S100A9-In
