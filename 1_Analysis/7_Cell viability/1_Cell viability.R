library(tidyverse)
library(org.Hs.eg.db)


k <- keys(org.Hs.eg.db, keytype = "SYMBOL")
head(k)

# 
go_id = select(org.Hs.eg.db,
               keys = k,
               columns = c("GO", "ONTOLOGY"),
               keytype="SYMBOL")

dim(go_id)
head(go_id)


# 
# library(GO.db)
# goterms <- Term(GOTERM)
# GOlist=as.data.frame(goterms)
# head(GOlist)
# write.csv(GOlist,"GOlist.csv")

go_data = read.csv("GO_select.csv", header = T)
head(go_data)

# 

go_id_select = subset(go_id, GO %in% go_data$GO)
head(go_id_select)

data = merge(go_data,go_id_select,by = "GO")
head(data)

data = dplyr::select(data,c(goterms,SYMBOL))
head(data)

data_list = split(data$SYMBOL,data$goterms)
head(data_list)

# save(data_list, file = "go_data_list.rdata")
