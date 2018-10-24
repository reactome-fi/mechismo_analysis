
library(dplyr)
z <- read.delim(file="/Users/joshuaburkhart/interfacewise_de.tsv")
types <- as.character(z$Cancer.Type %>% unique())

for(i in 1:length(types)){
  type <- types[i]
  z %>%
    dplyr::filter(as.character(Cancer.Type) == type) %>%
    dplyr::arrange(DE.Gene.Wilcox.BH.Adj.p,
                   Num.Interface.Samples) %>%
    write.table(file=paste(type,"_interfacewise_de.tsv"))
}
