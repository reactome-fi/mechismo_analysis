library(magrittr)

read.delim("interfacewise_de.tsv") %>%
  dplyr::mutate(DE.Gene.Wilcox.BH.Adj.p = p.adjust(DE.Gene.Wilcox.p,method="BH")) %>%
  write.table("interfacewise_de.tsv",
    row.names = FALSE,
            sep="\t")
