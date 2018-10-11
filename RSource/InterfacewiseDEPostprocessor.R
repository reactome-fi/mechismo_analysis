library(magrittr)

read.delim("interfacewise_de.tsv") %>%
  dplyr::mutate(DE.Gene.Wilcox.BH.Adj.p = p.adjust(DE.Gene.Wilcox.p,method="BH")) %>%
  dplyr::select(Interface,
                DE.Gene,
                DE.Gene.Wilcox.p,
                DE.Gene.Wilcox.BH.Adj.p,
                Num.Interface.Samples,
                Num.NoInterface.Samples,
                Interface.NoInterface.Diff) %>%
  write.table("interfacewise_de.tsv",
    row.names = FALSE,
            sep="\t")
