library(magrittr)

LUSTRE_DIR <- "/home/users/burkhajo/results"

read.delim(paste(LUSTRE_DIR,
                 "/interfacewise_de.tsv",sep="")) %>%
  dplyr::group_by(Cancer.Type) %>%
  dplyr::mutate(DE.Gene.Wilcox.BH.Adj.p = p.adjust(DE.Gene.Wilcox.p,method="BH")) %>%
  dplyr::select(Cancer.Type,
                Interface,
                DE.Gene,
                DE.Gene.Wilcox.p,
                DE.Gene.Wilcox.BH.Adj.p,
                Num.Interface.Samples,
                Num.NoInterface.Samples,
                Interface.NoInterface.Diff) %>%
  write.table(paste(LUSTRE_DIR,"/interfacewise_de.tsv",sep=""),
    row.names = FALSE,
            sep="\t")
