library(magrittr)

LUSTRE_DIR <- "/home/exacloud/lustre1/WongLab/tmp"

read.delim(paste(LUSTRE_DIR,
                 "/pinterfacewise_de.tsv",sep="")) %>%
  dplyr::mutate(DE.Gene.Wilcox.BH.Adj.p = p.adjust(DE.Gene.Wilcox.p,method="BH")) %>%
  dplyr::select(Interface,
                DE.Gene,
                DE.Gene.Wilcox.p,
                DE.Gene.Wilcox.BH.Adj.p,
                Num.Interface.Samples,
                Num.NoInterface.Samples,
                Interface.NoInterface.Diff) %>%
  write.table(paste(LUSTRE_DIR,"/pinterfacewise_de.tsv",sep=""),
    row.names = FALSE,
            sep="\t")
