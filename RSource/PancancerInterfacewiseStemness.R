#libs
library(dplyr)
library(magrittr)
library(stringr)

#globs
MECH_INTERFACES <- "/Users/joshuaburkhart/Downloads/tcga_mechismo_stat_pancancer_significant.tsv"
DNA_METH_STEMNS <- "/Users/joshuaburkhart/Downloads/StemnessScores_DNAmeth.csv"
RNA_EXPR_STEMNS <- "/Users/joshuaburkhart/Downloads/StemnessScores_RNAexp.csv"
REACTOME_FIS <- "/Users/joshuaburkhart/Downloads/FIsInGene_071718_with_annotations.txt"
CANCER_CENSUS <- "/Users/joshuaburkhart/Downloads/Census_allWed Aug 22 22_25_41 2018.tsv"

#load data
mech_interfaces_df <- read.delim(MECH_INTERFACES,
                                 stringsAsFactors = FALSE,
                                 header = FALSE)
dna_meth_stemns_df <- read.delim(DNA_METH_STEMNS,
                                 stringsAsFactors = FALSE,
                                 sep = ',')
rna_expr_stemns_df <- read.delim(RNA_EXPR_STEMNS,
                                 stringsAsFactors = FALSE,
                                 sep = ',')
reactome_fis_df <- read.delim(REACTOME_FIS,
                              stringsAsFactors = FALSE)
cancer_census_df <- read.delim(CANCER_CENSUS,
                               stringsAsFactors = FALSE)

#wrangle
mech_interfaces_df2 <- mech_interfaces_df %>%
  dplyr::arrange(desc(V7)) %>%
  dplyr::filter(V7 >= 10)

dna_meth_stemns_df2 <- dna_meth_stemns_df %>%
  dplyr::mutate(
    TCGA.sample = stringr::str_extract(
      TCGAlong.id,
      "^TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}")) %>%
  dplyr::select(-TCGAlong.id)

rna_expr_stemns_df2 <- rna_expr_stemns_df %>%
  dplyr::mutate(
    TCGA.sample = stringr::str_extract(
      TCGAlong.id,
      "^TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}")) %>%
  dplyr::select(-TCGAlong.id)

stemness_df <- dna_meth_stemns_df2 %>%
  dplyr::full_join(rna_expr_stemns_df2)

reactome_fis_df2 <- reactome_fis_df %>%
  dplyr::mutate(fwd_fi = paste(Gene1,"-",Gene2,sep=""),
                rev_fi = paste(Gene2,"-",Gene1,sep="")) %>%
  dplyr::select(fwd_fi,rev_fi)

#results df
results_df <- data.frame(Interface = character(),
                         In.Reactome = logical(),
                         G1.In.Cancer.Census = logical(),
                         G2.In.Cancer.Census = logical(),
                         G1.Tier = numeric(),
                         G2.Tier = numeric(),
                         G1.Hallmark = logical(),
                         G2.Hallmark = logical(),
                         num_interface_samples = integer(),
                         mDNAsi.t.test = numeric(),
                         mDNAsi.wilcox = numeric(),
                         EREG.mDNAsi.t.test = numeric(),
                         EREG.mDNAsi.wilcox = numeric(),
                         DMPsi.t.test = numeric(),
                         DMPsi.wilcox = numeric(),
                         ENHsi.t.test = numeric(),
                         ENHsi.wilcox = numeric(),
                         mRNAsi.t.test = numeric(),
                         mRNAsi.wilcox = numeric(),
                         EREG.mRNAsi.t.test = numeric(),
                         EREG.mRNAsi.wilcox = numeric(),
                         interface_samples = character())

#calculate results
for(i in 1:nrow(mech_interfaces_df2)){
  gene1 <- mech_interfaces_df2[i,1]
  gene2 <- mech_interfaces_df2[i,2]
  interface <- paste(gene1,
                     "-",
                     gene2,
                     sep="")
  interface_samples <- stringr::str_extract_all(mech_interfaces_df2[i,15],
                                                "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %>%
    unlist()
  #mDNAsi distributions
  mDNAsi_mut <- stemness_df %>%
    dplyr::filter(TCGA.sample %in% interface_samples,
                  !is.na(mDNAsi)) %>%
    dplyr::select(mDNAsi)
  mDNAsi_wt <- stemness_df %>%
    dplyr::filter(!TCGA.sample %in% interface_samples,
                  !is.na(mDNAsi)) %>%
    dplyr::select(mDNAsi)
  
  #EREG.mDNAsi distributions
  EREG.mDNAsi_mut <- stemness_df %>%
    dplyr::filter(TCGA.sample %in% interface_samples,
                  !is.na(EREG.mDNAsi)) %>%
    dplyr::select(EREG.mDNAsi)
  EREG.mDNAsi_wt <- stemness_df %>%
    dplyr::filter(!TCGA.sample %in% interface_samples,
                  !is.na(EREG.mDNAsi)) %>%
    dplyr::select(EREG.mDNAsi)
  
  #DMPsi distributions
  DMPsi_mut <- stemness_df %>%
    dplyr::filter(TCGA.sample %in% interface_samples,
                  !is.na(DMPsi)) %>%
    dplyr::select(DMPsi)
  DMPsi_wt <- stemness_df %>%
    dplyr::filter(!TCGA.sample %in% interface_samples,
                  !is.na(DMPsi)) %>%
    dplyr::select(DMPsi)
  
  #ENHsi distributions
  ENHsi_mut <- stemness_df %>%
    dplyr::filter(TCGA.sample %in% interface_samples,
                  !is.na(ENHsi)) %>%
    dplyr::select(ENHsi)
  ENHsi_wt <- stemness_df %>%
    dplyr::filter(!TCGA.sample %in% interface_samples,
                  !is.na(ENHsi)) %>%
    dplyr::select(ENHsi)
  
  #mRNAsi distributions
  mRNAsi_mut <- stemness_df %>%
    dplyr::filter(TCGA.sample %in% interface_samples,
                  !is.na(mRNAsi)) %>%
    dplyr::select(mRNAsi)
  mRNAsi_wt <- stemness_df %>%
    dplyr::filter(!TCGA.sample %in% interface_samples,
                  !is.na(mRNAsi)) %>%
    dplyr::select(mRNAsi)
  
  #EREG.mRNAsi distributions
  EREG.mRNAsi_mut <- stemness_df %>%
    dplyr::filter(TCGA.sample %in% interface_samples,
                  !is.na(EREG.mRNAsi)) %>%
    dplyr::select(EREG.mRNAsi)
  EREG.mRNAsi_wt <- stemness_df %>%
    dplyr::filter(!TCGA.sample %in% interface_samples,
                  !is.na(EREG.mRNAsi)) %>%
    dplyr::select(EREG.mRNAsi)
  
  result_df <- data.frame(Interface = interface,
                          In.Reactome = (interface %in% reactome_fis_df2$fwd_fi |
                                           interface %in% reactome_fis_df2$rev_fi),
                          G1.In.Cancer.Census = gene1 %in% cancer_census_df$Gene.Symbol,
                          G2.In.Cancer.Census = gene2 %in% cancer_census_df$Gene.Symbol,
                          G1.Tier = ifelse(gene1 %in% cancer_census_df$Gene.Symbol,
                                           cancer_census_df %>%
                                             dplyr::filter(Gene.Symbol == gene1) %>%
                                             dplyr::select(Tier) %>%
                                             as.numeric(),
                                           0),
                          G2.Tier = ifelse(gene2 %in% cancer_census_df$Gene.Symbol,
                                           cancer_census_df %>%
                                             dplyr::filter(Gene.Symbol == gene2) %>%
                                             dplyr::select(Tier) %>%
                                             as.numeric(),
                                           0),
                          G1.Hallmark = ifelse(gene1 %in% cancer_census_df$Gene.Symbol,
                                               cancer_census_df %>%
                                                 dplyr::filter(Gene.Symbol == gene1) %>%
                                                 dplyr::select(Hallmark) %>%
                                                 as.logical(),
                                               FALSE),
                          G2.Hallmark = ifelse(gene2 %in% cancer_census_df$Gene.Symbol,
                                               cancer_census_df %>%
                                                 dplyr::filter(Gene.Symbol == gene2) %>%
                                                 dplyr::select(Hallmark) %>%
                                                 as.logical(),
                                               FALSE),
                          num_interface_samples = length(interface_samples),
                          mDNAsi.t.test = t.test(mDNAsi_mut$mDNAsi,
                                                 mDNAsi_wt$mDNAsi)$p.value,
                          mDNAsi.wilcox = wilcox.test(mDNAsi_mut$mDNAsi,
                                                      mDNAsi_wt$mDNAsi)$p.value,
                          EREG.mDNAsi.t.test = t.test(EREG.mDNAsi_mut$EREG.mDNAsi,
                                                      EREG.mDNAsi_wt$EREG.mDNAsi)$p.value,
                          EREG.mDNAsi.wilcox = wilcox.test(EREG.mDNAsi_mut$EREG.mDNAsi,
                                                           EREG.mDNAsi_wt$EREG.mDNAsi)$p.value,
                          DMPsi.t.test = t.test(DMPsi_mut$DMPsi,
                                                DMPsi_wt$DMPsi)$p.value,
                          DMPsi.wilcox = wilcox.test(DMPsi_mut$DMPsi,
                                                     DMPsi_wt$DMPsi)$p.value,
                          ENHsi.t.test = t.test(ENHsi_mut$ENHsi,
                                                ENHsi_wt$ENHsi)$p.value,
                          ENHsi.wilcox = wilcox.test(ENHsi_mut$ENHsi,
                                                     ENHsi_wt$ENHsi)$p.value,
                          mRNAsi.t.test = t.test(mRNAsi_mut$mRNAsi,
                                                 mRNAsi_wt$mRNAsi)$p.value,
                          mRNAsi.wilcox = wilcox.test(mRNAsi_mut$mRNAsi,
                                                      mRNAsi_wt$mRNAsi)$p.value,
                          EREG.mRNAsi.t.test = t.test(EREG.mRNAsi_mut$EREG.mRNAsi,
                                                      EREG.mRNAsi_wt$EREG.mRNAsi)$p.value,
                          EREG.mRNAsi.wilcox = wilcox.test(EREG.mRNAsi_mut$EREG.mRNAsi,
                                                           EREG.mRNAsi_wt$EREG.mRNAsi)$p.value,
                          interface_samples = paste(interface_samples,
                                                    collapse = ','))
  results_df <- results_df %>%
    rbind(result_df)
}

#adjust p-values
results_df2 <- results_df %>%
  dplyr::mutate(mDNAsi.t.test.BH = p.adjust(mDNAsi.t.test,method="BH"),
                mDNAsi.wilcox.BH = p.adjust(mDNAsi.wilcox,method="BH"),
                EREG.mDNAsi.t.test.BH = p.adjust(EREG.mDNAsi.t.test,method="BH"),
                EREG.mDNAsi.wilcox.BH = p.adjust(EREG.mDNAsi.wilcox,method="BH"),
                DMPsi.t.test.BH = p.adjust(DMPsi.t.test,method="BH"),
                DMPsi.wilcox.BH = p.adjust(DMPsi.wilcox,method="BH"),
                ENHsi.t.test.BH = p.adjust(ENHsi.t.test,method="BH"),
                ENHsi.wilcox.BH = p.adjust(ENHsi.wilcox,method="BH"),
                mRNAsi.t.test.BH = p.adjust(mRNAsi.t.test,method="BH"),
                mRNAsi.wilcox.BH = p.adjust(mRNAsi.wilcox,method="BH"),
                EREG.mRNAsi.t.test.BH = p.adjust(EREG.mRNAsi.t.test,method="BH"),
                EREG.mRNAsi.wilcox.BH = p.adjust(EREG.mRNAsi.wilcox,method="BH"))

write.table(results_df2,
            "PancancerInterfacewiseStemness.tsv",
            row.names = FALSE,
            sep = "\t")
