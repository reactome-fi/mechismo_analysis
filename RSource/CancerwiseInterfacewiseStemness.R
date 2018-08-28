#libs
library(dplyr)
library(magrittr)
library(ggplot2)
library(stringr)

#globs
MECH_INTERFACES <- "/Users/joshuaburkhart/Downloads/tcga_mechismo_stat_cancer_wise_undirected_significant.tsv"
DNA_METH_STEMNS <- "/Users/joshuaburkhart/Downloads/StemnessScores_DNAmeth.csv"
RNA_EXPR_STEMNS <- "/Users/joshuaburkhart/Downloads/StemnessScores_RNAexp.csv"
REACTOME_FIS <- "/Users/joshuaburkhart/Downloads/ProteinFIsInReactions_073118.txt"
CANCER_CENSUS <- "/Users/joshuaburkhart/Downloads/Census_allWed Aug 22 22_25_41 2018.tsv"
MECH_INPUT <- "/Users/joshuaburkhart/Downloads/TCGA_mech_input.tsv"

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
mech_input_df <- read.delim(MECH_INPUT,
                            stringsAsFactors = FALSE)

#wrangle
interface_cancer_types <- mech_interfaces_df[,1] %>%
  unique()

mech_interfaces_df2 <- mech_interfaces_df %>%
  dplyr::arrange(desc(V8)) %>%
  dplyr::filter(V8 >= 5)

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

cancer_census_df <- cancer_census_df %>%
  dplyr::mutate(Hallmark = ifelse(Hallmark == "Yes",TRUE,FALSE))

gene_samples_hash <- new.env()
for(i in 1:nrow(mech_input_df)){
  gene <- mech_input_df[i,"Gene.symbol"]
  samples_string <- mech_input_df[i,"Sample.ID"]
  samples_by_cancer_type_strings <- stringr::str_split(samples_string,
                                                           ";") %>%
    unlist()
  cancer_types_mtx <- stringr::str_match_all(samples_string,
                                           "([A-Z]+):")
  cancer_types <- cancer_types_mtx[[1]][,2] %>%
    unique()
  
  for(j in 1:length(cancer_types)){
    cancer_type <- cancer_types[j]
    cur_samples <- gene_samples_hash[[paste(gene,cancer_type,sep="-")]]
    if(is.null(cur_samples)){
      cur_samples <- c()
    }
    new_samples <- c()
    for(k in 1:length(samples_by_cancer_type_strings)){
      samples_by_cancer_type_string <- samples_by_cancer_type_strings[k]
      if(grepl(cancer_type,samples_by_cancer_type_string)){
        new_samples <- stringr::str_extract_all(samples_by_cancer_type_string,
              "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %>%
              unlist()
        samples <- union(cur_samples,new_samples)
      }
      gene_samples_hash[[paste(gene,cancer_type,sep="-")]] <- samples
    }
  }
}

#results df

results_df3 <- data.frame(Cancer.Type = character(),
                          Interface = character(),
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
                          interface_samples = character(),
                          mDNAsi.t.test.BH = numeric(),
                          mDNAsi.wilcox.BH = numeric(),
                          EREG.mDNAsi.t.test.BH = numeric(),
                          EREG.mDNAsi.wilcox.BH = numeric(),
                          DMPsi.t.test.BH = numeric(),
                          DMPsi.wilcox.BH = numeric(),
                          ENHsi.t.test.BH = numeric(),
                          ENHsi.wilcox.BH = numeric(),
                          mRNAsi.t.test.BH = numeric(),
                          mRNAsi.wilcox.BH = numeric(),
                          EREG.mRNAsi.t.test.BH = numeric(),
                          EREG.mRNAsi.wilcox.BH = numeric())

#calculate results
for(j in 1:length(interface_cancer_types)){
  print(paste(j," of ",length(interface_cancer_types),sep=""))
  cancer_type <- interface_cancer_types[j]
  mech_interfaces_df3 <- mech_interfaces_df2 %>%
    dplyr::filter(V1 == cancer_type)
  
  #results df
  results_df <- data.frame(Cancer.Type = character(),
                           Interface = character(),
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
  
  for(i in 1:nrow(mech_interfaces_df3)){
  gene1 <- mech_interfaces_df3[i,2]
  gene2 <- mech_interfaces_df3[i,3]
    interface <- paste(gene1,
                       "-",
                       gene2,
                       sep="")
    
    interface_samples <- stringr::str_extract_all(mech_interfaces_df3[i,16],
                                                  "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %>%
      unlist()
    gene_samples <- union(gene_samples_hash[[paste(gene1,cancer_type,sep="-")]],
                          gene_samples_hash[[paste(gene2,cancer_type,sep="-")]])
  
    gene_no_interface_samples <- setdiff(gene_samples,interface_samples)
  
    #mDNAsi distributions
    mDNAsi_mut <- stemness_df %>%
      dplyr::filter(TCGA.sample %in% interface_samples,
                    !is.na(mDNAsi)) %>%
      dplyr::select(mDNAsi)
    mDNAsi_wt <- stemness_df %>%
      dplyr::filter(TCGA.sample %in% gene_no_interface_samples,
                    !is.na(mDNAsi)) %>%
      dplyr::select(mDNAsi)
    
    #EREG.mDNAsi distributions
    EREG.mDNAsi_mut <- stemness_df %>%
      dplyr::filter(TCGA.sample %in% interface_samples,
                    !is.na(EREG.mDNAsi)) %>%
      dplyr::select(EREG.mDNAsi)
    EREG.mDNAsi_wt <- stemness_df %>%
      dplyr::filter(TCGA.sample %in% gene_no_interface_samples,
                    !is.na(EREG.mDNAsi)) %>%
      dplyr::select(EREG.mDNAsi)
    
    #DMPsi distributions
    DMPsi_mut <- stemness_df %>%
      dplyr::filter(TCGA.sample %in% interface_samples,
                    !is.na(DMPsi)) %>%
      dplyr::select(DMPsi)
    DMPsi_wt <- stemness_df %>%
      dplyr::filter(TCGA.sample %in% gene_no_interface_samples,
                    !is.na(DMPsi)) %>%
      dplyr::select(DMPsi)
    
    #ENHsi distributions
    ENHsi_mut <- stemness_df %>%
      dplyr::filter(TCGA.sample %in% interface_samples,
                    !is.na(ENHsi)) %>%
      dplyr::select(ENHsi)
    ENHsi_wt <- stemness_df %>%
      dplyr::filter(TCGA.sample %in% gene_no_interface_samples,
                    !is.na(ENHsi)) %>%
      dplyr::select(ENHsi)
    
    #mRNAsi distributions
    mRNAsi_mut <- stemness_df %>%
      dplyr::filter(TCGA.sample %in% interface_samples,
                    !is.na(mRNAsi)) %>%
      dplyr::select(mRNAsi)
    mRNAsi_wt <- stemness_df %>%
      dplyr::filter(TCGA.sample %in% gene_no_interface_samples,
                    !is.na(mRNAsi)) %>%
      dplyr::select(mRNAsi)
    
    #EREG.mRNAsi distributions
    EREG.mRNAsi_mut <- stemness_df %>%
      dplyr::filter(TCGA.sample %in% interface_samples,
                    !is.na(EREG.mRNAsi)) %>%
      dplyr::select(EREG.mRNAsi)
    EREG.mRNAsi_wt <- stemness_df %>%
      dplyr::filter(TCGA.sample %in% gene_no_interface_samples,
                    !is.na(EREG.mRNAsi)) %>%
      dplyr::select(EREG.mRNAsi)
    
    tryCatch(
      {
        result_df <- data.frame(Cancer.Type = cancer_type,
                                Interface = interface,
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
      },
      error=function(cond){
        #do nothing
      }
    )
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
  
  results_df3 <- results_df3 %>%
    rbind(results_df2) 
}

results_df3 <- results_df3 %>%
  dplyr::arrange(In.Reactome)

results_df3 %>%
  dplyr::select(Cancer.Type,
                Interface,
                In.Reactome,
                G1.In.Cancer.Census,
                G2.In.Cancer.Census,
                G1.Tier,
                G2.Tier,
                G1.Hallmark,
                G2.Hallmark,
                num_interface_samples,
                mDNAsi.wilcox,
                mDNAsi.wilcox.BH,
                mRNAsi.wilcox,
                mRNAsi.wilcox.BH,
                interface_samples) %>%
  write.table("CancerwiseInterfacewiseStemness.tsv",
            row.names = FALSE,
            sep = "\t")


results_df3 %>%
  ggplot2::ggplot(aes(-log10(mRNAsi.wilcox),fill=In.Reactome,color=In.Reactome)) +
  ggplot2::geom_density(alpha=0.1) +
  ggplot2::scale_color_brewer(palette="Set1")

ggsave("Cancerwise_mRNA_wilcox_density.png",width=10,height=10,dpi=600)

results_df3 %>%
  ggplot2::ggplot(aes(-log10(mDNAsi.wilcox),fill=In.Reactome,color=In.Reactome)) +
  ggplot2::geom_density(alpha=0.1) +
  ggplot2::scale_color_brewer(palette="Set1")

ggsave("Cancerwise_mDNAsi_wilcox_density.png",width=10,height=10,dpi=600)

results_df3 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.wilcox.BH),y=-log10(mDNAsi.wilcox.BH),color=In.Reactome,label=Interface)) +
  ggplot2::geom_point() +
  scale_color_brewer(palette="Set1")

ggsave("Cancerwise_interface_stemness_significance.png",width=10,height=10,dpi=600)
  
results_df3 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.wilcox.BH),y=-log10(mDNAsi.wilcox.BH),color=In.Reactome,label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")

ggsave("Cancerwise_interface_stemness_significance_labeled.png",width=10,height=10,dpi=600)

results_df3 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.wilcox.BH),y=-log10(mDNAsi.wilcox.BH),color=In.Reactome,label=Interface)) +
  ggplot2::geom_point() +
  scale_color_brewer(palette="Set1") +
  xlim(0,5) +
  ylim(0,5)

ggsave("Cancerwise_interface_stemness_significance_upr_bounded.png",width=10,height=10,dpi=600)
  
results_df3 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.wilcox.BH),y=-log10(mDNAsi.wilcox.BH),color=In.Reactome,label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")+
  xlim(0,5)+
  ylim(0,5)

ggsave("Cancerwise_interface_stemness_significance_upr_bounded_labeled.png",width=10,height=10,dpi=600)

results_df3 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.wilcox.BH),y=-log10(mDNAsi.wilcox.BH),color=In.Reactome,label=Interface)) +
  ggplot2::geom_point() +
  scale_color_brewer(palette="Set1") +
  xlim(2,NA) +
  ylim(2,NA)

ggsave("Cancerwise_interface_stemness_significance_lwr_bounded.png",width=10,height=10,dpi=600)
  
results_df3 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.wilcox.BH),y=-log10(mDNAsi.wilcox.BH),color=In.Reactome,label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")+
  xlim(2,NA)+
  ylim(2,NA)

ggsave("Cancerwise_interface_stemness_significance_lwr_bounded_labeled.png",width=10,height=10,dpi=600)

results_df3 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.wilcox.BH),y=-log10(mDNAsi.wilcox.BH),color=In.Reactome,label=Interface)) +
  ggplot2::geom_point() +
  scale_color_brewer(palette="Set1") +
  xlim(-log10(0.05),NA)

ggsave("Cancerwise_interface_stemness_significance_mRNAsi_lwr_bounded.png",width=10,height=10,dpi=600)
  
results_df3 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.wilcox.BH),y=-log10(mDNAsi.wilcox.BH),color=In.Reactome,label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")+
  xlim(-log10(0.05),NA)

ggsave("Cancerwise_interface_stemness_significance_mRNAsi_lwr_bounded_labeled.png",width=10,height=10,dpi=600)

results_df3 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.wilcox.BH),y=-log10(mDNAsi.wilcox.BH),color=In.Reactome,label=Interface)) +
  ggplot2::geom_point() +
  scale_color_brewer(palette="Set1") +
  ylim(-log10(0.05),NA)

ggsave("Cancerwise_interface_stemness_significance_mDNAsi_lwr_bounded.png",width=10,height=10,dpi=600)
  
results_df3 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.wilcox.BH),y=-log10(mDNAsi.wilcox.BH),color=In.Reactome,label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")+
  ylim(-log10(0.05),NA)

ggsave("Cancerwise_interface_stemness_significance_mDNAsi_lwr_bounded_labeled.png",width=10,height=10,dpi=600)

coi <- "LGG"

results_df3 %>%
  dplyr::filter(grepl(coi,Cancer.Type)) %>%
  dplyr::arrange(In.Reactome) %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.wilcox.BH),
                      y=-log10(mDNAsi.wilcox.BH),
                      color=In.Reactome,label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=.5,vjust=-.5)+
  scale_color_brewer(palette="Set1")

ggsave(paste(coi,"_interface_stemness_significance_labeled.png",sep=""),width=10,height=10,dpi=600)
