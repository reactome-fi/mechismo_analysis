#libs
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(stringr)

#globs
MECH_INTERFACES <- "/Users/joshuaburkhart/Downloads/tcga_mechismo_stat_pancancer_undirected_significant.tsv"
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
mech_interfaces_df2 <- mech_interfaces_df %>%
  dplyr::arrange(desc(V7))

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
  cur_samples <- gene_samples_hash[[gene]]
  if(is.null(cur_samples)){
    cur_samples <- c()
  }
  new_samples_string <- mech_input_df[i,"Sample.ID"]
  new_samples <- stringr::str_extract_all(new_samples_string,
                                          "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %>%
    unlist()
  samples <- union(cur_samples,new_samples)
  gene_samples_hash[[gene]] <- samples
}

save(gene_samples_hash,file="pancancer_gene_samples_hash.rda")

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
                         mDNAsi_interface_v_gene_no_interface_med_diff = numeric(),
                         mDNAsi.interface_v_all.wilcox = numeric(),
                         mDNAsi.interface_v_gene_no_interface.wilcox = numeric(),
                         mDNAsi.gene_v_all.wilcox = numeric(),
                         mDNAsi.gene_no_interface_v_all.wilcox = numeric(),
                         mRNAsi_interface_v_gene_no_interface_med_diff = numeric(),
                         mRNAsi.interface_v_all.wilcox = numeric(),
                         mRNAsi.interface_v_gene_no_interface.wilcox = numeric(),
                         mRNAsi.gene_v_all.wilcox = numeric(),
                         mRNAsi.gene_no_interface_v_all.wilcox = numeric(),
                         interface_samples = character())

#calculate results
for(i in 1:nrow(mech_interfaces_df2)){
  print(paste(i," of ",nrow(mech_interfaces_df2),sep=""))
  gene1 <- mech_interfaces_df2[i,1]
  gene2 <- mech_interfaces_df2[i,2]
  interface <- paste(gene1,
                     "-",
                     gene2,
                     sep="")
  interface_samples <- stringr::str_extract_all(mech_interfaces_df2[i,15],
                                                "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %>%
    unlist()
  
  gene_samples <- union(gene_samples_hash[[gene1]],gene_samples_hash[[gene2]])
  
  gene_no_interface_samples <- setdiff(gene_samples,interface_samples)
  
  #mDNAsi distributions
  mDNAsi_interface <- stemness_df %>%
    dplyr::filter(TCGA.sample %in% interface_samples,
                  !is.na(mDNAsi)) %>%
    dplyr::select(mDNAsi)
  mDNAsi_gene_no_interface <- stemness_df %>%
    dplyr::filter(TCGA.sample %in% gene_no_interface_samples,
                  !is.na(mDNAsi)) %>%
    dplyr::select(mDNAsi)
  mDNAsi_gene <- stemness_df %>%
    dplyr::filter(TCGA.sample %in% gene_samples,
                  !is.na(mDNAsi)) %>%
    dplyr::select(mDNAsi)
  mDNAsi_no_interface <- stemness_df %>%
    dplyr::filter(!(TCGA.sample %in% interface_samples),
                  !is.na(mDNAsi)) %>%
    dplyr::select(mDNAsi)
  mDNAsi_no_gene <- stemness_df %>%
    dplyr::filter(!(TCGA.sample %in% gene_samples),
                  !is.na(mDNAsi)) %>%
    dplyr::select(mDNAsi)
  
  #mRNAsi distributions
  mRNAsi_interface <- stemness_df %>%
    dplyr::filter(TCGA.sample %in% interface_samples,
                  !is.na(mRNAsi)) %>%
    dplyr::select(mRNAsi)
  mRNAsi_gene_no_interface <- stemness_df %>%
    dplyr::filter(TCGA.sample %in% gene_no_interface_samples,
                  !is.na(mRNAsi)) %>%
    dplyr::select(mRNAsi)
  mRNAsi_gene <- stemness_df %>%
    dplyr::filter(TCGA.sample %in% gene_samples,
                  !is.na(mRNAsi)) %>%
    dplyr::select(mRNAsi)
  mRNAsi_no_interface <- stemness_df %>%
    dplyr::filter(!(TCGA.sample %in% interface_samples),
                  !is.na(mRNAsi)) %>%
    dplyr::select(mRNAsi)
  mRNAsi_no_gene <- stemness_df %>%
    dplyr::filter(!(TCGA.sample %in% gene_samples),
                  !is.na(mRNAsi)) %>%
    dplyr::select(mRNAsi)
  
  #  top_stemness_interfaces <- c(
  #    "BRAF-MAP2K1",
  #    "PLCB3-GNAQ",
  #    "GNAQ-RGS2",
  #    "PLCB1-GNAQ",
  #    "GNA11-PLCB3",
  #    "GNAQ-RGS21",
  #    "GNAQ-RGS3",
  #    "GNB1-GNAQ",
  #    "GNAQ-PLCB2",
  #    "RGS18-GNAQ",
  #    "GNB2-GNAQ")
  #  
  #  if(interface %in% top_stemness_interfaces){
  if((interface %in% reactome_fis_df2$fwd_fi |
      interface %in% reactome_fis_df2$rev_fi) &
     grepl("GNAQ",interface)){  
    mDNAsi_distribution_df0 <- mDNAsi_interface %>%
      as.data.frame() %>%
      dplyr::mutate(mutation = paste("Within Interface (", nrow(mDNAsi_interface),")",sep=""))
    
    mDNAsi_distribution_df1 <- mDNAsi_gene_no_interface %>%
      as.data.frame() %>%
      dplyr::mutate(mutation = paste("Outside Interface (", nrow(mDNAsi_gene_no_interface),")",sep=""))
    
    mDNAsi_distribution_df2 <- mDNAsi_no_gene %>%
      as.data.frame() %>%
      dplyr::mutate(mutation = paste("No Mutation (", nrow(mDNAsi_no_gene), ")", sep=""))
    
    mDNAsi_distribution_df <- dplyr::bind_rows(mDNAsi_distribution_df0,
                                               mDNAsi_distribution_df1,
                                               mDNAsi_distribution_df2)
    
    ggpubr::compare_means(mDNAsi ~ mutation,data=mDNAsi_distribution_df)
    
    comparisons <- list(c(mDNAsi_distribution_df$mutation %>% unique() %>% .[1],
                          mDNAsi_distribution_df$mutation %>% unique() %>% .[2]),
                        c(mDNAsi_distribution_df$mutation %>% unique() %>% .[2],
                          mDNAsi_distribution_df$mutation %>% unique() %>% .[3]),
                        c(mDNAsi_distribution_df$mutation %>% unique() %>% .[1],
                          mDNAsi_distribution_df$mutation %>% unique() %>% .[3]))
    
    mDNAsi_distribution_df %>%
      ggplot(aes(x=mutation,y=mDNAsi,fill=mutation)) +
      geom_boxplot() +
      ggpubr::stat_compare_means(comparisons = comparisons) +
      xlab("Mutation") +
      ylab("mDNAsi") +
      ggtitle(interface)
    
    ggsave(paste(interface,"_","mDNAsi.png",sep=""),width=10,height=10,dpi=600)
    
    mRNAsi_distribution_df0 <- mRNAsi_interface %>%
      as.data.frame() %>%
      dplyr::mutate(mutation = paste("Within Interface (", nrow(mRNAsi_interface),")",sep=""))
    
    mRNAsi_distribution_df1 <- mRNAsi_gene_no_interface %>%
      as.data.frame() %>%
      dplyr::mutate(mutation = paste("Outside Interface (", nrow(mRNAsi_gene_no_interface),")",sep=""))
    
    mRNAsi_distribution_df2 <- mRNAsi_no_gene %>%
      as.data.frame() %>%
      dplyr::mutate(mutation = paste("No Mutation (", nrow(mRNAsi_no_gene), ")", sep=""))
    
    mRNAsi_distribution_df <- dplyr::bind_rows(mRNAsi_distribution_df0,
                                               mRNAsi_distribution_df1,
                                               mRNAsi_distribution_df2)
    
    ggpubr::compare_means(mRNAsi ~ mutation,data=mRNAsi_distribution_df)
    
    comparisons <- list(c(mRNAsi_distribution_df$mutation %>% unique() %>% .[1],
                          mRNAsi_distribution_df$mutation %>% unique() %>% .[2]),
                        c(mRNAsi_distribution_df$mutation %>% unique() %>% .[2],
                          mRNAsi_distribution_df$mutation %>% unique() %>% .[3]),
                        c(mRNAsi_distribution_df$mutation %>% unique() %>% .[1],
                          mRNAsi_distribution_df$mutation %>% unique() %>% .[3]))
    
    mRNAsi_distribution_df %>%
      ggplot(aes(x=mutation,y=mRNAsi,fill=mutation)) +
      geom_boxplot() +
      ggpubr::stat_compare_means(comparisons = comparisons) +
      xlab("Mutation") +
      ylab("mRNAsi") +
      ggtitle(interface)
    
    ggsave(paste(interface,"_","mRNAsi.png",sep=""),width=10,height=10,dpi=600)
  }
  
  tryCatch(
    {
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
                              mDNAsi_interface_v_gene_no_interface_med_diff =
                                median(mDNAsi_interface$mDNAsi) - median(mDNAsi_gene_no_interface$mDNAsi),
                              mDNAsi.interface_v_all.wilcox =
                                wilcox.test(mDNAsi_interface$mDNAsi,
                                            mDNAsi_no_interface$mDNAsi)$p.value,
                              mDNAsi.interface_v_gene_no_interface.wilcox =
                                wilcox.test(mDNAsi_interface$mDNAsi,
                                            mDNAsi_gene_no_interface$mDNAsi)$p.value,
                              mDNAsi.gene_v_all.wilcox =
                                wilcox.test(mDNAsi_gene$mDNAsi,
                                            mDNAsi_no_gene$mDNAsi)$p.value,
                              mDNAsi.gene_no_interface_v_all.wilcox =
                                wilcox.test(mDNAsi_gene_no_interface$mDNAsi,
                                            mDNAsi_no_gene$mDNAsi)$p.value,
                              mRNAsi_interface_v_gene_no_interface_med_diff =
                                median(mRNAsi_interface$mRNAsi) - median(mRNAsi_gene_no_interface$mRNAsi),
                              mRNAsi.interface_v_all.wilcox =
                                wilcox.test(mRNAsi_interface$mRNAsi,
                                            mRNAsi_no_interface$mRNAsi)$p.value,
                              mRNAsi.interface_v_gene_no_interface.wilcox =
                                wilcox.test(mRNAsi_interface$mRNAsi,
                                            mRNAsi_gene_no_interface$mRNAsi)$p.value,
                              mRNAsi.gene_v_all.wilcox =
                                wilcox.test(mRNAsi_gene$mRNAsi,
                                            mRNAsi_no_gene$mRNAsi)$p.value,
                              mRNAsi.gene_no_interface_v_all.wilcox =
                                wilcox.test(mRNAsi_gene_no_interface$mRNAsi,
                                            mRNAsi_no_gene$mRNAsi)$p.value,
                              interface_samples =
                                paste(interface_samples,
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
  dplyr::mutate(mDNAsi.interface_v_all.wilcox.BH = p.adjust(mDNAsi.interface_v_all.wilcox,method="BH"),
                mDNAsi.interface_v_gene_no_interface.wilcox.BH = p.adjust(mDNAsi.interface_v_gene_no_interface.wilcox,method="BH"),
                mDNAsi.gene_v_all.wilcox.BH = p.adjust(mDNAsi.gene_v_all.wilcox,method="BH"),
                mDNAsi.gene_no_interface_v_all.wilcox.BH = p.adjust(mDNAsi.gene_no_interface_v_all.wilcox,method="BH"),
                mRNAsi.interface_v_all.wilcox.BH = p.adjust(mRNAsi.interface_v_all.wilcox,method="BH"),
                mRNAsi.interface_v_gene_no_interface.wilcox.BH = p.adjust(mRNAsi.interface_v_gene_no_interface.wilcox,method="BH"),
                mRNAsi.gene_v_all.wilcox.BH = p.adjust(mRNAsi.gene_v_all.wilcox,method="BH"),
                mRNAsi.gene_no_interface_v_all.wilcox.BH = p.adjust(mRNAsi.gene_no_interface_v_all.wilcox,method="BH")) %>%
  dplyr::arrange(In.Reactome)

save(results_df2,file="pancancer_results.rda")

results_df2 %>%
  dplyr::select(Interface,
                In.Reactome,
                G1.In.Cancer.Census,
                G2.In.Cancer.Census,
                G1.Tier,
                G2.Tier,
                G1.Hallmark,
                G2.Hallmark,
                num_interface_samples,
                mDNAsi.interface_v_all.wilcox,
                mDNAsi.interface_v_all.wilcox.BH,
                mRNAsi.interface_v_all.wilcox,
                mRNAsi.interface_v_all.wilcox.BH,
                interface_samples) %>%
  write.table("Pancancer_interface_v_all_Stemness.tsv",
              row.names = FALSE,
              sep = "\t")
results_df2 %>%
  dplyr::select(Interface,
                In.Reactome,
                G1.In.Cancer.Census,
                G2.In.Cancer.Census,
                G1.Tier,
                G2.Tier,
                G1.Hallmark,
                G2.Hallmark,
                num_interface_samples,
                mDNAsi.gene_no_interface_v_all.wilcox,
                mDNAsi.gene_no_interface_v_all.wilcox.BH,
                mRNAsi.gene_no_interface_v_all.wilcox,
                mRNAsi.gene_no_interface_v_all.wilcox.BH,
                interface_samples) %>%
  write.table("Pancancer_gene_no_interface_v_all_Stemness.tsv",
              row.names = FALSE,
              sep = "\t")
results_df2 %>%
  dplyr::select(Interface,
                In.Reactome,
                G1.In.Cancer.Census,
                G2.In.Cancer.Census,
                G1.Tier,
                G2.Tier,
                G1.Hallmark,
                G2.Hallmark,
                num_interface_samples,
                mDNAsi_interface_v_gene_no_interface_med_diff,
                mDNAsi.interface_v_gene_no_interface.wilcox,
                mDNAsi.interface_v_gene_no_interface.wilcox.BH,
                mRNAsi_interface_v_gene_no_interface_med_diff,
                mRNAsi.interface_v_gene_no_interface.wilcox,
                mRNAsi.interface_v_gene_no_interface.wilcox.BH,
                interface_samples) %>%
  write.table("Pancancer_interface_v_gene_no_interface_Stemness.tsv",
              row.names = FALSE,
              sep = "\t")

#interface_v_all
results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.interface_v_all.wilcox.BH),
                      y=-log10(mDNAsi.interface_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")
ggsave("Pancancer_interface_v_all_stemness_significance_labeled.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.interface_v_all.wilcox.BH),
                      y=-log10(mDNAsi.interface_v_all.wilcox.BH),
                      color=In.Reactome)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1") +
  ggplot2::geom_hline(yintercept = -log10(0.05)) +
  ggplot2::geom_vline(xintercept = -log10(0.05))
ggsave("Pancancer_interface_v_all_stemness_significance_hvlines.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.interface_v_all.wilcox.BH),
                      y=-log10(mDNAsi.interface_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")+
  xlim(-log10(0.05),NA)+
  ylim(-log10(0.05),NA)
ggsave("Pancancer_interface_v_all_stemness_significance_lwr_bounded_labeled.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.interface_v_all.wilcox.BH),
                      y=-log10(mDNAsi.interface_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point() +
  scale_color_brewer(palette="Set1") +
  xlim(-log10(0.05),NA)
ggsave("Pancancer_interface_v_all_stemness_significance_mRNAsi_lwr_bounded.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.interface_v_all.wilcox.BH),
                      y=-log10(mDNAsi.interface_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")+
  xlim(-log10(0.05),NA)
ggsave("Pancancer_interface_v_all_stemness_significance_mRNAsi_lwr_bounded_labeled.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.interface_v_all.wilcox.BH),
                      y=-log10(mDNAsi.interface_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point() +
  scale_color_brewer(palette="Set1") +
  ylim(-log10(0.05),NA)
ggsave("Pancancer_interface_v_all_stemness_significance_mDNAsi_lwr_bounded.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.interface_v_all.wilcox.BH),
                      y=-log10(mDNAsi.interface_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")+
  ylim(-log10(0.05),NA)
ggsave("Pancancer_interface_v_all_stemness_significance_mDNAsi_lwr_bounded_labeled.png",width=10,height=10,dpi=600)

#interface_v_gene_no_interface
results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.interface_v_gene_no_interface.wilcox.BH),
                      y=-log10(mDNAsi.interface_v_gene_no_interface.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")
ggsave("Pancancer_interface_v_gene_no_interface_stemness_significance_labeled.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.interface_v_gene_no_interface.wilcox.BH),
                      y=-log10(mDNAsi.interface_v_gene_no_interface.wilcox.BH),
                      color=In.Reactome)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1") +
  ggplot2::geom_hline(yintercept = -log10(0.05)) +
  ggplot2::geom_vline(xintercept = -log10(0.05))
ggsave("Pancancer_interface_v_gene_no_interface_stemness_significance_hvlines.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.interface_v_gene_no_interface.wilcox.BH),
                      y=-log10(mDNAsi.interface_v_gene_no_interface.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")+
  xlim(-log10(0.05),NA)+
  ylim(-log10(0.05),NA)
ggsave("Pancancer_interface_v_gene_no_interface_stemness_significance_lwr_bounded_labeled.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.interface_v_gene_no_interface.wilcox.BH),
                      y=-log10(mDNAsi.interface_v_gene_no_interface.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point() +
  scale_color_brewer(palette="Set1") +
  xlim(-log10(0.05),NA)
ggsave("Pancancer_interface_v_gene_no_interface_stemness_significance_mRNAsi_lwr_bounded.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.interface_v_gene_no_interface.wilcox.BH),
                      y=-log10(mDNAsi.interface_v_gene_no_interface.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")+
  xlim(-log10(0.05),NA)
ggsave("Pancancer_interface_v_gene_no_interface_stemness_significance_mRNAsi_lwr_bounded_labeled.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.interface_v_gene_no_interface.wilcox.BH),
                      y=-log10(mDNAsi.interface_v_gene_no_interface.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point() +
  scale_color_brewer(palette="Set1") +
  ylim(-log10(0.05),NA)
ggsave("Pancancer_interface_v_gene_no_interface_stemness_significance_mDNAsi_lwr_bounded.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.interface_v_gene_no_interface.wilcox.BH),
                      y=-log10(mDNAsi.interface_v_gene_no_interface.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")+
  ylim(-log10(0.05),NA)
ggsave("Pancancer_interface_v_gene_no_interface_stemness_significance_mDNAsi_lwr_bounded_labeled.png",width=10,height=10,dpi=600)

#gene_v_all
results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.gene_v_all.wilcox.BH),
                      y=-log10(mDNAsi.gene_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")
ggsave("Pancancer_gene_v_all_stemness_significance_labeled.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.gene_v_all.wilcox.BH),
                      y=-log10(mDNAsi.gene_v_all.wilcox.BH),
                      color=In.Reactome)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1") +
  ggplot2::geom_hline(yintercept = -log10(0.05)) +
  ggplot2::geom_vline(xintercept = -log10(0.05))
ggsave("Pancancer_gene_v_all_stemness_significance_hvlines.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.gene_v_all.wilcox.BH),
                      y=-log10(mDNAsi.gene_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")+
  xlim(-log10(0.05),NA)+
  ylim(-log10(0.05),NA)
ggsave("Pancancer_gene_v_all_stemness_significance_lwr_bounded_labeled.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.gene_v_all.wilcox.BH),
                      y=-log10(mDNAsi.gene_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point() +
  scale_color_brewer(palette="Set1") +
  xlim(-log10(0.05),NA)
ggsave("Pancancer_gene_v_all_stemness_significance_mRNAsi_lwr_bounded.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.gene_v_all.wilcox.BH),
                      y=-log10(mDNAsi.gene_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")+
  xlim(-log10(0.05),NA)
ggsave("Pancancer_gene_v_all_stemness_significance_mRNAsi_lwr_bounded_labeled.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.gene_v_all.wilcox.BH),
                      y=-log10(mDNAsi.gene_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point() +
  scale_color_brewer(palette="Set1") +
  ylim(-log10(0.05),NA)
ggsave("Pancancer_gene_v_all_stemness_significance_mDNAsi_lwr_bounded.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.gene_v_all.wilcox.BH),
                      y=-log10(mDNAsi.gene_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")+
  ylim(-log10(0.05),NA)
ggsave("Pancancer_gene_v_all_stemness_significance_mDNAsi_lwr_bounded_labeled.png",width=10,height=10,dpi=600)

#gene_no_interface_v_all
results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.gene_no_interface_v_all.wilcox.BH),
                      y=-log10(mDNAsi.gene_no_interface_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")
ggsave("Pancancer_gene_no_interface_v_all_stemness_significance_labeled.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.gene_no_interface_v_all.wilcox.BH),
                      y=-log10(mDNAsi.gene_no_interface_v_all.wilcox.BH),
                      color=In.Reactome)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1") +
  ggplot2::geom_hline(yintercept = -log10(0.05)) +
  ggplot2::geom_vline(xintercept = -log10(0.05))
ggsave("Pancancer_gene_no_interface_v_all_stemness_significance_hvlines.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.gene_no_interface_v_all.wilcox.BH),
                      y=-log10(mDNAsi.gene_no_interface_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")+
  xlim(-log10(0.05),NA)+
  ylim(-log10(0.05),NA)
ggsave("Pancancer_gene_no_interface_v_all_stemness_significance_lwr_bounded_labeled.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.gene_no_interface_v_all.wilcox.BH),
                      y=-log10(mDNAsi.gene_no_interface_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point() +
  scale_color_brewer(palette="Set1") +
  xlim(-log10(0.05),NA)
ggsave("Pancancer_gene_no_interface_v_all_stemness_significance_mRNAsi_lwr_bounded.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.gene_no_interface_v_all.wilcox.BH),
                      y=-log10(mDNAsi.gene_no_interface_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")+
  xlim(-log10(0.05),NA)
ggsave("Pancancer_gene_no_interface_v_all_stemness_significance_mRNAsi_lwr_bounded_labeled.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.gene_no_interface_v_all.wilcox.BH),
                      y=-log10(mDNAsi.gene_no_interface_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point() +
  scale_color_brewer(palette="Set1") +
  ylim(-log10(0.05),NA)
ggsave("Pancancer_gene_no_interface_v_all_stemness_significance_mDNAsi_lwr_bounded.png",width=10,height=10,dpi=600)

results_df2 %>%
  ggplot2::ggplot(aes(x=-log10(mRNAsi.gene_no_interface_v_all.wilcox.BH),
                      y=-log10(mDNAsi.gene_no_interface_v_all.wilcox.BH),
                      color=In.Reactome,
                      label=Interface)) +
  ggplot2::geom_point()+
  ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
  scale_color_brewer(palette="Set1")+
  ylim(-log10(0.05),NA)
ggsave("Pancancer_gene_no_interface_v_all_stemness_significance_mDNAsi_lwr_bounded_labeled.png",width=10,height=10,dpi=600)

if(FALSE){
  goi <- "IDH1"
  
  results_df2 %>%
    dplyr::filter(grepl(goi,Interface)) %>%
    ggplot2::ggplot(aes(x=-log10(mRNAsi.interface_v_gene_no_interface.wilcox.BH),
                        y=-log10(mDNAsi.interface_v_gene_no_interface.wilcox.BH),
                        color=In.Reactome,
                        label=Interface)) +
    ggplot2::geom_point()+
    ggplot2::geom_text(aes(label=Interface),size=2.5,hjust=0,vjust=0)+
    scale_color_brewer(palette="Set1")
  
  ggsave(paste("Pancancer_",goi,"_interface_v_gene_no_interface_stemness_significance_labeled.png",sep=""),width=10,height=10,dpi=600)
}
