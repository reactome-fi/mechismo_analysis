#libs
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(foreach)
library(doParallel)

#globs
MECH_INTERFACES <- "/Users/joshuaburkhart/Downloads/tcga_mechismo_stat_pancancer_undirected_significant.tsv"
RNA_SEQ_DIR <- "/Users/joshuaburkhart/Downloads/stddata__2016_07_15/pancancer_rnaseq/"
REACTOME_FIS <- "/Users/joshuaburkhart/Downloads/ProteinFIsInReactions_073118.txt"
CANCER_CENSUS <- "/Users/joshuaburkhart/Downloads/Census_allWed Aug 22 22_25_41 2018.tsv"
MECH_INPUT <- "/Users/joshuaburkhart/Downloads/TCGA_mech_input.tsv"

cores=detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

#helper functions
load_data_dir_2_df <- function(data_dir){
  data_files <- dir(data_dir,
                    pattern='\\.data.txt',
                    full.names = TRUE)
  data_tables <- lapply(data_files,
                        read.delim,
                        check.names=FALSE)
  do.call(cbind,data_tables)
}

#load data
if(file.exists("rna_seq_df.rda")){
  load("rna_seq_df.rda")
}else{
  rna_seq_df <- load_data_dir_2_df(RNA_SEQ_DIR)
  rna_seq_df <- rna_seq_df[-1,] %>% #remove 'gene_id' and 'normalized_count' row
    t()
  col_names <- rna_seq_df[1,] #'gene name|entrez gene id'
  rna_seq_df <- rna_seq_df[-1,] #remove 'gene name|entrez gene id' row
  row_names <- rownames(rna_seq_df)
  rna_seq_df <- matrix(as.numeric(unlist(rna_seq_df)),nrow=nrow(rna_seq_df))
  rownames(rna_seq_df) <- lapply(row_names,
                                 str_extract,
                                 pattern="^TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %>%
    unlist()
  colnames(rna_seq_df) <- col_names
  save(file="rna_seq_df.rda",rna_seq_df)
}

mech_interfaces_df <- read.delim(MECH_INTERFACES,
                                 stringsAsFactors = FALSE,
                                 header = FALSE)
reactome_fis_df <- read.delim(REACTOME_FIS,
                              stringsAsFactors = FALSE)
cancer_census_df <- read.delim(CANCER_CENSUS,
                               stringsAsFactors = FALSE)
mech_input_df <- read.delim(MECH_INPUT,
                            stringsAsFactors = FALSE)
#wrangle
mech_interfaces_df2 <- mech_interfaces_df %>%
  dplyr::arrange(desc(V7))

reactome_fis_df2 <- reactome_fis_df %>%
  dplyr::mutate(fwd_fi = paste(Gene1,"-",Gene2,sep=""),
                rev_fi = paste(Gene2,"-",Gene1,sep="")) %>%
  dplyr::select(fwd_fi,rev_fi)

cancer_census_df <- cancer_census_df %>%
  dplyr::mutate(Hallmark = ifelse(Hallmark == "Yes",TRUE,FALSE))

if(file.exists("pancancer_gene_samples_hash.rda")){
  load("pancancer_gene_samples_hash.rda")  
}else{
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
}

#loop setup
results_df <- data.frame(Interface = character(),
                         DE.Gene = character(),
                         DE.Gene.Wilcox.p = numeric(),
                         In.Reactome = logical(),
                         G1.In.Cancer.Census = logical(),
                         G2.In.Cancer.Census = logical(),
                         G1.Tier = numeric(),
                         G2.Tier = numeric(),
                         G1.Hallmark = logical(),
                         G2.Hallmark = logical(),
                         N.Interface.Samples = integer(),
                         Interface.Samples = character())

de_gene_ids <- colnames(rna_seq_df)
num_interfaces <- nrow(mech_interfaces_df2)

#calculate results
results_df <- foreach(i=1:num_interfaces,.combine=rbind) %dopar% {
  print(paste(i," of ",nrow(mech_interfaces_df2),sep=""))
  gene1 <- mech_interfaces_df2[i,1]
  gene2 <- mech_interfaces_df2[i,2]
  interface <- paste(gene1,
                     "-",
                     gene2,
                     sep="")
  
  if(interface %in% reactome_fis_df2$fwd_fi |
     interface %in% reactome_fis_df2$rev_fi){
    print(".")
  }else{
    interface_samples <- stringr::str_extract_all(mech_interfaces_df2[i,15],
                                                  "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %>%
      unlist()
    interface_samples <- intersect(interface_samples,rownames(rna_seq_df))
    
    gene_samples <- union(gene_samples_hash[[gene1]],gene_samples_hash[[gene2]])
    gene_samples <- intersect(gene_samples,rownames(rna_seq_df))
    
    gene_no_interface_samples <- setdiff(gene_samples,interface_samples)
    gene_no_interface_samples <- intersect(gene_no_interface_samples,rownames(rna_seq_df))
    
    for(j in 1:ncol(rna_seq_df)){
      
      # plot genes of interest / most significant findingings
      if(FALSE){
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
      }
      
      tryCatch(
        {
          results_df <- data.frame(Interface = interface,
                                  DE.Gene = de_gene_ids[j],
                                  DE.Gene.Wilcox.p = wilcox.test(rna_seq_df[interface_samples,j],
                                                                 rna_seq_df[-interface_samples,j])$p,
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
                                  N.Interface.Samples = length(interface_samples),
                                  Interface.Samples =
                                    paste(interface_samples,
                                          collapse = ','))
          #results_df <- results_df %>%
          #  rbind(result_df)
        },
        error=function(cond){
          print(error)
        }
      )
    }
  }
  results_df
}
stopCluster(cl)

#adjust p-values
results_df2 <- results_df %>%
  dplyr::mutate(DE.Gene.Wilcox.BH.Adj.p = p.adjust(DE.Gene.Wilcox.p,method="BH")) %>%
  dplyr::arrange(In.Reactome)

save(results_df2,file="pancancer_de_results.rda")

results_df2 %>%
  dplyr::select(Interface,
                DE.Gene,
                DE.Gene.Wilcox.p,
                DE.Gene.Wilcox.BH.Adj.p,
                In.Reactome,
                G1.In.Cancer.Census,
                G2.In.Cancer.Census,
                G1.Tier,
                G2.Tier,
                G1.Hallmark,
                G2.Hallmark,
                N.Interface.Samples,
                Interface.Samples) %>%
  write.table("Pancancer_de_interface_v_all_Stemness.tsv",
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
