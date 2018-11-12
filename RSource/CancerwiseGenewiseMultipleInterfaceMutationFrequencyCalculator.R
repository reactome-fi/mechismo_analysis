library(magrittr)

MECH_INPUT <- "/Users/joshuaburkhart/Downloads/TCGA_mech_input.tsv"
MECH_OUTPUT <- "/Users/joshuaburkhart/Downloads/tcga_mechismo_stat_cancer_wise_undirected_significant.tsv"
REACTOME_FIS <- "/Users/joshuaburkhart/Downloads/ProteinFIsInReactions_073118.txt"

if(!file.exists("mech_input_df.rda")){
  mech_input_df <- read.delim(MECH_INPUT,
                              stringsAsFactors = FALSE)
  save(mech_input_df,file="mech_input_df.rda")
}else{
  load("mech_input_df.rda")
}

if(!file.exists("mech_output_df.rda")){
  mech_output_df <- read.delim(MECH_OUTPUT,
                               stringsAsFactors = FALSE,
                               header = FALSE)
  save(mech_output_df,file="mech_output_df.rda")
}else{
  load("mech_output_df.rda")
}

if(!file.exists("reactome_fis_df.rda")){
  reactome_fis_df <- read.delim(REACTOME_FIS,
                                stringsAsFactors = FALSE) %>%
    dplyr::mutate(fwd_fi = paste(Gene1,"-",Gene2,sep=""),
                  rev_fi = paste(Gene2,"-",Gene1,sep="")) %>%
    dplyr::select(fwd_fi,rev_fi)
  save(reactome_fis_df,file="reactome_fis_df.rda")
}else{
  load("reactome_fis_df.rda")
}

if(!file.exists("cancerwise_gene_interfaces.rda") | !file.exists("cancer_type_samples.rda")){
  gene_interfaces <- list()
  cancer_type_samples <- list()
  for(i in 1:nrow(mech_output_df)){
    interface_samples_string <- mech_output_df[i,16]
    interface_samples <- stringr::str_extract_all(interface_samples_string,
                                                  "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %>%
      unlist()
    interface_genes <- c(mech_output_df[i,2],
                         mech_output_df[i,3]) %>%
      sort()
    cancer_type <- mech_output_df[i,1]
    cur_1_interfaces <- list()
    if(interface_genes[1] %in% names(gene_interfaces[[cancer_type]])){
      cur_1_interfaces <- gene_interfaces[[cancer_type]][[interface_genes[1]]]
    }
    cur_2_interfaces <- list()
    if(interface_genes[2] %in% names(gene_interfaces[[cancer_type]])){
      cur_2_interfaces <- gene_interfaces[[cancer_type]][[interface_genes[2]]]
    }
    
    interface_name <- paste(interface_genes[1],"-",interface_genes[2],sep="")
    
    if(interface_name %in% reactome_fis_df$fwd_fi |
       interface_name %in% reactome_fis_df$rev_fi){
      
      if(interface_name %in% names(cur_1_interfaces) |
         interface_name %in% names(cur_2_interfaces)){
        #this shouldn't happen unless interface names are duplicated out-of-order
        print("ERROR: cur_#_interfaces[[interface_name]] should be empty.")
        break
      }
      cur_1_interfaces[[interface_name]] <- interface_samples
      cur_2_interfaces[[interface_name]] <- interface_samples
      cancer_type_samples[[cancer_type]] <- union(cancer_type_samples[[cancer_type]],interface_samples)
      
      gene_interfaces[[cancer_type]][[interface_genes[1]]] <- cur_1_interfaces
      gene_interfaces[[cancer_type]][[interface_genes[2]]] <- cur_2_interfaces
    }
  }
  save(gene_interfaces,file="cancerwise_gene_interfaces.rda")
  save(cancer_type_samples,file="cancer_type_samples.rda")
}else{
  load("cancerwise_gene_interfaces.rda")
  load("cancer_type_samples.rda")
}

# not sure how this needs to change but know we'll need to parse samples out like this
if(!file.exists("multi_gene_samples_hash.rda")){
  multi_gene_samples_hash <- new.env()
  for(i in 1:nrow(mech_input_df)){
    gene <- mech_input_df[i,"Gene.symbol"]
    
    new_samples_string <- mech_input_df[i,"Sample.ID"]
    new_samples <- stringr::str_extract_all(new_samples_string,
                                            "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %>%
      unlist()
    
    cur_samples <- multi_gene_samples_hash[[gene]]
    if(is.null(cur_samples)){
      cur_samples <- list()
    }else{
      new_multi_mutant_samples <- intersect(cur_samples[["single"]],new_samples)
      cur_samples[["multi"]] <- union(cur_samples[["multi"]],new_multi_mutant_samples)
    }
    cur_samples[["single"]] <- union(cur_samples[["single"]],new_samples)
    multi_gene_samples_hash[[gene]] <- cur_samples
  }
  save(multi_gene_samples_hash,file="multi_gene_samples_hash.rda")
}else{
  load("multi_gene_samples_hash.rda")
}

genes <- names(multi_gene_samples_hash)
gene_interface_mutation_distribution_df <- data.frame("Cancer.Type" = character(),
                                                      "Gene" = character(),
                                                      "Num.Singles.Samples" = numeric(),
                                                      "Num.Multis.Samples" = numeric(),
                                                      "Num.Interfaces" = numeric(),
                                                      "Single.v.Multi.Fisher" = numeric())
for(i in 1:length(genes)){
  gene <- genes[[i]]
  multi_samples <- multi_gene_samples_hash[[gene]][["multi"]]
  if(!is.null(multi_samples) & length(multi_samples) > 0){
    single_samples <- setdiff(multi_gene_samples_hash[[gene]][["single"]],
                              multi_gene_samples_hash[[gene]][["multi"]])
    if(length(single_samples) > 0){
      singles_inteface_counts <- list()
      multis_interface_counts <- list()
      for(k in 1:length(names(cancer_type_samples))){
        cancer_type <- names(cancer_type_samples)[k]
        interfaces <- names(gene_interfaces[[cancer_type]][[gene]])
        if(length(interfaces) > 1){
          for(j in 1:length(interfaces)){
            interface <- interfaces[j]
            interface_singles <- intersect(gene_interfaces[[cancer_type]][[gene]][[interface]],
                                           single_samples)
            interface_singles <- intersect(cancer_type_samples[[cancer_type]],
                                           single_samples)
            interface_multis <- intersect(gene_interfaces[[cancer_type]][[gene]][[interface]],
                                          multi_samples)
            interface_multis <- intersect(cancer_type_samples[[cancer_type]],
                                          multi_samples)
            singles_count <- length(interface_singles)
            multis_count <- length(interface_multis)#doubles aren't counted twice but should be
            
            singles_inteface_counts[[cancer_type]] <- c(singles_inteface_counts[[cancer_type]],
                                                        singles_count)
            multis_interface_counts[[cancer_type]] <- c(multis_interface_counts[[cancer_type]],
                                                        multis_count)
            if(length(singles_inteface_counts[[cancer_type]]) !=
               length(multis_interface_counts[[cancer_type]])){
              print("ERROR: interface counts should match")
              break
            }
          }
        }
      }
      
      cancer_types <- names(singles_inteface_counts)
      if(!is.null(cancer_types) & length(cancer_types) > 0){
        for(l in 1:length(cancer_types)){
          cancer_type <- names(singles_inteface_counts)[l]
          if(length(singles_inteface_counts[[cancer_type]]) > 1 &
             length(singles_inteface_counts[[cancer_type]]) == length(multis_interface_counts) &
             sum(singles_inteface_counts[[cancer_type]]) > 0 &
             sum(multis_interface_counts[[cancer_type]]) > 0){ # fisher test requires 2 or more non-zero row marginals
            fisher.prep <- rbind(singles_inteface_counts[[cancer_type]],
                                 multis_interface_counts[[cancer_type]])
            if(sum(as.numeric(colSums(fisher.prep) > 1)) > 1){
              
              tmp_df <- data.frame("Cancer.Type" = cancer_type,
                                   "Gene" = gene,
                                   "Num.Singles.Samples" = length(intersect(single_samples,cancer_type_samples[[cancer_type]])),
                                   "Num.Multis.Samples" = length(intersect(multi_samples,cancer_type_samples[[cancer_type]])),
                                   "Num.Interfaces" = length(gene_interfaces[[cancer_type]][[gene]]),
                                   "Single.v.Multi.Fisher" = fisher.test(fisher.prep,
                                                                         simulate.p.value = TRUE,
                                                                         B = 1e4)$p.value)
              gene_interface_mutation_distribution_df <-
                rbind(gene_interface_mutation_distribution_df,
                      tmp_df)
            }
          }
        }
      }
    }
  }
}

gene_interface_mutation_distribution_df %>%
  dplyr::mutate(Single.v.Multi.Fisher.BH.Adj.p = p.adjust(Single.v.Multi.Fisher,method = "BH")) %>%
  write.table("cancerwise_gene_interface_single_v_multi_mutation_distribution.tsv",
              row.names = FALSE,
              sep="\t")
