library(magrittr)

MECH_INPUT <- "/Users/joshuaburkhart/Downloads/TCGA_mech_input.tsv"
MECH_OUTPUT <- "/Users/joshuaburkhart/Downloads/tcga_mechismo_stat_pancancer_undirectesd.tsv"

if(!file.exists("mech_input_df.rda")){
  mech_input_df <- read.delim(MECH_INPUT,
                              stringsAsFactors = FALSE)
  save(mech_input_df,file="mech_input_df.rda")
}else{
  load("mech_input_df.rda")
}

if(!file.exists("mech_output_df.rda")){
  mech_output_df <- read.delim(MECH_OUTPUT,
                               stringsAsFactors = FALSE)
  save(mech_output_df,file="mech_output_df.rda")
}else{
  load("mech_output_df.rda")
}

if(!file.exists("gene_interfaces.rda")){
  gene_interfaces <- list()
  for(i in 1:nrow(mech_output_df)){
    interface_samples_string <- mech_output_df[i,"Unique.Samples"]
    interface_samples <- stringr::str_extract_all(interface_samples_string,
                                                  "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %>%
      unlist()
    interface_genes <- c(mech_output_df[i,"Partner.A"],
                         mech_output_df[i,"Partner.B"]) %>%
      sort()
    
    cur_1_interfaces <- list()
    if(interface_genes[1] %in% names(gene_interfaces)){
      cur_1_interfaces <- gene_interfaces[[interface_genes[1]]]
    }
    cur_2_interfaces <- list()
    if(interface_genes[2] %in% names(gene_interfaces)){
      cur_2_interfaces <- gene_interfaces[[interface_genes[2]]]
    }
    
    interface_name <- paste(interface_genes[1],"-",interface_genes[2],sep="")
    
    if(interface_name %in% names(cur_1_interfaces) |
       interface_name %in% names(cur_2_interfaces)){
      #this shouldn't happen unless interface names are duplicated out-of-order
      print("ERROR: cur_#_interfaces[[interface_name]] should be empty.")
      break
    }
    cur_1_interfaces[[interface_name]] <- interface_samples
    cur_2_interfaces[[interface_name]] <- interface_samples
    
    gene_interfaces[[interface_genes[1]]] <- cur_1_interfaces
    gene_interfaces[[interface_genes[2]]] <- cur_2_interfaces
  }
  save(gene_interfaces,file="gene_interfaces.rda")
}else{
  load("gene_interfaces.rda")
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
gene_interface_mutation_distribution_df <- data.frame("Gene" = character(),
                                                      "Num.Singles" = numeric(),
                                                      "Num.Multis" = numeric(),
                                                      "Num.Interfaces" = numeric(),
                                                      "Single.v.Multi.Wilcox" = numeric())
for(i in 1:length(genes)){
  gene <- genes[[i]]
  multi_samples <- multi_gene_samples_hash[[gene]][["multi"]]
  if(!is.null(multi_samples) & length(multi_samples) > 0){
    single_samples <- setdiff(multi_gene_samples_hash[[gene]][["single"]],
                              multi_gene_samples_hash[[gene]][["multi"]])
    if(length(single_samples) > 0){
      interfaces <- names(gene_interfaces[[gene]])
      if(length(interfaces) > 1){
        singles_inteface_ratios <- numeric()
        multis_interface_ratios <- numeric()
        for(j in 1:length(interfaces)){
          interface <- interfaces[j]
          interface_singles <- intersect(gene_interfaces[[gene]][[interface]],
                                         single_samples)
          interface_multis <- intersect(gene_interfaces[[gene]][[interface]],
                                        multi_samples)
          singles_ratio <- length(interface_singles) / length(single_samples)
          multis_ratio <- length(interface_multis) / length(multi_samples)
          singles_inteface_ratios <- c(singles_inteface_ratios,
                                           singles_ratio)
          multis_interface_ratios <- c(multis_interface_ratios,
                                           multis_ratio)
        }
        tmp_df <- data.frame("Gene" = gene,
                             "Num.Singles" = length(single_samples),
                             "Num.Multis" = length(multi_samples),
                             "Num.Interfaces" = length(interfaces),
                             "Single.v.Multi.Wilcox" = wilcox.test(x = singles_inteface_ratios,
                                                                   y = multis_interface_ratios,
                                                                   paired = TRUE,
                                                                   alternative = "two.sided")$p.value)
        gene_interface_mutation_distribution_df <-
          rbind(gene_interface_mutation_distribution_df,
                tmp_df)
      }
    }
  }
}

gene_interface_mutation_distribution_df %>%
  dplyr::mutate(Single.v.Multi.Wilcox.BH.Adj.p = p.adjust(Single.v.Multi.Wilcox,method = "BH")) %>%
  write.table("gene_interface_single_v_multi_mutation_distribution.tsv",
              row.names = FALSE,
              sep="\t")
