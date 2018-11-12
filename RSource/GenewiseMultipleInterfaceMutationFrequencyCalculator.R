library(magrittr)
library(ggplot2)

MECH_INPUT <- "/Users/joshuaburkhart/Downloads/TCGA_mech_input.tsv"
MECH_OUTPUT <- "/Users/joshuaburkhart/Downloads/tcga_mechismo_stat_pancancer_undirected_significant.tsv"
REACTOME_FIS <- "/Users/joshuaburkhart/Downloads/ProteinFIsInReactions_073118.txt"
MIN_SAMPLES <- 2

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

if(!file.exists("gene_interfaces.rda")){
  gene_interfaces <- list()
  for(i in 1:nrow(mech_output_df)){
    interface_samples_string <- mech_output_df[i,15]
    interface_samples <- stringr::str_extract_all(interface_samples_string,
                                                  "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %>%
      unlist()
    interface_genes <- c(mech_output_df[i,1],
                         mech_output_df[i,2]) %>%
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
    
    #if(interface_name %in% reactome_fis_df$fwd_fi |
    #   interface_name %in% reactome_fis_df$rev_fi){
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
    #}
  save(gene_interfaces,file="gene_interfaces.rda")
}else{
  load("gene_interfaces.rda")
}

# not sure how this needs to change but know we'll need to parse samples out like this
if(!file.exists("pan_multi_gene_samples_hash.rda")){
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
  save(multi_gene_samples_hash,file="pan_multi_gene_samples_hash.rda")
}else{
  load("pan_multi_gene_samples_hash.rda")
}

genes <- names(multi_gene_samples_hash)
gene_interface_mutation_distribution_df <- data.frame("Gene" = character(),
                                                      "Num.Singles.Samples" = numeric(),
                                                      "Num.Multis.Samples" = numeric(),
                                                      "Num.Interfaces" = numeric(),
                                                      "Single.v.Multi.Fisher" = numeric())
for(i in 1:length(genes)){
  gene <- genes[[i]]
  multi_samples <- multi_gene_samples_hash[[gene]][["multi"]]
  if(!is.null(multi_samples) & length(multi_samples) >= MIN_SAMPLES){
    single_samples <- setdiff(multi_gene_samples_hash[[gene]][["single"]],
                              multi_gene_samples_hash[[gene]][["multi"]])
    if(length(single_samples) >= MIN_SAMPLES){
      interfaces <- names(gene_interfaces[[gene]])
      if(length(interfaces) > 1){
        singles_inteface_counts <- numeric()
        multis_interface_counts <- numeric()
        for(j in 1:length(interfaces)){
          interface <- interfaces[j]
          interface_singles <- intersect(gene_interfaces[[gene]][[interface]],
                                         single_samples)
          interface_multis <- intersect(gene_interfaces[[gene]][[interface]],
                                        multi_samples)
          singles_count <- length(interface_singles)
          multis_count <- length(interface_multis)#doubles aren't counted twice but should be
          singles_inteface_counts <- c(singles_inteface_counts,
                                       singles_count)
          multis_interface_counts <- c(multis_interface_counts,
                                       multis_count)
          if(length(singles_inteface_counts) != length(multis_interface_counts)){
            print("ERROR: interface counts should match")
            break
          }
        }
        
        if(length(singles_inteface_counts) > 1 &
           length(singles_inteface_counts) == length(multis_interface_counts) &
           sum(singles_inteface_counts) > 0 &
           sum(multis_interface_counts) > 0){ # fisher test requires 2 or more non-zero row marginals
          fisher.prep <- rbind(singles_inteface_counts,
                               multis_interface_counts)
          if(sum(as.numeric(colSums(fisher.prep) > 1)) > 1){
          
          tmp_df <- data.frame("Gene" = gene,
                               "Num.Singles.Samples" = length(single_samples),
                               "Num.Multis.Samples" = length(multi_samples),
                               "Num.Interfaces" = length(interfaces),
                               "Single.v.Multi.Fisher" = fisher.test(t(fisher.prep),
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

gene_interface_mutation_distribution_df %>%
  dplyr::mutate(Single.v.Multi.Fisher.BH.Adj.p = p.adjust(Single.v.Multi.Fisher,method = "BH")) %>%
  write.table("gene_interface_single_v_multi_mutation_distribution.tsv",
              row.names = FALSE,
              sep="\t")

pmulti_samples <- multi_gene_samples_hash[["PIK3CA"]][["multi"]]
psingle_samples <- setdiff(multi_gene_samples_hash[["PIK3CA"]][["single"]],
                           multi_gene_samples_hash[["PIK3CA"]][["multi"]])
pinterfaces <- names(gene_interfaces[["PIK3CA"]])
pdists_df <- data.frame("Interface" = character(),
                        "Mutation.Group" = character(),
                        "Num.Samples" = numeric())
for(i in 1:length(pinterfaces)){
  tmp <- data.frame("Interface" = pinterfaces[i],
  "Mutation.Group" = "Single Mutation",
  "Num.Samples" = length(intersect(gene_interfaces[["PIK3CA"]][[pinterfaces[i]]],
                                 psingle_samples)))
  pdists_df <- rbind(pdists_df,tmp)
  
  tmp <- data.frame("Interface" = pinterfaces[i],
  "Mutation.Group" = "Multiple Mutations",
  "Num.Samples" = length(intersect(gene_interfaces[["PIK3CA"]][[pinterfaces[i]]],
                                  pmulti_samples)))
  pdists_df <- rbind(pdists_df,tmp)
}

pdists_df %>%
  ggplot(aes(x=Interface,y=Num.Samples,fill=Mutation.Group))+
      geom_col(position="stack")+
  theme(axis.text.x=element_text(angle = 290, hjust=0.0, vjust=1.0),
        plot.margin =unit( c(30,30,10,5),"points")
        )+
  scale_fill_discrete(guide = guide_legend(reverse=TRUE))
  

  
  
  
  
  
  
  
  


