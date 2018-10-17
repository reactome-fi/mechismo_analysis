#libs
library(dplyr)
library(magrittr)
library(stringr)
library(genefilter)

#globs
MECH_INTERFACES <- "/Users/joshuaburkhart/Downloads/tcga_mechismo_stat_pancancer_undirected_significant.tsv"
RNA_SEQ_DIR <- "/Users/joshuaburkhart/Downloads/stddata__2016_07_15/pancancer_rnaseq/"
REACTOME_FIS <- "/Users/joshuaburkhart/Downloads/ProteinFIsInReactions_073118.txt"
CANCER_CENSUS <- "/Users/joshuaburkhart/Downloads/Census_allWed Aug 22 22_25_41 2018.tsv"
MECH_INPUT <- "/Users/joshuaburkhart/Downloads/TCGA_mech_input.tsv"


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
if(!file.exists("prna_seq_df.rda")){
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
  rna_seq_df <- rna_seq_df[,-which(grepl("\\?",colnames(rna_seq_df)))]
  #rna_seq_df <- t(genefilter::varFilter(t(rna_seq_df))) # Variance filter fails with large data matrix
  save(file="prna_seq_df.rda",rna_seq_df)
}

if(!file.exists("pmech_interfaces_df.rda")){
  mech_interfaces_df <- read.delim(MECH_INTERFACES,
                                   stringsAsFactors = FALSE,
                                   header = FALSE)
  save(mech_interfaces_df,file="pmech_interfaces_df.rda")
}

if(!file.exists("preactome_fis_df.rda")){
  reactome_fis_df <- read.delim(REACTOME_FIS,
                                stringsAsFactors = FALSE) %>%
    dplyr::mutate(fwd_fi = paste(Gene1,"-",Gene2,sep=""),
                  rev_fi = paste(Gene2,"-",Gene1,sep="")) %>%
    dplyr::select(fwd_fi,rev_fi)
  save(reactome_fis_df,file="preactome_fis_df.rda")
}

if(!file.exists("pcancer_census_df.rda")){
  cancer_census_df <- read.delim(CANCER_CENSUS,
                                 stringsAsFactors = FALSE) %>%
    dplyr::mutate(Hallmark = ifelse(Hallmark == "Yes",TRUE,FALSE))
  save(cancer_census_df,file="pcancer_census_df.rda")
}

if(!file.exists("pmech_input_df.rda")){
  mech_input_df <- read.delim(MECH_INPUT,
                              stringsAsFactors = FALSE)
  save(mech_input_df,file="pmech_input_df.rda")
}else{
  load("pmech_input_df.rda")
}

if(!file.exists("pgene_samples_hash.rda")){
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
  save(gene_samples_hash,file="pgene_samples_hash.rda")
}