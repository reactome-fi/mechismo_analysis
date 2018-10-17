#libs
library(dplyr)
library(magrittr)
library(stringr)
library(foreach)
library(doParallel)

#globs
N_INTERFACE_SAMPLES_THRESH <- 5
N_UNUSED_CORES <- 0
SLURM_PROCID <- as.numeric(Sys.getenv("SLURM_PROCID"))
SLURM_JOB_NUM_NODES <- as.numeric(Sys.getenv("SLURM_NTASKS"))
LUSTRE_DIR <- "/home/exacloud/lustre1/WongLab/tmp"

if (!is.na(SLURM_PROCID)) {
  load("data/prna_seq_df.rda")
  load("data/pmech_interfaces_df.rda")
  load("data/preactome_fis_df.rda")
  load("data/pcancer_census_df.rda")
  load("data/pmech_input_df.rda")
  load("data/pgene_samples_hash.rda")
  
  cores = detectCores()
  cl <- makeCluster(cores[1] - N_UNUSED_CORES)
  registerDoParallel(cl)
  
  #loop setup
  results_df <- data.frame(
    Interface = character(),
    DE.Gene = character(),
    DE.Gene.Wilcox.p = numeric(),
    Num.Interface.Samples = numeric(),
    Num.NoInterface.Samples = numeric(),
    Interface.NoInterface.Diff = numeric()
  )
  
  de_gene_ids <- colnames(rna_seq_df)
  num_interfaces <- nrow(mech_interfaces_df)
  
  #calculate results
  results_df <- foreach(
    i = 1:num_interfaces,
    .combine = rbind,
    .packages = "magrittr"
  ) %dopar% {
    if (mod(i, SLURM_JOB_NUM_NODES) == SLURM_PROCID) {
      #split job
      gene1 <- mech_interfaces_df[i, 1] #1 for pancancer file/2 for cancerwise
      gene2 <- mech_interfaces_df[i, 2] #2 for pancancer file/3 for cancerwise
      interface <- paste(gene1,
                         "-",
                         gene2,
                         sep = "")
      
      if (interface %in% reactome_fis_df$fwd_fi |
          interface %in% reactome_fis_df$rev_fi) {
        interface_samples <-
          stringr::str_extract_all(mech_interfaces_df[i, 15], #15 for pancancer file/16 for cancerwise
                                   "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %>%
          unlist()
        interface_samples <- intersect(interface_samples,
                                       rownames(rna_seq_df))
        interface_samples <-
          interface_samples[!is.na(interface_samples)]
        
        li <- length(interface_samples)
        
        if (li > N_INTERFACE_SAMPLES_THRESH) {
          no_interface_samples <- setdiff(rownames(rna_seq_df),
                                          interface_samples)
          no_interface_samples <-
            no_interface_samples[!is.na(no_interface_samples)]
          
          lni <- length(no_interface_samples)
        
          if (lni > N_INTERFACE_SAMPLES_THRESH) {
            for (j in 1:ncol(rna_seq_df)) {
              result_df <- data.frame(
                Interface = interface,
                DE.Gene = de_gene_ids[j],
                DE.Gene.Wilcox.p = wilcox.test(rna_seq_df[interface_samples, j],
                                               rna_seq_df[no_interface_samples, j])$p.value,
                Num.Interface.Samples = li,
                Num.NoInterface.Samples = lni,
                Interface.NoInterface.Diff = mean(rna_seq_df[interface_samples,j]) - mean(rna_seq_df[no_interface_samples,j])
              )
              results_df <- results_df %>%
                rbind(result_df)
            }
          }
        }
      }
    }
    results_df
  }
  stopCluster(cl)
  
  results_df %>%
    dplyr::select(Interface,
                  DE.Gene,
                  DE.Gene.Wilcox.p,
                  Num.Interface.Samples,
                  Num.NoInterface.Samples,
                  Interface.NoInterface.Diff) %>%
    write.table(
      paste(LUSTRE_DIR,"/pslurm_proc_",
            SLURM_PROCID,
            "_de.tsv",
            sep = ""),
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t"
    )
}
