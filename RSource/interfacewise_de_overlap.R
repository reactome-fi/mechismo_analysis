
library(dplyr)

P.THRESH <- 0.05

CANCERTYPES <- c("CESC",
                 "TGCT",
                 "GBM",
                 "ESCA",
                 "PAAD",
                 "LIHC",
                 "LGG",
                 "BLCA",
                 "LUSC",
                 "GBMLGG",
                 "UCS",
                 "HNSC",
                 "BRCA",
                 "PCPG",
                 "STAD",
                 "STES",
                 "THCA",
                 "SKCM",
                 "UVM")

all_interfaces_mutex_df <- data.frame(Cancer.Type = character(),
                                      Interface.1 = character(),
                                      Interface.2 = character(),
                                      Both = numeric(),
                                      Interface.1.Only = numeric(),
                                      Interface.2.Only = numeric(),
                                      Neither = numeric(),
                                      Fisher.Test.p = numeric(),
                                      Fisher.Test.BH.Adj.p = numeric())

for(k in 1:length(CANCERTYPES)){
  
  CANCERTYPE <- CANCERTYPES[k]
  
  CANCERTYPE_INTERFACEWISE_DE <- paste("/Users/joshuaburkhart/",CANCERTYPE,"_interfacewise_de.tsv",sep="")
  
  cancertype_interfacewise_df <- read.delim(CANCERTYPE_INTERFACEWISE_DE)
  
  sig_gene_expr <- cancertype_interfacewise_df %>%
    dplyr::filter(DE.Gene.Wilcox.BH.Adj.p < P.THRESH)
  
  interface_mutex_df <- data.frame(Cancer.Type = character(),
                                   Interface.1 = character(),
                                   Interface.2 = character(),
                                   Both = numeric(),
                                   Interface.1.Only = numeric(),
                                   Interface.2.Only = numeric(),
                                   Neither = numeric(),
                                   Fisher.Test.p = numeric())
  
  interfaces <- as.character(unique(sig_gene_expr$Interface))
  
  if(length(interfaces > 1)){
    
    for(i in 1:length(interfaces)){
      interface.1 <- interfaces[i]
      
      for(j in 1:length(interfaces)){
        interface.2 <- interfaces[j]
        if(interface.1 != interface.2){
          
          interface.1.de <- sig_gene_expr %>%
            dplyr::filter(as.character(Interface) == interface.1)
          
          interface.2.de <- sig_gene_expr %>%
            dplyr::filter(as.character(Interface) == interface.2)
          
          both <- intersect(as.character(interface.1.de$DE.Gene),
                            as.character(interface.2.de$DE.Gene))
          
          interface.1.only <- setdiff(as.character(interface.1.de$DE.Gene),
                                      as.character(interface.2.de$DE.Gene))
          
          interface.2.only <- setdiff(as.character(interface.2.de$DE.Gene),
                                      as.character(interface.1.de$DE.Gene))
          
          neither <- as.character(cancertype_interfacewise_df$DE.Gene) %>%
            setdiff(.,as.character(interface.1.de$DE.Gene)) %>%
            setdiff(.,as.character(interface.2.de$DE.Gene))
          
          z <- matrix(c(length(neither),
                        length(interface.1.only),
                        length(interface.2.only),
                        length(both)),
                      nrow=2,
                      dimnames=list(c("no.interface.1",
                                      "interface.1"),
                                    c("no.interface.2",
                                      "interface.2")))
          
          z
          
          y <- fisher.test(z,alternative = "l")
          
          interface_mutex_row <- data.frame(Cancer.Type = CANCERTYPE,
                                            Interface.1 = interface.1,
                                            Interface.2 = interface.2,
                                            Both = length(both),
                                            Interface.1.Only = length(interface.1.only),
                                            Interface.2.Only = length(interface.2.only),
                                            Neither = length(neither),
                                            Fisher.Test.p = y$p.value)
          interface_mutex_df <- rbind(interface_mutex_df,
                                      interface_mutex_row)
        }
      }
    }
    
    interface_mutex_df <- interface_mutex_df %>%
      dplyr::mutate(Fisher.Test.BH.Adj.p = p.adjust(Fisher.Test.p,method="BH")) %>%
      dplyr::arrange(Fisher.Test.p)
    
    all_interfaces_mutex_df <- rbind(all_interfaces_mutex_df,
                                     interface_mutex_df)
  }
}
all_interfaces_mutex_df %>%
  write.table(file = "interfacewise_de_mutex.tsv",
              row.names = FALSE,
              sep="\t")
