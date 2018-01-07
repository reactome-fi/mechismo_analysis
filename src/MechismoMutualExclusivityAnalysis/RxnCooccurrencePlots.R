
library(magrittr)
library(EcoSimR)
library(dplyr)
library(plyr)
library(reshape2)
library(igraph)
library(pheatmap)

SIG_RXN_FILEPATH <- "/Users/joshuaburkhart/SoftwareProjects/Ogmios/datasets/MechismoSamplesToReactions_103017.txt"
OUT_DIR <- "/Users/joshuaburkhart/SoftwareProjects/Ogmios/results/Mechismo/"

sig_rxn_df <- read.csv(SIG_RXN_FILEPATH,header = TRUE,sep="\t")

sig_rxn_df <- sig_rxn_df %>%
  dplyr::select(-Cancer)

sig_rxn_df_t <- sig_rxn_df  %>%
  dplyr::select(-Sample) %>%
  t()

colnames(sig_rxn_df_t) <- sig_rxn_df$Sample

#from https://cran.r-project.org/web/packages/EcoSimR/vignettes/CoOccurrenceVignette.html
myModel <- cooc_null_model(speciesData=sig_rxn_df_t,suppressProg=TRUE)
summary(myModel)

plot(myModel,type="hist")

plot(myModel,type="cooc")

plot(myModel,type="burn_in")

#from https://stackoverflow.com/questions/13281303/creating-co-occurrence-matrix
dat2 <- reshape2::melt(sig_rxn_df_t)
w <- dcast(dat2, Var1~Var2)
x <- as.matrix(w[,-1])
x[is.na(x)] <- 0
x <- apply(x, 2,  function(x) as.numeric(x > 0))  #recode as 0/1
v <- x %*% t(x)                                   #the magic matrix 
diag(v) <- 0                                      #repalce diagonal
dimnames(v) <- list(w[, 1], w[,1])                #name the dimensions
g <- graph.adjacency(v, weighted=TRUE, mode ='undirected')
g <- simplify(g)
# set labels and degrees of vertices
V(g)$label <- V(g)$name
V(g)$degree <- degree(g)
plot(g)


svg(filename=paste(OUT_DIR,"SigRxnCooccurrenceAll.svg",sep=""),
    width=11,
    height=11,
    pointsize=10)
v %>% pheatmap::pheatmap()
dev.off()

system_command = paste("qlmanage -t -s 1366 -o ",
                       OUT_DIR,
                       " ",
                       OUT_DIR,
                       "*.svg", 
                       sep="")
system(paste(system_command))

cols_to_keep = colSums(v) > 1000
v <- v[cols_to_keep,cols_to_keep]

svg(filename=paste(OUT_DIR,"SigRxnCooccurrence.svg",sep=""),
    width=11,
    height=11,
    pointsize=10)
v %>% pheatmap::pheatmap()
dev.off()


system_command = paste("qlmanage -t -s 1366 -o ",
                       OUT_DIR,
                       " ",
                       OUT_DIR,
                       "*.svg", 
                       sep="")
system(paste(system_command))
