# libs
library(survival)
library(magrittr)
library(dplyr)

# global vars
#DATA_DIR <- "/home/burkhart/Software/Ogmios/results/Mechismo/"
DATA_DIR <- "/Users/joshuaburkhart/SoftwareProjects/Ogmios/results/Mechismo/"
OUT_DIR <- DATA_DIR
DIST_FILE <- paste(DATA_DIR,"PAAD_PatientPairDistances.csv",sep="")

# load patient distance data
dist_df <- read.csv(DIST_FILE)

min_dist_df <- dist_df %>%
  dplyr::select(-Min.Shortest.Path) %>%
  reshape(idvar = "Patient.1",
          timevar = "Patient.2",
          direction = "wide")

avg_dist_df <- dist_df %>%
  dplyr::select(-Average.Distance) %>%
  reshape(idvar = "Patient.1",
          timevar = "Patient.2",
          direction = "wide")

min_dist_mtx <- as.matrix(min_dist_df)
avg_dist_mtx <- as.matrix(avg_dist_df)

rownames(min_dist_mtx) <- min_dist_mtx[,1]
rownames(avg_dist_mtx) <- avg_dist_mtx[,1]

min_dist_mtx <- min_dist_mtx[,-1]
avg_dist_mtx <- avg_dist_mtx[,-1]

# cluster
min_clust <- dist(min_dist_mtx) %>%
  hclust()
avg_clust <- dist(avg_dist_mtx) %>%
  hclust()

# plot
svg(filename=paste(OUT_DIR,"min_clust.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
  plot(min_clust)
dev.off()
svg(filename=paste(OUT_DIR,"avg_clust.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
  plot(avg_clust)
dev.off()
