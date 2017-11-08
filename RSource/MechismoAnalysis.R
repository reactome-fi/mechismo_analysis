plot.tcga.cancer <- function(tcga.dist.file, main = "TCGA Cancers") {
    # result.dir <- "/Users/wug/git/Ogmios/results"
    # tcga.dist.file <- paste(result.dir, "TCGACancerPairWiseReactionNetworkDist_092617.txt", sep="/")
    # tcga.dist.file <- paste(result.dir, "TCGACancerPairWiseReactionNetworkDist_100217.txt", sep="/")
    tcga.dist.df <-
        read.delim(
            tcga.dist.file,
            header = T,
            sep = "\t",
            check.names = FALSE
        )
    # print(names(tcga.dist.df))
    # convert the data frame into a matrix
    tcga.matrix <- as.matrix(tcga.dist.df[2:length(tcga.dist.df)])
    # print(tcga.matrix)
    tcga.dist <- as.dist(tcga.matrix)
    # print(tcga.dist)
    tcga.clust <- hclust(tcga.dist,
                         method = "ward.D")
    # tcga.clust <- hclust(reaction.dist, method = "ward.D")
    
    plot(
        tcga.clust,
        # main = "TCGA Cancers (29 in Total)",
        # main = "TCGA Cancers (32 in Total)",
        main = main,
        xlab = "Cancer Samples",
        ylab = "Network Distance",
        sub = "",
        cex = 0.40
    )
    
    # Get clusters manually
    clusters <- (identify(tcga.clust))
    clusters
}

# Use this function to create a data frame for survival analysis
create.survival.data <- function(clin.file, clusters, sample.names, exclude.no.clustered.samples = TRUE) {
    clin <- read.delim(clin.file,
                       header = TRUE,
                       sep = "\t",
                       check.names = FALSE)
    clin.samples <- clin[, 1]
    # Column for cluster
    # -1 for other samples are not subject to clustering
    if (!exclude.no.clustered.samples) {
        sample.clusters <- rep(-1, length(clin.samples))
    }
    else {
        # NA for other samples are not subject to clustering
        sample.clusters <- rep(NA, length(clin.samples))
    }
    # Column for duration
    duration <- rep(NA, length(clin.samples))
    # Parse the patient barcodes
    for (i in 1 : length(sample.names)) {
        sample.names[i] <- substr(sample.names[i], 0, 12)
    }
    clin.names <- as.character(clin.samples)
    for (i in 1:length(clin.samples)) {
        # This is just an index array
        sample <- clin.names[i]
        # Look for clusters for this sample
        for (j in 1:length(clusters)) {
            clusterSamples <- clusters[[j]]
            if (is.element(sample, sample.names[clusterSamples])) {
                sample.clusters[i] <- j
                break;
            }
        }
        if (!is.na(clin$CLI_days_to_death[i])) {
            duration[i] = clin$CLI_days_to_death[i]
        }else if (!is.na(clin$CLI_days_to_last_followup[i])) {
            duration[i] = clin$CLI_days_to_last_followup[i]
        }
    }
    # print(sample.clusters)
    # print(duration)
    samples <- clin$tcga_participant_barcode
    vital_status <- clin$CLI_vital_statu
    surv.clin <- cbind(samples, vital_status, duration, sample.clusters)
    # print(surv.clin)
    surv.clin
}

# Run survival analysis
run.survival.analysis <- function(clin.file,
                                  sample.clusters,
                                  sample.names,
                                  cancer,
                                  exclude.no.clustered.samples = TRUE) {
    require(survival);
    surv.clin <- create.survival.data(clin.file, sample.clusters, sample.names, exclude.no.clustered.samples)
    # print(surv.clin)
    surv <- Surv(surv.clin[, 3], surv.clin[, 2] == 1);
    # print(surv);
    surv.bycluster <- survfit(surv ~ surv.clin[, 4]);
    par(mfrow = c(1, 1));
    # old.par <- par(mar = c(5, 5, 1, 1));
    colors <- c("red", "green", "blue", "cyan", "yellow", "purple");
    plot(surv.bycluster, col=colors, xlab="Survival Time (Day)", ylab="Percent", main = cancer);
    labels <- names(surv.bycluster$strata);
    labels <- paste("Cluster", 1:length(labels), sep=" ");
    print(labels)
    # Removing any prefix in the format of clinPlusClusters$CLUSTER=.
    labels <- sub("((\\w)+\\[, (\\w)+\\.(\\w)+]=)", "", labels);
    legend("topright", labels, lty=1, col=colors, cex=0.5);
    surv.diff <- survdiff(surv ~ surv.clin[, 4]);
    print(surv.diff);
    coxph <- coxph(surv ~ surv.clin[, 4]);
    print(summary(coxph));
    # par(old.par);
}

calculate.sample.dist <- function(file.name, cancer) {
    sample.to.reactions <- read.delim(file.name, header = TRUE, sep = "\t", check.names = FALSE)
    which <- as.character(sample.to.reactions[, 1]) == cancer
    sample.to.reactions <- sample.to.reactions[which, ]
    sample.to.reactions.matrix <- as.matrix(sample.to.reactions[, 3 : length(sample.to.reactions)]);
    rtn <- dist(sample.to.reactions.matrix, method="binary")
    rtn
}


# Plot the cancer types
result.dir <- "/Users/wug/git/Ogmios/results"
sample.to.reaction.file <- paste(result.dir, "MechismoSamplesToReactions_103017.txt", sep = "/")
reaction.dist <- calculate.sample.dist(sample.to.reaction.file, "UCEC")
# print(summary(dist))

# tcga.dist.file <- paste(result.dir, "TCGACancerPairWiseReactionNetworkDist_092617.txt", sep="/")
# tcga.dist.file <- paste(result.dir, "TCGACancerPairWiseReactionNetworkDist_100217.txt", sep="/")
# plot.tcga.cancer(tcga.dist.file)

# Plot for individual samples
# cancer = "TCGA PAAD"
# tcga.sample.dist.file <- paste(result.dir,
#                               "MechismoSamplePairWiseReactionNetworkDist_PAAD_103117.txt",
#                                sep = "/")
# clin.file <- "/Users/wug/datasets/TCGA/clinical/PAAD/PAAD-TP.samplefeatures.txt"

# cancer = "TCGA SKCM"
# tcga.sample.dist.file <- paste(result.dir,
#                                "MechismoSamplePairWiseReactionNetworkDist_SKCM_103117.txt",
#                                sep = "/")
# clin.file <- "/Users/wug/datasets/TCGA/clinical/SKCM/SKCM-TM.samplefeatures.txt"

# cancer = "TCGA LGG"
# tcga.sample.dist.file <- paste(result.dir,
#                                "MechismoSamplePairWiseReactionNetworkDist_LGG_103117.txt",
#                                sep = "/")
# clin.file <- "/Users/wug/datasets/TCGA/clinical/LGG/LGG-TP.samplefeatures.txt"

cancer = "TCGA UCEC"
tcga.sample.dist.file <- paste(result.dir,
                               "MechismoSamplePairWiseReactionNetworkDist_UCEC_103117.txt",
                               sep = "/")
clin.file <- "/Users/wug/datasets/TCGA/clinical/UCEC/UCEC-TP.samplefeatures.txt"

sample.clusters <- plot.tcga.cancer(tcga.sample.dist.file, cancer)
print(sample.clusters)
sample.names <- read.delim(tcga.sample.dist.file, header = TRUE, sep = "\t", check.names = FALSE)[, 1]
sample.names <- as.character(sample.names)
run.survival.analysis(clin.file, sample.clusters, sample.names, cancer, exclude.no.clustered.samples = TRUE)
