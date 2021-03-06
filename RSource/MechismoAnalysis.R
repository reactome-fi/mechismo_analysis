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

output.clusters <- function(clusters, names) {
    cat(paste("Total clusters:", length(clusters), "\n", sep = " "))
    for (i in 1 : length(clusters)) {
        cluster <- clusters[[i]]
        cat(paste("Cluster", i, ":", length(cluster), "\n", sep=" "))
        for (j in 1 : length(cluster)) {
            cat(paste(names[cluster[j]], "\t", sep=""))
        }
        cat("\n")
    }
}

analyze.mechiso.fi <- function(file, use.reactome.only = FALSE) {
    mechismo.fi <- read.delim(file, sep = "\t", header = T)
    # filter to use Reactome FI only
    # if (use.reactome.only) {
    #     which <- 
    # }
    # Caclualte pair-wise correlation
    names <- c('SignificantCancers', 'RankScore', 'PanCancer_PScore')
    par(mfrow = c(2, 2));
    for (i in 1 : (length(names) - 1)) {
        for (j in (i + 1) : length(names)) {
            print(paste(names[i], "~", names[j], sep = " "))
            cor.result <- cor.test(mechismo.fi[, names[i]], mechismo.fi[, names[j]])
            print(cor.result)
            plot(mechismo.fi[, names[i]], mechismo.fi[, names[j]],
                 xlab = names[i],
                 ylab = names[j])
        }
    }
    # Try to merge three dimentional information into 2D using categories
    total.fis <- dim(mechismo.fi)[1]
    disc.sig.cancer <- rep(0, total.fis)
    for (i in 1 : total.fis) {
        counter <- mechismo.fi[i, 'SignificantCancers']
        if (counter > 5) {
            disc.sig.cancer[i] <- 6
        }else {
            disc.sig.cancer[i] <- counter
        }
    }
    print("Counters of SignificantCancers:")
    print(with(mechismo.fi, table(SignificantCancers)))
    print("Counters after discretizing SignificantCancers:")
    print(table(disc.sig.cancer))
    with(mechismo.fi, plot(RankScore, PanCancer_PScore, pch = disc.sig.cancer))
    disc.sig.cancer.sort <- sort(unique(disc.sig.cancer))
    legend('bottomright', c('0','1','2','3','4','5','>6'), pch = disc.sig.cancer.sort)
}

plot.mechiso.fi <- function(file) {
    results <- read.delim(file, sep="\t", header = T)
    par(mfrow = c(2, 2))
    main <- "Histogram of Significant Interactions"
    hist(results$SignificantCancers, breaks = 20, main = main, xlab = "Significant Cancers")
    hist(results$SignificantCancers, breaks = 20, main = main, ylim = c(1, 100), xlab = "Significant Cancers")
    # Filter to Reactome FIs only
    which <- results$InReactome == 1
    reactome.results <- results[which, ]
    main <- "Histogram of Significant Reactome FIs"
    hist(reactome.results$SignificantCancers, main = main, breaks = 20, xlab = "Significant Cancers")
    hist(reactome.results$SignificantCancers, main = main, breaks = 20, ylim = c(1, 100), xlab = "Significant Cancers")
}

# Plot the cancer types
result.dir <- "/Users/wug/git/Ogmios/results"
# mechismo.fi.file <- "FI_Cancer_FDR_Filtered_051018.txt";
mechismo.fi.file <- "FI_Cancer_FDR_Filtered_082118.txt";
analyze.mechiso.fi(paste(result.dir, mechismo.fi.file, sep = "/"))
stop("Done!")

sample.to.reaction.file <- paste(result.dir, "MechismoSamplesToReactions_103017.txt", sep = "/")
reaction.dist <- calculate.sample.dist(sample.to.reaction.file, "UCEC")
# print(summary(dist))

# tcga.dist.file <- paste(result.dir, "TCGACancerPairWiseReactionNetworkDist_092617.txt", sep="/")
# tcga.dist.file <- paste(result.dir, "TCGACancerPairWiseReactionNetworkDist_100217.txt", sep="/")
# plot.tcga.cancer(tcga.dist.file)

# Plot for individual samples
# cancer = "TCGA PAAD"
# tcga.sample.dist.file <- paste(result.dir,
#                               "MechismoSamplePairWiseReactionNetworkDist_PAAD_111617.txt",
#                                sep = "/")
# clin.file <- "/Users/wug/datasets/TCGA/clinical/PAAD/PAAD-TP.samplefeatures.txt"

cancer = "TCGA SKCM"
tcga.sample.dist.file <- paste(result.dir,
                               "MechismoSamplePairWiseReactionNetworkDist_SKCM_112817.txt",
                               sep = "/")
clin.file <- "/Users/wug/datasets/TCGA/clinical/SKCM/SKCM-TM.samplefeatures.txt"

cancer = "TCGA LGG"
tcga.sample.dist.file <- paste(result.dir,
                               "MechismoSamplePairWiseReactionNetworkDist_LGG_111617.txt",
                               sep = "/")
clin.file <- "/Users/wug/datasets/TCGA/clinical/LGG/LGG-TP.samplefeatures.txt"

cancer = "TCGA UCEC"
tcga.sample.dist.file <- paste(result.dir,
                               "MechismoSamplePairWiseReactionNetworkDist_UCEC_111317.txt",
                               sep = "/")
clin.file <- "/Users/wug/datasets/TCGA/clinical/UCEC/UCEC-TP.samplefeatures.txt"
# 
# cancer = "TCGA COADREAD"
# tcga.sample.dist.file <- paste(result.dir,
#                                "MechismoSamplePairWiseReactionNetworkDist_COADREAD_112817.txt",
#                                sep = "/")
# clin.file <- "/Users/wug/datasets/TCGA/clinical/COADREAD/COADREAD-TP.samplefeatures.txt"
# 
# cancer = "TCGA THCA"
# tcga.sample.dist.file <- paste(result.dir,
#                                "MechismoSamplePairWiseReactionNetworkDist_THCA_111617.txt",
#                                sep = "/")
# clin.file <- "/Users/wug/datasets/TCGA/clinical/THCA/THCA-TP.samplefeatures.txt"
# 
# cancer = "TCGA STAD"
# tcga.sample.dist.file <- paste(result.dir,
#                                "MechismoSamplePairWiseReactionNetworkDist_STAD_111617.txt",
#                                sep = "/")
# clin.file <- "/Users/wug/datasets/TCGA/clinical/STAD/STAD-TP.samplefeatures.txt"
# 
# cancer = "TCGA LUAD"
# tcga.sample.dist.file <- paste(result.dir,
#                                "MechismoSamplePairWiseReactionNetworkDist_LUAD_111617.txt",
#                                sep = "/")
# clin.file <- "/Users/wug/datasets/TCGA/clinical/LUAD/LUAD-TP.samplefeatures.txt"
# 
cancer = "TCGA HNSC"
tcga.sample.dist.file <- paste(result.dir,
                               "MechismoSamplePairWiseReactionNetworkDist_HNSC_111617.txt",
                               sep = "/")
clin.file <- "/Users/wug/datasets/TCGA/clinical/HNSC/HNSC-TP.samplefeatures.txt"

sample.clusters <- plot.tcga.cancer(tcga.sample.dist.file, cancer)
# print(sample.clusters)
sample.names <- read.delim(tcga.sample.dist.file, header = TRUE, sep = "\t", check.names = FALSE)[, 1]
sample.names <- as.character(sample.names)
output.clusters(sample.clusters, sample.names)
cat("\n") # Give it an extra line
run.survival.analysis(clin.file, sample.clusters, sample.names, cancer, exclude.no.clustered.samples = T)
