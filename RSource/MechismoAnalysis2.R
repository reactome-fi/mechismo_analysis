
library(survival)
library(magrittr)
library(dplyr)
library(stringr)
library(qvalue)
plot.tcga.cancer <- function(tcga.dist.file,
                             main = "TCGA Cancers") {
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
create.survival.data <- function(clin.file,
                                 clusters,
                                 sample.names,
                                 exclude.no.clustered.samples = TRUE) {
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

# Use this function to create a data frame for survival analysis
create.survival.data <- function(clin,
                                 affected.samples) {
  clin2 <- clin %>%
    dplyr::select(tcga_participant_barcode,
                  CLI_days_to_last_followup,
                  CLI_days_to_death,
                  CLI_vital_status) %>%
    dplyr::filter(grepl("TCGA",as.character(tcga_participant_barcode))) %>%
    dplyr::mutate(tcga_participant_barcode = as.character(tcga_participant_barcode),
                  CLI_days_to_last_followup = as.numeric(CLI_days_to_last_followup),
                  CLI_days_to_death = as.numeric(CLI_days_to_death),
                  CLI_vital_status = as.numeric(as.character(CLI_vital_status))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(CLI_vital_status = ifelse(CLI_vital_status == 1,
                                            1,
                                            0)) %>%
    dplyr::mutate(duration = ifelse(CLI_vital_status == 1,
                                    CLI_days_to_death,
                                    CLI_days_to_last_followup)) %>%
    dplyr::mutate(group = ifelse(tcga_participant_barcode %in% affected.samples$Patient,
                                 1,
                                 0))
    
    # dplyr::mutate(group = ifelse(tcga_participant_barcode %in% affected.samples$X1,
    #                              1,
    #                              ifelse(tcga_participant_barcode %in% affected.samples$X2,
    #                                     1,
    #                                     ifelse(tcga_participant_barcode %in% affected.samples$X3,
    #                                            1,
    #                                            ifelse(tcga_participant_barcode %in% affected.samples$X4,
    #                                                   2,
    #                                                   0)))))
  clin2
}

# Run survival analysis
run.survival.analysis <- function(clin,
                                  cancer,
                                  affected.samples) {
  
  surv.clin <- create.survival.data(clin,
                                    affected.samples)
  
  surv <- survival::Surv(time=surv.clin$duration,
                         event=surv.clin$CLI_vital_status,
                         type="right");
  surv.bycluster <- survival::survfit(surv ~ surv.clin$group);
  png(paste("/home/burkhart/Software/Ogmios/results/Mechismo/Figs/",
            cancer,
            ".png",
            sep=""))
  par(mfrow = c(1, 1));
  colors <- c("red", "green", "blue", "cyan", "orange", "purple");
  plot(surv.bycluster,
       col=colors,
       xlab="Survival Time (Days)",
       ylab="Percent",
       main = cancer);
  labels <- names(surv.bycluster$strata);
  labels <- paste("Group",0:(length(labels) - 1), sep=" ");
  legend("topright", labels, lty=1, col=colors, cex=0.5);
  
  surv.diff <- survival::survdiff(surv ~ surv.clin$group);
  coxphresult <- survival::coxph(surv ~ surv.clin$group);
  
  cSum <- summary(coxphresult)
  
  mtext(paste("Score p-value = ",
                   round(cSum$sctest[3],digits=6),
                   sep=""),
        side=3)
  dev.off()
  
  # thisRow <- data.frame(File=cancer,
  #                       logtest=cSum$logtest[3],
  #                       waldtest=cSum$waldtest[3],
  #                       sctest=cSum$sctest[3])
  # all_pvals <<- rbind(all_pvals,thisRow)
  
  print(surv.diff);
  print(summary(coxphresult));
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

result.dir <- "/home/burkhart/Software/Ogmios/results/Mechismo/"
data.dir <- "/home/burkhart/Software/Ogmios/datasets/FirehoseClinical/"
clin.file <- paste(data.dir,"pancancer.samplefeatures.txt",sep="")
cancer = "pancancer"

affected.samples <- read.csv(paste(result.dir,"cluster0UnionDistances.csv",sep=""),
                              header=TRUE,
                              sep=",",
                              check.names=FALSE) %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(Patient = ifelse(is.na(Patient),
                                 NA,
                                 stringr::str_extract(Patient,
                                                      "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")))

clin <- read.csv(clin.file,
                 header = TRUE,
                 sep = "\t",
                 check.names = FALSE)

run.survival.analysis(clin,
                      cancer,
                      affected.samples)

#filenames <- list.files(result.dir,pattern="COADgroup.*\\.csv",full.names = TRUE)
#all_pvals <<- data.frame(File=character(),
                         # logtest=numeric(),
                         # waldtest=numeric(),
                         # sctest=numeric())
#for(f in filenames){
#  affected.samples <- read.csv(f,
#                               header=TRUE,
#                               sep=",",
#                               stringsAsFactors = FALSE,
#                                na.strings = "null") %>%
#     dplyr::rowwise() %>% 
#     dplyr::mutate(X1 = ifelse(is.na(X1),NA,stringr::str_extract(X1,"TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")),
#                   X2 = ifelse(is.na(X2),NA,stringr::str_extract(X2,"TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")),
#                   X3 = ifelse(is.na(X3),NA,stringr::str_extract(X3,"TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")),
#                   X4 = ifelse(is.na(X4),NA,stringr::str_extract(X4,"TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")))
#   
#   cancer2 <- paste(cancer,basename(f),sep=" ")
#   print(paste("Generating survival plot for",cancer2,sep=" "))
#   run.survival.analysis(clin,
#                         cancer2,
#                         affected.samples)
# }
# 
# all_pvals <- cbind(all_pvals,qvalue::qvalue(all_pvals$sctest)$qvalue)
# 
# colnames(all_pvals) <- c("Group",
#                          "CoxPH Log Test pvalue",
#                          "CoxPH Wald Test pvalue",
#                          "CoxPH Score Test pvalue",
#                          "qvalue")
# all_pvals %>%
#   write.table(paste(result.dir,
#                   cancer,
#                   "coxph_survival_significance.tsv",
#                   sep=""),
#             sep="\t",
#             row.names = FALSE)