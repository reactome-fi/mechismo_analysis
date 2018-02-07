# libs
library(survival)
library(magrittr)
library(dplyr)

# global vars
K=4
CANCER_CODE <- "PAAD"
#DATA_DIR <- "/home/burkhart/Software/Ogmios/results/Mechismo/"
DATA_DIR <- "/Users/joshuaburkhart/SoftwareProjects/Ogmios/results/Mechismo/"
CLIN_DIR <- "/Users/joshuaburkhart/SoftwareProjects/Ogmios/datasets/FirehoseClinical/"
FIG_DIR <- paste(DATA_DIR,"Figs/",sep="")
DIST_FILE <- paste(DATA_DIR,
                   "MechismoSamplePairwiseFINetworkDist_",
                  CANCER_CODE,
                  ".tsv",
                   sep="")
CLIN_FILE <- paste(CLIN_DIR,
                   CANCER_CODE,
                   "-TP.samplefeatures.txt",
                   sep="")

######################
## cluster patients ##
######################

# load patient distance data
dist_df <- read.csv(DIST_FILE,
                    sep="\t",
                    header = TRUE,
                    row.names = 1)

# cluster
h_clustering <- dist(dist_df) %>%
  hclust(method="ward.D")

# plot
svg(filename=paste(FIG_DIR,CANCER_CODE,"_hclust.svg",sep=""),
    width=8,
    height=8,
    pointsize=10)
plot(h_clustering,
     main=paste(CANCER_CODE," Patient FI Distances",sep=""),
     ylab = "Patients",
     xlab = "FI Distances")
dev.off()

# extract top K clusters
top_groups <- cutree(h_clustering,
                         k=K)

# extract common patient barcode substring from top_groups
# transform like "TCGA-AY-4070-01" -> "TCGA-AY-4070"
names(top_groups) <- names(top_groups) %>%
    stringr::str_extract("TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")

###############################
## perform survival analysis ##
###############################

# load clinical data
clin_df <- read.csv(
  CLIN_FILE,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# extract survival info from clinical data
clin_df <- clin_df %>%
  dplyr::mutate(
    days_to_last_followup = as.numeric(CLI_days_to_last_followup),
    days_to_death = as.numeric(CLI_days_to_death),
    vital_status = as.numeric(CLI_vital_status)
  ) %>%
  dplyr::select(tcga_participant_barcode,
                days_to_last_followup,
                days_to_death,
                vital_status) %>%
  dplyr::filter(tcga_participant_barcode %in%
                  names(top_groups))

top_groups <- top_groups[clin_df$tcga_participant_barcode]

# calculate duration
surv_df <- clin_df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(duration = ifelse(
    vital_status == 1,
    as.numeric(days_to_death),
    as.numeric(days_to_last_followup)
  ))

# calculate Cox PH
coxph_surv <- survival::Surv(
  time = surv_df$duration,
  event = surv_df$vital_status,
  type = "right"
)

# report Cox PH summary
coxph_summary <-
  survival::coxph(coxph_surv ~ surv_df$duration) %>%
  summary()
                 
print(coxph_summary)

logrank.pvals <- coxph_summary$logtest[3]
wald.pvals <- coxph_summary$waldtest[3]
score.pvals <- coxph_summary$sctest[3]

# calculate survival
surv.fit <- survival::survfit(coxph_surv ~ top_groups)

# plot
par(mfrow = c(1, 1))
svg(filename=paste(FIG_DIR,CANCER_CODE,"_top_",K,"_survival.svg",sep=""),
    width=8,
    height=8,
    pointsize=10)
plot(
  surv.fit,
  col = 3:6,
  lty = 2:5,
  lwd = 2,
  xlab = "Survival Time (Days)",
  ylab = "Ratio",
  main = paste("Top ",
               K,
               " ",
               CANCER_CODE,
               " Patient Clusters",
               sep = "")
)
legend(
  "topright",
  legend=names(surv.fit$strata),
  col=3:6,
  lty=2:5,
  lwd = 2,
  horiz=FALSE,
  bty='n')
surv.diff <- survival::survdiff(coxph_surv ~ top_groups)
chisq.pval <- pchisq(surv.diff$chisq, 3, lower.tail = FALSE)
mtext(paste("Chi squared p = ",
            round(chisq.pval, digits = 6),
            sep = ""),
      side = 3)
print(surv.diff)
dev.off()
