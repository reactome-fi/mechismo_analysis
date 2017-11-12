
library(survival)
library(magrittr)
library(dplyr)
library(stringr)
library(qvalue)

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
    dplyr::mutate(group = ifelse(tcga_participant_barcode %in% affected.samples$X1,
                                 1,
                                 ifelse(tcga_participant_barcode %in% affected.samples$X2,
                                        1,
                                        ifelse(tcga_participant_barcode %in% affected.samples$X3,
                                               1,
                                               ifelse(tcga_participant_barcode %in% affected.samples$X4,
                                                      2,
                                                      0)))))
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
  
  thisRow <- data.frame(File=cancer,
                        logtest=cSum$logtest[3],
                        waldtest=cSum$waldtest[3],
                        sctest=cSum$sctest[3])
  all_pvals <<- rbind(all_pvals,thisRow)
  
  print(surv.diff);
  print(summary(coxphresult));
}

result.dir <- "/home/burkhart/Software/Ogmios/results/Mechismo/"
data.dir <- "/home/burkhart/Software/Ogmios/datasets/FirehoseClinical/"
clin.file <- paste(data.dir,"THCA-TP.samplefeatures.txt",sep="")
cancer = "THCA"
clin <- read.csv(clin.file,
                 header = TRUE,
                 sep = "\t",
                 check.names = FALSE)

filenames <- list.files(result.dir,pattern="THCAgroup.*\\.csv",full.names = TRUE)
all_pvals <<- data.frame(File=character(),
                         logtest=numeric(),
                         waldtest=numeric(),
                         sctest=numeric())
for(f in filenames){
  affected.samples <- read.csv(f,
                               header=TRUE,
                               sep=",",
                               stringsAsFactors = FALSE,
                               na.strings = "null") %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(X1 = ifelse(is.na(X1),NA,stringr::str_extract(X1,"TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")),
                  X2 = ifelse(is.na(X2),NA,stringr::str_extract(X2,"TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")),
                  X3 = ifelse(is.na(X3),NA,stringr::str_extract(X3,"TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")),
                  X4 = ifelse(is.na(X4),NA,stringr::str_extract(X4,"TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")))
  
  cancer2 <- paste(cancer,basename(f),sep=" ")
  print(paste("Generating survival plot for",cancer2,sep=" "))
  run.survival.analysis(clin,
                        cancer2,
                        affected.samples)
}

all_pvals <- cbind(all_pvals,qvalue::qvalue(all_pvals$sctest)$qvalue)

colnames(all_pvals) <- c("Group",
                         "CoxPH Log Test pvalue",
                         "CoxPH Wald Test pvalue",
                         "CoxPH Score Test pvalue",
                         "qvalue")
all_pvals %>%
  write.table(paste(result.dir,
                  cancer,
                  "coxph_survival_significance.tsv",
                  sep=""),
            sep="\t",
            row.names = FALSE)






