

library(survival)
library(magrittr)
library(dplyr)
library(stringr)
library(qvalue)

P_VAL_THRESH <- 0.05

# Use this function to create a data frame for survival analysis
create.survival.data <- function(cancer.dist.df,
                                 cancer.clin.df) {
  
  lt.median.1 <- cancer.dist.df %>%
    dplyr::filter(`Distance to 1` < median(cancer.dist.df$`Distance to 1`))
  lt.median.2 <- cancer.dist.df %>%
    dplyr::filter(`Distance to 2` < median(cancer.dist.df$`Distance to 2`))
  lt.median.3 <- cancer.dist.df %>%
  dplyr::filter(`Distance to 3` < median(cancer.dist.df$`Distance to 3`))
  lt.median.T <- cancer.dist.df %>%
    dplyr::filter(`Distance to T` < median(cancer.dist.df$`Distance to T`))
  lt.median.All <- cancer.dist.df %>%
    dplyr::filter(`Distance to All` < median(cancer.dist.df$`Distance to All`))
  
  cancer.clin.df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(duration = ifelse(
      vital_status == 1,
      as.numeric(days_to_death),
      as.numeric(days_to_last_followup)
    )) %>%
    dplyr::filter(!is.na(duration)) %>%
    dplyr::mutate(
      group1 = ifelse(tcga_participant_barcode %in% lt.median.1$Patient,
                      1,
                      0),
      group2 = ifelse(tcga_participant_barcode %in% lt.median.2$Patient,
                      1,
                      0),
      group3 = ifelse(tcga_participant_barcode %in% lt.median.3$Patient,
                      1,
                      0),
      groupT = ifelse(tcga_participant_barcode %in% lt.median.T$Patient,
                      1,
                      0),
      groupAll = ifelse(tcga_participant_barcode %in% lt.median.All$Patient,
                        1,
                        0)
    )
}

# Run survival analysis
analyze.survival <- function(cancer.code,
                             cancer.dist.df,
                             cancer.clin.df) {
  
  surv.clin <- create.survival.data(cancer.dist.df,
                                    cancer.clin.df)
  
  if(length(surv.clin$duration) > 1 &
     length(surv.clin$vital_status) > 1){
    
  surv <- survival::Surv(
    time = surv.clin$duration,
    event = surv.clin$vital_status,
    type = "right"
  )
  
  generate.survival.plot(cancer.code, surv, surv.clin$group1, "1 upstream rxn")
  generate.survival.plot(cancer.code, surv, surv.clin$group2, "2 upstream rxns")
  generate.survival.plot(cancer.code, surv, surv.clin$group3, "3+ upstream rxns")
  generate.survival.plot(cancer.code, surv, surv.clin$groupT, "target rxn")
  generate.survival.plot(cancer.code,
                         surv,
                         surv.clin$groupAll,
                         "all target-related rxns")
  }
}

generate.survival.plot <- function(cancer.code, surv, group, group.name) {
  
  surv.fit <- survival::survfit(surv ~ group)
  
  coxph.result <- survival::coxph(surv ~ group)
  
  coxph.summary <- summary(coxph.result)
  
  coxph.sctest.p <- coxph.summary$sctest[3]
  
  if(coxph.sctest.p < P_VAL_THRESH){
  png(
    paste(
      "/home/burkhart/Software/Ogmios/results/Mechismo/Figs/",
      cancer.code,
      group.name,
      ".png",
      sep = ""
    )
  )
  par(mfrow = c(1, 1))
  colors <- c("red", "green", "blue", "cyan", "orange", "purple")
  plot(
    surv.fit,
    col = colors,
    xlab = "Survival Time (Days)",
    ylab = "Percent",
    main = paste(cancer.code,group.name,sep="")
  )
  labels <- c(paste("Patients far from from ",group.name,sep=""),
              paste("Patients close to from ",group.name,sep=""))
  legend("topright",
         labels,
         lty = 1,
         col = colors,
         cex = 0.5)
  
  surv.diff <- survival::survdiff(surv ~ group)
  
  mtext(paste(
    "Score p-value = ",
    round(coxph.sctest.p, digits = 6),
    sep = ""
  ),
  side = 3)
  dev.off()
  
  print(surv.diff)
  print(coxph.summary)
  }
}

cancer.codes <- list("COAD",
                     "LGG",
                     "BRCA",
                     "PAAD",
                     "LUAD",
                     "LUSC",
                     "GBM",
                     "OV",
                     "KIRC",
                     "KIRP",
                     "PRAD",
                     "SARC",
                     "STAD",
                     "THCA",
                     "THYM",
                     "UVM",
                     "ACC",
                     "BLCA",
                     "HNSC",
                     "UCEC")
result.dir <- "/home/burkhart/Software/Ogmios/results/Mechismo/"
data.dir <-
  "/home/burkhart/Software/Ogmios/datasets/FirehoseClinical/"

dist.file <- paste(result.dir, "patientDistances.csv", sep = "")

dist.df <- read.csv(dist.file,
                    header = TRUE,
                    sep = ",",
                    check.names = FALSE)

max.dist <- max(dist.df$`Distance to All`)

for (cancer in 1:length(cancer.codes)) {
  print(paste("Analyzing ",cancer.codes[[cancer]],"...",sep=""))
  
clin.file <- paste(data.dir,cancer.codes[[cancer]],"-TP.samplefeatures.txt", sep = "")

clin.df <- read.csv(
  clin.file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE)
  
  cancer.dist.df <- dist.df %>%
    dplyr::filter(grepl(cancer.codes[[cancer]], `Cancer Type`),
                  `Distance to All` < max.dist) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Patient = ifelse(
      is.na(Patient),
      NA,
      stringr::str_extract(Patient,
                           "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")))
  
  cancer.clin.df <- clin.df %>%
    dplyr::mutate(
      days_to_last_followup = as.numeric(CLI_days_to_last_followup),
      days_to_death = as.numeric(CLI_days_to_death),
      vital_status = as.numeric(CLI_vital_status)
    ) %>%
    dplyr::select(tcga_participant_barcode,
                  days_to_last_followup,
                  days_to_death,
                  vital_status)
  
  analyze.survival(cancer.codes[[cancer]],
                   cancer.dist.df,
                   cancer.clin.df)
}
