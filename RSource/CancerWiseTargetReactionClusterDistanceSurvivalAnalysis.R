







library(survival)
library(magrittr)
library(dplyr)
library(stringr)
library(qvalue)

chisq.pvals <<- numeric()
logrank.pvals <<- numeric()
wald.pvals <<- numeric()
score.pvals <<- numeric()

# Use this function to create a data frame for survival analysis
hinge.group.numeric.survival.data <- function(cancer.dist.column,
                                              cancer.clin.df) {
  quarterNumber <- round(nrow(cancer.dist.column) / 4)
  lower.hinge <- cancer.dist.column %>%
    .[, 3] %>%
    as.numeric() %>%
    sort() %>%
    head(quarterNumber) %>%
    max()
  lower.hinge.patients <-
    cancer.dist.column[cancer.dist.column[, 3] <= lower.hinge, ]
  upper.hinge <- cancer.dist.column %>%
    .[, 3] %>%
    as.numeric() %>%
    sort() %>%
    tail(quarterNumber) %>%
    min()
  upper.hinge.patients <-
    cancer.dist.column[cancer.dist.column[, 3] >= upper.hinge, ]
  return.df <- data.frame
  if (is.matrix(lower.hinge.patients) &
      is.matrix(upper.hinge.patients) &&
      nrow(lower.hinge.patients > 1) &&
      nrow(upper.hinge.patients > 1)) {
    return.df <- cancer.clin.df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(duration = ifelse(
        vital_status == 1,
        as.numeric(days_to_death),
        as.numeric(days_to_last_followup)
      )) %>%
      dplyr::filter(!is.na(duration)) %>%
      dplyr::mutate(
        group = ifelse(
          tcga_participant_barcode %in% lower.hinge.patients[, 1],
          "L",
          "X"
        ),
        group = ifelse(
          tcga_participant_barcode %in% upper.hinge.patients[, 1],
          "U",
          group
        )
      ) %>%
      dplyr::filter(group != "X") %>%
      dplyr::mutate(group = as.factor(group))
  }
  return.df
}



generate.survival.plot <-
  function(cancer.code, surv, group, cancer.type, group.name) {
    surv.fit <- survival::survfit(surv ~ group)# + cancer.type)
    
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
    colors <- rainbow(length(unique(cancer.type)) * 2)
    plot(
      surv.fit,
      col = colors,
      xlab = "Survival Time (Days)",
      ylab = "Percent",
      main = paste(cancer.code, group.name, sep = "")
    )
    #legend("topright",
    #       lty = 1,
    #       col = colors,
    #       cex = 0.5)
    
    surv.diff <- survival::survdiff(surv ~ group)# + cancer.type)
    chisq.pval <- pchisq(surv.diff$chisq, 1, lower.tail = FALSE)
    chisq.pvals <<- c(chisq.pval, chisq.pvals)
    mtext(paste("Chi squared p = ",
                round(chisq.pval, digits = 6),
                sep = ""),
          side = 3)
    dev.off()
    
    print(surv.diff)
  }

result.dir <- "/home/burkhart/Software/Ogmios/results/Mechismo/"
data.dir <-
  "/home/burkhart/Software/Ogmios/datasets/FirehoseClinical/"

dist.file <- paste(result.dir,
                   "BLCATargetReactionGroupPatientDistances.csv",
                   sep = "")

dist.df <- read.csv(dist.file,
                    header = TRUE,
                    sep = ",",
                    check.names = FALSE)

#loop through each column

for (column in 3:ncol(dist.df)) {
  print(paste("Processing column ", column, "...", sep = ""))
  clin.file <-
    paste(data.dir, "BLCA-TP.samplefeatures.txt", sep = "")
  
  clin.df <- read.csv(
    clin.file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  cancer.dist.df <- dist.df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(`Patient Barcode` = ifelse(
      is.na(`Patient Barcode`),
      NA,
      stringr::str_extract(`Patient Barcode`,
                           "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")
    ))
  
  cancer.clin.df <- clin.df %>%
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
                    cancer.dist.df$`Patient Barcode`)
  
  cancer.code <- "BLCA"
  cancer.dist.column <- as.matrix(cancer.dist.df[, c(1,2, column)])
  
    if (max(as.numeric(cancer.dist.column[, 3])) > 1) {
      coxph.surv.clin <- cancer.clin.df %>%
        dplyr::rowwise() %>%
        dplyr::mutate(duration = ifelse(
          vital_status == 1,
          as.numeric(days_to_death),
          as.numeric(days_to_last_followup)
        )) %>%
        dplyr::filter(!is.na(duration)) %>%
        dplyr::filter(tcga_participant_barcode %in% cancer.dist.column[, 1]) %>%
        dplyr::full_join(
          as.data.frame(cancer.dist.column),
          by = c("tcga_participant_barcode" = "Patient Barcode")
        )
      if (length(coxph.surv.clin$duration) > 1 &
          length(coxph.surv.clin$vital_status) > 1) {
        coxph.surv <- survival::Surv(
          time = coxph.surv.clin$duration,
          event = coxph.surv.clin$vital_status,
          type = "right"
        )
        
        coxph.result <-
          survival::coxph(coxph.surv ~
                            as.numeric(as.matrix(coxph.surv.clin[, 7])))# +
                            #as.factor(as.matrix(coxph.surv.clin[, 6]))) 
        
        coxph.summary <- summary(coxph.result)
        
        print(coxph.summary)
        
        logrank.pvals <<- coxph.summary$logtest[3]
        wald.pvals <<- coxph.summary$waldtest[3]
        score.pvals <<- coxph.summary$sctest[3]
      }
      
      cancer.clin.df <- cancer.clin.df %>%
        dplyr::full_join(as.data.frame(cancer.dist.column),
                         by=c("tcga_participant_barcode" = "Patient Barcode"))
      
      surv.clin <-
        hinge.group.numeric.survival.data(cancer.dist.column,
                                          cancer.clin.df)
      if (!is.null(dim(surv.clin))) {
        if (length(surv.clin$duration) > 1 &
            length(surv.clin$vital_status) > 1) {
          surv <- survival::Surv(
            time = surv.clin$duration,
            event = surv.clin$vital_status,
            type = "right"
          )
          
          generate.survival.plot(cancer.code,
                                 surv,
                                 surv.clin$group,
                                 surv.clin$`Cancer Type`,
                                 colnames(cancer.dist.column)[3])
        }
      }
    }
    
    if (max(as.numeric(cancer.dist.column[, 3])) == 1) {
      included.patients <-
        cancer.dist.column[cancer.dist.column[, 3] == 1, ]
      
      cancer.clin.df <- cancer.clin.df %>%
        dplyr::full_join(as.data.frame(cancer.dist.column),
                         by=c("tcga_participant_barcode" = "Patient Barcode"))
      
      surv.clin <- cancer.clin.df %>%
        dplyr::rowwise() %>%
        dplyr::mutate(duration = ifelse(
          vital_status == 1,
          as.numeric(days_to_death),
          as.numeric(days_to_last_followup)
        )) %>%
        dplyr::filter(!is.na(duration)) %>%
        dplyr::mutate(group = ifelse(tcga_participant_barcode %in% included.patients[, 1],
                                     1,
                                     0))
      if (length(surv.clin$duration) > 1 &
          length(surv.clin$vital_status) > 1) {
        surv <- survival::Surv(
          time = surv.clin$duration,
          event = surv.clin$vital_status,
          type = "right"
        )
        
        generate.survival.plot(cancer.code,
                               surv,
                               surv.clin$group,
                               surv.clin$`Cancer Type`,
                               colnames(cancer.dist.column)[3])
      }
    }
  }

bh.chisq.pvals <- p.adjust(chisq.pvals, method = "BH")
bh.logrank.pvals <- p.adjust(logrank.pvals, method = "BH")
bh.wald.pvals <- p.adjust(wald.pvals, method = "BH")
bh.score.pvals <- p.adjust(score.pvals, method = "BH")

pvals.df <- data.frame(
  chisq.p = chisq.pvals,
  chisq.bh = bh.chisq.pvals,
  logrank.p = logrank.pvals,
  logrank.bh = bh.logrank.pvals,
  wald.p = wald.pvals,
  wald.bh = bh.wald.pvals,
  score.p = score.pvals,
  score.bh = bh.score.pvals
)

write.csv(pvals.df,
          file = paste(result.dir, "significance.csv", sep = ""),
          sep = ",")
