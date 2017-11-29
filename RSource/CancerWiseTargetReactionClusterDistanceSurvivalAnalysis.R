





library(survival)
library(magrittr)
library(dplyr)
library(stringr)
library(qvalue)

pvals <<- numeric()

# Use this function to create a data frame for survival analysis
hinge.group.numeric.survival.data <- function(cancer.dist.column,
                                              cancer.clin.df) {
  print("D")
  quarterNumber <- round(nrow(cancer.dist.column) / 3)
  lower.hinge <- cancer.dist.column %>%
    .[, 2] %>%
    as.numeric() %>%
    sort() %>%
    head(quarterNumber) %>%
    max()
  print("D0")
  lower.hinge.patients <-
    cancer.dist.column[cancer.dist.column[, 2] <= lower.hinge, ]
  print("D1")
  upper.hinge <- cancer.dist.column %>%
    .[, 2] %>%
    as.numeric() %>%
    sort() %>%
    tail(quarterNumber) %>%
    min()
  print("D2")
  upper.hinge.patients <-
    cancer.dist.column[cancer.dist.column[, 2] >= upper.hinge, ]
  print("D3")
  print(lower.hinge.patients)
  print("D4")
  print(upper.hinge.patients)
  print("D5")
  print(lower.hinge)
  print(upper.hinge)
  print(class(lower.hinge.patients))
  print(class(upper.hinge.patients))
  print("D6")
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

# Use this function to create a data frame for survival analysis
create.numeric.survival.data <- function(cancer.dist.column,
                                         cancer.clin.df) {
  print("E")
  cancer.clin.df %>%
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
}

# Use this function to create a data frame for survival analysis
logical.survival.data <- function(cancer.dist.column,
                                  cancer.clin.df) {
  print("F")
  included.patients <-
    cancer.dist.column[cancer.dist.column[, 2] == 1, ]
  
  cancer.clin.df %>%
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
}

# Run survival analysis
analyze.survival <- function(cancer.code,
                             cancer.dist.column,
                             cancer.clin.df) {
  z <- nrow(cancer.dist.column)
  print("C")
  print(nrow(cancer.dist.column))
  print(ncol(cancer.dist.column))
  if (z > 0) {
    print("G")
    if (max(as.numeric(cancer.dist.column[, 2])) > 1) {
      print("H")
      coxph.surv.clin <-
        create.numeric.survival.data(cancer.dist.column,
                                     cancer.clin.df)
      if (length(coxph.surv.clin$duration) > 1 &
          length(coxph.surv.clin$vital_status) > 1) {
        print("J")
        coxph.surv <- survival::Surv(
          time = coxph.surv.clin$duration,
          event = coxph.surv.clin$vital_status,
          type = "right"
        )
        
        coxph.result <-
          survival::coxph(coxph.surv ~ as.numeric(as.matrix(coxph.surv.clin[, 6])))
        
        print(summary(coxph.result))
      }
      print("DDD")
      surv.clin <-
        hinge.group.numeric.survival.data(cancer.dist.column,
                                          cancer.clin.df)
      if (!is.null(dim(surv.clin))) {
        print("ZZZ")
        if (length(surv.clin$duration) > 1 &
            length(surv.clin$vital_status) > 1) {
          print("J")
          surv <- survival::Surv(
            time = surv.clin$duration,
            event = surv.clin$vital_status,
            type = "right"
          )
          
          generate.survival.plot(cancer.code,
                                 surv,
                                 surv.clin$group,
                                 colnames(cancer.dist.column)[2])
        }
      }
    }
    
    # if (max(as.numeric(cancer.dist.column[, 2])) == 1) {
    #   print("K")
    #   surv.clin <- logical.survival.data(cancer.dist.column,
    #                                      cancer.clin.df)
    #   if (length(surv.clin$duration) > 1 &
    #       length(surv.clin$vital_status) > 1) {
    #     print("L")
    #     surv <- survival::Surv(
    #       time = surv.clin$duration,
    #       event = surv.clin$vital_status,
    #       type = "right"
    #     )
    #     
    #     generate.survival.plot(cancer.code,
    #                            surv,
    #                            surv.clin$group,
    #                            colnames(cancer.dist.column)[2])
    #   }
    # }
  }
}

generate.survival.plot <-
  function(cancer.code, surv, group, group.name) {
    surv.fit <- survival::survfit(surv ~ group)
    
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
    colors <- c("red", "green")
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
    
    surv.diff <- survival::survdiff(surv ~ group)
    pval <- pchisq(surv.diff$chisq, 0, lower.tail = FALSE)
    pvals <<- c(pval,pvals)
    mtext(paste("Chi squared p = ",
                round(
                  pval, digits = 6
                ),
                sep = ""),
          side = 3)
    dev.off()
    
    print(surv.diff)
  }

result.dir <- "/home/burkhart/Software/Ogmios/results/Mechismo/"
data.dir <-
  "/home/burkhart/Software/Ogmios/datasets/FirehoseClinical/"

dist.file <- paste(
  result.dir,
  "LargestComponentOnly_FDR0p01/PAADTargetReactionGroupPatientDistances.csv",
  sep = ""
)

dist.df <- read.csv(dist.file,
                    header = TRUE,
                    sep = ",",
                    check.names = FALSE)

#loop through each column

for (column in 3:ncol(dist.df)) {
  clin.file <- paste(data.dir, "PAAD-TP.samplefeatures.txt", sep = "")
  
  clin.df <- read.csv(
    clin.file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  max.dist <- max(dist.df[, column])
  
  print(max.dist)
  
  cancer.dist.df <- dist.df %>%
    dplyr::filter(.[[column]] < max.dist) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(`Patient Barcode` = ifelse(
      is.na(`Patient Barcode`),
      NA,
      stringr::str_extract(`Patient Barcode`,
                           "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")
    ))
  
  print("A")
  print(dim(cancer.dist.df))
  
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
  
  print("B")
  print(dim(cancer.dist.df))
  
  analyze.survival("PAAD",
                   as.matrix(cancer.dist.df[, c(1, column)]),
                   cancer.clin.df)
}

bh.pvals <- p.adjust(pvals,method="BH")
print(bh.pvals)
