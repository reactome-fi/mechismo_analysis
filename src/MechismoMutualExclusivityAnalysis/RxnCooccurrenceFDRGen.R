
library(magrittr)
library(dplyr)
library(ggplot2)

REAL_RXN_COOCCURR_PATH <- "/Users/joshuaburkhart/SoftwareProjects/Ogmios/datasets/rxnCooccurrence.csv"
RAND_RXN_COOCCURR_PATH <- "/Users/joshuaburkhart/SoftwareProjects/Ogmios/datasets/random_rxnCooccurrence.csv"

realRxnCooccurr <- read.table(REAL_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value),
                rxn = as.factor("real"))

randRxnCooccurr <- read.table(RAND_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value),
                rxn = as.factor("random"))

z <- rbind(realRxnCooccurr,randRxnCooccurr)

z %>%
  ggplot(aes(x=negLog10FDR,color=rxn)) +
  geom_step(aes(y=1 - (..y..)),
            stat="ecdf")

nrow(realRxnCooccurr)
nrow(randRxnCooccurr)

z %>% ggplot(aes(x=negLog10FDR,
                 color=rxn)) +
  stat_bin(data=subset(z,as.character(rxn)=="real"),
           aes(y=34368-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="random"),
           aes(y=34002-cumsum(..count..)),
           geom="step")

z %>% ggplot(aes(x=negLog10FDR,
                 color=rxn)) +
  stat_bin(data=subset(z,as.character(rxn)=="real"),
           aes(y=34368-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="random"),
           aes(y=34002-cumsum(..count..)),
           geom="step") +
  ylim(low = 0, high = 3000)
           