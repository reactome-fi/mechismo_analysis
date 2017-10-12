
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
            stat="ecdf") +
  xlab("-log10(FDR)") +
  ylab("relative density")

# https://stackoverflow.com/questions/17832608/how-to-use-earlier-declared-variables-within-aes-in-ggplot-with-special-operator
ding <- nrow(realRxnCooccurr)
dong <- nrow(randRxnCooccurr)

# https://stackoverflow.com/questions/18379933/plotting-cumulative-counts-in-ggplot2
z %>% ggplot(aes(x=negLog10FDR,
                 color=rxn)) +
  stat_bin(data=subset(z,as.character(rxn)=="real"),
           aes(ding = ding, y=ding-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="random"),
           aes(dong = dong,y=dong-cumsum(..count..)),
           geom="step") +
  xlab("-log10(FDR)") +
  ylab("co-occurring reaction pairs")

z %>% ggplot(aes(x=negLog10FDR,
                 color=rxn)) +
  stat_bin(data=subset(z,as.character(rxn)=="real"),
           aes(ding = ding, y=ding-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="random"),
           aes(dong = dong, y=dong-cumsum(..count..)),
           geom="step") +
  ylim(low = 0, high = 3000) +
  xlab("-log10(FDR)") +
  ylab("co-occurring reaction pairs")
           