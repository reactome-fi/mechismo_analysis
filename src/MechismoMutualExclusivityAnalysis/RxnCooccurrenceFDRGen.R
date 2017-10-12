
library(magrittr)
library(dplyr)
library(ggplot2)

REAL_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rxnCooccurrence.csv"
RAND_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/random_rxnCooccurrence.csv"
REW1_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_1rxnCooccurrence.csv"
REW2_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_2rxnCooccurrence.csv"
REW3_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_3rxnCooccurrence.csv"
REW4_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_4rxnCooccurrence.csv"
REW5_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_5rxnCooccurrence.csv"
REW6_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_6rxnCooccurrence.csv"
REW7_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_7rxnCooccurrence.csv"
REW8_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_8rxnCooccurrence.csv"
REW9_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_9rxnCooccurrence.csv"
REW10_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_10rxnCooccurrence.csv"

IMAGES_AS_SVG = TRUE

FIGURES_DIR <- "/home/burkhart/Software/Ogmios/results/Mechismo/Figs/"

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

rew1RxnCooccurr <- read.table(REW1_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value),
                rxn = as.factor("rew1"))


rew2RxnCooccurr <- read.table(REW2_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value),
                rxn = as.factor("rew2"))

rew3RxnCooccurr <- read.table(REW3_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value),
                rxn = as.factor("rew3"))

rew4RxnCooccurr <- read.table(REW4_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value),
                rxn = as.factor("rew4"))

rew5RxnCooccurr <- read.table(REW5_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value),
                rxn = as.factor("rew5"))

rew6RxnCooccurr <- read.table(REW6_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value),
                rxn = as.factor("rew6"))


rew7RxnCooccurr <- read.table(REW7_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value),
                rxn = as.factor("rew7"))

rew8RxnCooccurr <- read.table(REW8_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value),
                rxn = as.factor("rew8"))

rew9RxnCooccurr <- read.table(REW9_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value),
                rxn = as.factor("rew9"))

rew10RxnCooccurr <- read.table(REW10_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value),
                rxn = as.factor("rew10"))

z <- rbind(realRxnCooccurr,
           randRxnCooccurr,
           rew1RxnCooccurr,
           rew2RxnCooccurr,
           rew3RxnCooccurr,
           rew4RxnCooccurr,
           rew5RxnCooccurr,
           rew6RxnCooccurr,
           rew7RxnCooccurr,
           rew8RxnCooccurr,
           rew9RxnCooccurr,
           rew10RxnCooccurr
           )

if(IMAGES_AS_SVG){
  svg(filename=paste(FIGURES_DIR,"relative_density.svg",sep=""),
      width=8,
      height=8,
      pointsize=20)
}
z %>%
  ggplot(aes(x=negLog10FDR,color=rxn)) +
  geom_step(aes(y=1 - (..y..)),
            stat="ecdf") +
  xlab("-log10(FDR)") +
  ylab("relative density")
if(IMAGES_AS_SVG){
  dev.off()
}

# https://stackoverflow.com/questions/17832608/how-to-use-earlier-declared-variables-within-aes-in-ggplot-with-special-operator
ding <- nrow(realRxnCooccurr)
dong <- nrow(randRxnCooccurr)
thiz <- nrow(rew1RxnCooccurr)
b2 <- nrow(rew2RxnCooccurr)
b3 <- nrow(rew3RxnCooccurr)
b4 <- nrow(rew4RxnCooccurr)
b5 <- nrow(rew5RxnCooccurr)
b6 <- nrow(rew6RxnCooccurr)
b7 <- nrow(rew7RxnCooccurr)
b8 <- nrow(rew8RxnCooccurr)
b9 <- nrow(rew9RxnCooccurr)
b10 <- nrow(rew10RxnCooccurr)

# https://stackoverflow.com/questions/18379933/plotting-cumulative-counts-in-ggplot2
if(IMAGES_AS_SVG){
  svg(filename=paste(FIGURES_DIR,"rxnCooccurrence.svg",sep=""),
      width=8,
      height=8,
      pointsize=20)
}
z %>% ggplot(aes(x=negLog10FDR,
                 color=rxn)) +
  stat_bin(data=subset(z,as.character(rxn)=="real"),
           aes(ding = ding, y=ding-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="random"),
           aes(dong = dong,y=dong-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew1"),
           aes(thiz = thiz,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew2"),
           aes(thiz = b2,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew3"),
           aes(thiz = b3,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew4"),
           aes(thiz = b4,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew5"),
           aes(thiz = b5,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew6"),
           aes(thiz = b6,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew7"),
           aes(thiz = b7,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew8"),
           aes(thiz = b8,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew9"),
           aes(thiz = b9,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew10"),
           aes(thiz = b10,y=thiz-cumsum(..count..)),
           geom="step") +
  xlab("-log10(FDR)") +
  ylab("co-occurring reaction pairs")
if(IMAGES_AS_SVG){
  dev.off()
}

if(IMAGES_AS_SVG){
  svg(filename=paste(FIGURES_DIR,"rxnCooccurrenceZoom.svg",sep=""),
      width=8,
      height=8,
      pointsize=20)
}
z %>% ggplot(aes(x=negLog10FDR,
                 color=rxn)) +
  stat_bin(data=subset(z,as.character(rxn)=="real"),
           aes(ding = ding, y=ding-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="random"),
           aes(dong = dong,y=dong-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew1"),
           aes(thiz = thiz,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew2"),
           aes(thiz = b2,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew3"),
           aes(thiz = b3,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew4"),
           aes(thiz = b4,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew5"),
           aes(thiz = b5,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew6"),
           aes(thiz = b6,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew7"),
           aes(thiz = b7,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew8"),
           aes(thiz = b8,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew9"),
           aes(thiz = b9,y=thiz-cumsum(..count..)),
           geom="step") +
  stat_bin(data=subset(z,as.character(rxn)=="rew10"),
           aes(thiz = b10,y=thiz-cumsum(..count..)),
           geom="step") +
  ylim(low = 0, high = 2200) +
  #xlim(low = 0, high = 125) +
  xlab("-log10(FDR)") +
  ylab("co-occurring reaction pairs")
if(IMAGES_AS_SVG){
  dev.off()
}

sum(realRxnCooccurr$negLog10FDR) / (ding * max(realRxnCooccurr$negLog10FDR))
sum(randRxnCooccurr$negLog10FDR) / (ding * max(realRxnCooccurr$negLog10FDR))
sum(rew1RxnCooccurr$negLog10FDR) / (ding * max(realRxnCooccurr$negLog10FDR))
sum(rew2RxnCooccurr$negLog10FDR) / (ding * max(realRxnCooccurr$negLog10FDR))
sum(rew3RxnCooccurr$negLog10FDR) / (ding * max(realRxnCooccurr$negLog10FDR))
sum(rew4RxnCooccurr$negLog10FDR) / (ding * max(realRxnCooccurr$negLog10FDR))
sum(rew5RxnCooccurr$negLog10FDR) / (ding * max(realRxnCooccurr$negLog10FDR))
sum(rew6RxnCooccurr$negLog10FDR) / (ding * max(realRxnCooccurr$negLog10FDR))
sum(rew7RxnCooccurr$negLog10FDR) / (ding * max(realRxnCooccurr$negLog10FDR))
sum(rew8RxnCooccurr$negLog10FDR) / (ding * max(realRxnCooccurr$negLog10FDR))
sum(rew9RxnCooccurr$negLog10FDR) / (ding * max(realRxnCooccurr$negLog10FDR))
sum(rew10RxnCooccurr$negLog10FDR) / (ding * max(realRxnCooccurr$negLog10FDR))

system_command = paste("qlmanage -t -s 1000 -o ", FIGURES_DIR, " ",FIGURES_DIR,"*.svg",sep="")
system(paste(system_command))
           