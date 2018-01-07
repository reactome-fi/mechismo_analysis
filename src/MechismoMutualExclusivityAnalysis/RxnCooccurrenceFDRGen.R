
library(magrittr)
library(dplyr)
library(ggplot2)
library(pander)

panderOptions('table.alignment.default',function(df) ifelse(sapply(df,is.numeric),'right','left'))

REAL_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/INrxnCooccurrence.csv"
RAND_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/random_INrxnCooccurrence.csv"
REW1_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_1INrxnCooccurrence.csv"
REW2_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_2INrxnCooccurrence.csv"
REW3_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_3INrxnCooccurrence.csv"
REW4_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_4INrxnCooccurrence.csv"
REW5_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_5INrxnCooccurrence.csv"
REW6_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_6INrxnCooccurrence.csv"
REW7_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_7INrxnCooccurrence.csv"
REW8_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_8INrxnCooccurrence.csv"
REW9_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_9INrxnCooccurrence.csv"
REW10_RXN_COOCCURR_PATH <- "/home/burkhart/Software/Ogmios/results/Mechismo/rewired_10INrxnCooccurrence.csv"

IMAGES_AS_SVG = TRUE

FIGURES_DIR <- "/home/burkhart/Software/Ogmios/results/Mechismo/Figs/"

realRxnCooccurr <- read.csv(REAL_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::filter(Co.occurrence.BH.Adjusted.P.value < 1) %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value + .Machine$double.xmin),
                rxn = as.factor("real"),
                CooccurrenceDataSource = "Real Data")

randRxnCooccurr <- read.csv(RAND_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::filter(Co.occurrence.BH.Adjusted.P.value < 1) %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value + .Machine$double.xmin),
                rxn = as.factor("random"),
                CooccurrenceDataSource = "Random Data")

rew1RxnCooccurr <- read.csv(REW1_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::filter(Co.occurrence.BH.Adjusted.P.value < 1) %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value + .Machine$double.xmin),
                rxn = as.factor("rew1"),
                CooccurrenceDataSource = "Random Rewiring")


rew2RxnCooccurr <- read.csv(REW2_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::filter(Co.occurrence.BH.Adjusted.P.value < 1) %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value + .Machine$double.xmin),
                rxn = as.factor("rew2"),
                CooccurrenceDataSource = "Random Rewiring")

rew3RxnCooccurr <- read.csv(REW3_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::filter(Co.occurrence.BH.Adjusted.P.value < 1) %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value + .Machine$double.xmin),
                rxn = as.factor("rew3"),
                CooccurrenceDataSource = "Random Rewiring")

rew4RxnCooccurr <- read.csv(REW4_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::filter(Co.occurrence.BH.Adjusted.P.value < 1) %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value + .Machine$double.xmin),
                rxn = as.factor("rew4"),
                CooccurrenceDataSource = "Random Rewiring")

rew5RxnCooccurr <- read.csv(REW5_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::filter(Co.occurrence.BH.Adjusted.P.value < 1) %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value + .Machine$double.xmin),
                rxn = as.factor("rew5"),
                CooccurrenceDataSource = "Random Rewiring")

rew6RxnCooccurr <- read.csv(REW6_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::filter(Co.occurrence.BH.Adjusted.P.value < 1) %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value + .Machine$double.xmin),
                rxn = as.factor("rew6"),
                CooccurrenceDataSource = "Random Rewiring")


rew7RxnCooccurr <- read.csv(REW7_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::filter(Co.occurrence.BH.Adjusted.P.value < 1) %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value + .Machine$double.xmin),
                rxn = as.factor("rew7"),
                CooccurrenceDataSource = "Random Rewiring")

rew8RxnCooccurr <- read.csv(REW8_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::filter(Co.occurrence.BH.Adjusted.P.value < 1) %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value + .Machine$double.xmin),
                rxn = as.factor("rew8"),
                CooccurrenceDataSource = "Random Rewiring")

rew9RxnCooccurr <- read.csv(REW9_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::filter(Co.occurrence.BH.Adjusted.P.value < 1) %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value + .Machine$double.xmin),
                rxn = as.factor("rew9"),
                CooccurrenceDataSource = "Random Rewiring")

rew10RxnCooccurr <- read.csv(REW10_RXN_COOCCURR_PATH,
                              header = TRUE,
                              sep = ",") %>%
  dplyr::filter(Co.occurrence.BH.Adjusted.P.value < 1) %>%
  dplyr::mutate(negLog10FDR = -log10(Co.occurrence.BH.Adjusted.P.value + .Machine$double.xmin),
                rxn = as.factor("rew10"),
                CooccurrenceDataSource = "Random Rewiring")

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
  ggplot(aes(x=negLog10FDR,color=CooccurrenceDataSource)) +
  geom_step(aes(y=1 - (..y..)),
            stat="ecdf") +
  xlab("-log10(FDR)") +
  ylab("relative density") +
  scale_color_manual(values=c("red","black","darkgreen"))
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
                 color=CooccurrenceDataSource)) +
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
  ylab("co-occurring reaction pairs") +
  scale_color_manual(values=c("red","black","darkgreen"))
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
                 color=CooccurrenceDataSource)) +
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
  ylim(low = 0, high = 100) +
  xlab("-log10(FDR)") +
  ylab("co-occurring reaction pairs") +
  scale_color_manual(values=c("red","black","darkgreen"))
if(IMAGES_AS_SVG){
  dev.off()
}

analysisNames <- c("Real",
                   "Random",
                   "Rewire 1",
                   "Rewire 2",
                   "Rewire 3",
                   "Rewire 4",
                   "Rewire 5",
                   "Rewire 6",
                   "Rewire 7",
                   "Rewire 8",
                   "Rewire 9",
                   "Rewire 10")

analysisFDRs <- c(
sum(realRxnCooccurr$negLog10FDR) / (ding * max(realRxnCooccurr$negLog10FDR)),
sum(randRxnCooccurr$negLog10FDR) / (dong * max(randRxnCooccurr$negLog10FDR)),
sum(rew1RxnCooccurr$negLog10FDR) / (thiz * max(rew1RxnCooccurr$negLog10FDR)),
sum(rew2RxnCooccurr$negLog10FDR) / (b2 * max(rew2RxnCooccurr$negLog10FDR)),
sum(rew3RxnCooccurr$negLog10FDR) / (b3 * max(rew3RxnCooccurr$negLog10FDR)),
sum(rew4RxnCooccurr$negLog10FDR) / (b4 * max(rew4RxnCooccurr$negLog10FDR)),
sum(rew5RxnCooccurr$negLog10FDR) / (b5 * max(rew5RxnCooccurr$negLog10FDR)),
sum(rew6RxnCooccurr$negLog10FDR) / (b6 * max(rew6RxnCooccurr$negLog10FDR)),
sum(rew7RxnCooccurr$negLog10FDR) / (b7 * max(rew7RxnCooccurr$negLog10FDR)),
sum(rew8RxnCooccurr$negLog10FDR) / (b8 * max(rew8RxnCooccurr$negLog10FDR)),
sum(rew9RxnCooccurr$negLog10FDR) / (b9 * max(rew9RxnCooccurr$negLog10FDR)),
sum(rew10RxnCooccurr$negLog10FDR) / (b10 * max(rew10RxnCooccurr$negLog10FDR)))

data.frame(Analysis = analysisNames,
           FDR = analysisFDRs) %>%
  pander(split.cells = 50, split.table = Inf)

system_command = paste("qlmanage -t -s 1000 -o ", FIGURES_DIR, " ",FIGURES_DIR,"*.svg",sep="")
system(paste(system_command))
           