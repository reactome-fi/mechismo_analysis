library(knitr)
library(markdown)
library(rmarkdown)

input_csvs <- list.files("/media/burkhart/Media/Software/Ogmios/results/heatmapData/",pattern = "*interactions_with_mutated_interfaces.csv")
for(input_csv in input_csvs){
  input_csv_path <- paste("/media/burkhart/Media/Software/Ogmios/results/heatmapData/",input_csv,sep="")
  rmarkdown::render("/media/burkhart/Media/Software/Ogmios/src/MechismoMutualExclusivityAnalysis/FirehoseInteractionHeatmaps.Rmd",
                    output_format = "pdf_document",
                    output_file = paste("FirehoseInteractionHeatmaps_",gsub("firehose_([A-Z]+).*","\\1",input_csv),".pdf",sep=""),
                    output_dir = "/media/burkhart/Media/Software/Ogmios/src/MechismoMutualExclusivityAnalysis/")
}
