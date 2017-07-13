library(knitr)
library(markdown)
library(rmarkdown)

cancer_types <- list.files("/media/burkhart/Media/Software/Ogmios/datasets/Mechismo/cancer_types/")
for(cancer_type in cancer_types){
  input_txt_path <- paste("/media/burkhart/Media/Software/Ogmios/datasets/Mechismo/cancer_types/",
                          cancer_type,
                          "/sample_id.txt",
                          sep="")
  if(length(readLines(input_txt_path)) > 1000){
    rmarkdown::render("/media/burkhart/Media/Software/Ogmios/src/MechismoMutualExclusivityAnalysis/MechismoHeatmaps.Rmd",
                      output_format = "pdf_document",
                      output_file = paste("MechismoHeatmaps_",cancer_type,".pdf",sep=""),
                      output_dir = "/media/burkhart/Media/Software/Ogmios/src/MechismoMutualExclusivityAnalysis/MechismoOutput/")
  }
}