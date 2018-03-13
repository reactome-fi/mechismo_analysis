#! /bin/bash

wget \
    -r \
    -nd \
    -P ~/Software/Ogmios/datasets/FirehoseClinical/ \
    --no-parent \
    -A '*.samplefeatures.txt' \
    http://gdac.broadinstitute.org/runs/analyses__latest/reports/cancer/{COADREAD,LAML,ACC,BLCA,LGG,BRCA,CESC,CHOL,LCML,COAD,CNTL,ESCA,FPPP,GBM,HNSC,KICH,KIRC,KIRP,LIHC,LUAD,LUSC,DLBC,MESO,MISC,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THYM,THCA,UCS,ACC,BLCA,LGG,BRCA,CESC,CHOL,LCML,COAD,CNTL,ESCA,FPPP,GBM,HNSC,KICH,KIRC,KIRP,LIHC,LUAD,LUSC,DLBC,MESO,MISC,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THYM,THCA,UCS,UCEC,UVM,UCEC,UVM}/Aggregate_AnalysisFeatures/
