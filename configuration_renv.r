setwd("/home/benjamin/PROJETS/5_Developpement/RNA-Seq_Analyses-0322/")

renv::init()

renv::install('readr')
renv::install('reshape')
renv::install('ggplot2')
renv::install('RColorBrewer')
renv::install('pheatmap')
renv::install('ggrepel')
renv::install('rmarkdown')
BiocManager::install(c("DESeq2", "tximport","limma","biomaRt"))

renv::snapshot()
