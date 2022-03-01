# RNA-Seq analysis with DESeq2

## Setup

To setup the pipeline, you need to restore the R environment. Start by cloning this repository, then open RStudio and set the working directory. Run the `renv::restore` command to verify the absence of errors. If there is none, the R environment is now set.
If you have any errors, please follow the R instruction to update or install the recommended packages.

You certainly shouldn't move or modify any file in the git directory before the install.

### The configuration file

The script is based on two tabulated configuration files, as PROJET.conf & PROJET.contrastes. The first one associate the RSEM count file for each sample to sample name and the different groups it belongs to. The second one lists all the contrast DESeq2 will perform.

PROJET.conf
```
File	Name	Condition	Preparation
/path/to/te/RSEM/file sample_name Group Batch
```

You can add many conditions or batch in this file, but you will need to configure the RNA-Seq_pipeline.2112.Rmd in consequences.

PROJET.contrast
```
Condition Patients  Controls
```
As the script uses the `contrast()` option in the `DESeq2::results` function, this file will be parse and the script performs statistical analysis for each line. DESeq2 offers many options for the test, so if you do not want the default parameters, you will need to adapt the `DESeq2::results` command line.

### The script

All the pipeline is in the RNA-Seq_pipeline.2112.Rmd file. You need to modify some lines to adapt this script to your data and your project.

The header and first paragraph summarize the project name, the genome reference used, and the objectives of this project. Please be the more complete and comprehensive possible, as the report generated can be share with other teams.

#### biomaRt

On the plateform we are mainly using the Ensembl database, so the script uses the `biomaRt` package. As the command line for mouse and human are indicated, you may want to change the reference organism.

Once your new reference configurated, please be sure you ask for the exact same amount of information during the `getBM` command: 
```
getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol','description','go_id','name_1006'), filters = 'ensembl_gene_id', values = ensemblIDs, mart = ensembl,useCache = FALSE)
```
You can determine you new attributes by using the ListAttribute() on your ensembl object.

#### figures

All the commands and functions to generate the figures are listed in the GenomIC_RNASeq-pipeline_functions.r file. You can modify the setup of each figure here.

#### saving the data

At the end of this script, results are save in a results.rds file, with a particular configuration. Do not modify this file if you need to use our shiny application (currently in development) to generate new personnalized figures.

## Getting some results

Once the script and the configuration files are ready, you can generate the matrices, figures and report by simply knitting the script. The filename are obviously very long, but it summarizes the project, the database, the contrast and the nature of the file. 

# Sources

We based our pipeline on the recommended DESeq2 script, with the default figures for RNA-Seq data ggplot2 can do.

“Create Elegant Data Visualisations Using the Grammar of Graphics • Ggplot2.” n.d. Accessed April 16, 2021. https://ggplot2.tidyverse.org/.

Dobin, Alexander, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, and Thomas R. Gingeras. 2013. “STAR: Ultrafast Universal RNA-Seq Aligner.” Bioinformatics 29 (1): 15–21. https://doi.org/10.1093/bioinformatics/bts635.

Li, Bo, and Colin N. Dewey. 2011. “RSEM: Accurate Transcript Quantification from RNA-Seq Data with or Without a Reference Genome.” BMC Bioinformatics 12 (1): 323. https://doi.org/10.1186/1471-2105-12-323.

Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. “Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.” Genome Biology 15 (12): 550. https://doi.org/10.1186/s13059-014-0550-8.
