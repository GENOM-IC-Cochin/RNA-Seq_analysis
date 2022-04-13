# RNA-Seq_analysis


# Table of Contents

-   [Introduction](#orge3da5dc)
-   [Usage](#org7ffbcd4)
    -   [Prerequisite](#org520eba4)
    -   [Configuration](#orga7464b0)
-   [Contact](#org22656d3)
-   [Acknowledgement](#orgd45dd1a)



<a id="orge3da5dc"></a>

# Introduction

This is a R pipeline for the analyse of RNA-Seq data, from a raw count matrix to a list of differentially expressed genes, with figures and reports. The results.rds output is to be used with our SHARE web application.

It follows a STAR -> RSEM pipeline, available on the plateform.

<a id="org7ffbcd4"></a>

# Usage

<a id="org520eba4"></a>

## Prerequisite

To run the pipeline, a Renv restoration is needed

-   R version >=4: <https://www.r-project.org>
-   The corresponding version of the pipeline, depending on the date of your project
-   The .genes.results files from RSEM

<a id="orga7464b0"></a>

## Configuration

The .conf and .contrastes files must be updated with the right metadata and count files. These configuration files must be tabulated, the column names have to be conserved and one must verify the good execution of the read.table lines at the beginning of the RNA-Seq_pipeline.Rmd script.

Edition of lines in the RNA-Seq_pipeline.Rmd to correspond to the project:

-  the ensembl.database 
-  the configuration columns
-  the dds design


<a id="org22656d3"></a>

# Contact

If you have any comment, suggestion or question, feel free to post an issue.

