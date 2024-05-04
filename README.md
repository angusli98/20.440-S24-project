# Profiling Transcriptional Dysregulation of Tumor Stroma and Metastatic Lymph Node
## Overview
This repository contains files required for performing analysis and data visualization, as well as outputs thereof, for a term project for the Spring 2024 semester of 20.440: Analysis of Biological Networks at the Massachusetts Institute of Technology. A lung adenocarcinoma single-cell RNA sequencing dataset is downloaded using an automated tool, and gene expression analyses are performed on the stroma of primary tumor and metastatic lymph nodes.

This project uses the R Project (https://www.r-project.org/), Bioconductor (https://www.bioconductor.org/), TMExplorer (https://github.com/shooshtarilab/TMExplorer), Seurat 5 (https://satijalab.org/seurat/), rbioapi (https://cran.r-project.org/web/packages/rbioapi/index.html), and the PANTHER Classification System (https://www.pantherdb.org/).
### Full package literature citations
Gentleman, R. C., Carey, V. J., Bates, D. M., Bolstad, B., Dettling, M., Dudoit, S., Ellis, B., Gautier, L., Ge, Y., Gentry, J., Hornik, K., Hothorn, T., Huber, W., Iacus, S., Irizarry, R., Leisch, F., Li, C., Maechler, M., Rossini, A. J., Sawitzki, G., … Zhang, J. (2004). Bioconductor: open software development for computational biology and bioinformatics. Genome biology, 5(10), R80. https://doi.org/10.1186/gb-2004-5-10-r80

Christensen, E., Naidas, A., Chen, D., Husic, M., & Shooshtari, P. (2022). TMExplorer: A tumour microenvironment single-cell RNAseq database and search tool. PloS one, 17(9), e0272302. https://doi.org/10.1371/journal.pone.0272302

Hao, Y., Stuart, T., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Molla, G., Madad, S., Fernandez-Granda, C., & Satija, R. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature biotechnology, 42(2), 293–304. https://doi.org/10.1038/s41587-023-01767-y

Rezwani, M., Pourfathollah, A. A., & Noorbakhsh, F. (2022). rbioapi: user-friendly R interface to biologic web services' API. Bioinformatics (Oxford, England), 38(10), 2952–2953. https://doi.org/10.1093/bioinformatics/btac172

Thomas, P. D., Ebert, D., Muruganujan, A., Mushayahama, T., Albou, L. P., & Mi, H. (2022). PANTHER: Making genome-scale phylogenetics accessible to all. Protein science : a publication of the Protein Society, 31(1), 8–22. https://doi.org/10.1002/pro.4218

Mi, H., Muruganujan, A., Huang, X., Ebert, D., Mills, C., Guo, X., & Thomas, P. D. (2019). Protocol Update for large-scale genome and gene function analysis with the PANTHER classification system (v.14.0). Nature protocols, 14(3), 703–721. https://doi.org/10.1038/s41596-019-0128-8
## Data
Source data was generated via single-cell RNA sequencing of lung adenocarcinoma patient tissue samples using a 10x Genomics platform. Data files are too large to be uploaded directly to this repository but are downloaded directly during execution of the analysis code using TMExplorer (https://github.com/shooshtarilab/TMExplorer). Broadly, data are formatted as transcripts per million reads, per gene, per cell.
### Public availability of original dataset
Kim, N., Kim, H. K., Lee, K., Hong, Y., Cho, J. H., Choi, J. W., Lee, J. I., Suh, Y. L., Ku, B. M., Eum, H. H., Choi, S., Choi, Y. L., Joung, J. G., Park, W. Y., Jung, H. A., Sun, J. M., Lee, S. H., Ahn, J. S., Park, K., Ahn, M. J., … Lee, H. O. (2020). Single-cell RNA sequencing demonstrates the molecular and cellular reprogramming of metastatic lung adenocarcinoma. Nature communications, 11(1), 2285. https://doi.org/10.1038/s41467-020-16164-1

Accessed in repackaged format using TMExplorer: Christensen, E., Naidas, A., Chen, D., Husic, M., & Shooshtari, P. (2022). TMExplorer: A tumour microenvironment single-cell RNAseq database and search tool. PloS one, 17(9), e0272302. https://doi.org/10.1371/journal.pone.0272302
## Folder structure
The repository is divided into three subdirectories: code, data, and figures.

The code subdirectory contains an RStudio project file to initialize an appropriate environment for running analyses in RStudio, as well as R scripts for running analyses and generating figures.

The data subdirectory is empty because TMExplorer is used to download data files directly during analysis. As such, no local data files are required nor accessed.

The figures subdirectory contains the output figures and tables generated during data analysis performed by the scripts contained within the code subdirectory. Execution of scripts within the code subdirectory will recreate the files within the figures subdirectory.
## Installation
This project uses R, version 4.3.3 or later, which can be downloaded from https://cran.r-project.org/. The R packages Bioconductor, TMExplorer, Seurat 5, rbioapi, ggplot2, and here are also required.

Scripts for this project can optionally be run in the RStudio integrated development environment, available at https://posit.co/download/rstudio-desktop/.

An active Internet connection is required. 16 GB of RAM or greater are recommended. Due to default memory allocation behavior on Mac, the virtual memory may have to be increased for the dataset to be loaded. To do this in the terminal, execute
```
cd ~
touch .Renviron
open .Renviron
```
Then, add `R_MAX_VSIZE=100Gb` into .Renviron.

Packages required for this project are installed through Bioconductor, version 3.18 or later. To install Bioconductor, the following commands can be executed in the R console:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")
```
TMExplorer can be installed by executing the following commands in the R console:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TMExplorer")
```
Seurat 5 can be installed by executing the following command in the R console:
```
install.packages('Seurat')
```
rbioapi can be installed by executing the following command in the R console:
```
install.packages("rbioapi")
```
This project also uses the ggplot2 and here packages, which are included in standard installations of R.