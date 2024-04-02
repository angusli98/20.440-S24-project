# Profiling Transcriptional Dysregulation of Tumor Stroma and Metastatic Lymph Node
## Overview
This repository contains files required for performing analysis and data visualization, as well as outputs thereof, for a term project for the Spring 2024 semester of 20.440: Analysis of Biological Networks at the Massachusetts Institute of Technology. Currently, a lung adenocarcinoma single-cell RNA sequencing dataset is downloaded using an automated tool, principal components analysis is performed on the data, and a dimensionality reduction plot is produced based on the resultant principal components.

This project uses the R Project (https://www.r-project.org/), Bioconductor (https://www.bioconductor.org/), TMExplorer (https://github.com/shooshtarilab/TMExplorer), and Seurat 5 (https://satijalab.org/seurat/).
### Full package literature citations
Gentleman, R. C., Carey, V. J., Bates, D. M., Bolstad, B., Dettling, M., Dudoit, S., Ellis, B., Gautier, L., Ge, Y., Gentry, J., Hornik, K., Hothorn, T., Huber, W., Iacus, S., Irizarry, R., Leisch, F., Li, C., Maechler, M., Rossini, A. J., Sawitzki, G., … Zhang, J. (2004). Bioconductor: open software development for computational biology and bioinformatics. Genome biology, 5(10), R80. https://doi.org/10.1186/gb-2004-5-10-r80

Christensen, E., Naidas, A., Chen, D., Husic, M., & Shooshtari, P. (2022). TMExplorer: A tumour microenvironment single-cell RNAseq database and search tool. PloS one, 17(9), e0272302. https://doi.org/10.1371/journal.pone.0272302

Hao, Y., Stuart, T., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Molla, G., Madad, S., Fernandez-Granda, C., & Satija, R. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature biotechnology, 42(2), 293–304. https://doi.org/10.1038/s41587-023-01767-y
## Data
Source data was generated via single-cell RNA sequencing of lung adenocarcinoma patient tissue samples using a 10x Genomics platform. Data files are too large to be uploaded directly to this repository but are downloaded directly during execution of the analysis code using TMExplorer (https://github.com/shooshtarilab/TMExplorer). Broadly, data are formatted as transcripts per million reads, per gene, per cell.
### Public availability of original dataset
Kim, N., Kim, H. K., Lee, K., Hong, Y., Cho, J. H., Choi, J. W., Lee, J. I., Suh, Y. L., Ku, B. M., Eum, H. H., Choi, S., Choi, Y. L., Joung, J. G., Park, W. Y., Jung, H. A., Sun, J. M., Lee, S. H., Ahn, J. S., Park, K., Ahn, M. J., … Lee, H. O. (2020). Single-cell RNA sequencing demonstrates the molecular and cellular reprogramming of metastatic lung adenocarcinoma. Nature communications, 11(1), 2285. https://doi.org/10.1038/s41467-020-16164-1

Accessed in repackaged format using TMExplorer: Christensen, E., Naidas, A., Chen, D., Husic, M., & Shooshtari, P. (2022). TMExplorer: A tumour microenvironment single-cell RNAseq database and search tool. PloS one, 17(9), e0272302. https://doi.org/10.1371/journal.pone.0272302
## Folder structure
The repository is divided into three subdirectories: code, data, and figures.

The code subdirectory contains an RStudio project file to initialize an appropriate environment for running analyses in RStudio, as well as R scripts for running analyses and generating figures. There is presently a single R script, which downloads the dataset used for this project using TMExplorer, uses Seurat to normalize the data and perform principal components analysis, and creates a dimensionality reduction plot using the principal components.

The data subdirectory is empty because TMExplorer is used to download data files directly during analysis. As such, no local data files are required nor accessed.

The figures subdirectory contains the output figures generated during data analysis performed by the scripts contained within the code subdirectory. Currently, there is a single figure, which shows the results of dimensionality reduction via principal components analysis. Execution of scripts within the code subdirectory will recreate the figures within the figures subdirectory.
## Installation
This project uses R, version 4.3.3 or later, which can be downloaded from https://cran.r-project.org/. The R packages Bioconductor, TMExplorer, Seurat 5, ggplot2, and here are also required.

Scripts for this project can optionally be run in the RStudio integrated development environment, available at https://posit.co/download/rstudio-desktop/.

An active Internet connection is required. 16 GB of RAM or greater are recommended.

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
This project also uses the ggplot2 and here packages, which are included in standard installations of R.

To run the script for generating the first figure in the R console, set the working directory to the project code subdirectory using `setwd("PATH_TO_CODE_SUBDIRECTORY")`, then run
```
source("PCA_plot_fig1.R")
```
To run in RStudio, open code.Rproj within the code subdirectory. Then, open PCA_plot_fig1.R and click the Source button. Alternatively, enter `source("PCA_plot_fig1.R")` within the console window.