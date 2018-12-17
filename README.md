## Purpose
The R ['Script_NN-T61915'](Script_NN-T61915.R) outlines the analysis conducted for the study "Human microglia regional heterogeneity and phenotypes determined by multiplexed single-cell mass cytometry" and reproduces Figures 5_a_b_c_d_f, 6_b_c, 7_b_c_d_e, and S8_h of Böttcher et al. (Nature Neuroscience, 2018, https://doi.org/10.1038/s41593-018-0290-2). The accompanying R package ['NNhuMG'](NNhuMG) contains a number of utility functions required to walk through the analysis.

## Installation  
1) go to https://flowrepository.org/id/FR-FCM-ZYM6 and download zip folder with all fcs files

2) unzip and move folder (should be 'FlowRepository_FR-FCM-ZYM6_files') to your working directory

3) install required packages and dependencies (you need to have R (> 3.4.1), https://www.r-project.org/):  

```r
# dependencies
install.packages(c("ks", "dbscan", "feature", "shiny", "rmarkdown", "rgl", "knitr", "robust", "splancs", "flowCore", "lpSolve", "ggplot2", "gridExtra", "grid", "matrixcalc"))

# required to reproduce figures
install.packages(c("devtools", "R.utils", "lme4", "lmerTest", "flowFP", "vegan"))

# try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite(c("flowWorkspace", "CytoML", "flowType", "RchyOptimyx", "edgeR", "cydar"))

# use install.packages("emdist") if 0.3-2 doesn’t work 
library(devtools)
install_github("s-u/emdist") # version 0.3-2
install_github("steschlick/Boettcher_et_al_NN_2018/NNhuMG")

# set your working directory and move 'FlowRepository_FR-FCM-ZYM6_files’ folder here
# start walking through the analysis 'Script_NN-T61915'

```

## Disclaimer
The software in this package is publication version, in an early development stage and not fully tested nor documented. 
The authors maintain no responsibility for any possible code errors or bugs.
 
