################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

# LOAD PACKAGES
library(dlnm); library(mixmeta)
library(tsModel) ; library(splines) ; library(gnm) ; library(mgcv)
library(data.table)
library(dplyr) ; library(tidyr); library(lubridate)
library(foreach) ; library(doParallel)
library(MASS)
library(abind)
library(scales) ; library(corrplot)
library(tmap) ; library(ggplot2) ; library(sf)
library(grid) ; library(patchwork)
library(openxlsx) ; library(readxl) 