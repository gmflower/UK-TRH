################################################################################
# UK-TRH: SMALL-AREA ANALYSIS OF TEMPERATURE RELATED HOSPITALISATIONS IN ENGLAND
################################################################################

################################################################################
# LOAD PACKAGES AND SET PARAMETERS (REPRODUCABLE)
################################################################################

################################################################################
# LOAD PACKAGES

library(dlnm); library(mixmeta)
library(tsModel) ; library(splines) ; library(pbs)
library(gnm) ; library(mgcv)
library(data.table)
library(dplyr) ; library(tidyr); library(tibble) ; library(lubridate)
library(foreach) ; library(doParallel)
library(MASS)
library(abind)
library(scales) ; library(corrplot)
library(ggplot2) ; 
library(grid) ; library(patchwork)
library(openxlsx) ; library(readxl) 
library(Epi)
library(paletteer)
library(ggpubr)
library(cowplot)

################################################################################
# PARAMETERS

# PATH TO DIRECTORIES
lookuppath <- "V:/VolumeQ/AGteam/ONS/geography/lookup/"
hosppath <- "V:/VolumeQ/AGteam/HES"
envpath <- "V:/VolumeQ/AGteam/SmallArea/processed/"
holypath <- "V:/VolumeQ/AGteam/Holidays"

# CREATE FOLDERS (IF NEEDED)
if(!"temp" %in% list.files()) dir.create("temp")
if(!"figures" %in% list.files()) dir.create("figures")
if (!"figures/supplementary" %in% list.files("figures")) {
  dir.create("figures/supplementary", recursive = TRUE)}
if(!"tables" %in% list.files()) dir.create("tables")

# AGE GROUPS AND CUT-OFF POINTS
agelab <- c("18-64", "65-74", "75-84", "85 and above")
agevarlab <- c("age1864", "age6574", "age7584", "age85plus")
agecut <- c(18, as.numeric(substr(agelab[-1],1,2))-1, 150)

# SELECTION OF YEARS AND MONTHS (FOR SEASONAL ANALYSIS)
seqyear <- 2008:2019
seqmonth <- 6:8
  
# SEQUENCE OF PERCENTILES
predper <- c(seq(0,1,0.1), 2:98, seq(99,100,0.1))

# EXPOSURE-RESPONSE PARAMETERIZATION
varfun <- "ns"
varper <- c(50,90)

# LAG-REPONSE PARAMETERIZATION
maxlag <- 3
lagfun <- "ns"
lagknots <- logknots(maxlag, 1)

# DF/YEAR FOR TIME AND DF FOR DOY TO CONTROL OF TREND & SEASONALITY
dftrend <- 1
dfseas <- 4

# Reference list of hospital admission causes:
causelist <- list(
  list(name="Infectious and parasitic", shortname="infect", icd=c("A","B"), nchar=1),
  list(name="Bacterial diseases", shortname="bacterial", icd=c(paste0("A", 20:28), paste0("A",30:49)), nchar=3),
  list(name="Endocrine, nutritional, metabolic", shortname="endo", icd=c("E"), nchar=1),
  list(name="Diabetes mellitus", shortname="diabetes", icd=paste0("E",10:14), nchar=3),
  list(name="Metabolic disorders", shortname="metabolic", icd=c(paste0("E", 70:90)), nchar=3),
  list(name="Cardiovascular", shortname="cvd", icd="I", nchar=1), 
  list(name="Myocardial infarction", shortname="mi", icd=paste0("I",21:23), nchar=3),
  list(name="Heart failure", shortname="hf", icd=paste0("I",50), nchar=3),
  list(name="Hypotension", shortname="hypo", icd=paste0("I",95), nchar=3),
  list(name="Stroke", shortname="stroke", icd=paste0("I",60:69), nchar=3),
  list(name="Respiratory", shortname="resp", icd="J", nchar=1),
  list(name="Acute respiratory infection", shortname="ari", icd=c(paste0("J0", 0:6), paste0("J",20:22)), nchar=3),
  list(name="Pneumonia", shortname="pneumonia", icd=paste0("J",12:18), nchar=3),
  list(name="COPD", shortname="copd", icd=paste0("J",40:44), nchar=3),
  list(name="Asthma", shortname="asthma", icd=paste0("J",45:46), nchar=3),
  list(name="Genitourinary", shortname="genito", icd="N", nchar=1),
  list(name="Renal disease", shortname="renal", icd=paste0("N",0:3), nchar=2),  
  list(name="Acute renal failure", shortname="arf", icd=paste0("N",17), nchar=3)
)
