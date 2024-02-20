################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# PARAMETERS
################################################################################

# PATH TO TEMPERATURE DIRECTORY
tmeanpath <- "V:/VolumeQ/AGteam/ONSmortality/processed/LSOA19812018/"

# SEQUENCE OF DATES OF THE TMEAN PERIOD
date <- seq(as.Date("1981/1/1"), as.Date("2018/12/31"), "days")

# STARTING DATE
datestart <- as.Date("2000/01/01")

# SEQUENCE OF PERCENTILES
predper <- c(seq(0,1,0.1), 2:98, seq(99,100,0.1))

# EXPOSURE-RESPONSE PARAMETERIZATION
varfun <- "bs"
varper <- c(10,75,90)
vardegree <- 2

# LAG-REPONSE PARAMETERIZATION
maxlag <- 21
lagfun <- "ns"
lagknots <- logknots(maxlag, 3)

# DF/YEAR FOR CONTROL OF SEASONALITY
nkseas <- 1

# LISTS OF LSOA AND LAD
listlsoa <- unique(lookup$LSOA11CD)
listlad <- unique(lookup$LAD11CD)
