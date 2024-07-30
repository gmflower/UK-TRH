################################################################################
# UK-TRH: SMALL-AREA ANALYSIS OF TEMPERATURE RELATED HOSPITALISATIONS IN ENGLAND
################################################################################

################################################################################
# PARAMETERS
################################################################################

# PATH TO DIRECTORIES
lookuppath <- "V:/VolumeQ/AGteam/ONS/geography/lookup/"
#deathpath <- "V:/VolumeQ/AGteam/ONS/mortality/data"
hosppath <- "~/RProjects/UK-TRH/temp/COPD/"
tmeanpath <- "V:/VolumeQ/AGteam/LSOA_Level/processed/temperature/met_v12"
holypath <- "V:/VolumeQ/AGteam/Holidays"

# AGE GROUPS AND CUT-OFF POINTS
agelab <- c("0-64", "65-74", "75-84", "85 and above")
agevarlab <- c("age064", "age6574", "age7584", "age85plus")
agecut <- c(0, as.numeric(substr(agelab[-1],1,2))-1, 150)

# SELECTION OF YEARS AND MONTHS (FOR SEASONAL ANALYSIS)
seqyear <- 2008:2018
seqmonth <- 6:8
  
# SEQUENCE OF PERCENTILES
predper <- c(seq(0,1,0.1), 2:98, seq(99,100,0.1))

# EXPOSURE-RESPONSE PARAMETERIZATION
varfun <- "ns"
varper <- c(50,90)

# LAG-REPONSE PARAMETERIZATION
maxlag <- 2
lagfun <- "ns"
lagknots <- logknots(maxlag, 1)

# DF/YEAR FOR TIME AND DF FOR DOY TO CONTROL OF TREND & SEASONALITY
dftrend <- 1
dfseas <- 4
