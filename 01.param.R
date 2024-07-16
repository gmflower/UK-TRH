################################################################################
# UK-TRH: SMALL-AREA ANALYSIS OF TEMPERATURE RELATED HOSPITALISATIONS IN ENGLAND
################################################################################

################################################################################
# PARAMETERS
################################################################################

# PATH TO DIRECTORIES
lookuppath <- "V:/VolumeQ/AGteam/ONS/geography/lookup/"
#deathpath <- "V:/VolumeQ/AGteam/ONS/mortality/data"
hosppath <- "V:/VolumeQ/AGteam/HES/extracted"
tmeanpath <- "V:/VolumeQ/AGteam/LSOA_Level/processed/temperature"
holypath <- "V:/VolumeQ/AGteam/Holidays"


# AGE GROUPS AND CUT-OFF POINTS
agelab <- c("0-64", "65-74", "75-84", "85 and above")
agevarlab <- c("age064", "age6574", "age7584", "age85plus")
agecut <- c(0, as.numeric(substr(agelab[-1],1,2))-1, 150)

# Selection of years (and months for a summer-only analysis)
#seqyear <- 2009:2021
seqyear <- 2008:2019
smonth <- 6:9
  
# SEQUENCE OF PERCENTILES
predper <- c(seq(0,1,0.1), 2:98, seq(99,100,0.1))

# EXPOSURE-RESPONSE PARAMETERIZATION
#varfun <- "bs"
#varper <- c(10,75,90)
#vardegree <- 2
#vardf <- 5

# EXPOSURE-RESPONSE PARAMETERIZATION
# Seasonal analysis only: 
varfun <- "ns"
varper <- c(50,90)

# LAG-REPONSE PARAMETERIZATION
#maxlag <- 21
# Seasonal analysis only: define a short lag period
maxlag <- 3
lagfun <- "ns"
lagknots <- logknots(maxlag, 2)

# DF/YEAR FOR CONTROL OF SEASONALITY
nkseas <- 1

