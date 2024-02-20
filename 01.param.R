################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# PARAMETERS
################################################################################

# PATH TO DIRECTORIES
lookuppath <- "V:/VolumeQ/AGteam/ONS/geography/lookup/"
deathpath <- "V:/VolumeQ/AGteam/ONS/mortality/data"
tmeanpath <- "V:/VolumeQ/AGteam/MetData/Processed/LSOA_ukcp18_1kmgrid_v1030"

# AGE GROUPS AND CUT-OFF POINTS
agelab <- c("0-64", "65-74", "75-84", "85 and above")
agevarlab <- c("age064", "age6574", "age7584", "age85plus")
agecut <- c(0, as.numeric(substr(agelab[-1],1,2))-1, 150)

# SEQUENCE OF YEARS
seqyear <- 2000:2019

# SEQUENCE OF PERCENTILES
predper <- c(seq(0,1,0.1), 2:98, seq(99,100,0.1))

# EXPOSURE-RESPONSE PARAMETERIZATION
varfun <- "bs"
varper <- c(10,75,90)
vardegree <- 2
vardf <- 5

# LAG-REPONSE PARAMETERIZATION
maxlag <- 21
lagfun <- "ns"
lagknots <- logknots(maxlag, 3)

# DF/YEAR FOR CONTROL OF SEASONALITY
nkseas <- 1
