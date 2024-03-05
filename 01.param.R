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
tmeanpath <- "V:/VolumeQ/AGteam/MetData/Processed/LSOA_ukcp18_1kmgrid_v1030"

# AGE GROUPS AND CUT-OFF POINTS
agelab <- c("0-64", "65-74", "75-84", "85 and above")
agevarlab <- c("age064", "age6574", "age7584", "age85plus")
agecut <- c(0, as.numeric(substr(agelab[-1],1,2))-1, 150)

# SEQUENCE OF YEARS
seqyear <- 2009:2021
  
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


# CAUSE(S) TO SELECT: COMBINATION OF LETTERS AND LENGTH OF CHARACTERS
causelist <- list(
  #list(name="All", shortname="all", icd=LETTERS, nchar=1),
  list(name="Infectious and parasitic", shortname="infect", icd=c("A","B"), nchar=1),
  list(name="Lung cancer", shortname="lungcancer", icd=paste0("C",33:34), nchar=3),
  list(name="Diabetes mellitus", shortname="diabetes", icd=paste0("E",10:14), nchar=3),
  list(name="Mental and behavioural", shortname="mental", icd="F", nchar=1),
  list(name="Organic mental disorders", shortname="dementia", icd=paste0("F0",0:9), nchar=3),
  #list(name="Nervous system", shortname="nervous", icd="G", nchar=1),
  #list(name="Eye diseases", shortname="eye", icd=paste0("H",0:5), nchar=2),
  #list(name="Ear diseases", shortname="ear", icd=paste0("H",6:9), nchar=2),
  list(name="Cardiovascular", shortname="cvd", icd="I", nchar=1), 
  list(name="Ischaemic heart disease", shortname="ihd", icd=paste0("I",20:25), nchar=3),
  list(name="Myocardial infarction", shortname="mi", icd=paste0("I",21:23), nchar=3),
  list(name="Pulmonary heart disease", shortname="pulmheart", icd=paste0("I",26:28), nchar=3),
  list(name="Heart failure", shortname="hf", icd=paste0("I",50), nchar=3),
  list(name="Stroke", shortname="stroke", icd=paste0("I",60:69), nchar=3),
  list(name="Respiratory", shortname="resp", icd="J", nchar=1),
  list(name="Acute respiratory infection", shortname="ari", icd=c(paste0("J0", 0:6), paste0("J",20:22)), nchar=3),
  list(name="Influenza", shortname="flu", icd=c("J09","J10","J11"), nchar=3),
  list(name="Pneumonia", shortname="pneumonia", icd=paste0("J",12:18), nchar=3),
  list(name="COPD", shortname="copd", icd=paste0("J",40:44), nchar=3),
  list(name="Asthma", shortname="asthma", icd=paste0("J",45:46), nchar=3),
  list(name="Bronchiectasis", shortname="bronch", icd=c("J47","E84"), nchar=3),
  #list(name="Interstitial lung disease", shortname="ild", icd=paste0("J", c(60:84, 96:99)), nchar=3),
  #list(name="Digestive diseases", shortname="digestive", icd=paste0("K",0:9), nchar=2),
  #list(name="Skin diseases",  shortname="skin", icd="L", nchar=1),
  #list(name="Musculoskeletal", shortname="musculoskeletal", icd="M", nchar=1),
  #list(name="Genitourinary", shortname="genito", icd="N", nchar=1),
  #list(name="Pregnancy and perinatal", shortname="pregnancy", icd=c("O","P"), nchar=1),
  #list(name="Congenital", shortname="congenital", icd="Q", nchar=1),
  list(name="Accidents and injuries", shortname="accidents", icd=c(paste0("V", 0:9), paste0("W",0:9), paste0("X",0:5)), nchar=2),
  list(name="Intentional self harm", shortname="selfharm", icd=paste0("X",60:84), nchar=3)
)


