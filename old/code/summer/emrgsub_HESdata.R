################################################################################
#
# HES DATA
# CREATE CLEAN DATASETS WITH EMERGENCY ADMISSIONS
#
# OUTPUT: CLEAN YEARLY DATASETS WITH ONLY EMERGENCY ADMISSIONS, FIRST EPISODE
#   IN EACH SPELL, REMOVED SEQUENTIAL ADMISSION WITHIN A GAP PERIOD, WITH LSOA
#
################################################################################

# LOAD PACKAGES
library(data.table)

# SET THE FOLDERS
dir <- "V:/VolumeQ/AGteam/HES"

# LOOKUPS FOR LSOA
lookup1 <- read.csv(paste0("V:/VolumeQ/AGteam/ONS/geography/lookup/",
  "Lower_Layer_Super_Output_Area_(2001)_to_Lower_Layer_Super_Output_Area_",
  "(2011)_to_Local_Authority_District_(2011)_Lookup_in_England_and_Wales.csv"),
  fileEncoding="UTF-8-BOM")
lookup2 <- read.csv(paste0("V:/VolumeQ/AGteam/ONS/geography/lookup/",
  "Lower_Layer_Super_Output_Area__2011__to_",
  "Built-up_Area_Sub-division_to_Built-up_Area_to_Local_Authority_District_to_",
  "Region__December_2011__Lookup_in_England_and_Wales.csv"),
  fileEncoding="UTF-8-BOM")

# SELECT THE FILES
# NB: IN REVERSE ORDER SO TO KEEP OLD ADMISSIONS IN EARLY YEARS
files <- list.files(paste0(dir, "/original")) |> rev()
files <- files[substr(files, 1, 9)=="NIC329869"]
yearseq <- strsplit(files, "APC_", fixed=T) |> sapply("[[", 2) |> 
  substr(1, 4) |> as.numeric()

# SELECT THE VARIABLES TO EXTRACT
vars <- c("ADMIAGE","ADMIDATE","ADMIMETH",sprintf("DIAG_%02d",1:20),"DISDATE",
  "EPIKEY","LSOA01","LSOA11","SEX","SITETRET","Token_Person_ID","TOKEN_PERSON_ID")

# DEFINE ADMISSION METHOD (JUST 2* FOR EMERGENCY ADMISSIONS)
meth <- as.character(2)
methdig <- 1

# NUMBER OF DIAGNOSIS TO KEEP
diag <- sprintf("DIAG_%02d",1:10)

# MINIMUM PERIOD BETWEEN HOSPITALISATIONS
mingap <- 30

# OUTPUT FOLDER
outdir <- paste0(dir,"/emergency/")
if(!dir.exists(outdir)) dir.create(outdir)

# RUN THE LOOP (BACKWARD IN TIME)
for(i in seq(files)) {

  # SELECT
  zipfile <- files[i]
  cat(zipfile, "")
  
  # COPY AND UNZIP
  file.copy(from=paste0(dir, "/original/", zipfile), to=getwd(), overwrite=F)
  unzip(zipfile=zipfile, overwrite=F)
  
  # READ AND ERASE THE FILES
  file <- sub(".zip", ".txt", zipfile, fixed=T)
  data <- fread(file, select=vars)
  file.remove(zipfile, file)

  # REMOVE THE OTHER DIAGNOSIS
  data[, setdiff(sprintf("DIAG_%02d",1:20), diag):=NULL]
  
  # SET DISCHARGE DATES TO na WHERE NULL OR INVALID
  is.na(data$DISDATE) <- data$DISDATE %in% c("1801-01-01","1800-01-01")
  
  # REMOVE MISSING AND IMPLAUSIBLE ADMISSION DATES 
  # SEE DATA DICTIONARY FOR MISSING DATE CODES
  data <- data[!(is.na(ADMIDATE) | 
      ADMIDATE %in% c("1800-01-01", "1801-01-01", "1600-01-01", "1582-10-15"))]
  
  # MERGE WITH DATA FROM PREVIOUS ITERATION (IF NOT FIRST), AND REMOVE FILES
  # NB: SPELLS RECOMPUTED, KEEPING THOSE WITH DISCHARGE IN FOLLOWING YEARS
  if(i>1) {
    temp <- fread(paste0("temp", yearseq[i]+1,".txt"))
    temp$spell <- NULL
    data <- rbind(data, temp)
    rm(temp)
    file.remove(paste0("temp", yearseq[i]+1,".txt"))
  }

  # DEFINE SPELLS: ORDER (WITH NA FIRST) AND COUNT NON-MISSING DISCHARGE
  # NB: KIND OF BOTTLENECK IN TERMS OF COMPUTING TIME
  setorder(data, Token_Person_ID, ADMIDATE, DISDATE, na.last=FALSE)
  data[, spell:=cumsum(as.numeric((shift(!is.na(DISDATE),1,fill=TRUE)))),
    by=Token_Person_ID]
  
  # TAKE ONLY THE FIRST EPISODE WITHIN A SPELL
  # NB: KEEPING COPY OF DISCHARGE DATE IN THE LAST EPISODE OF THE SPELL
  last <- data[, last(.SD), by=c("Token_Person_ID","spell")] |>
    subset(select=c(Token_Person_ID, spell, DISDATE))
  data <- data[, first(.SD), by=c("Token_Person_ID","spell")] |>
    subset(select=-DISDATE) |> merge(last)
  rm(last)

  # SELECT EMERGENCY BY METHOD OF ADMISSION
  data <- data[substr(ADMIMETH,1,methdig)==meth]
  
  # REMOVE ADMISSIONS WITHIN MINIMUM-GAP DAYS FROM PREVIOUS ADMISSION
  # NB: CONSIDER THAT NON-EMERGENCY ADMISSIONS ARE ALREADY REMOVED
  # NB: THIS ALSO REMOVES SOME REPEATED RECORDS
  setkey(data, Token_Person_ID, ADMIDATE)
  data[, ind2:= c(mingap+1, diff(ADMIDATE)) >= mingap, by=Token_Person_ID]
  data <- data[(ind2)]
  data[, ind2:=NULL]

  # MANAGE LOCATION DATA (SET TO MISSING IF EMPTY OR NOT IN LOOKUPS)
  data[, `:=`(LSOA01=as.character(LSOA01), LSOA11=as.character(LSOA11))]
  data[, `:=`(LSOA01=fifelse(LSOA01=="", NA_character_, LSOA01),
    LSOA11=fifelse(LSOA11=="", NA_character_, LSOA11))]
  data[!LSOA01 %in% lookup1$LSOA01CD, LSOA01 := NA]
  data[!LSOA11 %in% lookup2$LSOA11CD, LSOA11 := NA]
  
  # REMOVE IF NO LOCATION DATA
  data <- data[!(is.na(LSOA01) & is.na(LSOA11))]
  
  # FILL MISSING LSOA11 WITH INFO ON LSOA01
  match <- lookup1$LSOA11CD[match(data$LSOA01, lookup1$LSOA01CD)]
  naind <- is.na(data$LSOA11)
  data$LSOA11[naind] <- match[naind]

  # DIVIDE BY YEAR (ONLY RECORDS WITH ADMISSIONS IN SELECTED YEAR KEPT)
  admicut <- as.Date(paste0(yearseq[i],"-12-31"))
  temp <- data[ADMIDATE <= admicut]
  data <- data[ADMIDATE > admicut]
  
  # SAVE TEMPORARY OF YEAR LATER IF NOT THE LAST ONE, SKIP THE REST IF FIRST
  if(i!=length(files)) fwrite(temp, file=paste0("temp", yearseq[i],".txt"))
  if(i==1) next
  
  # TRANSFORM IN A SIMPLE DATAFRAME
  data <- as.data.frame(data)
  
  # SAVE AS FILE FOR FOLLOWING YEAR (INCLUDING DATA FROM PREVIOUS ITERATION))
  saveRDS(data, file=paste0(outdir, "emrgHES", yearseq[i]+1, ".RDS"))
}
