################################################################################
#
# HES DATA
# EXTRACT SPECIFIC OUTCOMES FROM EMERGENCY ADMISSIONS
# AGGREGATE THE OUTPUT BY LSOA/DATE
#
################################################################################

# LOAD PACKAGES
library(data.table)

# SET THE FOLDERS
dir <- "V:/VolumeQ/AGteam/HES"

# SELECT THE FILES (IN REVERSE ORDER)
files <- list.files(paste0(dir, "/emergency"))
files <- files[substr(files, 1, 7)=="emrgHES"]
yearseq <- strsplit(files, "emrgHES", fixed=T) |> sapply("[[", 2) |> 
  substr(1, 4) |> as.numeric()

# REMOVE PRE-2008
files <- files[yearseq>=2008]
yearseq <- yearseq[yearseq>=2008]

# NUMBER OF DIAGNOSIS TO KEEP
diag <- sprintf("DIAG_%02d",1:1)

# AGE GROUPS AND CUT-OFF POINTS
agelab <- c("0-64", "65-74", "75-84", "85 and above")
agevarlab <- c("age064", "age6574", "age7584", "age85plus")
agecut <- c(0, as.numeric(substr(agelab[-1],1,2))-1, 150)

# AGE AND SEX
#agegr <- ...
#sex <- NULL

# CREATE EMPTY OBJECT
outdata <- NULL

# RUN THE LOOP (BACKWARD IN TIME)
for(i in seq(files)) {
  
  # PRINT
  cat("\n", files[i], "")
  
  # READ AND TRANSFORM IN DATA.TABLE
  data <- readRDS(paste0(dir,"/emergency/",files[i])) |> as.data.table()
  
    # SELECT COPD
    out <- "copd"
    icd <- paste0("J",40:44)
    icddig <- 3
    
    # PRINT
    cat(out,"")

    # FUNCTION TO SELECT ICD
    ficd <- function(x, icd, icddig) substr(x, 1, icddig) %in% icd
  
    # SELECT THE ICD CODES
    data[, ind1:=apply(sapply(.SD, ficd, icd, icddig), 1, any), .SDcols=diag]
    subset <- data[(ind1),]
  
    # CREATE AGE GROUP
    subset[, agegr:=cut(ADMIAGE, agecut, labels=agevarlab, include.lowest=T)]
    
    # AGGREGATE BY LSOA AND ADMISSION DATE, TRANSFORM IN DATAFRAME, RENAME
    subset <- subset[, list(sum(ind1)), by=c("LSOA11", "ADMIDATE","agegr")]
    #names(subset) <- c("LSOA11CD","date","agegr",out)
    names(subset) <- c("LSOA11CD","date","agegr","count")
    subset$cause <- out
    setkey(subset, LSOA11CD, date, cause, agegr)
    data$ind1 <- NULL
    
  # APPEND
  outdata <- rbind(outdata, subset)
}
  
# FILL NA
#outdata <- setnafill(outdata, type=c("const"), fill=0, nan=NA,
#  cols=4:ncol(outdata))

# ORDER AND TRANSFORM IN A SIMPLE DATAFRAME
setkey(outdata, cause, LSOA11CD, date, agegr)
outdata <- as.data.frame(outdata)
  
# Save merged file
saveRDS(outdata, file=paste0("~/RProjects/UK-TRH/temp/COPD/emrgcountHES.RDS"))


