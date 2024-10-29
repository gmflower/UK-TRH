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

# Keep 2002 onwards:
files <- files[yearseq>=2008]
yearseq <- yearseq[yearseq>=2008]

# Primary diagnosis only:
diag <- sprintf("DIAG_%02d",1:1)

# AGE GROUPS AND CUT-OFF POINTS
agelab <- c("0-64", "65-74", "75-84", "85 and above")
agevarlab <- c("age064", "age6574", "age7584", "age85plus")
agecut <- c(0, as.numeric(substr(agelab[-1],1,2))-1, 150)

# AGE AND SEX
#agegr <- ...
#sex <- NULL

# CAUSE(S) TO SELECT: COMBINATION OF LETTERS AND LENGTH OF CHARACTERS
causelist <- list(
  #list(name="All", shortname="all", icd=LETTERS, nchar=1),
  list(name="Infectious and parasitic", shortname="infect", icd=c("A","B"), nchar=1),
  list(name="Lung cancer", shortname="lungcancer", icd=paste0("C",33:34), nchar=3),
  list(name="Endocrine, nutritional, metabolic", shortname="endo", icd=c("E"), nchar=1),
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
  list(name="Hypotension", shortname="hypo", icd=paste0("I",95), nchar=3),
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
  list(name="Genitourinary", shortname="genito", icd="N", nchar=1),
  list(name="Renal disease", shortname="renal", icd=paste0("N",0:3), nchar=2), 
  list(name="Acute renal failure", shortname="arf", icd=paste0("N",17), nchar=3),
  list(name="Chronic kidney disease", shortname="ckd", icd=paste0("N",18), nchar=3),
  #list(name="Pregnancy and perinatal", shortname="pregnancy", icd=c("O","P"), nchar=1),
  #list(name="Congenital", shortname="congenital", icd="Q", nchar=1),
  list(name="Accidents and injuries", shortname="accidents", icd=c(paste0("V", 0:9), paste0("W",0:9), paste0("X",0:5)), nchar=2),
  list(name="Intentional self harm", shortname="selfharm", icd=paste0("X",60:84), nchar=3),
  list(name="Heat related illness", shortname="heat", icd=c(paste0("E", 86), paste0("T",67), paste0("X",30)), nchar=3)
)


# CREATE EMPTY OBJECT
outdata <- NULL

# RUN THE LOOP (BACKWARD IN TIME)
for(i in seq(files)) {
  
  # PRINT
  cat("\n", files[i], "")
  
  # READ AND TRANSFORM IN DATA.TABLE
  data <- readRDS(paste0(dir,"/emergency/",files[i])) |> as.data.table()
  
  for (j in seq(causelist)) {
    
    # SELECT CAUSE
    out <- causelist[[j]][[2]]
    icd <- causelist[[j]][[3]]
    icddig <- causelist[[j]][[4]]
    
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
    
    # MERGE
    if(j==1) yeardata <- subset else
      #yeardata <- merge(yeardata, subset, all=T, by=c("LSOA11CD", "date", "agegr"))
      yeardata <- rbind(yeardata, subset)
  }
    
  # APPEND
  outdata <- rbind(outdata, yeardata)
}
  
# FILL NA
#outdata <- setnafill(outdata, type=c("const"), fill=0, nan=NA,
#  cols=4:ncol(outdata))

# ORDER AND TRANSFORM IN A SIMPLE DATAFRAME
setkey(outdata, cause, LSOA11CD, date, agegr)
outdata <- as.data.frame(outdata)
  
# Save merged file
saveRDS(outdata, file=paste0("~/RProjects/UK-TRH/temp/simplified/emrgcountHES.RDS"))


