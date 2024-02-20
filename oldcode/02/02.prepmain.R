################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# PREPARE THE MAIN DATA
################################################################################

################################################################################
# LOOKUP

# LOAD LOOKUP (NOTE THE FILE ENCODING TO PREVENT WEIRD NAMES)
lookup <- read.csv(paste0(lookuppath, "Lower_Layer_Super_Output_Area__2011__to_",
  "Built-up_Area_Sub-division_to_Built-up_Area_to_Local_Authority_District_to_",
  "Region__December_2011__Lookup_in_England_and_Wales.csv"),
  fileEncoding="UTF-8-BOM")

# ORDER AND TRANSFORM
lookup <- lookup[with(lookup, order(LSOA11CD, LAD11CD)),]

# LISTS OF LSOA AND LAD (SORTED)
listlsoa <- unique(lookup$LSOA11CD)
listlad <- sort(unique(lookup$LAD11CD))

################################################################################
# MAIN DATASET

# LOAD MORTALITY DATA AND SELECT YEARS
onsdeath <- as.data.table(readRDS(paste0(deathpath, "/ONSmortality_20211111.RDS")))
onsdeath <- onsdeath[year(DOD) %in% seqyear,]
onsdeath[, agegr:=cut(ageinyrs, agecut, labels=agevarlab, include.lowest=T)]
setkey(onsdeath, lsoa, DOD)

# COLLAPSE AND RESHAPE BY AGE GROUP
onsdeath <- onsdeath[, list(d=length(DOD)),
  by=list(LSOA11CD=lsoa,date=DOD, age=agegr)]
onsdeath <- dcast(onsdeath, LSOA11CD+date~age, value.var="d", fill=0) 
setkey(onsdeath, LSOA11CD, date)

# LOAD THE TEMPERATURE DATA
listtmean <- lapply(seqyear, function(y) {
  cat(y, "")
  file <- paste0("tasmean_lsoa_", y, ".RDS")
  out <- data.table(readRDS(paste(tmeanpath, file, sep="/")))
  setkey(out, LSOA11CD, date)
  out
})
datatmean <- do.call(rbind, listtmean)
rm(listtmean)

# RENAME AND EXCLUDE NON-MATCHING LSOA
datatmean <- datatmean[LSOA11CD %in% listlsoa,]
setkey(datatmean, LSOA11CD, date)

# MERGE THE TWO, KEEPING ALL THE LATTER TO INCLUDE NO-COUNT DAYS
datafull <- merge(onsdeath, datatmean, all.y=T)

# MERGE LAD 
datafull <- merge(as.data.table(lookup[,c("LSOA11CD", "LAD11CD")]), datafull,
  by="LSOA11CD")

# FILL NO-COUNT
datafull[, (agevarlab):=lapply(.SD, nafill, fill=0), .SDcols=agevarlab]

# CREATE TIME VARS
datafull[, time:=as.numeric(date)]
datafull[, year:=year(date)]
datafull[, month:=month(date)]
datafull[, doy:=yday(date)]
datafull[, dow:=wday(date)]

# CLEAN
rm(onsdeath, datatmean)

################################################################################
# RURAL-URBAN CLASSIFICATION

dir <- "V:/VolumeQ/AGteam/ONS/rural_urban/"
file <- paste0("Rural_Urban_Classification_(2011)_of_Lower_Layer_",
  "Super_Output_Areas_in_England_and_Wales.csv")
ruralurban <- read.csv(paste0(dir, file))[c(1,3,4)]
names(ruralurban)[1] <- "LSOA11CD"
ruralurban <- merge(lookup["LSOA11CD"], ruralurban)

################################################################################
# POPULATION BY LSOA AND AGE GROUP

# LOAD FROM NOMIS
dir <- "V:/VolumeQ/AGteam/NOMIS/census2011/KeyStatistics/"
pop <- read.csv(paste0(dir, "KS102EW.csv"))[c(3,6:21)]

# COLLAPSE AGE GROUPS AND TRANSFORM IN LONG
pop <- cbind(pop[1], rowSums(pop[2:13]), pop[[14]], pop[[15]], 
  rowSums(pop[16:17]), rowSums(pop[2:17]))
names(pop) <- c("LSOA11CD", agevarlab, "all")
pop <- pop %>% 
  pivot_longer( -1) %>%
  as.data.frame()
names(pop) <- c("LSOA11CD", "agegr", "pop")
pop <- pop[order(pop$LSOA11CD, pop$agegr),]

# WEIGHTS EUROPEAN STANDARD POPULATION 2013
# NB: SEE V:\VolumeQ\AGteam\ONS\standardization
stdweight <- c(80500, 10500, 6500, 2500) / 100000
names(stdweight) <- unique(pop$agegr)[seq(length(stdweight))]

################################################################################
# SHAPEFILES

# LOAD LSOA SHAPEFILES (ROUGH SUPER GENERALISED)
source <- "V:/VolumeQ/AGteam/ONS/geography/shapefiles"
file <- "Lower_Layer_Super_Output_Areas_(December_2011)_Boundaries_Super_Generalised_Clipped_(BSC)_EW_V3"
file.copy(paste0(source, "/LSOA/", file, "-shp.zip"), getwd())
unzip(zipfile=paste0(file,"-shp.zip"), exdir=getwd())
lsoashp1 <- st_read(paste0(file, ".shp"))[2]
file.remove(list.files()[grep(file, list.files(), fixed=T)])
lsoashp1 <- lsoashp1[match(lookup$LSOA11CD, lsoashp1$LSOA11CD),]

# LOAD LSOA SHAPEFILES (ROUGH GENERALISED)
source <- "V:/VolumeQ/AGteam/ONS/geography/shapefiles"
file <- "Lower_Layer_Super_Output_Areas_(December_2011)_Boundaries_Generalised_Clipped_(BGC)_EW_V3"
file.copy(paste0(source, "/LSOA/", file, "-shp.zip"), getwd())
unzip(zipfile=paste0(file,"-shp.zip"), exdir=getwd())
lsoashp2 <- st_read(paste0(file, ".shp"))[2]
file.remove(list.files()[grep(file, list.files(), fixed=T)])
lsoashp2 <- lsoashp2[match(lookup$LSOA11CD, lsoashp2$LSOA11CD),]

# LOAD LAD SHAPEFILES
source <- "V:/VolumeQ/AGteam/ONS/geography/shapefiles"
file <- "Local_Authority_Districts_(December_2011)_Boundaries_EW_BGC"
file.copy(paste0(source, "/LAD/", file, ".zip"), getwd())
unzip(zipfile=paste0(file,".zip"), exdir=getwd())
ladshp <- st_read(paste0(file, ".shp"))[1]
file.remove(list.files()[grep(file, list.files(), fixed=T)])
names(ladshp)[1] <- "LAD11CD"
ladshp <- ladshp[match(listlad, ladshp$LAD11CD),]
