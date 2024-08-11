################################################################################
# UK-TRH: SMALL-AREA ANALYSIS OF TEMPERATURE RELATED HOSPITALISATIONS IN ENGLAND
################################################################################

################################################################################
# PREPARE THE MAIN AND OTHER DATASETS
################################################################################

################################################################################
# LOOKUP

# LOAD LOOKUP (NOTE THE FILE ENCODING TO PREVENT WEIRD NAMES)
lookup <- read.csv(paste0(lookuppath, "Lower_Layer_Super_Output_Area__2011__to_",
  "Built-up_Area_Sub-division_to_Built-up_Area_to_Local_Authority_District_to_",
  "Region__December_2011__Lookup_in_England_and_Wales.csv"),
  fileEncoding="UTF-8-BOM")

# EXCLUDE WALES, ISLES OF SCILLY, AND CITY OF LONDON
lookup <- subset(lookup, 
  substr(lookup$LAD11CD,1,1)!="W" & !LAD11CD %in% c("E06000053","E09000001"))

# ORDER AND TRANSFORM
lookup <- lookup[with(lookup, order(LSOA11CD, LAD11CD)),]

# LISTS OF LSOA AND LAD (SORTED)
listlsoa <- unique(lookup$LSOA11CD)
listlad <- sort(unique(lookup$LAD11CD))

################################################################################
# MAIN DATASETS

# LOAD HOSPITALISATIONS DATA, TRANFORM DATE 
hesdata <- as.data.table(readRDS(paste0(hosppath, "emrgcountHES.RDS")))
hesdata$date <- as.Date(hesdata$date)

# SELECT YEARS AND EXCLUDE NON-MATCHING LSOA
#hesdata <- hesdata[!is.na(agegr),]
hesdata <- hesdata[year(date) %in% seqyear & LSOA11CD %in% listlsoa,]

# Combine all age group counts (including NA):
hes_allages <- hesdata %>% 
  group_by(LSOA11CD, date, cause) %>% 
  summarise(count = sum(count)) %>% 
  as.data.table()
hes_allages$agegr <- "total"

# And append to age-specific data:
hesdata <- rbind(hesdata, hes_allages)

# Now remove NAs:
hesdata <- hesdata[!is.na(agegr),]

# SELECT ONLY SUMMER MONTHS
hesdata <- hesdata[month(date) %in% seqmonth,]

# CREATE A LIST OF CAUSES
setcause <- sort(unique(hesdata$cause))

# SET KEYS
setkey(hesdata, LSOA11CD, date)

# LOAD THE TEMPERATURE DATA
listtmean <- lapply(seqyear, function(y) {
  cat(y, "")
  file <- paste0("lsoa_daily_temp_", y, ".RDS")
  out <- data.table(readRDS(paste(tmeanpath, file, sep="/")))
  setkey(out, LSOA11CD, date)
  out
})
datatmean <- do.call(rbind, listtmean)
rm(listtmean)

# SELECT SUMMER MONTHS, EXCLUDE NON-MATCHING LSOA
datatmean <- datatmean[LSOA11CD %in% listlsoa,]

# SELECT SUMMER MONTHS
datatmean <- datatmean[month(date) %in% seqmonth,]

# SET KEYS
setkey(datatmean, LSOA11CD, date)

# MERGE LSOA AND LAD ONTO HES AND TMEAN DATA 
hesdata <- as.data.table(lookup[,c("LSOA11CD", "LAD11CD")]) |>
  merge(hesdata, by="LSOA11CD")
datatmean <- as.data.table(lookup[,c("LSOA11CD", "LAD11CD")]) |>
  merge(datatmean, by="LSOA11CD")

# HOLIDAYS FOR ENGLAND
holy <- readRDS(paste0(holypath, "/holidays_GB_1980_2025.RDS")) |>
  subset(select=c(date, GB_ENG)) |> rename(holy=GB_ENG) |> as.data.table()

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
# DERIVE STATS FOR DEFINING BASELINE EVENTS AT LSOA LEVEL

# MULTIPLICATIVE FACTOR FOR SEASONALITY OF EVENTS BY CAUSE AND AGE GROUP
seasevent <- lapply(setcause, function(x) {
  
  # RETRIEVE THE TOTAL COUNTS BY AGE GROUP AND DAY OF THE YEAR
  data <- hesdata[cause==x]
  data[, doy:=yday(date)]
  data <- data[, list(count=sum(count)), by=c("agegr","doy")]
  
  # FUNCTION TO COMPUTE THE SMOOTHED PROPORTION
  funbase <- function(count, doy) {
    mod <- glm(count ~ pbs(doy,4,Bound=c(1,366)), family=poisson())
    pred <- predict(mod, newdata=data.frame(doy=1:366), type="response")
    data.table(doy=1:366, mult=pred/mean(pred))
  }
  
  # APPLY THE FUNCTION, ADD THE CAUSE, REORDER
  data <- data[, funbase(count, doy), by="agegr"]
  data$cause <- x
  setcolorder(data, "cause")
  setkey(data, cause, agegr, doy)

  # RETURN
  data
}) |> Reduce(rbind, x=_)

# AVERAGE DAILY EVENTS BY LAD AND LSOA/LAD POP PROPORTION BY CAUSE AND AGE GROUP
ladevent <- hesdata[, list(count=sum(count)/length(seqyear)/365.25),
  by=c("LAD11CD","agegr","cause")]
lsoaprop <- as.data.table(lookup[,c("LSOA11CD", "LAD11CD")]) |>
  merge(pop, by="LSOA11CD")
lsoaprop[, prop:=pop/sum(pop), by=c("LAD11CD","agegr")]

################################################################################
# CENSUS VARIABLES

# % OF POP 65 AND OVER
dir <- "V:/VolumeQ/AGteam/NOMIS/census2011/KeyStatistics/"
KS102EW <- read.csv(paste0(dir, "KS102EW.csv"))
poptot65 <- data.frame(LSOA11CD=KS102EW$geography.code,
  poptot=KS102EW[[5]],
  pop65=rowSums(KS102EW[18:21])/KS102EW[[5]])
rm(KS102EW)

# POPULATION DENSITY
dir <- "V:/VolumeQ/AGteam/NOMIS/census2011/KeyStatistics/"
KS101EW <- read.csv(paste0(dir, "KS101EW.csv"))
popdens <- data.frame(LSOA11CD=KS101EW$geography.code, popdens=KS101EW[[12]])
rm(KS101EW)

# IMD ENGLAND (2015)
dir <- "V:/VolumeQ/AGteam/GOV/England/IMD/2015/"
engimd <- read_excel(paste0(dir, "File_2_ID_2015_Domains_of_deprivation.xlsx"),
  sheet=2)
engimd <- engimd[c(1,5,7,9,11,13,17,19)]
names(engimd) <- c("LSOA11CD", paste0("IMD", c("","income","employ","educ","health",
  "housing","environ")))

# IMD WALES (2014)
dir <- "V:/VolumeQ/AGteam/GOV/Wales/IMD/2014/"
walimd <- read_excel(paste0(dir, "170213 WIMD Ranks and Deciles 2014.xlsx"),
  sheet=2, skip=3, n_max=1909)
walimd <- walimd[c(1,6,8,10,14,12,22,20)]
names(walimd) <- names(engimd)

# IMD (STANDARDIZED AND APPENDED)
engimd[-1] <- lapply(engimd[-1], function(x) x/max(x))
walimd[-1] <- lapply(walimd[-1], function(x) x/max(x))
imd <- rbind(engimd, walimd)
rm(engimd,walimd)

# % HOUSE WITH CENTRAL HEATING
dir <- "V:/VolumeQ/AGteam/NOMIS/census2011/KeyStatistics/"
KS403EW <- read.csv(paste0(dir, "KS403EW.csv"))
heating <- data.frame(LSOA11CD=KS403EW$geography.code,
  heating=KS403EW[[7]]/KS403EW[[5]])
rm(KS403EW)

################################################################################
# SATELLITE DATA

# ALBEDO
albedo <- readRDS("data/ALBEDO_500m/BSA2011_lsoa.RDS")[,c(1,6)]
names(albedo) <- c("LSOA11CD", "albedo")

# EVI
evisum <- readRDS("data/EVI_250m/EVI_Jul2011_lsoa.RDS")[,c(1,6)]
names(evisum) <- c("LSOA11CD", "evisum")
eviwin <- readRDS("data/EVI_250m/EVI_Jan2011_lsoa.RDS")[,c(1,6)]
names(eviwin) <- c("LSOA11CD", "eviwin")
evi <- as.data.frame(merge(evisum, eviwin))
evi$evi <- rowMeans(evi[,-1])
rm(evisum, eviwin)

# IMPERVIOUS SURFACES
imperv <- readRDS("data/Imp_surface100m/IMP_SurF1002012_lsoa.RDS")[,c(1,6)]
names(imperv) <- c("LSOA11CD", "imperv")

# AGE BUILDINGS
building <- read.csv("data/agebuilding/voapropertyage.csv")
cutoffyear1 <- as.numeric(substr(colnames(building)[4:14],4,7))
cutoffyear2 <- as.numeric(substr(colnames(building)[4:14],9,13))
midyear <- c(1900, round((cutoffyear2 - cutoffyear1)/2 + cutoffyear1))
age2015 <- 2015-midyear
agebuild <- data.frame(LSOA11CD=building[[1]],
  agebuild=apply(building[3:14], 1, function(x) sum(x*age2015)/sum(x)))
rm(building)

################################################################################
# TEMPERATURE SUMMARIES

#tmeansumm <-  datafull[, list(meantmean=mean(tmean), 
#  rangetmean=diff(range(tmean))), by=LSOA11CD]
tmeansumm <-  datatmean[, list(meantmean=mean(tmean),
  rangetmean=diff(range(tmean))), by=LSOA11CD]

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
