################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# PREPARE THE DATA
################################################################################

################################################################################
# LOOKUP

# LOAD LOOKUP (NOTE THE FILE ENCODING TO PREVENT WEIRD NAMES)
dir <- "V:/VolumeQ/AGteam/ONS/geography/lookup/"
lookup <- read.csv(paste0(dir, "Lower_Layer_Super_Output_Area__2011__to_",
  "Built-up_Area_Sub-division_to_Built-up_Area_to_Local_Authority_District_to_",
  "Region__December_2011__Lookup_in_England_and_Wales.csv"),
  fileEncoding="UTF-8-BOM")

# ORDER
lookup <- lookup[with(lookup, order(LAD11CD, LSOA11CD)),]

################################################################################
# CENSUS VARIABLES

# % OF POP 65 AND OVER
dir <- "V:/VolumeQ/AGteam/NOMIS/census2011/KeyStatistics/"
KS102EW <- read.csv(paste0(dir, "KS102EW.csv"))
pop <- data.frame(LSOA11CD=KS102EW$geography.code,
  poptot=KS102EW[[5]],
  pop65=rowSums(KS102EW[18:21])/KS102EW[[5]])

# POPULATION DENSITY
dir <- "V:/VolumeQ/AGteam/NOMIS/census2011/KeyStatistics/"
KS101EW <- read.csv(paste0(dir, "KS101EW.csv"))
popdens <- data.frame(LSOA11CD=KS101EW$geography.code, popdens=KS101EW[[12]])

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

# % HOUSE WITH CENTRAL HEATING
dir <- "V:/VolumeQ/AGteam/NOMIS/census2011/KeyStatistics/"
KS403EW <- read.csv(paste0(dir, "KS403EW.csv"))
heating <- data.frame(LSOA11CD=KS403EW$geography.code,
  heating=KS403EW[[7]]/KS403EW[[5]])

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

################################################################################
# MERGE, AND AGGREGATE

# MERGE
lsoadata <- merge(lookup[c("LSOA11CD","LAD11CD")], pop, sort=F)
for(x in list(popdens, imd, heating, albedo, evi, imperv, agebuild))
  lsoadata <- merge(lsoadata, x, sort=F)

################################################################################
# RURAL-URBAN CLASSIFICATION

dir <- "V:/VolumeQ/AGteam/ONS/rural_urban/"
file <- paste0("Rural_Urban_Classification_(2011)_of_Lower_Layer_",
  "Super_Output_Areas_in_England_and_Wales.csv")
ruralurban <- read.csv(paste0(dir, file))[c(1,3,4)]
names(ruralurban)[1] <- "LSOA11CD"
ruralurban <- merge(lookup["LSOA11CD"], ruralurban, sort=F)

################################################################################
# SHAPEFILES

# LOAD LSOA SHAPEFILES
source <- "V:/VolumeQ/AGteam/ONS/geography/shapefiles"
file <- "Lower_Layer_Super_Output_Area__December_2011__EW_BSC_V2"
file.copy(paste0(source, "/LSOA_EnglandWales/", file, "-shp.zip"), getwd())
unzip(zipfile=paste0(file,"-shp.zip"), exdir=getwd())
lsoashp <- st_read(paste0(file, ".shp"))[2]
file.remove(list.files()[grep(file, list.files())])
lsoashp <- lsoashp[match(lookup$LSOA11CD, lsoashp$LSOA11CD),]

# LOAD LAD SHAPEFILES
source <- "V:/VolumeQ/AGteam/ONS/geography/shapefiles/LAD_EnglandWales/"
file <- "Local_Authority_Districts__December_2011__Boundaries_EW_BGC"
file.copy(paste0(source, file, "-shp.zip"), getwd())
unzip(zipfile=paste0(file,"-shp.zip"), exdir=getwd())
ladshp <- st_read(paste0(file, ".shp"))[1]
file.remove(list.files()[grep(file, list.files())])
names(ladshp)[1] <- "LAD11CD"
ladshp <- ladshp[match(unique(lookup$LAD11CD), ladshp$LAD11CD),]

################################################################################
# CLEAN, SAVE

# CLEAN
rm(list=ls()[!ls() %in% c("lookup","lsoadata","ruralurban","lsoashp","ladshp")])

# SAVE
save.image("temp/prep.RData")
