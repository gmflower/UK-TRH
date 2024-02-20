################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# PREPARE THE META-VARIABLES
################################################################################

################################################################################
# CENSUS VARIABLES

# % OF POP 65 AND OVER
dir <- "V:/VolumeQ/AGteam/NOMIS/census2011/KeyStatistics/"
KS102EW <- read.csv(paste0(dir, "KS102EW.csv"))
poptot65 <- data.frame(LSOA11CD=KS102EW$geography.code,
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

# MERGE THE META-VARIABLES ABOVE
lsoadata <- merge(lookup[,c("LSOA11CD","LAD11CD")], poptot65)
for(x in list(popdens, imd, heating, albedo, evi, imperv, agebuild))
  lsoadata <- merge(lsoadata, x)

# ADD LSOA TEMPERATURE DATA
lsoadata <-
  merge(lsoadata, do.call(rbind, lapply(stage1list, "[[", "tmeanlsoa")))

################################################################################
# PREPARE THE DATA AT LAD LEVEL

# AGGREAGATE BY LAD (SUM POP, AVERAGE THE REST)
temp <- lsoadata %>%
  group_by(LAD11CD) %>%
  summarise_at(vars(pop65:rangetmean), list(~weighted.mean(., poptot)))
laddata <- lsoadata %>%
  group_by(LAD11CD) %>%
  summarise(poptot=sum(poptot)) %>%
  merge(temp) %>% as.data.frame()

# CREATE A LIST OF SELECTED METAVAR
metanames <- c("pop65","popdens","IMDincome","IMDemploy","IMDeduc","IMDhealth",
  "IMDhousing","heating","agebuild","IMDenviron","albedo","imperv", "evi",
  "meantmean","rangetmean")
metalab <- c("Population\nabove 65","Population\ndensity","Income","Employment",
  "Education","Health &\ndisability score", "Access to\nhousing & services",
  "Heating\nat home","Age of\nbuildings","Living\nenvironment\nscore","Albedo",
  "Impervious\nsurfaces","Enhanced\nvegetation index","Average\ntemperature",
  "Temperature\nrange")
metavar <- laddata[metanames]

# CREATING GROUPING
dfmetagp <- data.frame(var=metanames, vargroup=c(rep("Demographic", 2), 
  rep("Socio-economic", 3), "Health & disability",
  rep("Housing & neighbourhood", 4), rep("Landscape", 3),
  rep("Climatological", 2)))
