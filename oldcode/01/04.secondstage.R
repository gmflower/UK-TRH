################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# SECOND STAGE
################################################################################

################################################################################
# PREPARE THE DATA AT LAD LEVEL

# ADD LSOA TEMPERATURE DATA
lsoadata <- merge(lsoadata, do.call(rbind, lapply(estlist, "[[", "tmeanlsoa")),
  sort=F)

# AGGREAGATE BY LAD (SUM POP, AVERAGE THE REST)
temp <- lsoadata %>%
  group_by(LAD11CD) %>%
  summarise_at(vars(pop65:rangetmean), list(~weighted.mean(., poptot)))
laddata <- lsoadata %>%
  group_by(LAD11CD) %>%
  summarise(poptot=sum(poptot)) %>%
  merge(temp) %>% as.data.frame()

# CREATE A LIST OF SELECTED METAVAR
metanames <- c("pop65","IMDincome","IMDemploy","IMDeduc","IMDhealth",
  "IMDhousing","IMDenviron","popdens","heating","agebuild","albedo","imperv",
  "evi","meantmean","rangetmean")
metalab <- c("Population\nabove 65","Income","Employment","Education","Health\n& disability",
  "Barriers\nto housing","Living\nenvironment","Population\ndensity","Heating\nat home","Age of\nbuildings",
  "Albedo","Impervious\nsurfaces","Enhanced\nvegetation index","Average\ntemperature",
  "Temperature\nrange")
metavar <- laddata[metanames]

################################################################################
# PCA AND META-ANALYSIS

# RUN PCA
pca <- prcomp(metavar, center=TRUE, scale=TRUE)
summary(pca)

# EXTRACT THE FIRST THREE COMPONENTS AS META-PREDICTORS
ladvar <- as.data.frame(pca$x[,1:3], row.names=listlad)

# EXTRACT THE PARAMETERS
coefall <- t(sapply(estlist, function(x) x$coefall))
vcovall <- lapply(estlist, function(x) x$vcovall)

# META-ANALYTICAL MODEL FOR OVERALL CUMULATIVE (WITH FIXD-EFFECTS)
metaall <- mixmeta(coefall~., S=vcovall, method="fixed", data=ladvar)
summary(metaall)

################################################################################
# SAVE

save.image("temp/secondstage.RData")
