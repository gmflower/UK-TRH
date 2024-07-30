################################################################################
# UK-TRH: SMALL-AREA ANALYSIS OF TEMPERATURE RELATED HOSPITALISATIONS IN ENGLAND
################################################################################

################################################################################
# PREPARE THE META-VARIABLES
################################################################################

################################################################################
# PREPARE THE DATA AT LSOA AND LAD LEVEL

# MERGE META-VARIABLES AND ERASE THEM
lsoadata <- list(lookup[,c("LSOA11CD","LAD11CD")], ruralurban, poptot65, 
  popdens, imd, heating, albedo, evi, imperv, agebuild, tmeansumm) |> 
  Reduce(merge, x=_)
rm(ruralurban, poptot65, popdens, imd, heating, albedo, evi, imperv, agebuild,
  tmeansumm)

# AGGREGATE BY LAD (EXCLUDE RURAL/URBAN, SUM POP, AVERAGE REST WEIGHTED BY POP)
temp <- lsoadata %>%
  group_by(LAD11CD) %>%
  summarise_at(vars(pop65:rangetmean), list(~weighted.mean(., poptot)))
laddata <- lsoadata %>%
  group_by(LAD11CD) %>%
  summarise(poptot=sum(poptot)) %>%
  merge(temp) %>% as.data.frame()
rm(temp)

################################################################################
# LIST OF META-VARIABLES USED FOR MODELLING

# CREATE A LIST OF METAVAR NAMES AND LABELS
metanames <- c("pop65","popdens","IMDincome","IMDemploy","IMDeduc","IMDhealth",
  "IMDhousing","heating","agebuild","IMDenviron","albedo","imperv", "evi",
  "meantmean","rangetmean")
metalab <- c("Population\nabove 65","Population\ndensity","Income","Employment",
  "Education","Health &\ndisability score", "Access to\nhousing & services",
  "Heating\nat home","Age of\nbuildings","Living\nenvironment\nscore","Albedo",
  "Impervious\nsurfaces","Enhanced\nvegetation index","Average\ntemperature",
  "Temperature\nrange")

# CREATING GROUPING
dfmetagp <- data.frame(var=metanames, vargroup=c(rep("Demographic", 2), 
  rep("Socio-economic", 3), "Health & disability",
  rep("Housing & neighbourhood", 4), rep("Landscape", 3),
  rep("Climatological", 2)))
