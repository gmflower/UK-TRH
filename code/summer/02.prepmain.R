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
hesdata <- as.data.table(readRDS(paste0(hosppath, "emrgcountHES_stacked.RDS")))
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
setcause <- setcause[c(1:6,8:11,13,15,18:19,21,23:25)]

# SET KEYS
setkey(hesdata, LSOA11CD, date)

# LOAD THE TEMPERATURE DATA
#listtmean <- lapply(seqyear, function(y) {
#  cat(y, "")
#  file <- paste0("lsoa_temp_", y, ".RDS")
#  out <- data.table(readRDS(paste(tmeanpath, file, sep="/")))
#  setkey(out, LSOA11CD, date)
#  out
#})
#datatmean <- do.call(rbind, listtmean)
#rm(listtmean)

# LOAD TEMPERATURE (ONLY MEAN), POLLUTION (ONLY PM2.5 AND NO2), AND MERGE
listenv <- lapply(seqyear, function(y) {
  cat(y, "")
  poll <- paste0(envpath, "pollutants/lsoa_no2pm10pm25_", y, ".RDS") |> 
    readRDS() |> as.data.table() |> subset(select=-pm10)
  tmean <- paste0(envpath, "temperature/met_v13/lsoa_temp_", y, ".RDS") |>
    readRDS() |> as.data.table() |> subset(select=LSOA11CD:tmean)
  merge(poll, tmean)
})
datatmean <- do.call(rbind, listenv)
rm(listenv)

# Rolling averages for NO2 and PM2.5:
datatmean$no2mean <- runMean(datatmean$no2, 0:1)
datatmean$pm25mean <- runMean(datatmean$pm25, 0:1)


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

