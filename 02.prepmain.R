################################################################################
# UK-TRH: SMALL-AREA ANALYSIS OF TEMPERATURE RELATED HOSPITALISATIONS IN ENGLAND
################################################################################

################################################################################
# PREPARE THE MAIN AND OTHER DATASETS (NON REPRODUCIBLE)
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

# LOAD HOSPITAL DATA AND TRANSFORM
listhes <- lapply(seqyear, function(y) {
  cat("\n", y, "")
  data <- readRDS(paste0(hosppath,"/emergency/emrgHES", y, ".RDS")) |> 
    as.data.table()
  # SELECT CAUSE
  yeardata <- rbindlist(lapply(seq_along(causelist), function(j) {
    out    <- causelist[[j]][[2]]
    icd    <- causelist[[j]][[3]]
    icddig <- causelist[[j]][[4]]
    subset <- data[substr(DIAG_01, 1, icddig) %in% icd]
    # CREATE AGE GROUPS AND REMOVE MISSING/UNDER 18
    subset[, agegr := cut(ADMIAGE, agecut, labels = agevarlab, include.lowest = TRUE)]
    subset <- subset[!is.na(agegr),]
    # AGGREGATE, RENAME AND ADD OUTCOME
    subset <- subset[, .(count = .N), 
                     by = .(LSOA11CD = LSOA11, date = ADMIDATE, agegr)]
    subset[, cause := out]
    # SELECT ONLY SUMMER MONTHS
    subset <- subset[month(date) %in% seqmonth,]
    # RETURN
    subset
  }))
  yeardata
})
hesdata <- rbindlist(listhes)
setkey(hesdata, cause, LSOA11CD, date, agegr)

# EXCLUDE NON-MATCHING LSOA
hesdata <- hesdata[LSOA11CD %in% listlsoa,]

# CREATE ALL AGE GROUP
hes_allages <- hesdata |>
  group_by(LSOA11CD, date, cause) |>
  summarise(count = sum(count)) |>
  as.data.table()
hes_allages$agegr <- "total"
hesdata <- rbind(hesdata, hes_allages)

# CREATE A LIST OF CAUSES
setcause <- unique(hesdata$cause)

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
datatmean <- datatmean[month(date) %in% seqmonth,]
datatmean <- datatmean[LSOA11CD %in% listlsoa,]

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
