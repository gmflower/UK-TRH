################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# OUTPUT
################################################################################

################################################################################
# TABLES OF RR AND MMT/MMP AT LSOA LEVEL

# RR
tabrrlsoa <- rbind(merge(unique(lookup[c("LSOA11CD","LSOA11NM")]), rrlsoa))
tabrrlsoa <- tabrrlsoa[with(tabrrlsoa, order(LSOA11CD, agegr, effect)),]
# for(i in 5:7) 
#   tabrrlsoa[,i] <- formatC(tabrrlsoa[,i], digits=3, format="f")
names(tabrrlsoa)[1:5] <- c("LSOA code", "LSOA name","Age group", "Effect", "RR")

# MMT/MMP
tabmmtmmplsoa <- rbind(merge(unique(lookup[c("LSOA11CD","LSOA11NM")]), mmtmmplsoa))
tabmmtmmplsoa <- tabmmtmmplsoa[with(tabmmtmmplsoa, order(LSOA11CD, agegr)),]
#tabmmtmmplsoa$mmt <- formatC(tabmmtmmplsoa$mmt, digits=1, format="f")
names(tabmmtmmplsoa) <- c("LSOA code", "LSOA name", "Age group", "MMT", "MMP")

# SAVE
tablist <- list("RR"=tabrrlsoa, "MMT&MMP"=tabmmtmmplsoa)
write.xlsx(tablist, file = "output/rrmmtmmp.xlsx", overwrite=T)

# CLEAN
rm(tabrrlsoa, tabmmtmmplsoa)

################################################################################
# TABLES OF EXCESS DEATHS FOR DIFFERENT AGGREGATION LEVELS

# WHOLE AREA
tabexcall <- as.data.frame(anall)
# for(i in c(4,6:8)-1) 
#   tabexcall[,i] <- formatC(tabexcall[,i], digits=1, format="f")

# REGION
tabexcreg <- merge(unique(lookup[c("RGN11CD","RGN11NM")]), anreg)
tabexcreg <- tabexcreg[with(tabexcreg, order(RGN11CD, agegr, effect)),]
# for(i in c(4,6:8)) 
#   tabexcreg[,i] <- formatC(tabexcreg[,i], digits=1, format="f")
names(tabexcreg)[1:2] <- c("Region code", "Region name")

# LAD
tabexclad <- rbind(merge(unique(lookup[c("LAD11CD","LAD11NM")]), anlad))
tabexclad <- tabexclad[with(tabexclad, order(LAD11CD, agegr, effect)),]
# for(i in c(4,6:8)+1) 
#   tabexclad[,i] <- formatC(tabexclad[,i], digits=1, format="f")
names(tabexclad)[1:2] <- c("LAD code", "LAD name")

# LSOA
tabexclsoa <- rbind(merge(unique(lookup[c("LSOA11CD","LSOA11NM")]), anlsoa))
tabexclsoa <- tabexclsoa[with(tabexclsoa, order(LSOA11CD, agegr, effect)),]
# for(i in c(4,6:8)+1) 
#   tabexclsoa[,i] <- formatC(tabexclsoa[,i], digits=1, format="f")
names(tabexclsoa)[1:2] <- c("LSOA code", "LSOA name")

# IMD QUINTILES
tabexcimd <- as.data.frame(animd)
# for(i in c(4,6:8)) 
#   tabexcimd[,i] <- formatC(tabexcimd[,i], digits=1, format="f")
names(tabexcimd)[1] <- "IMD quintile"

# RURAL/URBAN
tabexcrururb <- as.data.frame(anrururb)
# for(i in c(4,6:8)) 
#   tabexcrururb[,i] <- formatC(tabexcrururb[,i], digits=1, format="f")
names(tabexcrururb)[1] <- "Rural/urban classification"

# RENAME COLUMNS
names(tabexcall)[2:6-1] <- names(tabexcreg)[2:6+1] <- names(tabexclad)[2:6+1] <-
  names(tabexclsoa)[2:6+1] <- names(tabexcimd)[2:6] <- names(tabexcrururb)[2:6] <- 
  c("Age group", "Effect", "Total deaths", "Population", "Annual excess deaths")

# SAVE
tablist <- list("EnglandWales"=tabexcall, Region=tabexcreg, LAD=tabexclad,
  LSOA=tabexclsoa, IMD=tabexcimd, "RuralUrban"=tabexcrururb)
write.xlsx(tablist, file = "output/excdeath.xlsx", overwrite=T)

# CLEAN
rm(tabexcall, tabexcreg, tabexclad, tabexclsoa, tabexcrururb, tabexcimd)

################################################################################
# TABLES OF STANDARDISED RATE FOR DIFFERENT AGGREGATION LEVELS

# WHOLE AREA
tabrateall <- as.data.frame(stdrateall)
# for(i in c(3:5)-1) 
#   tabrateall[,i] <- formatC(tabrateall[,i], digits=1, format="f")

# REGION
tabratereg <- merge(unique(lookup[c("RGN11CD","RGN11NM")]), stdratereg)
tabratereg <- tabratereg[with(tabratereg, order(RGN11CD, effect)),]
# for(i in c(3:5)) 
#   tabratereg[,i] <- formatC(tabratereg[,i], digits=1, format="f")
names(tabratereg)[1:2] <- c("Region code", "Region name")

# LAD
tabratelad <- rbind(merge(unique(lookup[c("LAD11CD","LAD11NM")]), stdratelad))
tabratelad <- tabratelad[with(tabratelad, order(LAD11CD, effect)),]
# for(i in c(3:5)+1) 
#   tabratelad[,i] <- formatC(tabratelad[,i], digits=1, format="f")
names(tabratelad)[1:2] <- c("LAD code", "LAD name")

# LSOA
tabratelsoa <- rbind(merge(unique(lookup[c("LSOA11CD","LSOA11NM")]), stdratelsoa))
tabratelsoa <- tabratelsoa[with(tabratelsoa, order(LSOA11CD, effect)),]
# for(i in c(3:5)+1) 
#   tabratelsoa[,i] <- formatC(tabratelsoa[,i], digits=1, format="f")
names(tabratelsoa)[1:2] <- c("LSOA code", "LSOA name")

# IMD QUINTILES
tabrateimd <- as.data.frame(stdrateimd)
# for(i in c(3:5)) 
#   tabrateimd[,i] <- formatC(tabrateimd[,i], digits=1, format="f")
names(tabrateimd)[1] <- "IMD quintile"

# RURAL/URBAN
tabraterururb <- as.data.frame(stdraterururb)
# for(i in c(3:5)) 
#   tabraterururb[,i] <- formatC(tabraterururb[,i], digits=1, format="f")
names(tabraterururb)[1] <- "Rural/urban classification"

# RENAME COLUMNS
names(tabrateall)[2:3-1] <- names(tabratereg)[2:3+1] <- names(tabratelad)[2:3+1] <-
  names(tabratelsoa)[2:3+1] <- names(tabrateimd)[2:3] <- names(tabraterururb)[2:3] <- 
  c("Effect", "Excess std rate (x 100,000)")

# SAVE
tablist <- list("EnglandWales"=tabrateall, Region=tabratereg, LAD=tabratelad,
  LSOA=tabratelsoa, IMD=tabrateimd, "RuralUrban"=tabraterururb)
write.xlsx(tablist, file = "output/excstdrate.xlsx", overwrite=T)

# CLEAN
rm(tabrateall, tabratereg, tabratelad, tabratelsoa, tabraterururb, tabrateimd)

