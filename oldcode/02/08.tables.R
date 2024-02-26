################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# TABLES
################################################################################

# FUNCTION TO FORMAT VARIABLES WITH INTERVALS
funformat <- function(mat, digits=0) {
  mat <- as.matrix(mat)
  paste0(formatC(mat[,1, drop=T], digits=digits, format="f", big.mark=","), " (",
    formatC(mat[,2, drop=T], digits=digits, format="f", big.mark=",") , " to ",
    formatC(mat[,3, drop=T], digits=digits, format="f", big.mark=","), ")")
}

################################################################################
# DESCRIPTIVE TABLE

# TEMPERATURE STATISTICS
tmeanlad <- t(sapply(stage1list, function(x) c(mean(x$tmeanlsoa$meantmean),
  range(x$tmeanper))))
tmeanlad <- as.data.frame(tmeanlad)
names(tmeanlad) <- c("meantmean", "mintmean", "maxtmean")
tmeanlad$LAD11CD <- names(stage1list)

# TABLE WITH SUMMARY STATISTICS PER REGION
tabdesc <- lookup %>% 
  merge(tmeanlad) %>% merge(anlsoa[agegr=="all" & effect=="tot"]) %>%
  group_by(RGN11CD,RGN11NM) %>%
  summarise(LAD=length(unique(LAD11CD)), 
    LSOA=formatC(length(unique(LSOA11CD)), digits=0, format="f"),
    totdeath=formatC(round(sum(totdeath)), digits=0, format="f"),
    meantmean=mean(meantmean), mintmean=min(mintmean), maxtmean=max(maxtmean)) %>%
  as.data.frame()

# ADD TOTAL
totstat <- lookup %>% 
  merge(tmeanlad) %>% merge(anlsoa[agegr=="all" & effect=="tot"]) %>%
  summarise(LAD=length(unique(LAD11CD)),
    LSOA=formatC(length(unique(LSOA11CD)), digits=0, format="f"),
    totdeath=formatC(round(sum(totdeath)), digits=0, format="f"),
    meantmean=mean(meantmean), mintmean=min(mintmean), maxtmean=max(maxtmean))
tabdesc <- rbind(tabdesc[-1], cbind(RGN11NM="England & Wales", totstat))

# PRINT
names(tabdesc) <- c(" ", "LAD (N)", "LSOA (N)", "Annual deaths")
ftmean <- funformat(tabdesc[,5:7], digits=1)
write.csv(cbind(tabdesc[1:4], "Temperature (C)"=ftmean), row.names=F,
  file="tables/tabdesc.csv")

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

################################################################################
# TABLE OF EXCESS MORTALITY/STANDARDISED RATE BY REGION

# EXCESS MORTALITY AND STANDARDISED RATE
tabexcreg <- cbind(anreg[agegr=="all"&effect!="tot", c(1,3)],
  exc=funformat(anreg[agegr=="all"&effect!="tot",6:8], digits=0),
  stdrate=funformat(stdratereg[effect!="tot",3:5], digits=2)) %>%
  pivot_wider(names_from=2, values_from=3:4)
tabexcreg <- merge(unique(lookup[c("RGN11CD","RGN11NM")]), tabexcreg)[-1]

# ADD ALL
tabexcreg <- rbind(tabexcreg, c("England & Wales", 
  funformat(anall[agegr=="all"&effect!="tot", 5:7], digits=0),
  funformat(stdrateall[, 2:4], digits=2)))

# SAVE
write.csv(tabexcreg, row.names=F, file="tables/tabexcreg.csv")

################################################################################
# TABLE OF SUMMARY DISTRIBUTION OF EFFECT SUMMARIES BY LSOA

rrtab <- rrlsoa %>% 
  group_by(effect, agegr) %>% 
  summarise(as.data.frame(t(c(mean(rr), quantile(rr)))))
rrtab[3:8] <- formatC(as.matrix(rrtab[3:8]), digits=3, format="f")

antab <- anlsoa %>% 
  filter(totdeath>0 & pop>0 & effect!="tot") %>%
  mutate(af=excdeath/totdeath*100) %>%
  group_by(effect, agegr) %>% 
  summarise(as.data.frame(t(c(mean(af), quantile(af)))))
antab[3:8] <- formatC(as.matrix(antab[3:8]), digits=2, format="f")

mmtmmptab <- mmtmmplsoa %>% 
  pivot_longer(3:4, names_to="effect") %>%
  group_by(effect, agegr) %>% 
  summarise(as.data.frame(t(c(mean(value), quantile(value)))))
mmtmmptab[3:8] <- formatC(as.matrix(mmtmmptab[3:8]), digits=1, format="f")

stdratetab <- stdratelsoa %>% 
  filter(effect!="tot") %>%
  group_by(effect) %>% 
  summarise(as.data.frame(t(c(mean(stdrate), quantile(stdrate)))))
stdratetab[2:7] <- formatC(as.matrix(stdratetab[2:7]), digits=2, format="f")
stdratetab$agegr <- ""

tabeffsum <- cbind(stat=rep(c("RR","AF","MMT/MMP","Rate"),
  c(nrow(rrtab), nrow(antab), nrow(mmtmmptab), nrow(stdratetab))),
  rbind(rrtab, antab, mmtmmptab, stdratetab))
names(tabeffsum)[4:9] <- c("Mean", "Min", "25th", "Median", "75th", "Max")

# SAVE
write.csv(tabeffsum, row.names=F, file="tables/tabeffsum.csv")

# CLEAN
rm(rrtab, antab, mmtmmptab, stdratetab)

################################################################################
# TABLE OF EFFECT MODIFICATION (RATIO OF RELATIVE RISKS)

tabrrr <- pivot_wider(cbind(rrrvar[2:3], rrr=funformat(rrrvar[4:6], 3)),
  names_from=2, values_from=3)

# SAVE
write.csv(tabrrr, row.names=F, file="tables/tabrrr.csv")
