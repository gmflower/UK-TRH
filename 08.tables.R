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

# TABLE WITH SUMMARY STATISTICS PER REGION
tabdesc <- lookup |>
  merge(lsoadata[c("LSOA11CD","meantmean")]) |>
  merge(lsoatmeanper[c("LSOA11CD", "0.0%","100.0%")]) |>
  rename(mintmean="0.0%", maxtmean="100.0%") |>
  merge(anlsoa[agegr=="all" & effect=="tot"]) |>
  group_by(RGN11CD,RGN11NM) %>%
  summarise(LAD=length(unique(LAD11CD)), 
    LSOA=formatC(length(unique(LSOA11CD)), digits=0, format="f"),
    totdeath=formatC(round(sum(totdeath)), digits=0, format="f"),
    meantmean=mean(meantmean), mintmean=min(mintmean), maxtmean=max(maxtmean)) |>
  as.data.frame()

# ADD TOTAL
totstat <- lookup |>
  merge(lsoadata[c("LSOA11CD","meantmean")]) |>
  merge(lsoatmeanper[c("LSOA11CD", "0.0%","100.0%")]) |>
  rename(mintmean="0.0%", maxtmean="100.0%") |>
  merge(anlsoa[agegr=="all" & effect=="tot"]) |>
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
  funformat(stdrateall[effect!="tot", 2:4], digits=2)))

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
