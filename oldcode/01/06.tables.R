################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# TABLES
################################################################################

################################################################################
# DESCRIPTIVE TABLE

# TEMPERATURE STATISTICS
tmeanlad <- t(sapply(estlist, function(x) c(mean(x$tmeanlsoa$meantmean),
  range(x$tmeanper))))
tmeanlad <- as.data.frame(tmeanlad)
names(tmeanlad) <- c("meantmean", "mintmean", "maxtmean")
tmeanlad$LAD11CD <- names(estlist)

# TABLE WITH SUMMARY STATISTICS PER REGION
tabstat <- lookup %>% 
  merge(tmeanlad) %>% merge(lsoares) %>%
  group_by(RGN11CD,RGN11NM) %>%
  summarise(LAD=length(unique(LAD11CD)), 
    LSOA=formatC(length(unique(LSOA11CD)), digits=0, format="f"),
    totdeath=formatC(round(sum(totdeath)), digits=0, format="f"),
    meantmean=mean(meantmean), mintmean=min(mintmean), maxtmean=max(maxtmean)) %>%
  as.data.frame()

# ADD TOTAL
totstat <- lookup %>% 
  merge(tmeanlad) %>% merge(lsoares) %>%
  summarise(LAD=length(unique(LAD11CD)),
    LSOA=formatC(length(unique(LSOA11CD)), digits=0, format="f"),
    totdeath=formatC(round(sum(totdeath)), digits=0, format="f"),
    meantmean=mean(meantmean), mintmean=min(mintmean), maxtmean=max(maxtmean))
tabstat <- rbind(tabstat[-1], cbind(RGN11NM="England & Wales", totstat))

# PRINT
names(tabstat) <- c(" ", "LAD (N)", "LSOA (N)", "Deaths (per year)")
ftmean <- paste0(formatC(tabstat[,5], digits=1, format="f"), " (",
  formatC(tabstat[,6], digits=1, format="f") , " to ",
  formatC(tabstat[,7], digits=1, format="f"), ")")
write.csv(cbind(tabstat[1:4], "Temperature (C)"=ftmean), row.names=F,
  file="tables/tabstat.csv")
  
################################################################################
# TABLES OF EXCESS MORTALITY FOR DIFFERENT AGGREGATION LEVELS

tabreg <- rbind(regres[,-1], cbind(RGN11NM="England & Wales", allres))
names(tabreg)[1:2] <- c(" ","Deaths (per year)")
for(i in 2:ncol(tabreg)) 
  tabreg[,i] <- formatC(tabreg[,i], digits=1, format="f")
write.csv(tabreg, row.names=F, file="tables/tabreg.csv")

tabimd <- imdres
names(tabimd)[1:2] <- c("IMD (quintiles)","Deaths (per year)")
tabimd[,1] <- paste("Quintile", 1:5)
for(i in 2:ncol(tabimd))
  tabimd[,i] <- formatC(tabimd[,i], digits=1, format="f")
write.csv(tabimd, row.names=F, file="tables/tabimd.csv")

tabrururb <- rururbres
names(tabrururb)[1:2] <- c("Rural/urban classification","Deaths (per year)")
for(i in 2:ncol(tabrururb))
  tabrururb[,i] <- formatC(tabrururb[,i], digits=1, format="f")
write.csv(tabrururb, row.names=F, file="tables/tabrururb.csv")

################################################################################
# ALL THE RESULTS

# NAMES AND WRITE
tablist <- list("EnglandWales"=allres, Region=regres, LAD=ladres, LSOA=lsoares,
  IMD=imdres, "RuralUrban"=rururbres)
write.xlsx(tablist, file = "output/output.xlsx")
