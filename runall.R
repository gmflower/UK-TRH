################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# RUN ALL THE SCRIPTS
################################################################################

listscript <- c("00.pkg.R", "01.param.R", "02.prepmain.R", "03.firststage.R",
  "04.prepmeta.R", "05.secondstage.R", "06.effects.R", "07.plot.R", 
  "08.tables.R", "09.maps.R", "10.resadd.R")
sapply(listscript, source, echo=T)
