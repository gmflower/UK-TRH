################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# RESULTS
################################################################################

# TOTAL DEATHS INCLUDED IN THE ANALYSIS AND PER YEAR
round(sum(subset(anlsoa, effect=="tot" & agegr=="all")$totdeath) *
    (length(seqyear)))
round(sum(subset(anlsoa, effect=="tot" & agegr=="all")$totdeath))

# MMP ACROSS AGE GROUPS
for(i in seq(predparage)) print(tmeanper[which.min(bvar%*%predparage[[i]]$fit)])

# RR AT 1ST AND 99TH PERCENTILES FOR 0-64 AND 85PLUS
cp064 <- crosspred(bvar, coef=predparage[[1]]$fit, vcov=predparage[[1]]$vcov,
  model.link="log", cen=tmeanper[which.min(bvar%*%predparage[[1]]$fit)],
  at=tmeanper[c("1.0%","99.0%")])
cp85plus <- crosspred(bvar, coef=predparage[[4]]$fit, vcov=predparage[[4]]$vcov,
  model.link="log", cen=tmeanper[which.min(bvar%*%predparage[[4]]$fit)],
  at=tmeanper[c("1.0%","99.0%")])
cp064[c("allRRfit","allRRlow","allRRhigh")]
cp85plus[c("allRRfit","allRRlow","allRRhigh")]

# TESTS ON THE META-PREDICTORS AND ON HETEROGENEITY
drop1(metaall, test="Chisq")
qtest(metaall)
print(summary(metaall), digits=4)

# RANGE OF CORRELATION BETWEEN LSOA-LEVEL VARIABLES (AT LAD AND LSOA LEVEL)
range(corvarlad[lower.tri(corvarlad)])
range(corvarlsoa[lower.tri(corvarlsoa)])

# VARIANCE EXPLAINED BY THE FIRST THREE COMPONENTS
cumsum(pca$sdev^2 / sum(pca$sdev^2))[ncomp]

# LR TEST OF PCA INDICATORS
# NB: TRICK TO TEST ALL INDICATORS TOGETHER (anova STILL NOT AVAILABLE)
drop1(update(metaall, .~agegr+as.matrix(ladvar[,3:5])), test="Chisq")

# NUMBER OF EXCESS DEATHS AND RATE
round(subset(anall, agegr=="all")[,5:7])
round(as.matrix(stdrateall[,-1]), 2)

# HEAT-RELATED EXCESS MORTALITY IN LONDON
subset(stdratereg, effect=="heat" & RGN11CD=="E12000007")

# RANGE AND IQR OF RATES AND MMT/MMP
quantile(subset(stdratelsoa, effect=="heat")$stdrate)
quantile(subset(stdratelsoa, effect=="cold")$stdrate)
quantile(mmtmmplsoa$mmt)
quantile(mmtmmplsoa$mmp)

