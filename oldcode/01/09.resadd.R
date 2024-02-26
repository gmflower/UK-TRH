################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# RESULTS
################################################################################

# TOTAL DEATHS INCLUDED IN THE ANALYSIS AND PER YEAR
round(sum(lsoares$totdeath)*(sum(date>=datestart)/365.25))
round(sum(lsoares$totdeath))

# RANGE OF CORRELATION BETWEEN LSOA-LEVEL VARIABLES (AT LAD AND LSOA LEVEL)
range(corvarlad[lower.tri(corvarlad)])
range(corvarlsoa[lower.tri(corvarlsoa)])

# VARIANCE EXPLAINED BY THE FIRST THREE COMPONENTS
cumsum(pca$sdev^2 / sum(pca$sdev^2))[3]

# MMP
tmeanper[which.min(bvar%*%predparall$fit)]

# TESTS ON THE META-PREDICTORS AND ON HETEROGENEITY
drop1(metaall, test="Chisq")
qtest(metaall)
summary(metaall)

# RANGE OF RR AND MMT
range(lsoares$RR99)
range(lsoares$RR01)
range(lsoares$mmt)

# NUMBER OF EXCESS DEATHS AND FRACTION
round(allres)
round(allres[,-1]/allres[,1]*100, 2)

# HEAT-RELATED EXCESS MORTALITY IN LONDON
regres
regres[7,7:9] / regres[7,3] * 100
