################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# SECOND STAGE
################################################################################

################################################################################
# PCA AND META-ANALYSIS

# RUN THE PCA
pca <- prcomp(metavar, center=T, scale=T)
summary(pca)

# EXTRACT THE COMPONENTS (EXPANDED BY AGE GROUPS)
pcavar <- pca$x[rep(seq(listlad), each=length(agelab)),]

# EXTRACT THE PARAMETERS (BIND AGE-SPECIFIC FOR COEF, UNLIST 1 LEVEL FOR VCOV)
coefall <- do.call(rbind, lapply(stage1list, function(x) 
  t(sapply(x$estlist, "[[", "coefall"))))
vcovall <- unlist(lapply(stage1list, function(x) 
  lapply(x$estlist, "[[", "vcovall")), recursive=F)

# FIND THE OPTIMAL NUMBER OF COMPONENTS FOR THE PCA
# NB: PRE-SET TO 3
ncomp <- which.min(sapply(1:10, function(j) {
  ladvar <- cbind(LAD11CD=rep(listlad, each=length(agelab)),
    agegr=rep(agelab, length(listlad)), as.data.frame(pcavar[,seq(j)]))
  AIC(mixmeta(coefall~., S=vcovall, method="fixed", data=ladvar[,-1]))
}))
ncomp <- 3

# DEFINE THE META-PREDICTORS (EXPANDED BY AGE GROUP)
ladvar <- cbind(LAD11CD=rep(listlad, each=length(agevarlab)),
  agegr=rep(agevarlab, length(listlad)), as.data.frame(pcavar[, seq(ncomp)]))

# META-ANALYTICAL MODEL FOR OVERALL CUMULATIVE (WITH FIXD-EFFECTS)
metaall <- mixmeta(coefall~., S=vcovall, method="fixed", data=ladvar[,-1])
summary(metaall)

