################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# EFFECTS
################################################################################

# PREPARE THE PARALLELIZATION
ncores <- detectCores()
cl <- parallel::makeCluster(max(1,ncores-1))
registerDoParallel(cl)

# PACKAGE LIST FOR PARALLELIZATION
pack <- c("dlnm", "data.table", "tsModel", "MASS")

################################################################################
# BY LAD

# EXTRACT COEF/VCOV AT LAD LEVEL 
ladpar <- predict(metaall, vcov=T)

# RUN THE LOOP
ladres <- foreach(i=seq(length(listlad)), .packages=pack, .combine=rbind) %dopar% {

  # COMPUTE TEMPERATURE PERCENTILES
  tmeanper <- estlist[[i]]$tmeanper

  # DEFINE PARAMETERS OF THE EXPOSURE-RESPONSE FUNCTION
  argvar <- list(fun=varfun, knots=tmeanper[c("10.0%","75.0%","90.0%")],
    Bound=tmeanper[c("0.0%","100.0%")])
  argvar$degree <- vardegree
  
  # IDENTIFY MMT/MMP (BETWEEN 1ST AND 99TH)
  indextr <- which(predper<1 | predper>99)
  bvar <- do.call(onebasis, c(list(x=tmeanper[-indextr]), argvar))
  mmt <- tmeanper[-indextr][[which.min(bvar%*%ladpar[[i]]$fit)]]
  mmp <- predper[[which(tmeanper%in%mmt)]] 

  # COMPUTE RR
  cp <- crosspred(bvar, coef=ladpar[[i]]$fit, vcov=ladpar[[i]]$vcov,
    model.link="log", cen=mmt, at=tmeanper[c("1.0%","99.0%")])
  rr <- c(rbind(cp$allRRfit, cp$allRRlow, cp$allRRhigh))
  names(rr) <- c(t(outer(c("RR01","RR99"), c("","_95%CIlow","_95%CIhigh"),
    paste0)))

  # RESULTS
  cbind(data.frame(LAD11CD=listlad[i], LAD11NM=unique(lookup$LAD11NM[i]), mmt=mmt,
    mmp=mmp), t(rr))
}

################################################################################
# BY LSOA

# DEFINE THE PREDICTORS AT LSOA LEVEL
lsoavar <- as.data.frame(predict(pca, newdata=lsoadata[metanames])[,1:3],
  row.names=listlsoa)

# EXTRACT COEF/VCOV AT LSOA LEVEL (USE newdata TO ASSIGN NAMES)
# NB: USED TO QUICKLY COMPUTE RR AND POINT ESTIMATE OF EXCESS DEATHS
lsoapar <- predict(metaall, newdata=lsoavar, vcov=T)

# NUMBER OF SIMULATION FOR eCI OF EXCESS DEATHS
nsim <- 100

# TO COMPUTE eCI OF EXCESS DEATHS ACCOUNTING FOR CORRELATIONS ACROSS LSOA:
# - RECONSTRUCT THE MODEL MATRIX OF THE META-REGRESSION AT LSOA LEVEL
# - SAMPLE THE COEF OF THE META-REGRESSION
Xdeslsoa <- model.matrix(formula(metaall)[c(1,3)], data=lsoavar)
set.seed(13041975)
mvcoefsim <- mvrnorm(nsim, coef(metaall), vcov(metaall))

# RUN THE LOOP
lsoareslist <- foreach(i=seq(length(listlsoa)), .packages=pack) %dopar% {
  
  # EXTRACT THE LSOA DATA, SELECT FROM YEAR 2000 ON
  data <- data.table(readRDS(paste0(tmeanpath, listlsoa[i], ".Rds")))
  data$date <- date
  data <- data[data$date>=datestart,]
  
  # COMPUTE TEMPERATURE PERCENTILES
  tmeanper <- quantile(data$tmean, predper/100, na.rm=T)
  
  # DEFINE PARAMETERS OF THE EXPOSURE-RESPONSE FUNCTION
  argvar <- list(fun=varfun, knots=tmeanper[c("10.0%","75.0%","90.0%")],
    Bound=tmeanper[c("0.0%","100.0%")])
  argvar$degree <- vardegree
  
  # IDENTIFY MMT/MMP (BETWEEN 1ST AND 99TH)
  indextr <- which(predper<1 | predper>99)
  bvar <- do.call(onebasis, c(list(x=tmeanper[-indextr]), argvar))
  mmt <- tmeanper[-indextr][[which.min(bvar%*%lsoapar[[i]]$fit)]]
  mmp <- predper[[which(tmeanper%in%mmt)]] 
  
  # COMPUTE RR
  cp <- crosspred(bvar, coef=lsoapar[[i]]$fit, vcov=lsoapar[[i]]$vcov,
    model.link="log", cen=mmt, at=tmeanper[c("1.0%","99.0%")])
  rr <- c(rbind(cp$allRRfit, cp$allRRlow, cp$allRRhigh))
  names(rr) <- c(t(outer(c("RR01","RR99"), c("","_95%CIlow","_95%CIhigh"),
    paste0)))
  
  # DERIVE THE CENTERED BASIS
  bvar <- do.call(onebasis, c(list(x=data$tmean), argvar))
  cenvec <- do.call(onebasis, c(list(x=mmt), argvar))
  bvarcen <- scale(bvar, center=cenvec, scale=F) 
  
  # FORWARD MOVING AVERAGE OF DEATHS
  deaths <- rowMeans(as.matrix(Lag(data$count, -21:0)))
  
  # INDICATOR FOR COLD/HEAT DAYS
  indheat <- data$tmean>=mmt
  
  # COMPUTE THE YEARLY EXCESS DEATHS
  anday <- (1-exp(-bvarcen%*%lsoapar[[i]]$fit))*deaths
  an <- c(ANtot=sum(anday, na.rm=T), ANheat=sum(anday[indheat], na.rm=T),
    ANcold=sum(anday[!indheat], na.rm=T)) / (sum(date>=datestart)/365.25)
  
  # SIMULATED DISTRIBUTION OF YEARLY EXCESS DEATHS
  # - EXPAND THE DESING MATRIX AND MULTIPLE WITH SAMPLES COEF
  # - COMPUTE THE EXCESS DEATHS
  set.seed(13041975+i)
  ansim <- sapply(seq(nsim), function(j) {
    coef <- drop((Xdeslsoa[i,,drop=FALSE] %x% diag(5)) %*% mvcoefsim[j,])
    anday <- (1-exp(-bvarcen%*%coef))*deaths
    c(tot=sum(anday, na.rm=T), heat=sum(anday[indheat], na.rm=T), 
      cold=sum(anday[!indheat], na.rm=T))
  }) / (sum(date>=datestart)/365.25)
  colnames(ansim) <- paste0("sim", seq(nsim))

  # TOTAL YEARLY DEATHS
  totdeath <- sum(deaths, na.rm=T) / (sum(date>=datestart)/365.25)
  
  # MAIN RESULTS
  res <- cbind(data.frame(LSOA11CD=listlsoa[i], LSOA11NM=lookup$LSOA11NM[i],
    mmt=mmt, mmp=mmp), t(rr), totdeath=totdeath, t(an))

  # RETURN
  list(res=res, ansim=ansim)
}

# REMOVE PARALLELIZATION
stopCluster(cl)

################################################################################
# TABLE BY LSOA

# EXTRACT lsoares
lsoares <- do.call(rbind, lapply(lsoareslist, "[[", "res"))

# EXTRACT ansim, BIND, AND PERMUTE ITS DIMENSIONS
ansim <- aperm(abind(lapply(lsoareslist, "[[", "ansim"), along=3), c(3,1,2))
dimnames(ansim)[[1]] <- listlsoa

# COMPUTE eCI
anci <- matrix(c(aperm(apply(ansim, 1:2, quantile, c(2.5,97.5)/100), c(1,3,2))),
  nrow=length(listlsoa), byrow=T)
colnames(anci) <- c(t(outer(c("ANtot","ANheat","ANcold"),
  c("_95%eCIlow" ,"_95%eCIhigh"), paste0)))

# PUT TOGETHER AND REORDER
lsoares <- cbind(lsoares, anci) %>%
  dplyr::select(-starts_with("AN"), starts_with("ANtot"), starts_with("ANheat"),
    starts_with("ANcold"))

# ERASE OBJ
rm(lsoareslist, anci)

################################################################################
# TABLE BY LAD

# AGGREGATE AN BY LAD AND MERGE
anlad <-  cbind(lsoares, LAD11CD=lookup$LAD11CD) %>%
  group_by(LAD11CD) %>% 
  summarise_at(vars(totdeath,ANtot, ANheat, ANcold), sum)
ladres <- merge(ladres, anlad)

# COMPUTE eCI
anci <- do.call(cbind, lapply(c("tot","heat","cold"), function(x) {
  aggr <- apply(ansim[,x,], 2, tapply, lookup$LAD11CD, sum)
  t(apply(aggr[,-1], 1, quantile, c(2.5,97.5)/100))
}))
dimnames(anci) <- list(NULL, c(t(outer(c("ANtot","ANheat","ANcold"),
  c("_95%eCIlow" ,"_95%eCIhigh"), paste0))))

# PUT TOGETHER AND REORDER
ladres <- cbind(ladres, anci) %>%
  dplyr::select(-starts_with("AN"), starts_with("ANtot"), starts_with("ANheat"),
    starts_with("ANcold"))

# ERASE OBJ
rm(anlad, anci)

################################################################################
# TABLE BY REGION

# AGGREGATE AN BY REG
anreg <-  cbind(lsoares, lookup[c("RGN11CD","RGN11NM")]) %>%
  group_by(RGN11CD, RGN11NM) %>% 
  summarise_at(vars(totdeath,ANtot, ANheat, ANcold), sum) %>%
  as.data.frame()

# COMPUTE eCI
anci <- do.call(cbind, lapply(c("tot","heat","cold"), function(x) {
  aggr <- apply(ansim[,x,], 2, tapply, lookup$RGN11CD, sum)
  t(apply(aggr[,-1], 1, quantile, c(2.5,97.5)/100))
}))
dimnames(anci) <- list(NULL, c(t(outer(c("ANtot","ANheat","ANcold"),
  c("_95%eCIlow" ,"_95%eCIhigh"), paste0))))

# PUT TOGETHER AND REORDER
regres <- cbind(anreg, anci) %>%
  dplyr::select(-starts_with("AN"), starts_with("ANtot"), starts_with("ANheat"),
    starts_with("ANcold"))

# ERASE OBJ
rm(anreg, anci)

################################################################################
# TABLE BY DEPRIVATION QUINTILES

# FACTOR FOR IMD QUINTILES
fimd <- cut(lsoadata$IMD, 5, labels=paste0("IMDquint",1:5))

# AGGREGATE AN BY IMD
animd <-  cbind(lsoares, IMD=fimd) %>%
  group_by(IMD) %>% 
  summarise_at(vars(totdeath,ANtot, ANheat, ANcold), sum) %>%
  as.data.frame()

# COMPUTE eCI
anci <- do.call(cbind, lapply(c("tot","heat","cold"), function(x) {
  aggr <- apply(ansim[,x,], 2, tapply, fimd, sum)
  t(apply(aggr[,-1], 1, quantile, c(2.5,97.5)/100))
}))
dimnames(anci) <- list(NULL, c(t(outer(c("ANtot","ANheat","ANcold"),
  c("_95%eCIlow" ,"_95%eCIhigh"), paste0))))

# PUT TOGETHER AND REORDER
imdres <- cbind(animd, anci) %>%
  dplyr::select(-starts_with("AN"), starts_with("ANtot"), starts_with("ANheat"),
    starts_with("ANcold"))

# ERASE OBJ
rm(animd, anci)

################################################################################
# TABLE BY RURAL-URBAN CLASSIFICATION

# DEFINE THE CATEGORIES
ruralurban$frururb <- factor(substr(ruralurban$RUC11CD, 1, 1), 
  labels=c("Major conurbation", "Minor conurbation", "Urban city and town",
    "Rural town", "Rural village"))

# AGGREGATE AN BY RURAL-URBAN
anrururb <-  merge(lsoares, ruralurban) %>%
  group_by(RurUrb=frururb) %>% 
  summarise_at(vars(totdeath,ANtot, ANheat, ANcold), sum) %>%
  as.data.frame()

# COMPUTE eCI
anci <- do.call(cbind, lapply(c("tot","heat","cold"), function(x) {
  aggr <- apply(ansim[,x,], 2, tapply, ruralurban$frururb, sum)
  t(apply(aggr[,-1], 1, quantile, c(2.5,97.5)/100))
}))
dimnames(anci) <- list(NULL, c(t(outer(c("ANtot","ANheat","ANcold"),
  c("_95%eCIlow" ,"_95%eCIhigh"), paste0))))

# PUT TOGETHER AND REORDER
rururbres <- cbind(anrururb, anci) %>%
  dplyr::select(-starts_with("AN"), starts_with("ANtot"), starts_with("ANheat"),
    starts_with("ANcold"))

# ERASE OBJ
rm(anrururb, anci)

################################################################################
# TABLE OVERALL

# AGGREGATE ALL
anall <-  lsoares %>%
  summarise_at(vars(totdeath,ANtot, ANheat, ANcold), sum) %>%
  as.data.frame()

# COMPUTE eCI
anci <- t(c(apply(apply(ansim, 2:3, sum), 1, quantile, c(2.5,97.5)/100)))
dimnames(anci) <- list(NULL, c(t(outer(c("ANtot","ANheat","ANcold"),
  c("_95%eCIlow" ,"_95%eCIhigh"), paste0))))

# PUT TOGETHER AND REORDER
allres <- cbind(anall, anci) %>%
  dplyr::select(-starts_with("AN"), starts_with("ANtot"), starts_with("ANheat"),
    starts_with("ANcold"))

# ERASE OBJ
rm(anall, anci)

################################################################################
# CLEAN AND SAVE

# REMOVE ansim
rm(ansim)

# SAVE
save.image("temp/effects.RData")
