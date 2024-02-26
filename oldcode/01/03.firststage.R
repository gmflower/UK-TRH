################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# FIRST STAGE
################################################################################

# PREPARE THE PARALLELIZATION
ncores <- detectCores()
cl <- parallel::makeCluster(max(1,ncores-1))
registerDoParallel(cl)

# PACKAGE LIST FOR PARALLELIZATION
pack <- c("dlnm", "data.table", "gnm", "tsModel", "splines")

# WRITE A TEXT FILE TO TRACE ITERATIONS
writeLines(c(""), "temp/logstage1.txt")
cat(as.character(as.POSIXct(Sys.time())),file="temp/logstage1.txt",append=T)

# RUN THE LOOP
# (APPROX 9 MIN)
estlist <- foreach(i=seq(length(listlad)), .packages=pack) %dopar% {
  
  # STORE ITERATION
  cat("\n", "iter=",i, as.character(Sys.time()), "\n",
    file="temp/logstage1.txt", append=T)
  
  # EXTRACT THE SELECTED LSOAs BELONGING TO THIS LAD
  lsoasel <- lsoadata$LSOA11CD[lsoadata$LAD11CD==listlad[i]]
  
  # EXTRACT THE LSOA DATA AND APPEND THEM
  datalist <- lapply(seq(length(lsoasel)), function(j) {
    temp <- data.table(readRDS(paste0(tmeanpath, lsoasel[j], ".Rds")))
    cbind(data.table(LSOA11CD=lsoasel[j]), date=date, temp)
  })
  data <- do.call(rbind, datalist)
  
  # SELECT FROM YEAR 2000 ON
  data <- data[data$date>=datestart,]
  
  # SUMMARISE TEMPERATURE DATA BY LSOA (USED LATER AS META-PREDICTORS)
  tmeanlsoa <-  data[, list(meantmean=mean(tmean), 
    rangetmean=diff(range(tmean))), by=LSOA11CD]
  
  # GENERATE VARIABLES FOR THE ANALYSIS
  data$time <- as.numeric(data$date)
  data$year <- year(data$date)
  data$month <- month(data$date)
  data$doy <- yday(data$date)
  data$dow <- factor(wday(data$date))
  
  # COMPUTE TEMPERATURE PERCENTILES
  tmeanper <- quantile(data$tmean, predper/100, na.rm=T)
  
  # CREATE STRATUM VARIABLE
  data$stratum <- with(data, factor(paste(LSOA11CD,year,month,sep=":")))
  
  # PARAMETERIZE THE CB of temperature
  argvar <- list(fun=varfun, knots=quantile(data$tmean, varper/100, na.rm=T))
  argvar$degree <- vardegree
  arglag <- list(fun=lagfun, knots=lagknots)
  
  # CREATE THE CB of temperature
  cbtemp <- crossbasis(data$tmean, lag=maxlag, argvar=argvar, arglag=arglag,
    group=factor(data$LSOA11CD))
  
  # SPECIFY KNOTS OF SPLINES OF TIME 
  tknots <- equalknots(data$time, df=nkseas*length(unique(data$year)))
  
  # RUN MODEL THE MODEL ON NON-EMPTY STRATA
  sub <- tapply(data$count, data$stratum, sum)[data$stratum]>0
  mod <- gnm(count ~ cbtemp + ns(time,knots=tknots) + dow, eliminate=stratum,
    family=quasipoisson(), data=data, na.action="na.exclude", subset=sub)
  
  # PREDICT
  redall <- crossreduce(cbtemp, mod, cen=15)
  
  # RETURN
  list(coefall=coef(redall), vcovall=vcov(redall),
    conv=mod$converged,
    disp=sum(residuals(mod,type="pearson")^2, na.rm=T)/mod$df.residual,
    tmeanper=tmeanper, tmeanlsoa=tmeanlsoa, ndeaths=sum(data$count,na.rm=T))
  
}
names(estlist) <- listlad

# REMOVE PARALLELIZATION
stopCluster(cl)
file.remove("temp/logstage1.txt")

################################################################################
# CHECKS AND SAVE

# CHECK CONVERGENCE
all(sapply(estlist, function(x) x$conv))

# SAVE
save.image("temp/firststage.RData")
