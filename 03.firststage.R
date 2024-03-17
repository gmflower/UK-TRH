################################################################################
# UK-TRH: SMALL-AREA ANALYSIS OF TEMPERATURE RELATED HOSPITALISATIONS IN ENGLAND
################################################################################

################################################################################
# FIRST STAGE
################################################################################

# PREPARE THE PARALLELIZATION
ncores <- detectCores()
cl <- parallel::makeCluster(max(1,ncores-2))
registerDoParallel(cl)

# PACKAGE LIST FOR PARALLELIZATION
pack <- c("dlnm", "data.table", "gnm", "tsModel", "splines")

# WRITE A TEXT FILE TO TRACE ITERATIONS
writeLines(c(""), "temp/logstage1.txt")
cat(as.character(as.POSIXct(Sys.time())),file="temp/logstage1.txt",append=T)

# RUN THE LOOP
# NB: split AUTOMATICALLY RE-ORDER BY THE SPLITTING VAR
#stage1list <- foreach(data=split(datafull, datafull$LAD11CD), i=seq(listlad),
#  .packages=pack) %dopar% {
stage1list <- foreach(hes=split(hesdata, hesdata$LAD11CD), i=seq(listlad),
  .packages=pack) %dopar% {
    
  # STORE ITERATION (1 EVERY 10)
  if(i%%10==0) cat("\n", "iter=",i, as.character(Sys.time()), "\n",
    file="temp/logstage1.txt", append=T)
 
    clist <- lapply(c(listcause[1],listcause[2],listcause[5],listcause[6],listcause[11]), function(k) {
    
    # Merge datatmean and hesdata(single LAD)   
    data <- merge(hes[cause==k], datatmean[LAD11CD==listlad[i],], all.y=T, by.x=c("LAD11CD", "LSOA11CD", "date"), by.y=c("LAD11CD", "LSOA11CD", "date")) 
    
    # Fill in missings:
    data[, (agevarlab):=lapply(.SD, nafill, fill=0), .SDcols=agevarlab]
    data$cause <- k
      
    # CREATE TIME VARS
    data[, time:=as.numeric(date)]
    data[, year:=year(date)]
    data[, month:=month(date)]
    data[, doy:=yday(date)]
    data[, dow:=wday(date)]
    
    # COMPUTE TEMPERATURE PERCENTILES AT LAD LEVEL
    ladtmeanper <- quantile(data[LAD11CD==listlad[i], tmean], predper/100, na.rm=T)
    
    # COMPUTE TEMPERATURE PERCENTILES AT LSOA LEVEL (TO BE STORED)
    lsoatmeanper <-  data[, lapply(c(0,varper,100), function(x) 
      quantile(tmean, x/100, na.rm=T)), by=LSOA11CD] |> as.data.frame()
    names(lsoatmeanper)[-1] <- paste0(c(0,varper,100), ".0%")
      
    # CREATE STRATUM VARIABLE
    data$stratum <- with(data, factor(paste(LSOA11CD,year,month,sep=":")))
    
    # PARAMETERIZE THE CB of temperature
    argvar <- list(fun=varfun, knots=ladtmeanper[paste0(varper, ".0%")])
    argvar$degree <- vardegree
    arglag <- list(fun=lagfun, knots=lagknots)

    # CREATE THE CB of temperature
    cbtemp <- crossbasis(data$tmean, lag=maxlag, argvar=argvar, arglag=arglag,
      group=factor(data$LSOA11CD))
    
    # SPECIFY KNOTS OF SPLINES OF TIME 
    tknots <- equalknots(data$time, df=nkseas*length(unique(data$year)))
    
      # LOOP ACROSS AGE GROUPS
      estlist <- lapply(seq(agevarlab), function(j) {
        
        # RUN MODEL THE MODEL ON NON-EMPTY STRATA
        data$count <- data[[agevarlab[j]]]
        data[, sub:=sum(count)>0, by=list(stratum)]
        mod <- gnm(count ~ cbtemp + ns(time,knots=tknots) + factor(dow),
          eliminate=stratum, family=quasipoisson(), data=data,
          na.action="na.exclude", subset=sub)
    
        # PREDICT
        redall <- crossreduce(cbtemp, mod, cen=15)
        
        # RETURN COEF/VCOV, CONVERGENCE, DISPERSION, TOTAL DEATHS
        list(coefall=coef(redall), vcovall=vcov(redall), conv=mod$converged,
          disp=sum(residuals(mod,type="pearson")^2, na.rm=T)/mod$df.residual,
          ndeaths=sum(data$count,na.rm=T))
      })
  
  # RENAME
  names(estlist) <- agevarlab
  
  # Create a list for results by cause:
  list(estlist,ladtmeanper=ladtmeanper, lsoatmeanper=lsoatmeanper)
  
  # Close the cause loop:
  })
  
  names(clist) <- c(listcause[1],listcause[2],listcause[5],listcause[6],listcause[11])
  
  # RETURN ESTIMATES ABOVE, LAD TMEAN DISTRIBUTUON, LSOA TMEAN AVERAGE AND RANGE,
  #  AND LSOA-SPECIFIC PERCENTILES
  #list(estlist=estlist, ladtmeanper=ladtmeanper, lsoatmeanper=lsoatmeanper)
  list(clist)
  
}
names(stage1list) <- listlad

# REMOVE PARALLELIZATION
stopCluster(cl)

################################################################################
# CHECKS, CLEAN, AND SAVE

# CHECK CONVERGENCE AND DISPERSION
all(unlist(lapply(stage1list, function(x) sapply(x$estlist, "[[", "conv"))))
plot(unlist(lapply(stage1list, function(x) sapply(x$estlist, "[[", "disp"))))

# CLEAN
#file.remove("temp/logstage1.txt")
