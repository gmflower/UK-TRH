################################################################################
# UK-TRH: SMALL-AREA ANALYSIS OF TEMPERATURE RELATED HOSPITALISATIONS IN ENGLAND
################################################################################

################################################################################
# FIRST STAGE
################################################################################

# Create amended agevarlab:
agevarlab <- c("age064", "age6574", "age7584", "age85plus", "total")

# PREPARE THE PARALLELIZATION
ncores <- detectCores()
cl <- parallel::makeCluster(max(1,ncores-2))
registerDoParallel(cl)

# PACKAGE LIST FOR PARALLELIZATION
pack <- c("dlnm", "data.table", "gnm", "tsModel", "splines")

# WRITE A TEXT FILE TO TRACE ITERATIONS
writeLines(c(""), "temp/logstage1.txt")
cat(as.character(as.POSIXct(Sys.time())),file="temp/logstage1.txt",append=T)

# RUN THE LOOP, FIRST BY LAD
# NB: split AUTOMATICALLY RE-ORDER BY THE SPLITTING VAR
stage1list <- foreach(hes=split(hesdata, hesdata$LAD11CD), 
  dtmean=split(datatmean,datatmean$LAD11CD), i=seq(listlad),
  .packages=pack) %dopar% {
    
    # STORE ITERATION (1 EVERY 10)
    if(i%%10==0) cat("\n", "iter=",i, as.character(Sys.time()), "\n",
      file="temp/logstage1.txt", append=T)
    
    # CREATE TIME VARS
    dtmean[, time:=as.numeric(date)]
    dtmean[, year:=year(date)]
    dtmean[, month:=month(date)]
    dtmean[, doy:=yday(date)]
    dtmean[, dow:=wday(date)]
    
    # MERGE HOLIDAYS (ALSO FILLING MISSING) AND ORDER
    dtmean <- merge(dtmean, holy, by="date", all.x=T)
    dtmean[, holy:=nafill(as.numeric(holy), fill=0)]
    setkey(dtmean, LSOA11CD, date)
    
    # COMPUTE TEMPERATURE PERCENTILES AT LAD LEVEL
    ladtmeanper <- quantile(dtmean$tmean, predper/100, na.rm=T)
    
    # COMPUTE TEMPERATURE PERCENTILES AT LSOA LEVEL (TO BE STORED)
    lsoatmeanper <-  dtmean[, lapply(c(0,varper,100), function(x) 
      quantile(tmean, x/100, na.rm=T)), by=LSOA11CD] |> as.data.frame()
    names(lsoatmeanper)[-1] <- paste0(c(0,varper,100), ".0%")
    
    # CREATE STRATUM VARIABLE
    dtmean[, stratum:=factor(paste(LSOA11CD,year,month,sep=":"))]

    # PARAMETERIZE THE CB OF TEMPERATURE
    #argvar <- list("thr", thr=16)
    argvar <- list(fun=varfun, knots=ladtmeanper[paste0(varper, ".0%")])
    #arglag <- list("strata", df=1)
    arglag <- list(fun=lagfun, knots=lagknots)
    
    # CREATE THE CB of temperature
    cbtemp <- crossbasis(dtmean$tmean, lag=maxlag, argvar=argvar, arglag=arglag,
      group=paste0(dtmean$LSOA11CD, dtmean$year)) 
    
    # KNOTS OF SPLINE OF DAY OF THE YEAR
    kseas <- equalknots(dtmean$doy, df=dfseas)

    # LOOP ACROSS CAUSES
    clist <- lapply(seq(setcause), function(k) {
      
      # SELECT HES
      hescause <- hes[cause==setcause[k],]
      
      # RESHAPE COUNTS BY AGE AS COLUMNS
      hescause <- dcast(hescause, cause+LSOA11CD+date~agegr, value.var="count",
        fill=0) 
      
      # MERGE HES AND TMEAN DATA (SINGLE LAD, SINGLE CAUSE)   
      # NB: KEEP THE CTS STRUCTURE BY KEEPING ALL TMEAN DATA
      # THEN FILL THE MISSING COUNTS
      data <- merge(hescause, dtmean, all.y=T, by.x=c("LSOA11CD", "date"),
        by.y=c("LSOA11CD", "date"))
      
      # Check a column for each age group exists:
      for (a in seq(agevarlab)) {
        if(agevarlab[a] %in% colnames(data)) {
          print(paste0("age group ",agevarlab[a]," already in dataset"))
        } else
          as.data.table(data[, (agevarlab[a]):=as.numeric(NA)])
      }
      
      data[, (agevarlab):=lapply(.SD, nafill, fill=0), .SDcols=agevarlab]
      #data[, (unique(hesdata$agegr)):=lapply(.SD, nafill, fill=0), .SDcols=unique(hesdata$agegr)]
      #data$total[is.na(data$total)] <- 0
      
      # LOOP ACROSS AGE GROUPS
      estlist <- lapply(seq(agevarlab), function(j) {
      #estlist <- lapply(seq(unique(hesdata$agegr)), function(j) {

        # Create counts for the specific age group:
        data$count <- data[[agevarlab[j]]]
        #data$count <- data[[unique(hesdata$agegr)[j]]]
        
        # Check sufficient counts:
        if (sum(data$count)>=15) {
          
          print(paste0("Section 1 ","age group ",j))
          # Run the model on non-empty stratum:
          data[, sub:=sum(count)>0, by=list(stratum)]

          # Seasonal analysis only:
          mod <- gnm(count ~ cbtemp + ns(doy,knots=kseas) + factor(dow) + holy,
          eliminate=stratum, family=quasipoisson(), data=data,
          na.action="na.exclude", subset=sub)        
        
          # check NA for cbtemp
          if (is.na(mod[["coefficients"]][["cbtempv1.l1"]])) {
            
            print(paste0("Section 2 ","age group ",j))
            list(coefall=coef(list(coefficients=list(cbtemp=NA))), vcovall=NA, conv=mod$converged,
            disp=sum(residuals(mod,type="pearson")^2, na.rm=T)/mod$df.residual,
            nevent=sum(data$count,na.rm=T))
            
          } else {
            
            print(paste0("Section 3 ","age group ",j))
            redall <- crossreduce(cbtemp, mod, cen=15)
            list(coefall=coef(redall), vcovall=vcov(redall), conv=mod$converged,
            disp=sum(residuals(mod,type="pearson")^2, na.rm=T)/mod$df.residual,
            nevent=sum(data$count,na.rm=T))
          }
        } else {
          
          #dont run the model and save empty list
          print(paste0("Section 4 ","age group ",j))
          list(coefall=coef(list(coefficients=list(cbtemp=NA))), vcovall=NA, conv=NA,
            disp=as.numeric(NA),
            nevent=sum(data$count,na.rm=T))
          
        }
        
      })
      
      # RENAME AND RETURN
      names(estlist) <- agevarlab
      estlist
    })

    # NAME
    names(clist) <- setcause
    
    # RETURN ESTIMATES ABOVE, LAD TMEAN DISTRIBUTUON, LSOA TMEAN AVERAGE AND RANGE,
    #  AND LSOA-SPECIFIC PERCENTILES
    list(clist=clist, ladtmeanper=ladtmeanper, lsoatmeanper=lsoatmeanper)
    
  }
names(stage1list) <- listlad

# REMOVE PARALLELIZATION
stopCluster(cl)

################################################################################
# CHECKS, CLEAN, AND SAVE

# CHECK CONVERGENCE AND DISPERSION
all(unlist(lapply(stage1list, function(y)
  lapply(y$clist, function(x) sapply(x, "[[", "conv")))))
plot(unlist(lapply(stage1list, function(y)
  lapply(y$clist, function(x) sapply(x, "[[", "disp")))))

# CLEAN
#file.remove("temp/logstage1.txt")

# SAVE FIRST-STAGE OBJECT
saveRDS(stage1list, "C:/Users/LSHGF3/Documents/RProjects/UK-TRH/temp/asthma/stage1list.RDS")
