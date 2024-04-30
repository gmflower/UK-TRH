################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# EFFECTS
################################################################################

# PREPARE THE PARALLELIZATION
ncores <- detectCores()
cl <- parallel::makeCluster(max(1,ncores-2))
registerDoParallel(cl)

# PACKAGE LIST FOR PARALLELIZATION
pack <- c("dlnm", "data.table", "tsModel", "MASS", "mixmeta", "abind")

# WRITE A TEXT FILE TO TRACE ITERATIONS
writeLines(c(""), "temp/effects.txt")
cat(as.character(as.POSIXct(Sys.time())),file="temp/effects.txt",append=T)

################################################################################
# LOOPS BY CAUSE/LSOA/AGE

# NUMBER OF SIMULATION FOR eCI OF EXCESS EVENTS
nsim <- 500

# LOOP BY CAUSE
for(k in seq(setcause)) {
  
  # PRINT
  cat(setcause[k], "")
  
  # ETRACT COEF/VCOV OF THE META-REGRESSION AND SAMPLE
  coefmeta <- coefmetalist[[k]]
  vcovmeta <- vcovmetalist[[k]]
  set.seed(13041975)
  mvcoefsim <- mvrnorm(nsim, coefmeta, vcovmeta)
  
  # RUN THE LOOP BY LSOA (PARALLELISED)
  # NB: split AUTOMATICALLY RE-ORDER BY THE SPLITTING VAR
  effectlist <- foreach(dtmean=split(datatmean,datatmean$LSOA11CD),
    i=seq(listlsoa), .packages=pack) %dopar% {
    
    # STORE ITERATION (1 EVERY 100)
    if(i%%100==0) cat("\n", "iter=",i, as.character(Sys.time()), "\n",
      file="temp/effects.txt", append=T)
    
    # COMPUTE TEMPERATURE PERCENTILES
    tmeanper <- quantile(dtmean$tmean, predper/100, na.rm=T)
      
    # DEFINE PARAMETERS OF THE EXPOSURE-RESPONSE FUNCTION
    argvar <- list(fun=varfun, knots=tmeanper[paste0(varper, ".0%")],
      Bound=tmeanper[c("0.0%","100.0%")])
    argvar$degree <- vardegree
    
    # LOOP ACROSS AGE GROUPS
    estlist <- lapply(seq(agevarlab), function(j) {
      
      # RECONSTRUCT THE MODEL MATRIX OF THE META-REGRESSION AT LSOA LEVEL
      lsoavar <- cbind(agegr=agevarlab[j], lsoacomp[i,])
      Xdeslsoa <- paste0("~", paste0(c("agegr", names(lsoacomp)[-1]), collapse="+")) |>
        formula() |> delete.response() |> 
        model.matrix(data=lsoavar, xlev=list(agegr=agevarlab))
      
      # PREDICT COEF/VCOV FOR LSOA/AGE
      fit <- (Xdeslsoa %x% diag(vardf)) %*% coefmeta |> drop()
      vcov <- (Xdeslsoa %x% diag(vardf)) %*% vcovmeta %*% 
        t(Xdeslsoa %x% diag(vardf))
      lsoapar <- list(fit=fit, vcov=vcov)
      
      # IDENTIFY MMT/MMP (BETWEEN 1ST AND 99TH)
      indextr <- which(predper<1 | predper>99)
      bvar <- do.call(onebasis, c(list(x=tmeanper[-indextr]), argvar))
      mmt <- tmeanper[-indextr][[which.min(bvar%*%lsoapar$fit)]]
      mmp <- predper[[which(tmeanper%in%mmt)]]
      mmtmmp <- data.frame(LSOA11CD=listlsoa[i], agegr=agevarlab[j], 
        mmt=mmt,mmp=mmp)
      
      # COMPUTE RR
      cp <- crosspred(bvar, coef=lsoapar$fit, vcov=lsoapar$vcov, model.link="log",
        cen=mmt, at=tmeanper[c("1.0%","99.0%")])
      rr <- t(rbind(cp$allRRfit, cp$allRRlow, cp$allRRhigh))
      dimnames(rr) <- list(NULL, c("rr", "95%eCIlow" ,"95%eCIhigh"))
      rr <- cbind(data.frame(LSOA11CD=listlsoa[i], agegr=agevarlab[j], 
        effect=c("cold","heat")), rr)
      
      # DERIVE THE CENTERED BASIS
      bvar <- do.call(onebasis, c(list(x=dtmean$tmean), argvar))
      cenvec <- do.call(onebasis, c(list(x=mmt), argvar))
      bvarcen <- scale(bvar, center=cenvec, scale=F) 
      
      # FORWARD MOVING AVERAGE OF EVENTS, RECONSTRUCTED AS BASELINE COUNTS
      prop <- subset(lsoaprop, LSOA11CD==listlsoa[i] & agegr==agevarlab[j])$prop
      count <- subset(ladevent, LAD11CD==subset(lookup, LSOA11CD==listlsoa[i])$LAD11CD & 
          agegr==agevarlab[j] & cause==setcause[k])$count
      seas <- subset(seasevent, agegr==agevarlab[j] & cause==setcause[k])$mult
      events <- (count*prop*seas[yday(dtmean$date)]) |> Lag(-maxlag:0) |>
        rowMeans()
      
      # INDICATOR FOR COLD/HEAT DAYS
      indheat <- dtmean$tmean>=mmt
      
      # COMPUTE THE YEARLY EXCESS EVENTS
      anday <- (1-exp(-bvarcen%*%lsoapar$fit))*events
      an <- c(sum(anday[!indheat], na.rm=T), sum(anday[indheat], na.rm=T),
        sum(anday, na.rm=T)) / (length(anday)/365.25)
      
      # SIMULATED DISTRIBUTION OF YEARLY EXCESS EVENTS
      ansim <- sapply(seq(nsim), function(s) {
        coef <- drop((Xdeslsoa %x% diag(vardf)) %*% mvcoefsim[s,])
        anday <- (1-exp(-bvarcen%*%coef))*events
        c(sum(anday[!indheat], na.rm=T), sum(anday[indheat], na.rm=T),
          sum(anday, na.rm=T))
      }) / (length(anday)/365.25)
      ansim <- cbind(an, ansim)
      dimnames(ansim) <- list(c("cold","heat","tot"), 
        c("est", paste0("sim", seq(nsim))))
      
      # TOTAL YEARLY EVENTS
      totevent <- data.frame(LSOA11CD=listlsoa[i], agegr=agevarlab[j],
        totevent=sum(events, na.rm=T) / (length(anday)/365.25))
      
      # RETURN
      list(mmtmmp=mmtmmp, rr=rr, ansim=ansim, totevent=totevent)
    })
    
    # PUT TOGETHER BY AGE
    mmtmmp <- do.call(rbind, lapply(estlist, "[[", "mmtmmp"))
    rr <- do.call(rbind, lapply(estlist, "[[", "rr"))
    ansim <- abind(lapply(estlist, "[[", "ansim"), along=3)
    totevent <- do.call(rbind, lapply(estlist, "[[", "totevent"))
    
    # COMPUTE ALL-AGE OF EXCESS EVENTS AND ADD IT
    ansim <- abind(ansim, apply(ansim, 1:2, sum))
    dimnames(ansim)[[3]] <- c(agevarlab, "all")
    totevent <- rbind(totevent, data.frame(LSOA11CD=listlsoa[i], agegr="all",
      totevent=sum(totevent$totevent)))
    
    # RETURN
    list(mmtmmp=mmtmmp, rr=rr, ansim=ansim, totevent=totevent)
  }
  
  # CLEAN (ALSO FULL DATA TO FREE MEMORY)
  rm(coefmeta, vcovmeta, mvcoefsim)
  #file.remove("temp/effects.txt")
  
  ################################################################################
  # TABLES BY LSOA
  
  # EXTRACT TABLES OTHER THAN ansim
  mmtmmplsoa <- do.call(rbind, lapply(effectlist, "[[", "mmtmmp"))
  rrlsoa <- do.call(rbind, lapply(effectlist, "[[", "rr"))
  toteventlsoa <- do.call(rbind, lapply(effectlist, "[[", "totevent"))
  
  # EXTRACT ansim, BIND, AND PERMUTE ITS DIMENSIONS
  ansim <- aperm(abind(lapply(effectlist, "[[", "ansim"), along=4), c(4,3,1,2))
  dimnames(ansim)[[1]] <- listlsoa
  
  # EXTRACT AN FROM ANSIM AND TRANSFORM BOTH IN DATA.TABLES
  ansim <- as.data.table(ansim)
  names(ansim)[1:4] <- c(names(rrlsoa)[1:3], "sim")
  anlsoa <- ansim[sim=="est",c(1:3,5)]
  names(anlsoa)[4] <- "excevent"
  ansim <- ansim[sim!="est",]
  
  # MERGE TOTAL EVENTS AND POPULATION
  anlsoa <- merge(anlsoa, toteventlsoa, by=c("LSOA11CD", "agegr"))
  anlsoa <- merge(anlsoa, pop, by=c("LSOA11CD", "agegr"))
  setcolorder(anlsoa, c(1:3,5:6,4))
  ansim <- merge(ansim, pop, by=c("LSOA11CD", "agegr"))
  
  # COMPUTE AN AND eCI
  anlsoa <- merge(anlsoa, ansim[, list("95%eCIlow"=quantile(value, 0.025),
    "95%eCIhigh"=quantile(value, 0.975)), by=LSOA11CD:effect],
    by=c("LSOA11CD", "agegr", "effect"))
  
  # DERIVE LAD-SPECIFIC MULTIPLIERS TO COMPUTE LSOA-SPECIFIC STANDARDISED RATES
  # NB: AVOID ISSUES WITH SMALL NUMBER OF EVENTS/POP
  stdratelsoa <- merge(anlsoa, lookup[c("LSOA11CD","LAD11CD")], by="LSOA11CD") 
  stdratelsoa[, rate:=sum(totevent)/sum(pop)*100000, by=c("LAD11CD", "agegr")]
  stdratelsoa[, mult:=rate/pmax(totevent,1)] 
  
  # ADD mult TO ansim (ONLY FOR LSOA TABLES)
  temp <- merge(ansim, stdratelsoa[,c (1:3,11)],
    by=c("LSOA11CD", "agegr", "effect"))
  
  # COMPUTE STD RATE
  stdratelsoa <- stdratelsoa[agegr!="all", list(stdrate=sum(excevent*mult*stdweight)),
    by=c("LSOA11CD", "effect")]
  temp <- temp[agegr!="all", list(value=sum(value*mult*stdweight)),
    by=c("LSOA11CD", "effect", "sim")]
  stdratelsoa <- merge(stdratelsoa, temp[, list("95%eCIlow"=quantile(value, 0.025),
    "95%eCIhigh"=quantile(value, 0.975)), by=LSOA11CD:effect],
    by=c("LSOA11CD", "effect"))
  
  ################################################################################
  # TABLE BY LAD
  
  # AGGREGATE BY LAD
  anlad <- merge(anlsoa[,c(1:6)], lookup[c("LSOA11CD","LAD11CD")], by="LSOA11CD")
  anlad <- anlad[, list(totevent=sum(totevent), pop=sum(pop), excevent=sum(excevent)),
    by=c("LAD11CD", "agegr", "effect")]
  temp <- merge(ansim, lookup[c("LSOA11CD","LAD11CD")], by="LSOA11CD")
  temp <- temp[, list(value=sum(value), pop=sum(pop)),
    by=c("LAD11CD", "agegr", "effect", "sim")]
  
  # COMPUTE AN AND eCI
  anlad <- merge(anlad, temp[, list("95%eCIlow"=quantile(value, 0.025),
    "95%eCIhigh"=quantile(value, 0.975)), by=LAD11CD:effect],
    by=c("LAD11CD", "agegr", "effect"))
  
  # COMPUTE STD RATE
  stdratelad <- anlad[agegr!="all", list(stdrate=sum(excevent/pop*100000*stdweight)),
    by=c("LAD11CD", "effect")]
  temp <- temp[agegr!="all", list(value=sum(value/pop*100000*stdweight)),
    by=c("LAD11CD", "effect", "sim")]
  stdratelad <- merge(stdratelad, temp[, list("95%eCIlow"=quantile(value, 0.025),
    "95%eCIhigh"=quantile(value, 0.975)), by=LAD11CD:effect],
    by=c("LAD11CD", "effect"))
  
  ################################################################################
  # TABLE BY REGION
  
  # AGGREGATE BY REG
  anreg <- merge(anlsoa[,c(1:6)], lookup[c("LSOA11CD","RGN11CD")], by="LSOA11CD")
  anreg <- anreg[, list(totevent=sum(totevent), pop=sum(pop), excevent=sum(excevent)),
    by=c("RGN11CD", "agegr", "effect")]
  temp <- merge(ansim, lookup[c("LSOA11CD","RGN11CD")], by="LSOA11CD")
  temp <- temp[, list(value=sum(value), pop=sum(pop)),
    by=c("RGN11CD", "agegr", "effect", "sim")]
  
  # COMPUTE AN AND eCI
  anreg <- merge(anreg, temp[, list("95%eCIlow"=quantile(value, 0.025),
    "95%eCIhigh"=quantile(value, 0.975)), by=RGN11CD:effect],
    by=c("RGN11CD", "agegr", "effect"))
  
  # COMPUTE STD RATE
  stdratereg <- anreg[agegr!="all", list(stdrate=sum(excevent/pop*100000*stdweight)),
    by=c("RGN11CD", "effect")]
  temp <- temp[agegr!="all", list(value=sum(value/pop*100000*stdweight)),
    by=c("RGN11CD", "effect", "sim")]
  stdratereg <- merge(stdratereg, temp[, list("95%eCIlow"=quantile(value, 0.025),
    "95%eCIhigh"=quantile(value, 0.975)), by=RGN11CD:effect],
    by=c("RGN11CD", "effect"))
  
  ################################################################################
  # TABLE BY DEPRIVATION QUINTILES
  
  # FACTOR FOR IMD QUINTILES
  fimd <- lsoadata[c("LSOA11CD", "IMD")]
  fimd$fimd <- cut(fimd$IMD, 5, labels=paste0("IMDquint",1:5))
  
  # AGGREGATE BY IMD QUINTILES
  animd <- merge(anlsoa[,c(1:6)], fimd[c("LSOA11CD","fimd")], by="LSOA11CD")
  animd <- animd[, list(totevent=sum(totevent), pop=sum(pop), excevent=sum(excevent)),
    by=c("fimd", "agegr", "effect")]
  temp <- merge(ansim, fimd[c("LSOA11CD","fimd")], by="LSOA11CD")
  temp <- temp[, list(value=sum(value), pop=sum(pop)),
    by=c("fimd", "agegr", "effect", "sim")]
  
  # COMPUTE AN AND eCI
  animd <- merge(animd, temp[, list("95%eCIlow"=quantile(value, 0.025),
    "95%eCIhigh"=quantile(value, 0.975)), by=fimd:effect],
    by=c("fimd", "agegr", "effect"))
  
  # COMPUTE STD RATE
  stdrateimd <- animd[agegr!="all", list(stdrate=sum(excevent/pop*100000*stdweight)),
    by=c("fimd", "effect")]
  temp <- temp[agegr!="all", list(value=sum(value/pop*100000*stdweight)),
    by=c("fimd", "effect", "sim")]
  stdrateimd <- merge(stdrateimd, temp[, list("95%eCIlow"=quantile(value, 0.025),
    "95%eCIhigh"=quantile(value, 0.975)), by=fimd:effect],
    by=c("fimd", "effect"))
  
  # ERASE OBJ
  rm(fimd)
  
  ################################################################################
  # TABLE OVERALL
  
  # AGGREGATE ALL
  anall <- anlsoa[, list(totevent=sum(totevent), pop=sum(pop), excevent=sum(excevent)),
    by=c("agegr", "effect")]
  temp <- ansim[, list(value=sum(value), pop=sum(pop)),
    by=c("agegr", "effect", "sim")]
  
  # COMPUTE AN AND eCI
  anall <- merge(anall, temp[, list("95%eCIlow"=quantile(value, 0.025),
    "95%eCIhigh"=quantile(value, 0.975)), by=agegr:effect],
    by=c("agegr", "effect"))
  
  # COMPUTE STD RATE
  stdrateall <- anall[agegr!="all", list(stdrate=sum(excevent/pop*100000*stdweight)),
    by=c("effect")]
  temp <- temp[agegr!="all", list(value=sum(value/pop*100000*stdweight)),
    by=c("effect", "sim")]
  stdrateall <- merge(stdrateall, temp[, list("95%eCIlow"=quantile(value, 0.025),
    "95%eCIhigh"=quantile(value, 0.975)), by=effect], by=c("effect"))
  
  # SAVE THE OUTPUT
  save(anlsoa, stdratelsoa, anlad, stdratelad, anreg, stdratereg, animd,
    stdrateimd, anall, stdrateall, mmtmmplsoa, rrlsoa, toteventlsoa,
    file=paste0("output/output_",setcause[k],".RData"))
  
  # REMOVE OBJECTS
  rm(effectlist, ansim, temp)
  rm(mmtmmplsoa, rrlsoa, toteventlsoa)
  rm(anlsoa, stdratelsoa, anlad, stdratelad, anreg, stdratereg, animd,
    stdrateimd, anall, stdrateall)
}

# REMOVE PARALLELIZATION
stopCluster(cl)
