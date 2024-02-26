################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# FORECASTED EXCESS MORTALITY DURING THE HEATWAVE OF JULY 2022
################################################################################

# LOAD THE FORECASTED TEMPERATURE DURING THE JULY 2022 HEATWAVE
tmeanfcast <- read.csv("data/tmeanforecast.csv") |> arrange(LSOA11CD)
all(lookup$LSOA11CD %in% tmeanfcast$LSOA11CD)
fcastdates <- format(as.Date(names(tmeanfcast)[-1], format="%B%d_%Y"),
  format="%d %B %Y")
names(tmeanfcast)[-1] <- fcastdates

# PREPARE THE PARALLELIZATION
ncores <- detectCores()
cl <- parallel::makeCluster(max(1,ncores-2))
registerDoParallel(cl)

# WRITE A TEXT FILE TO TRACE ITERATIONS
writeLines(c(""), "temp/fcast.txt")
cat(as.character(as.POSIXct(Sys.time())),file="temp/fcast.txt",append=T)

# SAMPLE THE COEF OF THE META-REGRESSION
set.seed(13041975)
mvcoefsim <- mvrnorm(nsim, coef(metaall), vcov(metaall))

# RUN THE LOOP
fcastreslist <- foreach(i=seq(listlsoa), .packages=pack) %dopar% {
  
  # STORE ITERATION (1 EVERY 100)
  if(i%%100==0) cat("\n", "iter=",i, as.character(Sys.time()), "\n",
    file="temp/fcast.txt", append=T)
  
  # EXTRACT TEMPERATURE PERCENTILES (CONVERT TO VECTOR WITH NAMES)
  tmeanper <- lsoatmeanper[i,-1] |> as.matrix() |> drop()
  
  # DEFINE PARAMETERS OF THE EXPOSURE-RESPONSE FUNCTION
  argvar <- list(fun=varfun, knots=tmeanper[seq(varper)+1],
    Bound=tmeanper[-(seq(varper)+1)])
  argvar$degree <- vardegree
  
  # LOOP ACROSS AGE GROUPS
  estlist <- lapply(seq(agevarlab), function(j) {
    
    # DEFINE THE PREDICTORS AND THEN COEF/VCOV FOR LSOA/AGE
    lsoavar <-  cbind(agegr=agevarlab[j], 
      as.data.frame(predict(pca, newdata=lsoadata[i, metanames]))[,seq(ncomp)])
    lsoapar <- predict(metaall, newdata=lsoavar, vcov=T)
    
    # IDENTIFY THE MMT AND SET FORECASTED TEMPERATURES (ABOVE MMT ONLY) 
    mmt <- subset(mmtmmplsoa, LSOA11CD==lookup$LSOA11CD[i] & 
        agegr==agevarlab[j])$mmt
    tmean <- pmax(as.numeric(tmeanfcast[i,-1]), mmt)
    
    # DERIVE THE CENTERED BASIS (SUPPRESS WARNINGS DUE TO BOUNDARIES)
    bvar <- suppressWarnings(do.call(onebasis, c(list(x=tmean), argvar)))
    cenvec <- do.call(onebasis, c(list(x=mmt), argvar))
    bvarcen <- scale(bvar, center=cenvec, scale=F) 
    
    # DEATHS (FROM YEARLY AVERAGES)
    deaths <- subset(totdeathlsoa, LSOA11CD==lookup$LSOA11CD[i] & 
        agegr==agevarlab[j])$totdeath / 365.25

    # COMPUTE THE DAILY EXCESS DEATHS
    anday <- drop((1-exp(-bvarcen%*%lsoapar$fit))*deaths)
    
    # RECONSTRUCT THE MODEL MATRIX OF THE META-REGRESSION AT LSOA LEVEL
    Xdeslsoa <- model.matrix(delete.response(terms(metaall)), data=lsoavar,
      contrasts.arg=metaall$contrasts, xlev=metaall$xlevels)

    # SIMULATED DISTRIBUTION OF DAILY EXCESS DEATHS
    andaysim <- sapply(seq(nsim), function(s) {
      coef <- drop((Xdeslsoa %x% diag(vardf)) %*% mvcoefsim[s,])
      drop((1-exp(-bvarcen%*%coef))*deaths)
    })
    anday <- cbind(anday, andaysim)

    # RETURN
    return(anday)
  })
  
  # PUT TOGETHER BY AGE, THEN PERMUTE
  est <- abind(estlist, along=3) |> aperm(c(3,1,2))
  
  # ADD ALL-AGE
  all <- array(apply(est, 2:3, sum), dim=c(1,dim(est)[-1]))
  est <- abind(est,all, along=1)
  
  # ADD NAMES
  dimnames(est) <- list(c(agevarlab, "all"), fcastdates,
    c("est",paste0("sim", seq(nsim))))

  # RETURN
  return(est)
}

# REMOVE PARALLELIZATION
stopCluster(cl)

# CLEAN (ALSO FULL DATA TO FREE MEMORY)
#file.remove("temp/fcast.txt")

################################################################################
# EXTRACT LSOA-SPECIFIC, WITHOUT SIM, WITH POP AND TEMPERATURE DIFF WITH MAX

# EXTRACT POINT ESTIMATES ONLY
fcastlsoa <- lapply(fcastreslist, function(x) x[,,"est"]) |>
  abind(along=length(dim(fcastreslist[[1]]))) |> aperm(c(3,1:2))
dimnames(fcastlsoa)[[1]] <- listlsoa
fcastlsoa <- as.data.table(fcastlsoa)
names(fcastlsoa) <- c("LSOA11CD", "agegr", "date", "an")
fcastlsoa <- merge(fcastlsoa, pop, by=c("LSOA11CD","agegr"))

# ADD TEMPERATURE DIFFERENCE WITH MAX
lsoatmax <- data.table(lsoatmeanper[,c("LSOA11CD","100.0%")])
names(lsoatmax)[2] <- "tmeanmax"
lsoamaxfcast <- melt(data.table(tmeanfcast), 1, variable.factor=F)
names(lsoamaxfcast)[2:3] <- c("date","tmeanfcast")
lsoatmax <- merge(lsoatmax, lsoamaxfcast)
lsoatmax$diff <- lsoatmax$tmeanfcast - lsoatmax$tmeanmax
fcastlsoa <- merge(fcastlsoa, lsoatmax[,c("LSOA11CD","date","diff")],
  by=c("LSOA11CD","date"))
rm(lsoatmax,lsoamaxfcast)

################################################################################
# EXTRACT REGION-SPECIFIC, WITH SIM AND WITHOUT POP, ONLY HW PERIOD

# AGGREGATE BY REGION, KEEP ONLY HEATWAVE DAYS
fcastaggr <- lapply(unique(lookup$RGN11NM), function(reg) 
  Reduce('+', fcastreslist[lookup$RGN11NM==reg])) |> abind(along=4) |>
  aperm(c(4,1:3))
dimnames(fcastaggr)[[1]] <- unique(lookup$RGN11NM)
fcastaggr <- as.data.table(fcastaggr)
names(fcastaggr) <- c("Region", "agegr", "date", "type", "an")
fcastaggr <- fcastaggr[date %in% paste(17:19, "July 2022")]

# ADD HW PERIOD
temp <- fcastaggr[, list(an=sum(an)), by=c("Region","agegr","type")]
temp$date <- "17-19 July 2022"
fcastaggr <- rbind(fcastaggr, temp)
fcastaggr$date <- factor(fcastaggr$date, levels=unique(fcastaggr$date))

# REMOVE BIG OBJECT
#rm(fcastreslist, temp)

################################################################################
# TABLES

# BY REGION
temp <- fcastaggr[agegr=="all"&date=="17-19 July 2022", list(an=sum(an)),
  by=c("Region","type")]
fcastreg <- temp[type=="est"] |> 
  merge(temp[type!="est", list(anlow=quantile(an, 0.025)),  by=Region]) |> 
  merge(temp[type!="est", list(anhigh=quantile(an, 0.975)),  by=Region])
popreg <- merge(lookup[c("LSOA11CD","RGN11NM")], subset(pop, agegr=="all")) |>
  group_by(RGN11NM) |> summarise(pop=sum(pop)) |> rename(Region=RGN11NM)
fcastreg <- merge(fcastreg, popreg)
fcastreg[, paste0("rate",c("","low","high")):=(lapply(.SD, function(x) x/pop*10^6)),
  .SDcols=an:anhigh]
rm(temp, popreg)

# BY AGE
temp <- fcastaggr[agegr!="all"&date=="17-19 July 2022", list(an=sum(an)),
  by=c("agegr","type")]
fcastage <- temp[type=="est"] |> 
  merge(temp[type!="est", list(anlow=quantile(an, 0.025)),  by=agegr]) |> 
  merge(temp[type!="est", list(anhigh=quantile(an, 0.975)),  by=agegr])
popage <- merge(lookup[c("LSOA11CD","RGN11NM")], subset(pop, agegr!="all")) |>
  group_by(agegr) |> summarise(pop=sum(pop))
fcastage <- merge(fcastage, popage)
fcastage[, paste0("rate",c("","low","high")):=(lapply(.SD, function(x) x/pop*10^6)),
  .SDcols=an:anhigh]
rm(temp, popage)

# BY DATE
temp <- fcastaggr[agegr=="all"&date!="17-19 July 2022", list(an=sum(an)),
  by=c("date","type")]
fcastdate <- temp[type=="est"] |> 
  merge(temp[type!="est", list(anlow=quantile(an, 0.025)),  by=date]) |> 
  merge(temp[type!="est", list(anhigh=quantile(an, 0.975)),  by=date])
fcastdate$pop <- sum(subset(pop, agegr=="all")$pop)
fcastdate[, paste0("rate",c("","low","high")):=(lapply(.SD, function(x) x/pop*10^6)),
  .SDcols=an:anhigh]
rm(temp)

# TOT
temp <- fcastaggr[agegr=="all"&date=="17-19 July 2022", list(an=sum(an)),
  by=c("type")]
fcastall <- cbind(data.table(name="all"), temp[type=="est"],
  temp[type!="est", list(anlow=quantile(an, 0.025))],
  temp[type!="est", list(anhigh=quantile(an, 0.975))])
fcastall$pop <- sum(subset(pop, agegr=="all")$pop)
fcastall[, paste0("rate",c("","low","high")):=(lapply(.SD, function(x) x/pop*10^6)),
  .SDcols=an:anhigh]
rm(temp)

# PRODUCE TABLE
fcasttab <- rename(fcastreg, name=Region) |>
  rbind(rename(fcastage, name=agegr)) |> rbind(rename(fcastdate, name=date)) |>
  rbind(fcastall)
fcasttab <- cbind(fcasttab[,c("name","pop")],
  fcasttab[, lapply(.SD, function(x) formatC(x, format="f", digits=0, big.mark=",")),
    .SDcols=an:anhigh],
  fcasttab[, lapply(.SD, function(x) formatC(x, format="f", digits=1, big.mark=",")),
    .SDcols=rate:ratehigh]
)
fcasttab$pop <- formatC(fcasttab$pop, format="f", digits=0, big.mark=",")
fcasttab$an <- paste0(fcasttab$an, " (", fcasttab$anlow, " to ",
  fcasttab$anhigh, ")")
fcasttab$rate <- paste0(fcasttab$rate, " (", fcasttab$ratelow, " to ",
  fcasttab$ratehigh, ")")
fcasttab <- fcasttab[, c("name","pop","an","rate")]

# SAVE
write.csv(fcasttab, row.names=F, file="tables/fcasttab.csv")

################################################################################
# MAPS OF EXCESS DEATH RATE AND TEMPERATURE DIFFERENCE WITH MAX

fcastmap <- merge(lsoashp1, lookup[c("LSOA11CD","LAD11CD")]) |>
  merge(fcastlsoa[agegr=="all",c("LSOA11CD","date","an","pop","diff")]) |>
  mutate(rate=an/pop*10^6)
ladmap <- merge(ladshp, unique(lookup[c("LAD11CD","RGN11CD")]), sort=F)

colval <- c("white", colorRampPalette(c("yellow","orange","red","red4"))(11))
cutrate <-  c(-1000,0:10,1000)
ggplot(data=subset(fcastmap, date=="19 July 2022")) + 
  geom_sf(aes(fill=cut(rate,  cutrate, inc=T)), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.2) +
  scale_fill_manual(values=colval, drop=F,
    name="Heat-related mortality rate (x 1,000,000)") + 
  guides(fill=guide_colorsteps(barwidth=18, barheight=0.6, title.position="top",
    title.hjust=0.5)) +
  theme_void() +
  #facet_wrap(vars(date), labeller=function(x) format(x, format="%d %B %Y")) +
  theme(legend.position="bottom") 

ukfcastdeath <- ggplot(data=fcastmap) + 
  geom_sf(aes(fill=cut(rate,  cutrate, inc=T)), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.2) +
  scale_fill_manual(values=colval, drop=F,
    name="Heat-related mortality rate (x 1,000,000)") + 
  guides(fill=guide_colorsteps(barwidth=18, barheight=0.6, title.position="top",
    title.hjust=0.5)) +
  theme_void() +
  facet_wrap(vars(date)) +
  theme(legend.position="bottom") 

png("figures/fcastdeath.png", height=3400*0.7, width=4000*0.7, res=288)
ukfcastdeath
dev.off()

pdf("figures/fcastdeath.pdf", height=11.2*0.7, width=13.20*0.7)
ukfcastdeath
dev.off()

################################################################################
# MAP OF TEMPERATURE ABOVE MAX

ggplot(data=subset(fcastmap, date=="19 July 2022")) + 
  geom_sf(aes(fill=diff), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.2) +
  scale_fill_steps2(low="green", high="deeppink4", n.breaks=15,
    name="Temperature difference from recorded maximum (Celsius)") + 
  guides(fill=guide_colorsteps(barwidth=18, barheight=0.6, title.position="top",
    title.hjust=0.5)) +
  theme_void() +
  #facet_wrap(vars(date), labeller=function(x) format(x, format="%d %B %Y")) +
  theme(legend.position="bottom")

ukfcastdiff <- ggplot(data=fcastmap) + 
  geom_sf(aes(fill=diff), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.2) +
  scale_fill_steps2(low="green", high="deeppink4", n.breaks=15,
    name="Temperature difference from maximum recorded in 2000-2019") + 
  guides(fill=guide_colorsteps(barwidth=18, barheight=0.6, title.position="top",
    title.hjust=0.5)) +
  theme_void() +
  facet_wrap(vars(date)) +
  theme(legend.position="bottom") 

png("figures/fcastdiff.png", height=3400*0.7, width=4000*0.7, res=288)
ukfcastdiff
dev.off()