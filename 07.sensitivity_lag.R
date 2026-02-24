################################################################################
# UK-TRH: SMALL-AREA ANALYSIS OF TEMPERATURE RELATED HOSPITALISATIONS IN ENGLAND
################################################################################

################################################################################
# SENSITIVITY ANALYSIS - EXTENDED LAG (5 DAYS) (REPRODUCIBLE FROM LINE 166)
################################################################################

# EXPOSURE-RESPONSE PARAMETERIZATION
varfun <- "ns"
varper <- c(50,90)

# LAG-REPONSE PARAMETERIZATION
maxlag <- 5
lagfun <- "ns"
lagknots <- logknots(maxlag, 1)

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
                        argvar <- list(fun=varfun, knots=ladtmeanper[paste0(varper, ".0%")])
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
                          
                          # ENSURE ALL AGE GROUP COLUMNS EXIST AND FILL NON-CASE DAYS
                          for (a in seq(agevarlab)) {
                            if(agevarlab[a] %in% colnames(data)) {
                              print(paste0("age group ",agevarlab[a]," already in dataset"))
                            } else
                              as.data.table(data[, (agevarlab[a]):=as.numeric(NA)])
                          }
                          data[, (agevarlab):=lapply(.SD, nafill, fill=0), .SDcols=agevarlab]
                          
                          # LOOP ACROSS AGE GROUPS
                          estlist <- lapply(seq(agevarlab), function(j) {
                            # CREATE COUNTS
                            data$count <- data[[agevarlab[j]]]
                            
                            # CHECK SUFFICIENT COUNTS
                            if (sum(data$count)>=15) {
                              
                              print(paste0("Section 1 ","age group ",j))
                              
                              # RUN THE MODEL ON NON-EMPTY STRATA
                              data[, sub:=sum(count)>0, by=list(stratum)]
                              mod <- gnm(count ~ cbtemp + ns(doy,knots=kseas) + factor(dow) + holy + no2mean + pm25mean,
                                         eliminate=stratum, family=quasipoisson(), data=data,
                                         na.action="na.exclude", subset=sub)        
                              
                              # RETURN MODEL COEFFICIENTS
                              if (is.na(mod[["coefficients"]][["cbtempv1.l1"]])) {
                                print(paste0("Section 2 ","age group ",j))
                                list(coefall=coef(list(coefficients=list(cbtemp=NA))), vcovall=NA, conv=mod$converged,
                                     disp=sum(residuals(mod,type="pearson")^2, na.rm=T)/mod$df.residual,
                                     nevent=sum(data$count,na.rm=T))
                              } else {
                                print(paste0("Section 3 ","age group ",j))
                                redall <- crossreduce(cbtemp, mod, cen=ladtmeanper[["50.0%"]])
                                list(coefall=coef(redall), vcovall=vcov(redall), conv=mod$converged,
                                     disp=sum(residuals(mod,type="pearson")^2, na.rm=T)/mod$df.residual,
                                     nevent=sum(data$count,na.rm=T))
                              }
                            } else {
                              # DON'T RUN THE MODEL FOR LOW COUNTS
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
                        names(clist) <- setcause
                        # RETURN ESTIMATES ABOVE, LAD TMEAN DISTRIBUTUON, LSOA TMEAN AVERAGE AND RANGE,
                        #  AND LSOA-SPECIFIC PERCENTILES
                        list(clist=clist, ladtmeanper=ladtmeanper, lsoatmeanper=lsoatmeanper)
                      }
names(stage1list) <- listlad

# REMOVE PARALLELIZATION
stopCluster(cl)

################################################################################
# CHECKS, CLEAN AND SAVE

# CHECK CONVERGENCE AND DISPERSION
all(unlist(lapply(stage1list, function(y)
  lapply(y$clist, function(x) sapply(x, "[[", "conv")))))
plot(unlist(lapply(stage1list, function(y)
  lapply(y$clist, function(x) sapply(x, "[[", "disp")))))

# CLEAN
file.remove("temp/logstage1.txt")

# SAVE
saveRDS(stage1list, "./data/stage1list_senslag.RDS")

################################################################################
# SECOND STAGE (REPRODUCIBLE)

# READ DATA

stage1list <- readRDS("./data/stage1list_senslag.RDS")

# EXTRACT THE RESULTS

# INITIALISE OBJECTS
coefs_list <- vcovs_list <- pooled_results_list <- 
  coefs_list_c <- vcovs_list_c <- pooled_results_list_c <- NULL

for (c in seq(setcause)) {
  print(setcause[c])
  for (a in seq(agevarlab)) {
    print(agevarlab[a])
    
    # COLLATE COEFFICIENTS, SKIPPING NAS
    coefs<-lapply(seq(stage1list), function(i) unlist(
      stage1list[[i]][["clist"]][[setcause[c]]][[agevarlab[a]]][["coefall"]] ))
    lad_coef <- NULL
    coefs_all <- NULL
    for (i in 1:length(listlad)) {
      lad_coef <- unlist(coefs[i])
      if (is.na(lad_coef[1])) {
      } else {
        coefs_all <- rbind(coefs_all, lad_coef)
      }
    }
    coefs_list[a]<-list(coefs_all)
    
    # COLLACE VARIANCE INFORMATION
    vcovs<-lapply(seq(stage1list), function(i) stage1list[[i]][["clist"]][[setcause[c]]][[agevarlab[a]]][["vcovall"]])
    vcovs<-Filter(function(a) any(!is.na(a)), vcovs) 
    vcovs_list[a]<-list(vcovs)
    
  }
  
  names(coefs_list) <- agevarlab
  names(vcovs_list) <- agevarlab
  
  coefs_list_c[c] <- list(coefs_list)
  vcovs_list_c[c] <- list(vcovs_list)
  
}

names(coefs_list_c) <- setcause
names(vcovs_list_c) <- setcause

# POOL THE RESULTS IN THE SECOND STAGE

# RECREATE THE ERF
argvar <- list(
  fun = varfun,
  knots = quantile(datatmean$tmean, varper/100, na.rm = T),
  Bound = range(datatmean$tmean, na.rm = T)
)
bvar <- do.call(onebasis,c(list(x=datatmean$tmean), argvar))

# IMPLEMENT THE META ANALYSIS FOR EACH CAUSE AND AGE GROUP
for (c in seq(setcause)) {
  
  for (a in seq(agevarlab)) {
    mix<-mixmeta(coefs_list_c[[setcause[c]]][[agevarlab[a]]]~1, 
                 vcovs_list_c[[setcause[c]]][[agevarlab[a]]], method="reml", 
                 na.action = "na.omit")
    s <- summary(mix)
    
    # CALCULATE PREDICTIONS
    pooled_results_list[a] <- 
      list(crosspred(bvar, coef=coef(mix), vcov=vcov(mix), model.link="log",
                     by=0.1, cen=quantile(datatmean$tmean, 50/100, na.rm = T)))
  }
  
  names(pooled_results_list) <- agevarlab
  pooled_results_list_c[c] <- list(pooled_results_list)
}

names(pooled_results_list_c) <- setcause

saveRDS(pooled_results_list_c, "./temp/pooled_results_senslag.RDS")

################################################################################
# FOREST PLOT

# SET TEMPERATURE REF POINTS
heat_temp <- round(quantile(datatmean$tmean, 99/100, na.rm = T),1)

names <- unlist(lapply(seq(causelist), function(x) causelist[[x]][["name"]]))
shortnames <- unlist(lapply(seq(causelist), function(x) causelist[[x]][["shortname"]]))
cause_ref_list <- data.frame(names, shortnames)

rr <- rr_lowci <- rr_upci <- cause <- all_plot_data <- NULL
plot_data <- NULL

for (c in seq(setcause)) {
  rr <- rr_lowci <- rr_upci <- cause <- NULL
  for (a in seq(ages)) {
    
    rr <- c(rr,pooled_results_list_c[[setcause[c]]][[ages[a]]][["allRRfit"]][[paste0(heat_temp)]])
    rr_lowci <- c(rr_lowci,pooled_results_list_c[[setcause[c]]][[ages[a]]][["allRRlow"]][[paste0(heat_temp)]])
    rr_upci <- c(rr_upci,pooled_results_list_c[[setcause[c]]][[ages[a]]][["allRRhigh"]][[paste0(heat_temp)]])
    cause <- setcause[c]
  }
  
  plot_data <- data.frame(rr, rr_lowci, rr_upci, ages, cause)
  all_plot_data <- rbind(all_plot_data,plot_data)
}

all_plot_data <- left_join(all_plot_data,cause_ref_list,  by = join_by(cause == shortnames))


all_plot_data$grouping <- NULL
all_plot_data$grouping[all_plot_data$cause %in% c("resp","ari","pneumonia","copd","asthma")] <- "Respiratory" 
all_plot_data$grouping[all_plot_data$cause %in% c("cvd","mi","pulmheart","stroke","ihd","hf", "hypo")] <- "Cardiovascular" 
all_plot_data$grouping[all_plot_data$cause %in% c("genito","renal","arf","ckd")] <- "Genitourinary"
all_plot_data$grouping[all_plot_data$cause %in% c("infect", "bacterial")] <- "Infectious and parasitic" 
all_plot_data$grouping[all_plot_data$cause %in% c("endo","diabetes", "metabolic")] <- "Endocrine, nutritional, metabolic"

# SET ORDER
all_plot_data <-all_plot_data[with(all_plot_data, 
                                   order(grouping, names, cause, ages)),]
all_plot_data <- all_plot_data[, c("grouping", "names", "cause", "ages", "rr", "rr_lowci", "rr_upci")]
all_plot_data$id <- as.numeric(factor(all_plot_data$cause, 
                                      levels = c("cvd","hf","hypo","mi","stroke",
                                                 "endo","diabetes","metabolic","genito","arf","renal","infect","bacterial",
                                                 "resp","ari","asthma","copd","pneumonia"))) +
  as.numeric(factor(all_plot_data$grouping, 
                    levels = unique(all_plot_data$grouping)))
all_plot_data$ages <- factor(all_plot_data$ages, levels = ages)
all_plot_data$age_lab <- "first"
all_plot_data$age_lab[all_plot_data$ages=="total"] <- "second"


# SET GRID
unid <- unique(all_plot_data$id)
bglines <- data.frame(pos = unid[-1][diff(unid)==1] - .5)

# AGE GROUP SHAPES
all_plot_data$shape_group <- ifelse(all_plot_data$ages == "total", 19, 5)

# ICD POSITIONING
icd_pos <- aggregate(id ~ grouping, data = all_plot_data, mean)

# PLOT
bgplot <- ggplot(all_plot_data,
                 aes(x=id, group = ages, col = ages)) +
  theme_classic() +
  scale_x_reverse() +
  scale_x_reverse(name = "",
                  breaks = unique(all_plot_data$id),
                  labels = unique(all_plot_data$names), 
                  sec.axis = sec_axis(trans = ~., name = "", breaks = icd_pos$id, 
                                      labels = icd_pos$grouping)) +
  geom_hline(yintercept = 1) +
  theme(axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.text.x.bottom = element_text(angle = 90, vjust = 0, hjust = 1,
                                          size = 8), plot.title = element_text(hjust = 0.5, size = 9),
        axis.text.x.top = element_text(size = 8),
        panel.grid.major.y = element_line(linetype = 3, colour = "grey")) +
  geom_pointrange(position = position_dodge(0.8), size = 0.3) +
  scale_y_continuous(limits = range(with(all_plot_data, c(rr_lowci, rr_upci)))) +
  ggtitle("Relative Risk of Hospital Admission\nEngland, 2008 - 2019") 

# CONFIDENCE INTERVALS
heatplot <- bgplot +
  aes(y = rr, ymin = rr_lowci, ymax = rr_upci, shape = factor(shape_group)) +  # Mapping shape to the new column
  scale_color_manual(guide = "none",
                     values = c("black",paletteer_c("ggthemes::Red", 4))) +
  scale_shape_manual(guide = "none",
                     values = c(19, 5)) +  # 1 for circles, 5 for diamonds
  scale_linetype_manual(values = c("solid","solid","solid","solid","dotted")) +
  ylab(sprintf("RR at 99th\ntemperature percentile")) 

# CREATE LEGEND
legplot <- ggplot(all_plot_data, aes(x = id, group = ages, col = ages, shape = factor(shape_group))) +
  theme_void() + ylim(c(0, 0)) +
  geom_pointrange(aes(y = rr, ymin = rr_lowci, ymax = rr_upci, shape = factor(shape_group))) +
  scale_color_manual(name = "Age group",
                     values = c("black",paletteer_c("ggthemes::Red", 4))) + 
  scale_shape_manual(guide = "none",
                     values = c(19,5)) +  # 1 for circles, 5 for diamonds
  theme(legend.position = "top", legend.direction = "horizontal") 

pdf("./figures/supplementary/combined_forest_plot_senslag.pdf", height=6.5, width=16)

# COMBINE PLOT AND LEGEND
heatplot / legplot + 
  plot_layout(heights = c(1, .1))

dev.off()

png_plot <- heatplot / legplot + 
  plot_layout(heights = c(1, .1))

# SAVE
ggsave("./figures/supplementary/combined_forest_plot_senslag.png", png_plot, height=7, width=14, dpi=1200)
write.csv(all_plot_data, file="./tables/forest_plot_data_senslag.csv")

