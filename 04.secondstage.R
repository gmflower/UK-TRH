################################################################################
# UK-TRH: SMALL-AREA ANALYSIS OF TEMPERATURE RELATED HOSPITALISATIONS IN ENGLAND
################################################################################

################################################################################
# SECOND STAGE (REPRODUCIBLE)
################################################################################

################################################################################
# READ DATA

stage1list <- readRDS("./data/stage1list.RDS")

################################################################################
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

################################################################################
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
