################################################################################
# UK-TRH: SMALL-AREA ANALYSIS OF TEMPERATURE RELATED HOSPITALISATIONS IN ENGLAND
################################################################################

################################################################################
# SECOND STAGE
################################################################################


# (1) Extract results from the first stage model:

coefs_list <- vcovs_list <- pooled_results_list <- coefs_list_c <- vcovs_list_c <- pooled_results_list_c <- NULL

for (c in seq(setcause)) {
    print(setcause[c])
  for (a in seq(agevarlab)) {
    print(agevarlab[a])
  
    # Collate coefficients, excluding NAs:
    coefs<-lapply(seq(stage1list), function(i) unlist(stage1list[[i]][["clist"]][[setcause[c]]][[agevarlab[a]]][["coefall"]] ))
    lad_coef <- NULL
    coefs_all <- NULL
    for (i in 1:length(listlad)) {
      lad_coef<-unlist(coefs[i])
      if (is.na(lad_coef[1])) {
        #nothing
      } else {
      coefs_all<-rbind(coefs_all,lad_coef)
      }
    }
    
    coefs_list[a]<-list(coefs_all)
    
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


# (2) Pool the results of the first stage model across LADs:

# Re-create the exposure response relationship:
argvar <- list(
  fun = "ns",
  knots = quantile(datatmean$tmean, c(50, 90) / 100, na.rm = T),
  Bound = range(datatmean$tmean, na.rm = T)
)
bvar <- do.call(onebasis,c(list(x=datatmean$tmean), argvar))

heterogeneity_results <- list()

# Loop over ages and causes to pool results, and store model predictions:
for (c in seq(setcause)) {
  
  for (a in seq(agevarlab)) {
  # Pool the results for a particular age group and cause/diagnosis:
  mix<-mixmeta(coefs_list_c[[setcause[c]]][[agevarlab[a]]]~1, vcovs_list_c[[setcause[c]]][[agevarlab[a]]], method="reml", na.action = "na.omit")
  
  s <- summary(mix)
  
  heterogeneity_results[[setcause[c]]][[agevarlab[a]]] <- list(
    tau2_diag = round(diag(mix$Psi), 5),
    tau2_scalar = round(mean(diag(mix$Psi)), 5),
    I2_vector = round(s$i2, 1),
    I2_mean = round(mean(s$i2), 1),
    Q = round(s$qstat$Q[1], 2),
    Q_df = s$qstat$df[1],
    Q_p = ifelse(s$qstat$p[1] < 0.001, "<0.001", round(s$qstat$p[1], 3))
  )
  
  # Calculate predictions for the pooled model using meta-analysis coefficients:
  pooled_results_list[a] <- list(crosspred(bvar, coef=coef(mix), vcov=vcov(mix), model.link="log",
    by=0.1, cen=quantile(datatmean$tmean, 50/100, na.rm = T)))
  
  }

  names(pooled_results_list) <- agevarlab
  pooled_results_list_c[c] <- list(pooled_results_list)
}

names(pooled_results_list_c) <- setcause


saveRDS(pooled_results_list_c, "./temp/pooled_results_list_c_extended.RDS")


# All ages:
for (c in setcause) {
  
  pdf(paste0("./figures/supplementary/",c,"_exposure_response.pdf"), height=9, width=9.5)
  
  layout(matrix(1:6,3, byrow=T))
  par(mar=c(4,4,2,0.5), las=1, mgp=c(2.5,1,0))
  
  for (a in c(5,1,2,3,4)) {
    plot(pooled_results_list_c[[c]][[a]], type="l", col=a, ylab="RR", ylim=c(.6,1.6), lwd=2,
      xlab="Temperature (C)", main=paste(agevarlab[a]))
  }
  
  dev.off()
}

# Just total age group:
for (c in setcause) {
  
  pdf(paste("./figures/supplementary/total_age_",c,"_exposure_response.pdf"), height=9, width=9.5)
  
  layout(matrix(1:1,1, byrow=T))
  par(mar=c(4,4,2,0.5), las=1, mgp=c(2.5,1,0))

  plot(pooled_results_list_c[[c]][[5]], type="l", col=5, ylab="RR", ylim=c(.6,1.6), lwd=2,
    xlab="Temperature (C)")
  
  dev.off()
}


# Overlayed:
layout(matrix(1, 1, byrow=T))
plot(pooled_results_list_c[[c]][[1]], type="l", col=1, ylab="RR", ylim=c(.9,1.2), lwd=2,
    xlab="Temperature (C)", main=paste(c))
legend(20, 1.15, legend=c(agevarlab),  
       lty=1, lwd=2, col=c(1,2,3,4,5)
)
for (a in 2:5) {
lines(pooled_results_list_c[[c]][[a]], type="l", col=a, ylab="RR", ylim=c(.9,1.2), lwd=2,
    xlab="Temperature (C)", main=paste(c))
}

