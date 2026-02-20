################################################################################
# UK-TRH: SMALL-AREA ANALYSIS OF TEMPERATURE RELATED HOSPITALISATIONS IN ENGLAND
################################################################################

################################################################################
# SECOND STAGE
################################################################################

################################################################################
# RUN PCA, PREPARE META-VARIABLES, AND EXTRACT PARAMETERS

# RUN THE PCA
pca <- prcomp(laddata[,metanames], center=T, scale=T)
summary(pca)

# EXPANSION INDEX OF LADLSOA BY AGE GROUP
explad <- rep(seq(listlad), each=length(agelab))
explsoa <- rep(seq(listlsoa), each=length(agelab))

# EXTRACT THE COEF/VCOV AT LAD LEVEL
# NB: LIST OF DATAFRAMES BY CAUSE
ladcoeflist <- lapply(seq(setcause), function(k) {
  lapply(seq(stage1list), function(i) 
    t(sapply(stage1list[[i]]$clist[[k]], "[[", "coefall"))) |> 
    Reduce(rbind, x=_) |> as.data.frame() |> 
    cbind(LAD11CD=listlad[explad], agegr=agevarlab) |>
    relocate(LAD11CD, agegr) |> remove_rownames()
})
ladvcovlist <- lapply(seq(setcause), function(k) {
  lapply(seq(stage1list), function(i) 
    lapply(stage1list[[i]]$clist[[k]], "[[", "vcovall")) |> 
    unlist(recursive=F) |> lapply(vechMat) |> Reduce(rbind, x=_) |> 
    as.data.frame() |> cbind(LAD11CD=listlad[explad], agegr=agevarlab) |>
    relocate(LAD11CD, agegr) |> remove_rownames()
})
names(ladcoeflist) <- names(ladvcovlist) <- setcause

# FIND THE OPTIMAL NUMBER OF COMPONENTS FOR THE PCA
# PRE-SET NUMBER OF COMPONENT TO 3
ncomp <- 3

# SAVE THE COMPONENTS AT LAD LEVEL
ladcomp <- cbind(LAD11CD=listlad, as.data.frame(pca$x[, seq(ncomp)]))

# PREDICT THE COMPONENTS AT LSOA LEVEL
lsoacomp <- predict(pca, newdata=lsoadata[, metanames])[,seq(ncomp)] |>
  as.data.frame() |> cbind(LSOA11CD=listlsoa) |> relocate(LSOA11CD)

################################################################################
# META-ANALYTICAL MODEL FOR OVERALL CUMULATIVE (WITH FIXD-EFFECTS)

# EMPTY OBJECTS
coefmetalist <- vcovmetalist <- vector(mode="list", length=length(setcause))
names(coefmetalist) <- names(vcovmetalist) <- setcause
hetmeta <- matrix(NA, nrow=length(setcause), ncol=2, 
  dimnames=list(setcause, c("qstat","i2stat")))
#convmeta <- rep(NA, length(setcause))
#names(convmeta) <- setcause

# DEFINE FORMULA (NB: IDENTICAL FOR ALL CAUSES)
fmeta <- paste0("cbind(", paste(names(ladcoeflist[[1]][-c(1:2)]), collapse=", "),
  ") ~ ", paste(c("agegr", names(ladcomp)[-1]), collapse=" + ")) |>
  as.formula()

# RUN THE META-ANALYTICAL MODEL LOOPING ACROSS CAUSES
for(k in seq(setcause)) {
  
  # PRINT
  cat(setcause[k], "")
  
  # RUN THE META-REGRESSION
  metaall <- mixmeta(fmeta, S=ladvcovlist[[k]][-c(1:2)], method="fixed", 
    data=merge(ladcoeflist[[k]], ladcomp))
  
  # SUMMARY
  summeta <- summary(metaall)
  
  # STORE THE RESULTS
  coefmetalist[[k]] <- coef(metaall)
  vcovmetalist[[k]] <- vcov(metaall)
  #convmeta[i] <- metaall$converged
  hetmeta[k,] <- c(summeta$qstat$pvalue[1], summeta$i2stat[1])
}
(rm(metaall, summeta))

# CHECK CONVERGENCE
#all(convmeta)
