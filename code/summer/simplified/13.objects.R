################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# OBJECTS
################################################################################

# EXTRACT AND SAVE TEMPERATURE PERCENTILES AT LAD AND LSOA LEVEL
ladtmeanper <- lapply(stage1list, '[[', "ladtmeanper") |> lapply(t) |>
  Reduce(rbind, x=_) |> as.data.frame() |> cbind(LAD11CD=listlad) |> 
  relocate(LAD11CD)
lsoatmeanper <- lapply(stage1list, '[[', "lsoatmeanper") |> 
  Reduce(rbind, x=_) |> arrange(LSOA11CD)
saveRDS(ladtmeanper, file="objects/ladtmeanper.RDS")
saveRDS(lsoatmeanper, file="objects/lsoatmeanper.RDS")

# SAVE DATA AT LAD AND LSOA LEVEL
saveRDS(laddata, file="objects/laddata.RDS")
saveRDS(lsoadata, file="objects/lsoadata.RDS")

# SAVE COEF/VCOV AT LAD LEVEL
saveRDS(ladcoef, file="objects/ladcoef.RDS")
saveRDS(ladvcov, file="objects/ladvcov.RDS")

# SAVE LIST OF META-VARIABLES AND LABELS
saveRDS(metanames, file="objects/metanames.RDS")
saveRDS(metalab, file="objects/metalab.RDS")

# SAVE PCA COMPONENTS AT LAD AND LSOA LEVEL
saveRDS(ladcomp, file="objects/ladcomp.RDS")
saveRDS(lsoacomp, file="objects/lsoacomp.RDS")

# SAVE META-ANALYTICAL COEF/VCOV
saveRDS(coefmeta, file="objects/coefmeta.RDS")
saveRDS(vcovmeta, file="objects/vcovmeta.RDS")

# SAVE MMT/MMP AT LSOA LEVEL
saveRDS(mmtmmplsoa, file="objects/lsoammtmmp.RDS")
