################################################################################
# UK-TRH: SMALL-AREA ANALYSIS OF TEMPERATURE RELATED HOSPITALISATIONS IN ENGLAND
################################################################################

################################################################################
# SECOND STAGE
################################################################################

################################################################################
# PCA AND META-ANALYSIS

# RUN THE PCA
pca <- prcomp(laddata[,metanames], center=T, scale=T)
summary(pca)

# EXPANSION INDEX OF LADLSOA BY AGE GROUP
explad <- rep(seq(listlad), each=length(agelab))
explsoa <- rep(seq(listlsoa), each=length(agelab))

# EXTRACT THE COEF/VCOV AT LAD LEVEL
# NB:BIND AGE-SPECIFIC FOR COEF, UNLIST 1 LEVEL FOR VCOV
#ladcoef <- lapply(stage1list, function(x) 
#  t(sapply(x$estlist, "[[", "coefall"))) |> Reduce(rbind, x=_) |> 
#  as.data.frame() |> cbind(LAD11CD=listlad[explad], agegr=agevarlab) |>
#  relocate(LAD11CD, agegr) |> remove_rownames()

ladcoef <- lapply(stage1list,function(y)
  lapply(y$clist, function(x) t(sapply(x$estlist, "[[", "coefall"))))

#ladvcov <- unlist(lapply(stage1list, function(x) 
#  lapply(x$estlist, "[[", "vcovall")), recursive=F) |> lapply(vechMat) |>
#  lapply(t) |> Reduce(rbind, x=_) |> as.data.frame() |>
#  cbind(LAD11CD=listlad[explad], agegr=agevarlab) |> relocate(LAD11CD, agegr)

ladvcov <- lapply(stage1list,function(y)
  lapply(y$clist, function(x) lapply(x$estlist, "[[", "vcovall")))


# FIND THE OPTIMAL NUMBER OF COMPONENTS FOR THE PCA
# PRE-SET NUMBER OF COMPONENT TO 3
ncomp <- 3

# SAVE THE COMPONENTS AT LAD LEVEL
ladcomp <- cbind(LAD11CD=listlad, as.data.frame(pca$x[, seq(ncomp)]))

# PREDICT THE COMPONENTS AT LSOA LEVEL
lsoacomp <- predict(pca, newdata=lsoadata[, metanames])[,seq(ncomp)] |>
  as.data.frame() |> cbind(LSOA11CD=listlsoa) |> relocate(LSOA11CD)

# META-ANALYTICAL MODEL FOR OVERALL CUMULATIVE (WITH FIXD-EFFECTS)

# loop across causes to pick up right coefs (start cause loop)
coefmeta <- NULL
vcovmeta <- NULL
cmetaall <- NULL

for (i in listcause[c(1,2,5,6,11)]) {

print(i)
  
coefcause <- lapply(ladcoef, "[[", i) |> Reduce(rbind, x=_) |> 
  as.data.frame() |> cbind(LAD11CD=listlad[explad], agegr=agevarlab) |>
  relocate(LAD11CD, agegr) |> remove_rownames()

vcovcause <- unlist(lapply(ladvcov, "[[", i),recursive=F)  |> lapply(vechMat) |>
  lapply(t) |> Reduce(rbind, x=_) |> as.data.frame() |>
  cbind(LAD11CD=listlad[explad], agegr=agevarlab) |> relocate(LAD11CD, agegr)

fmeta <- paste0("cbind(", paste(names(coefcause[-c(1:2)]), collapse=", "),
  ") ~ ", paste(c("agegr", names(ladcomp)[-1]), collapse=" + ")) |>
  as.formula()

metaall <- mixmeta(fmeta, S=vcovcause[-c(1:2)], method="fixed", 
  data=merge(coefcause, ladcomp))
summary(metaall)

# EXTRACT META-ANALYTICAL COEF/VCOV
cmetaall[[i]] <- metaall
coefmeta[[i]] <- coef(metaall)
vcovmeta[[i]] <- vcov(metaall)
}

# end cause loop
