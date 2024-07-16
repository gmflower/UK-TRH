################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# PLOTS
################################################################################

################################################################################
# EXPOSURE-RESPONSES

# EXTRACT AVERAGE DISTRIBUTION OF TEMPERATURES
tmeanper <- sapply(stage1list, "[[", "ladtmeanper") |> rowMeans()
per <- names(tmeanper)

# DEFINE EXPOSURE-RESPONSE
argvar <- list(fun=varfun, knots=tmeanper[paste0(varper, ".0%")],
  Bound=tmeanper[c("0.0%","100.0%")])
argvar$degree <- vardegree
bvar <- do.call(onebasis, c(list(x=tmeanper), argvar))

# DEFINE META-PREDICTORS FOR AVERAGE CURVES BY AGE
# NB: ALL THE PCA COMPONENTS ARE 0 BY DEF
preddataage <- ladcomp[seq(agevarlab),-1]
preddataage[,] <- 0
preddataage <- cbind(agegr=agevarlab, preddataage)

# LABELS FOR PERCENTILES
labperc <- paste(c(0,1,25,50,75,99,100),sep="")


# LOOP ACROSS CAUSES
for(k in seq(setcause)) {

  # MULTI-PANEL PLOT
  pdf(paste0("figures/",setcause[k],"expresp.pdf"), height=6, width=8.5)
  layout(matrix(1:4, 2, byrow=T))
  par(mar=c(4,4,2,0.5), las=1, mgp=c(2.5,1,0))

    
  # LOOP
  for(j in seq(agelab)) {
    
    # EXTRACT COEF/VCOV OF META-ANALYTICAL MODEL
    coefmeta <- coefmetalist[[k]]
    vcovmeta <- vcovmetalist[[k]]
    
    # RECONSTRUCT THE MODEL MATRIX OF THE META-REGRESSION AT LSOA LEVEL
    Xdeslsoa <- paste0("~", paste0(c("agegr", names(lsoacomp)[-1]), collapse="+")) |>
      formula() |> delete.response() |> 
      model.matrix(data=preddataage[j,], xlev=list(agegr=agevarlab))
    
    # PREDICT PARAMETERS FOR AVERAGE CURVES BY AGE
    fit <- (Xdeslsoa %x% diag(vardf)) %*% coefmeta |> drop()
    vcov <- (Xdeslsoa %x% diag(vardf)) %*% vcovmeta %*% 
      t(Xdeslsoa %x% diag(vardf))
    
    #cen <- tmeanper[[which.min(bvar%*%fit)]]
    cen <- tmeanper[["90.0%"]]
    cpall <- crosspred(bvar, coef=fit, vcov=vcov, at=tmeanper, model.link="log",
      cen=cen)
    plot(cpall, type="n", ci="n", ylim=c(0.5,2.5), ylab="RR",
      xlab="Temperature percentile", xaxt="n", yaxt="n")
    axis(1, at=tmeanper[paste(labperc, ".0%", sep="")], 
      labels=paste(labperc, "%", sep=""))
    axis(2, at=-1:7/5+1)
    abline(v=tmeanper[paste(labperc, ".0%", sep="")], col=grey(0.85), lty=3)
    abline(h=-1:7/5+1, col=grey(0.85), lty=3)
    lines(cpall, ci="area", col=j, ci.arg=list(col=alpha(j,0.2)))
    title(setcause[k])
    mtext(paste("Aged", agelab[j]))
  }
dev.off()
}



################################################################################
# META-VARIABLES AND PCA

# CORRELATION MATRIX PLOT
corvarlad <- cor(laddata[metanames])
corvarlsoa <- cor(lsoadata[metanames])
dimnames(corvarlad) <- dimnames(corvarlsoa) <- list(metalab, metalab)
fcol <- colorRampPalette(c(muted("blue"), "white", muted("red")))
png("figures/corvar.png", height=800, width=1600, pointsize=24)
layout(t(1:2))
corrplot(corvarlad, tl.cex=0.5, mar=c(0,1,2,1), title="LAD level",
  col=fcol(200), tl.col="black")
corrplot(corvarlsoa, tl.cex=0.5, mar=c(0,1,2,1), title="LSOA level",
  col=fcol(200), tl.col="black")
dev.off()

# PCA VARIANCE PLOTS
png("figures/pcaplot.png", height=800, width=1400, pointsize=24)
layout(t(1:2))
par(mar=c(5,4,3,0.5), las=1, mgp=c(2.5,1,0))
screeplot(pca, type="l", npcs=10, pch=19, col=2, cex=1.5,
  main="Screeplot of the first 10 PCs")
title(xlab="PC#")
box()
cumpro <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
plot(cumpro[0:15], type="o", pch=19, col=2, cex=1.5, ylim=c(0,1), xlab="PC #", 
  ylab="Amount of explained variance", main="Cumulative variance plot")
abline(h=cumpro[ncomp], lty=5)
dev.off()

# EXTRACT COORDINATES (CORRELATIONS)
pcacoord <- t(apply(pca$rotation,1, "*", pca$sdev))
rownames(pcacoord) <- metalab

# PLOTS COORDINATES
png("figures/pcascore1.png", height=2000, width=2000, res=288)
layout(1:6, heights=c(rep(1,5),0.7))
par(mar=c(0,4,4,1), las=1, mgp=c(2.5,1,0))
for(i in seq(1:5)) {
  x <- pcacoord[,i]
  col <- c(grey(0.9),"red")[seq(x)%in%order(abs(x),decreasing=T)[1:4]+1]
  bar <- barplot(x, names.arg=F, las=2, ylim=c(-1,1), col=col, 
    main=paste("Principal component", i))
  abline(h=0)
}
par(mar=c(0,4,0,1))
plot(c(bar), c(bar), type="n", ylim=c(0,1), axes=F, ylab="", xlab="", bty="n",
  xlim=par()$usr[1:2], xaxs="i")
text(c(bar), 1, srt=50, adj=c(1,1), xpd=TRUE, labels=names(x), cex=1)
layout(1)
dev.off()

png("figures/pcascore2.png", height=1600, width=2400, res=288)
as.data.frame(pcacoord[,1:6]) %>%
  cbind(lab=factor(metalab, levels=metalab)) %>%
  pivot_longer(1:6) %>%
  ggplot() +
  geom_tile(aes(x=lab, y=name, fill=value), col="white") +
  geom_rect(aes(xmin=0.5, xmax=length(metalab)+0.5, ymin=6.5-ncomp, ymax=6.5), 
    fill=NA, col=1, linewidth=1.3) + 
  scale_fill_gradient2(low=muted("blue"), mid="white", high=muted("red"),
    midpoint=0, name="Correlation") +
  guides(fill=guide_colorsteps(barwidth=10, barheight=0.6, title.position="top",
    title.hjust=0.5)) + 
  scale_y_discrete(limits = rev) + 
  xlab(NULL) + ylab(NULL) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top")
dev.off()

################################################################################
# STANDARDISED RATES BY REGION, IMD QUINTILES, RURAL/URBAN CLASSIFICATION

# EXTRACT STANDARDISED RATES FOR DIFFERENT CATEGORIES
eff <- as.data.frame(rbind(unlist(stdrateall[,-1]),
  pivot_wider(stdratereg, names_from=2, values_from=3:5)[-1],
  pivot_wider(stdrateimd, names_from=2, values_from=3:5)[-1],
  pivot_wider(stdraterururb, names_from=2, values_from=3:5)[-1]))

# DEFINE THE X-AXIS POINTS
xpoint <- c(1, 4:13, 16:20, 23:27)

# LAYOUT
pdf("figures/rateplot.pdf", height=5.5, width=7)
layout(1:3, heights=c(1,1,0.5))
par(mar=c(0,5,3,1), las=1, mgp=c(2.5,1,0))

# HEAT
plot(stdrate_heat~xpoint, eff, type="n", axes=F, xlab="", ylim=c(-0.5,4),
  ylab="Heat-related std mortality\nrate (x 100,000)")
abline(h=-1:5, col=grey(0.8), lty=3)
arrows(x0=xpoint, x1=xpoint, y0=eff[,"95%eCIlow_heat"],
  y1=eff[,"95%eCIhigh_heat"], length=0.05, angle=90, code=3)
points(stdrate_heat~xpoint, eff, pch=16, col="red3", cex=1.2)
axis(side=2, cex.axis=0.8)

# ADD TITLES
titles <- c("Overall","Region","IMD quintiles","Rural/Urban\nclassification")
mtext(titles, side=3, at=c(1,8.5,18,25), font=2, cex=0.7, line=1)

# COLD
par(mar=c(1,5,1,1))
plot(stdrate_cold~xpoint, eff, type="n", axes=F, xlab="", ylim=c(80,160),
  ylab="Cold-related std mortality\nrate (x 100,000)")
abline(h=8:16*10, col=grey(0.8), lty=3)
arrows(x0=xpoint, x1=xpoint, y0=eff[,"95%eCIlow_cold"],
  y1=eff[,"95%eCIhigh_cold"], length=0.05, angle=90, code=3)
points(stdrate_cold~xpoint, eff, pch=16, col="navy", cex=1.3)
axis(side=2, cex.axis=0.8)

# DEFINE LABELS
labels <- c(arrange(unique(lookup[c("RGN11CD","RGN11NM")]), RGN11CD)[[2]],
  "1 - most deprived","2","3","4","5 - least deprived",
  levels(anrururb$frururb))

# ADD LABELS
par(mar=c(0,4,0,1))
plot(xpoint[-1], xpoint[-1], type="n", ylim=c(0,1), axes=F, ylab="", xlab="", bty="n",
  xlim=par()$usr[1:2], xaxs="i")
text(xpoint[-1], 1, srt=50, adj=c(1,1), xpd=TRUE, labels=labels, cex=0.9)

layout(1)
dev.off()
