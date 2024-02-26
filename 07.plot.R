################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# PLOTS
################################################################################

################################################################################
# EXPOSURE-RESPONSES

# EXTRACT AVERAGE DISTRIBUTION OF TEMPERATURES
tmeanper <- colMeans(ladtmeanper[,-1])
per <- names(tmeanper)

# DEFINE EXPOSURE-RESPONSE
argvar <- list(fun=varfun, knots=tmeanper[paste0(varper, ".0%")],
  Bound=tmeanper[c("0.0%","100.0%")])
argvar$degree <- vardegree
bvar <- do.call(onebasis, c(list(x=tmeanper), argvar))

# PREDICT PARAMETERS FOR AVERAGE CURVES BY AGE
# NB: ALL THE PCA COMPONENTS ARE 0 BY DEF
preddataage <- ladcomp[seq(agevarlab),-1]
preddataage[,] <- 0
predparage <- predict(metaall, newdata=cbind(agegr=agevarlab, preddataage),
  vcov=T)

# LABELS FOR PERCENTILES
labperc <- paste(c(0,1,25,50,75,99,100),sep="")

# MULTI-PANEL PLOT
pdf("figures/expresp.pdf", height=6, width=8.5)
layout(matrix(1:4, 2, byrow=T))
par(mar=c(4,4,2,0.5), las=1, mgp=c(2.5,1,0))

# LOOP
for(i in seq(agelab)) {
  cen <- tmeanper[[which.min(bvar%*%predparage[[i]]$fit)]]
  cpall <- crosspred(bvar, coef=predparage[[i]]$fit, vcov=predparage[[i]]$vcov,
    at=tmeanper, model.link="log", cen=cen)
  plot(cpall, type="n", ci="n", ylim=c(0.8,2.5), ylab="RR",
    xlab="Temperature percentile", xaxt="n", yaxt="n")
  axis(1, at=tmeanper[paste(labperc, ".0%", sep="")], 
    labels=paste(labperc, "%", sep=""))
  axis(2, at=-1:7/5+1)
  abline(v=tmeanper[paste(labperc, ".0%", sep="")], col=grey(0.85), lty=3)
  abline(h=-1:7/5+1, col=grey(0.85), lty=3)
  lines(cpall, ci="area", col=i, ci.arg=list(col=alpha(i,0.2)))
  title(paste("Aged", agelab[i]))
}

dev.off()

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
# EFFECTS OF META-VARIABLES

# BACK-TRANSFORM COEF/VCOV IN THE ORIGINAL SCALE
# NB: USED FOR TESTING EACH ORIGINAL VARIABLE
L <- pca$rotation[,seq(ncomp)]
indcoef <- grep("PC", metaall$lab$p)
backcoef <- L %*% coef(metaall, format = "matrix")[indcoef,]
expL <- L %x% diag(vardf)
indvcov <- grep("PC", rownames(vcov(metaall)))
backvcov <- expL %*% vcov(metaall)[indvcov,indvcov] %*%  t(expL)

# DEFINE DISTRIBUTION OF AVERAGE VALUES OF CHARACTERISTICS (TWO ROWS)
avgmetavar <- as.data.frame(t(pca$center))[c(1,1),]

# MULTI-PANEL PLOT OF EXPOSURE-RESPONSES FOR DIFFERENT LEVELS
png("figures/effmod1.png", height=4500, width=3000, res=288)
layout(matrix(1:15, nrow=5, byrow=T))
par(mar=c(5,4,3,0.5), las=1, mgp=c(2.5,1,0))

# LOOP ACROSS CHARACTERISTICS
for(i in seq(ncol(avgmetavar))) {
  
  # SET CHARACTERISTIC
  var <- names(avgmetavar)[i]
  
  # CHANGE CHARACTERISTICS TO 10TH/90TH PERCENTILE VALUES
  metavarpred <- avgmetavar
  metavarpred[,var] <- quantile(lsoadata[,var], c(10,90)/100)
  
  # PREDICT THE VALUES OF COMPONENTS AND THEN COEF/VCOV
  pcavarpred <- as.data.frame(predict(pca, newdata=metavarpred)[,seq(ncomp)])
  ladvarpred <- cbind(agegr=agevarlab[3], pcavarpred)
  predpar <- predict(metaall, newdata=ladvarpred, vcov=T)
  
  # DEFINE CENTERING POINT
  cen <- tmeanper[[which.min(bvar%*%predpar[[1]]$fit)]]

  # CROSS-PREDICT
  cpall1 <- crosspred(bvar, coef=predpar[[1]]$fit, vcov=predpar[[1]]$vcov,
    at=tmeanper, model.link="log", cen=cen)
  cpall2 <- crosspred(bvar, coef=predpar[[2]]$fit, vcov=predpar[[2]]$vcov,
    at=tmeanper, model.link="log", cen=cen)
  
  # PLOT
  plot(cpall1, type="n", ci="n", ylim=c(0.8,2.5), ylab="RR",
    xlab="Temperature percentile", xaxt="n", yaxt="n", main=metalab[i])
  axis(1, at=tmeanper[paste(labperc, ".0%", sep="")], 
    labels=paste(labperc, "%", sep=""))
  axis(2, at=-1:7/5+1, las=1)
  abline(v=tmeanper[paste(labperc, ".0%", sep="")], col=grey(0.85), lty=3)
  abline(h=-1:7/5+1, col=grey(0.85), lty=3)
  lines(cpall1, ci="area", col=2, ci.arg=list(col=alpha(2,0.2)))
  lines(cpall2, ci="area", col=4, ci.arg=list(col=alpha(4,0.2)))
  legend("top", c("1st percentile", "99th pecentile"), lty=1, col=c(2,4),
    bty="n", ncol=2, inset=0.05, cex=1)
  
  # COMPUTE AND ADD THE TEST P-VALUE
  # NB: EXTRACT RELATED PARTS FROM ORIGINAL COEF/VCOV
  inds <- (1:vardf) + (i - 1) * vardf
  waldstat <- backcoef[i,,drop=F] %*% solve(backvcov[inds,inds]) %*%
    backcoef[i,]
  pval <- 1 - pchisq(waldstat, vardf)
  psign <- ifelse(pval<0.0005, " < ", " = ")
  text(mean(range(tmeanper)), 2.2, adj=c(0.5, 0.5), cex=1,
    labels=paste("p-value", formatC(pval, digits=3, format="f"), sep=psign))
}
dev.off()

# MULTI-PANEL PLOT OF EXPOSURE-RESPONSES FOR EFFECT MODIFICATION
png("figures/effmod2.png", height=4500, width=3000, res=288)
layout(matrix(1:15, nrow=5, byrow=T))
par(mar=c(5,4,3,0.5), las=1, mgp=c(2.5,1,0))

# LOOP ACROSS CHARACTERISTICS
for(i in seq(ncol(avgmetavar))) {
  
  # DEFINE CENTERING POINT
  cen <- tmeanper[["95.0%"]]
  
  # EXTRACT BACK-TRANSFORMED COEF-VCOV AND MULTIPLY BY SCALING OF 10TH VS 90TH
  var <- metanames[i]
  fscale <- diff((quantile(lsoadata[,var], c(10,90)/100) - pca$center[i]) /
      pca$scale[i])
  varcoef <- backcoef[i,] * fscale
  inds <- (1:vardf) + (i - 1) * vardf
  varvcov <- backvcov[inds,inds] * fscale

  # CROSS-PREDICT
  cpall <- crosspred(bvar, coef=varcoef, vcov=varvcov, at=tmeanper,
    model.link="log", cen=cen)

  # PLOT
  plot(cpall, type="n", ci="area", ylim=c(0.96,1.20), ylab="Ratio of RR",
    xlab="Temperature percentile", xaxt="n", yaxt="n", main=metalab[i])
  axis(1, at=tmeanper[paste(labperc, ".0%", sep="")], 
    labels=paste(labperc, "%", sep=""))
  axis(2, at=-1:4/20+1, las=1)
  abline(v=tmeanper[paste(labperc, ".0%", sep="")], col=grey(0.85), lty=3)
  abline(h=-1:4/20+1, col=grey(0.85), lty=3)
  lines(cpall)

  # COMPUTE AND ADD THE TEST P-VALUE
  # NB: EXTRACT RELATED PARTS FROM ORIGINAL COEF/VCOV
  waldstat <- backcoef[i,,drop=F] %*% solve(backvcov[inds,inds]) %*%
    backcoef[i,]
  pval <- 1 - pchisq(waldstat, vardf)
  psign <- ifelse(pval<0.0005, " < ", " = ")
  text(mean(range(tmeanper)), 2.2, adj=c(0.5, 0.5), cex=1,
    labels=paste("p-value", formatC(pval, digits=3, format="f"), sep=psign))
}
dev.off()

# LOOP ACROSS CHARACTERISTICS
rrrvar <- do.call(rbind, lapply(seq(metanames), function(i) {
  
  # EXTRACT BACK-TRANSFORMED COEF-VCOV AND MULTIPLY BY SCALING OF 10TH VS 90TH
  var <- metanames[i]
  fscale <- diff((quantile(lsoadata[,var], c(10,90)/100) - pca$center[i]) /
      pca$scale[i])
  varcoef <- backcoef[i,] * fscale
  inds <- (1:vardf) + (i - 1) * vardf
  varvcov <- backvcov[inds,inds] * fscale
  
  # COMPUTE EFFECT MODIFICATION FOR HEAT/COLD USING 95TH AS CENTERING
  cpall <- crosspred(bvar, coef=varcoef, vcov=varvcov,
    at=tmeanper[c("1.0%","99.0%")], model.link="log", cen=tmeanper[c("95.0%")])

  # STORE THE RESULTS (INVERT ORDER HEAT/COLD)
  df <- data.frame(var=var, lab=metalab[i])
  df <- cbind(df, effect=c("Heat", "Cold"), rrr=rev(cpall$allRRfit),
    rrr.low=rev(cpall$allRRlow), rrr.high=rev(cpall$allRRhigh))
  rownames(df) <- NULL
  df
}))

# ADD GROUPING AND CREATE LABEL ORDERING
rrrvar <- merge(rrrvar, dfmetagp, sort=F)
for(j in c(1:3, 7)) rrrvar[[j]] <- 
  factor(rrrvar[[j]], levels=unique(rrrvar[[j]]))

# POLAR PLOT OF EFFECT MODIFICATION
pdf("figures/effmodpolar.pdf", height=7, width=10)

ggplot(rrrvar, aes(x=factor(lab, levels=unique(lab)))) +
  geom_errorbar(aes(ymin=rrr.low, ymax=rrr.high), width=0.2) +
  geom_hline(yintercept=1, lty=2) +
  geom_point(aes(y=rrr, fill=vargroup), shape=21, size=3.5) +
  scale_fill_manual(name="Categories of small-area characteristics",
    values=hue_pal()(6)[c(1:2,6,4,3,5)]) +
  ylim(0.96, 1.03) +
  coord_polar(start=0.37) +
  facet_wrap(~effect) +
  theme_minimal() +
  theme(legend.position="bottom", strip.text=element_text(size=14),
    axis.title=element_blank(), axis.text.y=element_blank(),
    axis.text.x=element_text(size=8, vjust=-0.5)) + 
  guides(fill=guide_legend(title.position="top", title.hjust=0.5, nrow=1))

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

################################################################################
# ATTRIBUTABLE FRACTION BY REGION, IMD QUINTILES, RURAL/URBAN CLASSIFICATION

# EXTRACT ATTRIBUTABLE FRACTION FOR DIFFERENT CATEGORIES
eff <- as.data.frame(rbind(anall[agegr=="all"], anreg[agegr=="all",-1], 
  animd[agegr=="all",-1], anrururb[agegr=="all",-1]))
for(i in 5:7) eff[,i] <- eff[,i]/eff[,3]*100
eff <- as.data.frame(pivot_wider(eff, names_from=2, values_from=5:7))

# DEFINE THE X-AXIS POINTS
xpoint <- c(1, 4:13, 16:20, 23:27)

# LAYOUT
png("figures/afplot.png", height=1500, width=1920, pointsize=48)
layout(1:3, heights=c(1,1,0.5))
par(mar=c(0,5,3,1), las=1, mgp=c(2.5,1,0))

# HEAT
plot(excdeath_heat~xpoint, eff, type="n", axes=F, xlab="", ylim=c(-0.1,0.5),
  ylab="Heat-related mortality\nattributable fraction (%)")
abline(h=-1:5/10, col=grey(0.8), lty=3)
arrows(x0=xpoint, x1=xpoint, y0=eff[,"95%eCIlow_heat"],
  y1=eff[,"95%eCIhigh_heat"], length=0.05, angle=90, code=3)
points(excdeath_heat~xpoint, eff, pch=16, col="red3", cex=1.2)
axis(side=2, cex.axis=0.8)

# ADD TITLES
titles <- c("Overall","Region","IMD quintiles","Rural/Urban\nclassification")
mtext(titles, side=3, at=c(1,8.5,18,25), font=2, cex=0.7, line=1)

# COLD
par(mar=c(1,5,1,1))
plot(excdeath_cold~xpoint, eff, type="n", axes=F, xlab="", ylim=c(6,16),
  ylab="Cold-related mortality\nattributable fraction (%)")
abline(h=3:8*2, col=grey(0.8), lty=3)
arrows(x0=xpoint, x1=xpoint, y0=eff[,"95%eCIlow_cold"],
  y1=eff[,"95%eCIhigh_cold"], length=0.05, angle=90, code=3)
points(excdeath_cold~xpoint, eff, pch=16, col="navy", cex=1.3)
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
