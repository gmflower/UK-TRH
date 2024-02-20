################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# PLOTS
################################################################################

################################################################################
# EXPOSURE-RESPONSES

# EXTRACT DISTRIBUTION OF TEMPERATURES
tmeandist <- t(sapply(estlist, function(x) x$tmeanper))
summary(tmeandist[,c(1,ncol(tmeandist))])
tmeanper <- colMeans(tmeandist)
per <- names(tmeanper)

# DEFINE EXPOSURE-RESPONSE
argvar <- list(fun=varfun, knots=tmeanper[c("10.0%","75.0%","90.0%")],
  Bound=tmeanper[c("0.0%","100.0%")])
argvar$degree <- vardegree
bvar <- do.call(onebasis, c(list(x=tmeanper), argvar))

# PREDICT PARAMETERS FOR AVERAGE CURVES (NB: ALL THE PCA COMPONENTS ARE 0 BY DEF)
meanladvar <- as.data.frame(t(colMeans(ladvar)))
predparall <- predict(metaall, newdata=meanladvar, vcov=T)

# MULTI-PANEL PLOT
pdf("figures/expresp.pdf", height=5, width=7)
par(mar=c(5,4,1,0.5), las=1, mgp=c(2.5,1,0))

# OBTAIN PREDICTIONS
cen <- tmeanper[[which.min(bvar%*%predparall$fit)]]
cpall <- crosspred(bvar, coef=predparall$fit, vcov=predparall$vcov, at=tmeanper,
  model.link="log", cen=cen)
cpallcold <- crosspred(bvar, coef=predparall$fit, vcov=predparall$vcov,
  at=tmeanper[tmeanper<=cen], model.link="log", cen=cen)
cpallheat <- crosspred(bvar, coef=predparall$fit, vcov=predparall$vcov,
  at=tmeanper[tmeanper>=cen], model.link="log", cen=cen)

# LAD-SPECIFIC AND POOLED
plot(cpall, type="n", ci="n", ylim=c(0.8,2.5), ylab="RR", xlab="Temperature percentiles", xaxt="n")
for(i in seq(listlad))
  lines(crosspred(bvar, coef=coefall[i,], vcov=vcovall[[i]],
    model.link="log", cen=cen, at=tmeanper), col=grey(0.85), lty=5)
abline(h=1)
labperc <- paste(c(0,1,25,50,75,99,100),sep="")
axis(1, at=tmeanper[paste(labperc, ".0%", sep="")], 
  labels=paste(labperc, "%", sep=""))
abline(v=tmeanper["1.0%"], lty="dashed", col="blue")
abline(v=tmeanper["99.0%"], lty="dashed", col="red")
lines(cpallcold, col=4, ci="area", ci.arg=list(col=alpha("blue", 0.2)))
lines(cpallheat, col=2, ci="area", ci.arg=list(col=alpha("red", 0.2)))

dev.off()

################################################################################
# META-VARIABLES AND PCA

# CORRELATION MATRIX PLOT
corvarlad <- cor(metavar)
corvarlsoa <- cor(lsoadata[metanames])
dimnames(corvarlad) <- dimnames(corvarlsoa) <- list(metalab, metalab)
png("figures/corvar.png", height=800, width=1600, pointsize=24)
layout(t(1:2))
corrplot(corvarlad, tl.cex=0.5, mar=c(0,1,2,1), title="LAD level")
corrplot(corvarlsoa, tl.cex=0.5, mar=c(0,1,2,1), title="LSOA level")
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
abline(h=cumpro[3], lty=5)
dev.off()

# EXTRACT COORDINATES (CORRELATIONS)
pcacoord <- t(apply(pca$rotation,1, "*", pca$sdev))
rownames(pcacoord) <- metalab

# PLOTS COORDINATES
pdf("figures/pcascore1.pdf", height=5, width=8)
layout(1:4, heights=c(1,1,1,0.7))
par(mar=c(0,4,4,1), las=1, mgp=c(2.5,1,0))
for(i in seq(1:3)) {
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

pdf("figures/pcascore2.pdf", height=3.7, width=8.5)
as.data.frame(pcacoord[,1:3]) %>%
  cbind(lab=factor(metalab, levels=metalab)) %>%
  pivot_longer(1:3) %>%
  ggplot() +
  geom_tile(aes(x=lab, y=name, fill=value), col="white") +
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
# EXCESS DEATHS

# COMPUTE ATTRIBUTABLE FRACTION FOR DIFFERENT CATEGORIES
af <- rbind(allres, regres[,-c(1,2)], imdres[,-1], rururbres[,-1])
for(i in 2:ncol(af)) af[i] <- af[i]/af[1]*100
af <- af[,-1]

# DEFINE THE X-AXIS POINTS
xpoint <- c(1, 4:13, 16:20, 23:27)

# LAYOUT
pdf("figures/afplot.pdf", height=5.5, width=7)
layout(1:3, heights=c(1,1,0.5))
par(mar=c(0,4,3,1), las=1, mgp=c(2.5,1,0))

# HEAT
plot(ANheat~xpoint, af, type="n", axes=F, xlab="", ylim=c(-0.1,0.5),
  ylab="Heat-related excess mortality (%)")
abline(h=-1:5/10, col=grey(0.8), lty=3)
arrows(x0=xpoint, x1=xpoint, y0=af[,"ANheat_95%eCIlow"],
  y1=af[,"ANheat_95%eCIhigh"], length=0.05, angle=90, code=3)
points(ANheat~xpoint, af, pch=16, col="red3", cex=1.2)
axis(side=2, cex.axis=0.8)

# ADD TITLES
titles <- c("Overall","Region","IMD quintiles","Rural/Urban\nclassification")
mtext(titles, side=3, at=c(1,8.5,18,25), font=2, cex=0.7, line=1)

# COLD
par(mar=c(1,4,0,1))
plot(ANcold~xpoint, af, type="n", axes=F, xlab="", ylim=c(6,17),
  ylab="Cold-related excess mortality (%)")
abline(h=3:8*2, col=grey(0.8), lty=3)
arrows(x0=xpoint, x1=xpoint, y0=af[,"ANcold_95%eCIlow"],
  y1=af[,"ANcold_95%eCIhigh"], length=0.05, angle=90, code=3)
points(ANcold~xpoint, af, pch=16, col="navy", cex=1.3)
axis(side=2, cex.axis=0.8)

# DEFINE LABELS
labels <- c(regres$RGN11NM,"1 - most deprived","2","3","4","5 - least deprived",
  levels(rururbres$RurUrb))

# ADD LABELS
par(mar=c(0,4,0,1))
plot(xpoint[-1], xpoint[-1], type="n", ylim=c(0,1), axes=F, ylab="", xlab="", bty="n",
  xlim=par()$usr[1:2], xaxs="i")
text(xpoint[-1], 1, srt=50, adj=1, xpd=TRUE, labels=labels, cex=0.9)

layout(1)
dev.off()

# SAVE AS 5X6 PDF
