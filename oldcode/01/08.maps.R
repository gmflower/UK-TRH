################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# MAPS
################################################################################

# MERGE SHAPEFILES AND RESULTS
lsoamap <- merge(lsoashp, lsoares, sort=F)
ladmap <- merge(ladshp, ladres, sort=F)

# ADD LAD AND REGION CODES
lsoamap <- merge(lsoamap, lookup[c("LSOA11CD","LAD11CD","RGN11CD")], sort=F)
ladmap <- merge(ladmap, unique(lookup[c("LAD11CD","RGN11CD")]), sort=F)

# BOX OF LONDON AND CAMDEN
#lndbox <- st_as_sfc(st_bbox(subset(lsoamap, lookup$RGN11CD=="E12000007")))
#cmdbox <- st_as_sfc(st_bbox(subset(lsoamap, lookup$LAD11CD=="E09000007")))

################################################################################
# TMAP

# ukheat <- tm_shape(lsoamap) + 
#   tm_fill("RR99", palette="YlOrRd", title="RR Heat", style="quantile", n=10) + 
#   tm_shape(ladmap) + tm_borders(col="black") +
#   tm_layout(frame=F)
# 
# ukmmt <- tm_shape(lsoamap) + 
#   tm_fill("mmt", palette="PiYG", title="MMT (C)", style="quantile", n=10) + 
#   tm_shape(ladmap) + tm_borders(col="black") +
#   tm_layout(frame=F)
# 
# ukcold <- tm_shape(lsoamap) + 
#   tm_fill("RR01", palette="Blues", title="RR cold", style="quantile", n=10) + 
#   tm_shape(ladmap) + tm_borders(col="black") +
#   tm_layout(frame=F)
# 
# lndheat <- tm_shape(subset(lsoamap, lookup$RGN11CD=="E12000007")) + 
#   tm_fill("RR99", palette="YlOrRd", title="RR Heat", style="quantile", n=10) + 
#   tm_shape(subset(ladmap, lookup$RGN11CD=="E12000007")) + tm_borders(col="black") +
#   tm_layout(frame=F)
# 
# lndmmt <- tm_shape(subset(lsoamap, lookup$RGN11CD=="E12000007")) + 
#   tm_fill("mmt", palette="PiYG", title="MMT (C)", style="quantile", n=10) + 
#   tm_shape(subset(ladmap, lookup$RGN11CD=="E12000007")) + tm_borders(col="black") +
#   tm_layout(frame=F)
# 
# lndcold <- tm_shape(subset(lsoamap, lookup$RGN11CD=="E12000007")) + 
#   tm_fill("RR01", palette="Blues", title="RR cold", style="quantile", n=10) + 
#   tm_shape(subset(ladmap, lookup$RGN11CD=="E12000007")) + tm_borders(col="black") +
#   tm_layout(frame=F)
# 
# 
# pdf("maptmap.pdf", height=10, width=16)
# grid.newpage()
# pushViewport(viewport(layout=grid.layout(2, 3, height=c(7,3))))
# print(ukheat, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(ukmmt, vp=viewport(layout.pos.row=1, layout.pos.col=2))
# print(ukcold, vp=viewport(layout.pos.row=1, layout.pos.col=3))
# print(lndheat, vp=viewport(layout.pos.row=2, layout.pos.col=1))
# print(lndmmt, vp=viewport(layout.pos.row=2, layout.pos.col=2))
# print(lndcold, vp=viewport(layout.pos.row=2, layout.pos.col=3))
# dev.off()

################################################################################
# RRs AND MMT BY LSOA

# DEFINE THE CUT-OFF POINTS
percut <- seq(0,100,length=10)

# HEAT
cutheat <- quantile(lsoamap$RR99, percut/100)
ukheat <- ggplot(data=lsoamap) + 
  geom_sf(aes(fill=cut(RR99, cutheat, inc=T)), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.4) +
  scale_fill_brewer(palette="YlOrRd", name="RR heat") + 
  theme_void() +
  theme(legend.position="none", plot.title=element_text(hjust=1)) 
lndheat <- ggplot(data=subset(lsoamap, RGN11CD=="E12000007")) + 
  geom_sf(aes(fill=cut(RR99, cutheat, inc=T)), col=NA) +
  geom_sf(data=subset(ladmap, RGN11CD=="E12000007"), col=1, fill=NA, size=0.4) +
  scale_fill_brewer(palette="YlOrRd", name="RR heat", drop=F) + 
  theme_void() +
  theme(legend.position="none", plot.title=element_text(hjust=1)) 
cmdheat <- ggplot(data=subset(lsoamap, LAD11CD=="E09000007")) + 
  geom_sf(aes(fill=cut(RR99, cutheat, inc=T)), col=1, size=0.4) +
  guides(fill=guide_colorsteps(barwidth=18, barheight=0.6, title.position="top",
    title.hjust=0.5)) + 
  scale_fill_brewer(palette="YlOrRd", name="RR heat", drop=F) + 
  theme_void() +
  theme(legend.position="bottom", plot.title=element_text(hjust=1)) 

# MMT
cutmmt <- quantile(lsoamap$mmt, percut/100)
ukmmt <- ggplot(data=lsoamap) + 
  geom_sf(aes(fill=cut(mmt, cutmmt, inc=T)), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.4) +
  scale_fill_brewer(palette="PiYG", direction=-1, name="MMT (C)") + 
  theme_void() + labs(title="England & Wales") +
  theme(legend.position="none", plot.title=element_text(hjust=1)) 
lndmmt <- ggplot(data=subset(lsoamap, RGN11CD=="E12000007")) + 
  geom_sf(aes(fill=cut(mmt, cutmmt, inc=T)), col=NA) +
  geom_sf(data=subset(ladmap, RGN11CD=="E12000007"), col=1, fill=NA, size=0.4) +
  scale_fill_brewer(palette="PiYG", direction=-1, name="MMT (C)", drop=F) + 
  theme_void() + labs(title="Greater London") +
  theme(legend.position="none", plot.title=element_text(hjust=1)) 
cmdmmt <- ggplot(data=subset(lsoamap, LAD11CD=="E09000007")) + 
  geom_sf(aes(fill=cut(mmt, cutmmt, inc=T)), col=1, size=0.4) +
  guides(fill=guide_colorsteps(barwidth=18, barheight=0.6, title.position="top",
    title.hjust=0.5)) + 
  scale_fill_brewer(palette="PiYG", direction=-1, name="MMT (C)", drop=F) + 
  theme_void() + labs(title="Camden") +
  theme(legend.position="bottom", plot.title=element_text(hjust=1)) 

# COLD
cutcold <- quantile(lsoamap$RR01, percut/100)
ukcold <- ggplot(data=lsoamap) + 
  geom_sf(aes(fill=cut(RR01, cutcold, inc=T)), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.4) +
  scale_fill_brewer(palette="Blues", name="RR cold") + 
  theme_void() +
  theme(legend.position="none", plot.title=element_text(hjust=1)) 
lndcold <- ggplot(data=subset(lsoamap, RGN11CD=="E12000007")) + 
  geom_sf(aes(fill=cut(RR01, cutcold, inc=T)), col=NA) +
  geom_sf(data=subset(ladmap, RGN11CD=="E12000007"), col=1, fill=NA, size=0.4) +
  scale_fill_brewer(palette="Blues", name="RR cold", drop=F) + 
  theme_void() +
  theme(legend.position="none", plot.title=element_text(hjust=1)) 
cmdcold <- ggplot(data=subset(lsoamap, LAD11CD=="E09000007")) + 
  geom_sf(aes(fill=cut(RR01, cutcold, inc=T)), col=1, size=0.4) +
  guides(fill=guide_colorsteps(barwidth=18, barheight=0.6, title.position="top",
    title.hjust=0.5)) + 
  scale_fill_brewer(palette="Blues", name="RR cold", drop=F) + 
  theme_void() +
  theme(legend.position="bottom", plot.title=element_text(hjust=1)) 

# SAVE
# pdf("figures/lsoarrmap.pdf", height=20, width=10)
# ukheat + lndheat + cmdheat + 
#   ukmmt + lndmmt + cmdmmt + 
#   ukcold + lndcold + cmdcold +
#   plot_layout(nrow=3, ncol=3, byrow=F, height=c(1,0.5,0.4))
# dev.off()

png("figures/lsoarrmap.png", height=1700, width=2000, res=144)
ukheat + lndheat + cmdheat + 
  ukmmt + lndmmt + cmdmmt + 
  ukcold + lndcold + cmdcold +
  plot_layout(nrow=3, ncol=3, byrow=F, height=c(1,0.5,0.4))
dev.off()

################################################################################
# RRs AND MMT BY LAD

# HEAT
ukheat <- ggplot(data=ladmap) + 
  geom_sf(aes(fill=cut(RR99, cutheat, inc=T)), col=1, size=0.4) +
  scale_fill_brewer(palette="YlOrRd", name="RR heat") + 
  theme_void() +
  theme(legend.position="none", plot.title=element_text(hjust=1)) 
lndheat <- ggplot(data=subset(ladmap, RGN11CD=="E12000007")) + 
  geom_sf(aes(fill=cut(RR99, cutheat, inc=T)), col=1, size=0.4) +
  scale_fill_brewer(palette="YlOrRd", name="RR heat", drop=F) + 
  theme_void() +
  theme(legend.position="none", plot.title=element_text(hjust=1)) 
cmdheat <- ggplot(data=subset(ladmap, LAD11CD=="E09000007")) + 
  geom_sf(aes(fill=cut(RR99, cutheat, inc=T)), col=1, size=0.4) +
  guides(fill=guide_colorsteps(barwidth=18, barheight=0.6, title.position="top",
    title.hjust=0.5)) + 
  scale_fill_brewer(palette="YlOrRd", name="RR heat", drop=F) + 
  theme_void() +
  theme(legend.position="bottom", plot.title=element_text(hjust=1)) 

# MMT
ukmmt <- ggplot(data=ladmap) + 
  geom_sf(aes(fill=cut(mmt, cutmmt, inc=T)), col=1, size=0.4) +
  scale_fill_brewer(palette="PiYG", direction=-1, name="MMT (C)") + 
  theme_void() + labs(title="England & Wales") +
  theme(legend.position="none", plot.title=element_text(hjust=1)) 
lndmmt <- ggplot(data=subset(ladmap, RGN11CD=="E12000007")) + 
  geom_sf(aes(fill=cut(mmt, cutmmt, inc=T)), col=1, size=0.4) +
  scale_fill_brewer(palette="PiYG", direction=-1, name="MMT (C)", drop=F) + 
  theme_void() + labs(title="Greater London") +
  theme(legend.position="none", plot.title=element_text(hjust=1))
cmdmmt <- ggplot(data=subset(ladmap, LAD11CD=="E09000007")) + 
  geom_sf(aes(fill=cut(mmt, cutmmt, inc=T)), col=1, size=0.4) +
  guides(fill=guide_colorsteps(barwidth=18, barheight=0.6, title.position="top",
    title.hjust=0.5)) + 
  scale_fill_brewer(palette="PiYG", direction=-1, name="MMT (C)", drop=F) + 
  theme_void() + labs(title="Camden") +
  theme(legend.position="bottom", plot.title=element_text(hjust=1))

# COLD
ukcold <- ggplot(data=ladmap) + 
    geom_sf(aes(fill=cut(RR01, cutcold, inc=T)), col=1, size=0.4) +
    scale_fill_brewer(palette="Blues", name="RR cold") + 
    theme_void() +
  theme(legend.position="none", plot.title=element_text(hjust=1)) 
lndcold <- ggplot(data=subset(ladmap, RGN11CD=="E12000007")) + 
  geom_sf(aes(fill=cut(RR01, cutcold, inc=T)), col=1, size=0.4) +
  scale_fill_brewer(palette="Blues", name="RR cold", drop=F) + 
  theme_void() +
  theme(legend.position="none", plot.title=element_text(hjust=1)) 
cmdcold <- ggplot(data=subset(ladmap, LAD11CD=="E09000007")) + 
  geom_sf(aes(fill=cut(RR01, cutcold, inc=T)), col=1, size=0.4) +
  guides(fill=guide_colorsteps(barwidth=18, barheight=0.6, title.position="top",
    title.hjust=0.5)) + 
  scale_fill_brewer(palette="Blues", name="RR cold", drop=F) + 
  theme_void() +
  theme(legend.position="bottom", plot.title=element_text(hjust=1)) 

# SAVE
# pdf("figures/ladrrmap.pdf", height=20, width=10)
# ukheat + lndheat + cmdheat + 
#   ukmmt + lndmmt + cmdmmt + 
#   ukcold + lndcold + cmdcold +
#   plot_layout(nrow=3, ncol=3, byrow=F, height=c(1,0.5,0.4))
# dev.off()

png("figures/ladrrmap.png", height=1700, width=2000, res=144)
ukheat + lndheat + cmdheat + 
  ukmmt + lndmmt + cmdmmt + 
  ukcold + lndcold + cmdcold +
  plot_layout(nrow=3, ncol=3, byrow=F, height=c(1,0.5,0.4))
dev.off()

################################################################################
# ORIGINAL VARIABLES

png("figures/varmap.png", height=2500, width=1400, res=144)

row <- rep(1:5, each=3)
col <- rep(1:3, 5)

grid.newpage()
pushViewport(viewport(layout=grid.layout(5, 3)))

for(i in seq(metanames)) {
  cut <- unique(round(quantile(lsoadata[[metanames[i]]], percut/100), 3))
  map <- ggplot(data=merge(lsoamap, lsoadata[c("LSOA11CD", metanames[i])])) + 
    geom_sf(aes(fill=cut(get(metanames[i]), cut, inc=T)), col=NA) +
    geom_sf(data=ladmap, col=1, fill=NA, size=0.4) +
    guides(fill=guide_colorsteps(barwidth=0.6, barheight=10)) + 
    scale_fill_brewer(palette=i, name="") + 
    theme_void() + labs(title=metalab[i]) +
    theme(plot.title=element_text(hjust=1))
  print(map, vp=viewport(layout.pos.row=row[i], layout.pos.col=col[i]))
}

dev.off()

################################################################################
# PRINCIPAL COMPONENT

png("figures/pcamap.png", height=1400, width=1400, res=144)

row <- c(1,1,2)
col <- list(1:2,3:4,2:3)

grid.newpage()
pushViewport(viewport(layout=grid.layout(2, 4)))

for(i in seq(ncol(lsoavar))) {
  cut <- unique(round(quantile(lsoavar[,i], percut/100), 3))
  map <- ggplot(data=lsoamap) + 
    geom_sf(aes(fill=cut(lsoavar[,i], cut, inc=T)), col=NA) +
    geom_sf(data=ladmap, col=1, fill=NA, size=0.4) +
    guides(fill=guide_colorsteps(barwidth=0.6, barheight=10)) + 
    scale_fill_brewer(palette=i, name="") + 
    theme_void() + labs(title=paste("Principal component", i)) +
    theme(plot.title=element_text(hjust=1))
  print(map, vp=viewport(layout.pos.row=row[[i]], layout.pos.col=col[[i]]))
}

dev.off()



