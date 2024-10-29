################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# MAPS
################################################################################

# ADD LAD AND REGION CODES
lsoamap1 <- merge(lsoashp1, lookup[c("LSOA11CD","LAD11CD","RGN11CD")], sort=F)
lsoamap2 <- merge(lsoashp2, lookup[c("LSOA11CD","LAD11CD","RGN11CD")], sort=F)
ladmap <- merge(ladshp, unique(lookup[c("LAD11CD","RGN11CD")]), sort=F)

# MERGE STATS TO LSOA MAPS: RR, MMT/MMP, STD RATE
# NB: THE FIRST TWO BY AGE
lsoamap1 <- merge(lsoamap1, 
  pivot_wider(rrlsoa[,1:4], names_from=2:3, values_from=4, names_prefix="rr_"),
  sort=F)
lsoamap1 <- merge(lsoamap1, 
  pivot_wider(mmtmmplsoa[,1:4], names_from=2, values_from=3:4),
  sort=F)
lsoamap1 <- merge(lsoamap1, 
  pivot_wider(stdratelsoa[,1:3], names_from=2, values_from=3,
    names_prefix="stdrate_"),
  sort=F)

# AVERAGE MMT MMP
lsoamap1$mmt <- rowMeans(st_drop_geometry(lsoamap1[,grep("mmt",names(lsoamap1))]))
lsoamap1$mmp <- rowMeans(st_drop_geometry(lsoamap1[,grep("mmp",names(lsoamap1))]))

# REPLICATE IN HIGH-RESOLUTION LSOA MAP
lsoamap2 <- merge(lsoamap2, st_drop_geometry(lsoamap1), sort=F) 

# MERGE STATS TO LAD: ONLY STDRATE
ladmap <- merge(ladmap, 
  pivot_wider(stdratelad[,1:3], names_from=2, values_from=3,
    names_prefix="stdrate_"),
  sort=F)

# AGGREGATE TO REGIONS, SIMPLIFYING THE SHAPES
regmap <- ladmap %>% 
  merge(unique(lookup[c("RGN11CD","RGN11NM")])) %>%
  group_by(RGN11CD, RGN11NM) %>%  summarise() %>%
  ms_simplify(keep=0.025, keep_shapes=T)

# BOX OF LONDON AND CAMDEN
#lndbox <- st_as_sfc(st_bbox(subset(lsoamap, lookup$RGN11CD=="E12000007")))
#cmdbox <- st_as_sfc(st_bbox(subset(lsoamap, lookup$LAD11CD=="E09000007")))

################################################################################
# STDRATEs AND MMT BY LSOA

# DEFINE THE CUT-OFF POINTS
percut <- seq(0,100,length=10)

# HEAT
cutheat <- round(quantile(lsoamap1$stdrate_heat, percut/100), 1)
ukheat <- ggplot(data=lsoamap1) + 
  geom_sf(aes(fill=cut(stdrate_heat, cutheat, inc=T)), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.2) +
  scale_fill_brewer(palette="YlOrRd") + 
  theme_void() +
  theme(legend.position="none") 
lndheat <- ggplot(data=subset(lsoamap2, RGN11CD=="E12000007")) + 
  geom_sf(aes(fill=cut(stdrate_heat, cutheat, inc=T)), col=NA) +
  geom_sf(data=subset(ladmap, RGN11CD=="E12000007"), col=1, fill=NA, size=0.2) +
  scale_fill_brewer(palette="YlOrRd", drop=F) + 
  theme_void() +
  theme(legend.position="none") 
cmdheat <- ggplot(data=subset(lsoamap2, LAD11CD=="E09000007")) + 
  geom_sf(aes(fill=cut(stdrate_heat, cutheat, inc=T)), col=1, size=0.2) +
  guides(fill=guide_colorsteps(barwidth=18, barheight=0.6, title.position="top",
    title.hjust=0.5)) + 
  scale_fill_brewer(palette="YlOrRd", drop=F,
    name="Heat-related standardised mortality rate (x 100,000)") + 
  theme_void() +
  theme(legend.position="bottom") 

cutmmt <- round(quantile(lsoamap1$mmt, percut/100), 1)
ukmmt <- ggplot(data=lsoamap1) + 
  geom_sf(aes(fill=cut(mmt, cutmmt, inc=T)), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.2) +
  scale_fill_brewer(palette="BrBG", direction=1) + 
  theme_void() + labs(title="England & Wales") +
  theme(legend.position="none", plot.title=element_text(hjust=0)) +
  annotation_scale(style="ticks", width_hint=.3, location="tr",
    pad_x=unit(1, "cm"), pad_y=unit(1, "cm"))
lndmmt <- ggplot(data=subset(lsoamap2, RGN11CD=="E12000007")) + 
  geom_sf(aes(fill=cut(mmt, cutmmt, inc=T)), col=NA) +
  geom_sf(data=subset(ladmap, RGN11CD=="E12000007"), col=1, fill=NA, size=0.2) +
  scale_fill_brewer(palette="BrBG", direction=1, drop=F) + 
  theme_void() + labs(title="Greater London") +
  theme(legend.position="none", plot.title=element_text(hjust=0)) +
  annotation_scale(style="ticks", width_hint=.3, location="tr") 
cmdmmt <- ggplot(data=subset(lsoamap2, LAD11CD=="E09000007")) + 
  geom_sf(aes(fill=cut(mmt, cutmmt, inc=T)), col=1, size=0.2) +
  guides(fill=guide_colorsteps(barwidth=18, barheight=0.6, title.position="top",
    title.hjust=0.5)) +
  scale_fill_brewer(palette="BrBG", direction=1, drop=F,
    name="Minimum mortality temperature (MMT)") + 
  theme_void() + labs(title="Borough of Camden") +
  theme(legend.position="bottom", plot.title=element_text(hjust=0)) +
  annotation_scale(style="ticks", width_hint=.3, location="tr",
    pad_x=unit(0, "cm"), pad_y=unit(0, "cm")) 

# COLD
cutcold <- round(quantile(lsoamap1$stdrate_cold, percut/100))
ukcold <- ggplot(data=lsoamap1) + 
  geom_sf(aes(fill=cut(stdrate_cold, cutcold, inc=T)), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.2) +
  scale_fill_brewer(palette="Blues") + 
  theme_void() +
  theme(legend.position="none") 
lndcold <- ggplot(data=subset(lsoamap2, RGN11CD=="E12000007")) + 
  geom_sf(aes(fill=cut(stdrate_cold, cutcold, inc=T)), col=NA) +
  geom_sf(data=subset(ladmap, RGN11CD=="E12000007"), col=1, fill=NA, size=0.2) +
  scale_fill_brewer(palette="Blues", drop=F) + 
  theme_void() +
  theme(legend.position="none") 
cmdcold <- ggplot(data=subset(lsoamap2, LAD11CD=="E09000007")) + 
  geom_sf(aes(fill=cut(stdrate_cold, cutcold, inc=T)), col=1, size=0.2) +
  guides(fill=guide_colorsteps(barwidth=18, barheight=0.6, title.position="top",
    title.hjust=0.5)) + 
  scale_fill_brewer(palette="Blues", drop=F,
    name="Cold-related standardised mortality rate (x 100,000)") + 
  theme_void() +
  theme(legend.position="bottom") 

png("figures/lsoaratemap.png", height=3400, width=4000, res=288)
ukheat + lndheat + cmdheat + 
  ukmmt + lndmmt + cmdmmt + 
  ukcold + lndcold + cmdcold +
  plot_layout(nrow=3, ncol=3, byrow=F, height=c(1,0.5,0.4))
dev.off()

pdf("figures/lsoaratemap.pdf", height=11.2, width=13.20)
ukheat + lndheat + cmdheat + 
  ukmmt + lndmmt + cmdmmt + 
  ukcold + lndcold + cmdcold +
  plot_layout(nrow=3, ncol=3, byrow=F, height=c(1,0.5,0.4))
dev.off()

################################################################################
# AGE-SPECIFIC RR BY LSOA

# HEAT
cutheat <- round(quantile(subset(rrlsoa, effect=="heat")$rr, percut[-c(1:2)]/100), 3)
heat064 <- ggplot(data=lsoamap1) + 
  geom_sf(aes(fill=cut(`rr_age064_heat`, cutheat, inc=T)), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.2) +
  scale_fill_brewer(palette="YlOrRd", drop=F) + 
  theme_void() + labs(title=paste("Aged", agelab[1])) +
  theme(legend.position="none", plot.title=element_text(hjust=1)) 
heat6574 <- ggplot(data=lsoamap1) + 
  geom_sf(aes(fill=cut(`rr_age6574_heat`, cutheat, inc=T)), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.2) +
  scale_fill_brewer(palette="YlOrRd", drop=F) + 
  theme_void() + labs(title=paste("Aged", agelab[2])) +
  theme(legend.position="none", plot.title=element_text(hjust=1)) 
heat7584 <- ggplot(data=lsoamap1) + 
  geom_sf(aes(fill=cut(`rr_age7584_heat`, cutheat, inc=T)), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.2) +
  scale_fill_brewer(palette="YlOrRd", drop=F) + 
  theme_void() + labs(title=paste("Aged", agelab[3])) +
  theme(legend.position="none", plot.title=element_text(hjust=1)) 
heat85plus <- ggplot(data=lsoamap1) + 
  geom_sf(aes(fill=cut(`rr_age85plus_heat`, cutheat, inc=T)), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.2) +
  guides(fill=guide_colorsteps(barwidth=18, barheight=0.6, title.position="top",
    title.hjust=0.5)) + 
  scale_fill_brewer(palette="YlOrRd", drop=F,
    name="Heat-related relative risk (RR)") + 
  theme_void() + labs(title=paste("Aged", agelab[4])) +
  theme(legend.position="bottom", plot.title=element_text(hjust=1)) 

# COLD
cutcold <- round(quantile(subset(rrlsoa, effect=="cold")$rr, percut[-c(1:2)]/100), 2)
cold064 <- ggplot(data=lsoamap1) + 
  geom_sf(aes(fill=cut(`rr_age064_cold`, cutcold, inc=T)), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.2) +
  scale_fill_brewer(palette="Blues", drop=F) + 
  theme_void() + 
  theme(legend.position="none", plot.title=element_text(hjust=1))
cold6574 <- ggplot(data=lsoamap1) + 
  geom_sf(aes(fill=cut(`rr_age6574_cold`, cutcold, inc=T)), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.2) +
  scale_fill_brewer(palette="Blues", drop=F) + 
  theme_void() + 
  theme(legend.position="none", plot.title=element_text(hjust=1))
cold7584 <- ggplot(data=lsoamap1) + 
  geom_sf(aes(fill=cut(`rr_age7584_cold`, cutcold, inc=T)), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.2) +
  scale_fill_brewer(palette="Blues", drop=F) + 
  theme_void() + 
  theme(legend.position="none", plot.title=element_text(hjust=1))
cold85plus <- ggplot(data=lsoamap1) + 
  geom_sf(aes(fill=cut(`rr_age85plus_cold`, cutcold, inc=T)), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.2) +
  guides(fill=guide_colorsteps(barwidth=18, barheight=0.6, title.position="top",
    title.hjust=0.5)) + 
  scale_fill_brewer(palette="Blues", drop=F,
    name="Cold-related relative risk (RR)") + 
  theme_void() + 
  theme(legend.position="bottom", plot.title=element_text(hjust=1))

png("figures/lsoarrmap.png", height=4534, width=2666, res=288)
heat064 + heat6574 + heat7584 + heat85plus +
  cold064 + cold6574 + cold7584 + cold85plus +
  plot_layout(nrow=4, ncol=2, byrow=F)
dev.off()

################################################################################
# MMP

# MMP
cutmmp <- round(quantile(lsoamap1$mmp, percut/100), 1)
ukmmp <- ggplot(data=lsoamap1) + 
  geom_sf(aes(fill=cut(mmp, cutmmp, inc=T)), col=NA) +
  geom_sf(data=ladmap, col=1, fill=NA, size=0.2) +
  scale_fill_brewer(palette="PiYG", direction=1) + 
  theme_void() + labs(title="England & Wales") +
  theme(legend.position="none", plot.title=element_text(hjust=0.5)) +
  annotation_scale(style="ticks", width_hint=.3, location="tr", 
    pad_y=unit(3, "cm"))
lndmmp <- ggplot(data=subset(lsoamap2, RGN11CD=="E12000007")) + 
  geom_sf(aes(fill=cut(mmp, cutmmp, inc=T)), col=NA) +
  geom_sf(data=subset(ladmap, RGN11CD=="E12000007"), col=1, fill=NA, size=0.2) +
  guides(fill=guide_colorsteps(barwidth=18, barheight=0.6, title.position="top",
    title.hjust=0.5)) + 
  scale_fill_brewer(palette="PiYG", direction=1, drop=F,
    name="Minimum mortality temperature percentile (MMP)") + 
  theme_void() + labs(title="Greater London") +
  theme(legend.position="bottom", plot.title=element_text(hjust=0.5)) +
  annotation_scale(style="ticks", width_hint=.3, location="tr")
cmdmmp <- ggplot(data=subset(lsoamap2, LAD11CD=="E09000007")) + 
  geom_sf(aes(fill=cut(mmp, cutmmp, inc=T)), col=1, size=0.2) +
  scale_fill_brewer(palette="PiYG", direction=1, drop=F) + 
  theme_void() + labs(title="Borough of Camden") +
  theme(legend.position="none", plot.title=element_text(hjust=0.5))  +
  annotation_scale(style="ticks", width_hint=.3, location="tr")

png("figures/lsoammpmap.png", height=2000, width=4000, res=288)
ukmmp + lndmmp + cmdmmp + 
  plot_layout(nrow=1, ncol=3, byrow=T, widths=c(1,0.7,0.6))
dev.off()

################################################################################
# ORIGINAL VARIABLES

png("figures/varmap.png", height=5000, width=3400, res=288)

row <- rep(1:5, each=3)
col <- rep(1:3, 5)

grid.newpage()
pushViewport(viewport(layout=grid.layout(5, 3)))

for(i in seq(metanames)) {
  cut <- unique(round(quantile(lsoadata[[metanames[i]]], percut/100), 3))
  map <- ggplot(data=merge(lsoamap1, lsoadata[c("LSOA11CD", metanames[i])])) + 
    geom_sf(aes(fill=cut(get(metanames[i]), cut, inc=T)), col=NA) +
    geom_sf(data=ladmap, col=1, fill=NA, size=0.2) +
    guides(fill=guide_colorsteps(barwidth=0.6, barheight=10)) + 
    scale_fill_brewer(palette=i, name="") + 
    theme_void() + labs(title=metalab[i]) +
    theme(plot.title=element_text(hjust=1))
  print(map, vp=viewport(layout.pos.row=row[i], layout.pos.col=col[i]))
}

dev.off()

################################################################################
# PRINCIPAL COMPONENT

# DERIVE PCA COMPONENTS AT LSOA LEVEL
lsoavar <- predict(pca, newdata=lsoadata)[,seq(ncomp)]

png("figures/pcamap.png", height=2800, width=2800, res=288)

row <- c(1,1,2)
col <- list(1:2,3:4,2:3)

grid.newpage()
pushViewport(viewport(layout=grid.layout(2, 4)))

for(i in seq(ncol(lsoavar))) {
  cut <- unique(round(quantile(lsoavar[,i], percut/100), 3))
  map <- ggplot(data=lsoamap1) + 
    geom_sf(aes(fill=cut(lsoavar[,i], cut, inc=T)), col=NA) +
    geom_sf(data=ladmap, col=1, fill=NA, size=0.4) +
    guides(fill=guide_colorsteps(barwidth=0.6, barheight=10)) + 
    scale_fill_brewer(palette=i, name="") + 
    theme_void() + labs(title=paste("Principal component", i)) +
    theme(plot.title=element_text(hjust=1))
  print(map, vp=viewport(layout.pos.row=row[[i]], layout.pos.col=col[[i]]))
}

dev.off()

################################################################################
# MAP OF REGIONS

png("figures/regmap.png", height=2000, width=2000, res=288)

ggplot(data=regmap) + 
  geom_sf(aes(fill=factor(RGN11NM, levels=RGN11NM)), col=1) +
  scale_fill_brewer(palette="Paired", name="Region") +
  theme_void() 

dev.off()