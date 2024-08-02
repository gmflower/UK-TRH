################################################################################
# EXPOSURE-RESPONSES

# EXTRACT AVERAGE DISTRIBUTION OF TEMPERATURES
tmeanper <- sapply(stage1list, "[[", "ladtmeanper") |> rowMeans()
per <- names(tmeanper)

# DEFINE EXPOSURE-RESPONSE
argvar <- list(fun=varfun, knots=tmeanper[paste0(varper, ".0%")],
  Bound=tmeanper[c("0.0%","100.0%")])
bvar <- do.call(onebasis, c(list(x=tmeanper), argvar))

#------ Plot seasonality and temperature
tmeanseas <- datatmean[, .(tmean = mean(tasmean)), by = .(doy = yday(date))]

seasevent2 <- seasevent[doy %in% unique(tmeanseas$doy),] 
seasevent2[, mult := scale(mult), by = .(cause, agegr)]

ggplot(seasevent2, aes(x = doy)) + 
  theme_classic() +
  geom_line(aes(y = scale(mult), col = agegr)) + 
  facet_wrap(~ cause) + 
  geom_line(aes(y = scale(tmean)), data = tmeanseas, col = "black") + 
  labs(y = "Scaled seasonaility") + 
  geom_vline(xintercept = yday(seq.Date(as.Date("1989-06-01"), 
    as.Date("1989-09-01"), by = "month")), linetype = 2)
ggsave("figures/seasevent.pdf")


#----- Plot all ERFs
allfserf <- rbindlist(ladcoeflist, idcol = "cause")
allfserf <- cbind(allfserf[, .(cause, LAD11CD, agegr)], 
  t(exp(bvar %*% t(as.matrix(allfserf[
    , .SD, .SDcols = patterns("b[[:digit:]]")]))))) |> 
  melt(measure.vars = patterns("%$"), variable.name = "per", value.name = "rr")
allfserf[, tmean := tmeanper[per]]

ggplot(allfserf[cause == "cvd"]) + 
  geom_line(aes(x = tmean, y = rr, group = LAD11CD), col = "grey") + 
  facet_wrap(cause ~ agegr) + 
  geom_hline(yintercept = 1)
