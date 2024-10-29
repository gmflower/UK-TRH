################################################################################
# UK-TRH: SMALL-AREA ANALYSIS OF TEMPERATURE RELATED HOSPITALISATIONS IN ENGLAND
################################################################################

#### Forest plots of second stage results

# Temperature reference points:
ref_temp <- quantile(datatmean$tmean, 50/100, na.rm = T)
heat_temp <- quantile(datatmean$tmean, 99/100, na.rm = T)

# Define groups of diagnoses:
cause_groups <- list(
  list(group_name="Cardiovascular", causes=c("cvd","mi","pulmheart","stroke")),
  list(group_name="Respiratory", causes=c("resp","ari","pneumonia","copd","asthma")),
  list(group_name="Genitourinary", causes=c("genito","renal","arf"))
)
names(cause_groups) <- c("Cardiovascular","Respiratory","Genitourinary")
  

names(cause_groups)[1]
causes <- cause_groups[[1]][["causes"]]

Temperature <- c(-5,-5,-5)
RR  <- c(cp_allcause$allRRfit["-5"],cp_cvd$allRRfit["-5"],cp_resp$allRRfit["-5"]) 
RR_lowCI <- c(cp_allcause$allRRlow["-5"],cp_cvd$allRRlow["-5"],cp_resp$allRRlow["-5"])
RR_upperCI <- c(cp_allcause$allRRhigh["-5"],cp_cvd$allRRhigh["-5"],cp_resp$allRRhigh["-5"])
Group <- c("All Cause","Cardiovascular","Respiratory")

cold <- data.frame(Temperature, RR, RR_lowCI, RR_upperCI, Group)

Temperature <- c(25,25,25)
RR  <- c(cp_allcause$allRRfit["25"],cp_cvd$allRRfit["25"],cp_resp$allRRfit["25"]) 
RR_lowCI <- c(cp_allcause$allRRlow["25"],cp_cvd$allRRlow["25"],cp_resp$allRRlow["25"])
RR_upperCI <- c(cp_allcause$allRRhigh["25"],cp_cvd$allRRhigh["25"],cp_resp$allRRhigh["25"])
Group <- c("All Cause","Cardiovascular","Respiratory")

hot <- data.frame(Temperature, RR, RR_lowCI, RR_upperCI, Group)

cold_fp <- ggplot(data=cold, aes(x=Group, y=RR, ymin=RR_lowCI, ymax=RR_upperCI, 
                                  group=Group , color=Group)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  xlab("") + ylab("Relative Risk") +
  ylim(0.95,2.50) +
  ggtitle("Extreme Cold") +
  scale_colour_manual(values = c("Black","red","#4472C4")) +
  theme_bw() + # use a white background
  theme(plot.title = element_text(hjust = 0.5)) 

hot_fp <- ggplot(data=hot, aes(x=Group, y=RR, ymin=RR_lowCI, ymax=RR_upperCI, group=Group , color=Group)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  xlab("") + ylab("Relative Risk") +
  ylim(0.95,2.50) +
  ggtitle("Extreme Heat") +
  scale_colour_manual(values = c("Black","red","#4472C4")) +
  theme_bw() + # use a white background
  theme(plot.title = element_text(hjust = 0.5)) 

#grid.arrange(cold_fp, hot_fp, ncol=2, top = textGrob("England - Temperature Extremes", gp=gpar(fontsize=12)))
combined_plot <- ggarrange(cold_fp, hot_fp, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
combined_plot
annotate_figure(combined_plot, top = text_grob("England", 
                                      color = "black", face = "bold", size = 14))
