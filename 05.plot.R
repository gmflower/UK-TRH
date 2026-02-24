################################################################################
# UK-TRM: SMALL-AREA ANALYSIS OF TEMPERATURE-MORTALITY IN ENGLAND & WALES
################################################################################

################################################################################
# PLOTS (REPRODUCIBLE)
################################################################################

################################################################################
# HOSPITAL ADMISSIONS BY AGE

names <- unlist(lapply(seq(causelist), function(x) causelist[[x]][["name"]]))
shortnames <- unlist(lapply(seq(causelist), function(x) causelist[[x]][["shortname"]]))
cause_ref_list <- data.frame(names, shortnames)

cause_ref_list$grouping[cause_ref_list$shortnames %in% c("resp","ari","pneumonia","copd","asthma")] <- 5 
cause_ref_list$grouping[cause_ref_list$shortnames %in% c("cvd","mi","stroke","hf", "hypo")] <- 4 
cause_ref_list$grouping[cause_ref_list$shortnames %in% c("genito","renal","arf","ckd")] <- 3
cause_ref_list$grouping[cause_ref_list$shortnames %in% c("infect","bacterial")] <- 2
cause_ref_list$grouping[cause_ref_list$shortnames %in% c("endo","diabetes", "metabolic")] <- 1 

cause_ref_list$names[names=="Respiratory"] <- "Respiratory (All)"
cause_ref_list$names[names=="Cardiovascular"] <- "Cardiovascular (All)"
cause_ref_list$names[names=="Genitourinary"] <- "Genitourinary (All)"
cause_ref_list$names[names=="Infectious and parasitic"] <- "Infectious and parasitic (All)"
cause_ref_list$names[names=="Endocrine, nutritional, metabolic"] <- "Endocrine, nutritional, metabolic (All)"

hes_aggregated <- hesdata |> 
  filter(agegr=="total" & cause %in% setcause) |> 
  group_by(cause) |> 
  summarise(count = sum(count)) |> 
  left_join(cause_ref_list, by = join_by(cause==shortnames)) |> 
  as.data.frame()
hes_aggregated <- arrange(hes_aggregated,count)
leadingcauses <- unique(hes_aggregated$cause)

hes_byage <- hesdata |> 
  filter(agegr!="total" & cause %in% setcause) |> 
  group_by(cause, agegr) |> 
  summarise(count = sum(count)) |> 
  left_join(cause_ref_list, by = join_by(cause==shortnames)) |> 
  as.data.table()

# Keep only leading causes:
hes_byage <- hes_byage[cause %in% leadingcauses,]

pdf("./figures/total_admission_burden.pdf", height=9, width=9.5)

fig1 <- ggplot(data=hes_byage, aes(fill=agegr, x=count/1000, y=reorder(reorder(names,count),grouping))) +
  # First geom_bar: original bars (no outlines)
  geom_bar(position="stack", stat="identity", color=NA) +  
  # Second geom_bar: bars with black outlines
  geom_bar(data = hes_byage |> filter(names %in% c("Respiratory (All)", "Cardiovascular (All)", 
                                                   "Genitourinary (All)", "Infectious and parasitic (All)", 
                                                   "Endocrine, nutritional, metabolic (All)")), 
           aes(fill=agegr, x=count/1000, y=reorder(reorder(names,count),grouping)), 
           position="stack", stat="identity", color="black", size=0.5, width=0.9) +  # Add black outline, same width as others
  scale_fill_manual(values=c(paletteer_c("ggthemes::Red", 4))) + 
  labs(
    x="Thousand admissions (x 1,000)",
    y="") +
  theme_bw() +
  theme(legend.position = "bottom",plot.margin = unit(c(1,10,1,1), "lines")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_continuous(labels = scales::comma_format(big.mark = ","), breaks = c(0:7*200)) +  
  ggtitle("Unplanned Hospital Admissions", 
          subtitle="England from 2008 to 2019 (June to August)") +
  guides(fill=guide_legend(title="Age group")) +
  annotation_custom(
    grob = textGrob("Respiratory", gp = gpar(fontsize = 8)), 
    xmin = 3150, xmax = 7.5, ymin = 31, ymax = 1   # Position of the label for "Group 1"
  ) +
  annotation_custom(
    grob = textGrob("Cardiovascular", gp = gpar(fontsize = 8)),
    xmin = 3150, xmax = 7.5, ymin = 20, ymax = 2  # Position of the label for "Group 2"
  ) +
  annotation_custom(
    grob = textGrob("Genitourinary", gp = gpar(fontsize = 8)),
    xmin = 3150, xmax = 7.5, ymin = 11, ymax = 3  # Position of the label for "Group 3"
  ) +
  annotation_custom(
    grob = textGrob("Infectious &\nParasitic", gp = gpar(fontsize = 8)),
    xmin = 3150, xmax = 7.5, ymin = 5, ymax = 4  # Position of the label for "Group 4"
  ) +
  annotation_custom(
    grob = textGrob("Endocrine, Nutritional &\nMetabolic", gp = gpar(fontsize = 8)),
    xmin = 3150, xmax = 7.5, ymin = 0, ymax = 4  # Position of the label for "Group 5"
  ) +
  coord_cartesian(clip = 'off')
fig1 

dev.off()

################################################################################
# FOREST PLOT

# SET TEMPERATURE REF POINTS
heat_temp <- round(quantile(datatmean$tmean, 99/100, na.rm = T),1)

names <- unlist(lapply(seq(causelist), function(x) causelist[[x]][["name"]]))
shortnames <- unlist(lapply(seq(causelist), function(x) causelist[[x]][["shortname"]]))
cause_ref_list <- data.frame(names, shortnames)

rr <- rr_lowci <- rr_upci <- cause <- all_plot_data <- NULL
ages <- c("total","age1864","age6574","age7584","age85plus") # PLOT ORDER
plot_data <- NULL

for (c in seq(setcause)) {
  rr <- rr_lowci <- rr_upci <- cause <- NULL
  for (a in seq(ages)) {
    
    rr <- c(rr,pooled_results_list_c[[setcause[c]]][[ages[a]]][["allRRfit"]][[paste0(heat_temp)]])
    rr_lowci <- c(rr_lowci,pooled_results_list_c[[setcause[c]]][[ages[a]]][["allRRlow"]][[paste0(heat_temp)]])
    rr_upci <- c(rr_upci,pooled_results_list_c[[setcause[c]]][[ages[a]]][["allRRhigh"]][[paste0(heat_temp)]])
    cause <- setcause[c]
  }
  
  plot_data <- data.frame(rr, rr_lowci, rr_upci, ages, cause)
  all_plot_data <- rbind(all_plot_data,plot_data)
}

all_plot_data <- left_join(all_plot_data,cause_ref_list,  by = join_by(cause == shortnames))


all_plot_data$grouping <- NULL
all_plot_data$grouping[all_plot_data$cause %in% c("resp","ari","pneumonia","copd","asthma")] <- "Respiratory" 
all_plot_data$grouping[all_plot_data$cause %in% c("cvd","mi","pulmheart","stroke","ihd","hf", "hypo")] <- "Cardiovascular" 
all_plot_data$grouping[all_plot_data$cause %in% c("genito","renal","arf","ckd")] <- "Genitourinary"
all_plot_data$grouping[all_plot_data$cause %in% c("infect", "bacterial")] <- "Infectious and parasitic" 
all_plot_data$grouping[all_plot_data$cause %in% c("endo","diabetes", "metabolic")] <- "Endocrine, nutritional, metabolic"

# SET ORDER
all_plot_data <-all_plot_data[with(all_plot_data, 
                                   order(grouping, names, cause, ages)),]
all_plot_data <- all_plot_data[, c("grouping", "names", "cause", "ages", "rr", "rr_lowci", "rr_upci")]
all_plot_data$id <- as.numeric(factor(all_plot_data$cause, 
                                      levels = c("cvd","hf","hypo","mi","stroke",
                                                 "endo","diabetes","metabolic","genito","arf","renal","infect","bacterial",
                                                 "resp","ari","asthma","copd","pneumonia"))) +
  as.numeric(factor(all_plot_data$grouping, 
                    levels = unique(all_plot_data$grouping)))
all_plot_data$ages <- factor(all_plot_data$ages, levels = ages)
all_plot_data$age_lab <- "first"
all_plot_data$age_lab[all_plot_data$ages=="total"] <- "second"


# SET GRID
unid <- unique(all_plot_data$id)
bglines <- data.frame(pos = unid[-1][diff(unid)==1] - .5)

# AGE GROUP SHAPES
all_plot_data$shape_group <- ifelse(all_plot_data$ages == "total", 19, 5)

# ICD POSITIONING
icd_pos <- aggregate(id ~ grouping, data = all_plot_data, mean)

# PLOT
bgplot <- ggplot(all_plot_data,
                 aes(x=id, group = ages, col = ages)) +
  theme_classic() +
  scale_x_reverse() +
  scale_x_reverse(name = "",
                  breaks = unique(all_plot_data$id),
                  labels = unique(all_plot_data$names), 
                  sec.axis = sec_axis(trans = ~., name = "", breaks = icd_pos$id, 
                                      labels = icd_pos$grouping)) +
  geom_hline(yintercept = 1) +
  theme(axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.text.x.bottom = element_text(angle = 90, vjust = 0, hjust = 1,
                                          size = 8), plot.title = element_text(hjust = 0.5, size = 9),
        axis.text.x.top = element_text(size = 8),
        panel.grid.major.y = element_line(linetype = 3, colour = "grey")) +
  geom_pointrange(position = position_dodge(0.8), size = 0.3) +
  scale_y_continuous(limits = range(with(all_plot_data, c(rr_lowci, rr_upci)))) +
  ggtitle("Relative Risk of Hospital Admission\nEngland, 2008 - 2019") 

# CONFIDENCE INTERVALS
heatplot <- bgplot +
  aes(y = rr, ymin = rr_lowci, ymax = rr_upci, shape = factor(shape_group)) +  # Mapping shape to the new column
  scale_color_manual(guide = "none",
                     values = c("black",paletteer_c("ggthemes::Red", 4))) +
  scale_shape_manual(guide = "none",
                     values = c(19, 5)) +  # 1 for circles, 5 for diamonds
  scale_linetype_manual(values = c("solid","solid","solid","solid","dotted")) +
  ylab(sprintf("RR at 99th\ntemperature percentile")) 

# CREATE LEGEND
legplot <- ggplot(all_plot_data, aes(x = id, group = ages, col = ages, shape = factor(shape_group))) +
  theme_void() + ylim(c(0, 0)) +
  geom_pointrange(aes(y = rr, ymin = rr_lowci, ymax = rr_upci, shape = factor(shape_group))) +
  scale_color_manual(name = "Age group",
                     values = c("black",paletteer_c("ggthemes::Red", 4))) + 
  scale_shape_manual(guide = "none",
                     values = c(19,5)) +  # 1 for circles, 5 for diamonds
  theme(legend.position = "top", legend.direction = "horizontal") 

pdf("./figures/combined_forest_plot.pdf", height=6.5, width=16)

# COMBINE PLOT AND LEGEND
heatplot / legplot + 
  plot_layout(heights = c(1, .1))

dev.off()

png_plot <- heatplot / legplot + 
  plot_layout(heights = c(1, .1))

# SAVE
ggsave("./figures/fig2.png", png_plot, height=7, width=14, dpi=1200)
write.csv(all_plot_data, file="./tables/forest_plot_data.csv")

################################################################################
# SUPPLEMENTARY

# ALL AGES
for (c in setcause) {
  
  pdf(paste0("./figures/supplementary/",c,"_exposure_response.pdf"), height=9, width=9.5)
  
  layout(matrix(1:6,3, byrow=T))
  par(mar=c(4,4,2,0.5), las=1, mgp=c(2.5,1,0))
  
  for (a in c(5,1,2,3,4)) {
    plot(pooled_results_list_c[[c]][[a]], type="l", col=a, ylab="RR", ylim=c(.6,1.6), lwd=2,
         xlab="Temperature (C)", main=paste(agevarlab[a]))
  }
  
  dev.off()
}

# TOTAL
for (c in setcause) {
  
  pdf(paste("./figures/supplementary/total_age_",c,"_exposure_response.pdf"), height=9, width=9.5)
  
  layout(matrix(1:1,1, byrow=T))
  par(mar=c(4,4,2,0.5), las=1, mgp=c(2.5,1,0))
  
  plot(pooled_results_list_c[[c]][[5]], type="l", col=5, ylab="RR", ylim=c(.6,1.6), lwd=2,
       xlab="Temperature (C)")
  
  dev.off()
}