################################################################################
# UK-TRH: SMALL-AREA ANALYSIS OF TEMPERATURE RELATED HOSPITALISATIONS IN ENGLAND
################################################################################

################################################################################
# Forest plots of second stage results
################################################################################

# Temperature reference points:
heat_temp <- round(quantile(datatmean$tmean, 99/100, na.rm = T),1)

names <- unlist(lapply(seq(causelist), function(x) causelist[[x]][["name"]]))
shortnames <- unlist(lapply(seq(causelist), function(x) causelist[[x]][["shortname"]]))
cause_ref_list <- data.frame(names, shortnames)

rr <- rr_lowci <- rr_upci <- cause <- all_plot_data <- NULL
ages <- c("total","age1864","age6574","age7584","age85plus")
#ages <- c("total","age1834","age3564","age6574","age7584","age85plus")
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
#all_plot_data$grouping[all_plot_data$cause %in% c("mental","dementia")] <- "Mental and behavioural" 
all_plot_data$grouping[all_plot_data$cause %in% c("endo","diabetes", "metabolic")] <- "Endocrine, nutritional, metabolic"

# Order by ICD group, diagnosis and age
all_plot_data <-all_plot_data[with(all_plot_data, 
  order(grouping, names, cause, ages)),]
# Order columns left to right:
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


# Background lines:
unid <- unique(all_plot_data$id)
bglines <- data.frame(pos = unid[-1][diff(unid)==1] - .5)

# Create a new column to specify shapes based on age group
all_plot_data$shape_group <- ifelse(all_plot_data$ages == "total", 19, 5)  # 5 is diamond, 1 is circle

# ICD grouping label position:
icd_pos <- aggregate(id ~ grouping, data = all_plot_data, mean)

# Create background plot:
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
  #coord_flip()

# Add confidence intervals as vertical bars for the 'total' group
heatplot <- bgplot +
  aes(y = rr, ymin = rr_lowci, ymax = rr_upci, shape = factor(shape_group)) +  # Mapping shape to the new column
  scale_color_manual(guide = "none",
    values = c("black",paletteer_c("ggthemes::Red", 4))) +
  scale_shape_manual(guide = "none",
    values = c(19, 5)) +  # 1 for circles, 5 for diamonds
  scale_linetype_manual(values = c("solid","solid","solid","solid","dotted")) +
  ylab(sprintf("RR at 99th\ntemperature percentile")) 
 # coord_flip() 
  
# # Add a dashed line for the 'total' group
# heatplot <- heatplot + 
#   geom_line(data = subset(all_plot_data, ages == "total"),
#             aes(x = rr, y = id), linetype = "dashed", color = "black", size = 1.2)

# Create legend plot
legplot <- ggplot(all_plot_data, aes(x = id, group = ages, col = ages, shape = factor(shape_group))) +
  theme_void() + ylim(c(0, 0)) +
  geom_pointrange(aes(y = rr, ymin = rr_lowci, ymax = rr_upci, shape = factor(shape_group))) +
  scale_color_manual(name = "Age group",
    values = c("black",paletteer_c("ggthemes::Red", 4))) + 
  scale_shape_manual(guide = "none",
    values = c(19,5)) +  # 1 for circles, 5 for diamonds
  theme(legend.position = "top", legend.direction = "horizontal") 

pdf("./figures/combined_forest_plot.pdf", height=6.5, width=16)

# Combine heatplot and legplot using patchwork
heatplot / legplot + 
  plot_layout(heights = c(1, .1))

dev.off()

png_plot <- heatplot / legplot + 
  plot_layout(heights = c(1, .1))

ggsave("./figures/fig2.png", png_plot, height=8, width=16, dpi=1200)
ggsave("./figures/fig2.png", png_plot, height=7, width=14, dpi=1200)

write.csv(all_plot_data, file="./tables/forest_plot_data.csv")