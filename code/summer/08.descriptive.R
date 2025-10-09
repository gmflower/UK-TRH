################################################################################
# UK-TRH: SMALL-AREA ANALYSIS OF TEMPERATURE RELATED HOSPITALISATIONS IN ENGLAND
################################################################################

# Descriptive analysis of hospital admissions dataset

names <- unlist(lapply(seq(causelist), function(x) causelist[[x]][["name"]]))
shortnames <- unlist(lapply(seq(causelist), function(x) causelist[[x]][["shortname"]]))
cause_ref_list <- data.frame(names, shortnames)

cause_ref_list$grouping[cause_ref_list$shortnames %in% c("resp","ari","pneumonia","copd","asthma")] <- 5 
cause_ref_list$grouping[cause_ref_list$shortnames %in% c("cvd","mi","stroke","hf", "hypo")] <- 4 
cause_ref_list$grouping[cause_ref_list$shortnames %in% c("genito","renal","arf","ckd")] <- 3
cause_ref_list$grouping[cause_ref_list$shortnames %in% c("infect","bacterial")] <- 2
#cause_ref_list$grouping[cause_ref_list$shortnames %in% c("mental","dementia")] <- 2 
cause_ref_list$grouping[cause_ref_list$shortnames %in% c("endo","diabetes", "metabolic")] <- 1 

cause_ref_list$names[names=="Respiratory"] <- "Respiratory (All)"
cause_ref_list$names[names=="Cardiovascular"] <- "Cardiovascular (All)"
cause_ref_list$names[names=="Genitourinary"] <- "Genitourinary (All)"
cause_ref_list$names[names=="Infectious and parasitic"] <- "Infectious and parasitic (All)"
#cause_ref_list$names[names=="Mental and behavioural"] <- "Mental and behavioural (All)"
cause_ref_list$names[names=="Endocrine, nutritional, metabolic"] <- "Endocrine, nutritional, metabolic (All)"

hes_aggregated <- hesdata %>% 
  filter(agegr=="total" & cause %in% setcause) %>% 
  group_by(cause) %>% 
  summarise(count = sum(count)) %>% 
  left_join(cause_ref_list, by = join_by(cause==shortnames)) %>% 
  as.data.frame()
hes_aggregated <- arrange(hes_aggregated,count)
#hes_aggregated <- tail(hes_aggregated,20)
leadingcauses <- unique(hes_aggregated$cause)

hes_byage <- hesdata %>% 
  filter(agegr!="total" & cause %in% setcause) %>% 
  group_by(cause, agegr) %>% 
  summarise(count = sum(count)) %>% 
  left_join(cause_ref_list, by = join_by(cause==shortnames)) %>% 
  as.data.table()

# Keep only leading causes:
hes_byage <- hes_byage[cause %in% leadingcauses,]

#hes_byage <- hes_byage %>% 
#  mutate(names = factor(names, levels = names))

pdf("./figures/total_admission_burden.pdf", height=9, width=9.5)

fig1 <- ggplot(data=hes_byage, aes(fill=agegr, x=count/1000, y=reorder(reorder(names,count),grouping))) +
  # First geom_bar: original bars (no outlines)
  geom_bar(position="stack", stat="identity", color=NA) +  
  # Second geom_bar: bars with black outlines
  geom_bar(data = hes_byage %>% filter(names %in% c("Respiratory (All)", "Cardiovascular (All)", 
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
  
write.csv(hes_byage, file="./tables/hes_byage.csv")

ggsave("./figures/fig1.png", fig1, height=9, width=9.5, dpi=1200)

  