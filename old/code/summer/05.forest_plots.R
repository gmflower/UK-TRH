################################################################################
# UK-TRH: SMALL-AREA ANALYSIS OF TEMPERATURE RELATED HOSPITALISATIONS IN ENGLAND
################################################################################

#### Forest plots of second stage results

# Temperature reference points:
ref_temp <- round(quantile(datatmean$tmean, 50/100, na.rm = T),1)
heat_temp <- round(quantile(datatmean$tmean, 99/100, na.rm = T),1)

# Define groups of diagnoses:
cause_groups <- list(
  list(group_name="Cardiovascular", causes=c("cvd","mi","pulmheart","stroke","ihd","hf", "hypo")),
  list(group_name="Respiratory", causes=c("resp","ari","pneumonia","copd","asthma")),
  list(group_name="Genitourinary", causes=c("genito","renal","arf")),
  list(group_name="Endocrine, nutritional, metabolic", causes=c("endo","diabetes")),
  list(group_name="Mental and behavioural", causes=c("mental","dementia")),
  list(group_name="Infectious and parasitic", causes=c("infect"))
)
names(cause_groups) <- c("Cardiovascular","Respiratory","Genitourinary",
  "Endocrine, nutritional, metabolic", "Mental and behavioural", 
  "Infectious and parasitic")
  

#names(cause_groups)[1]
#causes <- cause_groups[[1]][["causes"]]

Temperature <- c(rep(heat_temp,5))
plots <- NULL

for (c in (c(cause_groups[["Respiratory"]][["causes"]],
            cause_groups[["Cardiovascular"]][["causes"]],
            cause_groups[["Genitourinary"]][["causes"]],
            cause_groups[["Endocrine, nutritional, metabolic"]][["causes"]],
            cause_groups[["Mental and behavioural"]][["causes"]],
            cause_groups[["Infectious and parasitic"]][["causes"]]))) {

  print(c)
  RR  <- c(pooled_results_list_c[[c]][["total"]][["allRRfit"]][[paste0(heat_temp)]],
    pooled_results_list_c[[c]][["age064"]][["allRRfit"]][[paste0(heat_temp)]],
    pooled_results_list_c[[c]][["age6574"]][["allRRfit"]][[paste0(heat_temp)]],
    pooled_results_list_c[[c]][["age7584"]][["allRRfit"]][[paste0(heat_temp)]],
    pooled_results_list_c[[c]][["age85plus"]][["allRRfit"]][[paste0(heat_temp)]]) 
  
  RR_lowCI <- c(pooled_results_list_c[[c]][["total"]][["allRRlow"]][[paste0(heat_temp)]],
    pooled_results_list_c[[c]][["age064"]][["allRRlow"]][[paste0(heat_temp)]],
    pooled_results_list_c[[c]][["age6574"]][["allRRlow"]][[paste0(heat_temp)]],
    pooled_results_list_c[[c]][["age7584"]][["allRRlow"]][[paste0(heat_temp)]],
    pooled_results_list_c[[c]][["age85plus"]][["allRRlow"]][[paste0(heat_temp)]])
  
  RR_upperCI <- c(pooled_results_list_c[[c]][["total"]][["allRRhigh"]][[paste0(heat_temp)]],
    pooled_results_list_c[[c]][["age064"]][["allRRhigh"]][[paste0(heat_temp)]],
    pooled_results_list_c[[c]][["age6574"]][["allRRhigh"]][[paste0(heat_temp)]],
    pooled_results_list_c[[c]][["age7584"]][["allRRhigh"]][[paste0(heat_temp)]],
    pooled_results_list_c[[c]][["age85plus"]][["allRRhigh"]][[paste0(heat_temp)]])
  
  Group <- c("total", "age064", "age6574", "age7584", "age85plus")
  
  plot_data <- list(data.frame(Temperature, RR, RR_lowCI, RR_upperCI, Group))
  
  plots[[c]] <- ggplot(data=as.data.frame(plot_data), aes(x=Group, y=RR, ymin=RR_lowCI, ymax=RR_upperCI, 
                                  group=Group , color=Group)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 
  xlab("") + ylab("Relative Risk") +
  ylim(0.55,1.6) +
  ggtitle(paste0(c)) +
  scale_colour_manual(values = c(paletteer_c("ggthemes::Blue", 5))) +
  theme_bw() + # use a white background
  theme(plot.title = element_text(hjust = 0.5, size=10),
    axis.text.x=element_blank(), legend.box = "horizontal",
    plot.margin = unit(c(0.4,0.2,0,0.2), 'lines'),
    legend.margin=margin(t = 0, unit='cm')) +
  guides(color = guide_legend(nrow = 1))
}


# Respiratory plot:
combined_respiratory <- annotate_figure(ggarrange(plots[[1]],plots[[2]],plots[[3]],
  plots[[4]],plots[[5]],
  ncol=7, nrow=1, legend=FALSE), fig.lab = "Respiratory", fig.lab.pos = c("top.left"), 
  fig.lab.face = "bold", fig.lab.size = 10)
# Cardiovascular plot:
combined_cardiovascular <- annotate_figure(ggarrange(plots[[6]],
  plots[[7]],plots[[8]],plots[[9]],plots[[10]],plots[[11]],plots[[12]],
  ncol=7, nrow=1, legend=FALSE), fig.lab = "Cardiovascular", fig.lab.pos = c("top.left"), 
  fig.lab.face = "bold", fig.lab.size = 10)
# Genitourinary plot:
combined_genitourinary <- annotate_figure(ggarrange(plots[[13]],plots[[14]],
  plots[[15]],
  ncol=7, nrow=1, legend=FALSE), fig.lab = "Genitourinary", fig.lab.pos = c("top.left"),
  fig.lab.face = "bold", fig.lab.size = 10)
# Endocrine plot:
combined_endocrine <- annotate_figure(ggarrange(plots[[16]],plots[[17]],
  ncol=7, nrow=1, legend=FALSE), fig.lab = "Endocrine, nutritional, metabolic", 
  fig.lab.pos = c("top.left"), fig.lab.face = "bold", fig.lab.size = 10)
# Mental & behavioral plot:
combined_mental <- annotate_figure(ggarrange(plots[[18]],plots[[19]],
  ncol=7, nrow=1, legend=FALSE), fig.lab = "Mental and behavioural", 
  fig.lab.pos = c("top.left"), fig.lab.face = "bold", fig.lab.size = 10)
combined_infect <- annotate_figure(ggarrange(plots[[20]],
  ncol=7, nrow=1, legend=FALSE), fig.lab = "Infectious and parasitic", 
  fig.lab.pos = c("top.left"), fig.lab.face = "bold", fig.lab.size = 10)
test_legend<-get_legend(plots[[1]])
as_ggplot(test_legend)

pdf(paste0("figures/forest_plot.pdf"), height=10, width=12)
  
ggarrange(combined_respiratory, combined_cardiovascular, combined_genitourinary, 
  combined_endocrine, combined_mental, combined_infect, ncol=1, nrow=6, common.legend = TRUE, 
  legend="bottom",legend.grob=get_legend(plots[[1]]))

dev.off()