################################################################################
# Script: Jaccard Similarity Analysis and Visualization Across DEA Methods
#
# Description:
# This script parses and visualizes intra-method Jaccard similarity indices
# generated from the DE-Signal-10 simulated datasets.
#
# Input:
#   - CSV file containing Jaccard similarity summaries.
#
# Output:
#   - Faceted ggplot figure displaying average Jaccard similarity indices across datasets and DEA methods.
#
# Author: Erda Qorri
# Date: 16-05-2026


# Load libraries
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

##### Read Jaccard summary #####
jac_dist <- read.csv("path_to_jaccard_summary.csv")
jac_dist <- read.csv("Downloads/in_method_jaccard_summary.csv")


##### Parse dataset, filtering strategy, and method #####
jac_dist_parsed <- jac_dist %>%
  extract(
    sample,
    into = c("Dataset", "Filtering", "Method"),
    regex = "^([^_]+)_([^_]+)_DE0\\.1_(.*)_DE_results_all_reps$",
    remove = FALSE
  ) %>%
  select(Dataset, Filtering, Method, avg_jaccard) %>%
  mutate(
    Method = recode(
      Method,
      "DESeq2_BetaPrior" = "DESeq2-BP",
      "DESeq2_LRT" = "DESeq2-LRT",
      "DESeq2_WaldTest" = "DESeq2-Wald",
      "Limma-voom_default_TMMwsp" = "LV-Def-TMMwsp",
      "Limma-voom_DFRBT_TMM" = "LV-DefRBT-TMM",
      "Limma-voom_DFRBT_TMMwsp" = "LV-DefRBT-TMMwsp",
      "Limma-voom_LmFit_TMM" = "VoomLmFit-TMM",
      "Limma-voom_LmFit_TMMwsp" = "VoomLmFit-TMMwsp",
      "Limma-voom_default_TMM" = "LV-Def-TMM",
      "Limma-voom_DFRBT_quantile" = "LV-DefRBT-QT"
    ),
    Filtering = recode(
      Filtering,
      "autofilter" = "Autofilter",
      "min5" = "Min 5",
      "min1" = "Min 1"
    ),
    Filtering = factor(Filtering, levels = c("Autofilter", "Min 5", "Min 1"))
  )

##### Plot avaerage Jaccard similarity #####
ggplot(jac_dist_parsed, aes(x = Method, y = avg_jaccard, color = Filtering)) +
  geom_segment(
    aes(x = Method, xend = Method, y = 0, yend = avg_jaccard),
    position = position_dodge(0.5),
    show.legend = FALSE
  ) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  facet_wrap(~ Dataset) +
  labs(
    title = "Jaccard Similarity Index Across Methods",
    x = "",
    y = "Average Jaccard Index"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
# facet_wrap(~ dataset) +
theme_classic(base_size = 17.5) +
  # facet_wrap(~ dataset, ncol = 2, scales = "free_y") +
  theme(
    strip.background = element_rect(fill = "#cfe8ff", color = "#cfe8ff", linewidth = 1.2),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "top",
    legend.text = element_text(size = 18),
    legend.spacing.y = unit(2, "cm")
  ) 




