################################################################################
#' Visualizing output of nontarget on yeast data.
#'
#' Author: Kevin Mildau
#' Date: 2022 June
################################################################################
library(homologueDiscoverer)
prefix <- "output/"
source("R/00-utils-validation-and-comparison.R")

yeast_nt <- readRDS(file = paste0(prefix, "yeast_out_nontarget_decrement.RDS"))

# Multiple Assignment Histogram ------------------------------------------------
yeast_nt %>%
  mutate(., n_assignments = as.factor(n_assignments)) %>%
  select(., peak_id, n_assignments) %>%
  unique() %>%
  ggplot(., aes(x = n_assignments)) +
  geom_bar(stat = "count",) + theme_minimal() +
  geom_text(stat='count', aes(label=..count..), vjust=-1, size = 2) +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(family="mono")) +
  xlab("number of series assignments per peak") +
  xlab("number of series assignments per peak") +
  ylab("Counts (log10 scale)") +
  ggtitle("Nontarget Multiple Assignment Overview Yeast Data") +
  scale_y_log10()
ggsave(paste0(prefix, "yeast-multiple-assignments.pdf"), device = "pdf",
       width = 20, height = 12, units = "cm", dpi = 300)

# Static annotated peak table plot ---------------------------------------------
p <- plotAnnotatedStatic(yeast_nt) +
  ggtitle("Nontarget Annotated Peak Table - Decrementing Series in Yeast Data")
ggsave(plot = p, filename = paste0(prefix, "yeast_out_nontarget_aptb.pdf"),
       device = "pdf", width = 20, height = 12, units = "cm", dpi = 300)

# Interactive Plots ------------------------------------------------------------
# Highlight the number of multiple series assignments per peak
plotNontarget(yeast_nt)
# Interactive annotated peak table plot
ptbPlotInteractive(yeast_nt)
sdbPlotInteractive(sdbCreate(yeast_nt))
