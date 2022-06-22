################################################################################
#' Plotting homologueDiscoverer & nontarget output on mtbls1358 data.
#'
#' Author: Kevin Mildau, 2022 June
#'
#' Description:
################################################################################

# Load Packages & data #########################################################
library(devtools)
library(homologueDiscoverer)
prefix <- "output/"
source("R/00-utils-validation-and-comparison.R")
# Load annotated peak tables
mtbls1358_out_hd <-
  readRDS(paste0(prefix, "mtbls1358_out_homologueDiscoverer.RDS"))
mtbls1358_out_nt <-
  readRDS(paste0(prefix, "mtbls1358_out_nontarget.RDS"))

# Plot Annotated Peak Table ####################################################

p <- plotAnnotatedStatic(mtbls1358_out_hd) +
  ggtitle("homologueDiscoverer Annotated Peak Table - MTBLS1358 Data")
ggsave(plot = p,
       filename = paste0(prefix, "mtbls1358_out_homologueDiscoverer_aptb.pdf"),
       device = "pdf", width = 20, height = 12, units = "cm", dpi = 300)

p <- plotAnnotatedStatic(mtbls1358_out_nt) +
  ggtitle("nontarget Annotated Peak Table - MTBLS1358 Data")
ggsave(plot = p, filename = paste0(prefix, "mtbls1358_out_nontarget_aptb.pdf"),
       device = "pdf", width = 20, height = 12, units = "cm", dpi = 300)

# Multiple assignment plotting #################################################
mtbls1358_out_nt %>%
  mutate(., n_assignments = as.factor(n_assignments)) %>%
  select(., peak_id, n_assignments) %>%
  unique() %>%
  ggplot(., aes(x = n_assignments)) +
  geom_bar(stat = "count",) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(family="mono")) +
  xlab("number of series assignments per peak") +
  xlab("number of series assignments per peak") +
  ylab("Counts (log10 scale)") +
  ggtitle("Nontarget Multiple Assignment Overview MTBLS1358 Data") +
  scale_y_log10() +
  geom_text(stat='count', aes(label=..count..), vjust=-1, size = 2)
ggsave(paste0(prefix, "mtbls1358-multiple-assignments.pdf"), device = "pdf",
       width = 30, height = 16, units = "cm", dpi = 300)

# Interactive Plotting #########################################################
# homologueDiscoverer
ptbPlotInteractive(mtbls1358_out_hd)
sdbPlotInteractive(sdbCreate(mtbls1358_out_hd))

# Nontarget
plotNontarget(mtbls1358_out_nt)
ptbPlotInteractive(mtbls1358_out_nt)








