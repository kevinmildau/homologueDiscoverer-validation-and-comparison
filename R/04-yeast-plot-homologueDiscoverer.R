################################################################################
#' Visualizing output of homologueDiscoverer on yeast data.
#'
#' Author: Kevin Mildau
#' Date: 2022 June
################################################################################
library(homologueDiscoverer)
prefix <- "output/"
source("R/00-utils-validation-and-comparison.R")

yeast_01 <- readRDS(file = paste0(prefix, "yeast_out_homologueDiscoverer_decrement.RDS"))
yeast_02 <- readRDS(file = paste0(prefix, "yeast_out_homologueDiscoverer_increment.RDS"))
yeast <- readRDS(file = paste0(prefix, "yeast_hd.RDS"))

# Combine runs into single series_db -------------------------------------------
sdbrs <- sdbCreate(yeast_01)
sdbrs <- sdbPush(sdbrs, yeast_02)

# Plot Static ------------------------------------------------------------------
# Run annotation run on raw yeast data to generate annotated peak table
yeast_anno <- annotateHomologues(yeast, sdbrs,
                                 rt_tolerance = 5,
                                 ppm_matching_tolerance = 5,
                                 ppm_step_tolerance =  5,
                                 min_match_length = 4)
p <- plotAnnotatedStatic(yeast_anno) +
  ggtitle("homologueDiscoverer Annotated Peak Table - Yeast Data")
ggsave(plot = p, filename = paste0(prefix,
                                   "yeast_out_homologueDiscoverer_aptb.pdf"),
       device = "pdf", width = 20, height = 12, units = "cm", dpi = 300)

# Interactive Plotting ---------------------------------------------------------
sdbPlotInteractive(sdbrs)
ptbPlotInteractive(yeast_anno)
