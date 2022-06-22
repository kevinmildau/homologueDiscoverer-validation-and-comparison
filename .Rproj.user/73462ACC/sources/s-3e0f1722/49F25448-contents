################################################################################
#' Running homologue identifiers packages on yeast data to find decrementing
#' homologue series.
#'
#' Author: Kevin Mildau, 2022 June
#'
#' Description:
#' --> Run nontarget::homol search on yeast data.
#' --> Run homologueDiscoverer on yeast data.
################################################################################

# Load Packages & data #########################################################
library(devtools)
library(nontarget)
library(enviPat)
library(homologueDiscoverer)
prefix <- "output/"
source("R/00-utils-validation-and-comparison.R")
data("isotopes") # required for nontarget, from enviPat

# Load peak tables
yeast_hd <- readRDS(paste0(prefix, "yeast_hd.RDS"))
yeast_nt <- readRDS(paste0(prefix, "yeast_nt.RDS"))
yeast_index <- readRDS(paste0(prefix, "yeast_index.RDS"))

# Run nontarget for decrementing series ########################################
outnt01 <- nontarget::homol.search(peaklist = yeast_nt, isotopes, FALSE,
                                   minmz = 10, maxmz = 50, mztol = 5,
                                   minrt = -100, maxrt = 0.1, # for decrement
                                   ppm = TRUE,
                                   minlength = 6,
                                   rttol = 50,
                                   R2 = 0.98) # default value
outnt01_ptb <- createPseudoPeakTableFromNonTarget(outnt01)
# Add sample peak_ids to nontarget output
outnt01_ptb <- left_join(outnt01_ptb, yeast_index, by = c("mz", "rt"))

saveRDS(object = outnt01_ptb,
        file = paste0(prefix, "yeast_out_nontarget_decrement.RDS"))

# Run homologueDiscoverer for decrementing and incrementing series each ########
outrs01 <- detectHomologues(yeast_hd,
                            mz_min = 10, mz_max = 50,
                            rt_min = 0.1, rt_max = 100,
                            min_series_length = 6,
                            ppm_tolerance = 5,
                            search_mode = "untargeted",
                            step_mode = "decrement")
outrs02 <- detectHomologues(yeast_hd,
                            mz_min = 10, mz_max = 50,
                            rt_min = 0.1, rt_max = 100,
                            min_series_length = 6,
                            ppm_tolerance = 5,
                            search_mode = "untargeted",
                            step_mode = "increment")
saveRDS(object = outrs01,
        file = paste0(prefix, "yeast_out_homologueDiscoverer_decrement.RDS"))
saveRDS(object = outrs02,
        file = paste0(prefix, "yeast_out_homologueDiscoverer_increment.RDS"))

