################################################################################
#' Running homologue identifiers packages on mtbls1358 data for output
#' comparison.
#'
#' Author: Kevin Mildau, 2022 June
#'
#' Description:
################################################################################

# Load Packages & data #########################################################
library(devtools)
library(homologueDiscoverer)
library(nontarget)
library(enviPat)
prefix <- "output/"
source("R/00-utils-validation-and-comparison.R")
data("isotopes") # required for nontarget, from enviPat

# Load peak tables
mtbls1358_hd <- readRDS(paste0(prefix, "mtbls1358_hd.RDS"))
mtbls1358_nt <- readRDS(paste0(prefix, "mtbls1358_nt.RDS"))
mtbls1358_index <- readRDS(paste0(prefix, "mtbls1358_index.RDS"))

# Run nontarget ----------------------------------------------------------------
out_nt <- nontarget::homol.search(peaklist = mtbls1358_nt, isotopes, FALSE,
                                  minmz = 13, maxmz = 15, mztol = 10,
                                  minrt = 0.1, maxrt = 100, ppm = TRUE,
                                  minlength = 4, # default value
                                  rttol = 50,
                                  R2 = 0.98) # default value
out_nt_ptb <- createPseudoPeakTableFromNonTarget(out_nt)
out_nt_ptb <- left_join(out_nt_ptb, mtbls1358_index, by = c("mz", "rt"))

# Run homologueDiscoverer ------------------------------------------------------
out_hd <- detectHomologues(mtbls1358_hd,
                           mz_min = 13, mz_max = 15,
                           rt_min = 0.1, rt_max = 100,
                           min_series_length = 4,
                           ppm_tolerance = 10,
                           search_mode = "untargeted",
                           step_mode = "increment")

# Save output ##################################################################
saveRDS(object = out_nt_ptb,
        file = paste0(prefix, "mtbls1358_out_nontarget.RDS"))
saveRDS(object = out_hd,
        file = paste0(prefix, "mtbls1358_out_homologueDiscoverer.RDS"))


# Additional nontarget runs for comparisons ####################################
outnt_rttol_1 <-
  nontarget::homol.search(peaklist = mtbls1358_nt, isotopes, FALSE,
                          minmz = 13, maxmz = 15, mztol = 10,
                          minrt = 0.1, maxrt = 100, ppm = TRUE,
                          minlength = 4, # default value
                          rttol = 0.5, # <- default value
                          R2 = 0.98) # default value
outnt_rttol_2 <-
  nontarget::homol.search(peaklist = mtbls1358_nt, isotopes, FALSE,
                          minmz = 13, maxmz = 15, mztol = 10,
                          minrt = 0.1, maxrt = 100, ppm = TRUE,
                          minlength = 4, # default value
                          rttol = 1, # <- default value
                          R2 = 0.98) # default value
outnt_rttol_3 <-
  nontarget::homol.search(peaklist = mtbls1358_nt, isotopes, FALSE,
                          minmz = 13, maxmz = 15, mztol = 10,
                          minrt = 0.1, maxrt = 100, ppm = TRUE,
                          minlength = 4, # default value
                          rttol = 5, # <- default value
                          R2 = 0.98) # default value
outnt_rttol_1 <- createPseudoPeakTableFromNonTarget(outnt_rttol_1)
outnt_rttol_1 <- left_join(outnt_rttol_1, mtbls1358_index, by = c("mz", "rt"))
outnt_rttol_2 <- createPseudoPeakTableFromNonTarget(outnt_rttol_2)
outnt_rttol_2 <- left_join(outnt_rttol_2, mtbls1358_index, by = c("mz", "rt"))
outnt_rttol_3 <- createPseudoPeakTableFromNonTarget(outnt_rttol_3)
outnt_rttol_3 <- left_join(outnt_rttol_3, mtbls1358_index, by = c("mz", "rt"))

saveRDS(object = outnt_rttol_1,
        file = paste0(prefix, "mtbls1358_out_nontarget_01.RDS"))
saveRDS(object = outnt_rttol_2,
        file = paste0(prefix, "mtbls1358_out_nontarget_02.RDS"))
saveRDS(object = outnt_rttol_3,
        file = paste0(prefix, "mtbls1358_out_nontarget_03.RDS"))
