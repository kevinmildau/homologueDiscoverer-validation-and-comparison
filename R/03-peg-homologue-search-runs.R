################################################################################
#' Running homologue identifiers packages on peg data for validation and
#' output comparison.
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
peg17_joined <- readRDS(paste0(prefix, "peg17_joined.RDS"))
peg17_joined_hd <- readRDS(paste0(prefix, "peg17_joined_hd.RDS"))
peg17_joined_nt <- readRDS(paste0(prefix, "peg17_joined_nt.RDS"))

peg70_joined <- readRDS(paste0(prefix, "peg70_joined.RDS"))
peg70_joined_hd <- readRDS(paste0(prefix, "peg70_joined_hd.RDS"))
peg70_joined_nt <- readRDS(paste0(prefix, "peg70_joined_nt.RDS"))

# Homologue Identifier Runs ####################################################
# rsidentifier
out_rs_peg17 <- rsidentifier::detectHomologues(peg17_joined_hd,
                                                 mz_min = 5, mz_max = 100,
                                                 rt_min = 0.1, rt_max = 100,
                                                 min_series_length = 4,
                                                 ppm_tolerance = 20,
                                                 search_mode = "untargeted",
                                                 step_mode = "increment")

out_rs_peg17 <-
  left_join(out_rs_peg17,
            select(peg17_joined, c("peak_id", "peg", "mix", "plasma")),
            by = c("peak_id"))

out_rs_peg70 <- rsidentifier::detectHomologues(peg70_joined_hd,
                                                 mz_min = 5, mz_max = 100,
                                                 rt_min = 0.1, rt_max = 100,
                                                 min_series_length = 4,
                                                 ppm_tolerance = 20,
                                                 search_mode = "untargeted",
                                                 step_mode = "increment")
out_rs_peg70 <-
  left_join(out_rs_peg70,
            select(peg70_joined, c("peak_id", "peg", "mix", "plasma")),
            by = c("peak_id"))

# nontarget
out_nt_peg17 <- nontarget::homol.search(peaklist = peg17_joined_nt,
                                        isotopes, FALSE, minmz = 5, maxmz = 100,
                                        mztol = 20, minrt = 0.1, maxrt = 100,
                                        ppm = TRUE, minlength = 4, rttol = 50,
                                        R2 = 0.98)
out_nt_ptb_peg17 <- createPseudoPeakTableFromNonTarget(out_nt_peg17)
out_nt_ptb_peg17 <-
  left_join(out_nt_ptb_peg17,
            select(peg17_joined,
                   c("mz", "rt", "peak_id", "peg", "mix", "plasma")),
            by = c("mz", "rt"))
out_nt_peg70 <- nontarget::homol.search(peaklist = peg70_joined_nt,
                                        isotopes, FALSE, minmz = 5, maxmz = 100,
                                        mztol = 20, minrt = 0.1, maxrt = 100,
                                        ppm = TRUE, minlength = 4, rttol = 50,
                                        R2 = 0.98)
out_nt_ptb_peg70 <- createPseudoPeakTableFromNonTarget(out_nt_peg70)
out_nt_ptb_peg70 <-
  left_join(out_nt_ptb_peg70,
            select(peg70_joined,
                   c("mz", "rt", "peak_id", "peg", "mix", "plasma")),
            by = c("mz", "rt"))

# save run output ##############################################################
saveRDS(object = out_rs_peg17 , file = paste0(prefix, "peg17_out_hd.RDS"))
saveRDS(object = out_rs_peg70 , file = paste0(prefix, "peg70_out_hd.RDS"))
saveRDS(object = out_nt_ptb_peg17 , file = paste0(prefix, "peg17_out_nt.RDS"))
saveRDS(object = out_nt_ptb_peg70 , file = paste0(prefix, "peg70_out_nt.RDS"))
