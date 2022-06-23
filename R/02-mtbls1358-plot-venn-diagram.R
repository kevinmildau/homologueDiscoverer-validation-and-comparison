################################################################################
#' Plotting comparative venn diagrams for nontarget and homologueDiscoverer on
#' mtbls1358 data.
#'
#' Author: Kevin Mildau
#' Date: 2022 June
################################################################################

# Load Packages & data ---------------------------------------------------------
library(devtools)
library(homologueDiscoverer)
prefix <- "output/"
source("R/00-utils-validation-and-comparison.R")

# Load annotated peak tables
out_hd <-
  readRDS(paste0(prefix, "mtbls1358_out_homologueDiscoverer.RDS"))
out_nt <- readRDS(paste0(prefix, "mtbls1358_out_nontarget.RDS"))
out_nt_01 <- readRDS(paste0(prefix, "mtbls1358_out_nontarget_01.RDS"))
out_nt_02 <- readRDS(paste0(prefix, "mtbls1358_out_nontarget_02.RDS"))
out_nt_03 <- readRDS(paste0(prefix, "mtbls1358_out_nontarget_03.RDS"))

# Exact Series Overlap between nontarget runs and homologueDiscoverer ----------
pdf(file = paste0(prefix, "mtbls1358_out_venn_series.pdf"), family = "Helvetica",
    title = "R Graphics Output", fonts = NULL, width = 10, height = 6)
display_venn(list("HD" = createSeriesStrings(out_hd),
                  "RTTOL 50" = createSeriesStrings(out_nt),
                  "RTTOL DEF" = createSeriesStrings(out_nt_01),
                  "RTTOL 1" = createSeriesStrings(out_nt_02),
                  "RTTOL 5" = createSeriesStrings(out_nt_03)),
             lwd = 2, lty = 'blank',
             fill = c("#FF00FF", "#999999", "#E69F00", "#56B4E9", "#009E73"),
             cex = .9, fontface = "italic",
             cat.cex = 1, cat.fontface = "bold", cat.default.pos = "outer",
             cat.dist = c(0.055, 0.1, 0.1, 0.1,0.1))
dev.off()

# Annotated peak overlap between nontarget runs and homologueDiscoverer --------
pdf(file = paste0(prefix, "mtbls1358_out_venn_peaks.pdf"), family = "Helvetica",
    title = "R Graphics Output", fonts = NULL, width = 10, height = 6)
display_venn(list(HD = createAnnotedPeakIdList(out_hd),
                  "RTTOL 50" = createAnnotedPeakIdList(out_nt),
                  "RTTOL DEF" = createAnnotedPeakIdList(out_nt_01),
                  "RTTOL 1" = createAnnotedPeakIdList(out_nt_02),
                  "RTTOL 5" = createAnnotedPeakIdList(out_nt_03)),
             lwd = 2, lty = 'blank',
             fill = c("#FF00FF", "#999999", "#E69F00", "#56B4E9", "#009E73"),
             cex = .9, fontface = "italic",
             cat.cex = 1, cat.fontface = "bold", cat.default.pos = "outer",
             cat.dist = c(0.055, 0.1, 0.1, 0.1,0.1))
dev.off()


# Plot multiple assignments for different runs ---------------------------------
out_nt_01 %>%
  mutate(., n_assignments = as.factor(n_assignments)) %>%
  select(., peak_id, n_assignments) %>%
  unique() %>%
  ggplot(., aes(x = n_assignments)) +
  geom_bar(stat = "count",) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(family="mono")) +
  xlab("number of series assignments per peak") +
  ylab("Counts (log10 scale)") +
  ggtitle("Nontarget Multiple Assignment Overview MTBLS1358 Data - Default RT Tolerance") +
  scale_y_log10() +
  geom_text(stat='count', aes(label=..count..), vjust=-1, size = 2)
ggsave(paste0(prefix, "mtbls1358-multiple-assignments-rttol-def.pdf"), device = "pdf",
       width = 30, height = 16, units = "cm", dpi = 300)

out_nt_02 %>%
  mutate(., n_assignments = as.factor(n_assignments)) %>%
  select(., peak_id, n_assignments) %>%
  unique() %>%
  ggplot(., aes(x = n_assignments)) +
  geom_bar(stat = "count",) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(family="mono")) +
  xlab("number of series assignments per peak") +
  ylab("Counts (log10 scale)") +
  ggtitle("Nontarget Multiple Assignment Overview MTBLS1358 Data - RT Tolerance of 1") +
  scale_y_log10() +
  geom_text(stat='count', aes(label=..count..), vjust=-1, size = 2)
ggsave(paste0(prefix, "mtbls1358-multiple-assignments-rttol-1.pdf"), device = "pdf",
       width = 30, height = 16, units = "cm", dpi = 300)

out_nt_03 %>%
  mutate(., n_assignments = as.factor(n_assignments)) %>%
  select(., peak_id, n_assignments) %>%
  unique() %>%
  ggplot(., aes(x = n_assignments)) +
  scale_y_log10() +
  geom_bar(stat = "count") +
  geom_text(stat='count', aes(label=..count..), vjust=-1, size = 2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(family="mono")) +
  xlab("number of series assignments per peak") +
  ylab("Counts (log10 scale)") +
  ggtitle("Nontarget Multiple Assignment Overview MTBLS1358 Data - RT Tolerance of 5")
ggsave(paste0(prefix, "mtbls1358-multiple-assignments-rttol-5.pdf"), device = "pdf",
       width = 30, height = 16, units = "cm", dpi = 300)
