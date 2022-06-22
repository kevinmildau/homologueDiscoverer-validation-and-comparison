################################################################################
#' Plotting homologueDiscoverer & nontarget output on peg data.
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
peg17_hd <- readRDS(paste0(prefix, "peg17_out_hd.RDS"))
peg70_hd <- readRDS(paste0(prefix, "peg70_out_hd.RDS"))
peg17_nt <- readRDS(paste0(prefix, "peg17_out_nt.RDS"))
peg70_nt <- readRDS(paste0(prefix, "peg70_out_nt.RDS"))


# Plot Annotated Peak Table ####################################################

p <- plotAnnotatedStatic(peg17_hd, max_legend = 0) +
  ggtitle("homologueDiscoverer Annotated Peak Table - PEG17 Data")
ggsave(plot = p, filename = paste0(prefix, "peg17_out_homologueDiscoverer_aptb.pdf"),
       device = "pdf", width = 20, height = 12, units = "cm", dpi = 300)

p <- plotAnnotatedStatic(peg70_hd, max_legend = 0) +
  ggtitle("homologueDiscoverer Annotated Peak Table - PEG70 Data")
ggsave(plot = p, filename = paste0(prefix, "peg70_out_homologueDiscoverer_aptb.pdf"),
       device = "pdf", width = 20, height = 12, units = "cm", dpi = 300)

p <- plotAnnotatedStatic(peg17_nt, max_legend = 0) +
  ggtitle("Nontarget Annotated Peak Table - PEG17 Data")
ggsave(plot = p, filename = paste0(prefix, "peg17_out_nontarget_aptb.pdf"),
       device = "pdf", width = 20, height = 12, units = "cm", dpi = 300)

p <- plotAnnotatedStatic(peg70_nt, max_legend = 0) +
  ggtitle("Nontarget Annotated Peak Table - PEG70 Data")
ggsave(plot = p, filename = paste0(prefix, "peg70_out_nontarget_aptb.pdf"),
       device = "pdf", width = 20, height = 12, units = "cm", dpi = 300)

# Venn Plotting ################################################################
ids_1 <- peg17_hd %>% filter(., !is.na(homologue_id) & plasma ) %>% pull(., peak_id) # plasma homologues
ids_2 <- peg17_hd %>% filter(., !is.na(homologue_id) & peg ) %>% pull(., peak_id)    # peg homologues
ids_3 <- peg17_hd %>% filter(., !is.na(homologue_id) & mix ) %>% pull(., peak_id)    # mix homologues
ids_4 <- peg17_hd %>% filter(., is.na(homologue_id) & peg ) %>% pull(., peak_id) # missed peg homologues
x <- list(
  "Plasma Peaks Detected" = ids_1,
  "PEG Peaks Detected" = ids_2,
  "PEG Peaks Missed" = ids_4
)
pdf(file = paste0(prefix, "peg17_out_validation_hd.pdf"), family = "Helvetica",
    title = "R Graphics Output", fonts = NULL, width = 5, height = 3)
display_venn(x,              lwd = 2, lty = 'blank',
             fill = c( "#E69F00", "#56B4E9", "#009E73"),
             cex = .9, fontface = "italic",
             cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "text",
             cat.dist = c(-0.05, 0.1, 0.1))
dev.off()

ids_1 <- peg70_hd %>% filter(., !is.na(homologue_id) & plasma ) %>% pull(., peak_id) # plasma homologues
ids_2 <- peg70_hd %>% filter(., !is.na(homologue_id) & peg ) %>% pull(., peak_id)    # peg homologues
ids_3 <- peg70_hd %>% filter(., !is.na(homologue_id) & mix ) %>% pull(., peak_id)    # mix homologues
ids_4 <- peg70_hd %>% filter(., is.na(homologue_id) & peg ) %>% pull(., peak_id) # missed peg homologues
x <- list(
  "Plasma Peaks Detected" = ids_1,
  "PEG Peaks Detected" = ids_2,
  "PEG Peaks Missed" = ids_4
)
pdf(file = paste0(prefix, "peg70_out_validation_hd.pdf"), family = "Helvetica",
    title = "R Graphics Output", fonts = NULL, width = 5, height = 3)
display_venn(x,              lwd = 2, lty = 'blank',
             fill = c( "#E69F00", "#56B4E9", "#009E73"),
             cex = .9, fontface = "italic",
             cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "text",
             cat.dist = c(-0.05, 0.1, 0.1))
dev.off()


ids_1 <- peg17_nt %>% filter(., !is.na(homologue_id) & plasma ) %>% pull(., peak_id) # plasma homologues
ids_2 <- peg17_nt %>% filter(., !is.na(homologue_id) & peg ) %>% pull(., peak_id)    # peg homologues
ids_3 <- peg17_nt %>% filter(., !is.na(homologue_id) & mix ) %>% pull(., peak_id)    # mix homologues
ids_4 <- peg17_nt %>% filter(., is.na(homologue_id) & peg ) %>% pull(., peak_id) # missed peg homologues
x <- list(
  "Plasma Peaks Detected" = ids_1,
  "PEG Peaks Detected" = ids_2,
  "PEG Peaks Missed" = ids_4
)
pdf(file = paste0(prefix, "peg17_out_validation_nt.pdf"), family = "Helvetica",
    title = "R Graphics Output", fonts = NULL, width = 5, height = 3)
display_venn(x,              lwd = 2, lty = 'blank',
             fill = c( "#E69F00", "#56B4E9", "#009E73"),
             cex = .9, fontface = "italic",
             cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "text",
             cat.dist = c(-0.05, 0.1, 0.1))
dev.off()

ids_1 <- peg70_nt %>% filter(., !is.na(homologue_id) & plasma ) %>% pull(., peak_id) # plasma homologues
ids_2 <- peg70_nt %>% filter(., !is.na(homologue_id) & peg ) %>% pull(., peak_id)    # peg homologues
ids_3 <- peg70_nt %>% filter(., !is.na(homologue_id) & mix ) %>% pull(., peak_id)    # mix homologues
ids_4 <- peg70_nt %>% filter(., is.na(homologue_id) & peg ) %>% pull(., peak_id) # missed peg homologues
x <- list(
  "Plasma Peaks Detected" = ids_1,
  "PEG Peaks Detected" = ids_2,
  "PEG Peaks Missed" = ids_4
)
pdf(file = paste0(prefix, "peg70_out_validation_nt.pdf"), family = "Helvetica",
    title = "R Graphics Output", fonts = NULL, width = 5, height = 3)
display_venn(x,              lwd = 2, lty = 'blank',
             fill = c( "#E69F00", "#56B4E9", "#009E73"),
             cex = .9, fontface = "italic",
             cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "text",
             cat.dist = c(-0.05, 0.1, 0.1))
dev.off()

# Comparing HD and NT ##########################################################
# Peaks
ids_1 <- peg17_hd %>% filter(., !is.na(homologue_id)) %>% pull(., peak_id)
ids_2 <- peg17_nt %>% filter(., !is.na(homologue_id)) %>% pull(., peak_id)
x <- list(
  HD = ids_1,
  NT = ids_2
)
pdf(file = paste0(prefix, "peg17_out_comparison_peaks.pdf"), family = "Helvetica",
    title = "R Graphics Output", fonts = NULL, width = 5, height = 3)
display_venn(x, fill = c( "#E69F00", "#56B4E9"), cat.default.pos = "text")
dev.off()

ids_1 <- peg70_hd %>% filter(., !is.na(homologue_id)) %>% pull(., peak_id)
ids_2 <- peg70_nt %>% filter(., !is.na(homologue_id)) %>% pull(., peak_id)
x <- list(
  HD = ids_1,
  NT = ids_2
)
pdf(file = paste0(prefix, "peg70_out_comparison_peaks.pdf"), family = "Helvetica",
    title = "R Graphics Output", fonts = NULL, width = 5, height = 3)
display_venn(x, fill = c( "#E69F00", "#56B4E9"), cat.default.pos = "text")
dev.off()

# Series
x <- list(
  HD = createSeriesStrings(peg17_hd),
  NT = createSeriesStrings(peg17_nt)
)
pdf(file = paste0(prefix, "peg17_out_comparison_series.pdf"), family = "Helvetica",
    title = "R Graphics Output", fonts = NULL, width = 5, height = 3)
display_venn(x, fill = c( "#E69F00", "#56B4E9"))
dev.off()

x <- list(
  HD = createSeriesStrings(peg70_hd),
  NT = createSeriesStrings(peg70_nt)
)
pdf(file = paste0(prefix, "peg70_out_comparison_series.pdf"), family = "Helvetica",
    title = "R Graphics Output", fonts = NULL, width = 5, height = 3)
display_venn(x, fill = c( "#E69F00", "#56B4E9"))
dev.off()

# Visualize multiple assignment of nontarget ###################################
nass <- peg17_nt %>%
  select(., peak_id, n_assignments) %>%
  unique() %>%
  pull(n_assignments)

pdf(file = paste0(prefix, "peg17_hist.pdf"), family = "Helvetica",
    fonts = NULL, width = 10, height = 6)
hist(nass,
     main = "Histogramm of number of series assignments per peak by Nontarget on PEG17 Data")
dev.off()

pdf(file = paste0(prefix, "peg17_hist_limited.pdf"), family = "Helvetica",
    fonts = NULL, width = 10, height = 6)
hist(nass, ylim = c(0,50),
     main = "Histogramm of number of series assignments per peak by Nontarget on Peg17 Data (y-axis limited)")
dev.off()

nass <- peg70_nt %>%
  select(., peak_id, n_assignments) %>%
  unique() %>%
  pull(n_assignments)

pdf(file = paste0(prefix, "peg70_hist.pdf"), family = "Helvetica",
    fonts = NULL, width = 10, height = 6)
hist(nass,
     main = "Histogramm of number of series assignments per peak by Nontarget on PEG70 Data")
dev.off()

pdf(file = paste0(prefix, "peg70_hist_limited.pdf"), family = "Helvetica",
    fonts = NULL, width = 10, height = 6)
hist(nass, ylim = c(0,50),
     main = "Histogramm of number of series assignments per peak by Nontarget on Peg70 Data (y-axis limited)")
dev.off()

peg17_nt %>%
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
  ggtitle("Nontarget Multiple Assignment Overview PEG17 Data") +
  scale_y_log10()
ggsave(paste0(prefix, "peg17-multiple-assignments.pdf"), device = "pdf",
       width = 20, height = 12, units = "cm", dpi = 300)

peg70_nt %>%
  mutate(., n_assignments = as.factor(n_assignments)) %>%
  select(., peak_id, n_assignments) %>%
  unique() %>%
  ggplot(., aes(x = n_assignments)) +
  geom_bar(stat = "count") + theme_minimal() +
  geom_text(stat='count', aes(label=..count..), vjust=-1, size = 2) +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(family="mono")) +
  xlab("number of series assignments per peak") +
  xlab("number of series assignments per peak") +
  ylab("Counts (log10 scale)") +
  ggtitle("Nontarget Multiple Assignment Overview PEG70 Data") +
  scale_y_log10()
ggsave(paste0(prefix, "peg70-multiple-assignments.pdf"), device = "pdf",
       width = 20, height = 12, units = "cm", dpi = 300)


