################################################################################
#' Loading and pre-processing data for validation and comparison.
#'
#' Author: Kevin Mildau, 2022 June
#'
#' Description:
#' --> Load PEG Data
#' --> Load mtbls1358 data
#' --> Load yeast data
#' For each dataset, create a peak table as expected by homologueDiscoverer and
#' nontarget.
################################################################################
library("dplyr")
library("tibble")
library("readr")
library("homologueDiscoverer")
prefix <- "output/"
source("R/00-utils-validation-and-comparison.R")


# MTBLS Data Loading ###########################################################
# Read and convert to tibble
mtbls1358 <- tibble(read_tsv("raw-data/mtbls1358-msdial.txt", skip = 4))
# Fix column names
names(mtbls1358) <- stringr::str_replace_all(names(mtbls1358),
                                             "\\(|\\)| |/", "_")
# Rename columns & post-process peak table
mtbls1358 <- mtbls1358 %>%
  mutate(id = as.integer(Alignment_ID),
         peak_id = as.integer(Alignment_ID),
         mz = Average_Mz,
         rt = Average_Rt_min_ * 60,
         msms_assigned = MS_MS_assigned,
         msms_matched = MS_MS_matched) %>%
  rowwise() %>%
  # Generate median intensity column for peaks
  mutate(intensity =  median(c(!!! rlang::syms(grep('HEP', names(.),
                                                    value=TRUE))),
                             na.rm = TRUE)) %>%
  select(id, peak_id, mz, rt, intensity, msms_assigned, msms_matched) %>%
  ungroup()

# generate mtbls1358 input for nontarget
mtbls1358_nt <- select(mtbls1358, mz, intensity, rt) %>% as.data.frame(.)

# generate mtbls1358 input for homologueDiscoverer
mtbls1358_hd <- select(mtbls1358, peak_id, mz, rt, intensity)

# generate mtbls1358 index tibble to allow for mt and rt pair id assignment in
# nontarget output
mtbls1358_index <- select(mtbls1358, peak_id, mz,rt)


saveRDS(object = mtbls1358_nt, file = paste0(prefix, "mtbls1358_nt.RDS"))
saveRDS(object = mtbls1358_hd, file = paste0(prefix, "mtbls1358_hd.RDS"))
saveRDS(object = mtbls1358_index, file = paste0(prefix, "mtbls1358_index.RDS"))

# Load yeast data ##############################################################
# Read and convert to tibble
yeast <- tibble(read_tsv("raw-data/12c13c-yeast-msdial.txt", skip = 4))
# Fix column names
names(yeast) <- stringr::str_replace_all(names(yeast),
                                             "\\(|\\)| |/", "_")
# Rename columns & post-process peak table
yeast <- yeast %>%
  mutate(id = as.integer(Alignment_ID),
         peak_id = as.integer(Alignment_ID),
         mz = Average_Mz,
         rt = Average_Rt_min_ * 60) %>%
  rowwise() %>%
  # Generate median intensity column for peaks
  mutate(intensity =  median(c(!!! rlang::syms(grep('12C', names(.),
                                                    value=TRUE))),
                             na.rm = TRUE)) %>%
  select(id, peak_id, mz, rt, intensity) %>%
  ungroup()
# Remove all peaks below 3rd quartile of intensity over all data from table.
# This is done to reduce the size of this dataset for illustration purposes.
intensity_threshold <- quantile(yeast$intensity, probs = 0.75)
yeast <- filter(yeast, intensity >= intensity_threshold)

# generate yeast input for nontarget
yeast_nt <- select(yeast, mz, intensity, rt) %>% as.data.frame(.)

# generate yeast input for rsidentifier
yeast_hd <- select(yeast, peak_id, mz, rt, intensity)

# generate yeast index tibble to allow for mt and rt pair id assignment in
# nontarget output
yeast_index <- select(yeast, peak_id, mz,rt)

saveRDS(object = yeast_nt, file = paste0(prefix, "yeast_nt.RDS"))
saveRDS(object = yeast_hd, file = paste0(prefix, "yeast_hd.RDS"))
saveRDS(object = yeast_index, file = paste0(prefix, "yeast_index.RDS"))

# Load PEG17 and PEG70 data ####################################################

# PEG17
peg17 <- read.csv( "raw-data/all17K/tab17PEG.csv", header=T, check.names=F)
peg17[1:3, 1:3]
tmp <- do.call(rbind,
               lapply(strsplit(colnames(peg17)[2:ncol(peg17)],";"),
                      matrix,ncol=2,byrow=TRUE))
class(tmp) <- "numeric"
colnames(tmp)<- c("mz","rt")
head(tmp)
tmp <- as_tibble(tmp)
tmp$intensity <- colMeans(peg17[2:dim(peg17)[2]])
peg17 <- tmp %>% rowid_to_column(., "peak_id") %>% mutate(., rt = rt * 60)

# PEG17_plasma
peg17_plasma <- read.csv( "raw-data/all17K/tab17plasmaspikedPEG.csv",
                   header=T, check.names=F)
peg17_plasma[1:3, 1:3]
tmp <- do.call(rbind,
               lapply(strsplit(colnames(peg17_plasma)[2:ncol(peg17_plasma)],
                               ";"),
                      matrix,ncol=2,byrow=TRUE))
class(tmp) <- "numeric"
colnames(tmp)<- c("mz","rt")
tmp <- as_tibble(tmp)
tmp$intensity <- colMeans(peg17_plasma[2:dim(peg17_plasma)[2]])
peg17_plasma <- tmp %>% rowid_to_column(., "peak_id") %>%
  mutate(., rt = rt * 60)

# PEG70
peg70 <- read.csv( "raw-data/all70K/tab70PEG.csv", header=T, check.names=F)
tmp <- do.call(rbind,
               lapply(strsplit(colnames(peg70)[2:ncol(peg70)],";"),
                      matrix,ncol=2,byrow=TRUE))
class(tmp) <- "numeric"
colnames(tmp)<- c("mz","rt")
head(tmp)
tmp <- as_tibble(tmp)
tmp$intensity <- colMeans(peg70[2:dim(peg70)[2]])
peg70 <- tmp %>% rowid_to_column(., "peak_id") %>% mutate(., rt = rt * 60)

# PEG70_plasma
peg70_plasma <- read.csv( "raw-data/all70K/tab70plasmaspikedPEG.csv",
                          header=T, check.names=F)
peg70_plasma[1:3, 1:3]
tmp <- do.call(rbind,
               lapply(strsplit(colnames(peg70_plasma)[2:ncol(peg70_plasma)],
                               ";"),
                      matrix,ncol=2,byrow=TRUE))
class(tmp) <- "numeric"
colnames(tmp)<- c("mz","rt")
tmp <- as_tibble(tmp)
tmp$intensity <- colMeans(peg70_plasma[2:dim(peg70_plasma)[2]])
peg70_plasma <- tmp %>% rowid_to_column(., "peak_id") %>%
  mutate(., rt = rt * 60)


peg17_nt <- select(peg17, c("mz", "intensity", "rt")) %>% as.data.frame(.)
peg70_nt <- select(peg70, c("mz", "intensity", "rt")) %>% as.data.frame(.)
peg17_plasma_nt<- select(peg70_plasma, c("mz", "intensity", "rt")) %>%
  as.data.frame(.)
peg70_plasma_nt <- select(peg70_plasma, c("mz", "intensity", "rt")) %>%
  as.data.frame(.)

peg17_index <- select(peg17, c("mz", "peak_id", "rt")) %>% as.data.frame(.)
peg70_index <- select(peg70, c("mz", "peak_id", "rt")) %>% as.data.frame(.)
peg17_plasma_index <- select(peg70_plasma, c("mz", "peak_id", "rt")) %>%
  as.data.frame(.)
peg70_plasma_index <- select(peg70_plasma, c("mz", "peak_id", "rt")) %>%
  as.data.frame(.)

# Joining peg data for easier venn diagram creation ############################
peg17_joined <-
  full_join(select(peg17, c("peak_id", "mz", "rt")),
            select(peg17_plasma, c("peak_id", "mz", "rt", "intensity")),
            by = c("mz", "rt")) %>%
  # Create matching row identifiers between samples with only peg and
  # those with also plasma
  rowid_to_column(., var = "peak_id") %>%
  mutate(., peg = if_else(!is.na(peak_id.x), TRUE, FALSE), # peg if from peg sample
         plasma = if_else(is.na(peak_id.x) & !is.na(peak_id.y), T, F), # exclusively found in plasma
         mix = if_else(!is.na(peak_id.y), TRUE, FALSE)) # mix if from peg+plasma sample

peg70_joined <-
  full_join(select(peg70, c("peak_id", "mz", "rt")),
            select(peg70_plasma, c("peak_id", "mz", "rt", "intensity")),
            by = c("mz", "rt")) %>%
  # Create matching row identifiers between samples with only peg and
  # those with also plasma
  rowid_to_column(., var = "peak_id") %>%
  mutate(., peg = if_else(!is.na(peak_id.x), TRUE, FALSE), # peg if from peg sample
         plasma = if_else(is.na(peak_id.x) & !is.na(peak_id.y), T, F), # exclusively found in plasma
         mix = if_else(!is.na(peak_id.y), TRUE, FALSE)) # mix if from peg+plasma sample

peg17_joined_nt <- peg17_joined %>%
  filter(., mix == TRUE) %>%
  select(., c("mz", "intensity", "rt")) %>%
  as.data.frame(.)
peg17_joined_index <- peg17_joined %>% select(., c("peak_id", "mz", "rt"))

peg70_joined_nt <- peg70_joined %>%
  filter(., mix == TRUE) %>%
  select(., c("mz", "intensity", "rt")) %>%
  as.data.frame(.)
peg70_joined_index <- peg70_joined %>% select(., c("peak_id", "mz", "rt"))

peg17_joined_hd <- filter(peg17_joined, mix == T) %>%
  select(., c("peak_id", "mz", "rt", "intensity"))
peg70_joined_hd <- filter(peg70_joined, mix == T) %>%
  select(., c("peak_id", "mz", "rt", "intensity"))


saveRDS(object = peg17_joined, file = paste0(prefix, "peg17_joined.RDS"))
saveRDS(object = peg70_joined, file = paste0(prefix, "peg70_joined.RDS"))

saveRDS(object = peg17_joined_nt, file = paste0(prefix, "peg17_joined_nt.RDS"))
saveRDS(object = peg70_joined_nt, file = paste0(prefix, "peg70_joined_nt.RDS"))

saveRDS(object = peg17_joined_hd, file = paste0(prefix, "peg17_joined_hd.RDS"))
saveRDS(object = peg70_joined_hd, file = paste0(prefix, "peg70_joined_hd.RDS"))
