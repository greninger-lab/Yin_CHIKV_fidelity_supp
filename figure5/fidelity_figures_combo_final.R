library(tidyverse)
library(readxl)
library(stringr)
library(ggvenn)
library(cowplot)
library(reshape2)
library(diverse)
library(gridExtra)

# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

'%!in%' <- function(x,y)!('%in%'(x,y))
#######################################################
# Read in readcounts files generated from VarScan
#######################################################
setwd("/Volumes/lizso_backup_drive/BaseSpace_231030/231018_Einstein_PY_AViDD_CHIKV_fidelity_repeat-401900794/20231027_SURV_MPXV_BEN_and_SLU_LS_TA_KK-701001306/231130_CHIKV_fidelity_NGS/fidelity_experiment_1")
readcounts_experiment1 <- list.files(".", pattern="_interest.csv") %>% 
  map_df(~read_csv(.))
readcounts_experiment1$EXPERIMENT <- "exp1"
setwd("/Volumes/lizso_backup_drive/BaseSpace_231030/231018_Einstein_PY_AViDD_CHIKV_fidelity_repeat-401900794/20231027_SURV_MPXV_BEN_and_SLU_LS_TA_KK-701001306/231130_CHIKV_fidelity_NGS/fidelity_experiment_2")
readcounts_experiment2 <- list.files(".", pattern="_interest.csv") %>% 
  map_df(~read_csv(.))
readcounts_experiment2$EXPERIMENT <- "exp2"
readcounts_all <- rbind(readcounts_experiment1, readcounts_experiment2)

readcounts_all$SampleID <- str_replace_all(readcounts_all$SampleID, "-", "_")

readcounts_all$SampleID[readcounts_all$SampleID == "230808_Einstein_PY_025_WT_plasmid_merge_ds_309k"] <- "230808_Einstein_PY_025_plasmid_WT_1_merge_ds_309k"
readcounts_all$SampleID[readcounts_all$SampleID == "230808_Einstein_PY_026_Q192L_plasmid_merge_ds_309k"] <- "230808_Einstein_PY_026_plasmid_Q192L_1_merge_ds_309k"
readcounts_all$SampleID[readcounts_all$SampleID == "230808_Einstein_PY_027_C483Y_plasmid_merge_ds_309k"] <- "230808_Einstein_PY_027_plasmid_C483Y_1_merge_ds_309k"

readcounts_all <- readcounts_all %>% filter(SampleID %!in% c(
                                                             "230515_Einstein_PY_003_P0stock_Q192L_1_readcounts",
                                                             "231018_Einstein_PY_022_RNA_WT_2_readcounts",
                                                             "231018_Einstein_PY_023_RNA_Q192L_2_readcounts",
                                                             "231018_Einstein_PY_024_RNA_C483Y_2_readcounts"
                                                             ))
CDS_annotations <- read_csv("CDS_annotations.csv")

setwd("/Volumes/lizso_backup_drive/BaseSpace_231030/231018_Einstein_PY_AViDD_CHIKV_fidelity_repeat-401900794/20231027_SURV_MPXV_BEN_and_SLU_LS_TA_KK-701001306/231130_CHIKV_fidelity_NGS/figures/")
########################################
# Make dataframes from readcounts files
########################################

# clean up csv file -- make "None" into 0 and interpret counts as numeric instead of characters
readcounts_all <- readcounts_all %>%
  mutate_all(~ ifelse(. == "None", 0, .)) %>%
  mutate(count_A = as.numeric(count_A)) %>%
  mutate(count_C = as.numeric(count_C)) %>%
  mutate(count_G = as.numeric(count_G)) %>%
  mutate(count_T = as.numeric(count_T))

# total mismatched (sum depth of all alt bases)
readcounts_all <- readcounts_all %>% mutate(count_mismatch = count_A * (ref_base != "A") +
                                                count_C * (ref_base != "C") +
                                                count_G * (ref_base != "G") +
                                                count_T * (ref_base != "T")) %>%
  mutate(freq_mismatch = count_mismatch/depth*100) %>%
  mutate(freq_mismatch.squared = (freq_mismatch/100)^2)

# unique mismatched: count number of alternative alleles (0 to 3)
readcounts_all <- readcounts_all %>% mutate(count_unique_mismatch = replace_na((count_A/count_A) * (ref_base != "A"),0) +
                                              replace_na((count_C/count_C) * (ref_base != "C"),0) +
                                              replace_na((count_G/count_G) * (ref_base != "G"),0) +
                                              replace_na((count_T/count_T) * (ref_base != "T"),0))

# total transitions, individual SNP combos, and unique individual SNP combos
readcounts_all <- readcounts_all %>% mutate(count_transitions = count_A * (ref_base == "G") +
                                              count_C * (ref_base == "T") +
                                              count_G * (ref_base == "A") +
                                              count_T * (ref_base == "C")) %>%
                                              mutate(count_AtoC = count_C * (ref_base == "A")) %>%
                                              mutate(count_AtoG = count_G * (ref_base == "A")) %>%
                                              mutate(count_AtoT = count_T * (ref_base == "A")) %>%
                                              mutate(count_CtoA = count_A * (ref_base == "C")) %>%
                                              mutate(count_CtoG = count_G * (ref_base == "C")) %>%
                                              mutate(count_CtoT = count_T * (ref_base == "C")) %>%
                                              mutate(count_GtoA = count_A * (ref_base == "G")) %>%
                                              mutate(count_GtoC = count_C * (ref_base == "G")) %>%
                                              mutate(count_GtoT = count_T * (ref_base == "G")) %>%
                                              mutate(count_TtoA = count_A * (ref_base == "T")) %>%
                                              mutate(count_TtoC = count_C * (ref_base == "T")) %>%
                                              mutate(count_TtoG = count_G * (ref_base == "T")) %>%
                                              mutate(pos_AtoC = replace_na((count_C/count_C) * (ref_base == "A"),0)) %>%
                                              mutate(pos_AtoG = replace_na((count_G/count_G) * (ref_base == "A"),0)) %>%
                                              mutate(pos_AtoT = replace_na((count_T/count_T) * (ref_base == "A"),0)) %>%
                                              mutate(pos_CtoA = replace_na((count_A/count_A) * (ref_base == "C"),0)) %>%
                                              mutate(pos_CtoG = replace_na((count_G/count_G) * (ref_base == "C"),0)) %>%
                                              mutate(pos_CtoT = replace_na((count_T/count_T) * (ref_base == "C"),0)) %>%
                                              mutate(pos_GtoA = replace_na((count_A/count_A) * (ref_base == "G"),0)) %>%
                                              mutate(pos_GtoC = replace_na((count_C/count_C) * (ref_base == "G"),0)) %>%
                                              mutate(pos_GtoT = replace_na((count_T/count_T) * (ref_base == "G"),0)) %>%
                                              mutate(pos_TtoA = replace_na((count_A/count_A) * (ref_base == "T"),0)) %>%
                                              mutate(pos_TtoC = replace_na((count_C/count_C) * (ref_base == "T"),0)) %>%
                                              mutate(pos_TtoG = replace_na((count_G/count_G) * (ref_base == "T"),0))
readcounts_all$SampleID <- str_replace_all(readcounts_all$SampleID, "-", "_")

# Get the mutation rate per 10k bases
per10k_mismatched <- readcounts_all %>% group_by(SampleID) %>% summarise(total_depth = sum(depth),
                                                                          total_mismatches = sum(count_mismatch),
                                                                          total_unique_mismatched = sum(count_unique_mismatch),
                                                                          per10k_mismatched = total_mismatches/total_depth *10000,
                                                                          total_positions = length(count_mismatch),
                                                                          positions_with_mismatch = sum(count_mismatch != 0, na.rm = TRUE),
                                                                          percent_positions_with_mismatch = positions_with_mismatch/total_positions * 100,
                                                                          total_transitions = sum(count_transitions),
                                                                          total_transversions = total_mismatches - total_transitions,
                                                                          it_ver_ratio = total_transitions/total_transversions,
                                                                          mean_depth = mean(depth),
                                                                          min_depth = min(depth),
                                                                          rmsd = sqrt(mean(freq_mismatch.squared)),
                                                                         EXPERIMENT = EXPERIMENT
                                                                        ) %>% unique()

# Split sampleid into useful categories: celltype, genotype, and biological replicate
per10k_mismatched <- per10k_mismatched %>% mutate(celltype = str_split_i(SampleID, "_", 5)) %>%
  mutate(genotype = str_split_i(SampleID, "_", 6)) %>%
  mutate(replicate = str_split_i(SampleID, "_", 7))

## MEAN DEPTH WHOLE GENOME
mean(per10k_mismatched$mean_depth)

# readcounts_all.plot <- readcounts_all %>% mutate(celltype = str_split_i(SampleID, "_", 5)) %>%
#   mutate(genotype = str_split_i(SampleID, "_", 6)) %>%
#   mutate(replicate = str_split_i(SampleID, "_", 7))
# depth.plot <- ggplot(readcounts_all.plot, aes(x=position, y = depth, color = celltype)) + geom_line()
# depth.plot

## MIN DEPTH WHOLE GENOME
min(per10k_mismatched$min_depth)

# Get frequency of each possible substitution over total nucleotides sequenced
SNP_stats <- readcounts_all %>% group_by(SampleID, EXPERIMENT) %>% summarise(total_depth = sum(depth),
  AtoC.freq = sum(count_AtoC)/total_depth * 10000,
  AtoG.freq = sum(count_AtoG)/total_depth * 10000,
  AtoT.freq = sum(count_AtoT)/total_depth * 10000,
  CtoA.freq = sum(count_CtoA)/total_depth * 10000,
  CtoG.freq = sum(count_CtoG)/total_depth * 10000,
  CtoT.freq = sum(count_CtoT)/total_depth * 10000,
  GtoA.freq = sum(count_GtoA)/total_depth * 10000,
  GtoC.freq = sum(count_GtoC)/total_depth * 10000,
  GtoT.freq = sum(count_GtoT)/total_depth * 10000,
  TtoA.freq = sum(count_TtoA)/total_depth * 10000,
  TtoC.freq = sum(count_TtoC)/total_depth * 10000,
  TtoG.freq = sum(count_TtoG)/total_depth * 10000) %>% pivot_longer(cols = c(AtoG.freq,
                                                          AtoC.freq,
                                                          AtoT.freq,
                                                          CtoA.freq,
                                                          CtoG.freq,
                                                          CtoT.freq,
                                                          GtoA.freq,
                                                          GtoC.freq,
                                                          GtoT.freq,
                                                          TtoA.freq,
                                                          TtoC.freq,
                                                          TtoG.freq), names_to = "substitution")
SNP_stats$SampleID <- str_replace_all(SNP_stats$SampleID, "-", "_")
SNP_stats <- SNP_stats %>% mutate(celltype = str_split_i(SampleID, "_", 5)) %>%
  mutate(genotype = str_split_i(SampleID, "_", 6)) %>%
  mutate(replicate = str_split_i(SampleID, "_", 7))

# Get frequency of each possible substitution over total number of mismatched positions
SNP_stats_pos <- readcounts_all %>% group_by(SampleID, EXPERIMENT) %>% summarise(total_positions_A = sum(count_mismatch > 0),
                                                                     total_positions_C = sum(count_mismatch > 0),
                                                                     total_positions_G = sum(count_mismatch > 0),
                                                                     total_positions_T = sum(count_mismatch > 0),
                                                                     AtoC = sum(pos_AtoC)/total_positions_A,
                                                                     AtoG = sum(pos_AtoG)/total_positions_A,
                                                                     AtoU = sum(pos_AtoT)/total_positions_A,
                                                                     CtoA = sum(pos_CtoA)/total_positions_C,
                                                                     CtoG = sum(pos_CtoG)/total_positions_C,
                                                                     CtoU = sum(pos_CtoT)/total_positions_C,
                                                                     GtoA = sum(pos_GtoA)/total_positions_G,
                                                                     GtoC = sum(pos_GtoC)/total_positions_G,
                                                                     GtoU = sum(pos_GtoT)/total_positions_G,
                                                                     UtoA = sum(pos_TtoA)/total_positions_T,
                                                                     UtoC = sum(pos_TtoC)/total_positions_T,
                                                                     UtoG = sum(pos_TtoG/total_positions_T)) %>% pivot_longer(cols = c(AtoG,
                                                                                                                               AtoC,
                                                                                                                               AtoU,
                                                                                                                               CtoA,
                                                                                                                               CtoG,
                                                                                                                               CtoU,
                                                                                                                               GtoA,
                                                                                                                               GtoC,
                                                                                                                               GtoU,
                                                                                                                               UtoA,
                                                                                                                               UtoC,
                                                                                                                               UtoG), names_to = "substitution")
test <- readcounts_all %>% group_by(SampleID) %>% summarise(total_positions_A = sum(count_mismatch > 0),
                                                                     total_positions_C = sum(count_mismatch > 0),
                                                                     total_positions_G = sum(count_mismatch > 0),
                                                                     total_positions_T = sum(count_mismatch > 0),
                                                                     AtoC = sum(pos_AtoC),
                                                                     AtoG = sum(pos_AtoG),
                                                                     AtoU = sum(pos_AtoT),
                                                                     CtoA = sum(pos_CtoA),
                                                                     CtoG = sum(pos_CtoG),
                                                                     CtoU = sum(pos_CtoT),
                                                                     GtoA = sum(pos_GtoA),
                                                                     GtoC = sum(pos_GtoC),
                                                                     GtoU = sum(pos_GtoT),
                                                                     UtoA = sum(pos_TtoA),
                                                                     UtoC = sum(pos_TtoC),
                                                                     UtoG = sum(pos_TtoG)) %>% pivot_longer(cols = c(AtoG,
                                                                                                                                       AtoC,
                                                                                                                                       AtoU,
                                                                                                                                       CtoA,
                                                                                                                                       CtoG,
                                                                                                                                       CtoU,
                                                                                                                                       GtoA,
                                                                                                                                       GtoC,
                                                                                                                                       GtoU,
                                                                                                                                       UtoA,
                                                                                                                                       UtoC,
                                                                                                                                       UtoG), names_to = "substitution")
test2 <- test %>% group_by(SampleID) %>% summarise(total = total_positions_A,
                                                   total2 = sum(value)) %>% unique()
SNP_stats_pos <- SNP_stats_pos %>% mutate(celltype = str_split_i(SampleID, "_", 5)) %>%
  mutate(genotype = str_split_i(SampleID, "_", 6)) %>%
  mutate(replicate = str_split_i(SampleID, "_", 7))
#########################################################################################################
# Plot Overall Frequency and position frequency of all possible substitutions for U2OS and C636
#########################################################################################################
SNP_stats$x <- paste0(SNP_stats$substitution, "_", SNP_stats$genotype)
x_breaks <- SNP_stats %>% filter(genotype == "Q192L") %>% select(x) %>% unique()
x_breaks_final <- x_breaks$x
x_labels <- str_split_i(x_breaks_final, ".freq", 1)


substitution_plot_U2OS <- ggplot(SNP_stats %>% filter(celltype %in% c("U2OS")), aes(x=x, y=value, color = genotype)) +
  geom_point(size=4, alpha=0.3) + 
  theme_classic(base_size = 24) +
  ylab("mutations per 10k") +
  xlab("") +
  labs(title= "U2OS") + 
  scale_x_discrete(breaks = x_breaks_final, labels = x_labels) +
  scale_y_continuous(limits = c(0,0.5)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 2, vjust = 0.5, size = 24))
substitution_plot_U2OS

substitution_plot_C636 <- ggplot(SNP_stats %>% filter(celltype %in% c("C636")), aes(x=x, y=value, color = genotype)) +
  geom_point(size=4, alpha=0.3) + 
  theme_classic(base_size = 24) +
  ylab("mutations per 10k") +
  xlab("") +
  labs(title= "C636") + 
  scale_x_discrete(breaks = x_breaks_final, labels = x_labels) +
  scale_y_continuous(limits = c(0,0.5)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 2, vjust = 0.5, size = 24))
substitution_plot_C636

SNP_stats_pos$x <- paste0(SNP_stats_pos$substitution, "_", SNP_stats_pos$genotype)
x_breaks <- SNP_stats_pos %>% filter(genotype == "Q192L") %>% select(x) %>% unique()
x_breaks_final <- x_breaks$x
x_labels <- str_split_i(x_breaks_final, "_", 1)

substitution_pos_plot_U2OS <- ggplot(SNP_stats_pos %>% filter(celltype %in% c("U2OS")), aes(x=x, y=value, color = genotype)) +
  geom_point(size=4, alpha=0.3) + 
  theme_classic(base_size = 24) +
  ylab("mutated site frequency") +
  xlab("") + 
  labs(title= "U2OS") + 
  scale_x_discrete(breaks = x_breaks_final, labels = x_labels) +
  scale_y_continuous(limits = c(0,0.4)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 2, vjust = 0.5, size = 24))
substitution_pos_plot_U2OS

ggsave(filename = "substitution_pos_U2OS.jpeg", plot = substitution_pos_plot_U2OS, width = 6, height= 6, device = "jpeg")

substitution_pos_plot_C636 <- ggplot(SNP_stats_pos %>% filter(celltype %in% c("C636")), aes(x=x, y=value, color = genotype)) +
  geom_point(size=4, alpha=0.3) + 
  theme_classic(base_size = 24) +
  ylab("mutated site frequency") +
  xlab("") + 
  labs(title= "C6/36") + 
  scale_x_discrete(breaks = x_breaks_final, labels = x_labels) +
  scale_y_continuous(limits = c(0,0.4)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 2, vjust = 0.5, size = 24))
substitution_pos_plot_C636
ggsave(filename = "substitution_pos_C636.jpeg", plot = substitution_pos_plot_C636, width = 6, height= 6, device = "jpeg")

substitution_pos_plot_P0stock <- ggplot(SNP_stats_pos %>% filter(celltype %in% c("P0stock")), aes(x=x, y=value, color = genotype)) +
  geom_point(size=4, alpha=0.3) + 
  theme_classic(base_size = 24) +
  ylab("mutated site frequency") +
  xlab("") + 
  labs(title= "P0 stock") + 
  scale_x_discrete(breaks = x_breaks_final, labels = x_labels) +
  scale_y_continuous(limits = c(0,0.4)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 2, vjust = 0.5, size = 24))
substitution_pos_plot_P0stock

substitution_pos_plot_RNA <- ggplot(SNP_stats_pos %>% filter(celltype %in% c("RNA")), aes(x=x, y=value, color = genotype)) +
  geom_point(size=4, alpha=0.3) + 
  theme_classic(base_size = 24) +
  ylab("mutated site frequency") +
  xlab("") + 
  labs(title= "RNA") + 
  scale_x_discrete(breaks = x_breaks_final, labels = x_labels) +
  scale_y_continuous(limits = c(0,0.4)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 2, vjust = 0.5, size = 24))
substitution_pos_plot_RNA
ggsave(filename = "substitution_pos_RNA.jpeg", plot = substitution_pos_plot_RNA, width = 6, height= 6, device = "jpeg")
res.aov.U2OS <- aov(value ~ genotype, data = SNP_stats_pos %>% filter(celltype == "C636", substitution == "GtoA"))
tukey_HSD.U2OS <- TukeyHSD(res.aov.U2OS)
############################################################
# get shannon entropy
############################################################
shannon_entropy <- per10k_mismatched %>% select(SampleID)
shannon_entropy$mean <- 69
shannon_entropy$sd <- 69

for(sample in shannon_entropy$SampleID){
  print(sample)
  shannon_matrix <- readcounts_all %>% filter(SampleID == sample) %>% select(count_A, count_C, count_G, count_T)
  shannon_matrix <- as.matrix(shannon_matrix)
  shannon_out <- diversity(shannon_matrix, type="entropy")
  shannon_entropy$mean[shannon_entropy$SampleID == sample] <- mean(shannon_out$entropy)
  shannon_entropy$sd[shannon_entropy$SampleID == sample] <- sd(shannon_out$entropy)
}
shannon_entropy <- shannon_entropy %>% mutate(celltype = str_split_i(SampleID, "_", 5)) %>%
  mutate(genotype = str_split_i(SampleID, "_", 6)) %>%
  mutate(replicate = str_split_i(SampleID, "_", 7))

shannon_entropy.summary <- shannon_entropy %>% group_by(celltype, genotype) %>% summarise(entropy.geom_mean = exp(mean(log(mean))),
                                                                                          entropy.geom_sd = exp(sd(log(mean))))
shannon_entropy.summary$x <- paste0(shannon_entropy.summary$celltype, "_",shannon_entropy.summary$genotype)
shannon_entropy.summary <- shannon_entropy.summary %>% mutate(x = factor(x, levels = c(shannon_entropy.summary$x[shannon_entropy.summary$celltype == "plasmid"],
                                                                                       shannon_entropy.summary$x[shannon_entropy.summary$celltype == "RNA"],
                                                                                       shannon_entropy.summary$x[shannon_entropy.summary$celltype == "P0stock"],
                                                                                       shannon_entropy.summary$x[shannon_entropy.summary$celltype == "U2OS"],
                                                                                       shannon_entropy.summary$x[shannon_entropy.summary$celltype == "C636"])))

shannon_entropy.plot <- ggplot(shannon_entropy.summary, aes(x = x, fill = genotype,y = entropy.geom_mean)) +
  geom_bar(stat = 'identity') + 
  ylab("shannon entropy") +
  labs(title = "Shannon entropy") +
  geom_errorbar(aes(x=x, ymin = entropy.geom_mean/entropy.geom_sd, ymax = entropy.geom_mean*entropy.geom_sd), width = 0.4, colour = "black") +
  theme_classic(base_size = 18) + theme(axis.text.x = element_text(angle=45, hjust=1))

shannon_entropy.plot
ggsave(filename="shannon_entropy_combo.jpeg", plot=shannon_entropy.plot, device="jpeg", width = 12, height = 6)

res.aov.U2OS <- aov(mean ~ genotype, data = shannon_entropy %>% filter(celltype == "U2OS"))
tukey_HSD.U2OS <- TukeyHSD(res.aov.U2OS)
############################################################
# Plot rmsd
############################################################
res.aov.U2OS <- aov(rmsd ~ genotype, data = per10k_mismatched %>% filter(celltype == "U2OS"))
tukey_HSD.U2OS <- TukeyHSD(res.aov.U2OS)

rmsd.summary <- per10k_mismatched %>% group_by(celltype, genotype) %>% summarise(rmsd.mean = exp(mean(log(rmsd))),
                                                                                rmsd.sd = exp(sd(log(rmsd))))

rmsd.summary$x <- paste0(rmsd.summary$celltype, "_",rmsd.summary$genotype)
rmsd.summary <- rmsd.summary %>% mutate(x = factor(x, levels = c(rmsd.summary$x[rmsd.summary$celltype == "plasmid"],
                                                                 rmsd.summary$x[rmsd.summary$celltype == "RNA"],
                                                                 rmsd.summary$x[rmsd.summary$celltype == "P0stock"],
                                                                 rmsd.summary$x[rmsd.summary$celltype == "U2OS"],
                                                                 rmsd.summary$x[rmsd.summary$celltype == "C636"])))
rmsd.plot <- ggplot(rmsd.summary, aes(x = x, fill = genotype, y = rmsd.mean)) +
  geom_bar(stat = 'identity') +
  ylab("rmsd") +
  labs(title = "Root mean squared deviation (rmsd)") +
  geom_errorbar(aes(x=x, ymin = rmsd.mean / rmsd.sd, ymax = rmsd.mean * rmsd.sd), width = 0.4, colour = "black") +
  # scale_y_continuous(limits = c(0,.015)) +
  theme_classic(base_size = 18) + theme(axis.text.x = element_text(angle=45, hjust=1))
rmsd.plot

ggsave(filename="rmsd_combo.jpeg", plot=rmsd.plot, device="jpeg", width = 12, height = 6)

#############################################
# Plot Frequency of SNPs for U2OS and C636
#############################################
min_mismatches <- 0
max_mismatches <- per10k_mismatched %>% ungroup() %>% select(per10k_mismatched) %>% max() * 1.1
min_positions <- 0
max_positions <- per10k_mismatched %>% ungroup() %>% select(percent_positions_with_mismatch) %>% max() * 1.1
min_transitions <- 0
max_transitions <- per10k_mismatched %>% ungroup() %>% select(it_ver_ratio) %>% max() * 1.1

mismatched.means <- per10k_mismatched %>% group_by(celltype, genotype) %>% summarise(per10k_mismatched.mean = mean(per10k_mismatched),
                                                                                     per10k_mismatched.sd = sd(per10k_mismatched),
                                                                                     percent_positions_with_mismatch.mean = mean(percent_positions_with_mismatch),
                                                                                     percent_positions_with_mismatch.sd = sd(percent_positions_with_mismatch),
                                                                                     it_ver_ratio.mean = mean(it_ver_ratio),
                                                                                     it_ver_ratio.sd = sd(it_ver_ratio))
mismatched.means$x <- paste0(mismatched.means$celltype, "_",mismatched.means$genotype)
mismatched.means <- mismatched.means %>% mutate(x = factor(x, levels = c(mismatched.means$x[mismatched.means$celltype == "plasmid"],
                                                                         mismatched.means$x[mismatched.means$celltype == "RNA"],
                                                                         mismatched.means$x[mismatched.means$celltype == "P0stock"],
                                                                         mismatched.means$x[mismatched.means$celltype == "U2OS"],
                                                                         mismatched.means$x[mismatched.means$celltype == "C636"])))

plot_percent_mismatches <- ggplot(mismatched.means,aes(x=x, y=per10k_mismatched.mean, fill=genotype)) + 
  geom_bar(stat='identity') + 
  geom_errorbar(aes(x=x, ymin=per10k_mismatched.mean-per10k_mismatched.sd, ymax=per10k_mismatched.mean+per10k_mismatched.sd), width=0.4, colour="black") +
  theme_classic(base_size = 24) +
  ylab("mutations per 10k bases") +
  xlab("") +
  labs(title = "Overall mutation frequency") + 
  scale_y_continuous(limits = c(min_mismatches,ceiling(max_mismatches))) +
  theme(legend.position = "bottom", 
        axis.text.x= element_text(angle = 45, hjust =1), 
        axis.text.y = element_text(size = 30))

plot_percent_mismatches
ggsave(filename = "mutations_per_10k_combo.jpeg", plot = plot_percent_mismatches, device = "jpeg", width = 10, height = 8)

plot_percent_positions_with_mismatch <- ggplot(mismatched.means, aes(x=x, y=percent_positions_with_mismatch.mean, fill=genotype)) + 
  geom_bar(stat='identity') + 
  geom_errorbar(aes(x=x, ymin=percent_positions_with_mismatch.mean-percent_positions_with_mismatch.sd, ymax=percent_positions_with_mismatch.mean+percent_positions_with_mismatch.sd), width=0.4, colour="black") +
  theme_classic(base_size = 24) +
  ylab("% positions with mutation") +
  labs(title = "Mutations by position") + 
  xlab("") +
  scale_y_continuous(limits = c(min_positions,max_positions)) +
  theme(legend.position = "bottom",
        axis.text.x= element_text(angle = 45, hjust =1),
        axis.text.y = element_text(size = 30))
plot_percent_positions_with_mismatch
ggsave(filename = "position_mismatches_combo.jpeg", plot = plot_percent_positions_with_mismatch, device = "jpeg", width = 10, height = 8)

res.aov.U2OS <- aov(per10k_mismatched ~ genotype, data = per10k_mismatched %>% filter(celltype == "U2OS"))
tukey_HSD.U2OS <- TukeyHSD(res.aov.U2OS)

#############################################
# Plot Frequency of SNPs for U2OS and C636 SPLIT BY EXPERIMENT
#############################################
min_mismatches <- 0
max_mismatches <- per10k_mismatched %>% ungroup() %>% select(per10k_mismatched) %>% max() * 1.1
min_positions <- 0
max_positions <- per10k_mismatched %>% ungroup() %>% select(percent_positions_with_mismatch) %>% max() * 1.1
min_transitions <- 0
max_transitions <- per10k_mismatched %>% ungroup() %>% select(it_ver_ratio) %>% max() * 1.1

mismatched.means <- per10k_mismatched %>% ungroup() %>% group_by(celltype, genotype, EXPERIMENT) %>% summarise(per10k_mismatched.mean = mean(per10k_mismatched),
                                                                                     per10k_mismatched.sd = sd(per10k_mismatched),
                                                                                     percent_positions_with_mismatch.mean = mean(percent_positions_with_mismatch),
                                                                                     percent_positions_with_mismatch.sd = sd(percent_positions_with_mismatch),
                                                                                     it_ver_ratio.mean = mean(it_ver_ratio),
                                                                                     it_ver_ratio.sd = sd(it_ver_ratio))
mismatched.means$x <- paste0(mismatched.means$celltype, "_",mismatched.means$genotype, "_", mismatched.means$EXPERIMENT)
mismatched.means <- mismatched.means %>% mutate(x = factor(x, levels = c(mismatched.means$x[mismatched.means$celltype == "plasmid"],
                                                                         mismatched.means$x[mismatched.means$celltype == "RNA"],
                                                                         mismatched.means$x[mismatched.means$celltype == "P0stock"],
                                                                         mismatched.means$x[mismatched.means$celltype == "U2OS"],
                                                                         mismatched.means$x[mismatched.means$celltype == "C636"])))

plot_percent_mismatches <- ggplot(mismatched.means,aes(x=x, y=per10k_mismatched.mean, fill=genotype)) + 
  geom_bar(stat='identity') + 
  geom_errorbar(aes(x=x, ymin=per10k_mismatched.mean-per10k_mismatched.sd, ymax=per10k_mismatched.mean+per10k_mismatched.sd), width=0.4, colour="black") +
  theme_classic(base_size = 24) +
  ylab("mutations per 10k bases") +
  xlab("") +
  labs(title = "Overall mutation frequency") + 
  scale_y_continuous(limits = c(min_mismatches,ceiling(max_mismatches))) +
  theme(legend.position = "bottom", 
        axis.text.x= element_text(angle = 45, hjust =1), 
        axis.text.y = element_text(size = 30))

plot_percent_mismatches
ggsave(filename = "mutations_per_10k_combo_split.jpeg", plot = plot_percent_mismatches, device = "jpeg", width = 10, height = 8)

plot_percent_positions_with_mismatch <- ggplot(mismatched.means, aes(x=x, y=percent_positions_with_mismatch.mean, fill=genotype)) + 
  geom_bar(stat='identity') + 
  geom_errorbar(aes(x=x, ymin=percent_positions_with_mismatch.mean-percent_positions_with_mismatch.sd, ymax=percent_positions_with_mismatch.mean+percent_positions_with_mismatch.sd), width=0.4, colour="black") +
  theme_classic(base_size = 24) +
  ylab("% positions with mutation") +
  labs(title = "Mutations by position") + 
  xlab("") +
  scale_y_continuous(limits = c(min_positions,max_positions)) +
  theme(legend.position = "bottom",
        axis.text.x= element_text(angle = 45, hjust =1),
        axis.text.y = element_text(size = 30))
plot_percent_positions_with_mismatch
ggsave(filename = "position_mismatches_combo_split.jpeg", plot = plot_percent_positions_with_mismatch, device = "jpeg", width = 10, height = 8)

##########################################################################################
# Get mean # SNPs per allele frequency bin 
##########################################################################################

SNPs_per_AF <- readcounts_all %>% group_by(SampleID) %>% summarise(less_1 = sum(freq_mismatch > 0 & freq_mismatch < 1),
                                                                          one_to_five = sum(freq_mismatch >=1 & freq_mismatch < 5),
                                                                          greater_5 = sum(freq_mismatch >= 5 & freq_mismatch < 10),
                                                                          greater_10 = sum(freq_mismatch >= 10))

SNPs_per_AF$SampleID <- str_replace_all(SNPs_per_AF$SampleID, "-", "_")

SNPs_per_AF <- SNPs_per_AF %>% mutate(celltype = str_split_i(SampleID, "_", 5)) %>%
  mutate(genotype = str_split_i(SampleID, "_", 6)) %>%
  mutate(replicate = str_split_i(SampleID, "_", 7))

SNPs_per_AF.mean <- SNPs_per_AF %>% 
  group_by(celltype, genotype) %>%
  summarise(a_less_1 = mean(less_1),
            b_one_to_five = mean(one_to_five),
            c_greater_5 = mean(greater_5),
            d_greater_10 = mean(greater_10),
            n = n()) %>%
             pivot_longer(cols = c(a_less_1, b_one_to_five, c_greater_5, d_greater_10), names_to = "freq_bin", values_to = "value")

SNPs_per_AF.sd <- SNPs_per_AF %>% 
  group_by(celltype, genotype) %>%
  summarise(a_less_1 = sd(less_1),
            b_one_to_five = sd(one_to_five),
            c_greater_5 = sd(greater_5),
            d_greater_10 = sd(greater_10),
            n = n()) %>%
       pivot_longer(cols = c(a_less_1, b_one_to_five, c_greater_5, d_greater_10), names_to = "freq_bin", values_to = "value.sd")

SNPs_per_AF <- merge(SNPs_per_AF.mean, SNPs_per_AF.sd, by=c("celltype", "genotype", "freq_bin"))

SNPs_per_AF$x <- paste0(SNPs_per_AF$celltype, "_",SNPs_per_AF$genotype)
SNPs_per_AF <- SNPs_per_AF %>% mutate(x = factor(x, levels = c(unique(SNPs_per_AF$x[SNPs_per_AF$celltype == "plasmid"]),
                                                                             unique(SNPs_per_AF$x[SNPs_per_AF$celltype == "RNA"]),
                                                                             unique(SNPs_per_AF$x[SNPs_per_AF$celltype == "P0stock"]),
                                                                             unique(SNPs_per_AF$x[SNPs_per_AF$celltype == "U2OS"]),
                                                                             unique(SNPs_per_AF$x[SNPs_per_AF$celltype == "C636"]))))

SNPs_per_AF.plot <- ggplot(SNPs_per_AF, aes(fill=freq_bin, y=value, x=x)) + 
  geom_bar(position= position_dodge(),stat="identity") +
  geom_errorbar(position = position_dodge(width=0.9),stat="identity",aes(x=x, ymin=value-value.sd, ymax=value+value.sd), width=0.4, colour="black") +
  theme_classic(base_size = 18) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_fill_manual(name = "Allele frequency", values = c("d_greater_10"="black",
                                                          "c_greater_5"="red",
                                                          "b_one_to_five"="blue",
                                                          "a_less_1"="green3"), 
                    labels = c("d_greater_10"=">10%",
                                                          "c_greater_5"=">5%",
                                                          "b_one_to_five"=">1%",
                                                          "a_less_1"="<1%")) +
  # scale_y_continuous(limits = c(0,30), breaks = seq(0,30,by=5)) +
  labs(title = "Mean SNPs colored by allele frequency") +
  ylab("mean # SNPs")

SNPs_per_AF.plot
ggsave(filename="SNPs_per_AF_combo.jpeg", plot=SNPs_per_AF.plot, device="jpeg", width = 12, height = 6)

SNPs_per_AF_over1 <- SNPs_per_AF %>% filter(freq_bin != "a_less_1")

SNPs_per_AF_over1.plot <- ggplot(SNPs_per_AF_over1, aes(fill=freq_bin, y=value, x=x)) + 
  geom_bar(position= position_dodge(),stat="identity") +
  geom_errorbar(position = position_dodge(width=0.9),stat="identity",aes(x=x, ymin=value-value.sd, ymax=value+value.sd), width=0.4, colour="black") +
  theme_classic(base_size = 18) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_fill_manual(name = "Allele frequency", values = c("d_greater_10"="black",
                                                          "c_greater_5"="red",
                                                          "b_one_to_five"="blue",
                                                          "a_less_1"="green3"), 
                    labels = c("d_greater_10"=">10%",
                               "c_greater_5"=">5%",
                               "b_one_to_five"=">1%",
                               "a_less_1"="<1%")) +
  # scale_y_continuous(limits = c(0,30), breaks = seq(0,30,by=5)) +
  labs(title = "Mean SNPs colored by allele frequency") +
  ylab("mean # SNPs")

SNPs_per_AF_over1.plot
write_csv(SNPs_per_AF, file = "Fig5D_SNVs_allele_frequency_bins.csv")
ggsave(filename="SNPs_per_AF_over1_combo.jpeg", plot=SNPs_per_AF_over1.plot, device="jpeg", width = 12, height = 6)
##########################################################################################
# Get mean # SNPs per allele frequency bin split by EXPERIMENT
##########################################################################################

SNPs_per_AF <- readcounts_all %>% group_by(SampleID, EXPERIMENT) %>% summarise(less_1 = sum(freq_mismatch > 0 & freq_mismatch < 1),
                                                                   one_to_five = sum(freq_mismatch >=1 & freq_mismatch < 5),
                                                                   greater_5 = sum(freq_mismatch >= 5 & freq_mismatch < 10),
                                                                   greater_10 = sum(freq_mismatch >= 10))
                                                            

SNPs_per_AF$SampleID <- str_replace_all(SNPs_per_AF$SampleID, "-", "_")

SNPs_per_AF <- SNPs_per_AF %>% mutate(celltype = str_split_i(SampleID, "_", 5)) %>%
  mutate(genotype = str_split_i(SampleID, "_", 6)) %>%
  mutate(replicate = str_split_i(SampleID, "_", 7))

SNPs_per_AF.mean <- SNPs_per_AF %>% 
  group_by(celltype, genotype, EXPERIMENT) %>%
  summarise(a_less_1 = mean(less_1),
            b_one_to_five = mean(one_to_five),
            c_greater_5 = mean(greater_5),
            d_greater_10 = mean(greater_10),
            n = n()) %>%
  pivot_longer(cols = c(a_less_1, b_one_to_five, c_greater_5, d_greater_10), names_to = "freq_bin", values_to = "value")

SNPs_per_AF.sd <- SNPs_per_AF %>% 
  group_by(celltype, genotype, EXPERIMENT) %>%
  summarise(a_less_1 = sd(less_1),
            b_one_to_five = sd(one_to_five),
            c_greater_5 = sd(greater_5),
            d_greater_10 = sd(greater_10),
            n = n()) %>%
  pivot_longer(cols = c(a_less_1, b_one_to_five, c_greater_5, d_greater_10), names_to = "freq_bin", values_to = "value.sd")

SNPs_per_AF <- merge(SNPs_per_AF.mean, SNPs_per_AF.sd, by=c("celltype", "genotype", "EXPERIMENT", "freq_bin"))

SNPs_per_AF$x <- paste0(SNPs_per_AF$celltype, "_",SNPs_per_AF$genotype,"_",SNPs_per_AF$EXPERIMENT)
SNPs_per_AF <- SNPs_per_AF %>% mutate(x = factor(x, levels = c(unique(SNPs_per_AF$x[SNPs_per_AF$celltype == "plasmid"]),
                                                               unique(SNPs_per_AF$x[SNPs_per_AF$celltype == "RNA"]),
                                                               unique(SNPs_per_AF$x[SNPs_per_AF$celltype == "P0stock"]),
                                                               unique(SNPs_per_AF$x[SNPs_per_AF$celltype == "U2OS"]),
                                                               unique(SNPs_per_AF$x[SNPs_per_AF$celltype == "C636"]))))

SNPs_per_AF.plot <- ggplot(SNPs_per_AF, aes(fill=freq_bin, y=value, x=x)) + 
  geom_bar(position= position_dodge(),stat="identity") +
  geom_errorbar(position = position_dodge(width=0.9),stat="identity",aes(x=x, ymin=value-value.sd, ymax=value+value.sd), width=0.4, colour="black") +
  theme_classic(base_size = 18) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_fill_manual(name = "Allele frequency", values = c("d_greater_10"="black",
                                                          "c_greater_5"="red",
                                                          "b_one_to_five"="blue",
                                                          "a_less_1"="green3"), 
                    labels = c("d_greater_10"=">10%",
                               "c_greater_5"=">5%",
                               "b_one_to_five"=">1%",
                               "a_less_1"="<1%")) +
  # scale_y_continuous(limits = c(0,30), breaks = seq(0,30,by=5)) +
  labs(title = "Mean SNPs colored by allele frequency") +
  ylab("mean # SNPs")

SNPs_per_AF.plot
ggsave(filename="SNPs_per_AF_combo_split.jpeg", plot=SNPs_per_AF.plot, device="jpeg", width = 12, height = 6)

SNPs_per_AF_over1 <- SNPs_per_AF %>% filter(freq_bin != "a_less_1")

SNPs_per_AF_over1.plot <- ggplot(SNPs_per_AF_over1, aes(fill=freq_bin, y=value, x=x)) + 
  geom_bar(position= position_dodge(),stat="identity") +
  geom_errorbar(position = position_dodge(width=0.9),stat="identity",aes(x=x, ymin=value-value.sd, ymax=value+value.sd), width=0.4, colour="black") +
  theme_classic(base_size = 18) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_fill_manual(name = "Allele frequency", values = c("d_greater_10"="black",
                                                          "c_greater_5"="red",
                                                          "b_one_to_five"="blue",
                                                          "a_less_1"="green3"), 
                    labels = c("d_greater_10"=">10%",
                               "c_greater_5"=">5%",
                               "b_one_to_five"=">1%",
                               "a_less_1"="<1%")) +
  # scale_y_continuous(limits = c(0,30), breaks = seq(0,30,by=5)) +
  labs(title = "Mean SNPs colored by allele frequency") +
  ylab("mean # SNPs")

SNPs_per_AF_over1.plot
write_csv(SNPs_per_AF, file = "Fig5D_SNVs_allele_frequency_bins_per_experiment.csv")
ggsave(filename="SNPs_per_AF_over1_combo_split.jpeg", plot=SNPs_per_AF_over1.plot, device="jpeg", width = 12, height = 6)




#######################################################
# make tables for peiqi to use in prism
#######################################################
# Overall mutation frequency 

PY_mut_per10k <- per10k_mismatched %>% select(SampleID, EXPERIMENT, celltype, genotype, replicate, total_depth, total_mismatches, per10k_mismatched)
PY_mut_per10k.summary.per.exp <- PY_mut_per10k %>% group_by(EXPERIMENT, celltype, genotype) %>% summarise(n= n(),
                                                                                                          total_depth.mean = mean(total_depth),
                                                                                                          total_mismatches.mean = mean(total_mismatches),
                                                                                                          per10k_mismatched.mean = mean(per10k_mismatched),
                                                                                                          per10k_mismatched.sd = sd(per10k_mismatched))
PY_mut_per10k.summary <- PY_mut_per10k %>% group_by(celltype, genotype) %>% summarise(n= n(),
                                                                                               total_depth.mean = mean(total_depth),
                                                                                               total_mismatches.mean = mean(total_mismatches),
                                                                                               per10k_mismatched.mean = mean(per10k_mismatched),
                                                                                               per10k_mismatched.sd = sd(per10k_mismatched))

write_csv(PY_mut_per10k, file = "Fig5A_mutations_per_10k_per_sample.csv")
write_csv(PY_mut_per10k.summary.per.exp, file = "Fig5A_mutations_per_10k_per_experiment.csv")
write_csv(PY_mut_per10k.summary, file = "Fig5A_mutations_per_10k.csv")

#######################################################
# % Positions with mutation
#######################################################

PY_mut_pos <- per10k_mismatched %>% select(SampleID, EXPERIMENT, celltype, genotype, replicate, total_positions, positions_with_mismatch, percent_positions_with_mismatch)
PY_mut_pos.summary.per.exp <- PY_mut_pos %>% group_by(EXPERIMENT, celltype, genotype) %>% summarise(n= n(),
                                                                                                          total_positions.mean = mean(total_positions),
                                                                                                          positions_with_mismatch.mean = mean(positions_with_mismatch),
                                                                                                          percent_positions_with_mismatch.mean = mean(percent_positions_with_mismatch),
                                                                                                          percent_positions_with_mismatch.sd = sd(percent_positions_with_mismatch))
PY_mut_pos.summary <- PY_mut_pos %>% group_by(celltype, genotype) %>% summarise(n= n(),
                                                                                               total_positions.mean = mean(total_positions),
                                                                                               positions_with_mismatch.mean = mean(positions_with_mismatch),
                                                                                               percent_positions_with_mismatch.mean = mean(percent_positions_with_mismatch),
                                                                                               percent_positions_with_mismatch.sd = sd(percent_positions_with_mismatch))

write_csv(PY_mut_pos, file = "Fig5B_percent_positions_with_mutation_per_sample.csv")
write_csv(PY_mut_pos.summary.per.exp, file = "Fig5B_percent_positions_with_mutation_per_experiment.csv")
write_csv(PY_mut_pos.summary, file = "Fig5B_percent_positions_with_mutation.csv")

#######################################################
# Shannon entropy
#######################################################
shannon_entropy <- merge(shannon_entropy, PY_mut_pos)
PY_shannon_entropy.summary.per.exp <- shannon_entropy %>% group_by(EXPERIMENT, celltype, genotype) %>% summarise(n= n(),
                                                                                                                 shannon_entropy.geom_mean = exp(mean(log(mean))),
                                                                                                                 shannon_entropy.geom_sd = exp(sd(log(mean))))
PY_shannon_entropy.summary <- shannon_entropy %>% group_by(celltype, genotype) %>% summarise(n= n(),
                                                                                                         shannon_entropy.geom_mean = exp(mean(log(mean))),
                                                                                                         shannon_entropy.geom_sd = exp(sd(log(mean))))
write_csv(shannon_entropy, file = "Fig5C_shannon_entropy_per_sample.csv")
write_csv(PY_shannon_entropy.summary.per.exp, file = "Fig5C_shannon_entropy_per_experiment.csv")
write_CSV(PY_shannon_entropy.summary, file="Fig5C_shannon_entropy.csv")

#######################################################
# mutation bias overall 
#######################################################
PY_SNP_stats.summary.per.experiment <- SNP_stats %>% group_by(EXPERIMENT, celltype, genotype,substitution) %>% summarise(n = n(),
                                                                                                            freq.mean = mean(value),
                                                                                                            freq.sd = sd(value),
                                                                                                            total_depth.mean = mean(total_depth))

PY_SNP_stats.summary <- SNP_stats %>% group_by(celltype, genotype,substitution) %>% summarise(n = n(),
                                                                                              freq.mean = mean(value),
                                                                                              freq.sd = sd(value),
                                                                                              total_depth.mean = mean(total_depth))

write_csv(SNP_stats, file = "Fig5EandF_mutation_bias_per_10k_per_sample.csv")
write_csv(PY_SNP_stats.summary.per.experiment, file = "Fig5EandF_mutation_bias_per_10k_per_experiment.csv")
write_csv(PY_SNP_stats.summary, file = "Fig5EandF_mutation_bias_per_10k.csv")


#######################################################
# mutation bias by positions
#######################################################
PY_SNP_stats_pos <- SNP_stats_pos %>% mutate(total_mutated_positions = total_positions_A) %>% select(SampleID, EXPERIMENT, celltype, genotype, replicate, substitution, value, total_mutated_positions, x)

PY_SNP_stats_pos.summary.per.experiment <- PY_SNP_stats_pos %>% group_by(EXPERIMENT, celltype, genotype,substitution) %>% summarise(n = n(),
                                                                                                                         freq.mean = mean(value),
                                                                                                                         freq.sd = sd(value),
                                                                                                                         total_mutated_positions.mean = mean(total_mutated_positions))

PY_SNP_stats_pos.summary <- PY_SNP_stats_pos %>% group_by(celltype, genotype,substitution) %>% summarise(n = n(),
                                                                                              freq.mean = mean(value),
                                                                                              freq.sd = sd(value),
                                                                                              total_mutated_positions.mean = mean(total_mutated_positions))

write_csv(PY_SNP_stats_pos, file = "Fig5EandF_mutation_bias_per_mutated-positions_per_sample.csv")
write_csv(PY_SNP_stats_pos.summary.per.experiment, file = "Fig5EandF_mutation_bias_per_mutated-positions_per_experiment.csv")
write_csv(PY_SNP_stats_pos.summary, file = "Fig5EandF_mutation_bias_per_mutated-positions.csv")

