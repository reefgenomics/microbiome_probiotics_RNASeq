#!/usr/bin/Rscript

# Load required libraries
library(pacman)
p_load(phyloseq, microbiome, httpgd, ggrepel, ggplot2, ggpubr, dplyr, readr)

# Load count data and preprocess
cnt <- read.table("/home/colinl/Proj/microbiome_probiotics_RNASeq/results/star_salmon/salmon.merged.gene_counts.tsv", 
          header = TRUE, sep = "\t", row.names = 1)
cnt <- cnt[, -1]  # Remove the first column (unnecessary)
cnt <- round(cnt[, order(colnames(cnt))], digits = 0)  # Round counts and reorder columns alphabetically

# Fix column names by removing prefixes
cnt <- cnt %>% rename_all(~sub("X", "", .))  # Remove "X" prefix
cnt <- cnt %>% rename_all(~sub("2024_Pver_", "", .))  # Remove "2024_Pver_" prefix # for cleaner names

# Load and clean tx2gene mapping
tx2gene <- read_tsv("/home/colinl/Proj/microbiome_probiotics_RNASeq/results/star_salmon/tx2gene.tsv")

# Handle improperly formatted pseudogenes
tx2gene <- tx2gene %>%
  mutate(gene_id = ifelse(grepl("agat-pseudogene-", gene_id), transcript_id, gene_id),
     transcript_id = ifelse(grepl("agat-pseudogene-", gene_name), gene_name, transcript_id))

# Remove "gene-" prefix from gene IDs
tx2gene$gene_id <- sub("gene-", "", tx2gene$gene_id)

# Match row names of count data to gene IDs in tx2gene
rownames(cnt) <- tx2gene$gene_id[match(rownames(cnt), tx2gene$gene_name)]

# Load sample metadata
samples <- read.csv("/home/colinl/Proj/microbiome_probiotics_RNASeq/samplesheet.csv", header = TRUE)
samples <- samples$sample

# Create metadata dataframe
metadata <- data.frame(
  sample = samples,
  sample_clean = sub("2024_Pver_", "", samples),  # Clean sample names
  treatment = factor(sub(".*_(CC|CP|VD|VL)_.*", "\\1", samples), levels = c("CC", "CP", "VD", "VL")),
  timepoint = factor(sub(".*_T(\\d+)_.*", "T\\1", samples), levels = c("T0", "T1", "T2")),
  colony = factor(sub(".*_T\\d+_(\\d+)_\\d+", "\\1", samples))
)

# Ensure metadata columns are factors
metadata$timepoint <- as.factor(metadata$timepoint)
metadata$treatment <- as.factor(metadata$treatment)
metadata$colony <- as.factor(metadata$colony)

# Set row names of metadata to match cleaned sample names
rownames(metadata) <- metadata$sample_clean

# Create phyloseq object
cnt.p <- otu_table(cnt, taxa_are_rows = TRUE)
sam.p <- sample_data(metadata)
phy <- phyloseq(cnt.p, sam.p)

# Transform data using log10 transformation
phy.t <- microbiome::transform(phy, transform = "log10", target = "OTU", shift = 0, scale = 1)

# Perform ordination analysis using RDA method and Euclidean distance
ordi <- ordinate(phy.t, method = "RDA", distance = "euclidean")

# Define color palettes for treatments and timepoints
col_treatment <- c("CC" = "#065dd6", "CP" = "#06d6a0", "VD" = "#ffd166", "VL" = "#ef476f")
col_timepoint <- c("T0" = "#06d6a0", "T1" = "#ffd166", "T2" = "#ef476f")

# Plot overall ordination
plot_all <- plot_ordination(phy.t, ordi, color = "treatment", shape = "timepoint") +
  geom_point(size = 3, alpha = 1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(aes(label = sample_clean)) +
  theme(legend.position = "left") +
  scale_color_manual(values = col_treatment) +
  labs(title = "Ordination plot")

# Generate plots for each treatment
plot_treatment_list <- list()
for (i in c("CC", "CP", "VD", "VL")) {
  col1 <- subset_samples(phy.t, treatment == i)
  ordi1 <- ordinate(col1, method = "RDA", distance = "euclidean")
  
  plot1 <- plot_ordination(col1, ordi1, color = "timepoint") +
  geom_point(size = 3, alpha = 1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(aes(label = sample_clean)) +
  theme(legend.position = "left") +
  scale_color_manual(values = col_timepoint) +
  labs(title = paste("Treatment ", i, sep = ""))
  
  plot_treatment_list[[i]] <- plot1
}

# Generate plots for each timepoint
plot_time_list <- list()
for (i in c("T0", "T1", "T2")) {
  col1 <- subset_samples(phy.t, timepoint == i)
  ordi1 <- ordinate(col1, method = "RDA", distance = "euclidean")
  
  plot1 <- plot_ordination(col1, ordi1, color = "treatment") +
  geom_point(size = 3, alpha = 1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(aes(label = sample_clean)) +
  theme(legend.position = "left") +
  scale_color_manual(values = col_treatment) +
  labs(title = paste("Timepoint ", i, sep = ""))
  
  plot_time_list[[i]] <- plot1
}

# Activate httpgd graphics device
# httpgd::hgd()

# Arrange and display treatment plots
ggarrange(plot_treatment_list[[1]], plot_treatment_list[[2]],
      plot_treatment_list[[3]], plot_treatment_list[[4]],
      ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

# Arrange and display combined plots
ggarrange(plot_all,
      ggarrange(plot_time_list[[1]], plot_time_list[[2]], plot_time_list[[3]],
          ncol = 3, nrow = 1, common.legend = TRUE, legend = "none"),
      ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")

# Save the combined plots as a PDF
pdf(file = "/home/colinl/Proj/microbiome_probiotics_RNASeq/results/deseq2_R/plots/ordination_plot.pdf", 
  width = 11, height = 8.5, pointsize = 6)

ggarrange(plot_all,
      ggarrange(plot_time_list[[1]], plot_time_list[[2]], plot_time_list[[3]],
          ncol = 3, nrow = 1, common.legend = TRUE, legend = "none"),
      ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")

dev.off()  # Close the PDF device
