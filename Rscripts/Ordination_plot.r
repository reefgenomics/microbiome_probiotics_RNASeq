#!/usr/bin/Rscript
# Load required libraries
library(pacman)
p_load(phyloseq, microbiome)
# graphics libraries
p_load(httpgd, ggrepel, ggplot2, ggpubr)
# library(patchwork)

cnt <- read.table("/home/colinl/Proj/microbiome_probiotics_RNASeq/results/star_salmon/salmon.merged.gene_counts.tsv", header = TRUE, sep = "\t", row.names = 1)
cnt <- cnt[, -1]
cnt <- round(cnt[,order(colnames(cnt))], digits = 0)
# fix colnames names remove the X beofre the 2024
cnt <- cnt %>% rename_all(~sub("X", "", .))
cnt <- cnt %>% rename_all(~sub("2024_Pver_", "", .))

tx2gene <- read_tsv("/home/colinl/Proj/microbiome_probiotics_RNASeq/results/star_salmon/tx2gene.tsv")
# clean up inproper formatting of pseudogenes agat converison 
# when the gene_id start with "agat-pseudogene-" invert gene_id and transcript_id columns
tx2gene <- tx2gene %>%
  mutate(gene_id = ifelse(grepl("agat-pseudogene-", gene_id), transcript_id, gene_id),
         transcript_id = ifelse(grepl("agat-pseudogene-", gene_name), gene_name, transcript_id))
#remove gene- prefix
tx2gene$gene_id <- sub("gene-", "", tx2gene$gene_id)

#fix cnt rownames matchin tx2gene gene_name to the colname and writing gene_id
rownames(cnt) <- tx2gene$gene_id[match(rownames(cnt), tx2gene$gene_name)]

samples <- read.csv("/home/colinl/Proj/microbiome_probiotics_RNASeq/samplesheet.csv", header = TRUE)
samples <- samples$sample

#make sample dataframe with colnames from cnt
metadata <- data.frame(
  sample = samples,
  sample_clean = sub("2024_Pver_","", samples),
  treatment = factor(sub(".*_(CC|CP|VD|VL)_.*", "\\1", samples),
                     levels = c("CC", "CP", "VD", "VL")),
  timepoint = factor(sub(".*_T(\\d+)_.*", "T\\1", samples),
                     levels = c("T0", "T1", "T2")),
  colony = factor(sub(".*_T\\d+_(\\d+)_\\d+", "\\1", samples))
)

#set day, treatment and replicate as factor
metadata$timepoint <- as.factor(metadata$timepoint)
metadata$Treatment <- as.factor(metadata$Treatment)
metadata$colony <- as.factor(metadata$colony)

rownames(metadata) <- metadata$sample_clean

# Create phyloseq object
cnt.p <- otu_table(cnt, taxa_are_rows = TRUE)
sam.p <- sample_data(metadata)
phy <- phyloseq(cnt.p, sam.p)

# Transform data using log10 transformation
phy.t <- microbiome::transform(phy, transform = "log10", target = "OTU", shift = 0, scale = 1)

# Perform ordination analysis using RDA method and euclidean distance
ordi <- ordinate(phy.t, method = "RDA", distance = "euclidean")

# Plot ordination results with color and shape based on metadata variables
plot_ordination(phy.t, ordi, color = "treatment", shape = "timepoint") +
  geom_point(size = 3, alpha = 1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

plot_all <- plot_ordination(phy.t, ordi, color = "treatment", shape = "timepoint") +
    geom_point(size = 3, alpha = 1) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_text_repel(aes(label = sample_clean)) +
    theme(legend.position = "left") +
    scale_color_manual(values = col_treatment) +
    labs(title = paste("Ordination plot", sep = ""))

col_treatment <- c("CC" = "#065dd6", "CP" = "#06d6a0" ,"VD" = "#ffd166","VL" = "#ef476f")
col_timepoint <- c("T0" = "#06d6a0","T1" = "#ffd166","T2" = "#ef476f")

plot_treatment_list <- list()
x <- 1
for (i in c("CC","CP","VD","VL")) {
  # Subset data for Colony 1 and perform ordination analysis
  col1 <- subset_samples(phy.t, treatment == i)
  ordi1 <- ordinate(col1, method = "RDA", distance = "euclidean")

  # Plot ordination results for Colony 1 with color based on Treatment2 variable, and add text labels
  plot1 <- plot_ordination(col1, ordi1, color = "timepoint") +
    geom_point(size = 3, alpha = 1) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_text_repel(aes(label = sample_clean)) +
    theme(legend.position = "left") +
    scale_color_manual(values = col_timepoint) +
    labs(title = paste("Treatment ", i, sep = ""))

  # Save the plot to the list
  plot_treatment_list[[x]] <- plot1
  x <- x + 1
  # Print the plot
  # print(plot1)
}

plot_time_list <- list()
x <- 1
# Create a list of colors for the treatments
for (i in c("T0", "T1", "T2")) {
  # Subset data for Colony 1 and perform ordination analysis
  col1 <- subset_samples(phy.t, timepoint == i)
  ordi1 <- ordinate(col1, method = "RDA", distance = "euclidean")

  # Define the shape type based on the treatment variable. this provide consistents shapes between differnt plots
  shape_type <- ifelse(i == "T0", 16,
                        ifelse(i == "T1", 17,
                              ifelse(i == "T2", 15, 13)))

  # Plot ordination results for Colony 1 with color based on Treatment2 variable, and add text labels
  plot1 <- plot_ordination(col1, ordi1, color = "treatment") +
    geom_point(size = 3, alpha = 1, shape = shape_type) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_text_repel(aes(label = sample_clean)) +
    theme(legend.position = "left") +
    scale_color_manual(values = col_treatment) +
    labs(title = paste("Timepoint ", i, sep = ""))

  # Save the plot to the list
  plot_time_list[[x]] <- plot1
  x <- x + 1
  # Print the plot
  # print(plot1)
}

# Activate httpgd graphics device
httpgd::hgd()
# Arrange and plot the figures to show all treatments by specific timepoint
ggarrange(plot_treatment_list[[1]],
          plot_treatment_list[[2]],
          plot_treatment_list[[3]],
          plot_treatment_list[[4]],
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")

# Arrange and plot the figures
ggarrange(plot_all,
      ggarrange(plot_time_list[[1]],
      plot_time_list[[2]],
      plot_time_list[[3]],
      ncol = 3, nrow = 1,
      common.legend = TRUE,
      legend = "none"),
      ncol = 1, nrow = 2,
      common.legend = TRUE,
      legend = "bottom")

# Save the combined plots as a PDF in landscape orientation
pdf(file = "/home/colinl/Proj/microbiome_probiotics_RNASeq/results/deseq2_R/plots/ordination_plot.pdf", width = 11, height = 8.5, pointsize = 6)

# Arrange and plot the figures
ggarrange(plot_all, 
      ggarrange(plot_time_list[[1]],
      plot_time_list[[2]],
      plot_time_list[[3]],
      ncol = 3, nrow = 1,
      common.legend = TRUE,
      legend = "none"),
      ncol = 1, nrow = 2,
      common.legend = TRUE,
      legend = "bottom")

# Close the PDF device
dev.off()
