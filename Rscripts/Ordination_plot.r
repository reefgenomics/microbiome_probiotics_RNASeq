#!/usr/bin/Rscript

# Load required libraries
library(pacman)
p_load(phyloseq, microbiome, httpgd, ggrepel, ggplot2, ggpubr, dplyr, readr, DESeq2)

# # set output directory variable
out_dir <- "/home/colinl/Proj/microbiome_probiotics_RNASeq/results/deseq2_R/Ordination_plot"
# Create output directory if it doesn't exist
if (!dir.exists(out_dir)) {
        dir.create(paste(out_dir,"plot", "",  sep = "/"), recursive = TRUE)
}
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
shape_time <- c("T0" = 16, "T1" = 15, "T2" = 17)

# Plot overall ordination
plot_all <- plot_ordination(phy.t, ordi, color = "treatment", shape = "timepoint") +
  geom_point(size = 3, alpha = 1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(aes(label = sample_clean)) +
  theme(legend.position = "left") +
  scale_shape_manual(values = shape_time) +
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
  # scale_shape_manual(values = shape_time) +
  labs(title = paste("Treatment ", i, sep = ""))
  
  plot_treatment_list[[i]] <- plot1
}

# Generate plots for each timepoint
plot_time_list <- list()
for (i in c("T0", "T1", "T2")) {
  col1 <- subset_samples(phy.t, timepoint == i)
  ordi1 <- ordinate(col1, method = "RDA", distance = "euclidean")
  
  plot1 <- plot_ordination(col1, ordi1, color = "treatment", shape = "timepoint") +
  geom_point(size = 3, alpha = 1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(aes(label = sample_clean)) +
  theme(legend.position = "left") +
  scale_shape_manual(values = c("T0" = 16, "T1" = 15, "T2" = 17)) +
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
pdf(file = paste0(out_dir, "/plot/ordination_plot.pdf"),
  width = 11, height = 8.5, pointsize = 6)

  ggarrange(plot_all,
        ggarrange(plot_time_list[[1]], plot_time_list[[2]], plot_time_list[[3]],
            ncol = 3, nrow = 1, common.legend = TRUE, legend = "none"),
        ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")

  ggarrange(plot_treatment_list[[1]], plot_treatment_list[[2]],
        plot_treatment_list[[3]], plot_treatment_list[[4]],
      ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

dev.off()  # Close the PDF device

### --- New section for PCA plot --- ###
# load the dds objects
object_list <- list.files(path = "/home/colinl/Proj/microbiome_probiotics_RNASeq/results/deseq2_R/dds_objects/", pattern = "*.rds", full.names = TRUE)

# httpgd::hgd()
# Loop through each DESeq2 object and perform PCA analysis
pca_plot_list <- list()
for (i in c(1:length(object_list))) {
    comp_name <- gsub(".rds", "",gsub("DESeq_data_", "", basename(object_list[i])))
    dds_n <- readRDS(object_list[i])
    # universe <- rownames(dds_n)#[rowSums(counts(dds_n, normalized=TRUE) >= 10) >= 6]  
    
    vsd <- vst(dds_n)
    pca_plot <-  plotPCA(vsd, intgroup = c( "timepoint", "treatment")) +
      geom_point(aes(color = treatment, shape = timepoint), size = 4) +
      ggtitle(paste0("PCA DEseq  - ", comp_name)) +
      theme_bw() +
      scale_color_manual(values = col_treatment) + 
      scale_shape_manual(values = shape_time)
    
    pca_plot_list[[comp_name]] <- pca_plot

    # print(pca_plot)
}

#custom function to extract legend from ggplot object
extract_legend <- function(plot) {
    g <- ggplotGrob(plot)
    legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
    return(legend)
}

legend <- extract_legend(pca_plot_list[["C0"]] + theme(legend.position = "bottom"))

pdf(file = paste0(out_dir, "/plot/PCA_plot.pdf"),
      width = 11, height = 8.5, pointsize = 6)

# [1] "T0-ALL_TREATMENTS"
# [1] "T1-ALL_TREATMENTS"
# [1] "T2-ALL_TREATMENTS"
ggarrange(
  ggarrange(pca_plot_list[["C0"]] + ggtitle("PCA - DEseq2 all"), pca_plot_list[["T0-ALL_TREATMENTS"]],
          pca_plot_list[["T0-ALL_TREATMENTS"]], pca_plot_list[["T0-ALL_TREATMENTS"]],
            ncol = 2, nrow = 2, common.legend = TRUE, legend = "none"),
            legend,
            ncol = 1, nrow = 2, heights = c(2, 0.1))

# [1] "CC_T0-1-2"
# [1] "CP_T0-1-2"
# [1] "VD_T0-1-2"
# [1] "VL_T0-1-2"
ggarrange(
  ggarrange(pca_plot_list[["CC_T0-1-2"]], pca_plot_list[["CP_T0-1-2"]],
      pca_plot_list[["VD_T0-1-2"]], pca_plot_list[["VL_T0-1-2"]],
      ncol = 2, nrow = 2, common.legend = TRUE, legend = "none"),
      legend,
      ncol = 1, nrow = 2, heights = c(2, 0.1))

# [1] "CCvCP-T0"
# [1] "CCvVD-T0"
# [1] "CCvVL-T0"
# [1] "CPvVD-T0"
# [1] "CPvVL-T0"
# [1] "VDvVL-T0"

ggarrange(
  ggarrange(pca_plot_list[["CCvCP-T0"]], pca_plot_list[["CCvVD-T0"]],
        pca_plot_list[["CCvVL-T0"]], pca_plot_list[["CPvVD-T0"]],
        pca_plot_list[["CPvVL-T0"]], pca_plot_list[["VDvVL-T0"]],
        ncol = 2, nrow = 3, common.legend = TRUE, legend = "none"),
        legend,
        ncol = 1, nrow = 2, heights = c(2, 0.1))

# [1] "CCvCP-T1"
# [1] "CCvVD-T1"
# [1] "CCvVL-T1"
# [1] "CPvVD-T1"
# [1] "CPvVL-T1"
# [1] "VDvVL-T1"

ggarrange(
  ggarrange(pca_plot_list[["CCvCP-T1"]], pca_plot_list[["CCvVD-T1"]],
        pca_plot_list[["CCvVL-T1"]], pca_plot_list[["CPvVD-T1"]],
        pca_plot_list[["CPvVL-T1"]], pca_plot_list[["VDvVL-T1"]],
        ncol = 2, nrow = 3, common.legend = TRUE, legend = "none"),
        legend,
        ncol = 1, nrow = 2, heights = c(2, 0.1))

# [1] "CCvCP-T2"
# [1] "CCvVD-T2"
# [1] "CPvVL-T2"
# [1] "CPvVD-T2"
# [1] "CCvVL-T2"
# [1] "VDvVL-T2"

ggarrange(
  ggarrange(pca_plot_list[["CCvCP-T2"]], pca_plot_list[["CCvVD-T2"]],
        pca_plot_list[["CCvVL-T2"]], pca_plot_list[["CPvVD-T2"]],
        pca_plot_list[["CPvVL-T2"]], pca_plot_list[["VDvVL-T2"]],
        ncol = 2, nrow = 3, common.legend = TRUE, legend = "none"),
        legend,
        ncol = 1, nrow = 2, heights = c(2, 0.1))

# [1] "CCvCP_T0-1-2"
# [1] "CCvVD_T0-1-2"
# [1] "CCvVL_T0-1-2"
# [1] "CPvVD_T0-1-2"
# [1] "CPvVL_T0-1-2"
# [1] "VDvVL_T0-1-2"

ggarrange(
  ggarrange(pca_plot_list[["CCvCP_T0-1-2"]], pca_plot_list[["CCvVD_T0-1-2"]],
        pca_plot_list[["CCvVL_T0-1-2"]], pca_plot_list[["CPvVD_T0-1-2"]],
        pca_plot_list[["CPvVL_T0-1-2"]], pca_plot_list[["VDvVL_T0-1-2"]],
        ncol = 2, nrow = 3, common.legend = TRUE, legend = "none"),
        legend,
        ncol = 1, nrow = 2, heights = c(2, 0.1))

dev.off()  # Close the PDF device
