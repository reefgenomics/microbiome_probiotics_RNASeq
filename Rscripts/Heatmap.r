#!/usr/bin/Rscript
library(pacman) # use pacman as a failsafe in case of missing packages.
p_load(readr,
        DESeq2,
        ellipse,
        ggplot2,
        dplyr,
        tidyr,
        data.table,
        tximport,
        apeglm,
        pheatmap)

# set output directory variable
out_dir <- "/home/colinl/Proj/microbiome_probiotics_RNASeq/results/deseq2_R/heatmap"
# Create output directory if it doesn't exist
if (!dir.exists(out_dir)) {
        dir.create(paste(out_dir, "plot", sep = "/"), recursive = TRUE)
        dir.create(paste(out_dir, "matrixes", sep = "/"), recursive = TRUE)
}

# Preparing count data from nf-core/rnaseq pipeline star_salmon output.
counts <- read.table("/home/colinl/Proj/microbiome_probiotics_RNASeq/results/star_salmon/salmon.merged.gene_counts.tsv", header = TRUE, sep = "\t", row.names = 1)
counts <- counts[, -1]
counts <- round(counts[, order(colnames(counts))], digits = 0)
# fix colnames names removing the X before the 2024. #sample set specific
counts <- counts %>% rename_all( ~ sub("X2024_Pver_", "", .))
colnames(counts)
# Convert counts to a matrix for expression analysis
expr_mat <- as.matrix(counts)
# Assume 'expr_mat' is your expression matrix
# Filter to keep genes with a sum above a threshold (e.g., 10)
filtered_expr <- expr_mat[rowSums(expr_mat) > 10, ]

# Scale the rows of the matrix
scaled_expr <- t(scale(t(filtered_expr)))
# pdf("/home/colinl/Proj/microbiome_probiotics_RNASeq/heatmap_expressed_genes.pdf")

# Save the heatmap to a PNG file
png(paste0(out_dir,"/plot/heatmap_expressed_genes.png"), width = 2000, height = 1500, res = 300)
pheatmap(scaled_expr,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        show_rownames = FALSE,    # Hide gene names if there are too many
        show_colnames = TRUE,
        main = "Heatmap of Expressed Genes")
dev.off()
# save the scaled expression matrix
write.table(scaled_expr, file = paste0(out_dir,"/matrixes/scaled_expression_matrix.txt"), sep = "\t", quote = FALSE, col.names = NA)

### --- subset by timepoint --- ###
# Subset the expression matrix by timepoint
# Define the timepoints
timepoints <- c("T0", "T1", "T2")
# Create a list to store the expression matrices for each timepoint
expr_matrices <- list()
# Loop through each timepoint and create a subset of the expression matrix
for (timepoint in timepoints) {
        # Create a pattern to match the timepoint
        pattern <- paste0("[CV][CPDL]_", timepoint, "_[0-9][0-9]_[0-9][0-9]")
        # Subset the expression matrix
        subset_counts <- counts[, grepl(pattern, colnames(counts))]
        length(colnames(subset_counts))

        expr_mat <- as.matrix(subset_counts)
        # Filter to keep genes with a sum above a threshold (e.g., 10)
        filtered_expr <- expr_mat[rowSums(expr_mat) > 10, ]

        # Scale the rows of the matrix
        scaled_expr <- t(scale(t(filtered_expr)))

        # Store the subset in the list
        expr_matrices[[timepoint]] <- scaled_expr
        # Save the heatmap for each timepoint
        png(paste0(out_dir,"/plot/heatmap_expressed_genes_timepoint_", timepoint, ".png"), width = 2000, height = 1500, res = 300)
        pheatmap(scaled_expr,
                color = colorRampPalette(c("blue", "white", "red"))(100),
                clustering_distance_rows = "euclidean",
                clustering_distance_cols = "euclidean",
                show_rownames = FALSE,
                show_colnames = TRUE,
                main = paste("Heatmap of Expressed Genes at Timepoint", timepoint))
        dev.off()
}

# save the scaled expression matrix for each timepoint
for (timepoint in timepoints) {
        write.table(expr_matrices[[timepoint]], file = paste0(out_dir, "/matrixes/scaled_expression_matrix_", timepoint, ".txt"), sep = "\t", quote = FALSE, col.names = NA)
}
