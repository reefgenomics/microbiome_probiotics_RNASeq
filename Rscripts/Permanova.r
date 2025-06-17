
# PERMANOVA analysis of RNA-Seq data using DESeq2-transformed counts.
# Loads sample metadata and DESeq2 objects, computes distance matrices,
# runs PERMANOVA and pairwise tests to assess effects of treatment and timepoint,
# and saves results to output files.

# Clean up the environment
rm(list = ls())

# Load required libraries
library(vegan)
library(pairwiseAdonis)
library(DESeq2)

# prepare metadata dataframe. Extracted form samplesheet.csv sample names in out case.
samples <- read.csv("/home/colinl/Proj/microbiome_probiotics_RNASeq/samplesheet.csv", header = TRUE)
samples <- samples$sample

#make sample dataframe with colnames from counts
metadata <- data.frame(
        sample = samples,
        treatment = factor(sub(".*_(CC|CP|VD|VL)_.*", "\\1", samples), # extract treatment from sample name
                        levels = c("CC", "CP", "VD", "VL")), # set levels
        timepoint = factor(sub(".*_T(\\d+)_.*", "T\\1", samples), # extract timepoint from sample name
                        levels = c("T0", "T1", "T2")), # set levels
        colony = factor(sub(".*_T\\d+_(\\d+)_\\d+", "\\1", samples)) # extract colony from sample name
)
rownames(metadata) <- metadata$sample # set sample as rownames

#set timepoint, treatment and colony as factor
metadata$timepoint <- as.factor(metadata$timepoint)
metadata$treatment <- as.factor(metadata$treatment)
metadata$colony <- as.factor(metadata$colony) # colony will not be used as factor in the design formula for DEseq2 but it is useful for the metadata

colnames(metadata)
metadata <- metadata[, c("treatment", "timepoint", "colony")] # keep only relevant columns

# Create DESeq2 object and run variance stabilizing transformation
object_list <- list.files(path = "/home/colinl/Proj/microbiome_probiotics_RNASeq/results/deseq2_R/dds_objects/", pattern = "*.rds", full.names = TRUE)
object_list <- object_list[order(nchar(object_list))] # sort by length of file name to clean up multi output text file later

# Load the first DESeq2 object to get the matrix
dds <- readRDS(object_list[1])
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

# Calculate distance matrix (e.g., Euclidean)
dist_mat <- dist(t(vsd_mat))

# Run PERMANOVA
# Since samples from the same colony may not be independent, it's sensible to use colony as a strata (blocking factor).
# This controls for the non-independence of samples from the same colony.
perm_control <- with(metadata, how(
    within = Within(type = "series"),   # preserves order of timepoints sequence
    plots  = Plots(strata = colony),    # treats each colony as a ‘plot’ of repeated measures
    blocks = colony,                    # no swaps between different colonies
    nperm  = 99999
))

# set.seed(123) # only for reproducibility
adonis_result <- adonis2(
            dist_mat ~ treatment + timepoint,
            data = metadata,
            permutations = perm_control,
            by = "margin"
        )
print(adonis_result)

# set.seed(123) # only for reproducibility
# Perform pairwise comparisons using pairwise.adonis2
results_pairwise <- pairwise.adonis2(
            dist_mat ~ treatment + timepoint,
            data = metadata,
            nperm = 99999,
            plots = Plots(strata = metadata$colony),
            strata = metadata$colony,
            p.adjust = "fdr",
            by = "margin"
        )
# explore with e.g.: "results_pairwise$VD_vs_VL"

# Save the results to a file
out_dir <- "/home/colinl/Proj/microbiome_probiotics_RNASeq/results/deseq2_R/permanova_results/joined-matrix/"
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}
# Save the adonis results to a text file
out_file <- paste0(out_dir, "permanova_results.txt")
cat("PERMANOVA Results:\n", file = out_file)
capture.output(print(adonis_result), file = out_file, append = TRUE)
# Save the pairwise results to a text file
out_file_pairwise <- paste0(out_dir, "pairwise_permanova_results.txt")
cat("Pairwise PERMANOVA Results:\n", file = out_file_pairwise)
capture.output(print(results_pairwise), file = out_file_pairwise, append = TRUE)

### run pairwise adonis2 for with blocks ###
## loop though each treatment pair and timepoint
tpair <- combn(levels(metadata$treatment), 2, simplify = FALSE)
for (i in tpair) {
    meta_tpair <- metadata[metadata$treatment %in% i, ]
    print(i)

    #Pairwise Permanova for each timepoint
    for (tp in unique(meta_tpair$timepoint)) {
        # Subset metadata for the current timepoint
        meta_tp <- meta_tpair[meta_tpair$timepoint == tp, ]
        
        # Subset distance matrix for the current timepoint
        vsd_mat_tp <- vsd_mat[, rownames(meta_tp)]
        dist_mat_tp <- dist(t(vsd_mat_tp))
        
        perm_control <- with(meta_tp, how(
        within = Within(type = "series"), # preserves order of timepoints sequence
        plots  = Plots(strata = colony),    # treats each colony as a ‘plot’ of repeated measures
        blocks = colony,                    # no swaps between different colonies
        nperm  = 999999
        ))

        # set.seed(123) # only for reproducibility
        pairwise_adonis_result_tp <- adonis2(
                    dist_mat_tp ~ treatment,
                    data = meta_tp,
                    permutations = perm_control,
                    by = "margin"
        )

        print(paste("Results for", i[1], "vs", i[2], "timepoint", tp))
        print(pairwise_adonis_result_tp)

    # Save the adonis results to a text file
    out_file <- paste0(out_dir, i[1], "v", i[2], "_permanova_results.txt")

    cat(paste("Results for", i[1], "vs", i[2], "timepoint", tp, ":\n"), file = out_file, append = TRUE)
    capture.output(print(pairwise_adonis_result_tp), file = out_file, append = TRUE)
    cat("\n", file = out_file, append = TRUE)
    }
}

### -- invert loop --- ### 
## loop though each timepoint pair and treatment
tpair <- combn(levels(metadata$timepoint), 2, simplify = FALSE)
for (i in tpair) {
    meta_tpair <- metadata[metadata$timepoint %in% i, ]
    print(i)

    for (tp in unique(meta_tpair$treatment)) {
    #Pairwise Permanova for each timepoint
        # Subset metadata for the current timepoint
        meta_tp <- meta_tpair[meta_tpair$treatment == tp, ]
        
        # Subset distance matrix for the current timepoint
        vsd_mat_tp <- vsd_mat[, rownames(meta_tp)]
        dist_mat_tp <- dist(t(vsd_mat_tp))
        
        perm_control <- with(meta_tp, how(
        within = Within(type = "series"), # preserves order of timepoints sequence
        plots  = Plots(strata = colony),    # treats each colony as a ‘plot’ of repeated measures
        blocks = colony,                    # no swaps between different colonies
        nperm  = 99999
        ))

        # set.seed(123) # only for reproducibility
        pairwise_adonis_result_tp <- adonis2(
                    dist_mat_tp ~ timepoint,
                    data = meta_tp,
                    permutations = perm_control,
                    by = "margin"
        )

        print(paste("Results for", i[1], "vs", i[2], "treatment", tp))
        print(pairwise_adonis_result_tp)

    # Save the adonis results to a text file
    out_file <- paste0(out_dir, i[1], "v", i[2], "_permanova_results.txt")

    cat(paste("Results for", i[1], "vs", i[2], "treatment", tp, ":\n"), file = out_file, append = TRUE)
    capture.output(print(pairwise_adonis_result_tp), file = out_file, append = TRUE)
    cat("\n", file = out_file, append = TRUE)
    }
}


### run pairwise adonis2 for without blocks ###
out_dir <- "/home/colinl/Proj/microbiome_probiotics_RNASeq/results/deseq2_R/permanova_results_no_blocks/joined-matrix/"
## loop though each treatment pair and timepoint
tpair <- combn(levels(metadata$treatment), 2, simplify = FALSE)
for (i in tpair) {
    meta_tpair <- metadata[metadata$treatment %in% i, ]
    print(i)

    #Pairwise Permanova for each timepoint
    for (tp in unique(meta_tpair$timepoint)) {
        # Subset metadata for the current timepoint
        meta_tp <- meta_tpair[meta_tpair$timepoint == tp, ]
        
        # Subset distance matrix for the current timepoint
        vsd_mat_tp <- vsd_mat[, rownames(meta_tp)]
        dist_mat_tp <- dist(t(vsd_mat_tp))
        
        perm_control <- with(meta_tp, how(
        # within = Within(type = "series"), # preserves order of timepoints sequence
        # plots  = Plots(strata = colony),    # treats each colony as a ‘plot’ of repeated measures
        blocks = colony,                    # no swaps between different colonies
        nperm  = 999999
        ))

        # set.seed(123) # only for reproducibility
        pairwise_adonis_result_tp <- adonis2(
                    dist_mat_tp ~ treatment,
                    data = meta_tp,
                    permutations = perm_control,
                    by = "margin"
        )

        print(paste("Results for", i[1], "vs", i[2], "timepoint", tp))
        print(pairwise_adonis_result_tp)

    # Save the adonis results to a text file
    out_file <- paste0(out_dir, i[1], "v", i[2], "_permanova_results.txt")

    cat(paste("Results for", i[1], "vs", i[2], "timepoint", tp, ":\n"), file = out_file, append = TRUE)
    capture.output(print(pairwise_adonis_result_tp), file = out_file, append = TRUE)
    cat("\n", file = out_file, append = TRUE)
    }
}

### -- invert loop --- ###
## loop though each timepoint pair and treatment 
tpair <- combn(levels(metadata$timepoint), 2, simplify = FALSE)

for (i in tpair) {
    meta_tpair <- metadata[metadata$timepoint %in% i, ]
    print(i)

    for (tp in unique(meta_tpair$treatment)) {
    #Pairwise Permanova for each timepoint
        # Subset metadata for the current timepoint
        meta_tp <- meta_tpair[meta_tpair$treatment == tp, ]
        
        # Subset distance matrix for the current timepoint
        vsd_mat_tp <- vsd_mat[, rownames(meta_tp)]
        dist_mat_tp <- dist(t(vsd_mat_tp))
        
        perm_control <- with(meta_tp, how(
        # within = Within(type = "series"), # preserves order of timepoints sequence
        # plots  = Plots(strata = colony),    # treats each colony as a ‘plot’ of repeated measures
        # blocks = colony,                    # no swaps between different colonies
        nperm  = 99999
        ))

        # set.seed(123) # only for reproducibility
        pairwise_adonis_result_tp <- adonis2(
                    dist_mat_tp ~ timepoint,
                    data = meta_tp,
                    permutations = perm_control,
                    by = "margin"
        )

        print(paste("Results for", i[1], "vs", i[2], "treatment", tp))
        print(pairwise_adonis_result_tp)

    # Save the adonis results to a text file
    out_file <- paste0(out_dir, i[1], "v", i[2], "_permanova_results.txt")

    cat(paste("Results for", i[1], "vs", i[2], "treatment", tp, ":\n"), file = out_file, append = TRUE)
    capture.output(print(pairwise_adonis_result_tp), file = out_file, append = TRUE)
    cat("\n", file = out_file, append = TRUE)
    }
}

# Clean up object_list to remove objects with only one time point and 2 treatments.
# This is a filter to exclude objects that would be blocked by strata in the adonis test. # activate if you want to filter out objects with low number of permutations. 
# object_list_cut <- object_list[!grepl("[CV][CPDL]v[CV][CPDL]-T[012]\\.rds$", object_list)]
object_list_cut <- object_list # keep all objects for now

out_dir <- "/home/colinl/Proj/microbiome_probiotics_RNASeq/results/deseq2_R/permanova_results/individual-matrix/"
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

# Loop through each DESeq2 object in object_list_cut and perform PERMANOVA
for (i in c(2:length(object_list_cut))) {
    comp_name <- gsub(".rds", "",gsub("DESeq_data_", "", basename(object_list_cut[i])))
    dds_n <- readRDS(object_list_cut[i])

    vsd_n <- vst(dds_n, blind = TRUE)
    vsd_mat_n <- assay(vsd_n)

    dist_mat_n <- dist(t(vsd_mat_n))

    # metadata_n <- metadata[metadata$sample %in% colnames(vsd_mat_n), ]
    metadata_n <- metadata[rownames(metadata) %in% colnames(vsd_mat_n), ]

    perm_control <- with(metadata_n, how(
        within = Within(type = "series"),   # preserves order of timepoints sequence
        plots  = Plots(strata = colony),    # treats each colony as a ‘plot’ of repeated measures
        blocks = colony,                    # no swaps between different colonies
        nperm  = 99999
        ))

    result <- adonis2(dist_mat_n ~ treatment + timepoint, data = metadata_n,
                    permutations = perm_control, by = "margin")

    print(paste("Results for:", comp_name))
    print(result)

    # Generate grouping text based on the comparison name for output file naming
    if (grepl("^[CV][CPDL]_T0-1-2$", comp_name)) {
        grouping_text <- "ALL_TIMEPOINTS"
    } else if (grepl("^T[012]-ALL_TREATMENTS", comp_name)) {
        grouping_text <- "ALL_TREATMENTS"
    } else {
        grouping_text <- sub("([A-Za-z0-9]+v[A-Za-z0-9]+).*", "\\1", comp_name)
    }

    # Save the adonis results to a text file
    out_file <- paste0(out_dir, grouping_text, "_permanova_results.txt")

    cat(paste("comparison", comp_name, ":\n"), file = out_file, append = TRUE)
    capture.output(print(result), file = out_file, append = TRUE)
    cat("\n", file = out_file, append = TRUE)
}
