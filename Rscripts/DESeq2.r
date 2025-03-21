#!/usr/bin/Rscript
# DESeq2 analysis for the P.verrucosa RNASeq data
library(pacman) # use pacman as a failsafe in case of missing packages.
p_load(readr,
        DESeq2,
        ellipse,
        ggplot2,
        dplyr,
        tidyr,
        data.table,
        tximport,
        apeglm)

# set output directory variable
out_dir <- "/home/colinl/Proj/microbiome_probiotics_RNASeq/results/deseq2_R"
# Create output directory if it doesn't exist
if (!dir.exists(out_dir)) {
        dir.create(paste(out_dir,"dds_objects", sep = "/"), recursive = TRUE)
        dir.create(paste(out_dir,"summary", sep = "/"), recursive = TRUE)
}

# Preparing count data from nf-core/rnaseq pipeline star_salmon output.
counts <- read.table("/home/colinl/Proj/microbiome_probiotics_RNASeq/results/star_salmon/salmon.merged.gene_counts.tsv", header = TRUE, sep = "\t", row.names = 1)
counts <- counts[, -1]
counts <- round(counts[, order(colnames(counts))], digits = 0)
# fix colnames names removing the X before the 2024. #sample set specific
counts <- counts %>% rename_all( ~ sub("X", "", .))

# used to be sure gene names are consitent between annotation and counts
tx2gene <- read_tsv("/home/colinl/Proj/microbiome_probiotics_RNASeq/results/star_salmon/tx2gene.tsv")
# clean up improper formatting of pseudogenes due to agat conversion
# When the gene_id start with "agat-pseudogene-" invert gene_id and transcript_id columns. 
tx2gene <- tx2gene %>%
        mutate(gene_id = ifelse(grepl("agat-pseudogene-", gene_id), transcript_id, gene_id),
                transcript_id = ifelse(grepl("agat-pseudogene-", gene_name), gene_name, transcript_id))
# remove "gene-" prefix from gene_id
tx2gene$gene_id <- sub("gene-", "", tx2gene$gene_id)
# fix counts rownames matching tx2gene gene_name to the colname and writing gene_id
rownames(counts) <- tx2gene$gene_id[match(rownames(counts), tx2gene$gene_name)]

# Read in annotation data from NCBI 
all_annotations <- read_tsv("/home/colinl/Proj/microbiome_probiotics_RNASeq/Pver_ref/NCBIFilesDir/203993.gene.tsv")

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


# quick funcitons to print and save the number of DE genes
echo_message <- function(comp_title, Cn_res) {
        message(comp_title,
                "\nNumber of DE genes :", nrow(Cn_res),
                "\nNumber of up-regulated genes :", nrow(subset(Cn_res, log2FoldChange > 0)),
                "\nNumber of down-regulated genes: ", nrow(subset(Cn_res, log2FoldChange < 0)))
}
# echo_message(comp_title, Cn_res)

export_echo_message <- function(comp_title, log_file, Cn_res) {
        cat(paste(comp_title,
                "\nNumber of DE genes :", nrow(Cn_res),
                "\nNumber of up-regulated genes :", nrow(subset(Cn_res, log2FoldChange > 0)),
                "\nNumber of down-regulated genes: ", nrow(subset(Cn_res, log2FoldChange < 0)), "\n", sep = ""),
                file = paste(out_dir, "/summary/", log_file, sep = ""), append = TRUE)
}

# Function to write DEseq summary to dataframe
dds_counts_to_df <- function(comp_title, Cn_res) {
        data.frame(
                Comparison = gsub(" ", "_", sub("Comparison ", "", sub("for ", "", comp_title))),
                DE_genes = nrow(Cn_res),
                Up_regulated = nrow(subset(Cn_res, log2FoldChange > 0)),
                Down_regulated = nrow(subset(Cn_res, log2FoldChange < 0))
        )
}

# # Example usage
# dds_counts_to_df(comp_title, Cn_res)
# print(echo_df)

# echo_message(comp_title, Cn_res)

###––– DESeq2 analysis –––###
#comparison 0 all the samples all treatments vs all timepoints
C0 <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = timepoint ~ treatment) # set object with design formula
C0$treatment <- relevel(C0$treatment, ref = "CC") # set reference level for treatment
C0 <- DESeq(C0) # run DESeq2
C0_res <- as.data.frame(results(C0, # get results
                                lfcThreshold = 0,
                                alpha=0.05,
                                cooksCutoff=FALSE,
                                independentFiltering=FALSE)) %>%
                        filter(padj <= 0.05) # filter results for p-value adjusted of <= 0.05

C0_fin <- merge(C0_res[,c(2,6)], all_annotations, by.x="row.names", by.y = "Symbol", all.x=T) # merge results with annotation
colnames(C0_fin)[1] <- "Symbol" # correct column name
C0_fin$Comparison <- "Comparison 0" # add comparison column

# save DESeq(Cn) object to file # deactivate by default. Activate if re-running the DESeq() function becomes too slow
saveRDS(C0, paste(out_dir, "/dds_objects/DESeq_data_C0.rds", sep = ""))

# write message to console 
echo_message("Comparison 0 - All samples", C0_res)

# write message to text file
export_echo_message("Comparison 0 - All samples", "DESeq_summary_info.txt", C0_res)

# #write results to table
write.table(C0_fin, paste(out_dir, "/DESeq_Comparison-ALL.tsv", sep = ""), row.names = F, quote = F, sep = "\t")


###––– Other comparisons, same logic as above but within loops –––### 
#loop through all the Timepoints (T0, T1, T2)
rm(summary_dds_df) # remove summary_dds_df if it exists
for (i in c(0,1,2)) {
        countData <- (counts[,grepl(paste("2024_Pver_[CV][PCDL]_T",i,"_0[0-9]_0[1-9]", sep = ""), colnames(counts))]) # subset counts and metadata for each timepoint and count using regex
        sampleTable <- subset(metadata, rownames(metadata) %like% paste("2024_Pver_[CV][PCDL]_T",i,"_0[0-9]_0[1-9]", sep = "" ))

        comp_title <- paste("Comparison T", i, " all samples", sep = "")
        Cn <- DESeqDataSetFromMatrix(countData = countData, colData =  sampleTable, design = colony ~ treatment) 
        Cn$treatment <- relevel(Cn$treatment, ref = "CC")
        Cn <- DESeq(Cn)

        Cn_res <- as.data.frame(results(Cn, lfcThreshold = 0, alpha=0.05, cooksCutoff = FALSE, independentFiltering = FALSE)) %>% filter(padj <= 0.05)
        Cn_fin <- merge(Cn_res[,c(2,6)], all_annotations, by.x="row.names", by.y = "Symbol", all.x=T)
        colnames(Cn_fin)[1] <- "Symbol"
        Cn_fin$Comparison <- comp_title
        echo_message(comp_title, Cn_res)

        if (exists("summary_dds_df")) {
                summary_dds_df <- rbind(summary_dds_df, dds_counts_to_df(comp_title, Cn_res))
        } else {
                summary_dds_df <- dds_counts_to_df(comp_title, Cn_res)
        }
        # write summary to file
        write.table(summary_dds_df, paste(out_dir, "/summary/DESeq_summary_timepoint", ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
        # save DESeq(Cn) object to file # deactivate by default. Activate if re-running the DESeq() function becomes too slow
        saveRDS(Cn, paste(out_dir, "/dds_objects/DESeq_data_T",i,".rds", sep = ""))
        export_echo_message(comp_title, "DESeq_timepoints_summary_info.txt", Cn_res)
        write.table(Cn_fin, paste(out_dir, "/DESeq_Comparison-", comp_title, ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
}

# loop through all timpoints and treatment again CC. (CC vs CP , VD, VL) for T0,1,2.
rm(summary_dds_df) # remove summary_dds_df if it exists
for (j in c(0:2)) {
        for (i in c("CP", "VD", "VL")) {
                countData <- (counts[,grepl(paste("2024_Pver_CC_T",j,"_0[0-9]_0[1-9]|2024_Pver_",i,"_T",j,"_0[0-9]_0[1-9]", sep = ""), colnames(counts))])
                sampleTable <- subset(metadata, rownames(metadata) %like% paste("2024_Pver_CC_T",j,"_0[0-9]_0[1-9]|2024_Pver_",i,"_T",j,"_0[0-9]_0[1-9]", sep = ""))

                comp_title <- paste("Comparison CC vs ",i," for T",j, sep = "")
                Cn <- DESeqDataSetFromMatrix(countData = countData, colData =  sampleTable, design = colony ~ treatment) 
                Cn$treatment <- relevel(Cn$treatment, ref = "CC")
                Cn <- DESeq(Cn)

                Cn_res <- as.data.frame(results(Cn, lfcThreshold = 0, alpha=0.05, cooksCutoff = FALSE, independentFiltering = FALSE)) %>% filter(padj <= 0.05)
                Cn_fin <- merge(Cn_res[,c(2,6)], all_annotations, by.x="row.names", by.y = "Symbol", all.x=T)
                colnames(Cn_fin)[1] <- "Symbol"
                Cn_fin$Comparison <- comp_title
                echo_message(comp_title, Cn_res)

                if (exists("summary_dds_df")) {
                        summary_dds_df <- rbind(summary_dds_df, dds_counts_to_df(comp_title, Cn_res))
                } else {
                        summary_dds_df <- dds_counts_to_df(comp_title, Cn_res)
                }

                # write summary to file
                write.table(summary_dds_df, paste(out_dir, "/summary/DESeq_summary_Treatment", ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
                # save DESeq(Cn) object to file # deactivate by default. Activate if re-running the DESeq() function becomes too slow
                saveRDS(Cn, paste(out_dir, "/dds_objects/DESeq_data_CCvs",i,"_T",j,".rds", sep = ""))
                export_echo_message(comp_title, "DESeq_individual_comp_summary_info.txt", Cn_res)
                write.table(Cn_fin, paste(out_dir, "/DESeq_Comparison-", comp_title, ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
        }
}
# loop through all timpoints and treatment against CP (CP vs VD, VL) for T0,1,2.
# rm(summary_dds_df) # remove summary_dds_df if it exists
for (j in c(0:2)) {
        for (i in c("VD", "VL")) {
                countData <- (counts[,grepl(paste("2024_Pver_CP_T",j,"_0[0-9]_0[1-9]|2024_Pver_",i,"_T",j,"_0[0-9]_0[1-9]", sep = ""), colnames(counts))])
                sampleTable <- subset(metadata, rownames(metadata) %like% paste("2024_Pver_CP_T",j,"_0[0-9]_0[1-9]|2024_Pver_",i,"_T",j,"_0[0-9]_0[1-9]", sep = ""))

                comp_title <- paste("Comparison CP vs ",i," for T",j, sep = "")
                Cn <- DESeqDataSetFromMatrix(countData = countData, colData =  sampleTable, design = colony ~ treatment) 
                Cn$treatment <- relevel(Cn$treatment, ref = "CP")
                Cn <- DESeq(Cn)

                Cn_res <- as.data.frame(results(Cn, lfcThreshold = 0, alpha=0.05, cooksCutoff = FALSE, independentFiltering = FALSE)) %>% filter(padj <= 0.05)
                Cn_fin <- merge(Cn_res[,c(2,6)], all_annotations, by.x="row.names", by.y = "Symbol", all.x=T)
                colnames(Cn_fin)[1] <- "Symbol"
                Cn_fin$Comparison <- comp_title
                echo_message(comp_title, Cn_res)

                if (exists("summary_dds_df")) {
                        summary_dds_df <- rbind(summary_dds_df, dds_counts_to_df(comp_title, Cn_res))
                } else {
                        summary_dds_df <- dds_counts_to_df(comp_title, Cn_res)
                }

                # write summary to file
                write.table(summary_dds_df, paste(out_dir, "/summary/DESeq_summary_Treatment", ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
                # save DESeq(Cn) object to file # deactivate by default. Activate if re-running the DESeq() function becomes too slow
                saveRDS(Cn, paste(out_dir, "/dds_objects/DESeq_data_CPvs",i,"_T",j,".rds", sep = ""))
                export_echo_message(comp_title, "DESeq_individual_comp_summary_info.txt", Cn_res)
                write.table(Cn_fin, paste(out_dir, "/DESeq_Comparison-", comp_title, ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
        }
}

# loop through all timpoints and treatment: VD vs VL (for T0,1,2.)
# rm(summary_dds_df) # remove summary_dds_df if it exists
for (j in c(0:2)) {
        for (i in c("VL")) {
                countData <- (counts[,grepl(paste("2024_Pver_VD_T",j,"_0[0-9]_0[1-9]|2024_Pver_",i,"_T",j,"_0[0-9]_0[1-9]", sep = ""), colnames(counts))])
                sampleTable <- subset(metadata, rownames(metadata) %like% paste("2024_Pver_VD_T",j,"_0[0-9]_0[1-9]|2024_Pver_",i,"_T",j,"_0[0-9]_0[1-9]", sep = ""))

                comp_title <- paste("Comparison VD vs ",i," for T",j, sep = "")
                Cn <- DESeqDataSetFromMatrix(countData = countData, colData =  sampleTable, design = colony ~ treatment) 
                Cn$treatment <- relevel(Cn$treatment, ref = "VD")
                Cn <- DESeq(Cn)

                Cn_res <- as.data.frame(results(Cn, lfcThreshold = 0, alpha=0.05, cooksCutoff = FALSE, independentFiltering = FALSE)) %>% filter(padj <= 0.05)
                Cn_fin <- merge(Cn_res[,c(2,6)], all_annotations, by.x="row.names", by.y = "Symbol", all.x=T)
                colnames(Cn_fin)[1] <- "Symbol"
                Cn_fin$Comparison <- comp_title
                echo_message(comp_title, Cn_res)

                if (exists("summary_dds_df")) {
                        summary_dds_df <- rbind(summary_dds_df, dds_counts_to_df(comp_title, Cn_res))
                } else {
                        summary_dds_df <- dds_counts_to_df(comp_title, Cn_res)
                }

                # write summary to file
                write.table(summary_dds_df, paste(out_dir, "/summary/DESeq_summary_Treatment", ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
                # save DESeq(Cn) object to file # deactivate by default. Activate if re-running the DESeq() function becomes too slow
                saveRDS(Cn, paste(out_dir, "/dds_objects/DESeq_data_VDvs",i,"_T",j,".rds", sep = ""))
                export_echo_message(comp_title, "DESeq_individual_comp_summary_info.txt", Cn_res)
                write.table(Cn_fin, paste(out_dir, "/DESeq_Comparison-", comp_title, ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
        }
}

## DEseq2 within treatment at different timepoints
rm(summary_dds_df) # remove summary_dds_df if it exists
for (i in c("CC", "CP", "VD", "VL")) {
        # print(subset(metadata, rownames(metadata) %like% paste("2024_Pver_",i,"_T[0-2]_0[0-9]_0[1-9]", sep = "")))
        countData <- (counts[,grepl(paste("2024_Pver_",i,"_T[0-2]_0[0-9]_0[1-9]", sep = ""), colnames(counts))])
        sampleTable <- subset(metadata, rownames(metadata) %like% paste("2024_Pver_",i,"_T[0-2]_0[0-9]_0[1-9]", sep = ""))

        comp_title <- paste("Comparison for ", i," at T0-1-2", sep = "")
        Cn <- DESeqDataSetFromMatrix(countData = countData, colData =  sampleTable, design = ~ timepoint) 
        # Cn$treatment <- relevel(Cn$treatment, ref = "VD")
        Cn <- DESeq(Cn)

        Cn_res <- as.data.frame(results(Cn, lfcThreshold = 0, alpha=0.05, cooksCutoff = FALSE, independentFiltering = FALSE)) %>% filter(padj <= 0.05)
        Cn_fin <- merge(Cn_res[,c(2,6)], all_annotations, by.x="row.names", by.y = "Symbol", all.x=T)
        colnames(Cn_fin)[1] <- "Symbol"
        Cn_fin$Comparison <- comp_title
        echo_message(comp_title, Cn_res)

        if (exists("summary_dds_df")) {
                summary_dds_df <- rbind(summary_dds_df, dds_counts_to_df(comp_title, Cn_res))
        } else {
                summary_dds_df <- dds_counts_to_df(comp_title, Cn_res)
        }

        # write summary to file
        write.table(summary_dds_df, paste(out_dir, "/summary/DESeq_summary_Time", ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
        # save DESeq(Cn) object to file # deactivate by default. Activate if re-running the DESeq() function becomes too slow
        saveRDS(Cn, paste(out_dir, "/dds_objects/DESeq_data_",i,".rds", sep = ""))
        export_echo_message(comp_title, "DESeq_individual_comp_summary_info.txt", Cn_res)
        write.table(Cn_fin, paste(out_dir, "/DESeq_Comparison-", comp_title, ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
}


# loop through all treatment combinations again CC. (CC vs CP , VD, VL) for T0-1-2.
rm(summary_dds_df) # remove summary_dds_df if it exists
# j= "[0-2]"
for (j in c("[0-2]")) {
        for (i in c("CP", "VD", "VL")) {
                countData <- (counts[,grepl(paste("2024_Pver_CC_T",j,"_0[0-9]_0[1-9]|2024_Pver_",i,"_T",j,"_0[0-9]_0[1-9]", sep = ""), colnames(counts))])
                sampleTable <- subset(metadata, rownames(metadata) %like% paste("2024_Pver_CC_T",j,"_0[0-9]_0[1-9]|2024_Pver_",i,"_T",j,"_0[0-9]_0[1-9]", sep = ""))

                comp_title <- paste("Comparison CP vs ",i," for T", gsub("\\[0-2\\]", "0-1-2", j), sep = "")
                Cn <- DESeqDataSetFromMatrix(countData = countData, colData =  sampleTable, design = colony ~ treatment) 
                Cn$treatment <- relevel(Cn$treatment, ref = "CC")
                Cn <- DESeq(Cn)

                Cn_res <- as.data.frame(results(Cn, lfcThreshold = 0, alpha=0.05, cooksCutoff = FALSE, independentFiltering = FALSE)) %>% filter(padj <= 0.05)
                Cn_fin <- merge(Cn_res[,c(2,6)], all_annotations, by.x="row.names", by.y = "Symbol", all.x=T)
                colnames(Cn_fin)[1] <- "Symbol"
                Cn_fin$Comparison <- comp_title
                echo_message(comp_title, Cn_res)

                if (exists("summary_dds_df")) {
                        summary_dds_df <- rbind(summary_dds_df, dds_counts_to_df(comp_title, Cn_res))
                } else {
                        summary_dds_df <- dds_counts_to_df(comp_title, Cn_res)
                }

                # write summary to file
                write.table(summary_dds_df, paste(out_dir, "/summary/DESeq_summary_treatm_combination_times", ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
                # save DESeq(Cn) object to file # deactivate by default. Activate if re-running the DESeq() function becomes too slow
                saveRDS(Cn, paste(out_dir, "/dds_objects/DESeq_data_CCvs",i,"_T",j,".rds", sep = ""))
                export_echo_message(comp_title, "DESeq_individual_comp_summary_info.txt", Cn_res)
                write.table(Cn_fin, paste(out_dir, "/DESeq_Comparison-", comp_title, ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
        }
}
# loop through all timpoints and treatment against CP (CP vs VD, VL) for T0,1,2.
# rm(summary_dds_df) # remove summary_dds_df if it exists
for (j in c("[0-2]")) {
        for (i in c("VD", "VL")) {
                countData <- (counts[,grepl(paste("2024_Pver_CP_T",j,"_0[0-9]_0[1-9]|2024_Pver_",i,"_T",j,"_0[0-9]_0[1-9]", sep = ""), colnames(counts))])
                sampleTable <- subset(metadata, rownames(metadata) %like% paste("2024_Pver_CP_T",j,"_0[0-9]_0[1-9]|2024_Pver_",i,"_T",j,"_0[0-9]_0[1-9]", sep = ""))

                comp_title <- paste("Comparison CP vs ",i," for T", gsub("\\[0-2\\]", "0-1-2", j), sep = "")
                Cn <- DESeqDataSetFromMatrix(countData = countData, colData =  sampleTable, design = colony ~ treatment) 
                Cn$treatment <- relevel(Cn$treatment, ref = "CP")
                Cn <- DESeq(Cn)

                Cn_res <- as.data.frame(results(Cn, lfcThreshold = 0, alpha=0.05, cooksCutoff = FALSE, independentFiltering = FALSE)) %>% filter(padj <= 0.05)
                Cn_fin <- merge(Cn_res[,c(2,6)], all_annotations, by.x="row.names", by.y = "Symbol", all.x=T)
                colnames(Cn_fin)[1] <- "Symbol"
                Cn_fin$Comparison <- comp_title
                echo_message(comp_title, Cn_res)

                if (exists("summary_dds_df")) {
                        summary_dds_df <- rbind(summary_dds_df, dds_counts_to_df(comp_title, Cn_res))
                } else {
                        summary_dds_df <- dds_counts_to_df(comp_title, Cn_res)
                }

                # write summary to file
                write.table(summary_dds_df, paste(out_dir, "/summary/DESeq_summary_treatm_combination_times", ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
                # save DESeq(Cn) object to file # deactivate by default. Activate if re-running the DESeq() function becomes too slow
                saveRDS(Cn, paste(out_dir, "/dds_objects/DESeq_data_CPvs",i,"_T",j,".rds", sep = ""))
                export_echo_message(comp_title, "DESeq_individual_comp_summary_info.txt", Cn_res)
                write.table(Cn_fin, paste(out_dir, "/DESeq_Comparison-", comp_title, ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
        }
}

# loop through all timpoints and treatment: VD vs VL (for T0,1,2.)
# rm(summary_dds_df) # remove summary_dds_df if it exists
for (j in c("[0-2]")) {
        for (i in c("VL")) {
                countData <- (counts[,grepl(paste("2024_Pver_VD_T",j,"_0[0-9]_0[1-9]|2024_Pver_",i,"_T",j,"_0[0-9]_0[1-9]", sep = ""), colnames(counts))])
                sampleTable <- subset(metadata, rownames(metadata) %like% paste("2024_Pver_VD_T",j,"_0[0-9]_0[1-9]|2024_Pver_",i,"_T",j,"_0[0-9]_0[1-9]", sep = ""))

                comp_title <- paste("Comparison CP vs ",i," for T", gsub("\\[0-2\\]", "0-1-2", j), sep = "")
                Cn <- DESeqDataSetFromMatrix(countData = countData, colData =  sampleTable, design = colony ~ treatment) 
                Cn$treatment <- relevel(Cn$treatment, ref = "VD")
                Cn <- DESeq(Cn)

                Cn_res <- as.data.frame(results(Cn, lfcThreshold = 0, alpha=0.05, cooksCutoff = FALSE, independentFiltering = FALSE)) %>% filter(padj <= 0.05)
                Cn_fin <- merge(Cn_res[,c(2,6)], all_annotations, by.x="row.names", by.y = "Symbol", all.x=T)
                colnames(Cn_fin)[1] <- "Symbol"
                Cn_fin$Comparison <- comp_title
                echo_message(comp_title, Cn_res)

                if (exists("summary_dds_df")) {
                        summary_dds_df <- rbind(summary_dds_df, dds_counts_to_df(comp_title, Cn_res))
                } else {
                        summary_dds_df <- dds_counts_to_df(comp_title, Cn_res)
                }

                # write summary to file
                write.table(summary_dds_df, paste(out_dir, "/summary/DESeq_summary_treatm_combination_times", ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
                # save DESeq(Cn) object to file # deactivate by default. Activate if re-running the DESeq() function becomes too slow
                saveRDS(Cn, paste(out_dir, "/dds_objects/DESeq_data_VDvs",i,"_T",j,".rds", sep = ""))
                export_echo_message(comp_title, "DESeq_individual_comp_summary_info.txt", Cn_res)
                write.table(Cn_fin, paste(out_dir, "/DESeq_Comparison-", comp_title, ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
        }
}
