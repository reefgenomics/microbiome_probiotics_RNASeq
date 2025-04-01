
library(pacman)
p_load(clusterProfiler, AnnotationDbi, GenomicFeatures, AnnotationForge, org.Pverrucosa.eg.db, enrichplot, clusterProfiler, ggplot2, dplyr, tidyr, data.table, tximport, apeglm, topGO, deseq2)
# p_load(org.Hs.eg.db) # Substitute with coral-specific DB if available
library(org.Pverrucosa.eg.db) # See Pver_ref/README.md file for more info  # install.packages("/home/colinl/Proj/microbiome_probiotics_RNASeq/Pver_ref/org.Pverrucosa.eg.db", repos = NULL) 


# # set output directory variable
out_dir <- "/home/colinl/Proj/microbiome_probiotics_RNASeq/results/deseq2_R/GO_enrichment"
# Create output directory if it doesn't exist
if (!dir.exists(out_dir)) {
        dir.create(paste(out_dir,"plot", "cnetplot",  sep = "/"), recursive = TRUE)
        dir.create(paste(out_dir,"plot", "dotplot",  sep = "/"), recursive = TRUE)
        dir.create(paste(out_dir,"summary","fisher", sep = "/"), recursive = TRUE)
}

gff_file <- "/home/colinl/Proj/microbiome_probiotics_RNASeq/Pver_ref/genomic.simple.gtf.gz"

txdb <- makeTxDbFromGFF(gff_file, format = "gtf")
gene_map <- genes(txdb)
valid_columns <- columns(txdb)
# print(valid_columns) # Check for valid column names
gene_map$GENENAME <- mapIds(txdb, keys = gene_map$gene_id, 
                                            column = "GENEID", keytype = "GENEID", 
                                            multiVals = "first")
gene_map$GENENAME <- sub("gene-", "", gene_map$GENENAME)
gene_map$gene_id <- sub("gene-", "", gene_map$gene_id)

# load the dds objects
object_list <- list.files(path = "/home/colinl/Proj/microbiome_probiotics_RNASeq/results/deseq2_R/dds_objects/", pattern = "*.rds", full.names = TRUE)

httpgd::hgd()
# Loop through each DESeq2 object and perform individual GO enrichment analysis
for (i in c(1:length(object_list))) {
    comp_name <- gsub(".rds", "",gsub("DESeq_data_", "", basename(object_list[i])))
    dds_n <- readRDS(object_list[i])
    universe <- rownames(dds_n)#[rowSums(counts(dds_n, normalized=TRUE) >= 10) >= 6]

    res_dds_n <- results(dds_n,lfcThreshold = 0,alpha=0.05,cooksCutoff=FALSE,independentFiltering=FALSE)
    res_df <- data.frame(
      gene = gene_map$GENENAME[match(rownames(res_dds_n), gene_map$gene_id)],
      log2FC = res_dds_n$log2FoldChange,
      padj = res_dds_n$padj
    ) |> na.omit()

    sig_genes <- res_df$gene[res_df$padj < 0.05 & abs(res_df$log2FC) > 1]

    # GO Enrichment for VL vs VD
    ego_dds_n <- enrichGO(
      gene = sig_genes,
      universe = universe,
      OrgDb = org.Pverrucosa.eg.db, # Use coral annotation if possible
      keyType = "SYMBOL",
      ont = "BP", # Biological Process
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05
    )

    if (nrow(as.data.frame(ego_dds_n)) == 0) {
      message(paste("Error in processing", comp_name, ":", "no significant GO terms found."))
      cat(paste("Error in processing", comp_name, ":", "no significant GO terms found.", "\n"),
          file = paste0(out_dir, "/summary/GO_enrichment.log"), append = TRUE)
      next
    }
    # experiments # make object with ego result if not empty.
    assign(paste0("ego_", comp_name), ego_dds_n)
    message(paste("ego_", comp_name, " made!", sep = ""))
    cat(paste("ego_", comp_name, " made!", "\n", sep = ""), file = paste0(out_dir, "/summary/GO_enrichment.log"), append = TRUE)

      pdf(file = paste0(out_dir, "/plot/GO_enrichment_dotplot_", comp_name, ".pdf"))
      print(dotplot(ego_dds_n, showCategory = 20) + ggtitle(paste("GO Enrichment (", comp_name, ")", sep = "")))
      dev.off()

      pdf(file = paste0(out_dir, "/plot/GO_enrichment_cnetplot_", comp_name, ".pdf"))
      print(cnetplot(ego_dds_n, showCategory = 20) + ggtitle(paste("GO Enrichment (", comp_name, ")", sep = "")))
      dev.off()
      write.table(as.data.frame(ego_dds_n), file = paste0(out_dir, "/summary/GO_enrichment_", comp_name, ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}

#check the names of the dds object previously created
for (i in seq_along(object_list)) {
  comp_name <- gsub(".rds", "", gsub("DESeq_data_", "", basename(object_list[i])))
  print(comp_name)
}

# make gene list for each dds object
gene_lists <- list()
universe_all <- c()
for (i in seq_along(object_list)) {
  comp_name <- gsub(".rds", "", gsub("DESeq_data_", "", basename(object_list[i])))
  dds_n <- readRDS(object_list[i])
  universe <- rownames(dds_n)
  universe_all <- c(universe_all, universe)
  res_dds_n <- results(dds_n, lfcThreshold = 0, alpha = 0.05, cooksCutoff = FALSE, independentFiltering = FALSE)
  res_df <- data.frame(
    gene = gene_map$GENENAME[match(rownames(res_dds_n), gene_map$gene_id)],
    log2FC = res_dds_n$log2FoldChange,
    padj = res_dds_n$padj
  )
  res_df <- na.omit(res_df)
  
  sig_genes <- res_df$gene[res_df$padj < 0.05 & abs(res_df$log2FC) > 1]
  
  # Assign the vector of significant genes directly to the list element named by comp_name
  gene_lists[[comp_name]] <- sig_genes
}

# estimate universal gene list
universe_all <- (unique(universe_all))

# Run compareCluster using enrichGO for each gene list.
# Adjust keyType and OrgDb as needed for your data.
cc <- compareCluster(
    geneCluster = gene_lists,
    fun = "enrichGO",
    OrgDb = org.Pverrucosa.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    universe = universe_all
)

# Create a dot plot for the combined results.
dotplot(cc, showCategory = 20, size = "FoldEnrichment") + 
  ggtitle("GO Enrichment Comparison") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Save the results to a file
write.csv(as.data.frame(cc), file = paste0(out_dir, "/summary/GO_enrichment_results.csv"), row.names = FALSE)
# Save the compareCluster object

## other comparison to explore separately
pdf(paste0(out_dir,"/plot/GO_enrichment_comparisons.pdf"), width = 8, height = 12)
# [1] "T0-ALL_TREATMENTS"
# [1] "T1-ALL_TREATMENTS"
# [1] "T2-ALL_TREATMENTS"
# Extract gene lists containing "ALL_TREATMENTS"
all_treatments_gene_lists <- gene_lists[grep("ALL_TREATMENTS", names(gene_lists))]
names(all_treatments_gene_lists)
all_treatments <- compareCluster(
    geneCluster = all_treatments_gene_lists,
    fun = "enrichGO",
    OrgDb = org.Pverrucosa.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    universe = universe_all
)

dotplot(all_treatments, showCategory = 20, size = "FoldEnrichment") + 
  ggtitle("GO Enrichment Comparison") + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# [1] "CC_T0-1-2"
# [1] "CP_T0-1-2"
# [1] "VD_T0-1-2"
# [1] "VL_T0-1-2"
all_times_gene_lists <- gene_lists[grep("^[CPV][CPLD]_T0-1-2$", names(gene_lists))]
names(all_times_gene_lists)

all_times<- compareCluster(
    geneCluster = all_times_gene_lists,
    fun = "enrichGO",
    OrgDb = org.Pverrucosa.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    universe = universe_all
)

dotplot(all_times, showCategory = 20, size = "FoldEnrichment") + 
  ggtitle("GO Enrichment Comparison - Single treatment through timepoint") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# [1] "CCvCP-T0"
# [1] "CCvVD-T0"
# [1] "CCvVL-T0"
# [1] "CPvVD-T0"
# [1] "CPvVL-T0"
# [1] "VDvVL-T0"
treatment_t0_gene_lists <- gene_lists[grep("^[CV][CPD]v[CV][PDL]-T0$", names(gene_lists))]
names(treatment_t0_gene_lists)

t0_t<- compareCluster(
    geneCluster = treatment_t0_gene_lists,
    fun = "enrichGO",
    OrgDb = org.Pverrucosa.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    universe = universe_all
)

dotplot(t0_t, showCategory = 20, size = "FoldEnrichment") + 
  ggtitle("GO Enrichment Comparison - Treatments comparison at T0") +
      # hrbrthemes::theme_ipsum() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# [1] "CCvCP-T1"
# [1] "CCvVD-T1"
# [1] "CCvVL-T1"
# [1] "CPvVD-T1"
# [1] "CPvVL-T1"
# [1] "VDvVL-T1"
treatment_t1_gene_lists <- gene_lists[grep("^[CV][CPD]v[CV][PDL]-T1$", names(gene_lists))]
names(treatment_t1_gene_lists)

t1_t<- compareCluster(
    geneCluster = treatment_t1_gene_lists,
    fun = "enrichGO",
    OrgDb = org.Pverrucosa.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    universe = universe_all
)

dotplot(t1_t, showCategory = 20, size = "FoldEnrichment") + 
  ggtitle("GO Enrichment Comparison - Treatments comparison at T1") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# [1] "CCvCP-T2"
# [1] "CCvVD-T2"
# [1] "CPvVL-T2"
# [1] "CPvVD-T2"
# [1] "CCvVL-T2"
# [1] "VDvVL-T2"
treatment_t2_gene_lists <- gene_lists[grep("^[CV][CPD]v[CV][PDL]-T2$", names(gene_lists))]
names(treatment_t2_gene_lists)

t2_t<- compareCluster(
    geneCluster = treatment_t2_gene_lists,
    fun = "enrichGO",
    OrgDb = org.Pverrucosa.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    universe = universe_all
)

dotplot(t2_t, showCategory = 20, size = "FoldEnrichment") + 
  ggtitle("GO Enrichment Comparison - Treatments comparison at T2") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# [1] "CCvCP_T0-1-2"
# [1] "CCvVD_T0-1-2"
# [1] "CCvVL_T0-1-2"
# [1] "CPvVD_T0-1-2"
# [1] "CPvVL_T0-1-2"
# [1] "VDvVL_T0-1-2"
treatment_all_time_gene_lists <- gene_lists[grep("^[CV][CPD]v[CV][PDL]_T0-1-2$", names(gene_lists))]
names(treatment_all_time_gene_lists)

tt <- compareCluster(
    geneCluster = treatment_all_time_gene_lists,
    fun = "enrichGO",
    OrgDb = org.Pverrucosa.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    universe = universe_all
)

dotplot(tt, showCategory = 20, size = "FoldEnrichment") + 
  ggtitle("GO Enrichment Comparison - Treatments comparison at T0-1-2") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
dev.off()

###--- fisher test ---###
gene2GO <- read_tsv("/home/colinl/Proj/microbiome_probiotics_RNASeq/Pver_ref/NCBIFilesDir/gene2go.gz")
# head(gene2GO)
# Filter gene2GO for taxid = 203993
gene2GO_filtered <- gene2GO %>% filter(`#tax_id` == 203993)

# Create a mapping of gene IDs to GO terms
geneID2GO <- split(gene2GO_filtered$GO_ID, gene2GO_filtered$GeneID)

# Run Fisher test on all DESeq2 objects
fisher_results <- list()
for (i in seq_along(object_list)) {
  comp_name <- gsub(".rds", "", gsub("DESeq_data_", "", basename(object_list[i])))
  dds_n <- readRDS(object_list[i])
  universe <- rownames(dds_n)
  
  res_dds_n <- results(dds_n, lfcThreshold = 0, alpha = 0.05, cooksCutoff = FALSE, independentFiltering = FALSE)
  res_df <- data.frame(
    gene = gene_map$GENENAME[match(rownames(res_dds_n), gene_map$gene_id)],
    log2FC = res_dds_n$log2FoldChange,
    padj = res_dds_n$padj
  ) |> na.omit()
  
  sig_genes <- res_df$gene[res_df$padj < 0.05 & abs(res_df$log2FC) > 1]
  
  geneList <- factor(as.integer(universe %in% sig_genes))
  # names(geneList) <- universe
  names(geneList) <- sub("LOC", "", universe)

  myGOdata <- new("topGOdata", 
                  description = paste("Fisher Test for", comp_name), 
                  ontology = "BP", 
                  allGenes = geneList, 
                  annot = annFUN.gene2GO, 
                  gene2GO = geneID2GO)
  
  resultFisher <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
  
  fisher_results[[comp_name]] <- GenTable(myGOdata, 
                                          weightFisher = resultFisher, 
                                          orderBy = "weightFisher", 
                                          ranksOf = "weightFisher", 
                                          topNodes = 20)
  
  fisher_results[[comp_name]]$adjusted.p <- p.adjust(fisher_results[[comp_name]]$weightFisher, method = "bonferroni", n = length(fisher_results[[comp_name]]$weightFisher))
  fisher_results[[comp_name]]$q.value <- p.adjust(fisher_results[[comp_name]]$weightFisher, method = "fdr", n = length(fisher_results[[comp_name]]$weightFisher))

  message(paste("Fisher test completed for", comp_name))
}

# # Save Fisher test results
for (comp_name in names(fisher_results)) {
  write.table(fisher_results[[comp_name]],
              file = paste0(out_dir, "/summary/fisher/fisher_results_", comp_name,".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}



# counts <- read.table("/home/colinl/Proj/microbiome_probiotics_RNASeq/results/star_salmon/salmon.merged.gene_counts.tsv", header = TRUE, sep = "\t", row.names = 1)
# universe <- data.frame(Gene = rownames(counts))
# # Read in annotation data from NCBI 
# all_annotations <- read_tsv("/home/colinl/Proj/microbiome_probiotics_RNASeq/Pver_ref/NCBIFilesDir/203993.gene.tsv")
# gene2GO <- read_tsv("/home/colinl/Proj/microbiome_probiotics_RNASeq/Pver_ref/NCBIFilesDir/gene2go.gz")

# head(gene2GO)
# # Filter gene2GO for taxid = 203993
# gene2GO_filtered <- gene2GO %>% filter(`#tax_id` == 203993)

# # Create a mapping of gene IDs to GO terms
# geneID2GO <- split(gene2GO_filtered$GO_ID, gene2GO_filtered$GeneID)
# # universe=data.frame(Gene=rownames(cnt))
# # universe$Uniprot=sprot_annotations$V2[match(universe$Gene, sprot_annotations$V1)]
# # universe$Uniprot=ifelse(is.na(universe$Uniprot), trembl_annotations$V2[match(universe$Gene, trembl_annotations$V1)], as.character(universe$Uniprot))
# # universe$GO=unip_meta$Gene.ontology.IDs[match(universe$Uniprot, unip_meta$Entry)]
# # is.na(universe$GO) <- universe$GO == ""
# # universe_final=na.omit(universe[,-2])
# # write.table(universe_final, "./Input_files/GO_universe", quote = F, row.names = F, sep = "\t")

# geneList <- factor(as.integer(universe_all %in% gene_lists[[1]]))
# geneList <- sub("LOC", "", names(geneList))
# names(geneList) <- sub("LOC", "", universe_all)

# myGOdata <- new("topGOdata", 
#         description = "My project", 
#         ontology = "BP", 
#         allGenes = geneList, 
#         annot = annFUN.gene2GO, 
#         gene2GO = geneID2GO)

# resultFisher <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")

# p_load(gtsummary)

# ##GO enerichment
# p_load(topGO)
# #library(data.table)

# #1. get the universe of genes
# universe=data.frame(Gene=rownames(cnt))
# universe$Uniprot=sprot_annotations$V2[match(universe$Gene, sprot_annotations$V1)]
# universe$Uniprot=ifelse(is.na(universe$Uniprot), trembl_annotations$V2[match(universe$Gene, trembl_annotations$V1)], as.character(universe$Uniprot))
# universe$GO=unip_meta$Gene.ontology.IDs[match(universe$Uniprot, unip_meta$Entry)]
# is.na(universe$GO) <- universe$GO == ""
# universe_final=na.omit(universe[,-2])
# write.table(universe_final, "./Input_files/GO_universe", quote = F, row.names = F, sep = "\t")

# geneID2GO = readMappings(file = "Input_files/GO_universe")
# geneUniverse <- names(geneID2GO)
# #length(geneUniverse) #15388

# #2. Define your enriched genes
# enr1=read.table("outputs/DESeq_results.txt", header = T, row.names = 1, sep = "\t") 

# #3. do the GO enrichment
# geneList = factor(as.integer(geneUniverse %in% rownames(enr1)))
# names(geneList) <- geneUniverse
# myGOdata= new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO=geneID2GO)
# resultFisher = runTest(myGOdata, algorithm="weight01", statistic="fisher") # in the classic algorithm hierarchy isn't taken into account, so each GO term is tested independently.To take the GO hierarchy into account, we use algorithm='weight01'. 

# #adjusting p-values, using Bonferroni and FDRs (q-value)
# allGO = usedGO(object = myGOdata)
# allRes = GenTable(myGOdata, weightFisher = resultFisher, orderBy = "resultFisher", ranksOf = "weightFisher", topNodes = length(allGO))
# options(scipen = 999)
# allRes$adjusted.p = p.adjust(allRes$weightFisher, method = "bonferroni", n = length(allRes$weightFisher))
# allRes$q.value= p.adjust(allRes$weightFisher, method = "fdr", n = length(allRes$weightFisher))
# write.table(allRes, "TopGO_BP.txt", row.names = FALSE, sep = "\t", quote = FALSE)

browseURL("http://localhost:8787")
allRes %>%
  tbl_summary() %>% 
  gt::gtsave("output.png")
  # add_overall() %>%
  # add_p() #%>%
  #add_stat_label()
