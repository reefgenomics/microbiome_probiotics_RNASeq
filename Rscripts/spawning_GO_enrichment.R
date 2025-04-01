##GO enerichment
library(topGO)
#library(data.table)

#1. get the universe of genes
# universe=data.frame(Gene=rownames(cnt))
# universe$Uniprot=sprot_annotations$V2[match(universe$Gene, sprot_annotations$V1)]
# universe$Uniprot=ifelse(is.na(universe$Uniprot), trembl_annotations$V2[match(universe$Gene, trembl_annotations$V1)], as.character(universe$Uniprot))
# universe$GO=unip_meta$Gene.ontology.IDs[match(universe$Uniprot, unip_meta$Entry)]
# is.na(universe$GO) <- universe$GO == ""
# universe_final=na.omit(universe[,-2])
# write.table(universe_final, "./Input_files/GO_universe", quote = F, row.names = F, sep = "\t")

geneID2GO = readMappings(file = "Input_files/GO_universe")
geneUniverse <- names(geneID2GO)
#length(geneUniverse) #15388

#2. Define your enriched genes
enr1=read.table("outputs/DESeq_results.txt", header = T, row.names = 1, sep = "\t") 

#3. do the GO enrichment
geneList = factor(as.integer(geneUniverse %in% rownames(enr1)))
names(geneList) <- geneUniverse
myGOdata= new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO=geneID2GO)
resultFisher = runTest(myGOdata, algorithm="weight01", statistic="fisher") # in the classic algorithm hierarchy isn't taken into account, so each GO term is tested independently.To take the GO hierarchy into account, we use algorithm='weight01'. 

#adjusting p-values, using Bonferroni and FDRs (q-value)
allGO = usedGO(object = myGOdata)
allRes = GenTable(myGOdata, weightFisher = resultFisher, orderBy = "resultFisher", ranksOf = "weightFisher", topNodes = length(allGO))
options(scipen = 999)
allRes$adjusted.p= p.adjust(allRes$weightFisher, method = "bonferroni", n = length(allRes$weightFisher))
allRes$q.value= p.adjust(allRes$weightFisher, method = "fdr", n = length(allRes$weightFisher))
write.table(allRes, "TopGO_BP.txt", row.names = FALSE, sep = "\t", quote = FALSE)



###----- my functional erichment 

library(GenomicFeatures)
library(GenomicFeatures)
library(AnnotationDbi)
p_load(clusterProfiler)
# p_load(org.Hs.eg.db) # Substitute with coral-specific DB if available
# install.packages("/home/colinl/Proj/microbiome_probiotics_RNASeq/Pver_ref/org.Pverrucosa.eg.db", repos = NULL) 
library(OrgDb) # Substitute with coral-specific DB if available
library(AnnotationForge)

# makeOrgPackageFromNCBI(version = "0.1",
#                        author = "Luigi Colin <luigi.colin@uni-konstanz.de>",
#                        maintainer = "Luigi Colin <luigi.colin@uni-konstanz.de>",
#                        outputDir = "/home/colinl/Proj/microbiome_probiotics_RNASeq/Pver_ref/.",
#                        NCBIFilesDir = "/home/colinl/Proj/microbiome_probiotics_RNASeq/Pver_ref/NCBIFilesDir/.",
#                        tax_id = "203993",
#                        genus = "Pocillopora",
#                        species = "verrucosa")
library(org.Pverrucosa.eg.db)
gff_file <- "/home/colinl/Proj/microbiome_probiotics_RNASeq/Pver_ref/genomic.simple.gtf.gz"

txdb <- GenomicFeatures::makeTxDbFromGFF(gff_file, format = "gtf")
gene_map <- GenomicFeatures::genes(txdb)
valid_columns <- AnnotationDbi::columns(txdb)
print(valid_columns) # Check for valid column names
gene_map$GENENAME <- AnnotationDbi::mapIds(txdb, keys = gene_map$gene_id, 
                                           column = "GENEID", keytype = "GENEID", 
                                           multiVals = "first")
gene_map$GENENAME <- sub("gene-", "", gene_map$GENENAME)
gene_map$gene_id <- sub("gene-", "", gene_map$gene_id)


# select(org.Pverrucosa.eg.db, keys=ensids, columns=cols, keytype="ENSEMBL")

# # Extract mapping of NCBI gene symbols to UniProtKB/TrEMBL gene IDs
# ncbi_to_uniprot <- AnnotationDbi::select(
#   org.Pverrucosa.eg.db,
#   keys = keys(org.Pverrucosa.eg.db, keytype = "SYMBOL"),
#   columns = c("SYMBOL", "UNIPROT"),
#   keytype = "SYMBOL"
# )

# # Remove duplicates and NA values
# ncbi_to_uniprot <- na.omit(unique(ncbi_to_uniprot))
# head(ncbi_to_uniprot)
# # ##extract Deseq_resutls for CP vs CC
# # res_CP <- lfcShrink(C0, coef = "treatment_CP_vs_CC", type = "apeglm")
# # res_df_CP <- data.frame(
# #   gene = gene_map$GENENAME[match(rownames(res_CP), gene_map$gene_id)],
# #   log2FC = res_CP$log2FoldChange,
# #   padj = res_CP$padj
# # ) |> na.omit()
# # ids <- keys(org.Pverrucosa.eg.db) 
# # columns(org.Pverrucosa.eg.db)
# # cols <- c("SYMBOL", "GENENAME")
# select(org.Pverrucosa.eg.db, keys=ensids, columns=cols, keytype="ENSEMBL")
# # mapping <- select(org.Pverrucosa.eg.db, keys=ids, columns=c('ENTREZID','SYMBOL'), keytype='ENTREZID')
# # mapping[duplicated(mapping[,"SYMBOL"]),] 

# # select(org.Pverrucosa.eg.db, keys="HBD", columns=c('ENTREZID','SYMBOL'), keytype='SYMBOL')

# # select(org.Pverrucosa.eg.db, keys="MEMO1", columns=c('ENTREZID','SYMBOL'), keytype='SYMBOL')
# # select(org.Pverrucosa.eg.db, keys="TRNAV-CAC", columns=c('ENTREZID','SYMBOL'), keytype='SYMBOL')

object_list <- list.files(path = "/home/colinl/Proj/microbiome_probiotics_RNASeq/results/deseq2_R/dds_objects/", pattern = "*.rds", full.names = TRUE)

# # read RDS with DEseq2 results # this is worng, better to import the tsv tables. 
# dds_CP <- readRDS("/home/colinl/Proj/microbiome_probiotics_RNASeq/results/deseq2_R/dds_objects/DESeq_data_T2-ALL_TREATMENTS.rds")
# # All genes surviving initial filtering and included in DESeq2 analysis

# # read tsv as replacement for the above
# res_df_CP <- read.table("/home/colinl/Proj/microbiome_probiotics_RNASeq/results/deseq2_R/DESeq_Comparison-T2-ALL_TREATMENTS.tsv", header = T, sep = "\t")
# # gene        log2FC      padj
# head(res_df_CP)
# # Gene universe (background) for CP vs CC
# # universe_CP <- res_df_CP$Symbol
# universe_CP <- rownames(dds_CP)#[rowSums(counts(dds_CP, normalized=TRUE) >= 10) >= 6]

# # Significant genes (padj < 0.05, |log2FC| > 1) for CP vs CC
# sig_genes_CP <- res_df_CP$Symbol[res_df_CP$padj < 0.05 & abs(res_df_CP$log2FoldChange) > 1]

# # GO Enrichment for CP vs CC
# ego_CP <- enrichGO(
#   gene = sig_genes_CP,
#   universe = universe_CP,
#   OrgDb = org.Pverrucosa.eg.db, # Use coral annotation if possible
#   keyType = "SYMBOL",
#   ont = "BP", # Biological Process
#   pAdjustMethod = "BH",
#   qvalueCutoff = 0.05
# )

# ##extract Deseq_resutls for VL vs VD
# res_VL <- lfcShrink(dds, coef = "treatment_VL_vs_CC", type = "apeglm")
# res_df_VL <- data.frame(
#   gene = gene_map$GENENAME[match(rownames(res_VL), gene_map$gene_id)],
#   log2FC = res_VL$log2FoldChange,
#   padj = res_VL$padj
# ) |> na.omit()

# read tsv as replacement for the above
res_df_VL <- read.table("/home/colinl/Proj/microbiome_probiotics_RNASeq/results/deseq2_R/DESeq_Comparison-Comparison CC vs CP for T0.tsv", header = T, sep = "\t")


object_list <- list.files(path = "/home/colinl/Proj/microbiome_probiotics_RNASeq/results/deseq2_R/dds_objects/", pattern = "*.rds", full.names = TRUE)


httpgd::hgd()
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
      next
    }
    # experiments # make object with ego result if not empty.
    assign(paste0("ego_", comp_name), ego_dds_n)
    message(paste("ego_", comp_name, " made!", sep = ""))
    # Plot the results
      print(dotplot(ego_dds_n, showCategory = 20) + ggtitle(paste("GO Enrichment (", comp_name, ")", sep = "")))
}

names(as.data.frame(ego_dds_n))
cnetplot(ego_dds_n)
barplot(ego_dds_n)
dotplot(ego_dds_n, x = "FoldEnrichment")
dotplot(ego_dds_n, x = "RichFactor")
dotplot(ego_dds_n, x = "zScore", size = "qvalue")
dotplot(ego_dds_n, x = "FoldEnrichment", color = "qvalue")
dotplot(ego_dds_n, x = "FoldEnrichment")

# # Generate a heatmap from the ego object
# ego_matrix <- as.data.frame(ego_dds_n)[, c("ID", "Description", "p.adjust")]
# ego_matrix <- ego_matrix[order(ego_matrix$p.adjust), ]
# rownames(ego_matrix) <- ego_matrix$Description
# heatmap_data <- -log10(ego_matrix$p.adjust)
# heatmap_data <- as.matrix(heatmap_data)

# # Create the heatmap
# p_load("pheatmap")
# pheatmap(
#   heatmap_data,
#   cluster_rows = TRUE,
#   cluster_cols = FALSE,
#   display_numbers = TRUE,
#   main = "GO Enrichment Heatmap",
#   color = colorRampPalette(c("white", "blue"))(50)
# )
# heatmap(ego_dds_n)

for (i in seq_along(object_list)) {
  comp_name <- gsub(".rds", "", gsub("DESeq_data_", "", basename(object_list[i])))
  print(comp_name)
}

# Load necessary packages
library(clusterProfiler)
library(enrichplot)

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

universe_all <- (unique(universe_all))


# # Suppose you have several gene lists for different comparisons:
# gene_list1 <- c("geneA", "geneB", "geneC")
# gene_list2 <- c("geneD", "geneE", "geneF")
# gene_list3 <- c("geneG", "geneH", "geneI")

# Combine your gene lists into a list. Name each list for clarity.
gene_lists[[1]]
# Run compareCluster using enrichGO for each gene list.
# Adjust keyType and OrgDb as needed for your data.
cc <- compareCluster(
  geneCluster = gene_lists,
  fun = "enrichGO",
  OrgDb = org.Pverrucosa.eg.db,  # Your custom organism package
  keyType = "SYMBOL",           # Use the proper key type for your gene IDs
  ont = "BP",                    # You can choose "BP", "MF", or "CC"
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  universe = universe_all
)

# Create a dot plot for the combined results.
dotplot(cc, showCategory = 20, size = "FoldEnrichment") + 
  ggtitle("GO Enrichment Comparison") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# dds_VL <- readRDS("/home/colinl/Proj/microbiome_probiotics_RNASeq/results/Out_deseq/DESeq_data_CCvsVL_T2.rds")
# # Gene universe (background) for VL vs VD
# # universe_VL <- res_df_VL$Symbol
# universe_VL <- rownames(dds_VL)#[rowSums(counts(dds_CP, normalized=TRUE) >= 10) >= 6]

# # Significant genes (padj < 0.05, |log2FC| > 1) for VL vs VD
# sig_genes_VL <- res_df_VL$Symbol[res_df_VL$padj < 0.05 & abs(res_df_VL$log2FoldChange) > 1]

# # GO Enrichment for VL vs VD
# ego_VL <- enrichGO(
#   gene = sig_genes_VL,
#   universe = universe_VL,
#   OrgDb = org.Pverrucosa.eg.db, # Use coral annotation if possible
#   keyType = "SYMBOL",
#   ont = "BP", # Biological Process
#   pAdjustMethod = "BH",
#   qvalueCutoff = 0.05
# )

# ##extract Deseq_resutls for VL vs VD
# # res_VD <- lfcShrink(dds, coef = "treatment_VD_vs_CC", type = "apeglm")

# # res_df_VD <- data.frame(
# #   gene = gene_map$GENENAME[match(rownames(dds), gene_map$gene_id)],
# #   log2FC = res_VD$log2FoldChange,
# #   padj = res_VD$padj
# # ) |> na.omit()

# res_df_VD <- read.table("/home/colinl/Proj/microbiome_probiotics_RNASeq/results/Out_deseq/DESeq_Comparison-Comparison CC vs VD for T2.tsv", header = T, sep = "\t")

# # Gene universe (background) for VL vs VD
# # universe_VD <- res_df_VD$Symbol
# dds_VD <- readRDS("/home/colinl/Proj/microbiome_probiotics_RNASeq/results/Out_deseq/DESeq_data_CCvsVD_T2.rds")
# # Gene universe (background) for VL vs VD
# universe_VD <- rownames(dds_VD)#[rowSums(counts(dds_CP, normalized=TRUE) >= 10) >= 6]

# # Significant genes (padj < 0.05, |log2FC| > 1) for VL vs VD
# sig_genes_VD <- res_df_VD$Symbol[res_df_VD$padj < 0.05 & abs(res_df_VD$log2FoldChange) > 1]

# # GO Enrichment for VL vs VD
# ego_VD <- enrichGO(
#   gene = sig_genes_VD,
#   universe = universe_VD,
#   OrgDb = org.Pverrucosa.eg.db, # Use coral annotation if possible
#   keyType = "SYMBOL",
#   ont = "BP", # Biological Process
#   pAdjustMethod = "BH",
#   qvalueCutoff = 0.05
# )


# httpgd::hgd()
# # Combine the enrichment results
# # combined_ego <- list(CP_vs_CC = ego_CP, VL_vs_VD = ego_VL)

# # Plot the combined results
# dotplot(ego_CP, showCategory = 20) + ggtitle("GO Enrichment (CP vs CC)")
# goplot(ego_CP)
# dotplot(ego_VL, showCategory = 20) + ggtitle("GO Enrichment (VL vs CC)")
# dotplot(ego_VD, showCategory = 20) + ggtitle("GO Enrichment (VD vs CC)")


# res_df_VD <- read.table("/home/colinl/Proj/microbiome_probiotics_RNASeq/results/Out_deseq/DESeq_Comparison-ALL.tsv", header = T, sep = "\t")

# # Gene universe (background) for VL vs VD
# # universe_VD <- res_df_VD$Symbol
# # dds_VD <- readRDS("/home/colinl/Proj/microbiome_probiotics_RNASeq/results/Out_deseq/DESeq_data_CCvsVD_T2.rds")
# # Gene universe (background) for VL vs VD
# universe_VD <- rownames(C0)#[rowSums(counts(dds_CP, normalized=TRUE) >= 10) >= 6]

# # Significant genes (padj < 0.05, |log2FC| > 1) for VL vs VD
# sig_genes_VD <- res_df_VD$Symbol[res_df_VD$padj < 0.05 & abs(res_df_VD$log2FoldChange) > 1]

# # GO Enrichment for VL vs VD
# ego_VD <- enrichGO(
#   gene = sig_genes_VD,
#   universe = universe_VD,
#   OrgDb = org.Pverrucosa.eg.db, # Use coral annotation if possible
#   keyType = "SYMBOL",
#   ont = "BP", # Biological Process
#   pAdjustMethod = "BH",
#   qvalueCutoff = 0.05
# )

# ### -- exploraratory 

# # library(fenr) nice but not an option for coral data >> not enough refernce data
# # go_data <- fetch_go(species = "203993") 
# # go_terms <- prepare_for_enrichment(go_data$terms, go_data$mapping, universe)

# # enr <- functional_enrichment(
# #   background = universe,
# #   selected = sig_genes,
# #   terms = go_terms
# # ) 

# ## network aware 
# library(pathfindR)
# input_df <- res_df_CP[, c("Symbol", "log2FoldChange", "padj")]
# colnames(input_df) <- c("Gene.symbol", "logFC", "adj.P.Val")
# tail(input_df)

# gsets_list <- get_gene_sets_list(
#   source = "KEGG",
#   org_code = "mmu"
# )

# ### to do: create custom pocillopora pdam(?) kegg references

# output <- run_pathfindR(
#   input_df,
#   gene_sets = "GO-All", #'KEGG', 'Reactome', 'BioCarta','GO-All', 'GO-BP', 'GO-CC', 'GO-MF', 'cell_markers','mmu_KEGG' or 'Custom'. If 'Custom', the arguments ‘custom_genes’ and ‘custom_descriptions’
#   pin_name_path = "mmu_STRING", #c('Biogrid', 'STRING', 'GeneMania','IntAct', 'KEGG', 'mmu_STRING')
#   iterations = 10 # For stability
# )

# heatmap_plot <- enrichment_chart(output) # Customizable visualization
