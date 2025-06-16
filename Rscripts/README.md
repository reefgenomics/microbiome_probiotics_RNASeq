# Differential Expression Analysis with DESeq2

## overview

- Creating R environment for server. Optional: depends on work environment

```bash
  mamba create -n r_env_DESeq2 -c conda-forge -c bioconda r-base bioconductor-deseq2 r-ellipse r-ggplot2 r-dplyr r-tidyr \
  r-data.table bioconductor-phyloseq r-patchwork r-ggrepel r-pacman bioconductor-summarizedexperiment bioconductor-tximport bioconductor-apeglm bioconductor-GenomicFeatures \
  r-httpgd radian
```

- Import the count data into R.
- Perform differential expression analysis using DESeq2.
- Example R script:

  ```R
  library(DESeq2)
  countData <- read.csv("results/counts.csv", row.names = 1)
  colData <- read.csv("samplesheet.csv", row.names = 1)
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
  dds <- DESeq(dds)
  results <- results(dds)
  ```

## Interpretation and Visualization

- Interpret the results and visualize them using tools like ggplot2 or pheatmap.

# R script logic

- DESeq2.r: Main DESeq comparisons. Relies on nf-core/rnaseq, Output in [deseq2_R](results/deseq2_R)
- Ordination_plot.r: Main visualization. Relies on [DESeq2.r](Rscripts/DESeq2.r). Output in [deseq2_R/plot](results/deseq2_R/plots)
- Heatmap.r: Main visualization. Relies on [nf-core/rnaseq](https://github.com/nf-core/rnaseq/). Output in [deseq2_R/plot](results/deseq2_R/heatmap)
- spawning_GO_enrichment.R: enrichment analysis. Relies on [DESeq2.r](Rscripts/DESeq2.r).
- Permanova.r: PERMANOVA analysis. Relies on [DESeq2.r](Rscripts/DESeq2.r). Output in [deseq2_R/permanova_results](results/deseq2_R/permanova_results)

## References

- [DESeq2 documentation](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
