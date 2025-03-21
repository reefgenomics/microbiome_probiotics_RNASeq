
![Header Image](https://github.com/nf-core/rnaseq/raw/3.18.0/docs/images/nf-core-rnaseq_logo_light.png)

# nf-core/rnaseq pipeline results

this folder will contain the results of the nf-core/rnaseq pipeline and the output form the DEG analysis see: [DESeq2 README.md](Rscripts/README.md)
****
Expected output:

```text
├── deseq2_R
│   ├── plots
│   ├── summary
│   └── etc ... # Not complete currently.
├── fastp
│   ├── 2024_Pver_CC_T0_01_01.fastp.html
│   ├── 2024_Pver_CC_T0_01_01.fastp.json
│   ├── ...
│   ├── 2024_Pver_CC_T0_02_04.fastp.json
│   └── log
├── fastqc
│   ├── raw
│   └── trim
├── fq_lint
│   ├── raw
│   └── trimmed
├── multiqc
│   └── star_salmon
├── pipeline_info
│   ├── 'multiple logs'
├── salmon
│   ├── 2024_Pver_CC_T0_01_01
│   ├── 2024_Pver_CC_T0_02_04
│   ├── ...
│   ├── deseq2_qc
│   ├── salmon.merged.gene_counts_length_scaled.SummarizedExperiment.rds
│   ├── salmon.merged.gene_counts_length_scaled.tsv
│   ├── salmon.merged.gene_counts_scaled.SummarizedExperiment.rds
│   ├── salmon.merged.gene_counts_scaled.tsv
│   ├── salmon.merged.gene_counts.SummarizedExperiment.rds
│   ├── salmon.merged.gene_counts.tsv
│   ├── salmon.merged.gene_lengths.tsv
│   ├── salmon.merged.gene_tpm.tsv
│   ├── salmon.merged.transcript_counts.SummarizedExperiment.rds
│   ├── salmon.merged.transcript_counts.tsv
│   ├── salmon.merged.transcript_lengths.tsv
│   ├── salmon.merged.transcript_tpm.tsv
│   └── tx2gene.tsv
└── star_salmon
    ├── 2024_Pver_CC_T0_01_01
    ├── 2024_Pver_CC_T0_01_01.markdup.sorted.bam
    ├── 2024_Pver_CC_T0_01_01.markdup.sorted.bam.bai
    ├── ...
    ├── 2024_Pver_VL_T2_04_02.markdup.sorted.bam.bai
    ├── bigwig
    ├── deseq2_qc
    ├── dupradar
    ├── log
    ├── picard_metrics
    ├── qualimap
    ├── rseqc
    ├── salmon.merged.gene_counts_length_scaled.SummarizedExperiment.rds
    ├── salmon.merged.gene_counts_length_scaled.tsv
    ├── salmon.merged.gene_counts_scaled.SummarizedExperiment.rds
    ├── salmon.merged.gene_counts_scaled.tsv
    ├── salmon.merged.gene_counts.SummarizedExperiment.rds
    ├── salmon.merged.gene_counts.tsv
    ├── salmon.merged.gene_lengths.tsv
    ├── salmon.merged.gene_tpm.tsv
    ├── salmon.merged.transcript_counts.SummarizedExperiment.rds
    ├── salmon.merged.transcript_counts.tsv
    ├── salmon.merged.transcript_lengths.tsv
    ├── salmon.merged.transcript_tpm.tsv
    ├── samtools_stats
    ├── stringtie
    └── tx2gene.tsv
```

## References

- [nf-core/rnaseq documentation](https://nf-co.re/rnaseq)
