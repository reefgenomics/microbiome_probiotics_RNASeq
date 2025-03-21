# Reference Genome Information

The reference genome used is Pocillopora verrucosa [(ASM3666991v2)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036669915.1/). For more detailed information, please refer to the publication, [Genome-Wide Analysis of Cell Cycle-Regulating Genes in the Symbiotic Dinoflagellate Breviolum minutum](https://pubmed.ncbi.nlm.nih.gov/31551286/).

File availability at:
[NCBI genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036669915.1/)

## Get genome with NCBI dataset

```bash
# Install AGAT with mamba in the work or base environment 
mamba install -c conda-forge ncbi-datasets-cli

# Download genome and clean up
datasets download genome accession GCF_036669915.1 --include gff3,rna,cds,protein,genome,seq-report
unzip ncbi_dataset.zip
mv Pver_ref/download/ncbi_dataset/data/GCF_036669915.1/*.fna .
mv Pver_ref/download/ncbi_dataset/data/GCF_036669915.1/*.gff .
rm -r ncbi_dataset ncbi_dataset.zip
```

## Convert GFF3 to GTF

The nf-core/rnaseq 3.18.0 pipeline does not specifically require a GTF file, but I encountered issues with the automatic conversion from GFF to GTF.

Here are the steps I took for manual conversion:

```bash
# Install AGAT with mamba in the work or base environment 
mamba install agat -y

# Convert with AGAT  /home/colinl/Proj/microbiome_probiotics_RNASeq/Pver_ref/download/ncbi_dataset/data/GCF_036669915.1
agat_convert_sp_gff2gtf.pl --gff GCF_036669915.1/genomic.gff -o GCF_036669915.1/genomic.gtf ## --gtf_version relax  ?

# simplify gtf to include only gene id, transcript id and gene biotype. Fully formatted GTF from NCBI refseq not properly recognised by the pipeline.
python GTF_simplify.py genomic.gtf genomic.simple.gtf

# Compress the file
pigz *.gtf *.gff
```

## Simplify GTF

Fully formatted GTF from NCBI refseq not properly recognised by the pipeline. simplifying it to gtf to include only gene id, transcript id and gene biotype prevent incompatibility the issues with [nf-core/rnaseq](https://nf-co.re/rnaseq).
GTF_simplify.py is a custom python script in [Pver_ref](Pver_ref/GTF_simplify.py)

```bash
# simplify gtf to include only gene id, transcript id and gene biotype. Relies on "re" and "argparse"
python GTF_simplify.py genomic.gtf genomic.simple.gtf

# Compress the file
pigz *.gtf *.gff *.fna
```

<details>
<summary>Other refernce needed for Rscripts</summary>


### Make DB from NCBI for Gene ID for 203993

```R
library(AnnotationForge)

makeOrgPackageFromNCBI(version = "0.1",
                       author = "Luigi Colin <luigi.colin@uni-konstanz.de>",
                       maintainer = "Luigi Colin <luigi.colin@uni-konstanz.de>",
                       outputDir = "/home/colinl/Proj/microbiome_probiotics_RNASeq/Pver_ref/.",
                       NCBIFilesDir = "/home/colinl/Proj/microbiome_probiotics_RNASeq/Pver_ref/NCBIFilesDir/.", # only specified if pre-downloaded (ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/) or intend to keep them for further use
                       tax_id = "203993",
                       genus = "Pocillopora",
                       species = "verrucosa")

install.packages("/home/colinl/Proj/microbiome_probiotics_RNASeq/Pver_ref/org.Pverrucosa.eg.db", repos = NULL) 
```
Generated db files here: [org.Pverrucosa.eg.db](Pver_ref/org.Pverrucosa.eg.db.tar.gz)

### Bash Command to Extract Gene Info from NCBI

```bash
datasets summary gene taxon 203993 --as-json-lines | dataformat tsv gene --fields gene-id,gene-type,symbol,description,tax-id,tax-name > Pver_ref/NCBIFilesDir/203993.gene.tsv
```

- NCBI datasets and dataformat [documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/getting_started/)

</details>