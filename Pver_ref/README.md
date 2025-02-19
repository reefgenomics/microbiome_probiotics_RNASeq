# Reference Genome Information

The reference genome used is Pocillopora verrucosa genome v1.0. For more detailed information, please refer to the open access Genome Biology and Evolution publication, including the [main text](https://academic.oup.com/gbe/article/12/10/1911/5898631#209703315%22) and [supplemental information](https://academic.oup.com/gbe/article/12/10/1911/5898631).

File availability at:
[http://pver.reefgenomics.org/](http://pver.reefgenomics.org/)

![Header Image](http://pver.reefgenomics.org/img/title.png)

## Convert GFF3 to GTF

The nf-core/rnaseq 3.18.0 pipeline does not specifically require a GTF file, but I encountered issues with the automatic conversion from GFF to GTF.

Here are the steps I took for manual conversion:

```bash
# Install AGAT with mamba in the work or base environment 
mamba install agat -y

# Convert with AGAT 
agat_convert_sp_gff2gtf.pl --gff Pver_ref/Pver_genome_assembly_v1.0.gff.gz -o Pver_genome_assembly_v1.0.gtf

# Compress the file
pigz Pver_genome_assembly_v1.0.gtf
```
