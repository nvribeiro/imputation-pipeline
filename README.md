# Pipeline for preparing CeDNN genotypes for imputation
This pipeline was developed to convert a Final Report from Illumina Genome Studio (using chip GSA v3) to VCF files that can be used for imputation using the Michigan Server. The reference genome is GRCh37 and the referebce haplotype panel is the Haplotype Reference Consortium (HRC).

Note: this pipeline was developed to have imputated genotypes to demultiplex single-cell sequencing - it was not tested and validade for other uses (e.g. GWAS).

Created by: Nathan V Ribeiro <n.v.ribeiro@umcg.nl> - UMCG Immunogenetics Group

# How to use
Get the codes

Update the file config/config.yaml

Run