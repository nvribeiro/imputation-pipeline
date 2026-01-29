# Pipeline for preparing CeDNN genotypes for imputation
This pipeline was developed to convert a Final Report from Illumina Genome Studio (using chip GSA v3) to VCF files that can be used for imputation using the Michigan Server. The reference genome is GRCh37 and the referebce haplotype panel is the Haplotype Reference Consortium (HRC).

Note: this pipeline was developed to have imputated genotypes to demultiplex single-cell sequencing - it was not tested and validade for other uses (e.g. GWAS).

Created by: Nathan V Ribeiro <n.v.ribeiro@umcg.nl> - UMCG Immunogenetics Group

## How to use
### 1. Get the codes

### 2. Create the sampleinfo file
To create this file, you need to know the samples that are included in your FinalReport.txt. A quick way of doing that is by typing the following line of code in the terminal:
`cut -f2 FinalReport.txt | tail -n +2 | sort | uniq`

This will return the list of samples so you can look for the sex and phenotype in the CeDNN database. The `sample_info.csv` file should be a comma-separated text file with the following fields: IID (the CeDNN sample ID as specified in the FinalReport.txt), Sex (1 = male, 2 = female) and Phenotype (1 = control, 2 = celiac). If this information is not relevant for you, you can just fill Sex and Phenotype with 0s. Below is an example of how this file should look like.

IID,Sex,Phenotype
CeDNN_1234,1,1
CeDNN_1244,1,2
CeDNN_1254,2,2

Save this file as `sample_info.csv` in the folder where your pipeline is.

### 3. Update the file nextflow.config

### 4. Run