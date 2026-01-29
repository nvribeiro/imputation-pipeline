########################################################################
# Genotype imputation from CeDNN samples genotyped using Illumina GSA v3
# This pipeline takes a Final Report from Illumina GenomeStudio and creates the vcf (split per chromosome, ignores XY) ready for imputation
# Note: this pipeline was developed to have imputated genotypes to demultiplex single-cell sequencing - it was not tested and validade for other uses (e.g. GWAS)
# Created by: Nathan V Ribeiro <n.v.ribeiro@umcg.nl> - UMCG Immunogenetics Group
# Last modified: 
########################################################################

########################################
# 1. Global configuration & constants
########################################
configfile: "config/config.yaml"

DATASET = config["dataset"]
FINAL_REPORT = config["input"]["final_report"]
SAMPLE_INFO = config["input"]["sample_info"]
OUTPUT_DIR = config["output_dir"]
CHR = config["chromosomes"]

HRC_SCRIPT = config["HRC_script"]
HRC_REF = config["HRC_reference"]

PYTHON2_DIR = config["python2_path"]
OPENSSL_DIR = config["openssl_path"]
CHECKVCF_REF = config["checkvcf_ref"]

########################################
# 2. Final targets (what "done" means)
########################################

rule all:
    input:
        expand(f"{OUTPUT_DIR}/VCF_compressed/{DATASET}_final-chr{{chr}}.vcf.gz", chr=CHR)

########################################
# 3. Rules (one per pipeline step)
########################################
rule finalreport_to_lgen:
    input: 
        final_report = FINAL_REPORT,
        sample_info = SAMPLE_INFO,
    output:
        lgen = f"{OUTPUT_DIR}/{DATASET}.lgen",
        map = f"{OUTPUT_DIR}/{DATASET}.map",
        fam = f"{OUTPUT_DIR}/{DATASET}.fam"
    resources:
        mem="32gb"
    shell:
        """
        module load RPlus/4.2.1-foss-2022a-v22.12.1
        Rscript scripts/create_lgen.R \
            {input.final_report} \
            {input.sample_info} \
            {output.lgen} \
            {output.map} \
            {output.fam}
        """

rule lgen_to_plink:
    input:
        f"{OUTPUT_DIR}/{DATASET}"
    output:
        f"{OUTPUT_DIR}/{DATASET}"
    shell:
        """
        module load PLINK/1.9-beta6-20190617
        plink --lfile {input} --autosome \
        --make-bed --out {output}
        """

rule update_build:
    input:
        dataset = f"{OUTPUT_DIR}/{DATASET}",
        strand = "/groups/umcg-immunogenetics/tmp02/users/NathanRibeiro/tools/GSA_strand_v3/GSAMD-24v3-0-EA_20034606_A1-b37.strand"
    output:
        f"{OUTPUT_DIR}/GSA_updated/{DATASET}_GSA_updated"
    resources:
        mem="16gb"
    shell:
        """
        module load PLINK/1.9-beta6-20190617
        scripts/update_build.sh \
            {input.dataset} \
            {input.strand} \
            {output}
        """

rule plink_freq:
    input:
        f"{OUTPUT_DIR}/GSA_updated/{DATASET}_GSA_updated"
    output:
        f"{OUTPUT_DIR}/GSA_updated/{DATASET}_GSA_updated"
    resources:
        mem="16gb"
    shell:
    """
    module load PLINK/1.9-beta6-20190617
    plink --freq --bfile {input} --out {output}
    """

rule make_hrc_vcf:
    input:
        bim_file = f"{OUTPUT_DIR}/GSA_updated/{DATASET}_GSA_updated.bim",
        frq_file = f"{OUTPUT_DIR}/GSA_updated/{DATASET}_GSA_updated.frq",
        ref_file = HRC_REF
    output:
        f"{OUTPUT_DIR}/GSA_updated/Run-plink.sh"
    resources:
        mem="32gb"
    shell:
        """
        module load foss/2022a
        module load PerlPlus/5.34.1-GCCcore-11.3.0-v22.11.1
        module load PLINK/1.9-beta6-20190617

        perl {HRC_SCRIPT} \
            -b {input.bim_file} \
            -f {input.frq_file} \
            -r {input.ref_file} \
            -h
        
        sh {output}

        mkdir {OUTPUT_DIR}/HRC_formatted
        mv {OUTPUT_DIR}/GSA_updated/{DATASET}_GSA_updated-updated* {OUTPUT_DIR}/HRC_formatted
        mv {OUTPUT_DIR}/GSA_updated/*.txt {OUTPUT_DIR}/HRC_formatted
        mv {output} {OUTPUT_DIR}/HRC_formatted
        """

rule check_vcf:
    input:
        ref = CHECKVCF_REF,
        vcf = expand(f"{OUTPUT_DIR}/HRC_formatted/{DATASET}_GSA_updated-updated-chr{{chr}}.vcf", chr=CHR)
    output:
        expand(f"{OUTPUT_DIR}/checkVCF/check-chr{{chr}}", chr=CHR)
    shell:
        """
        export LD_LIBRARY_PATH={OPENSSL_DIR}:$LD_LIBRARY_PATH

        for chr in {CHR}; do \
            {PYTHON2_DIR}/python2.7 scripts/checkVCF.py \
            -r {input.ref} \
            -o {OUTPUT_DIR}/checkVCF/check-chr${{chr}} \
            {OUTPUT_DIR}/HRC_formatted/{DATASET}_GSA_updated-updated-chr${{chr}}.vcf
        done
        """

rule bgzip_vcf:
    input:
        expand(f"{OUTPUT_DIR}/HRC_formatted/{DATASET}_GSA_updated-updated-chr{{chr}}.vcf", chr=CHR)
    output:
        expand(f"{OUTPUT_DIR}/VCF_compressed/{DATASET}_final-chr{{chr}}.vcf.gz", chr=CHR)
    shell:
        """
        module load BCFtools/1.22-GCCcore-13.3.0

        for chr in {CHR}; do \
            bcftools view {input} -Oz -o {output}; \
        done
        """