nextflow.enable.dsl=2

Channel
    .of(params.chromosomes)
    .set { chr_ch }

workflow {

    /*
     * Step 1: finalreport → lgen/map/fam
     */
    lgen_files = CREATE_LGEN(
        file(params.input.final_report),
        file(params.input.sample_info)
    )

    /*
     * Step 2: lgen → plink bed
     */
    bed_prefix = LGEN_TO_PLINK(lgen_files)

    /*
     * Step 3: update build
     */
    updated_prefix = UPDATE_BUILD(bed_prefix)

    /*
     * Step 4: plink freq
     */
    freq_prefix = PLINK_FREQ(updated_prefix)

    /*
     * Step 5: HRC formatting
     */
    hrc_vcf = MAKE_HRC_VCF(freq_prefix)

    /*
     * Step 6: checkVCF (parallel per chr)
     */
    checked = CHECK_VCF(chr_ch, hrc_vcf.vcf)

    /*
     * Step 7: bgzip (parallel per chr)
     */
    BGZIP_VCF(chr_ch, hrc_vcf.vcf)
}

process CREATE_LGEN {

    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path final_report
    path sample_info

    output:
    path "${params.dataset}.lgen"
    path "${params.dataset}.map"
    path "${params.dataset}.fam"

    script:
    """
    module load RPlus/4.2.1-foss-2022a-v22.12.1

    Rscript scripts/create_lgen.R \
        $final_report \
        $sample_info \
        ${params.output_dir}/${params.dataset}.lgen \
        ${params.output_dir}/${params.dataset}.map \
        ${params.output_dir}/${params.dataset}.fam
    """
}

process LGEN_TO_PLINK {

    input:
    path lgen_files

    output:
    path "${params.dataset}.*"

    script:
    """
    module load PLINK/1.9-beta6-20190617

    plink --lfile ${params.dataset} \
          --autosome \
          --make-bed \
          --out ${params.dataset}
    """
}

process UPDATE_BUILD {

    publishDir "${params.output_dir}/GSA_updated", mode: 'copy'

    input:
    path bed_prefix

    output:
    path "${params.dataset}_GSA_updated.*"

    script:
    """
    module load PLINK/1.9-beta6-20190617

    scripts/update_build.sh \
        ${bed_prefix.baseName} \
        /groups/umcg-immunogenetics/tmp02/users/NathanRibeiro/tools/GSA_strand_v3/GSAMD-24v3-0-EA_20034606_A1-b37.strand \
        ${params.dataset}_GSA_updated
    """
}

process PLINK_FREQ {

    publishDir "${params.output_dir}/GSA_updated", mode: 'copy'

    input:
    path updated_prefix

    output:
    path "${params.dataset}_GSA_updated.frq*"

    script:
    """
    module load PLINK/1.9-beta6-20190617

    plink --freq \
          --bfile ${updated_prefix.baseName} \
          --out ${params.dataset}_GSA_updated
    """
}

process MAKE_HRC_VCF {

    publishDir "${params.output_dir}/HRC_formatted", mode: 'copy', pattern: "*-updated-chr*.vcf"

    input:
    path freq_prefix

    output:
    path "${params.dataset}_GSA_updated-updated-chr*.vcf", emit: vcf

    script:
    """
    module load foss/2022a
    module load PerlPlus/5.34.1-GCCcore-11.3.0-v22.11.1
    module load PLINK/1.9-beta6-20190617

    perl ${params.HRC_script} \
        -b ${freq_prefix.parent}/${params.dataset}_GSA_updated.bim \
        -f ${freq_prefix.parent}/${params.dataset}_GSA_updated.frq \
        -r ${params.HRC_reference} \
        -h

    sh Run-plink.sh
    """
}

process CHECK_VCF {

    publishDir "${params.output_dir}/checkVCF", mode: 'copy'

    input:
    val chr
    path hrc_vcf

    output:
    path "check-chr${chr}*"

    script:
    """
    export LD_LIBRARY_PATH=${params.openssl_path}:\$LD_LIBRARY_PATH

    ${params.python2_path}/python2.7 scripts/checkVCF.py \
        -r ${params.checkvcf_ref} \
        -o check-chr${chr} \
        ${hrc_vcf}
    """
}

process BGZIP_VCF {

    publishDir "${params.output_dir}/VCF_compressed", mode: 'copy'

    input:
    val chr
    path hrc_vcf

    output:
    path "${params.dataset}_final-chr${chr}.vcf.gz"

    script:
    """
    module load BCFtools/1.22-GCCcore-13.3.0

    bcftools view \
        ${hrc_vcf} \
        -Oz \
        -o ${params.dataset}_final-chr${chr}.vcf.gz
    """
}