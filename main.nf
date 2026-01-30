nextflow.enable.dsl=2

chr_ch = Channel.from(1..22)

workflow {

    /*
     * Step 1: finalreport → lgen/map/fam
     */
    lgen_ch = CREATE_LGEN(
        file(params.input.final_report),
        file(params.input.sample_info)
    )

    /*
     * Step 2: lgen → plink bed
     */
    plink_ch = LGEN_TO_PLINK(lgen_ch)

    /*
     * Step 3: update build
     */
    updated_ch = UPDATE_BUILD(plink_ch)

    /*
     * Step 4: plink freq
     */
    freq_ch = PLINK_FREQ(updated_ch)

    /*
     * Step 5: HRC formatting
     */
    hrc_ch = MAKE_HRC_VCF(freq_ch)

    /* # Create (prefix, chr, vcf) */
    hrc_by_chr = hrc_ch.map { prefix, vcf ->
        def m = (vcf.name =~ /chr(\d+)/)
        assert m
        tuple(prefix, m[0][1] as int, vcf)
    }

    /* # Restrict explicitly to chromosomes 1–22 */
    hrc_by_chr = hrc_by_chr.filter { prefix, chr, vcf ->
        chr in 1..22

    /*
     * Step 6: checkVCF (parallel per chr)
     */
    CHECK_VCF(hrc_by_chr)

    /*
     * Step 7: bgzip (parallel per chr)
     */
    BGZIP_VCF(hrc_by_chr)
}

process CREATE_LGEN {

    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path final_report
    path sample_info

    output:
    tuple val(params.dataset),
          path("${params.dataset}.lgen"),
          path("${params.dataset}.map"),
          path("${params.dataset}.fam")

    script:
    """
    module load RPlus/4.2.1-foss-2022a-v22.12.1

    Rscript ${projectDir}/scripts/create_lgen.R \
        $final_report \
        $sample_info \
        ${params.dataset}.lgen \
        ${params.dataset}.map \
        ${params.dataset}.fam
    """
}

process LGEN_TO_PLINK {

    memory '16 GB'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple val(prefix),
          path(lgen),
          path(map),
          path(fam)

    output:
    tuple val(prefix),
          path("${prefix}.bed"),
          path("${prefix}.bim"),
          path("${prefix}.fam")

    script:
    """
    module load PLINK/1.9-beta6-20190617

    plink --lfile ${prefix} \
          --make-bed \
          --out ${prefix} \
          --memory 16000
    """
}


process UPDATE_BUILD {
    publishDir "${params.output_dir}/GSA_updated", mode: 'copy'

    input:
    tuple val(prefix),
          path(bed),
          path(bim),
          path(fam)

    output:
    tuple val(prefix),
          path("${prefix}_GSA_updated.bed"),
          path("${prefix}_GSA_updated.bim"),
          path("${prefix}_GSA_updated.fam")

    script:
    """
    module load PLINK/1.9-beta6-20190617

    ${projectDir}/scripts/update_build.sh \
        ${prefix} \
        /groups/umcg-immunogenetics/tmp02/users/NathanRibeiro/tools/GSA_strand_v3/GSAMD-24v3-0-EA_20034606_A1-b37.strand \
        ${prefix}_GSA_updated
    """
}

process PLINK_FREQ {

    publishDir "${params.output_dir}/GSA_updated", mode: 'copy'

    input:
    tuple val(prefix),
          path(bed),
          path(bim),
          path(fam)

    output:
    tuple val(prefix),
          path(bim),
          path("${prefix}.frq"),
          path("${prefix}.log")

    script:
    """
    module load PLINK/1.9-beta6-20190617

    plink \
      --bfile ${prefix} \
      --freq \
      --out ${prefix}
    """
}

process MAKE_HRC_VCF {

    publishDir "${params.output_dir}/HRC_formatted", mode: 'copy'

    input:
    tuple val(prefix),
          path(bim),
          path(frq),
          path(log)

    output:
    tuple val(prefix),
          path("${prefix}-updated-chr*.vcf")

    script:
    """
    module load foss/2022a
    module load PerlPlus/5.34.1-GCCcore-11.3.0-v22.11.1
    module load PLINK/1.9-beta6-20190617

    perl ${projectDir}/scripts/HRC-1000G-check-bim-NoReadKey.pl \
        -b ${bim} \
        -f ${frq} \
        -r ${params.HRC_reference} \
        -h

    sh Run-plink.sh
    """
}

process CHECK_VCF {

    publishDir "${params.output_dir}/checkVCF", mode: 'copy'

    input:
    tuple val(prefix),
          val(chr),
          path(vcf)

    output:
    path "check-chr${chr}*"

    script:
    """
    export LD_LIBRARY_PATH=${params.openssl_path}:\$LD_LIBRARY_PATH

    ${params.python2_path}/python2.7 ${projectDir}/scripts/checkVCF.py \
        -r ${params.checkvcf_ref} \
        -o check-chr${chr} \
        ${vcf}
    """
}

process BGZIP_VCF {

    publishDir "${params.output_dir}/VCF_compressed", mode: 'copy'

    input:
    tuple val(prefix),
          val(chr),
          path(vcf)

    output:
    path "${prefix}_final-chr${chr}.vcf.gz"

    script:
    """
    module load BCFtools/1.22-GCCcore-13.3.0

    bcftools view \
        ${vcf} \
        -Oz \
        -o ${prefix}_final-chr${chr}.vcf.gz
    """
}