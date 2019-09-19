#!/usr/bin/env nextflow

Channel.fromPath(params.data)
    .ifEmpty { exit 1, "FAM file (w/ header) containing phenotype data not found: ${params.data}" }
    .set { data }
if (params.vcf) {
    Channel.fromPath(params.vcf)
    .ifEmpty { exit 1, "VCF file not found: ${params.vcf}" }
    .set { vcf }
}
if (params.bed) {
    Channel.fromPath(params.bed)
        .ifEmpty { exit 1, "PLINK binary pedigree file not found: ${params.bed}" }
        .set { bed }
}
if (params.bim) {
    Channel.fromPath(params.bim)
        .ifEmpty { exit 1, "PLINK BIM file not found: ${params.bim}" }
        .set { bim }
}
Channel.fromPath(params.snps)
    .ifEmpty { exit 1, "SNPs of interest file not found: ${params.snps}" }
    .set { snps }
Channel.fromPath(params.pheno)
    .ifEmpty { exit 1, "SNPs of interest file not found: ${params.pheno}" }
    .set { pheno }
Channel.fromPath(params.mapping)
    .ifEmpty { exit 1, "Mapping file not found: ${params.mapping}" }
    .set { mapping }

// Get as many processors as machine has
int threads = Runtime.getRuntime().availableProcessors()

if (params.data && params.vcf) {
    process vcf2plink {
    publishDir "${params.outdir}/vcf2plink", mode: 'copy'

    input:
    file vcf from vcf
    file fam from data

    output:
    set file('*.bed'), file('*.bim'), file('*.fam') into plink

    script:
    """
    sed '1d' $fam > tmpfile; mv tmpfile $fam
    plink --vcf $vcf
    rm plink.fam
    mv $fam plink.fam
    """
    }
} else if (params.bed && params.bim && params.data) {
    process preprocess_plink {

    input:
    file bed from bed
    file bim from bim
    file fam from data

    output:
    set file('*.bed'), file('*.bim'), file('*.fam') into plink

    script:
    """
    sed '1d' $fam > tmpfile; mv tmpfile $fam
    mv $fam plink.fam
    mv $bed plink.bed
    mv $bim plink.bim
    """
    }
} else {
    exit 1, "\nPlease specify either:\n1) `--vcf` AND `--data` inputs\nOR\n2) `--bed` AND `--bim` AND `--data` inputs"
}

// process filter {
//     publishDir "${params.output_dir}/filter", mode: 'copy'

//     input:
//     set file(bed), file(bim), file(fam) from plink

//     output:
//     file('*') into results

//     script:
//     """
//     plink --bfile plink --mind 0.1 --geno 0.1 --maf 0.05 --hwe 0.000001 --me 0.05 0.1 --tdt --ci 0.95 --out results1
//     """
// }


process plink {
    publishDir "${params.outdir}/plink", mode: 'copy'

    input:
    set file(bed), file(bim), file(fam) from plink
    file snps from snps

    output:
    file('*') into phewas

    script:
    """
    plink --recodeA --bfile plink --extract $snps --out r_genotypes
    """
}

process phewas {
    publishDir "${params.outdir}/phewas", mode: 'copy'
    cpus threads

    input:
    file genotypes from phewas
    file pheno from pheno
    file mapping from mapping

    output:
    set file("phewas_results.csv"), file("top_results.csv"), file("*.png") into plots

    script:
    """
    phewas.R $pheno ${task.cpus} $pheno_codes
    """
}


process visualisations {
    publishDir "${params.outdir}/Visualisations", mode: 'copy'

    container 'lifebitai/vizjson:latest'

    input:
    set file(phe), file(top), file(man) from plots

    output:
    file '.report.json' into viz

    script:
    """
    img2json.py "${params.outdir}/phewas/$man" "Phenotype Manhattan Plot" ${man}.json  
    csv2json.py $top "Top results from PheWAS by significance" ${top}.json
    combine_reports.py .
    """
}