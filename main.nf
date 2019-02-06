#!/usr/bin/env nextflow

Channel.fromPath(params.vcf)
    .ifEmpty { exit 1, "VCF file not found: ${params.vcf}" }
    .set { vcf }
Channel.fromPath(params.data)
    .ifEmpty { exit 1, "FAM file (w/ header) containing phenotype data not found: ${params.data}" }
    .set { data }
Channel.fromPath(params.snps)
    .ifEmpty { exit 1, "SNPs of interest file not found: ${params.snps}" }
    .set { snps }
Channel.fromPath(params.pheno)
    .ifEmpty { exit 1, "SNPs of interest file not found: ${params.pheno}" }
    .set { pheno }

// Get as many processors as machine has
int threads = Runtime.getRuntime().availableProcessors()

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

    output:
    file("*") into plots

    script:
    """
    phewas.R $pheno ${task.cpus}
    """
}


// process visualisations {
//     publishDir "${params.outdir}/Visualisations", mode: 'copy'

//     container 'lifebitai/vizjson:latest'

//     input:
//     set file(qq_dfam), file(qq_poo), file(qq_tdt), file(tdt) from plots

//     output:
//     file '.report.json' into viz

//     script:
//     """
//     for image in \$(ls *png); do
//         prefix="\${image%.*}"
//         if [[ \$prefix == "QQdfam" ]]; then
//             title="QQ Plot From The Sib-TDT Test To Show Family-based Disease Traits"
//         elif [ \$prefix == "QQpoo" ]; then
//             title="QQ Plot For The Parent of Orgin Analysis"
//         elif [ \$prefix == "QQtdt" ]; then
//             title="QQ Plot From The Transmission Disequilibrium Test To Show Family-based Associations"
//         elif [ \$prefix == "tdt" ]; then
//             title="Manhattan Plot From The Transmission Disequilibrium Test To Show Family-based Associations"
//         fi
//         img2json.py "${params.outdir}/plots/\$image" "\$title" \${prefix}.json
//     done
    
//     #table=\$(ls *.txt)
//     #prefix=\${table%.*}
//     #tsv2csv.py < \${prefix}.txt > \${prefix}.csv
//     #csv2json.py \${prefix}.csv "Combined Results" 'results-combined.json'
    
//     combine_reports.py .
//     """
// }