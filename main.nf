#!/usr/bin/env nextflow

if (params.data) {
    Channel.fromPath(params.data)
    .ifEmpty { exit 1, "FAM file (w/ header) containing phenotype data not found: ${params.data}" }
    .set { data }
}
if (params.vcf_file) {
    Channel.fromPath(params.vcf_file)
           .ifEmpty { exit 1, "VCF file containing  not found: ${params.vcf_file}" }
           .into { vcf_file; vcfs_to_split }
    vcfs_to_split
        .splitCsv(header: true)
        .map{ row -> [file(row.vcf)] }
        .set { vcfs }
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
    .ifEmpty { exit 1, "Phenotype file not found: ${params.pheno}" }
    .set { pheno }
Channel.fromPath(params.mapping)
    .ifEmpty { exit 1, "Mapping file not found: ${params.mapping}" }
    .set { mapping }

// Get as many processors as machine has
int threads = Runtime.getRuntime().availableProcessors()

if (params.vcf_file) {
    process file_preprocessing {
        publishDir 'results'
        container 'lifebitai/preprocess_gwas:latest'

        input:
        file vcfs from vcfs.collect()
        file vcf_file from vcf_file

        output:
        file 'merged.vcf' into vcf_plink
        file 'sample.phe' into data

        script:
        """
        # iterate through urls in csv replacing s3 path with the local one
        urls="\$(tail -n+2 $vcf_file | awk -F',' '{print \$2}')"
        for url in \$(echo \$urls); do
            vcf="\${url##*/}"
            sed -i -e "s~\$url~\$vcf~g" $vcf_file
        done
        # bgzip uncompressed vcfs
        for vcf in \$(tail -n+2 $vcf_file | awk -F',' '{print \$2}'); do
            if [ \${vcf: -4} == ".vcf" ]; then
                    bgzip -c \$vcf > \${vcf}.gz
                    sed -i "s/\$vcf/\${vcf}.gz/g" $vcf_file 
            fi
        done
        # remove any prexisting columns for sex 
        if grep -Fq "SEX" $vcf_file; then
            awk -F, -v OFS=, 'NR==1{for (i=1;i<=NF;i++)if (\$i=="SEX"){n=i-1;m=NF-(i==NF)}} {for(i=1;i<=NF;i+=1+(i==n))printf "%s%s",\$i,i==m?ORS:OFS}' $vcf_file > tmp.csv && mv tmp.csv $vcf_file
        fi
        # determine sex of each individual from VCF file & add to csv file
        echo 'SEX' > sex.txt
        for vcf in \$(tail -n+2 $vcf_file | awk -F',' '{print \$2}'); do
            bcftools index -f \$vcf
            SEX="\$(bcftools plugin vcf2sex \$vcf)"
            if [[ \$SEX == *M ]]; then
                    echo "1" >> sex.txt
            elif [ \$SEX == *F ]]; then
                    echo "2" >> sex.txt
            fi
        done
        # make fam file & merge vcfs
        paste -d, sex.txt $vcf_file > tmp.csv && mv tmp.csv $vcf_file
        make_fam2.py $vcf_file
        vcfs=\$(tail -n+2 $vcf_file | awk -F',' '{print \$3}')
        bcftools merge --force-samples \$vcfs > merged.vcf
        """
    }

    process plink {
        publishDir "${params.outdir}/plink", mode: 'copy'
        container 'alliecreason/plink:1.90'

        input:
        file vcf from vcf_plink
        file fam from data

        output:
        set file('*.bed'), file('*.bim'), file('*.fam') into plink

        script:
        """
        sed '1d' $fam > tmpfile; mv tmpfile $fam
        # remove contigs eg GL000229.1 to prevent errors
        sed -i '/^GL/ d' $vcf
        plink --vcf $vcf --make-bed
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
//     publishDir "${params.outdir}/filter", mode: 'copy'

//     input:
//     set file(bed), file(bim), file(fam) from plink

//     output:
//     file('*') into results

//     script:
//     """
//     plink --bfile plink --mind 0.1 --geno 0.1 --maf 0.05 --hwe 0.000001 --me 0.05 0.1 --tdt --ci 0.95 --out results1
//     """
// }

process recode {
    publishDir "${params.outdir}/plink", mode: 'copy'

    input:
    set file(bed), file(bim), file(fam) from plink
    file snps from snps

    output:
    file('*') into phewas

    script:
    extract = params.snps.contains("no_snps.txt") ? "" : "--extract $snps"
    """
    plink --recodeA --bfile ${bed.baseName} --out r_genotypes $extract
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
    phewas.R $pheno ${task.cpus} $params.pheno_codes
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