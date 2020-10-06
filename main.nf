#!/usr/bin/env nextflow

/*
========================================================================================
                         lifebit-ai/phewas
========================================================================================
 lifebit-ai/phewas pheWAS pipeline built for Genomics England and CloudOS
 #### Homepage / Documentation
 https://github.com/lifebit-ai/phewas
----------------------------------------------------------------------------------------
*/

/*--------------------------------------------------
  Channel setup
---------------------------------------------------*/

ch_input_cb_data = params.cohort_browser_phenofile ? Channel.value(params.cohort_browser_phenofile) : Channel.empty()
ch_input_meta_data = params.input_meta_data ? Channel.value(params.input_meta_data) : Channel.empty()



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
if (params.snps) {
    Channel.fromPath(params.snps)
    .ifEmpty { exit 1, "SNPs of interest file not found: ${params.snps}" }
    .set { snps }
}
if (params.pheno_file){
Channel.fromPath(params.pheno_file)
    .ifEmpty { exit 1, "Phenotype file not found: ${params.pheno_file}" }
    .set { pheno }
}
Channel.fromPath(params.mapping)
    .ifEmpty { exit 1, "Mapping file not found: ${params.mapping}" }
    .set { mapping }

// Get as many processors as machine has
int threads = Runtime.getRuntime().availableProcessors()

if (params.cohort_browser_phenofile && params.input_meta_data && params.vcf_file) {
    params.pheno = "PHE"
    process transform_cb_output{
    tag "$name"
    publishDir "${params.outdir}/design_matrix", mode: 'copy'

    input:
    val input_cb_data from ch_input_cb_data
    val input_meta_data from ch_input_meta_data
    val vcf_filename from vcf_file

    output:
    file("*.json") into ch_encoding_json
    file("*vcf_samples.csv") into vcf_file_prep
    file("*id_icd10_count.csv") into pheno

    script:
    """
    cp /opt/bin/* .

    mkdir -p ${params.outdir}/design_matrix
    
    transform_cb_output.R --input_cb_data "${params.cohort_browser_phenofile}" \
                          --input_meta_data "${params.input_meta_data}" \
                          --phenoCol "${params.phenoCol}" \
                          --vcf_list "${params.vcf_file}" \
                          --continuous_var_transformation "${params.continuous_var_transformation}" \
                          --outdir "." \
                          --outprefix "${params.output_tag}"
    """
    }
    if (params.cohort_browser_phenofile && params.case_group && params.mode == 'case_vs_control_contrast') {
    process add_design_matrix_case_vs_control_contrast{
        tag "$name"
        publishDir "${params.outdir}/contrasts", mode: 'copy'

        input:
        file(pheFile) from vcf_file_prep
        file(json) from ch_encoding_json

        output:
        file("${params.output_tag}_design_matrix_control_*.phe") into vcf_file_2

        script:
        """
        cp /opt/bin/* .

        mkdir -p ${params.outdir}/contrasts

        create_design.R --input_file ${pheFile} \
                        --mode "${params.mode}" \
                        --case_group "${params.case_group}" \
                        --outdir . \
                        --outprefix "${params.output_tag}" \
                        --phenoCol "${params.phenoCol}"
                        
        """
    }
    }

    if (params.cohort_browser_phenofile &&  params.case_group && params.mode == 'case_vs_groups_contrasts') {
    process add_design_matrix_case_vs_groups_contrasts{
        tag "$name"
        publishDir "${params.outdir}/contrasts", mode: 'copy'

        input:
        file(pheFile) from ch_transform_cb
        file(json) from ch_encoding_json

        output:
        file("${output_tag}_design_matrix_control_*.phe'") into vcf_file_2

        script:
        """
        cp /opt/bin/* .

        mkdir -p ${params.outdir}/contrasts

        create_design.R --input_file ${pheFile} \
                        --case_group "${params.case_group}" \
                        --outdir . \
                        --outprefix "${params.output_tag}" \
                        --phenoCol "${params.phenoCol}"
                        
        """
    }
    }

    if (params.cohort_browser_phenofile && params.mode == 'all_contrasts') {

    process add_design_matrix_all_contrasts{
        tag "$name"
        publishDir "${params.outdir}/contrasts", mode: 'copy'

        input:
        file(pheFile) from ch_transform_cb

        output:
        file("${output_tag}_design_matrix_control_*.phe'") into vcf_file_2

        script:
        """
        cp /opt/bin/* .

        mkdir -p ${params.outdir}/contrasts

        create_design.R --input_file ${pheFile} \
                        --mode ${params.mode}
                        --outdir . \
                        --outprefix "${params.output_tag}" \
                        --phenoCol "${params.phenoCol}"
                        
        """
    }
    }
    process file_preprocessing_cb {
        publishDir 'results'
        container 'lifebitai/preprocess_gwas:latest'

        input:
        file vcfs from vcfs.collect()
        file vcf_file from vcf_file_2

        output:
        file 'merged.vcf' into vcf_plink
        file 'sample.phe' into data, data2

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
        python make_fam.py $vcf_file
        vcfs=\$(tail -n+2 $vcf_file | awk -F',' '{print \$3}')
        bcftools merge --force-samples \$vcfs > merged.vcf
        """
    }

    process plink_cb {
        publishDir "${params.outdir}/plink", mode: 'copy'
        container 'alliecreason/plink:1.90'

        input:
        file vcf from vcf_plink
        file fam from data

        output:
        set file('*.bed'), file('*.bim'), file('*.fam') into plink, plink2

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
    if (!params.snps) {
    process get_snps_cb {
        publishDir 'results', mode: 'copy'
        container 'alliecreason/plink:1.90'

        input:
        set file(bed), file(bim), file(fam) from plink
        file pheno_file from data2

        output:
        file("snps.txt") into snps

        script:
        """
        plink --bed $bed --bim $bim --fam $fam --pheno $pheno_file --pheno-name PHE --threads ${task.cpus} --assoc --out out
        awk -F' ' '{if(\$9<${params.snp_threshold}) print \$2}' out.assoc > snps.txt
        """
    }
}
}

if (params.vcf_file && !params.cohort_browser_phenofile && !params.input_meta_data) {
    process file_preprocessing {
        publishDir 'results'
        container 'lifebitai/preprocess_gwas:latest'

        input:
        file vcfs from vcfs.collect()
        file vcf_file from vcf_file

        output:
        file 'merged.vcf' into vcf_plink
        file 'sample.phe' into data, data2

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
        python make_fam.py $vcf_file
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
        set file('*.bed'), file('*.bim'), file('*.fam') into plink, plink2

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
    if (!params.snps) {
    process get_snps {
        publishDir 'results', mode: 'copy'
        container 'alliecreason/plink:1.90'

        input:
        set file(bed), file(bim), file(fam) from plink
        file pheno_file from data2

        output:
        file("snps.txt") into snps

        script:
        """
        plink --bed $bed --bim $bim --fam $fam --pheno $pheno_file --pheno-name $params.pheno --threads ${task.cpus} --assoc --out out
        awk -F' ' '{if(\$9<${params.snp_threshold}) print \$2}' out.assoc > snps.txt
        """
    }
}
}

if (params.bed && params.bim && params.data) {
    process preprocess_plink {

        input:
        file bed from bed
        file bim from bim
        file fam from data

        output:
        set file('*.bed'), file('*.bim'), file('*.fam') into plink, plink2

        script:
        """
        sed '1d' $fam > tmpfile; mv tmpfile $fam
        mv $fam plink.fam
        mv $bed plink.bed
        mv $bim plink.bim
        """
    }
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
    set file(bed), file(bim), file(fam) from plink2
    file snps from snps

    output:
    file('*.raw') into phewas

    script:
    """
    plink --recodeA --bfile ${bed.baseName} --out r_genotypes --extract $snps
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
    mkdir -p assets/
    cp /assets/* assets/
    phewas.R --pheno_file "$pheno" --geno_file "$genotypes" --n_cpus ${task.cpus} --pheno_codes "$params.pheno_codes"
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
