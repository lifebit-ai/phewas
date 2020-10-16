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

ch_input_cb_data = params.phenofile ? Channel.value(params.phenofile) : Channel.empty()
ch_input_meta_data = params.metadata ? Channel.value(params.metadata) : Channel.empty()
gwas_input_ch = params.gwas_input ? Channel.value(params.gwas_input) : Channel.empty()

if (params.plink_input){
Channel
  .fromFilePairs("${params.plink_input}",size:3, flat : true)
  .ifEmpty { exit 1, "PLINK files not found: ${params.plink_input}" }
  .set { plinkCh }
}

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

if (params.mapping){
Channel.fromPath(params.mapping)
    .ifEmpty { exit 1, "Mapping file not found: ${params.mapping}" }
    .set { mapping }
}


// Get as many processors as machine has
int threads = Runtime.getRuntime().availableProcessors()

/*--------------------------------------------------
  Using files outside CB
---------------------------------------------------*/

if (params.vcf_file && !params.phenofile && !params.metadata) {
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
}

if (params.bed && params.bim && params.data && !params.phenofile && !params.metadata) {
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

if (!params.phenofile && !params.metadata){

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

    output:
    set file("*phewas_results.csv") into results_chr

    script:
    """
    mkdir -p assets/
    cp /assets/* assets/
    phewas.R --pheno_file "$pheno" --geno_file "$genotypes" --n_cpus ${task.cpus} --pheno_codes "$params.pheno_codes"
    """
    }
}


if (params.plink_input && params.phenofile && params.metadata) {

    process transform_cb_output {
    tag "$name"
    publishDir "${params.outdir}/design_matrix", mode: 'copy'

    input:
    val input_cb_data from ch_input_cb_data
    val input_meta_data from ch_input_meta_data

    output:
    file("*.json") into ch_encoding_json
    file("*id_code_count.csv") into codes_pheno
    file("*.phe") into pheno_file_prep

    script:
    """
    
    transform_cb_output.R --input_cb_data "${params.phenofile}" \
                          --input_meta_data "${params.metadata}" \
                          --phenoCol "${params.pheno_col}" \
                          --continuous_var_transformation "${params.continuous_var_transformation}" \
                          --continuous_var_aggregation "${params.continuous_var_aggregation}" \
                          --outdir "." \
                          --outprefix "${params.output_tag}"
    """
    }


    if (params.phenofile && params.case_group && params.design_mode == 'case_vs_control_contrast') {

        process add_design_matrix_case_vs_control_contrast {
        tag "$name"
        publishDir "${params.outdir}/contrasts", mode: 'copy'

        input:
        file(pheFile) from pheno_file_prep
        file(json) from ch_encoding_json

        output:
        file("${params.output_tag}_design_matrix_control_*.phe") into (pheno, pheno2)

        script:
        """

        create_design.R --input_file ${pheFile} \
                        --mode "${params.design_mode}" \
                        --case_group "${params.case_group}" \
                        --outdir . \
                        --outprefix "${params.output_tag}" \
                        --phenoCol "${params.pheno_col}"
                        
        """
        }
    }

    if (params.phenofile && params.case_group && params.design_mode == 'case_vs_groups_contrasts') {
        process add_design_matrix_case_vs_groups_contrasts {
            tag "$name"
            publishDir "${params.outdir}/contrasts", mode: 'copy'

            input:
            file(pheFile) from pheno_file_prep
            file(json) from ch_encoding_json

            output:
            file("${output_tag}_design_matrix_control_*.phe'") into (pheno, pheno2)

            script:
            """

            create_design.R --input_file ${pheFile} \
                            --case_group "${params.case_group}" \
                            --outdir . \
                            --outprefix "${params.output_tag}" \
                            --phenoCol "${params.pheno_col}"
                            
            """
        }
    }

    if (params.phenofile && params.design_mode == 'all_contrasts') {

        process add_design_matrix_all_contrasts {
            tag "$name"
            publishDir "${params.outdir}/contrasts", mode: 'copy'

            input:
            file(pheFile) from pheno_file_prep
            file(json) from ch_encoding_json

            output:
            file("${output_tag}_design_matrix_control_*.phe'") into (pheno, pheno2)

            script:
            """

            create_design.R --input_file ${pheFile} \
                            --mode ${params.design_mode}
                            --outdir . \
                            --outprefix "${params.output_tag}" \
                            --phenoCol "${params.pheno_col}"
                            
            """
        }
    }
    process preprocess_plink_cb {

        input:
        set val(name), file(bed), file(bim), file(fam) from plinkCh

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

    if (!params.snps) {
        process get_snps_cb {
            publishDir 'results', mode: 'copy'
            container 'alliecreason/plink:1.90'

            input:
            set file(bed), file(bim), file(fam) from plink
            file pheno_file from pheno2

            output:
            file("snps.txt") into snps

            script:
            """
            plink --bed $bed --bim $bim --fam $fam --pheno $pheno_file --pheno-name PHE --threads ${task.cpus} --assoc --allow-no-sex --out out
            awk -F' ' '{if(\$9<${params.snp_threshold}) print \$2}' out.assoc > snps.txt
            """
        }
    }

    process recode_cb {
    publishDir "${params.outdir}/plink", mode: 'copy'

    input:
    set file(bed), file(bim), file(fam) from plink2
    file snps from snps

    output:
    file('*.raw') into phewas

    script:
    """
    plink --recodeA --bfile ${bed.baseName} --out r_genotypes --extract $snps --allow-no-vars
    """
    }

    process phewas_cb {
    publishDir "${params.outdir}/phewas", mode: 'copy'
    cpus threads

    input:
    file genotypes from phewas
    file pheno from codes_pheno

    output:
    file("*phewas_results.csv") into results_chr

    script:
    """
    mkdir -p assets/
    cp /assets/* assets/
    phewas.R --pheno_file "$pheno" --geno_file "$genotypes" --n_cpus ${task.cpus} --pheno_codes "$params.pheno_codes"
    """
    }

    process merge_results {
    publishDir "${params.outdir}/merged_results", mode: 'copy'

    input:
    file("*phewas_result.csv") from results_chr.collect()

    output:
    set file("merged_results.csv"), file("merged_top_results.csv"), file("*png") into plots, plots2

    script:
    """
    plot_merged_results.R

    """
    }
}


/*---------------------------------
  Colocalization analysis
-----------------------------------*/

if (params.post_analysis == 'coloc'){

    process run_coloc {
        publishDir "${params.outdir}/colocalization", mode: "copy"
        
        input:
        file gwas_file from gwas_input_ch
        set file(merged_results), file(merged_top_results), file("*png") from plots2

        output:
        set file("*coloc_heatmap.png"), file("*coloc_results.csv") into coloc_results_ch

        script:
        """
        cp /opt/bin/* .
        run_coloc.R --phewas_summary "$merged_results" \
                    --gwas_summary "${params.gwas_input}" \
                    --gwas_trait_type "${params.gwas_trait_type}" \
                    --outprefix "${params.output_tag}"
        """
    }
    process build_report_coloc {
        tag "report"
        publishDir "${params.outdir}/MultiQC", mode: 'copy', pattern: '*.html'

        input:
        set file(coloc_plot), file(coloc_results) from coloc_results_ch
        set file(phewas_results), file(phewas_top_results), file(phewas_plot) from plots

        output:
        file("multiqc_report.html") into ch_report_outputs

        script:

        """
        mkdir assets/
        cp /assets/* assets/

        
        # Generates the report
        cp /opt/bin/phewas_report.Rmd .
        cp /opt/bin/DTable.R .
        cp /opt/bin/sanitise.R .
        cp /opt/bin/style.css .
        cp /opt/bin/logo.png .
        

        Rscript -e "rmarkdown::render('phewas_report.Rmd', params = list(phewas_manhattan='${phewas_plot}', phewas_results='${phewas_results}', coloc_results='${coloc_results}', coloc_heatmap='${coloc_plot}'))"
        mv phewas_report.html multiqc_report.html

        rm ./DTable.R
        rm ./sanitise.R
        rm ./style.css
        rm ./phewas_report.Rmd

        """
    }
}

if (!params.post_analysis){

    process build_report {
        tag "report"
        publishDir "${params.outdir}/MultiQC", mode: 'copy', pattern: '*.html'

        input:
        set file(phewas_results), file(phewas_top_results), file(phewas_plot) from plots

        output:
        file("multiqc_report.html") into ch_report_outputs

        script:

        """
        mkdir assets/
        cp /assets/* assets/

        
        # Generates the report
        cp /opt/bin/phewas_report.Rmd .
        cp /opt/bin/DTable.R .
        cp /opt/bin/sanitise.R .
        cp /opt/bin/style.css .
        cp /opt/bin/logo.png .
        

        Rscript -e "rmarkdown::render('phewas_report.Rmd', params = list(phewas_manhattan='${phewas_manhattan}', phewas_results='${phewas_results}', coloc_results='None', coloc_heatmap='None'))"
        mv phewas_report.html multiqc_report.html

        rm ./DTable.R
        rm ./sanitise.R
        rm ./style.css
        rm ./phewas_report.Rmd
        """
    }

}


