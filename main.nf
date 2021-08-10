#!/usr/bin/env nextflow

/*
========================================================================================
                         lifebit-ai/phewas
========================================================================================
 lifebit-ai/phewas pheWAS pipeline
 #### Homepage / Documentation
 https://github.com/lifebit-ai/phewas
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info """
    phewas - A pipeline for running pheWAS adapted for vcf and plink inputs.
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --plink_input /path/plink.{bed,bim,fam} –-input_phenofile pheno.phe  \
    –-input_id_code_count icd10_id_code_count.csv –-pheno_codes "icd10" [Options]
    
    Essential parameters:

    Genomic data

    –-plink_input : Path/URL to plink bim bed fam files, in format /path/plink.{bed,bim,fam}

    OR:

    –-bim : Path/URL to bim file.

    –-bed : Path/URL to bed file.

    –-fam : Path/URL to fam file.

    OR:

    –-agg_vcf_file : Path/URL to .csv file containing chr chunk information, path to aggregated VCFs, VCFs index. Columns must include chr,vcf,index.

    OR:

    –-individual_vcf_file : Path/URL to .csv file containing individual, path to individual VCFs.


    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

// Header log info

def summary = [:]

if (workflow.revision) summary['Pipeline Release'] = workflow.revision

summary['Output dir']                                  = params.outdir
summary['Launch dir']                                  = workflow.launchDir
summary['Working dir']                                 = workflow.workDir
summary['Script dir']                                  = workflow.projectDir
summary['User']                                        = workflow.userName

summary['input_phenofile']                             = params.input_phenofile
summary['input_id_code_count']                         = params.input_id_code_count

summary['individual_vcf_file']                         = params.individual_vcf_file
summary['agg_vcf_file']                = params.agg_vcf_file
summary['plink_input'] = params.plink_input
summary['fam']                 = params.fam
summary['bed']                 = params.bed
summary['bim']                 = params.bim


summary['snps'] =  params.snps
summary['snp_threshold'] = params.snp_threshold
summary['pheno_file'] = params.pheno_file
summary['pheno_codes'] = params.pheno_codes
summary['post_analysis'] = params.post_analysis
summary['gwas_input'] = params.gwas_input
summary['gwas_trait_type'] = params.gwas_trait_type




log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

/*--------------------------------------------------
  Channel setup
---------------------------------------------------*/

if (!params.agg_vcf_file && !params.individual_vcf_file && !params.plink_input && !params.fam && !params.bed && !params.bim) {
  exit 1, "You have not supplied any genotypic data.\
  \nPlease use either --plink_input/--bed,--bim,--fam OR --agg_vcf_file/--individual_vcf_file to provide genotypic data."
}

if (params.plink_input && params.fam && params.bed && params.bim) {
  exit 1, "Only one set of plink files is supported as input.\
  \nPlease use either --plink_input to specify path to a set of plink files or --bed,--bim,--fam to provide plink files individually."
}

if (params.agg_vcf_file && params.individual_vcf_file) {
    exit 1, "Only one set of vcf files is supported as input.\
    \nPlease use either --agg_vcf_file or --individual_vcf_file to supply multi-sample or per-sample VCFs respectively."
}

if (params.input_phenofile){
    Channel.fromPath(params.input_phenofile)
           .ifEmpty { "Phenotype file not found" }
           .into{ ch_pheno; ch_pheno2; ch_pheno3 }
}

ch_codes_pheno = params.input_id_code_count ? Channel.value(file(params.input_id_code_count)) : Channel.empty()
ch_gwas_input = params.gwas_input ? Channel.value(file(params.gwas_input)) : Channel.empty()

if (params.agg_vcf_file){
    Channel.fromPath(params.agg_vcf_file)
           .ifEmpty { exit 1, "VCF file containing  not found: ${params.agg_vcf_file}" }
           .into {ch_vcf_file; ch_vcfs_to_split; ch_index_to_split}
    ch_vcfs_to_split
        .splitCsv(header: true)
        .map{ row -> [file(row.vcf)] }
        .set { ch_vcfs }
    ch_index_to_split
        .splitCsv(header: true)
        .map{ row -> [file(row.index)] }
        .set { ch_indexes }
}

if (params.plink_input){
    Channel
    .fromFilePairs("${params.plink_input}",size:3, flat : true)
    .ifEmpty { exit 1, "PLINK files not found: ${params.plink_input}" }
    .set { plinkCh }
}

if (params.individual_vcf_file) {
    Channel.fromPath(params.individual_vcf_file)
           .ifEmpty { exit 1, "VCF file containing  not found: ${params.individual_vcf_file}" }
           .into { ch_vcf_file; ch_vcfs_to_split }
    ch_vcfs_to_split
        .splitCsv(header: true)
        .map{ row -> [file(row.vcf)] }
        .into { ch_vcfs; ch_vcf_ind }
}


if (params.fam) {
    Channel.fromPath(params.fam)
    .ifEmpty { exit 1, "FAM file (w/ header) containing phenotype info not found: ${params.fam}" }
    .set { ch_fam }
}
if (params.bed) {
    Channel.fromPath(params.bed)
        .ifEmpty { exit 1, "PLINK binary pedigree file not found: ${params.bed}" }
        .set { ch_bed }
}
if (params.bim) {
    Channel.fromPath(params.bim)
        .ifEmpty { exit 1, "PLINK BIM file not found: ${params.bim}" }
        .set { ch_bim }
}
if (params.snps) {
    Channel.fromPath(params.snps)
    .ifEmpty { exit 1, "SNPs of interest file not found: ${params.snps}" }
    .set { ch_snps }
}


// Get as many processors as machine has
int threads = Runtime.getRuntime().availableProcessors()

/*--------------------------------------------------
  Using different formats: vcf & plink
---------------------------------------------------*/
if (params.agg_vcf_file){
    process merge_agg_vcfs {

        input:
        file vcf_file from ch_vcf_file
        
 
        output:
        
        file 'vcf_files.txt' into ch_updated_vcf_list

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
        tail -n+2 $vcf_file| awk -F',' '{print \$2}' > vcf_files.txt
        """
    }
}

if (params.individual_vcf_file) {
    process merge_ind_vcfs {

        label 'file_preprocessing'

        input:
        file vcf_file from ch_vcf_file
        file vcfs from ch_vcf_ind.collect()

        output:
        file 'vcf_files.txt' into ch_updated_vcf_list

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
        paste -d, sex.txt $vcf_file > tmp.csv && mv tmp.csv $vcf_file
        tail -n+2 $vcf_file | awk -F',' '{print \$3}' > vcf_files.txt
        """
    }

}

if (params.agg_vcf_file || params.individual_vcf_file){

    process combine_vcfs {
        publishDir "${params.outdir}/vcf", mode: 'copy'

        input:
        file(vcfs) from ch_vcfs.collect()
        file vcf_list from  ch_updated_vcf_list
        file pheno_file from ch_pheno3

        output:
        file 'filtered_by_sample.vcf.gz' into ch_vcf_plink

        script:
        if ( params.concat_vcfs )
            """
            for vcf in ${vcfs}; do
                if [ \${vcf} == *vcf ]; then
                    bgzip -c \$vcf > \${vcf}.gz
                fi
            done
            for i in *.vcf.gz; do bcftools index \${i}; done
            vcfs_to_combine=\$(find . -name '*.vcf.gz'| paste -sd " ")
            sed '1d' $pheno_file | awk -F' ' '{print \$1}' > sample_file.txt
            bcftools concat \${vcfs_to_combine} -Oz -o merged.vcf.gz
            bcftools view -S sample_file.txt merged.vcf.gz -Oz -o filtered_by_sample.vcf.gz
            """
        else if ( !params.concat_vcfs )
            """
            for vcf in ${vcfs}; do
                if [ \${vcf} == *vcf ]; then
                    bgzip -c \$vcf > \${vcf}.gz
                fi
            done
            for i in *.vcf.gz; do bcftools index \${i}; done
            vcfs_to_combine=\$(find . -name '*.vcf.gz'| paste -sd " ")
            bcftools merge --force-samples \${vcfs_to_combine} -Oz -o merged.vcf.gz
            sed '1d' $pheno_file | awk -F' ' '{print \$1}' > sample_file.txt
            bcftools view -S sample_file.txt merged.vcf.gz -Oz -o filtered_by_sample.vcf.gz
            """
    }

    process vcf_2_plink {
        tag "plink"
        

        input:
        file vcf from ch_vcf_plink
        file fam from ch_pheno

        output:
        set file('*.bed'), file('*.bim'), file('*.fam') into ch_plink, ch_plink2

        script:
        """
        sed '1d' $fam > tmpfile; mv tmpfile $fam
        # remove contigs eg GL000229.1 to prevent errors
        gunzip $vcf -c > temp.vcf
        sed -i '/^GL/ d' temp.vcf
        plink --keep-allele-order \
        --vcf temp.vcf \
        --make-bed \
        --vcf-half-call m
        rm plink.fam
        mv $fam plink.fam
        """
    }
}

if (params.bed && params.bim && params.fam) {

    process preprocess_plink {

        input:
        file bed from ch_bed
        file bim from ch_bim
        file fam from ch_fam

        output:
        set file('*.bed'), file('*.bim'), file('*.fam') into ch_plink, ch_plink2

        script:
        """
        sed '1d' $fam > tmpfile; mv tmpfile $fam
        mv $fam plink.fam
        mv $bed plink.bed
        mv $bim plink.bim
        """
    }
}

if (params.plink_input) {

    process preprocess_plink_folder {

        input:
        set val(name), file(bed), file(bim), file(fam) from plinkCh

        output:
        set file('*.bed'), file('*.bim'), file('*.fam') into ch_plink, ch_plink2

        script:
        """
        sed '1d' $fam > tmpfile; mv tmpfile $fam
        mv $fam plink.fam
        mv $bed plink.bed
        mv $bim plink.bim
        """
    }
}

/*--------------------------------------------------
  Prepare resulting plink files
---------------------------------------------------*/

if (!params.snps) {
    process get_snps {
        tag "plink"

        input:
        set file(bed), file(bim), file(fam) from ch_plink
        file pheno_file from ch_pheno2

        output:
        file("snps.txt") into ch_snps

        script:
        """
        plink --keep-allele-order \
        --bed $bed \
        --bim $bim \
        --fam $fam \
        --pheno ${pheno_file} \
        --pheno-name PHE \
        --threads ${task.cpus} \
        --assoc \
        --out out
        awk -F' ' '{if(\$9<${params.snp_threshold}) print \$2}' out.assoc > snps.txt
        """
    }
}

process recode {
tag "plink"

input:
set file(bed), file(bim), file(fam) from ch_plink2
file snps from ch_snps

output:
file('*.raw') into phewas

script:
"""
plink --keep-allele-order \
--recodeA \
--bfile ${bed.baseName} \
--out r_genotypes \
--extract $snps
"""
}

/*--------------------------------------------------
  Run phewas
---------------------------------------------------*/

process phewas {
    cpus threads

    input:
    file genotypes from phewas
    file pheno from ch_codes_pheno

    output:
    file("*phewas_results.csv") into results_chr

    script:
    """
    mkdir -p assets/
    cp /assets/* assets/
    phewas.R --pheno_file "${pheno}" --geno_file "${genotypes}" --n_cpus ${task.cpus} --pheno_codes "$params.pheno_codes"
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
        

        Rscript -e "rmarkdown::render('phewas_report.Rmd', params = list(phewas_manhattan='${phewas_plot}', phewas_results='${phewas_results}', coloc_results='None', coloc_heatmap='None'))"
        mv phewas_report.html multiqc_report.html

        rm ./DTable.R
        rm ./sanitise.R
        rm ./style.css
        rm ./phewas_report.Rmd
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
        file gwas_file from ch_gwas_input
        set file(merged_results), file(merged_top_results), file("*png") from plots2

        output:
        set file("*coloc_heatmap.png"), file("*coloc_results.csv") into coloc_results_ch

        script:
        """
        run_coloc.R --phewas_summary "$merged_results" \
                    --gwas_summary "${gwas_file}" \
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

