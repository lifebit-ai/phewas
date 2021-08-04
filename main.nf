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
    nextflow run main.nf --input reads.csv [Options]
    
    Inputs Options:
    --input         A csv file with paths to fastq.gz or a list SRA accession.
                    [Required] (Default: None)
    Reference Options:
    --genome        A reference genome. Available - 'GRCh37' and 'GRCh38'
                    This options handles all other require reference and make sure the builds are up-to-date
                    [Required] (Default: 'GRCh38')
    --genome_fasta  A genome fasta file
                    [Optional if --genome provided]
    --anno_gtf      A annotation GTF file
                    [Optional if --genome provided]
    --star_index    STAR index. 
                    (Make sure the build done with STAR: v2.7.1a)
                    [Optional if --genome provided]
    --ctat_genome_lib CTAT genome library for STAR-Fusion
                    (Make sure the build done with STAR: v2.7.1a and STAR-Fusion: v1.7.0)
                    [Optional if --genome provided]
    
    Parabricks Options:
    --chim_enable_star_options  
                    All the options related to Chimeric junction out in STAR aligner
                    This also decides if the rna-fusion calling going to happen or not.
                    (Default: "--min-chim-segment 12 --min-chim-overhang 8 --out-chim-format 1")
    Additonal Options:
    --outdir        A output directory
                    [Default: 'results']
    --help          This help menu
    Resource Options:
    --max_cpus      Maximum number of CPUs (int)
                    (default: $params.max_cpus)  
    --max_memory    Maximum memory (memory unit)
                    (default: $params.max_memory)
    --max_time      Maximum time (time unit)
                    (default: $params.max_time)
    See here for more info: docs/usage.md
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
summary['container']                                   = params.container

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

ch_pheno = params.input_phenofile ? Channel.value(file(params.input_phenofile)) : Channel.empty()
ch_pheno2 = params.input_phenofile ? Channel.value(file(params.input_phenofile)) : Channel.empty()
ch_pheno3 = params.input_phenofile ? Channel.value(file(params.input_phenofile)) : Channel.empty()

codes_pheno = params.input_id_code_count ? Channel.value(file(params.input_id_code_count)) : Channel.empty()
gwas_input_ch = params.gwas_input ? Channel.value(file(params.gwas_input)) : Channel.empty()

if (params.agg_vcf_file){
    Channel.fromPath(params.agg_vcf_file)
           .ifEmpty { exit 1, "VCF file containing  not found: ${params.agg_vcf_file}" }
           .into {vcf_file; vcfs_to_split; index_to_split}
    vcfs_to_split
        .splitCsv(header: true)
        .map{ row -> [file(row.vcf)] }
        .set { vcfs }
    index_to_split
        .splitCsv(header: true)
        .map{ row -> [file(row.index)] }
        .set { indexes }
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
           .into { vcf_file; vcfs_to_split }
    vcfs_to_split
        .splitCsv(header: true)
        .map{ row -> [file(row.vcf)] }
        .set { vcfs }
}


if (params.fam) {
    Channel.fromPath(params.fam)
    .ifEmpty { exit 1, "FAM file (w/ header) containing phenotype info not found: ${params.fam}" }
    .set { fam }
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


// Get as many processors as machine has
int threads = Runtime.getRuntime().availableProcessors()

/*--------------------------------------------------
  Using different formats: vcf & plink
---------------------------------------------------*/
if (params.agg_vcf_file){
    process merge_agg_vcfs {
        publishDir 'results'
        container 'lifebitai/preprocess_gwas:latest'

        input:
        file(vcfs) from vcfs.collect()
        file(indexes) from indexes.collect()
        file vcf_file from vcf_file
        file(pheno_file) from ch_pheno3
 
        output:
        file 'filtered.vcf' into vcf_plink

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
        #merge vcfs & subset samples
        vcfs=\$(tail -n+2 $vcf_file | awk -F',' '{print \$2}')
        bcftools merge --force-samples \$vcfs > merged.vcf
        sed '1d' $pheno_file | awk -F' ' '{print \$1}' > sample_file.txt
        bcftools view -S sample_file.txt merged.vcf > filtered.vcf
        """
    }
}

if (params.individual_vcf_file) {
    process merge_ind_vcfs {
        publishDir 'results'
        container 'lifebitai/preprocess_gwas:latest'

        input:
        file vcfs from vcfs.collect()
        file vcf_file from vcf_file

        output:
        file 'merged.vcf' into vcf_plink

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
        #merge vcfs
        paste -d, sex.txt $vcf_file > tmp.csv && mv tmp.csv $vcf_file
        vcfs=\$(tail -n+2 $vcf_file | awk -F',' '{print \$3}')
        bcftools merge --force-samples \$vcfs > merged.vcf
        """
    }

}

if (params.agg_vcf_file || params.individual_vcf_file){
    process vcf_2_plink {
        tag "plink"
        publishDir "${params.outdir}/plink", mode: 'copy'
        

        input:
        file vcf from vcf_plink
        file fam from ch_pheno

        output:
        set file('*.bed'), file('*.bim'), file('*.fam') into plink, plink2

        script:
        """
        sed '1d' $fam > tmpfile; mv tmpfile $fam
        # remove contigs eg GL000229.1 to prevent errors
        sed -i '/^GL/ d' $vcf
        plink --keep-allele-order \
        --vcf $vcf \
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
        file bed from bed
        file bim from bim
        file fam from fam

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

if (params.plink_input) {

    process preprocess_plink_folder {

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
}

/*--------------------------------------------------
  Prepare resulting plink files
---------------------------------------------------*/

if (!params.snps) {
    process get_snps {
        tag "plink"
        publishDir 'results', mode: 'copy'

        input:
        set file(bed), file(bim), file(fam) from plink
        file pheno_file from ch_pheno2

        output:
        file("snps.txt") into snps

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
publishDir "${params.outdir}/plink", mode: 'copy'
tag "plink"

input:
set file(bed), file(bim), file(fam) from plink2
file snps from snps

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
    publishDir "${params.outdir}/phewas", mode: 'copy'
    container 'lifebitai/phewas:latest'
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
    phewas.R --pheno_file "${pheno}" --geno_file "${genotypes}" --n_cpus ${task.cpus} --pheno_codes "$params.pheno_codes"
    """
}

process merge_results {
    publishDir "${params.outdir}/merged_results", mode: 'copy'
    container 'lifebitai/phewas:latest'
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
        container 'lifebitai/phewas:latest'
        input:
        file gwas_file from gwas_input_ch
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

