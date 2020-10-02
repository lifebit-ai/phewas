# CB integration with GWAS pipeline

CB output needs to be prepared and fed into pheWAS pipeline so users have an end-to-end experience.

## 1. Requirements

According to https://lifebit.gitbook.io/cloudos/pipelines-documentations-and-examples-1/nextflow-pipelines/phewas and `main.nf` :

### 1.1 - Inputs

- **--vcf_file:** CSV file containing path to VCF files and the following columns:

    1. The `sampleid` column contains a unique sample ID for each individual. These IDs are the same as the IDs used in the file for the `pheno_file` parameter
    2. The `vcf` column contains the path to an input VCFs. Each VCF file contains one individual.
    3. The `PHENO` column contains a phenotype of interest such as SMOKER. This can be used for the `pheno` parameter to select to top SNPs of interest if the snps parameter is not provided. The phenotype values are `1` = control, `2` = case and `-0`/`0`/non-numeric = missing data if case/control
    4. The `SEX` column contains the values `1` = male,  `2` = female, `0` = unknown
    5. Other `PHENO` columns can also be used such as `AGE`

```
sampleid,vcf,SMOKER,SEX,AGE
1,s3://test-vcfs/mini_cohort/1_1.vcf.gz,1,2,39
2,s3://test-vcfs/mini_cohort/2_2.vcf.gz,1,1,52
3,s3://test-vcfs/mini_cohort/3_3.vcf.gz,2,1,42
4,s3://test-vcfs/mini_cohort/4_4.vcf.gz,2,2,44
5,s3://test-vcfs/mini_cohort/5_5.vcf.gz,2,2,23
```

OR

- **.bed .bim .fam files**

PLUS:

- **--pheno_file** CSV file containing phenotypes & the following columns:
    1. The `id` column contains a unique sample ID for each individual. These IDs are the same as the IDs used in the file for the `vcf_file` parameter.
    2. The second column can have the value `icd9` or `doid` depending on the disease classification codes used.

```
id,doid,count
1170,3055,4
1170,12750,4
1170,12385,4
1170,1498,4
1170,12385,9
1170,3055,10
1170,12385,5
1170,0060859,3
1170,12385,3
1170,14239,3
```

EITHER:
- **--pheno** Comma-separated list of phenotypes to be tested (i.e SMOKER)
- **--snps** File containing relevant SNPs to be tested. One SNP per row.
```
rs1002365
rs1023115
rs1041676
rs1048418
rs10489265
rs10493582
rs10493644
rs10494452
rs10494464
```

**Optional**

- **--pheno_codes** pheno codes used for `pheno_file`. Options are `icd9` or `doid` (default = doid)
- **--snp_threshold** used to extract the top SNPs (default = 0.05)

### 1.2 - Outputs
`phewas` folder containing:
- `phewas_man.png` Manhattan plot

`plink` folder containing:
- `r_genotypes.raw` contains information from the FAM file with the number of alleles for each given SNP of interest for all individuals (can be loaded into R)
- `r_genotypes.log` log file detailing the execution of the PLINK command

`vcf2plink` folder containing PLINK binaries:
- plink.bed binary version of ped or pedigree file
- plink.bim map file containing variant information
- plink.fam FAM files contain phenotypic information and how the individuals relate to one another

`Visualisations` directory containing:
- `.report.json` file used to render the visualisations within CloudOS

### 1.3 - Scripts


# Future work

