# README for the CNA Framework
Author: Guyuan Tang  
Date: 2023/10/16 - 2024/04/08

## 1. Description 
### 1.1 The project
This project is designed within a master thesis (BINP52). We aim to develop a pipeline to generate copy number profiles and detect copy number signatures from shallow whole genome sequening (sWGS) samples.  

### 1.2 The pipeline
The `Snakemake (v7.32.4)` pipeline will include steps from preprocessing raw reads to obtaining copy number signatures (HGSC CN signatures and pan-cancer CIN signatures). We designed two snakemake workflows, `Snakefile_solution.smk` and `Snakefile_CNsig.smk`, to include step 1-5 and 6-8, respectively. For step 9, we designed scripts `panConusig_pair_local.R` for detecting panConusig in local environment with required pair-sample sheet.  
#### Steps
The version of tools and packages to be used will be specified in each step (see chapter 3). The scripts within the pipeline are based on `Python (v3.11.6)` and `R (v4.3.2)`.
- (1) Preprocessing. This step includes quality assessment and quality trimming on the raw reads. (`Fastp` will be used for QC and trimming, together with `fastqc` and `multiQC` to generate the QC reports.)
- (2) Alignment. The human reference genome will be indexed. And the reads will be mapped to the reference genome. (`BWA` will be used for both indexing and alignment.)
- (3) Clean-up. After alignment, the SAM files will be sorted and the PCR duplicates will be marked and removed. Also, the .sorted.deduplicated.sam will be converted to BAM files. The BAM files will be indexed for later analysis. (`Picard` will be used for sorting SAM, marking duplicates, removing duplicates and converting SAM to BAM. `samtools` will be used for generating the clean_up stats and for indexing the BAM files.)
- (4) Relative copy number profile. The BAM files will be analyzed through fixed-size binning, filtering, correction, normalization to generate the read counts per bin. This data will then used for segmentation of bins and generating the relative copy number profile. (`QDNAseq` will be used for this step.)
- (5) Ploidy and cellularity solutions. The output file from `QDNAseq` contains relative copy number, and we need to estimate ploidy and cellularity in our samples to generate our final absolute copy number profile for comparison. (`Rascal` will be used for this step to find the solutions that best fit our study samples.)
- (6) Absolute copy number profile. We will further use other information (such as TP53 allele frequency) inferring the tumour fraction to select the best ploidy and cellularity solution. We apply this best solution to our relative copy number profile, and generate the final absolute copy number profile for each sample. (`Rascal` will be used for this step.)
- (7) Comparison with the recent HGSC signatures (n=7). The functions should be loaded from the github repository: https://bitbucket.org/britroc/cnsignatures.git .
- (8) Comparison with the Pan-Cancer signatures (n=17). The package `CINSignatureQuantification` will be used to generate the samply-by-component matrix for the Pan-Cancer chromosomal instability signatures.
- (9) Comparison with the panConusig signatures (n=25). Tools including `Battenberg` (`alleleCounter`, `impute2` and `beagle5` were included in this package), `ASCAT.sc` and `panConusig` will be used in this step.  
  
![pipeline](https://github.com/GuyuanTang/BINP52_CNA_Framework/blob/main/pipeline.jpg)

### 1.3 Structure and descriptions
`config` folder contains the `config.yaml` file used for the Snakemake workflow as well as the example input sample sheets.  
- `sample_1.tsv`, `sample_2.tsv` and `sample_pair_3.tsv` is the example input for *Part I*, *Part II* and *Part III*, respectively. All of them should be tab-separated files.  
  
`Data` folder can be a place to store the fastq files.  
`resources` folder contains all the reference files required for running the workflow.  
- `BinAnnotations` folder contains the modified QDNAseq bin annotations.  
- `feat_sig_mat.rds` is the signature-by-component definition matrix for HGSC CN signatures.  
- `PanCan.xlsx` and `PanCan_def.rds` are both the signature-by-component definition matrix for pan-cancer CIN signatures.  
- `Panconusig_id.txt` and `panConusig_def.rds` are both the signature-by-component definition matrix for panConusig signatures.  
- `battenberg` folder should contain the reference files for *Part III* to run. But since the sizes of these files are too large to be saved in github, we provide link to download them instead.  

`workflow` is the main folder containing Snakemake pipelines, environment setting yaml files, and scripts used by the workflow.  
`Other_scripts` contains scripts outside the workflow, including sample sheet generation, solution statistics, and signature analyses. Details as the followings:  
- `get_sample.py` was used to generate the input sample sheet for *Part I Solutions*.  
- `preprocess_stat.R` was used to analyze the performance of the preprocessing steps in *Part I Solutions*.  
- `All_Sample_Solutions.py` was used to extract the solutions (ploidy and cellularity) from *Part I Solutions* for each sample under different bin sizes.  
- `solution_stat.R` was used to analyze the solution outputs from *Part I Solutions*.  
- `get_solutions_samptab.R` was used to generate the input sample sheet containing selected bin sizes for *Part II Signatures*.  
- `select_tab.R` was used to extract the signature-by-component definition matrix for pan-cancer CIN signatures from the reference files they provided. We have stored the definition matrix as `resources/PanCan.xlsx` as well as `resources/PanCan_def.rds`.  


## 2. Operation guide
### 2.1 Workflow *Part I Solutions* 
```
# activate the conda environment installing snakemake
conda activate snakemake
# move into the working directory where you downloaded the github repository
cd <path to the downloaded github repository>
# run the pipeline, for example, 30 threads
snakemake --use-conda --configfile config/config.yaml --cores 30 --snakefile workflow/Snakefile_solution.smk
```
### 2.2 Workflow *Part II Signatures* 
Note: this step works for HGSC CN signatures as well as pan-cancer CIN signatures, but not includes the panConusig signatures.
```
# activate the conda environment installing snakemake
conda activate snakemake
# move into the working directory where you downloaded the github repository
cd <path to the downloaded github repository>
# run the pipeline, for example, 30 threads
snakemake --use-conda --configfile config/config.yaml --cores 30 --snakefile workflow/Snakefile_CNsig.smk
```
### 2.3 Workflow *Part III panConusig* 
Note: for panConusig, it requires "normal-tumor" pair, the sample sheet is different from other signatures (see example `config/sample_pair_3.tsv`). Furthermore, it starts again from the BAM files to generate allele-specific copy number profiles which are completely different from the other two signatures, so it is good to run separately. Besides, it would be better to run in local conda environment instead of within the Snakemake pipeline to aviod dependency conflicts (such as Java) as well as pathway errors of reference files.  
Please see suggestion on setting up the environment in section 3.3.
```
# Firstly, we need to add our local true working directory in the impute reference files
cat resources/battenberg/battenberg_impute_v3/impute_info00.txt | sed 's#<path_to_impute_reference_files>#resources/battenberg/battenberg_impute_v3#g' > resources/battenberg/battenberg_impute_v3/impute_info.txt

# Secondly, we need to build up the environment with required dependencies
Rscript workflow/scripts/panConusig_env.R

# Then, run panConusig generation step by step.
# Please remember to change the directories in the scripts to adapt to the user situation
Rscript workflow/scripts/panConusig_pair_local_1.R
Rscript workflow/scripts/panConusig_pair_local_2.R
Rscript workflow/scripts/panConusig_pair_local_3.R
```

## 3. Workflow Details
### 3.0 Sample tables generation
We include a python script `get_sample.py` to help generate the sample.tsv (see example `config/sample_1.tsv`). The final sample.tsv for each type will include the following columns:  
**sample_name**: the name of the sample (also the library)  
**fastq_1** and **fastq_2**: the paths of the sequencing reads  
The sample.tsv to be used in the workflow should be specified in the `config.yaml`.

### 3.1 Snakefile_solution
The `Snakefile_solution.smk` will finish the jobs in step 1-5, and the final output would be the solutions (ploidy and cellularity) for each sample.
```
rule all:
    input:
        expand(results + '{sample}/05_absolute_CN/solutions/{sample}_' + binsize + 'kb.solution.csv', sample=sample_df.sample_name)
```

#### 3.1.1 Preprocessing
This step includes quality assessment and quality trimming on the raw reads.  
*Tools, Packages and Dependencies*
```
  - fastp=0.23.4
  - multiqc=1.17
  - fastqc=0.12.1
  - typing_extension=4.8.0
```
Before executing the workflow, we specified the final output for the *Preprocessing* step.
```
rule preprocess:
    input:
        results + '01_preprocess/html/' + sample_group + '_multiqc_report.html'
```

Firstly, we used `fastp` to perform quality assessment and quality trimming on the raw reads. Each paired-end sample had a html file to view its sequencing stats before and after filtering and trimming.
```
rule fastp:
    input:
        R1 = lambda wildcards: sample_df.loc[wildcards.sample, 'fastq_1'],
        R2 = lambda wildcards: sample_df.loc[wildcards.sample, 'fastq_2']
    output:
        R1 = results + '01_preprocess/reads/{sample}_R1_preprocess.fastq.gz',
        html = results + '01_preprocess/html/{sample}_fastp.html',
        R2 = results + '01_preprocess/reads/{sample}_R2_preprocess.fastq.gz'
    log: 'log/fastp/{sample}_fastp.log'
    threads: 16
    params: json = results + '01_preprocess/html/{sample}_fastp.json'
    conda: "envs/preprocess_env.yaml"
    shell: """
    fastp --detect_adapter_for_pe \
        --correction --cut_right --thread {threads} \
        --html {output.html} --json {params.json} \
        --in1 {input.R1} --in2 {input.R2} \
        --out1 {output.R1} --out2 {output.R2} \
        2>{log}
    
    rm {params.json}
    """
```
Secondly, we generated the qc reports with `fastqc` to investigate the quality of the preprocessed reads.
```
rule fastqc:
    input:
        R1_seq = results + '01_preprocess/reads/{sample}_R1_preprocess.fastq.gz',
        R2_seq = results + '01_preprocess/reads/{sample}_R2_preprocess.fastq.gz'
    output:
        R1_html = results + '01_preprocess/html/{sample}_R1_preprocess_fastqc.html',
        R1_qc = results + '01_preprocess/reports/{sample}_R1_preprocess_fastqc.zip',
        R2_html = results + '01_preprocess/html/{sample}_R2_preprocess_fastqc.html',
        R2_qc = results + '01_preprocess/reports/{sample}_R2_preprocess_fastqc.zip'
    log: 'log/fastqc/{sample}.fastqc.log'
    params: 
        outdir = results + '01_preprocess/reports/',
        out_html = results + '01_preprocess/html/'
    threads: 2
    conda: 'envs/preprocess_env.yaml'
    shell: """
    fastqc -o {params.outdir} {input.R1_seq} {input.R2_seq} 2>{log}
    mv {params.outdir}*_fastqc.html {params.out_html}
    """
```
Thirdly, we used `multiqc` to combine all the qc reports to generate a clear html showing the preprocessed stats.  
*Note*: `multiqc` running in our workflow required the module `typing_extension`, otherwise it could not be run successfully in our local computer. More details could be found in the `envs/multiqc_env.yaml`.
```
rule multiqc:
    input:
        R1_qc = expand(results + '01_preprocess/reports/{sample}_R1_preprocess_fastqc.zip', sample=sample_df.sample_name),
        R2_qc= expand(results + '01_preprocess/reports/{sample}_R2_preprocess_fastqc.zip', sample=sample_df.sample_name)
    output:
        results + '01_preprocess/html/' + sample_group + '_multiqc_report.html'
    log:
        'log/multiqc.log'
    conda: 'envs/multiqc_env.yaml'
    params: 
        out_name = sample_group + '_multiqc_report.html',
        indir = results + '01_preprocess/reports',
        outdir = results + '01_preprocess/html/',
        group = sample_group
    shell: """
    multiqc -f -n {params.out_name} \
    -o {params.outdir} {params.indir} >{log} 2>{log}
    rm -r {params.outdir}/{params.group}_multiqc_report_data/
    """
```

#### 3.1.2 Alignment
This step includes downloading hg19 reference genome (if not prepared), indexing reference genome, and mapping preprocessed reads to the reference.  
*Tools, Packages and Dependencies*
```
 - bwa=0.7.17
```
Before we executing this step, we specified the final output of *Alignment*.
```
rule alignment:
    input:
        expand(results + '02_alignment/{sample}.unsorted.sam', sample=sample_df.sample_name)
```

Firstly, to prepare the hg19 reference genome, we provided a step to download the genome from UCSC. We downloaded the selected chromosomes (1-22, X) sequencing and combined them together as the hg19.ref.fa. If the genome has already been prepared, then it should be renamed into `hg19.ref.fa.gz`, which enables the automatic detection of the rule.
```
rule download_hg19:
    ### if the hg19 reference genome does not exist, this rule will execute to download and generate the hg19 reference genome
    output:
        genome = 'resources/genome/hg19.ref.fa.gz'
    log:
        'log/genome/download_hg19.log'
    shell: """
    for i in $(seq 1 22) X; do echo $i; wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr${i}.fa.gz -O resources/genome/chr${i}.fa.gz; done 2>{log}

    gunzip resources/genome/*.gz
    
    for a in $(seq 1 22) X; do cat resources/genome/chr${a}.fa >> resources/genome/hg19.ref.fa; done

    gzip resources/genome/hg19.ref.fa

    rm *.fa

    """
```

Secondly, we used `bwa` to index the hg19 reference genome. The option `-a bwtsw` is specifically used on human whole genome   
*Note*: This step will take a long time, because we could not specify the threads. The `ìndex` function in `bwa` would only use 1 core at a time, which will take approximately 9 hours to finish the indexing.
```
rule bwa_index:
    input:
        genome = 'resources/genome/hg19.ref.fa.gz'
    output:
        multiext('resources/genome/hg19', ".amb", ".ann", ".bwt", ".pac", ".sa")
    conda: 'envs/alignment.yaml'
    log: 'log/bwa/bwa_index.log'
    params: outdir = 'resources/genome/'
    threads: 20
    shell: """
    bwa index -p hg19 -a bwtsw {input.genome}
    mv hg19.* {params.outdir}
    """
```

Thirdly, `bwa` was used again to map the preprocessed reads to the indexed genome.
```
rule map_reads:
    input: 
        idx = rules.bwa_index.output,
        link_up = rules.preprocess.input,
        R1 = results + '01_preprocess/reads/{sample}_R1_preprocess.fastq.gz',
        R2 = results + '01_preprocess/reads/{sample}_R2_preprocess.fastq.gz'
    output:
        results + '02_alignment/{sample}.unsorted.sam'
    log: 'log/bwa_mapping/{sample}.log'
    params:
        index_ref = 'resources/genome/hg19'
    conda: 'envs/alignment.yaml'
    threads: config['bwa_mapping']['threads']
    shell: """
    bwa mem -M -t {threads} \
        {params.index_ref} {input.R1} {input.R2} > {output} \
        2>{log}
    """
```

#### 3.1.3 Clean-up
This step includes sorting the sam files, marking PCR duplicates, de-duplication, and indexing the final bam files.
*Tools, Packages and Dependencies*
```
  - picard=2.18.7
  - htslib=1.18
  - samtools=1.18
  - qualimap=2.2.2d
```

Before we executing this step, we specified the final output of *Clean-up*.
```
rule clean_up:
    input: 
        expand(results + '03_clean_up/{sample}/{sample}.sorted.dedup.bai', sample=sample_df.sample_name)
```

To avoid using out of space in the disk, we combined the steps using `picard` together (sorting the SAM files, marking and removing PCR dupicates). And we removed the unsorted and intermediate SAM files.
```
rule sort_sam_dedup: 
    ### using Picard to sort the sam files, to mark and remove the PCR duplicates, and to convert SAM into BAM
    input:
        sam = results + '02_alignment/{sample}.unsorted.sam'
    output:
        results + '03_clean_up/{sample}/{sample}.sorted.dedup.bam'
    log: 
        sort_sam = 'log/sort_sam/{sample}.log',
        de_duplicate = 'log/de_duplicate/{sample}.log'
    conda: 'envs/clean_up.yaml'
    threads: 2
    params: 
        metrix_file = results + '03_clean_up/{sample}/{sample}.metrics.txt',
        sorted_sam = results + '03_clean_up/{sample}/{sample}.sorted.sam'
    shell: """
    picard SortSam \
        INPUT={input.sam} \
        OUTPUT={params.sorted_sam} \
        SORT_ORDER=coordinate \
        2>{log.sort_sam}
    rm {input.sam}

    picard MarkDuplicates \
        INPUT={params.sorted_sam} \
        OUTPUT={output} \
        METRICS_FILE={params.metrix_file} \
        REMOVE_DUPLICATES=true \
        ASSUME_SORT_ORDER=coordinate \
        CLEAR_DT=false \
        2>{log.de_duplicate}
    rm {params.sorted_sam}

    """
```

Secondly, we use `samtools` to show the clean-up stats and to index the bam files. The clean-up stats will be saved to the log file named with the sample names.
```
rule index_bam:
    input:
        results + '{sample}/03_clean_up/{sample}.sorted.dedup.bam'
    output:
        bai = results + '{sample}/03_clean_up/{sample}.sorted.dedup.bai',
        stat = results + '{sample}/03_clean_up/{sample}.bamstat.txt'
    threads: 8
    conda: 'envs/clean_up.yaml'
    shell: """
    samtools flagstat {input} | tee {output.stat}

    samtools index -@ {threads} -o {output.bai} {input}
    """
```
Also, if required, we used `qualimap` to generate the stats of the BAM files: `qualimap bamqc -bam results/03_clean_up/{sample}/{sample}.sorted.dedup.bam --java-mem-size=4G`.

#### 3.1.4 Relative copy number profile
This step includes fix-sized binning, filtering, normalization, smoothening, and segmentation on the reads counts from our samples.  
*Tools, Packages and Dependencies*
```
  - bioconductor-qdnaseq=1.38.0
  - bioconductor-qdnaseq.hg19=1.32.0
```
Before we executing this step, we specified the final output of *relative_CN*.
```
rule relative_CN:
    input:
        rds = expand(results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.rds', sample=sample_df.sample_name),
        tsv = expand(results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.seg.tsv', sample=sample_df.sample_name)
```
The main script for finishing this step can be found as `scripts/runQDNAseq.R`. And the snakemake rule only describes the inputs, outputs, and required parameters.
```
rule QDNAseq:
    input:
        link_up = rules.clean_up.input,
        bamfile = results + '{sample}/03_clean_up/{sample}.sorted.dedup.bam'
    output:
        rds = results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.rds',
        igv = results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.igv',
        seg_tsv = results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.seg.tsv'
    params:
        sample = '{sample}',
        binsize = config['QDNAseq']['binsize'],
        outdir = results + '{sample}/04_relative_CN/' + binsize + 'kb/',
        maxSize = config['QDNAseq']['maxSize']
    threads: config['QDNAseq']['threads']
    conda: 'envs/QDNAseq.yaml'
    script: 'scripts/runQDNAseq.R'
```

#### 3.1.5 Ploidy and cellularity solutions
This step uses Rascal package to calculate the optimal solutions of ploidy and cellularity for each sample. This is also the last step of the `Snakefile_solution.smk`.  
*Tools, Packages and Dependencies*
```
 - r-base
 - r-devtools=2.4.5
 - r-ggplot2=3.4.4
 - r-dplyr=1.1.3
 - r-tidyverse=2.0.0
 - r-readr=2.1.4
 - bioconductor-qdnaseq=1.38.0
 - r-rascal=0.7.0 (this should be installed from github repository)
```
Before we executing this step, we specified the final output of CN_solutions.
```
rule CN_solution: # the rule is the same with rule all at this moment
    input:
        expand(results + '{sample}/05_absolute_CN/solutions/{sample}_' + binsize + 'kb.solution.csv', sample=sample_df.sample_name)
```

Firstly, we need to set the environment and install the rascal package from github. If the environment was settled, we would test it and output a log file.
```
rule rascal_env:
    output:
        "log/rascal_settle_info.txt"
    conda: 'envs/rascal.yaml'
    script: 'scripts/rascal_env.R'
```

Secondly, we used rascal to calculate the optimal solutions (ploidy and cellularity) of the samples.
```
rule rascal_solution:
    input:
        link_up = rules.relative_CN.input,
        env_set = 'log/rascal_settle_info.txt',
        rds = results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.rds'
    output:
        solution = results + '{sample}/05_absolute_CN/solutions/{sample}_' + binsize + 'kb.solution.csv'
    params:
        output_prefix = results + '{sample}/05_absolute_CN/solutions/{sample}_' + binsize + 'kb',
        min_cellularity = config['Rascal']['min_cellularity'],
        script = 'workflow/scripts/fit_CN_solution.R'
    threads: 5
    conda: 'envs/rascal.yaml'
    shell: '''
    Rscript {params.script} -i {input.rds} -o {params.output_prefix} --min-cellularity {params.min_cellularity}
    '''
```
### 3.2 Snakefile_CNsig
After deciding the bin size for each group (30kb for ffTumor, ArchivalVS and MaNiLaVS, 100kb for ffpe, and 50kb for other groups), and selecting the optimal ploidy and cellularity for each sample, the second snakemake `Snakefile_CNsig.smk` was designed for the following analyzation.  

The final outputs for this snakemake workflow would be the sample-by-signature matrix (calculated by cosine similarity based on the sample-by-component matrix) for each type of signature (HGSC CN signatures and pan-cancer CIN signatures). 
```
rule all:
    input:
        CN_sig_SS = results + 'signatures/CN_sig/CN_sig.SSmatrix.rds',
        PanCan_sig_SS = results + 'signatures/PanCan_sig/PanCan_sig.SSmatrix.rds'
```

*Note*: Since the output filenames will differ in the following steps which are difficult to be generated by the `extend` function, we designed functions to generate output file names, which were described in `workflow/scripts/common.py`.

#### 3.2.1 Absolute copy number profile
We then generate absolute copy number profile using `rascal`.  
*Tools, Packages and Dependencies*
```
 - r-base
 - r-devtools=2.4.5
 - r-ggplot2=3.4.4
 - r-dplyr=1.1.3
 - r-tidyverse=2.0.0
 - r-readr=2.1.4
 - bioconductor-qdnaseq=1.38.0
 - r-rascal=0.7.0 (this should be installed from github repository)
```

Before running the scripts, we first defined the final output of this step. 
```
rule absolute_seg_CN:
    input:
        files = get_output_absolute(sample_df, results)
```
Firstly, we need to check again if we set up the required environment for running rascal.
```
rule rascal_env:
    output:
        "log/rascal_settle_info2.txt"
    conda: 'envs/rascal.yaml'
    script: 'scripts/rascal_env.R'
```

Secondly, we run `rascal` to calculate the absolute copy numbers and to apply segmentation on the copy number profiles. In order to make the workflow look tidy, we did not descripe all the output files in the *output* section. This step will also save a tsv file `{sample}_{binsize}kb_CN.tsv` containing the absolute copy numbers (before segmentation).
```
rule rascal_absolute_CN:
    input:
        env_set = 'log/rascal_settle_info2.txt',
        sample_table = config['samp_solutions'],
        rds = lambda wildcards: sample_df.loc[wildcards.sample, 'rds']
    output:
        results + '{sample}/05_absolute_CN/{sample}_{binsize}kb_seg.tsv'
    wildcard_constraints:
        binsize="\d+"
    params:
        sample = '{sample}',
        outdir = results + '{sample}/05_absolute_CN/'
    threads: 5
    conda: 'envs/rascal.yaml'
    script: 'scripts/fit_absoluteCN.R'
```
The output tsv files followed the requirement on format (essential columns) for the following signature validation, including **chromosome**, **start**, **end**, **segVal**.

#### 3.2.2 CN signatures
This step will generate the sample-by-component matrix and the sample-by-signature matrix based on the recent published CN signatures (Macintyre *et al.,* 2018).  
*Tools, Packages and Dependencies*
```
 - r-base
 - r-nmf=0.21.0
 - r-flexmix=2.3_19
 - r-tidyverse=2.0.0
 - bioconductor-yapsa=1.28.0
 - bioconductor-qdnaseq=1.38.0
```
Firstly, we need to set up the required environment and install the analytic scripts from the github repository.
```
rule cn_sig_git:
    output:
        'workflow/scripts/cnsignatures/main_functions.R',
        'workflow/scripts/cnsignatures/helper_functions.R'
    shell: '''
    git clone https://bitbucket.org/britroc/cnsignatures.git 
    mv cnsignatures/ workflow/scripts/
    '''
```

Before running the analyzation, we check whether all the required files are well prepared.
```
rule cn_sig_check:
    input:
        link_up = rules.absolute_seg_CN.input,
        scripts = rules.cn_sig_git.output
    output:
        'log/cn_sig_settle.txt'
    shell: """
    echo 'Finished preparation for signature validation.' > {output}
    """
```

Secondly, we run the signature validation programme on our samples. We also save a RDS file of the sample-by-component matrix, txt files and simple heatmaps of both matrices for further analysis.
```
rule CN_signature:
    input:
        # check whether the essential git repository have been downloaded
        main = 'workflow/scripts/cnsignatures/main_functions.R',
        helper = 'workflow/scripts/cnsignatures/helper_functions.R',
        # check all the absolute copy number profiles have successfully been generated
        check_data = rules.cn_sig_check.output,
        # the input segment files
        seg_files = rules.absolute_seg_CN.input
    output:
        CN_sig = results + 'signatures/CN_sig/CN_sig.SSmatrix.rds'
    params:
        sample_info = config['samp_solutions'],
        def_SC = 'workflow/scripts/cnsignatures/data/feat_sig_mat.rds',
        indir = results,
        outdir = results + 'signatures/CN_sig/'
    threads: 20
    conda: 'envs/CNsig.yaml'
    script: 'scripts/CN_sig.R'
```

#### 3.2.3 Pan-Cancer signatures
This step will generate the sample-by-component matrix and the sample-by-signature matrix based on the recent published pan-cancer CIN signatures (Drews *et al.,* 2022).  
*Tools, Packages and Dependencies*
```
 - r-base
 - r-devtools=2.4.5
 - r-tidyverse=2.0.0
 - r-lsa=0.23.3
 - r-CINSignatureQuantification=1.2.0 (this should be installed from the github repository)
```

Firstly, we need to check if we have downloaded the definition signature-by-component table for PanCan signature. (This should be prepared in the `resources/` folder automatically.)
```
rule PanCan_def_tab:
    # user should download the supplements manually
    # website: https://www-nature-com.ludwig.lub.lu.se/articles/s41586-022-04789-9#Sec17
    input:
        config['supple_tab']
    output:
        'resources/PanCan.xlsx'
    params:
        select_tab = 'workflow/scripts/select_tab.R',
        folder = config['supple_folder']
    conda: 'envs/PanCanSig.yaml'
    shell:'''
        mv '{input}' resources/Supplementary_Tables_15-22.xlsx
        Rscript {select_tab} -i resources/Supplementary_Tables_15-22.xlsx -o resources/PanCan.xlsx
        rm -r {folder}
    '''
``` 

Secondly, we need to set up the environment and install the required scripts from github repository.
```
rule PanCan_sig_env:
    input:
        'resources/PanCan.xlsx'
    output:
        "log/PanCan_settle_info.txt"
    conda: 'envs/PanCanSig.yaml'
    script: 'scripts/PanCan_sig_env.R'
```

Thirdly, we run the signature validation programme on our samples. The actural outputs also include a RDS file containing all the matrix information (activities and weights) and two simple heatmaps (one for activities and one for the sample-by-component matrix) for each group.
```
rule PanCan_sig:
    input:
        env_set = "log/PanCan_settle_info.txt",
        seg_files = rules.absolute_seg_CN.input
    output:
        PanCan_sig_SS = results + 'signatures/PanCan_sig/PanCan_sig.SSmatrix.rds'
    params:
        sample_info = config['samp_solutions'],
        def_SC = 'resources/PanCan.xlsx',
        indir = results,
        outdir = results + 'signatures/PanCan_sig/'
    threads: 10
    conda: 'envs/PanCanSig.yaml'
    script: 'scripts/PanCan_sig.R'
```

### 3.3 panConusig signatures
This step will generate the sample-by-component matrix and the sample-by-signature matrix based on the recent published panConusig signatures (Steels *et al.,* 2022).  
*Tools, Packages and Dependencies*  
(Suggestions on installing the dependencies are also described below.)
```
# dependencies below are suggested to be installed in conda environment
 - r-base=4.3.3
 - r-tidyverse=2.0.0
 - r-biocmanager=1.30.22
 - cancerit-allelecount=4.3.0
 - r-devtools=2.4.5
 - r-gridextra=2.3
 - r-rcolorbrewer=1.1_3
 - r-splines2=0.5.1
 - r-readr=2.1.5
 - r-doparallel=1.0.17
 - r-ggplot2=3.4.4
 - r-gtools=3.9.5
 - parallel=20170422
 - bioconductor-variantannotation=1.48.1
 - r-xgboost=2.0.3
 - bioconductor-genomicranges=1.54.1
 - bioconductor-biostrings=2.70.1
 - bioconductor-dnacopy=1.76.0
 - curl=8.5.0
 - impute2=2.3.2
 - java-jdk=8.0.92
 - r-rjava=1.0_10
 - openjdk=21.0.2
 - r-lsa=0.73.3

 # dependencies below are suggested to be installed in R
 - r-minfi=1.48.0
 - r-conumee=1.36.0
 - r-Rsamtools=2.18.0
 - r-copynumber=1.29.0.9000
 - r-ASCAT=3.1.2
 - r-ASCAT.sc=0.1
 - r-battenberg=2.2.10
 - r-panConusig=0.1.0
```
Firstly, we apply preprocessing on BAM files (generated in *Part I Solutions*) to derive allele frequency files and phased haplotype data by using modified `Battenberg`. Instead of running all the remaining steps, the package was adapted to only producing the files we needed, which were the first two steps in the original function (see `workflow/scripts/usr_battenberg_pair.R`). Besides, the filter on read depth is also changed to 2. The default value was 10, which was too large for our sWGS samples.  
This step will take approximately an hour for each sample.  
Below is example for running on one sample: (see `workflow/scripts/panConusig_pair_local_1.R` for details)
```
# set the working directory
working_dir <- paste0("/home/researcher/TangGY/BINP52/Workflow/Draft/results/",sampleID,"/06_panConusig/")
dir.create(working_dir)
setwd(working_dir)
# load the parameters
TUMOURNAME = sampleID
NORMALNAME = sample_df[which(sample_df$Sample==sampleID),'Sample_ref']
results_output = "/home/researcher/TangGY/BINP52/Workflow/Draft/results/"
NORMALBAM = paste0(results_output,NORMALNAME,'/03_clean_up/',NORMALNAME,'.sorted.dedup.bam')
TUMOURBAM = paste0(results_output,TUMOURNAME,'/03_clean_up/',TUMOURNAME,'.sorted.dedup.bam')
IS.MALE = FALSE
SKIP_ALLELECOUNTING = FALSE
SKIP_PREPROCESSING = FALSE
SKIP_PHASING = FALSE
NTHREADS = 30
PRIOR_BREAKPOINTS_FILE = NULL
  
analysis = "paired"
  
JAVAJRE = "java"
ALLELECOUNTER = "alleleCounter"
IMPUTE_EXE = "impute2"
  
GENOMEBUILD = "hg19"
USEBEAGLE = T
  
# General static
# location of the reference files
if (GENOMEBUILD=="hg19") {
    impute_basedir = "/home/researcher/TangGY/BINP52/Workflow/Draft/resources/battenberg"
    IMPUTEINFOFILE = file.path(impute_basedir, "battenberg_impute_v3/impute_info.txt")
    G1000ALLELESPREFIX = file.path(impute_basedir, "battenberg_1000genomesloci2012_v3/1000genomesAlleles2012_chr")
    G1000LOCIPREFIX = file.path(impute_basedir, "battenberg_1000genomesloci2012_v3/1000genomesloci2012_chr")
    GCCORRECTPREFIX = file.path(impute_basedir, "battenberg_wgs_gc_correction_1000g_v3/1000_genomes_GC_corr_chr_")
    REPLICCORRECTPREFIX = file.path(impute_basedir, "battenberg_wgs_replic_correction_1000g_v3/1000_genomes_replication_timing_chr_")
    
    # WGS specific static
    PROBLEMLOCI = file.path(impute_basedir, "probloci_270415.txt.gz")
    GENOME_VERSION = "b37"
    GENOMEBUILD = "hg19"
    BEAGLE_BASEDIR = "/home/researcher/TangGY/BINP52/Workflow/Draft/resources/battenberg/beagle"
    BEAGLEJAR = file.path(BEAGLE_BASEDIR, "beagle.22Jul22.46e.jar")
    BEAGLEREF.template = file.path(BEAGLE_BASEDIR, GENOME_VERSION, "chrCHROMNAME.1kg.phase3.v5a.b37.bref3")
    BEAGLEPLINK.template = file.path(BEAGLE_BASEDIR, GENOME_VERSION, "plink.chrCHROMNAME.GRCh37.map")
    
    CHROM_COORD_FILE = file.path(impute_basedir, "gcCorrect_chromosome_coordinates_hg19.txt")
    
    
    PLATFORM_GAMMA = 1
    PHASING_GAMMA = 1
    SEGMENTATION_GAMMA = 10
    SEGMENTATIIN_KMIN = 3
    PHASING_KMIN = 1
    CLONALITY_DIST_METRIC = 0
    ASCAT_DIST_METRIC = 1
    MIN_PLOIDY = 1.6
    MAX_PLOIDY = 4.8
    MIN_RHO = 0.1
    MIN_GOODNESS_OF_FIT = 0.63
    BALANCED_THRESHOLD = 0.51
    MIN_NORMAL_DEPTH = 2 #default setting is 10 which is not suitable for sWGS, and for most positions, the good read depth is only 1
    MIN_BASE_QUAL = 20
    MIN_MAP_QUAL = 35
    CALC_SEG_BAF_OPTION = 1
  }
# run battenberg
battenberg_new(analysis=analysis,
                tumourname=TUMOURNAME, 
                normalname=NORMALNAME, 
                tumour_data_file=TUMOURBAM, 
                normal_data_file=NORMALBAM, 
                ismale=IS.MALE, 
                imputeinfofile=IMPUTEINFOFILE, 
                g1000prefix=G1000LOCIPREFIX, 
                g1000allelesprefix=G1000ALLELESPREFIX, 
                gccorrectprefix=GCCORRECTPREFIX, 
                repliccorrectprefix=REPLICCORRECTPREFIX, 
                problemloci=PROBLEMLOCI, 
                data_type="wgs",
                impute_exe=IMPUTE_EXE,
                allelecounter_exe=ALLELECOUNTER,
                usebeagle=USEBEAGLE,
                beaglejar=BEAGLEJAR,
                beagleref=BEAGLEREF.template,
                beagleplink=BEAGLEPLINK.template,
                beaglemaxmem=10,
                beaglenthreads=30,
                beaglewindow=40,
                beagleoverlap=4,
                javajre=JAVAJRE,
                nthreads=NTHREADS,
                platform_gamma=PLATFORM_GAMMA,
                phasing_gamma=PHASING_GAMMA,
                segmentation_gamma=SEGMENTATION_GAMMA,
                segmentation_kmin=SEGMENTATIIN_KMIN,
                phasing_kmin=PHASING_KMIN,
                clonality_dist_metric=CLONALITY_DIST_METRIC,
                ascat_dist_metric=ASCAT_DIST_METRIC,
                min_ploidy=MIN_PLOIDY,
                max_ploidy=MAX_PLOIDY,
                min_rho=MIN_RHO,
                min_goodness=MIN_GOODNESS_OF_FIT,
                uninformative_BAF_threshold=BALANCED_THRESHOLD,
                min_normal_depth=MIN_NORMAL_DEPTH,
                min_base_qual=MIN_BASE_QUAL,
                min_map_qual=MIN_MAP_QUAL,
                calc_seg_baf_option=CALC_SEG_BAF_OPTION,
                skip_allele_counting=SKIP_ALLELECOUNTING,
                skip_preprocessing=SKIP_PREPROCESSING,
                skip_phasing=SKIP_PHASING,
                prior_breakpoints_file=PRIOR_BREAKPOINTS_FILE,
                GENOMEBUILD=GENOMEBUILD,
                chrom_coord_file=CHROM_COORD_FILE)
```

Secondly, we run `ASCAT.sc` as suggested on the allele frequency files and phased haplotype data to calculate the segmented absolute allele-specific copy number profiles. The final output profiles will be adjusted to have the required column names for next step.  
This step will take approximately 2 min for each sample.  
Below is an example for running on one sample: (see `workflow/scripts/panConusig_pair_local_2.R` for details)
```
# set the working directory
working_dir <- paste0("/home/researcher/TangGY/BINP52/Workflow/Draft/results/",sampleID,"/06_panConusig/")
setwd(working_dir)

# run ASCAT.sc
ASCAT_out <- paste0("/home/researcher/TangGY/BINP52/Workflow/Draft/results/",sampleID,"/06_panConusig/ASCAT_out/")
dir.create(ASCAT_out)
results_output = paste0("/home/researcher/TangGY/BINP52/Workflow/Draft/results/",sampleID,'/06_panConusig/')
## path to required files
path_to_phases=lapply(sampleID, function(SAMPLENAME)
  {
    paste0(results_output, SAMPLENAME,
           "_beagle5_output_chr",c(1:22),".txt.vcf.gz")
  })
  
list_ac_counts_paths=lapply(sampleID, function(SAMPLENAME)
  {
    paste0(results_output, SAMPLENAME,
           "_alleleFrequencies_chr",c(1:22),".txt")
  })
bins <- sample_df[which(sample_df$Sample==sampleID),'Binsize']
bins <- as.numeric(bins) * 1000
TUMOURBAM = paste0("/home/researcher/TangGY/BINP52/Workflow/Draft/results/",sampleID,'/03_clean_up/',sampleID,'.sorted.dedup.bam')
## run the main programme
res <- run_sc_sequencing(tumour_bams=TUMOURBAM,
                        allchr=paste0("chr",c(1:22)),
                        sex='female',
                        binsize=bins,
                        chrstring_bam="chr",
                        purs = seq(0.01, 1, 0.001), #if start from 0, will affect the following operations
                        ploidies = seq(1.7, 5, 0.01),
                        maxtumourpsi=5,
                        build="hg19",
                        MC.CORES=30,
                        outdir=ASCAT_out,
                        projectname="ASCAT_CN",
                        segmentation_alpha=0.01,
                        path_to_phases=path_to_phases,
                        list_ac_counts_paths=list_ac_counts_paths,
                        predict_refit=TRUE,
                        multipcf=FALSE)
  
## alter the format of the final output
in_file <- paste0(ASCAT_out, 'as_cna_profile_', sampleID, '.sorted.dedup.bam_bam1.txt')
df <- read.csv(in_file, sep = '\t')
df$sample <- sampleID
df <- df %>% select('sample', 'chr', 'startpos', 'endpos', 'nA', 'nB')
colnames(df) <- c('sample', 'chr', 'startpos', 'endpos', 'nMajor', 'nMinor')
# print the file
write.table(df, file = paste0(ASCAT_out, sampleID, '_as_cna_profile.tsv'), quote = FALSE, sep = '\t', row.names = FALSE)
```

Finally, we used the outputs from the above step to extract the sample-by-component matrix as well as the samply-by-signature cosine similarity matrix (see `workflow/scripts/panConusig_pair_local_3.R` for details).

## Acknowledgement
Thanks a lot to all the support and advice from [Ingrid Lab](https://github.com/IngridHLab)!

## Other information
If you have any questions towards the workflow, such as handling the reference files, please feel free to contact the author through email: gu5747ta-s@student.lu.se.  
Any advice and questions are welcome. Thank you so much!