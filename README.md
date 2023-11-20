# README for the CNA Framework
Author: Guyuan Tang  
Date: 2023/10/16 - 

## 1. Description 
### 1.1 The project
This project is designed within a master thesis (BINP52). We aim to develop a pipeline to generate copy number profiles from shallow whole genome sequening (sWGS) samples.  

### 1.2 The pipeline
The `Snakemake (v7.32.4)` pipeline will include steps from preprocessing raw reads with QC to obtaining absolute copy number profiles.  
#### Steps
The version of tools and packages to be used will be specified in each step (see chapter 2). The scripts within the pipeline are based on `Python (v3.11.6)` and `R (v4.3.1)`.
- (1) Preprocessing. This step includes quality assessment and quality trimming on the raw reads. (`Fastp` will be used for QC and trimming, together with `fastqc` and `multiQC` to generate the QC reports.)
- (2) Alignment. The human reference genome will be indexed. And the reads will be mapped to the reference genome. (`BWA` will be used for both indexing and alignment.)
- (3) Clean-up. After alignment, the SAM files will be sorted and the PCR duplicates will be marked and removed. Also, the .sorted.deduplicated.sam will be converted to BAM files. The BAM files will be indexed for later analysis. (`Picard` will be used for sorting SAM, marking duplicates, removing duplicates and converting SAM to BAM. `samtools` will be used for generating the clean_up stats and for indexing the BAM files.)
- (4) Relative copy number profile. The BAM files will be analyzed through fixed-size binning, filtering, correction, normalization to generate the read counts per bin. This data will then used for segmentation of bins and generating the relative copy number profile. (`QDNAseq` will be used for this step.)
- (5) Ploidy and cellularity solutions. The output file from `QDNAseq` contains relative copy number, and we need to estimate ploidy and cellularity in our samples to generate our final absolute copy number profile for comparison. (`Rascal` will be used for this step to find the solutions that best fit our study samples.)
- (6) Absolute copy number profile. We will further use other information (such as TP53 allele frequency) inferring the tumour fraction to select the best ploidy and cellularity solution. We apply this best solution to our relative copy number profile, and generate the final absolute copy number profile for each sample. (`Rascal` will be used for this step.)
- (7) Comparison with the pan-cancer signatures.
- (8) Comparison with the recent HGSC signatures.

## 2. Workflow Details
### 2.0 Sample tables generation
We include a python script `get_sample.py` to help generate the sample.tsv for each type of samples from the provided xlsx file describing the samples. The final sample.tsv for each type will include the following columns:  
**sample_name**: the name of the sample (also the library)  
**patient**: the patient ID  
**fastq_1** and **fastq_2**: the paths of the sequencing reads  
The sample.tsv to be used in the workflow should be specified in the `config.yaml`.

### 2.1 Preprocessing
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
    threads: 10
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
    threads: 10
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

### 2.2 Alignment
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

Secondly, we used `bwa` to index the hg19 reference genome. The option `-a bwtsw`   
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
    threads: 10
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

### 2.3 Clean-up
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

Firstly, we used `picard` to sort the sam files with coordinate mode. And we removed the unsorted sam files to save space in the operational computer.
```
rule sort_sam: 
    input:
        sam = results + '02_alignment/{sample}.unsorted.sam',
        link_up = rules.alignment.input
    output:
        results + '03_clean_up/{sample}/{sample}.sorted.sam'
    log: 'log/sort_sam/{sample}.log'
    conda: 'envs/clean_up.yaml'
    shell: """
    picard SortSam \
        INPUT={input.sam} \
        OUTPUT={output} \
        SORT_ORDER=coordinate \
        2>{log}

    """
```

Secondly, we marked and removed the PCR duplicates detected by `picard`. Because we did not add the DT tag on our sorted sam files, so we decided to set `CLEAN_DT=false`.
```
rule de_duplicate:
    input: 
        results + '03_clean_up/{sample}/{sample}.sorted.sam'
    output:
        results + '03_clean_up/{sample}/{sample}.sorted.dedup.bam'
    log: 'log/de_duplicate/{sample}.log'
    threads:10
    params:
        metrix_file = results + '03_clean_up/{sample}/{sample}.metrics.txt'
    conda: 'envs/clean_up.yaml'
    shell: """
    picard MarkDuplicates \
        INPUT={input} \
        OUTPUT={output} \
        METRICS_FILE={params.metrix_file} \
        REMOVE_DUPLICATES=true \
        ASSUME_SORT_ORDER=coordinate \
        CLEAR_DT=false \
        2>{log}
    
    rm {input}
    """
```

Thirdly, we use `samtools` to show the clean-up stats and to index the bam files. The clean-up stats will be saved to the log file named with the sample names.
```
rule index_bam:
    input:
        results + '03_clean_up/{sample}/{sample}.sorted.dedup.bam'
    output:
        results + '03_clean_up/{sample}/{sample}.sorted.dedup.bai'
    log: 'log/bam_stat/{sample}.log'
    threads: 10
    conda: 'envs/clean_up.yaml'
    shell: """
    samtools flagstat {input} | tee {log}

    samtools index -@ {threads} -o {output} {input}
    """
```
Also, if required, we used `qualimap` to generate the stats of the BAM files: `qualimap bamqc -bam results/03_clean_up/{sample}/{sample}.sorted.dedup.bam --java-mem-size=4G`.

### 2.4 Relative copy number profile
This step includes fix-sized binning, filtering, normalization, smoothening, and segmentation on the reads counts from our samples.  
*Tools, Packages and Dependencies*
```
  - bioconductor-qdnaseq=1.36.0
  - bioconductor-qdnaseq.hg19=1.30.0
```
Before we executing this step, we specified the final output of *relative_CN*.
```
rule relative_CN:
    input:
        rds = expand(results + '04_relative_CN/{sample}/{sample}.rds', sample=sample_df.sample_name),
        tsv = expand(results + '04_relative_CN/{sample}/{sample}.seg.tsv', sample=sample_df.sample_name)
```
The main script for finishing this step can be found as `scripts/runQDNAseq.R`. And the snakemake rule only describes the inputs, outputs, and required parameters.
```
rule QDNAseq:
    input:
        link_up = rules.clean_up.input,
        bamfile = results + '03_clean_up/{sample}/{sample}.sorted.dedup.bam'
    output:
        rds = results + '04_relative_CN/{sample}/{sample}.rds',
        igv = results + '04_relative_CN/{sample}/{sample}.igv',
        seg_tsv = results + '04_relative_CN/{sample}/{sample}.seg.tsv'
    params:
        sample = '{sample}',
        binsize = config['QDNAseq']['binsize'],
        outdir = results + '04_relative_CN/{sample}/'
    threads: 10
    conda: 'envs/QDNAseq.yaml'
    script: 'scripts/runQDNAseq.R'
```

### 2.5 Ploidy and cellularity solutions
This step uses Rascal package to calculate the optimal solutions of ploidy and cellularity for each sample. This is also the last step of the `Snakefile_solution.smk`.  
*Tools, Packages and Dependencies*
```
 - r-base
 - r-devtools=2.4.5
 - r-ggplot2=3.4.4
 - r-dplyr=1.1.3
 - r-tidyverse=2.0.0
 - r-readr=2.1.4
```
Before we executing this step, we specified the final output of CN_solutions.
```
rule CN_solution: # the rule is the same with rule all at this moment
    input:
        expand(results + '05_absolute_CN/{sample}/{sample}.solution.csv', sample=sample_df.sample_name)
```

Firstly, we need to set the environment and download the rascal package from github. If the environment was settled, we would test it and output a log file.
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
        rds = results + '04_relative_CN/{sample}/{sample}.rds'
    output:
        solution = results + '05_absolute_CN/{sample}/{sample}.solution.csv'
    params:
        output_prefix = results + '05_absolute_CN/{sample}/{sample}',
        min_cellularity = config['Rascal']['min_cellularity'],
        script = 'workflow/scripts/fit_CN_solution.R'
    threads: 10
    conda: 'envs/rascal.yaml'
    shell: '''
    Rscript {params.script} -i {input.rds} -o {params.output_prefix} --min-cellularity {params.min_cellularity}
    '''
```
