# README for the CNA Framework
Author: Guyuan Tang  
Date: 2023/10/16 - 

## 1. Description 
### 1.1 The project
This project is designed within a master thesis (BINP52). We aim to develop a pipeline to generate copy number profiles from shallow whole genome sequening (sWGS) samples.  

### 1.2 The pipeline
The `Snakemake (v7.32.4)` pipeline will include steps from preprocessing raw reads with QC to obtaining absolute copy number profiles. We designed two snakemake workflows, `Snakefile_solution.smk` and `Snakefile_CNsig.smk`, to include step 1-5 and 6-8, respectively.  
#### Steps
The version of tools and packages to be used will be specified in each step (see chapter 2). The scripts within the pipeline are based on `Python (v3.11.6)` and `R (v4.3.2)`.
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

### 1.3 File structure and descriptions


## 2. Workflow Details
### 2.0 Sample tables generation
We include a python script `get_sample.py` to help generate the sample.tsv for each type of samples from the provided xlsx file describing the samples. The final sample.tsv for each type will include the following columns:  
**sample_name**: the name of the sample (also the library)  
**patient**: the patient ID  
**type**: the sample type  
**fastq_1** and **fastq_2**: the paths of the sequencing reads  
The sample.tsv to be used in the workflow should be specified in the `config.yaml`.

### 2.1 Snakefile_solution
The `Snakefile_solution.smk` will finish the jobs in step 1-5, and the final output would be the solutions (ploidy and cellularity) for each sample.
```
rule all:
    input:
        expand(results + '{sample}/05_absolute_CN/solutions/{sample}_' + binsize + 'kb.solution.csv', sample=sample_df.sample_name)
```

### 2.1.1 Preprocessing
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

### 2.1.2 Alignment
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

### 2.1.3 Clean-up
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

### 2.1.4 Relative copy number profile
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

### 2.1.5 Ploidy and cellularity solutions
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
### 2.2 Snakefile_CNsig
After deciding the bin size for each group (30kb for ffTumor, ArchivalVS and MaNiLaVS, 100kb for ffpe, and 50kb for other groups), and selecting the optimal ploidy and cellularity for each sample, the second snakemake `Snakefile_CNsig.smk` was designed for the following analyzation.  

The final outputs for this snakemake workflow would be the sample-by-signature matrix (calculated by cosine similarity based on the sample-by-component matrix) for each type of signature. 
```
rule all:
    input:
        CN_sig_SS = results + 'signatures/CN_sig/CN_sig.SSmatrix.rds',
        PanCan_sig_SS = results + 'signatures/PanCan_sig/PanCan_sig.SSmatrix.rds',
        panConusig_SS = results + 'signatures/panConusig/panConusig.SSmatrix.rds'
```

*Note*: Since the output filenames will differ in the following steps which are difficult to be generated by the `extend` function, we designed functions to generate output file names, which were described in `workflow/scripts/common.py`.

### 2.2.1 Absolute copy number profile
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

### 2.2.2 CN signatures
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

### 2.2.3 Pan-Cancer signatures
This step will generate the sample-by-component matrix and the sample-by-signature matrix based on the recent published CN signatures (Drews *et al.,* 2022).  
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

### 2.2.4 panConusig signatures
This step will generate the sample-by-component matrix and the sample-by-signature matrix (calculated by cosine similarity) of the panConusig signatures. (Steel *et al.,* 2022)  
*Tools, Packages and Dependencies*
```
 - r-base
 - r-devtools=2.4.5
 - r-tidyverse=2.0.0
 - r-lsa=0.73.3
 - r-biocmanager=1.30.22
 - cancerit-allelecount=4.3.0
 - impute2=2.3.2
 - r-gridextra=2.3
 - r-rcolorbrewer=1.1_3
 - r-splines2=0.5.1
 - r-readr=2.1.5
 - r-doParallel=1.0.17
 - r-ggplot2=3.4.4
 - r-gtools=3.9.5
 - parallel=20240122
 - java-jdk=8.0.92
 - r-rjava=1.0_10
 - bioconductor-variantannotation=1.48.1
 - r-xgboost=2.0.3
 - curl=8.5.0
 - bioconductor-genomicranges=1.54.1
 - bioconductor-biostrings=2.70.1
 - bioconductor-dnacopy=1.76.0
 - openjdk=21.0.2
```
Firstly, we need to prepare the environment for `Battenberg`, which should be installed in R.
```
rule panConusig_env:
    output:
        "log/panConusig_settle_info.txt"
    conda: 'envs/panConusig.yaml'
    script: 'scripts/panConusig_env.R'
```
Secondly, we need to prepare the reference files for `Battenberg` and `beagle5`.
```
rule panConusig_ref_prep:
    input:
        env_set = "log/panConusig_settle_info.txt",
        impute_00 = working_dir + 'resources/battenberg/impute_info00.txt'
    output:
        ref_battenberg = get_ref_battenberg(working_dir),
        ref_beagle = get_ref_beagle(working_dir)
    params:
        workdir = working_dir,
        impute_info = working_dir + 'resources/battenberg/impute_info.txt'
    shell: '''
    cat {input.impute_00} | sed 's#<path_to_impute_reference_files>#{params.workdir}#g' > {params.impute_info}
    '''
```
Thirdly, we run the `Battenberg` and `ASCAT.sc` to generate the allele frequency files and phased files which are required for the next step.
```
rule panConusig_preprocessing:
    input: 
        ref_files = rules.panConusig_ref_prep.output,
        usr_battenberg = "workflow/scripts/usr_battenberg.R",
        bam_file = results + '{sample}/03_clean_up/{sample}.sorted.dedup.bam'
    output:
        preprocess_out = get_panConusig_preprocess(results, wildcards.sample)
    params:
        sampleID = '{sample}',
        impute_ref_dir = working_dir + 'resources/battenberg',
        beagle_ref_dir = working_dir + 'resources/battenberg/beagle',
        result_dir = results + '{sample}/06_panConusig/',
        ASCAT_outdir = results + '{sample}/06_panConusig/ASCAT_out/',
        binsize = lambda wildcards: sample_df.loc[wildcards.sample, 'Binsize']
    threads: 30
    conda: 'envs/panConusig.yaml'
    script: 'scripts/panConusig_preprocess.R'
```
Finally, similar to CN_sig and PanCan_sig, we extract the sample-by-component matrix of panConusig and apply cosine similarity method to discover the signatures that share most components.
```
rule panConusig_sig:
    input:
        cna_profiles = get_cna_profile(sample_df, results)
    output:
        panConusig_SS = results + 'signatures/panConusig/panConusig.SSmatrix.rds'
    params:
        def_SC = 'resources/panConusig_id.txt',
        outdir = results + 'signatures/panConusig/'
    conda: 'envs/panConusig.yaml'
    script: 'scripts/panConusig.R'
```
