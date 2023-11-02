# Title: Snakefile_draft.smk
# Author: Guyuan TANG
# Date: 2023/10/16 - 

import pandas as pd

# specify the configuration file
configfile: "config/config.yaml"

# specify the samples
sample_df = (pd.read_csv(config['samples'], 
    sep='\t', 
    dtype={'sample_name':str, 'type_group':str, 'fastq_1':str, 'fastq_2':str})
    .set_index('sample_name', drop=False))


# specify the final output of the whole workflow
"""
the final output for the whole workflow should be the absolute copy number profile for each type/group of samples.
For example: archive (pre-diagnosis), diagnosis, tumor tissue.
"""
rule all:
    input:
        expand('results/05_absolute_CN/{sample}.absoluteCN.csv', sample=sample_df.sample_name)


########## 1 Preprocessing ####################
"""
The final output for the step preprocessing should be a multiqc report generating all the preprocessed QC statistics.
"""
rule preprocess:
    input:
        'results/01_preprocess/report/samples_report.html'

# 1.1 quality trimming on reads
rule fastp:
    ### we will use fastp here for trimming on adapter and quality.
    input:
        R1 = lambda wildcards: sample_df.loc[wildcards.sample, 'fastq_1']
        R2 = lambda wildcards: sample_df.loc[wildcards.sample, 'fastq_2']
    output: # do we need to store unpaired reads also?
        R1 = 'results/01_preprocess/reads/{sample}_{type_group}_R1.trimmed.fastq.gz'
        R1_html = 'results/01_preprocess/html/{sample}_R1.fastp.html'
        R1_json = 'results/01_preprocess/json/{sample}_R1.fastp.json'
        R2 = 'results/01_preprocess/reads/{sample}_{type_group}_R2.trimmed.fastq.gz'
        R2_html = 'results/01_preprocess/html/{sample}_R2.fastp.html'
        R2_json = 'results/01_preprocess/json/{sample}_R2.fastp.json'
    log:
        R1log = 
        R2log = 
    threads: 10
    conda: "envs/preprocess_env.yaml"
    shell: """
    fastp --detect_adadpter_for_pe \
        --correction --cut_right --thread {threads} \
        --html {output.R1_html} --json {output.R1_json} \
        --in1 {input.R1} --in2 {input.R2} \
        --out1 {output.R1} --out2 {output.R2} \
        2>{log.R1log}

    sed \
            's/{wildcards.sample}_{wildcards.unit}_R1/{wildcards.sample}_{wildcards.unit}_R2/g' \
            {log.R1log} > {log.R2log}
    """

# 1.2 quality assessment report for the reads
rule multiqc:
    ### we will use multiqc here to generate reports from the output of fastp.
    input:
        R1_json = expand('results/01_preprocess/json/{sample}_R1.fastp.json', sample=sample_df.sample_name),
        R2_json = expand('results/01_preprocess/json/{sample}_R2.fastp.json', sample=sample_df.sample_name)
    output:
        'results/01_preprocess/report/samples_report.html'
    log:
        multiqc_log = 
    shell: """
    multiqc -n samples_report.html \
        -o {output} \
        results/01_preprocess/json/ 2>{log}
    """


########## 2 Alignment ####################
"""
The final output for the alignment step would be the unsorted BAM files for all the included samples.
"""
rule alignment:
    input:
        expand('results/02_alignment/{sample}.unsorted.bam', sample=sample_df.sample_name)

# 2.1 downloading the human reference genome (GRCh37 - hg19)
## using hg19 because the QDNAseq in the later steps requires hg19 for generating CN profiles
rule download_hg19:
    ### if the hg19 reference genome does not exist, this rule will execute to download and generate the hg19 reference genome
    output:
        genome = 'resources/genome/hg19.ref.fa.gz'
    log:
        'resources/genome/download_hg19.log'
    shell: """
    wget 'https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.alt.fa.gz' -O hg19.ref.fa.gz 2>{log}
    """

# 2.2 indexing the hg19 reference genome
rule bwa_index:
    ### use bwa to index the reference genome
    input:
        genome = 'resources/genome/hg19.ref.fa.gz'
    output:
        multiext('hg19', ".amb", ".ann", ".bwt", ".pac", ".sa")
    conda: 'envs/alignment.yaml'
    threads: 
    shell: """
    bwa index -p hg19 {input.genome}
    """

# 2.3 mapping the reads to the indexed reference genome
rule map_reads:
    ### use bwa again for alignment
    input: 
        idx = rules.bwa_index.output,
        R1 = 'results/01_preprocess/reads/{sample}_{type_group}_R1.trimmed.fastq.gz',
        R2 = 'results/01_preprocess/reads/{sample}_{type_group}_R2.trimmed.fastq.gz'
    output:
        'results/02_alignment/{sample}.unsorted.sam'
    log: 
    params:
        index_ref = 'resources/genome/hg19.ref.fa.gz'
    conda: 'envs/alignment.yaml'
    threads: config['bwa_mapping']['threads']
    shell: """
    bwa mem -M -t {threads} \
        {params.idx} {input.R1} {input.R2} \
        -o {output} \
        2>{log}
    """


########## 3 Clean-up ####################
"""
The final output for the clean-up step should be the sorted, marked, and indexed BAM files.
"""
rule clean_up:
    input: 
        expand('results/03_clean_up/{sample}.clean.bam', sample=sample_df.sample_name)

# 3.1 sorting the SAM files
rule sort_sam: 
    ### using Picard to sort the sam files and remove the unsorted ones to save space
    input:
        'results/02_alignment/{sample}.unsorted.sam'
    output:
        'results/03_clean_up/{sample}.sorted.sam'
    log: 
    threads:
    conda: 'envs/clean_up.yaml'
    shell: """
    picard Sortsam \
        --INPUT {input} \
        --OUTPUT {output} \
        --SORTORDER coordinate \
        --VALIDATION_STRINGENCY SILENT
    
    rm {input}
    """

# 3.2 marking dupicates
rule de_duplicate:
    ### using Picard to remove PCR duplicates, and convert SAM file into BAM files
    input: 
        'results/03_clean_up/{sample}.sorted.sam'
    output:
        'results/03_clean_up/{sample}.sorted.dedup.bam'
    log:
    threads:
    params:
        metrix_file = 'results/03_clean_up/{sample}.metrics.txt'
    conda: 'envs/clean_up.yaml'
    shell: """
    picard MarkDuplicates \
        --INPUT {input} \
        --OUTPUT {output} \
        --METRICS_FILE {params.metrix_file} \
        --REMOVE_DUPLICATES true \
        --ASSUME_SORTED true \
        --VALIDATION_STRINGENCY SILENT
    
    rm {input}
    """


# 3.3 indexing the BAM files
rule index_bam:
    ### using samtools to index the bam files
    input:
        'results/03_clean_up/{sample}.sorted.dedup.bam'
    output:
        'results/03_clean_up/{sample}.clean.bam'
    log:
    threads:
    conda: 'envs/clean_up.yaml'
    shell: """
    samtools index -@ {threads} -o {output} {input}
    
    rm {input}
    """

########## 4 Relative CN profile ####################
"""
The final output for this step is the relative copy number (CN) profile generated by QDNAseq.
"""
rule relative_CN:
    input:
        rds = 'results/04_relative_CN/{sample}.rds',
        tsv = 'results/04_relative_CN/{sample}.tsv'
    
# 4.1 generating relative CN profile
rule QDNAseq:
    ### QDNAseq will be applied to generate relative copy number profile stored in RDS and tsv files for later analyses.
    input:
        bam = expand('results/03_clean_up/{sample}.clean.bam')
    output:
        rds = 'results/04_relative_CN/{sample}.rds',
        tsv = 'results/04_relative_CN/{sample}.tsv'
    params:
        bam_dir = 'results/03_clean_up/'
    threads:
    conda: 'envs/QDNAseq.yaml'
    script: 'script/QDNAseq.r'


########## 5 Absolute CN profile ####################
"""
The final output for this step and also the workflow would be the absolute copy number (CN) profile.
"""
rule absolute_CN:
    input:
        tsv = 'results/05_absoluteCN/{sample}.absoluteCN.tsv'

# 5.1 setting up the conda environment with rascal package downloaded from github
rule rascal_env:
    output:
        "other_info/rascal_settle_info.txt"
    conda: 'envs/rascal.yaml'
    script: 'scripts/rascal_env.R'


# 5.2 generating absolute copy number profiles based on the optimal solutions
rule rascal_absoluteCN:
    input:
        rds = 'results/04_relative_CN/{sample}.rds',
        env_set = "other_info/rascal_settle_info.txt"
    output:
        solution = 'results/05_absolute_CN/{sample}.solution.csv',
        solution_best = 'results/05_absolute_CN/{sample}.best_solution.csv',
        absolute_CN = 'results/05_absolute_CN/{sample}.absoluteCN.csv'
    params:
        output_prefix = 'results/05_absolute_CN/{sample}'
    threads:
    conda: 'envs/rascal.yaml'
    shell: '''
    Rscript /workflow/scripts/fit_absoluteCN.R -i {input.rds} -o {params.output_prefix} 
    '''



