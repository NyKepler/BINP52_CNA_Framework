# Title: Snakefile_draft.smk
# Author: Guyuan TANG
# Date: 2023/10/16 - 

import pandas as pd

# specify the configuration file
configfile: "config/config.yaml"

# specify the samples and their groups
sample_df = (pd.read_csv(config['samples'], 
    sep='\t', 
    dtype={'sample_name':str, 'patient':str, 'fastq_1':str, 'fastq_2':str})
    .set_index('sample_name', drop=False))

sample_group = config['type']

# specify the results location (output directory)
results = config['outputdir']

# specify the final output of the whole workflow
"""
the final output for the whole workflow should be the absolute copy number profile for each type/group of samples.
For example: archive (pre-diagnosis), diagnosis, tumor tissue.
"""
rule all:
    input:
        expand(results + '05_absolute_CN/{sample}/{sample}.solution.csv', sample=sample_df.sample_name)


########## 1 Preprocessing ####################
"""
The final output for the step preprocessing should be a multiqc report generating all the preprocessed QC statistics.
"""
rule preprocess:
    input:
        results + '01_preprocess/html/' + sample_group + '_multiqc_report.html'

# 1.1 quality trimming on reads
rule fastp:
    ### we will use fastp here for trimming on adapter and quality.
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

# 1.2 quality assessment of preprocessed reads with fastqc
rule fastqc:
    ### we will use fastqc to generate the quality control stats from the outputs of fastp
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

# 1.3 quality assessment report for the reads
rule multiqc:
    ### we will use multiqc here to generate reports from the output of fastqc.
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


########## 2 Alignment ####################
"""
The final output for the alignment step would be the unsorted BAM files for all the included samples.
"""
rule alignment:
    input:
        expand(results + '02_alignment/{sample}.unsorted.sam', sample=sample_df.sample_name)

# 2.1 downloading the human reference genome (GRCh37 - hg19)
## using hg19 because the QDNAseq in the later steps requires hg19 for generating CN profiles
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

# 2.2 indexing the hg19 reference genome
rule bwa_index:
    ### use bwa to index the reference genome
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


# 2.3 mapping the reads to the indexed reference genome
rule map_reads:
    ### use bwa again for alignment
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


########## 3 Clean-up ####################
"""
The final output for the clean-up step should be the sorted, marked, and indexed BAM files.
"""
rule clean_up:
    input: 
        expand(results + '03_clean_up/{sample}/{sample}.sorted.dedup.bai', sample=sample_df.sample_name)

# 3.1 sorting the SAM files
rule sort_sam: 
    ### using Picard to sort the sam files 
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

# 3.2 marking dupicates and de-duplication
rule de_duplicate:
    ### using Picard to remove PCR duplicates, and convert SAM file into BAM files
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


# 3.3 indexing the BAM files
rule index_bam:
    ### using samtools to show the stats of the sorted and deduplicates outputs and to index the BAM files
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

# to test the quality of the BAM files: using qualimap
# qualimap bamqc -bam results/03_clean_up/{sample}/{sample}.sorted.dedup.bam --java-mem-size=4G

########## 4 Relative CN profile ####################
"""
The final output for this step is the relative copy number (CN) profile generated by QDNAseq.
"""
rule relative_CN:
    input:
        rds = expand(results + '04_relative_CN/{sample}/{sample}.rds', sample=sample_df.sample_name),
        tsv = expand(results + '04_relative_CN/{sample}/{sample}.seg.tsv', sample=sample_df.sample_name)
    
# 4.1 generating relative CN profile
rule QDNAseq:
    ### QDNAseq will be applied to generate relative copy number profile stored in RDS and tsv files for later analyses.
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


########## 5 Absolute CN profile ####################
"""
The final output for this step and also the workflow would be the absolute copy number (CN) profile.
"""
rule absolute_CN:
    input:
        tsv = results + '05_absoluteCN/{sample}.absoluteCN.tsv'

# 5.1 setting up the conda environment with rascal package downloaded from github
rule rascal_env:
    input: rules.relative_CN.input
    output:
        results + "other_info/rascal_settle_info.txt"
    conda: 'envs/rascal.yaml'
    script: 'scripts/rascal_env.R'


# 5.2 generating absolute copy number profiles based on the optimal solutions
rule rascal_absoluteCN:
    input:
        rds = results + '04_relative_CN/{sample}.rds',
        env_set = results + "other_info/rascal_settle_info.txt"
    output:
        solution = results + '05_absolute_CN/{sample}.solution.csv',
        solution_best = results + '05_absolute_CN/{sample}.best_solution.csv',
        absolute_CN = results + '05_absolute_CN/{sample}.absoluteCN.csv'
        absolute_CN_seg = results + '05_absolute_CN/{sample}.absoluteCN.seg.csv'
    params:
        output_prefix = results + '05_absolute_CN/{sample}'
    threads:
    conda: 'envs/rascal.yaml'
    shell: '''
    Rscript /workflow/scripts/fit_absoluteCN.R -i {input.rds} -o {params.output_prefix} 
    '''



