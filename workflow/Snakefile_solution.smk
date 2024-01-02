# Title: Snakefile_solution.smk
# Author: Guyuan TANG
# Date: 2023/10/16 - 2023/11/30

import pandas as pd

# specify the configuration file
configfile: "config/config.yaml"

# specify the samples and their groups
sample_df = (pd.read_csv(config['samples'], 
    sep='\t', 
    dtype={'sample_name':str, 'patient':str, 'type':str, 'fastq_1':str, 'fastq_2':str})
    .set_index('sample_name', drop=False))

# specify the results location (output directory)
results = config['outputdir']

# specify the bin size used for annotation
binsize = str(config['QDNAseq']['binsize'])

# specify the final output of the whole workflow
"""
the final output for the whole workflow should be the estimated solutions of ploidy and cellularity for each sample.
"""
rule all:
    input:
        expand(results + '{sample}/05_absolute_CN/solutions/{sample}_' + binsize + 'kb.solution.csv', sample=sample_df.sample_name)


########## 1 Preprocessing ####################
"""
The final output for the step preprocessing should be a multiqc report generating all the preprocessed QC statistics.
"""
rule preprocess:
    input:
        results + 'fastqc/preprocess_multiqc_report.html'

# 1.1 quality assessment on raw reads
rule fastqc_0:
    ### we will use fastqc to assess the quality of raw reads first
    input:
        R1_seq = lambda wildcards: sample_df.loc[wildcards.sample, 'fastq_1'],
        R2_seq = lambda wildcards: sample_df.loc[wildcards.sample, 'fastq_2']
    output:
        R1_html = results + '{sample}/01_preprocess/html/{sample}_R1_raw_fastqc.html',
        R1_qc = results + 'fastqc/raw/{sample}_R1_raw_fastqc.zip',
        R2_html = results + '{sample}/01_preprocess/html/{sample}_R2_raw_fastqc.html',
        R2_qc = results + 'fastqc/raw/{sample}_R2_raw_fastqc.zip'
    log: 'log/fastqc/{sample}.raw.fastqc.log'
    params: 
        outdir = results + 'fastqc/raw/',
        R1_qc = results + 'fastqc/raw/{sample}_*R1_*fastqc.zip',
        R2_qc = results + 'fastqc/raw/{sample}_*R2_*fastqc.zip',
        R1_html = results + 'fastqc/raw/{sample}_*R1_*fastqc.html',
        R2_html = results + 'fastqc/raw/{sample}_*R2_*fastqc.html'
    threads: 2
    conda: 'envs/preprocess_env.yaml'
    shell: """
    fastqc -o {params.outdir} {input.R1_seq} {input.R2_seq} 2>{log}
    mv {params.R1_html} {output.R1_html}
    mv {params.R2_html} {output.R2_html}
    mv {params.R1_qc} {output.R1_qc}
    mv {params.R2_qc} {output.R2_qc}
    """

# 1.2 generate multiqc report on the raw reads
rule multiqc_0:
    ### we will use multiqc here to generate reports from the output of fastqc.
    input:
        R1_qc = expand(results + 'fastqc/raw/{sample}_R1_raw_fastqc.zip', sample=sample_df.sample_name),
        R2_qc = expand(results + 'fastqc/raw/{sample}_R2_raw_fastqc.zip', sample=sample_df.sample_name)
    output:
        results + 'fastqc/raw_multiqc_report.html'
    log:
        'log/raw_multiqc.log'
    conda: 'envs/multiqc_env.yaml'
    params: 
        out_name = 'raw_multiqc_report.html',
        indir = results + 'fastqc/raw/',
        outdir = results + 'fastqc/'
    shell: """
    multiqc -f -n {params.out_name} \
    -o {params.outdir} {params.indir} >{log} 2>{log}
    rm -r {params.outdir}/raw_multiqc_report_data/
    """


# 1.3 quality trimming on reads
rule fastp:
    ### we will use fastp here for trimming on adapter and quality.
    input:
        qc_raw = results + 'fastqc/raw_multiqc_report.html',
        R1 = lambda wildcards: sample_df.loc[wildcards.sample, 'fastq_1'],
        R2 = lambda wildcards: sample_df.loc[wildcards.sample, 'fastq_2']
    output:
        R1 = results + '{sample}/01_preprocess/reads/{sample}_R1_preprocess.fastq.gz',
        html = results + '{sample}/01_preprocess/html/{sample}_fastp.html',
        R2 = results + '{sample}/01_preprocess/reads/{sample}_R2_preprocess.fastq.gz'
    log: 'log/fastp/{sample}_fastp.log'
    threads: 16
    params: json = results + '{sample}/01_preprocess/html/{sample}_fastp.json'
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

# 1.4 quality assessment of preprocessed reads with fastqc
rule fastqc_1:
    ### we will use fastqc to generate the quality control stats from the outputs of fastp
    input:
        R1_seq = results + '{sample}/01_preprocess/reads/{sample}_R1_preprocess.fastq.gz',
        R2_seq = results + '{sample}/01_preprocess/reads/{sample}_R2_preprocess.fastq.gz'
    output:
        R1_html = results + '{sample}/01_preprocess/html/{sample}_R1_preprocess_fastqc.html',
        R1_qc = results + 'fastqc/preprocess/{sample}_R1_preprocess_fastqc.zip',
        R2_html = results + '{sample}/01_preprocess/html/{sample}_R2_preprocess_fastqc.html',
        R2_qc = results + 'fastqc/preprocess/{sample}_R2_preprocess_fastqc.zip'
    log: 'log/fastqc/{sample}.preprocess.fastqc.log'
    params: 
        outdir = results + 'fastqc/preprocess/',
        R1_html = results + 'fastqc/preprocess/{sample}_R1_preprocess_fastqc.html',
        R2_html = results + 'fastqc/preprocess/{sample}_R2_preprocess_fastqc.html'
    threads: 2
    conda: 'envs/preprocess_env.yaml'
    shell: """
    fastqc -o {params.outdir} {input.R1_seq} {input.R2_seq} 2>{log}
    mv {params.R1_html} {output.R1_html}
    mv {params.R2_html} {output.R2_html}
    """

# 1.5 quality assessment report for the reads
rule multiqc_1:
    ### we will use multiqc here to generate reports from the output of fastqc.
    input:
        R1_qc = expand(results + 'fastqc/preprocess/{sample}_R1_preprocess_fastqc.zip', sample=sample_df.sample_name),
        R2_qc = expand(results + 'fastqc/preprocess/{sample}_R2_preprocess_fastqc.zip', sample=sample_df.sample_name)
    output:
        results + 'fastqc/preprocess_multiqc_report.html'
    log:
        'log/preprocess_multiqc.log'
    conda: 'envs/multiqc_env.yaml'
    params: 
        out_name = 'preprocess_multiqc_report.html',
        indir = results + 'fastqc/preprocess/',
        outdir = results + 'fastqc/'
    shell: """
    multiqc -f -n {params.out_name} \
    -o {params.outdir} {params.indir} >{log} 2>{log}
    rm -r {params.outdir}/preprocess_multiqc_report_data/
    """


########## 2 Alignment ####################
"""
The final output for the alignment step would be the unsorted BAM files for all the included samples.
"""
rule alignment:
    input:
        expand(results + '{sample}/02_alignment/{sample}.unsorted.sam', sample=sample_df.sample_name)

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

    rm resources/genome/*.fa

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
    threads: 20
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
        R1 = results + '{sample}/01_preprocess/reads/{sample}_R1_preprocess.fastq.gz',
        R2 = results + '{sample}/01_preprocess/reads/{sample}_R2_preprocess.fastq.gz'
    output:
        results + '{sample}/02_alignment/{sample}.unsorted.sam'
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
        expand(results + '{sample}/03_clean_up/{sample}.sorted.dedup.bai', sample=sample_df.sample_name)

# 3.1 sorting the SAM files
rule sort_sam_dedup: 
    ### using Picard to sort the sam files, to mark and remove the PCR duplicates, and to convert SAM into BAM
    input:
        sam = results + '{sample}/02_alignment/{sample}.unsorted.sam'
    output:
        results + '{sample}/03_clean_up/{sample}.sorted.dedup.bam'
    log: 
        sort_sam = 'log/sort_sam/{sample}.log',
        de_duplicate = 'log/de_duplicate/{sample}.log'
    conda: 'envs/clean_up.yaml'
    threads: 2
    params: 
        metrix_file = results + '{sample}/03_clean_up/{sample}.metrics.txt',
        sorted_sam = results + '{sample}/03_clean_up/{sample}.sorted.sam'
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


# 3.2 indexing the BAM files
rule index_bam:
    ### using samtools to show the stats of the sorted and deduplicates outputs and to index the BAM files
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

# to test the quality of the BAM files: using qualimap
# qualimap bamqc -bam results/03_clean_up/{sample}/{sample}.sorted.dedup.bam --java-mem-size=4G

########## 4 Relative CN profile ####################
"""
The final output for this step is the relative copy number (CN) profile generated by QDNAseq.
"""
rule relative_CN:
    input:
        rds = expand(results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.rds', sample=sample_df.sample_name),
        tsv = expand(results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.seg.tsv', sample=sample_df.sample_name)
    
# 4.1 generating relative CN profile
rule QDNAseq:
    ### QDNAseq will be applied to generate relative copy number profile stored in RDS and tsv files for later analyses.
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


########## 5 Ploidy and cellularity solution ################
"""
The final output for this step would be the solutions of ploidy and tumour purity for each sample. It is also the last step in the first snakemake file (Snakefile_solution.smk)
"""
rule CN_solution: # the rule is the same with rule all at this moment
    input:
        expand(results + '{sample}/05_absolute_CN/solutions/{sample}_' + binsize + 'kb.solution.csv', sample=sample_df.sample_name)

# 5.1 setting up the conda environment with rascal package downloaded from github
rule rascal_env:
    output:
        "log/rascal_settle_info.txt"
    conda: 'envs/rascal.yaml'
    script: 'scripts/rascal_env.R'

# 5.2 calculate the optimal solutions (ploidy and cellularity) of the samples
rule rascal_solution:
    ### Rascal will be applied to calculate the optimal solutions for further deriving the absolute copy numbers
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



