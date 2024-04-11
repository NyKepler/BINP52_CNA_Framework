#!/usr/bin/bash

# Title: correct_ref_chr.sh
# Author: Guyuan TANG
# Date: 2024/1/21

# Description: this bash script is used for adding 'chr' to the chromosome names in the reference files for hg19. Notably, we only change the gc_correction files and replication files. Other reference files could be adjusted according to the method the developer provided online (https://github.com/cancerit/cgpBattenberg/wiki/Reference-name-conventions).

# Step:
## 1. unzip all the reference files
## 2. make adjustment on the files
## 3. gzip all the adjusted files

## the working directory should be ~/resources/battenberg (where you downloaded the references files)

# the gc correction files
tar -zxvf battenberg_wgs_gc_correction_1000g_v3.tar.gz
for i in {1..23}
do
    filename=battenberg_wgs_gc_correction_1000g_v3/1000_genomes_GC_corr_chr_${i}.txt
    gunzip ${filename}.gz
    # change the column orders
    echo -n "$(awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' ${filename})" > ${filename}
    # add 'chr'
    sed -i s'/^/chr&/g' ${filename}
    # change the headers
    sed -i s'/chrPosition chr/chr Position/' ${filename}
    gzip ${filename}
done

# the replication files
tar -zxvf battenberg_wgs_replic_correction_1000g_v3.tar.gz
for i in {1..23}
do
    filename=battenberg_wgs_replic_correction_1000g_v3/1000_genomes_replication_timing_chr_${i}.txt
    gunzip ${filename}.gz
    # add 'chr'
    sed -i s'/^/chr&/g' ${filename}
    # change the header
    sed -i s'/chrchr/chr/' ${filename}
    gzip ${filename}
done

