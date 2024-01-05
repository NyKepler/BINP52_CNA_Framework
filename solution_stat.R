# Title: solution_stat.R
# Author: Guyuan TANG
# Date: 2023/12/4

library(tidyverse)
library(readxl)
library(ggsignif)
library(ggpubr)

# select the binsize (unit: kb) 
binsize <- c('5kb','10kb','15kb','30kb','50kb','100kb')

# sample types
sample_groups <- c('ArchivalVS', 'Blood', 'Endome', 'ffpe', 'ffPlasma', 'ffTissue', 'ffTumor', 'MaNiLaVS')


# read the dataset
solution_table <- read_excel("MaNiLa_All_samplesheet_test.xlsx", sheet = 'rascal_5kb')


# generate the whole table
for (bin in binsize) {
  if (bin != '5kb'){
    df_bin <- read_excel("MaNiLa_All_samplesheet_test.xlsx", sheet = paste0('rascal_',bin))
    solution_table <- left_join(solution_table, df_bin, by=c('Library'='Library','Patient'='Patient', 'Type'='Type', 'Sample'='Sample', 'Fastq'='Fastq', 'Group'='Group', 'num_reads'='num_reads'))
  }
}

# stats of number of reads in different sample types
solution_table %>%
  group_by(Group) %>%
  summarise(
    mean_read = mean(num_reads/1000000),
    sd_read   = sd(num_reads/1000000)
  )
reads_plot <- ggplot(data=solution_table, aes(x=Group, y=num_reads/1000000)) +
  stat_boxplot(geom = "errorbar", width=0.1,linewidth=0.3) +
  geom_boxplot(fill="grey") +
  labs(x="Sample types",
       y="Number of reads (M)") +
  theme_classic()
ggsave('mean_reads.pdf', plot = reads_plot, dpi = 600, width = 10, height = 7, units = 'in')


# stats of number of segments for each sample types in different bin size
seg_df <- select(solution_table, "Group", "num_segments_5kb", "num_segments_10kb", "num_segments_15kb",
                 "num_segments_30kb", "num_segments_50kb", "num_segments_100kb") %>%
  rename(c('5'='num_segments_5kb', '10'='num_segments_10kb', '15'='num_segments_15kb',
               '30'='num_segments_30kb', '50'='num_segments_50kb', '100'='num_segments_100kb')) %>%
  pivot_longer(!Group, names_to = "binsize", values_to="segments")

a <- aggregate(seg_df$segments, by=list(seg_df$Group, seg_df$binsize), median) %>% 
  rename(c('Group' = 'Group.1', 'Binsize' = 'Group.2', 'Median' = 'x'))
b <- aggregate(seg_df$segments, by=list(seg_df$Group, seg_df$binsize), IQR) %>% 
  rename(c('Group' = 'Group.1', 'Binsize' = 'Group.2', 'IQR' = 'x'))
a <- left_join(a,b, by=c('Group'='Group', 'Binsize'='Binsize'))
print(a)

seg_plot <-  ggboxplot(seg_df, x='binsize', y='segments', 
                       color = 'black', fill = 'Group', palette = 'Paired', 
                       xlab = 'Bin size (kb)', ylab = 'Number of segments') +
  ylim(0,800) +
  theme_classic()
ggsave(paste0('median_segs.pdf'), plot = seg_plot, dpi = 600, width = 10, height = 7, units = 'in')


# stats of number of solutions in different sample types
solution_df <- select(solution_table, "Group", "rascal_5kb_num", "rascal_10kb_num", "rascal_15kb_num",
                      "rascal_30kb_num", "rascal_50kb_num", "rascal_100kb_num") %>%
  rename(c('5'='rascal_5kb_num', '10'='rascal_10kb_num', '15'='rascal_15kb_num',
           '30'='rascal_30kb_num', '50'='rascal_50kb_num', '100'='rascal_100kb_num')) %>%
  pivot_longer(!Group, names_to = "binsize", values_to="solutions")

a <- aggregate(solution_df$solutions, by=list(solution_df$Group, solution_df$binsize), median) %>% 
  rename(c('Group' = 'Group.1', 'Binsize' = 'Group.2', 'Median' = 'x'))
b <- aggregate(solution_df$solutions, by=list(solution_df$Group, solution_df$binsize), IQR) %>% 
  rename(c('Group' = 'Group.1', 'Binsize' = 'Group.2', 'IQR' = 'x'))
a <- left_join(a,b, by=c('Group'='Group', 'Binsize'='Binsize'))
print(a)

solution_plot <-  ggboxplot(solution_df, x='binsize', y='solutions', 
                       color = 'black', fill = 'Group', palette = 'Paired', 
                       xlab = 'Bin size (kb)', ylab = 'Number of solutions') +
  theme_classic()
ggsave('median_solutions.pdf', plot = solution_plot, dpi = 600, width = 10, height = 7, units = 'in')


# stats of median ploidy
ploidy_df <- select(solution_table,"Group", "rascal_5kb_ploidy_1", "rascal_10kb_ploidy_1", "rascal_15kb_ploidy_1",
                                   "rascal_30kb_ploidy_1", "rascal_50kb_ploidy_1", "rascal_100kb_ploidy_1") %>%
  rename(c('5'='rascal_5kb_ploidy_1', '10'='rascal_10kb_ploidy_1', '15'='rascal_15kb_ploidy_1',
           '30'='rascal_30kb_ploidy_1', '50'='rascal_50kb_ploidy_1', '100'='rascal_100kb_ploidy_1')) %>%
  pivot_longer(!Group, names_to = "binsize", values_to="ploidy")
ploidy_df$ploidy[which(ploidy_df$ploidy == -1)] <- NA

# calculate the number of zero solutions
bins <- c('5','10','15','30','50','100')
for (samples in sample_groups) {
  for (bin in bins) {
    a <- ploidy_df %>% filter(Group == samples, binsize == bin) %>% summarise_all(~sum(is.na(.)))
    cat('----- ',samples, '---',bin,'kb -----\n')
    print(a)
  }
}
sample_list <- c("ArchivalVS","ArchivalVS","ArchivalVS","ArchivalVS","ArchivalVS","ArchivalVS",
                 "Blood","Blood","Blood","Blood","Blood", "Blood",
                 "Endome","Endome","Endome","Endome","Endome","Endome",
                 "ffpe","ffpe","ffpe","ffpe","ffpe","ffpe",
                 "ffPlasma","ffPlasma","ffPlasma","ffPlasma","ffPlasma","ffPlasma",
                 "ffTissue","ffTissue","ffTissue","ffTissue","ffTissue","ffTissue",
                 "ffTumor","ffTumor","ffTumor","ffTumor","ffTumor","ffTumor",
                 "MaNiLaVS","MaNiLaVS","MaNiLaVS","MaNiLaVS","MaNiLaVS","MaNiLaVS")
bin_list <- c('5kb','10kb','15kb','30kb','50kb','100kb','5kb','10kb','15kb','30kb','50kb','100kb',
              '5kb','10kb','15kb','30kb','50kb','100kb','5kb','10kb','15kb','30kb','50kb','100kb',
              '5kb','10kb','15kb','30kb','50kb','100kb','5kb','10kb','15kb','30kb','50kb','100kb',
              '5kb','10kb','15kb','30kb','50kb','100kb','5kb','10kb','15kb','30kb','50kb','100kb')
num_list <- c(3,3,4,2,0,0,
              0,0,0,0,0,0,
              0,0,0,0,0,0,
              6,6,6,6,4,0,
              0,0,0,0,0,1,
              0,0,0,0,0,0,
              8,7,6,0,0,0,
              0,0,0,0,0,0)
test_df <- data.frame(sample_list,bin_list,num_list) %>% 
  rename(c('Group'='sample_list', 'Binsize'='bin_list', 'num'='num_list')) 

test_df$Binsize <- factor(test_df$Binsize, levels = c('5kb', '10kb', '15kb', '30kb', '50kb', '100kb'))

test_plot <- ggplot(test_df, aes(x = Binsize, y = num, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5)) +
  labs(x = "Binsize",
       y = "Number of Samples without Solutions") +
  scale_fill_brewer(palette = 'Paired') +
  theme_classic()
ggsave('0_solutions.pdf', plot = test_plot, dpi = 600, width = 10, height = 7, units = 'in')

# stats of ploidy
a <- ploidy_df %>% drop_na(ploidy)
b <- aggregate(a$ploidy, by=list(a$Group, a$binsize), median) %>% 
  rename(c('Group' = 'Group.1', 'Binsize' = 'Group.2', 'Median' = 'x'))
c <- aggregate(a$ploidy, by=list(a$Group, a$binsize), IQR) %>% 
  rename(c('Group' = 'Group.1', 'Binsize' = 'Group.2', 'IQR' = 'x'))
a <- left_join(b,c, by=c('Group' = 'Group', 'Binsize' = 'Binsize'))
print(a)
# plotting
ploidy_df$ploidy <- log2(ploidy_df$ploidy)
for (samples in sample_groups) {
  ploidy_plot <- ggboxplot(filter(ploidy_df, Group == samples), x='binsize', y='ploidy', 
                           color = 'black', fill = 'grey', palette = 'Paired') +
    labs(title = paste0(samples, ' Ploidy'), x = 'Bin size (kb)', y = 'Ploidy (log2)') +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(samples, '_ploidy.pdf'), plot = ploidy_plot, dpi = 600, width = 10, height = 7, units = 'in')
}




# Cellularity
cellularity_df <- select(solution_table,"Group", "rascal_5kb_cellularity_1", "rascal_10kb_cellularity_1", 
                         "rascal_15kb_cellularity_1","rascal_30kb_cellularity_1", "rascal_50kb_cellularity_1",
                         "rascal_100kb_cellularity_1") %>%
  rename(c('5'='rascal_5kb_cellularity_1', '10'='rascal_10kb_cellularity_1', 
           '15'='rascal_15kb_cellularity_1', '30'='rascal_30kb_cellularity_1', 
           '50'='rascal_50kb_cellularity_1', '100'='rascal_100kb_cellularity_1')) %>%
  pivot_longer(!Group, names_to = "binsize", values_to="cellularity")
cellularity_df$cellularity[which(cellularity_df$cellularity == -1)] <- NA
cellularity_df$cellularity[which(cellularity_df$cellularity == 1)] <- 0
cellularity_df$cellularity <- cellularity_df$cellularity * 100
# stats of cellularity
a <- cellularity_df %>% drop_na(cellularity)# stats of ploidy
a <- cellularity_df %>% drop_na(cellularity)
b <- aggregate(a$cellularity, by=list(a$Group, a$binsize), median) %>% 
  rename(c('Group' = 'Group.1', 'Binsize' = 'Group.2', 'Median' = 'x'))
c <- aggregate(a$cellularity, by=list(a$Group, a$binsize), IQR) %>% 
  rename(c('Group' = 'Group.1', 'Binsize' = 'Group.2', 'IQR' = 'x'))
a <- left_join(b,c, by=c('Group' = 'Group', 'Binsize' = 'Binsize'))
print(a)
# plotting
for (samples in sample_groups) {
  ploidy_plot <- ggboxplot(filter(cellularity_df, Group == samples), x='binsize', y='cellularity', 
                           color = 'black', fill = 'grey', palette = 'Paired') +
    labs(title = paste0(samples, ' cellularity'), x = 'Bin size (kb)', y = 'Cellularity (%)') +
    ylim(0,100) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0('cellularity/',samples, '_cellularity.pdf'), plot = ploidy_plot, dpi = 600, width = 10, height = 7, units = 'in')
}


# MAD
MAD_df <- select(solution_table,"Group", "rascal_5kb_MAD_1", "rascal_10kb_MAD_1", 
                         "rascal_15kb_MAD_1","rascal_30kb_MAD_1", "rascal_50kb_MAD_1",
                         "rascal_100kb_MAD_1") %>%
  rename(c('5'='rascal_5kb_MAD_1', '10'='rascal_10kb_MAD_1', 
           '15'='rascal_15kb_MAD_1', '30'='rascal_30kb_MAD_1', 
           '50'='rascal_50kb_MAD_1', '100'='rascal_100kb_MAD_1')) %>%
  pivot_longer(!Group, names_to = "binsize", values_to="MAD")
MAD_df$MAD[which(MAD_df$MAD == -1)] <- NA
MAD_df$MAD <- MAD_df$MAD * 100
# stats of MAD
a <- MAD_df %>% drop_na(MAD)
a <- MAD_df %>% drop_na(MAD)
b <- aggregate(a$MAD, by=list(a$Group, a$binsize), mean) %>% 
  rename(c('Group' = 'Group.1', 'Binsize' = 'Group.2', 'Mean' = 'x'))
c <- aggregate(a$MAD, by=list(a$Group, a$binsize), sd) %>% 
  rename(c('Group' = 'Group.1', 'Binsize' = 'Group.2', 'sd' = 'x'))
a <- left_join(b,c, by=c('Group' = 'Group', 'Binsize' = 'Binsize'))
print(a)
# plotting
for (samples in sample_groups) {
  MAD_plot <- ggboxplot(filter(MAD_df, Group == samples), x='binsize', y='MAD', 
                           color = 'black', fill = 'grey', palette = 'Paired') +
    labs(title = paste0(samples, ' Distance (MAD)'), x = 'Bin size (kb)', y = bquote('MAD '(10^-2))) +
    ylim(0,100) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0('MAD/',samples, '_MAD.pdf'), plot = MAD_plot, dpi = 600, width = 10, height = 7, units = 'in')
}


