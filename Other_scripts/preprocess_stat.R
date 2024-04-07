# Title: preprocess_stat.R
# Author: Guyuan TANG
# Date: 2024/03/18

# Description: this was designed for analysing the performance of preprocessing and clean-up steps.

library(tidyverse)
library(readxl)

# Load the sample preprocessing data
sample_df <- read_excel('preprocess_stats.xlsx',sheet = 'stat')

# number of reads
sample_df %>% 
  group_by(Group) %>%
  summarise(
    mean_read_0 = mean(num_read_0),
    sd_read_0 = sd(num_read_0),
    mean_read_1 = mean(num_read_1),
    sd_read_1 = sd(num_read_1),
    mean_read_2 = mean(num_read_2),
    sd_read_2 = sd(num_read_2)
  )
read_df <- sample_df %>% select('Sample','num_read_0', 'num_read_1', 'num_read_2') %>%
  rename(c('raw'='num_read_0','preprocess'='num_read_1','clean-up'='num_read_2')) %>%
  pivot_longer(!Sample, names_to = 'Status', values_to = 'Number_of_reads')
read_df$Group <- NA
for (i in read_df$Sample){
  read_df[which(read_df$Sample==i),'Group'] <- sample_df[which(sample_df$Sample==i),'Group']
}
read_df$Status <- factor(read_df$Status, levels = c('raw','preprocess','clean-up'))
read_df$Group <- factor(read_df$Group, levels = c('ArchivalVS', 'MaNiLaVS', 'ffTumor', 'ffTissue', 'Endome', 'ffpe', 'ffPlasma', 'Blood'))

reads_plot <- ggplot(data=read_df, aes(x=Status, y=Number_of_reads)) +
  stat_boxplot(geom = "errorbar", width=0.1,linewidth=0.3) +
  geom_boxplot(fill="grey") +
  labs(x="Status",
       y="Number of reads (M)") +
  theme_classic() +
  facet_grid(cols = vars(Group))
reads_plot <- reads_plot + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('figures/mean_reads_all.pdf', plot = reads_plot, dpi = 600, width = 10, height = 7, units = 'in')


# Q30
q30_df <- sample_df %>% select('Sample', 'Q30_base_0', 'Q30_base_1') %>%
  rename(c('raw'='Q30_base_0','preprocess'='Q30_base_1')) %>%
  pivot_longer(!Sample, names_to = 'Status', values_to = 'Percentage_of_Q30_bases')
q30_df$Group <- NA
for (i in q30_df$Sample){
  q30_df[which(q30_df$Sample==i),'Group'] <- sample_df[which(sample_df$Sample==i),'Group']
}
q30_df$Status <- factor(q30_df$Status, levels = c('raw','preprocess'))
q30_df$Group <- factor(q30_df$Group, levels = c('ArchivalVS', 'MaNiLaVS', 'ffTumor', 'ffTissue', 'Endome', 'ffpe', 'ffPlasma', 'Blood'))
q30_plot <- ggplot(data=q30_df, aes(x=Status, y=Percentage_of_Q30_bases)) +
  stat_boxplot(geom = "errorbar", width=0.1,linewidth=0.3) +
  geom_boxplot(fill="grey") +
  labs(x="Status",
       y="Percentage of Q30 bases (%)") +
  theme_classic() +
  facet_grid(cols = vars(Group))
q30_plot <- q30_plot + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('figures/mean_q30_all.pdf', plot = q30_plot, dpi = 600, width = 10, height = 7, units = 'in')


# mapping percentage
sample_df$Group <- factor(sample_df$Group, levels = c('ArchivalVS', 'MaNiLaVS', 'ffTumor', 'ffTissue', 'Endome', 'ffpe', 'ffPlasma', 'Blood'))
map_plot <- ggplot(data=sample_df, aes(x=Group, y=map_2)) +
  stat_boxplot(geom = "errorbar", width=0.1,linewidth=0.3) +
  geom_boxplot(fill="grey") +
  labs(x="Sample types",
       y="Percentage of mapped reads (%)") +
  theme_classic()
ggsave('figures/mean_map_all.pdf', plot = map_plot, dpi = 600, width = 10, height = 7, units = 'in')



# summarise of other information
sample_df %>% 
  group_by(Group) %>%
  summarise(
    length_0 = mean(mean_length_0),
    length_1 = mean(mean_length_1),
    dup_0 = mean(duplication_rate_0),
    mean_map = mean(map_2),
    sd_map = sd(map_2),
    Q30_0 = mean(Q30_base_0),
    Q30_1 = mean(Q30_base_1)
  )
