# Title: analysis_signatures.R
# Author: Guyuan TANG
# Date: 2024/1/10

# Description: the script was designed for analyzation on the extracted signatures.

setwd('/home/researcher/TangGY/BINP52/Workflow/Draft/')

########## Clustering ##########
library(readxl)
library(tidyverse)

sample_df <- read.table('config/solution_sample.tsv', sep = '\t', header = 1) %>% select('Sample', 'Patient', 'Type', 'Group')
sample_df['enrich_CN'] <- NA
sample_df['CN_similarity'] <- NA
sample_df['enrich_PanCan'] <- NA
sample_df['PanCan_similarity'] <- NA

# 1. CN signatures
## load the similarity data
CN_mat <- readRDS('results/signatures/CN_sig/CN_sig.SSmatrix.rds')
for (i in 1:nrow(sample_df)) {
  sampleID <- sample_df[i, 'Sample']
  enrich_sig <- CN_mat[which(CN_mat$sample==sampleID), 'enrich']
  similarity <- CN_mat[which(CN_mat$sample==sampleID), enrich_sig]
  # add to the whole sample table
  sample_df[i, 'enrich_CN'] <- enrich_sig
  sample_df[i, 'CN_similarity'] <- similarity
}
# 2. Pan-Cancer signatures
## load the similarity data
PC_mat <- readRDS('results/signatures/PanCan_sig/PanCan_sig.SSmatrix.rds')
for (i in 1:nrow(sample_df)) {
  sampleID <- sample_df[i, 'Sample']
  ## for pan-cancer signatures, some samples with low segment counts will not be computed
  if (sampleID %in% PC_mat$sample) {
    enrich_sig <- PC_mat[which(PC_mat$sample==sampleID), 'enrich']
    similarity <- PC_mat[which(PC_mat$sample==sampleID), enrich_sig]
    # add to the whole sample table
    sample_df[i, 'enrich_PanCan'] <- enrich_sig
    sample_df[i, 'PanCan_similarity'] <- similarity
  }
}
# save the clustering data into excel file
write.table(sample_df, file = 'results/signatures/All_sample_signature.csv', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ',')


########## Analyzation on clusters ##########
# sample_df <- read.csv('results/signatures/All_sample_signature.csv')

# 1. CN signatures
groups <- c('ArchivalVS', 'Blood', 'Endome', 'ffpe', 'ffPlasma', 'ffTissue', 'ffTumor', 'MaNiLaVS')
## if you only want to analyze the three groups, you can alter the group names in the above vector
## groups <- c('ArchivalVS','ffTumor','MaNiLaVS')
CN_df <- as.data.frame(matrix(ncol = 6, nrow = 0))
colnames(CN_df) <- c('Group', 'enrich_CN', 'mean_CN_similarity', 'sd_CN_similarity', 'sample_size', 'percentage')
for (groupID in groups) {
  info_df <- sample_df %>% 
            filter(Group == groupID) %>%
            group_by(enrich_CN) %>%
            summarise(mean_CN_similarity = mean(CN_similarity),
                      sd_CN_similarity = sd(CN_similarity),
                      sample_size = n())
  total_n <- sum(info_df$sample_size)
  for (i in 1:nrow(info_df)) {
    info_df[i, 'percentage'] = info_df[i, 'sample_size']/total_n*100
  }
  info_df['Group'] <- groupID
  info_df <- select(info_df, 'Group', 'enrich_CN':'percentage')
  cat(paste0(groupID, '\n'))
  print(info_df)
  CN_df <- rbind(CN_df, info_df)
}
## make the bar plot
p1 <- ggplot(data = CN_df, aes(x=Group, y=percentage, fill=enrich_CN)) +
      geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
      labs(x='Groups', y='Percentage (%)', fill='CN_signatures') +
      scale_fill_brewer(palette = 'Set2') +
      theme_classic()
ggsave('results/signatures/CN_sig/CN_sig_all.pdf', plot = p1, dpi = 600, width = 10, height = 7, units = 'in')

# 2. Pan-Cnacer signatures
PanCan_df <- as.data.frame(matrix(ncol = 6, nrow = 0))
colnames(PanCan_df) <- c('Group', 'enrich_PanCan', 'mean_PanCan_similarity', 'sd_PanCan_similarity', 'sample_size', 'percentage')
for (groupID in groups) {
  info_df <- sample_df %>%
    filter(Group == groupID) %>%
    filter(!is.na(enrich_PanCan)) %>%
    group_by(enrich_PanCan) %>%
    summarise(mean_PanCan_similarity = mean(PanCan_similarity),
              sd_PanCan_similarity = sd(PanCan_similarity),
              sample_size = n())
  total_n = sum(info_df$sample_size)
  for (i in 1:nrow(info_df)) {
    info_df[i, 'percentage'] = info_df[i, 'sample_size']/total_n*100
  }
  info_df['Group'] <- groupID
  info_df <- select(info_df, 'Group', 'enrich_PanCan':'percentage')
  print(info_df)
  PanCan_df <- rbind(PanCan_df, info_df)
}
## make the bar plot
PanCan_df$enrich_PanCan <- factor(PanCan_df$enrich_PanCan, levels = paste0('CX', c(1:17)))
p2 <- ggplot(data = PanCan_df, aes(x=Group, y=percentage, fill=enrich_PanCan)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='PanCan_signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic()
ggsave('results/signatures/PanCan_sig/PC_sig_all.pdf', plot = p2, dpi = 600, width = 10, height = 7, units = 'in')

######### Cellularity and ploidy differences ##########
# load the sample_signature data
# sample_df <- read.csv('results/signatures/All_sample_signature.csv', sep=',')
sample_df['Ploidy'] <- NA
sample_df['Cellularity'] <- NA
# load the solution table
solution_df <- read.csv('config/solution_sample.tsv', sep = '\t', row.names = 1)
# map the soluitons with samples
for (i in 1:nrow(sample_df)) {
  sampleID <- sample_df[i, 'Sample']
  sample_df[i, 'Ploidy'] <- solution_df[sampleID, 'Ploidy']
  sample_df[i, 'Cellularity'] <- solution_df[sampleID, 'Cellularity']
}

# calculate the mean cellularity and ploidy
## 1. CN signatures
for (groupID in groups) {
  info_df <- sample_df %>%
    filter(Group==groupID) %>%
    group_by(enrich_CN) %>%
    summarise(mean_CN_cellularity = mean(Cellularity),
              sd_CN_cellularity = sd(Cellularity),
              mean_CN_ploidy = mean(Ploidy),
              sd_CN_ploidy = sd(Ploidy))
  info_df['Group'] <- groupID
  info_df <- select(info_df, 'Group', 'enrich_CN', 'mean_CN_cellularity':'sd_CN_ploidy')
  cat(paste0(groupID, '\n'))
  print(info_df)
}
## 2. PanCan signatures
for (groupID in groups) {
  info_df <- sample_df %>%
    filter(Group==groupID) %>%
    group_by(enrich_PanCan) %>%
    summarise(mean_PanCan_cellularity = mean(Cellularity),
              sd_PanCan_cellularity = sd(Cellularity),
              mean_PanCan_ploidy = mean(Ploidy),
              sd_PanCan_ploidy = sd(Ploidy))
  info_df['Group'] <- groupID
  info_df <- select(info_df, 'Group', 'enrich_PanCan', 'mean_PanCan_cellularity':'sd_PanCan_ploidy')
  cat(paste0(groupID, '\n'))
  print(info_df)
}

######### Further analyzation #########
# samples and patients
patient_df <- sample_df %>% count(Patient) %>% filter(n>1)
for (patientID in patient_df$Patient) {
  a <- sample_df %>% filter(Patient==patientID) %>% count(Group)
  b <- sample_df %>% filter(Patient==patientID) %>% count(Type)
  cat(patientID)
  print(a)
  print(b)
}
## if the same patient has same signatures across samples
patient_df <- sample_df %>% count(Patient)
for (patientID in patient_df$Patient) {
  cat(paste0(patientID, '\n'))
  p <- sample_df %>% filter(Patient==patientID & Group!='ffpe')
  print('CN')
  CN_list <- unique(p$enrich_CN)
  print(CN_list)
  
  print('PanCan')
  PanCan_list <- unique(p$enrich_PanCan)
  print(PanCan_list)
  cat('\n')
}

## patients have both samples in ffTumor, ArchivalVS and MaNiLaVS
library(xlsx)
group_names <- c('ArchivalVS', 'ffTumor', 'MaNiLaVS')
patients <- unique(sample_df$Patient)
v_df <- as.data.frame(matrix(nrow = 0, ncol = 4))
for (patientID in patients) {
  n_ffTumor <- nrow(filter(sample_df, Patient == patientID & Group == 'ffTumor'))
  n_Archival <- nrow(filter(sample_df, Patient == patientID & Group == 'ArchivalVS'))
  n_MaNiLa <- nrow(filter(sample_df, Patient == patientID & Group == 'MaNiLaVS'))
  fill_content <- list(patientID, n_ffTumor, n_Archival, n_MaNiLa)
  v_df <- rbind(v_df, fill_content)
}
colnames(v_df) <- c('Patient', 'ffTumor', 'ArchivalVS', 'MaNiLaVS')

for (patientID in patients) {
  # tumor
  t_df <- sample_df %>% filter(Patient == patientID & Group == 'ffTumor')
  if (nrow(t_df)>0) {
    t_CN <- count(t_df, enrich_CN)
    t_CN <- t_CN[order(t_CN$n, decreasing = T),1][[1]]
    v_df[which(v_df$Patient==patientID),'ffTumor_CN'] <- t_CN
    
    t_PC <- count(t_df, enrich_PanCan)
    t_PC <- t_PC[order(t_PC$n, decreasing = T),1][[1]]
    v_df[which(v_df$Patient==patientID),'ffTumor_PC'] <- t_PC
  }
  else {
    v_df[which(v_df$Patient==patientID),'ffTumor_CN'] <- NA
    v_df[which(v_df$Patient==patientID),'ffTumor_PC'] <- NA
  }
  # archival
  a_df <- sample_df %>% filter(Patient == patientID & Group == 'ArchivalVS')
  if (nrow(a_df)>0) {
    a_CN <- count(a_df,enrich_CN)
    a_CN <- a_CN[order(a_CN$n, decreasing = T),1][[1]]
    v_df[which(v_df$Patient==patientID),'Archival_CN'] <- a_CN
    
    a_PC <- count(a_df,enrich_PanCan)
    a_PC <- a_PC[order(a_PC$n, decreasing = T),1][[1]]
    v_df[which(v_df$Patient==patientID),'Archival_PC'] <- a_PC
  }
  else {
    v_df[which(v_df$Patient==patientID),'Archival_CN'] <- NA
    v_df[which(v_df$Patient==patientID),'Archival_PC'] <- NA
  }
  # manila
  m_df <- sample_df %>% filter(Patient == patientID & Group == 'MaNiLaVS')
  if (nrow(m_df)>0) {
    m_CN <- count(m_df,enrich_CN)
    m_CN <- m_CN[order(m_CN$n, decreasing = T),1][[1]]
    v_df[which(v_df$Patient==patientID),'MaNiLa_CN'] <- m_CN
    
    m_PC <- count(m_df,enrich_PanCan)
    m_PC <- m_PC[order(m_PC$n, decreasing = T),1][[1]]
    v_df[which(v_df$Patient==patientID),'MaNiLa_PC'] <- m_PC
  }
  else {
    v_df[which(v_df$Patient==patientID),'MaNiLa_CN'] <- NA
    v_df[which(v_df$Patient==patientID),'MaNiLa_PC'] <- NA
  }
}

write.xlsx(v_df, file = 'results/signatures/reports/Reports.xlsx', sheetName = 'patients_3', append = TRUE, row.names = FALSE)

# whether there is correlation between CN_sig and PanCan_sig (chi-square test)
## in all samples
chi_df <- as.data.frame(matrix(nrow = 7, ncol = 5))
row.names(chi_df) <- c('s1', 's2', 's3', 's4', 's5', 's6', 's7')
colnames(chi_df) <- c('CX1', 'CX3', 'CX5', 'CX15', 'CX17')
for (CN_sig in row.names(chi_df)) {
  for (PC_sig in colnames(chi_df)) {
    sub_df <- sample_df %>% filter(enrich_CN == CN_sig & enrich_PanCan == PC_sig)
    content <- nrow(sub_df)
    chi_df[CN_sig, PC_sig] <- content
  }
}
chisq.test(chi_df)
# X-squared = 915.64, df = 24, p-value < 2.2e-16
write.xlsx(chi_df, file = 'results/signatures/reports/Reports.xlsx', sheetName = 'CN_PC_sig_chi', append = TRUE)


### test whether there are differences in ploidy and cellularity among signatures in each group ###
library(ggpubr)
# prepare the dataframe with log2 transformation on ploidy (better for plotting)
sample_df$enrich_CN <- factor(sample_df$enrich_CN, levels = c('s1', 's2', 's3', 's4', 's5', 's6', 's7'))
sample_df$enrich_PanCan <- factor(sample_df$enrich_PanCan, levels = c('CX1', 'CX3', 'CX5', 'CX15', 'CX17'))
lgt_df <- sample_df # %>% filter(Group=='ArchivalVS' | Group=='ffTumor' | Group=='MaNiLaVS')
lgt_df$Ploidy <- log2(lgt_df$Ploidy)
lgt_df_3 <- filter(lgt_df, Group=='ArchivalVS' | Group=='ffTumor' | Group=='MaNiLaVS')
## 1. CN signatures
# plotting
ploidy_plot <-  ggboxplot(lgt_df, x='Group', y='Ploidy', 
                          color = 'black', fill = 'enrich_CN', palette = 'Paired', 
                          xlab = 'Group', ylab = 'Ploidy (log2)') +
  theme_classic()
ggsave('results/signatures/reports/CN_ploidy_all.pdf', plot = ploidy_plot, dpi = 600, width = 10, height = 7, units = 'in')
cellularity_plot <-  ggboxplot(lgt_df, x='Group', y='Cellularity', 
                               color = 'black', fill = 'enrich_CN', palette = 'Paired', 
                               xlab = 'Group', ylab = 'Cellularity') +
  theme_classic()
ggsave('results/signatures/reports/CN_cellularity_all.pdf', plot = cellularity_plot, dpi = 600, width = 10, height = 7, units = 'in') 
## only 3 groups
ploidy_plot_3 <-  ggboxplot(lgt_df_3, x='Group', y='Ploidy', 
                            color = 'black', fill = 'enrich_CN', palette = 'Paired', 
                            xlab = 'Group', ylab = 'Ploidy (log2)') +
  theme_classic()
ggsave('results/signatures/reports/CN_ploidy_3.pdf', plot = ploidy_plot_3, dpi = 600, width = 10, height = 7, units = 'in')
cellularity_plot_3 <-  ggboxplot(lgt_df_3, x='Group', y='Cellularity', 
                                 color = 'black', fill = 'enrich_CN', palette = 'Paired', 
                                 xlab = 'Group', ylab = 'Cellularity') +
  theme_classic()
ggsave('results/signatures/reports/CN_cellularity_3.pdf', plot = cellularity_plot_3, dpi = 600, width = 10, height = 7, units = 'in')

# normal-distribution
library(car)
sub_df <- filter(sample_df, Group == 'ffTumor')
qqPlot(lm(sub_df$Ploidy~sub_df$enrich_CN), 
       data=sub_df,simulate=T,
       main="Q-Q Plot",labels=F) # ploidy does not follow the normal distribution
qqPlot(lm(sub_df$Cellularity~sub_df$enrich_CN), 
       data=sub_df,simulate=T,
       main="Q-Q Plot",labels=F) # cellularity does not follow the normal distribution, all of them have the same cellularity 0
# ploidy
kruskal.test(sub_df$Ploidy~sub_df$enrich_CN, data = sub_df)
# cellularity
kruskal.test(sub_df$Cellularity~sub_df$enrich_CN, data = sub_df)


# signature distribution
CN_df <- as.data.frame(matrix(ncol = 4, nrow = 0))
colnames(CN_df) <- c('enrich_CN', 'Group', 'sample_size', 'percentage')
CN_signatures <- c('s1', 's2', 's3', 's4', 's5', 's6', 's7')

for (CN_sig in CN_signatures) {
  info_df <- sample_df %>%
    # filter(Group=='ArchivalVS' | Group=='ffTumor' | Group=='MaNiLaVS') %>%
    filter(enrich_CN == CN_sig) %>%
    group_by(Group) %>% 
    summarise(sample_size = n())
  total_n = sum(info_df$sample_size)
  for (i in 1:nrow(info_df)) {
    info_df[i, 'percentage'] = info_df[i,'sample_size']/total_n*100
  }
  info_df['enrich_CN'] <- CN_sig
  info_df <- select(info_df, 'enrich_CN', 'Group', 'sample_size', 'percentage')
  cat(paste0(CN_sig, '\n'))
  print(info_df)
  CN_df <- rbind(CN_df, info_df)
}
## make the bar plot
p1 <- ggplot(data = CN_df, aes(x=enrich_CN, y=percentage, fill=Group)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='CN signatures', y='Percentage (%)', fill='Groups') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic()
ggsave('results/signatures/CN_sig_distribution.pdf', plot = p1, dpi = 600, width = 10, height = 7, units = 'in')
"
p1 <- ggplot(data = CN_df, aes(x=enrich_CN, y=percentage, fill=Group)) +
      geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
      labs(x='CN signatures', y='Percentage (%)', fill='Groups') +
      scale_fill_brewer(palette = 'Set2') +
      theme_classic()
ggsave('results/signatures/reports/CN_sig_distribution_3.pdf', plot = p1, dpi = 600, width = 10, height = 7, units = 'in')
"

## 2. Pan-Cancer signatures
# plotting
lgt_df <- filter(lgt_df, !is.na(enrich_PanCan))
lgt_df_3 <- filter(lgt_df_3, !is.na(enrich_PanCan))
ploidy_plot <-  ggboxplot(lgt_df, x='Group', y='Ploidy', 
                          color = 'black', fill = 'enrich_PanCan', palette = 'Paired', 
                          xlab = 'Group', ylab = 'Ploidy (log2)') +
  theme_classic()
ggsave('results/signatures/reports/PC_ploidy_all.pdf', plot = ploidy_plot, dpi = 600, width = 10, height = 7, units = 'in')
cellularity_plot <-  ggboxplot(lgt_df, x='Group', y='Cellularity', 
                               color = 'black', fill = 'enrich_PanCan', palette = 'Paired', 
                               xlab = 'Group', ylab = 'Cellularity') +
  theme_classic()
ggsave('results/signatures/reports/PC_cellularity_all.pdf', plot = cellularity_plot, dpi = 600, width = 10, height = 7, units = 'in') 
## only 3 groups
ploidy_plot_3 <-  ggboxplot(lgt_df_3, x='Group', y='Ploidy', 
                            color = 'black', fill = 'enrich_PanCan', palette = 'Paired', 
                            xlab = 'Group', ylab = 'Ploidy (log2)') +
  theme_classic()
ggsave('results/signatures/reports/PC_ploidy_3.pdf', plot = ploidy_plot_3, dpi = 600, width = 10, height = 7, units = 'in')
cellularity_plot_3 <-  ggboxplot(lgt_df_3, x='Group', y='Cellularity', 
                                 color = 'black', fill = 'enrich_PanCan', palette = 'Paired', 
                                 xlab = 'Group', ylab = 'Cellularity') +
  theme_classic()
ggsave('results/signatures/reports/PC_cellularity_3.pdf', plot = cellularity_plot_3, dpi = 600, width = 10, height = 7, units = 'in')

# normal-distribution
library(car)
sub_df <- filter(sample_df, Group == 'ffTumor') %>% filter(!is.na(enrich_PanCan))
qqPlot(lm(sub_df$Ploidy~sub_df$enrich_PanCan), 
       data=sub_df,simulate=T,
       main="Q-Q Plot",labels=F) # ploidy does not follow the normal distribution
qqPlot(lm(sub_df$Cellularity~sub_df$enrich_PanCan), 
       data=sub_df,simulate=T,
       main="Q-Q Plot",labels=F) # cellularity does not follow the normal distribution, all of them have the same cellularity 0
# ploidy
kruskal.test(sub_df$Ploidy~sub_df$enrich_PanCan, data = sub_df)
# cellularity
kruskal.test(sub_df$Cellularity~sub_df$enrich_PanCan, data = sub_df)


# signature distribution
PC_df <- as.data.frame(matrix(ncol = 4, nrow = 0))
colnames(PC_df) <- c('enrich_PanCan', 'Group', 'sample_size', 'percentage')
PanCan_signatures <- c('CX1', 'CX3', 'CX5', 'CX15', 'CX17')

for (PC_sig in PanCan_signatures) {
  info_df <- sample_df %>%
    #filter(Group=='ArchivalVS' | Group=='ffTumor' | Group=='MaNiLaVS') %>%
    filter(enrich_PanCan == PC_sig) %>%
    group_by(Group) %>% 
    summarise(sample_size = n())
  total_n = sum(info_df$sample_size)
  for (i in 1:nrow(info_df)) {
    info_df[i, 'percentage'] = info_df[i,'sample_size']/total_n*100
  }
  info_df['enrich_PanCan'] <- PC_sig
  info_df <- select(info_df, 'enrich_PanCan', 'Group', 'sample_size', 'percentage')
  cat(paste0(PC_sig, '\n'))
  print(info_df)
  PC_df <- rbind(PC_df, info_df)
}
PC_df$enrich_PanCan <- factor(PC_df$enrich_PanCan, levels = c('CX1', 'CX3', 'CX5', 'CX15', 'CX17'))
## make the bar plot
p2 <- ggplot(data = PC_df, aes(x=enrich_PanCan, y=percentage, fill=Group)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Pan-Cancer signatures', y='Percentage (%)', fill='Groups') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic()
ggsave('results/signatures/PC_sig_distribution.pdf', plot = p2, dpi = 600, width = 10, height = 7, units = 'in')
"
p2 <- ggplot(data = PC_df, aes(x=enrich_PanCan, y=percentage, fill=Group)) +
      geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
      labs(x='Pan-Cancer signatures', y='Percentage (%)', fill='Groups') +
      scale_fill_brewer(palette = 'Set2') +
      theme_classic()
ggsave('results/signatures/PC_sig_distribution_3.pdf', plot = p2, dpi = 600, width = 10, height = 7, units = 'in')
"