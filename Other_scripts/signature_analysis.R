# Title: signature_analysis.R
# Author: Guyuan TANG
# Date: 2024/4/20

# Description: the main script for analysis and plotting on signatures.


# set working directory
setwd('E:/1Lund_Lectures/BINP52_MasterProject/Manuscript/analysis/')

# load the required funcitons
source('function_analysis.R')

# load sample information
sample_df <- read.xlsx('sample_data.xlsx', sheetIndex = 1)

# group stratification
groups <- c('ArchivalVS','MaNiLaVS','ffTumor','ffTissue','Endome','FFPE','Blood','ffPlasma')
groups2 <- c('prediagnostic VS','diagnostic VS','ffTissue','Endome','FFPE','Blood','ffPlasma')

###### 1. CN signatures ######
# load the dataframe
CN_mat <- readRDS('Brenton/CN_sig.SSmatrix.rds')
## exclude the samples
for (sampleID in CN_mat$sample) {
  if (!(sampleID %in% sample_df$Sample)) {
    CN_mat <- filter(CN_mat, sample!=sampleID)
  }
}

# calculate the exposures
CN_exposure <- sig_exposure(SSmatrix = CN_mat,sig_type = 'CN',threshold=0.5)
saveRDS(CN_exposure, file = 'Brenton/CN_exposure.rds',compress = FALSE)
# output dataframe
CN_df <- select_top_sig(SSmatrix = CN_mat, number_of_top = 2, sig_type = 'CN')
# calculate the similarity differences between the top 2 signatures
CN_out <- delta_val_cal(input_df = CN_df, sig_type = 'CN', number_of_top = 2)
# add the sample information (patient, type, group, BH, BH_type)
CN_out <- add_df_info(SS_df = CN_out, sample_df = sample_df)
CN_out <- CN_out %>% select('sample', 'patient':'BH_type', 'enrich_CN_1':'delta_val')
# output the dataframe
write.xlsx(CN_out, file = 'All_sample_signature.xlsx', sheetName = 'CN', append = TRUE, row.names = FALSE)


# CN signature stats
CN_stat <- read.xlsx('All_sample_signature.xlsx',sheetName = 'CN')
## 1. distribution of the delta-value
CN_dist_delta_plot <- dis_plot_delta(CN_stat)
ggsave(filename = 'Brenton/Figures/delta_val_dist.pdf', plot = CN_dist_delta_plot, dpi = 600, width = 10, height = 7, units = 'in')
quantile(CN_stat$delta_val, 0.1) # 0.079 top 10%
# after cut-off on delta values
groups2 <- c('prediagnostic VS','diagnostic VS','ffTissue','Endome','FFPE','Blood','ffPlasma')
## 0.1 for ffTumor and FFPE
## 0.25 for all other groups
########## top 1 ##########
CN_stat_01 <- CN_stat[which(CN_stat$delta_val>=0.1),]
CN_stat_01 <- CN_stat_01 %>% filter(group=='ffTumor' | group=='FFPE')
CN_stat_25 <- CN_stat[which(CN_stat$delta_val>=0.25),]
CN_stat_25 <- CN_stat_25 %>% filter(group!='ffTumor' & group!='FFPE')
CN_stat_cut <- rbind(CN_stat_01, CN_stat_25)
#### top 2 enrich signature
CN_info_df_cut <- sig_stat_enrich(stat_df = CN_stat_cut, enrich_num = 1, sig_type = 'CN')
# plot the graph
CN_enrich_plot_1 <- ggplot(data = CN_info_df_cut, aes(x=BH, y=percentage, fill=enrich_CN_1)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='CN signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group_time))
CN_enrich_plot_1 <- CN_enrich_plot_1 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Brenton/Figures/CN_sig_BH_cut_1.pdf', plot = CN_enrich_plot_1, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
CN_info_df_cut <- as.data.frame(CN_info_df_cut)
write.xlsx(CN_info_df_cut, file = 'Brenton/CN_group.xlsx', sheetName = 'CN_cut_1', append = TRUE, row.names = FALSE)

########## top 2 ##########
CN_stat_01 <- CN_stat[which(CN_stat$delta_val<0.1),]
CN_stat_01 <- CN_stat_01 %>% filter(group=='ffTumor' | group=='FFPE')
CN_stat_25 <- CN_stat[which(CN_stat$delta_val<0.25),]
CN_stat_25 <- CN_stat_25 %>% filter(group!='ffTumor' & group!='FFPE')
CN_stat_cut <- rbind(CN_stat_01, CN_stat_25)
#### top 2 enrich signature
CN_info_df_cut <- sig_stat_enrich(stat_df = CN_stat_cut, enrich_num = 2, sig_type = 'CN')
# plot the graph
CN_enrich_plot_2 <- ggplot(data = CN_info_df_cut, aes(x=BH, y=percentage, fill=enrich_CN_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='CN signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group_time))
CN_enrich_plot_2 <- CN_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Brenton/Figures/CN_sig_BH_cut_2.pdf', plot = CN_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
CN_info_df_cut <- as.data.frame(CN_info_df_cut)
write.xlsx(CN_info_df_cut, file = 'Brenton/CN_group.xlsx', sheetName = 'CN_cut_2', append = TRUE, row.names = FALSE)

## 2. exposure distribution
### rank by group
CN_color_bar_group <- color_bar_prep(stat_df = CN_stat, exposure_df = CN_exposure, grouping_type = 'group_VS', sig_type = 'CN', sort = TRUE) + labs(fill = 'Types') 
CN_exposure_plot_group <- plot_exposure(stat_df = CN_stat, exposure_df = CN_exposure, grouping_type = 'group_VS', sig_type = 'CN', sort = TRUE)
CN_exposure_plot_group <- CN_exposure_plot_group %>% insert_bottom(CN_color_bar_group, height = 0.03)
ggsave(filename = 'Brenton/Figures/exposure_plot_groups_sort_v1.pdf', plot = CN_exposure_plot_group, dpi = 600, width = 40, height = 7, units = 'in')
### rank by BH
CN_normal_stat <- CN_stat %>% filter(group_VS == 'Normal tissues')
normal_sample <- CN_normal_stat$sample
CN_normal_exposure <- CN_exposure %>% filter(sample %in% normal_sample)
CN_color_bar_BH <- color_bar_prep(stat_df = CN_normal_stat, exposure_df = CN_normal_exposure, grouping_type = 'BH_group', sig_type = 'CN', sort = TRUE)
CN_exposure_plot_BH <- plot_exposure(stat_df = CN_normal_stat, exposure_df = CN_normal_exposure, grouping_type = 'BH_group', sig_type = 'CN', sort = TRUE)
CN_exposure_plot_BH <- CN_exposure_plot_BH %>% insert_bottom(CN_color_bar_BH, height = 0.03)
ggsave(filename = 'Brenton/Figures/exposure_plot_BH_normal.pdf', plot = CN_exposure_plot_BH, dpi = 600, width = 30, height = 7, units = 'in')

## 3. CN signatures enrichment
### in this part, combine ffTumor and ffTissue into ffTissue
groups2 <- c('prediagnostic VS','diagnostic VS','ffTissue','Endome','FFPE','Blood','ffPlasma')

#### top 1 enrich signature
CN_info_df <- sig_stat_enrich(stat_df = CN_stat, enrich_num = 1, sig_type = 'CN')
# plot the graph
CN_enrich_plot_1 <- ggplot(data = CN_info_df, aes(x=BH, y=percentage, fill=enrich_CN_1)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='CN signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group_time))
CN_enrich_plot_1 <- CN_enrich_plot_1 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Brenton/Figures/CN_sig_BH_all.pdf', plot = CN_enrich_plot_1, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
CN_info_df <- as.data.frame(CN_info_df)
write.xlsx(CN_info_df, file = 'Brenton/CN_group.xlsx', sheetName = 'CN_1', append = TRUE, row.names = FALSE)

#### top 2 enrich signature
CN_info_df <- sig_stat_enrich(stat_df = CN_stat, enrich_num = 2, sig_type = 'CN')
# plot the graph
CN_enrich_plot_2 <- ggplot(data = CN_info_df, aes(x=BH, y=percentage, fill=enrich_CN_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='CN signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group_time))
CN_enrich_plot_2 <- CN_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Brenton/Figures/CN_sig_BH_all_2.pdf', plot = CN_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
CN_info_df <- as.data.frame(CN_info_df)
write.xlsx(CN_info_df, file = 'Brenton/CN_group.xlsx', sheetName = 'CN_2', append = TRUE, row.names = FALSE)

## 4. VS samples time point exposure
VS_sheet <- read.xlsx(file = 'MaNiLa_All_samplesheet_240311_groups.xlsx', sheetName = 'VS_samples')
VS_stat <- CN_stat[which(CN_stat$group_VS=='Cervical'),]
for (sampleID in VS_stat$sample) {
  VS_stat[which(VS_stat$sample==sampleID), 'Diagnostic_time'] <- VS_sheet[which(VS_sheet$Library==sampleID), 'Diagnostic_time']
}

### extract the similarity matrix
indices <- rownames(VS_stat)
VS_mat <- CN_mat[indices,]
VS_exposure <- sig_exposure(SSmatrix = VS_mat,sig_type = 'CN',threshold = 0.5)
VS_color_bar_BH <- color_bar_prep(stat_df = VS_stat, exposure_df = VS_exposure, grouping_type = 'BH_group', sig_type = 'CN')
VS_exposure_plot_BH <- plot_exposure(stat_df = VS_stat, exposure_df = VS_exposure, grouping_type = 'BH_group', sig_type = 'CN')
VS_exposure_plot_BH <- VS_exposure_plot_BH %>% insert_bottom(VS_color_bar_BH, height = 0.03)
ggsave('Brenton/Figures/CN_VS_exposure_v1.pdf', plot = VS_exposure_plot_BH, dpi = 600, width = 10, height = 7, units = 'in')

### for Benign VS samples
Benign_VS_stat <- VS_stat[which(VS_stat$BH=='Benign'),]
unique(Benign_VS_stat$Diagnostic_time) # only "diagnostic"
indices <- rownames(Benign_VS_stat)
Benign_VS_mat <- VS_mat[indices,]
#### exposure
Benign_VS_exposure <- sig_exposure(SSmatrix = Benign_VS_mat,sig_type = 'CN',threshold = 0.5)
#### sort by s1
Benign_VS_exposure <- Benign_VS_exposure[order(-Benign_VS_exposure$s1),]
Benign_VS_exposure <- select(Benign_VS_exposure,!SS_sum)
#### plotting
Benign_VS_plot <- VS_exposure_plot(VS_exposure_df = Benign_VS_exposure, sig_type = 'CN', VS_stat_df = Benign_VS_stat) + labs(title = 'diagnostic') + theme(plot.title = element_text(size = 15, hjust = 0.5))
ggsave('Brenton/Figures/Benign_VS_exposure.pdf', plot = Benign_VS_plot, dpi = 600, width = 7, height = 5, units = 'in')

### for RRSO samples
RRSO_VS_stat <- VS_stat[which(VS_stat$BH=='RRSO'),]
unique(RRSO_VS_stat$Diagnostic_time) # "prediagnostic_6+", "diagnostic"
indices <- rownames(RRSO_VS_stat)
RRSO_VS_mat <- VS_mat[indices,]
#### exposure
RRSO_VS_exposure_in <- sig_exposure(SSmatrix = RRSO_VS_mat,sig_type = 'CN',threshold = 0.5)
#### sort by s1
# RRSO_VS_exposure <- RRSO_VS_exposure[order(-RRSO_VS_exposure$s1),]
RRSO_VS_exposure <- as.data.frame(matrix(nrow = 0, ncol = 9))
diag_group_num <- c()
for (diagnostic_group in c("prediagnostic_6+", "diagnostic")) {
  indices_val <- which(RRSO_VS_stat$Diagnostic_time==diagnostic_group)
  sort_df <- RRSO_VS_exposure_in[indices_val,]
  sort_df <- sort_df[order(-sort_df$s1),]
  RRSO_VS_exposure <- rbind(RRSO_VS_exposure, sort_df)
  diag_group_num <- c(diag_group_num, nrow(sort_df))
}
RRSO_VS_exposure <- select(RRSO_VS_exposure,!SS_sum)
# diag_group_num == c(13, 17)
pre_6_RRSO_exposure <- RRSO_VS_exposure[1:13,]
diag_RRSO_exposure <- RRSO_VS_exposure[14:30,]
#### plotting
pre_6_RRSO_exposure_plot <- VS_exposure_plot(VS_exposure_df = pre_6_RRSO_exposure, sig_type = 'CN', VS_stat_df = RRSO_VS_stat) + labs(title = 'prediagnostic > 6m') + theme(plot.title = element_text(size = 15, hjust = 0.5))
diag_RRSO_exposure_plot <- VS_exposure_plot(VS_exposure_df = diag_RRSO_exposure, sig_type = 'CN', VS_stat_df = RRSO_VS_stat)  + labs(title = 'diagnostic') + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(size = 15, hjust = 0.5))
RRSO_VS_plot <- ggarrange(pre_6_RRSO_exposure_plot, diag_RRSO_exposure_plot, ncol = 2, common.legend = 1, legend.grob = get_legend(pre_6_RRSO_exposure_plot), legend = 'right')
ggsave('Brenton/Figures/RRSO_VS_exposure.pdf', plot = RRSO_VS_plot, dpi = 600, width = 7, height = 5, units = 'in')

### for HGSC samples
HGSC_VS_stat <- VS_stat[which(VS_stat$BH=='HGSC'),]
unique(HGSC_VS_stat$Diagnostic_time) # "prediagnostic_0-6", "prediagnostic_6+", "diagnostic"
indices <- rownames(HGSC_VS_stat)
HGSC_VS_mat <- VS_mat[indices,]
#### exposure
HGSC_VS_exposure_in <- sig_exposure(SSmatrix = HGSC_VS_mat,sig_type = 'CN',threshold = 0.5)
#### sort by s1
HGSC_VS_exposure <- as.data.frame(matrix(nrow = 0, ncol = 9))
diag_group_num <- c()
for (diagnostic_group in c("prediagnostic_6+", "prediagnostic_0-6", "diagnostic")) {
  indices_val <- which(HGSC_VS_stat$Diagnostic_time==diagnostic_group)
  sort_df <- HGSC_VS_exposure_in[indices_val,]
  sort_df <- sort_df[order(-sort_df$s1),]
  HGSC_VS_exposure <- rbind(HGSC_VS_exposure, sort_df)
  diag_group_num <- c(diag_group_num, nrow(sort_df))
}
HGSC_VS_exposure <- select(HGSC_VS_exposure,!SS_sum)
# diag_group_num == c(62,30,36)
pre_6_HGSC_exposure <- HGSC_VS_exposure[1:62,]
pre_0_HGSC_exposure <- HGSC_VS_exposure[63:92,]
diag_HGSC_exposure <- HGSC_VS_exposure[93:128,]
#### plotting
#1
pre_6_HGSC_exposure_plot <- VS_exposure_plot(VS_exposure_df = pre_6_HGSC_exposure, sig_type = 'CN', VS_stat_df = HGSC_VS_stat) + labs(title = 'prediagnostic > 6m') + theme(plot.title = element_text(size = 15, hjust = 0.5))
#2
pre_0_HGSC_exposure_plot <- VS_exposure_plot(VS_exposure_df = pre_0_HGSC_exposure, sig_type = 'CN', VS_stat_df = HGSC_VS_stat)  + labs(title = 'prediagnostic 0-6m') + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(size = 15, hjust = 0.5))
#3
diag_HGSC_exposure_plot <- VS_exposure_plot(VS_exposure_df = diag_HGSC_exposure, sig_type = 'CN', VS_stat_df = HGSC_VS_stat)  + labs(title = 'diagnostic') + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(size = 15, hjust = 0.5))
#combine
HGSC_VS_plot <- ggarrange(pre_6_HGSC_exposure_plot, pre_0_HGSC_exposure_plot, diag_HGSC_exposure_plot, ncol = 3, common.legend = 1, legend.grob = get_legend(pre_0_HGSC_exposure_plot), legend = 'right')
ggsave('Brenton/Figures/HGSC_VS_exposure.pdf', plot = HGSC_VS_plot, dpi = 600, width = 10, height = 5, units = 'in')

## 6. reference exposure plots
#### HGSC ffTumor 
HGSC_tumor_indices <- which(CN_stat$group=='ffTumor' & CN_stat$BH=='HGSC')
HGSC_tumor_mat <- CN_mat[HGSC_tumor_indices,]
HGSC_tumor_stat <- CN_stat[HGSC_tumor_indices,]
HGSC_tumor_exposure <- sig_exposure(SSmatrix = HGSC_tumor_mat, sig_type = 'CN',threshold = 0.5)
#### sort by s1
HGSC_tumor_exposure <- HGSC_tumor_exposure[order(-HGSC_tumor_exposure$s1),]
HGSC_tumor_exposure <- select(HGSC_tumor_exposure,!SS_sum)
#### plotting
HGSC_tumor_exposure_plot <- ref_exposure_plot(ref_exposure_df = HGSC_tumor_exposure, sig_type = 'CN', ref_stat_df = HGSC_tumor_stat) + labs(title = 'HGSC ffTumor') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('Brenton/cosine_similarity/Figures/ref_HGSC_exposure.pdf', plot = HGSC_tumor_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### Benign ffTissue
Benign_tissue_indices <- which(CN_stat$group=='ffTissue' & CN_stat$BH=='Benign')
Benign_tissue_mat <- CN_mat[Benign_tissue_indices,]
Benign_tissue_stat <- CN_stat[Benign_tissue_indices,]
Benign_tissue_exposure <- sig_exposure(SSmatrix = Benign_tissue_mat, sig_type = 'CN',threshold = 0.5)
#### sort by s1
Benign_tissue_exposure <- Benign_tissue_exposure[order(-Benign_tissue_exposure$s1),]
Benign_tissue_exposure <- select(Benign_tissue_exposure, !SS_sum)
#### plotting
Benign_tissue_exposure_plot <- ref_exposure_plot(ref_exposure_df = Benign_tissue_exposure, sig_type = 'CN', ref_stat_df = Benign_tissue_stat) + labs(title = 'Benign ffTissue') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('Brenton/cosine_similarity/Figures/ref_Benign_exposure.pdf', plot = Benign_tissue_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### RRSO ffTissue+FFPE
RRSO_ref_indices <- which(CN_stat$BH=='RRSO' & (CN_stat$group=='ffTissue' | CN_stat$group=='FFPE'))
RRSO_ref_mat <- CN_mat[RRSO_ref_indices,]
RRSO_ref_stat <- CN_stat[RRSO_ref_indices,]
RRSO_ref_exposure <- sig_exposure(SSmatrix = RRSO_ref_mat, sig_type = 'CN',threshold = 0.5)
#### sort by s1
RRSO_ref_exposure <- RRSO_ref_exposure[order(-RRSO_ref_exposure$s1),]
RRSO_ref_exposure <- select(RRSO_ref_exposure, !SS_sum)
#### plotting
RRSO_ref_exposure_plot <- ref_exposure_plot(ref_exposure_df = RRSO_ref_exposure, sig_type = 'CN', ref_stat_df = RRSO_ref_stat) + labs(title = 'RRSO ffTissue&FFPE') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('Brenton/cosine_similarity/Figures/ref_RRSO_exposure.pdf', plot = RRSO_ref_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### combined
RRSO_ref_exposure_plot_1 <- RRSO_ref_exposure_plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
HGSC_tumor_exposure_plot_1 <- HGSC_tumor_exposure_plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
CN_ref_exposure_plot <- ggarrange(Benign_tissue_exposure_plot, RRSO_ref_exposure_plot_1, HGSC_tumor_exposure_plot_1, ncol = 3, common.legend = 1, legend.grob = get_legend(Benign_tissue_exposure_plot), legend = 'right')
ggsave('Brenton/cosine_similarity/Figures/ref_exposure_all.pdf', plot = CN_ref_exposure_plot, dpi = 600, width = 10, height = 5, units = 'in')


###### 2. Pan-Cancer signatures ######
# load the dataframe
PanCan_mat <- readRDS('Pan-Cancer/PanCan_sig.SSmatrix.rds')
## exclude the samples
for (sampleID in PanCan_mat$sample) {
  if (!(sampleID %in% sample_df$Sample)) {
    PanCan_mat <- filter(PanCan_mat, sample!=sampleID)
  }
}

# calculate the exposures
PanCan_exposure <- sig_exposure(SSmatrix = PanCan_mat,sig_type = 'PanCan',threshold = 0.5)
saveRDS(PanCan_exposure, file = 'Pan-Cancer/CIN_exposure.rds',compress = FALSE)

# output dataframe
PanCan_df <- select_top_sig(SSmatrix = PanCan_mat, number_of_top = 2, sig_type = 'PanCan')
# calculate the similarity differences between the top 2 signatures
PanCan_out <- delta_val_cal(input_df = PanCan_df, sig_type = 'PanCan', number_of_top = 2)
PanCan_out <- add_df_info(SS_df = PanCan_out, sample_df = sample_df)
PanCan_out <- PanCan_out %>% select('sample', 'patient':'BH_type', 'enrich_PanCan_1':'delta_val')
# output the dataframe
write.xlsx(PanCan_out, file = 'All_sample_signature.xlsx', sheetName = 'Pan_Cancer', append = TRUE, row.names = FALSE)


# PanCan signature stats
PanCan_stat <- read.xlsx('All_sample_signature.xlsx',sheetName = 'Pan_Cancer')
## 1. distribution of the delta-value
PanCan_dist_delta_plot <- dis_plot_delta(PanCan_stat)
ggsave(filename = 'Pan-Cancer/Figures/delta_val_dist.pdf', plot = PanCan_dist_delta_plot, dpi = 600, width = 10, height = 7, units = 'in')
quantile(PanCan_stat$delta_val, 0.1) # 0.016 top 10%
# delta value cut-off at 0.25
PanCan_stat_25 <- PanCan_stat[which(PanCan_stat$delta_val<0.25),]
#### top 2 enrich signature
PanCan_info_df_25 <- sig_stat_enrich(stat_df = PanCan_stat_25, enrich_num = 2, sig_type = 'PanCan')
PanCan_info_df_25$enrich_PanCan_2 <- factor(PanCan_info_df_25$enrich_PanCan_2, levels = paste0('CX', c(1:17)))
# save the info dataframe
PanCan_info_df_25 <- as.data.frame(PanCan_info_df_25)
write.xlsx(PanCan_info_df_25, file = 'Pan-Cancer/Pan-Cancer_group.xlsx', sheetName = 'PanCan_025_2', row.names = FALSE, append = TRUE)

# delta value cut-off at 0.1
PanCan_stat_01 <- PanCan_stat[which(PanCan_stat$delta_val<0.1),]
#### top 2 enrich signature
PanCan_info_df_01 <- sig_stat_enrich(stat_df = PanCan_stat_01, enrich_num = 2, sig_type = 'PanCan')
PanCan_info_df_01$enrich_PanCan_2 <- factor(PanCan_info_df_01$enrich_PanCan_2, levels = paste0('CX', c(1:17)))
# save the info dataframe
PanCan_info_df_01 <- as.data.frame(PanCan_info_df_01)
write.xlsx(PanCan_info_df_01, file = 'Pan-Cancer/Pan-Cancer_group.xlsx', sheetName = 'PanCan_01_2', row.names = FALSE, append = TRUE)

# after adjusting the delta-value cutoff
## 0.25 for ffTissue, Endome, and Blood
## 0.1 for all other groups
########## top 1 ##########
PanCan_stat_01 <- PanCan_stat[which(PanCan_stat$delta_val>=0.1),]
PanCan_stat_01 <- PanCan_stat_01 %>% filter(group_VS!='ffTissue' & group_VS!='Endome' & group_VS!='Blood')
PanCan_stat_25 <- PanCan_stat[which(PanCan_stat$delta_val>=0.25),]
PanCan_stat_25 <- PanCan_stat_25 %>% filter(group_VS=='ffTissue' | group_VS=='Endome' | group_VS=='Blood')
PanCan_stat_cut <- rbind(PanCan_stat_01, PanCan_stat_25)
#### top 1 enrich signature
PanCan_info_df_cut <- sig_stat_enrich(stat_df = PanCan_stat_cut, enrich_num = 1, sig_type = 'PanCan')
PanCan_info_df_cut$enrich_PanCan_1 <- factor(PanCan_info_df_cut$enrich_PanCan_1, levels = paste0('CX',c(1:17))) 
# plot the graph
PanCan_enrich_plot_1 <- ggplot(data = PanCan_info_df_cut, aes(x=BH, y=percentage, fill=enrich_PanCan_1)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='Pan-Cancer signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group_time))
PanCan_enrich_plot_1 <- PanCan_enrich_plot_1 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Pan-Cancer/Figures/PanCan_sig_BH_cut_1.pdf', plot = PanCan_enrich_plot_1, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
PanCan_info_df_cut <- as.data.frame(PanCan_info_df_cut)
write.xlsx(PanCan_info_df_cut, file = 'Pan-Cancer/Pan-Cancer_group.xlsx', sheetName = 'PanCan_cut_1', row.names = FALSE, append = TRUE)

########## top 2 ##########
PanCan_stat_01 <- PanCan_stat[which(PanCan_stat$delta_val<0.1),]
PanCan_stat_01 <- PanCan_stat_01 %>% filter(group_VS!='ffTissue' & group_VS!='Endome' & group_VS!='Blood')
PanCan_stat_25 <- PanCan_stat[which(PanCan_stat$delta_val<0.25),]
PanCan_stat_25 <- PanCan_stat_25 %>% filter(group_VS=='ffTissue' | group_VS=='Endome' | group_VS=='Blood')
PanCan_stat_cut <- rbind(PanCan_stat_01, PanCan_stat_25)
#### top 2 enrich signature
PanCan_info_df_cut <- sig_stat_enrich(stat_df = PanCan_stat_cut, enrich_num = 2, sig_type = 'PanCan')
PanCan_info_df_cut$enrich_PanCan_2 <- factor(PanCan_info_df_cut$enrich_PanCan_2, levels = paste0('CX',c(1:17)))
# plot the graph
PanCan_enrich_plot_2 <- ggplot(data = PanCan_info_df_cut, aes(x=BH, y=percentage, fill=enrich_PanCan_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='Pan-Cancer signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group_time))
PanCan_enrich_plot_2 <- PanCan_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Pan-Cancer/Figures/PanCan_sig_BH_cut_2.pdf', plot = PanCan_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
PanCan_info_df_cut <- as.data.frame(PanCan_info_df_cut)
write.xlsx(PanCan_info_df_cut, file = 'Pan-Cancer/Pan-Cancer_group.xlsx', sheetName = 'PanCan_cut_2', row.names = FALSE, append = TRUE)


## 2. exposure distribution
### rank by group
PanCan_color_bar_group <- color_bar_prep(stat_df = PanCan_stat, exposure_df = PanCan_exposure, grouping_type = 'group_VS', sig_type = 'Pan-Cancer', sort = TRUE) + labs(fill = 'Types')
PanCan_exposure_plot_group <- plot_exposure(stat_df = PanCan_stat, exposure_df = PanCan_exposure, grouping_type = 'group_VS', sig_type = 'Pan-Cancer', sort = TRUE)
PanCan_exposure_plot_group <- PanCan_exposure_plot_group %>% insert_bottom(PanCan_color_bar_group, height = 0.03)
ggsave(filename = 'Pan-Cancer/Figures/exposure_plot_groups_sort_v1.pdf', plot = PanCan_exposure_plot_group, dpi = 600, width = 40, height = 7, units = 'in')
### rank by BH
PanCan_normal_stat <- PanCan_stat %>% filter(group_VS == 'Normal tissues')
normal_sample <- PanCan_normal_stat$sample
PanCan_normal_exposure <- PanCan_exposure %>% filter(sample %in% normal_sample)
PanCan_color_bar_BH <- color_bar_prep(stat_df = PanCan_normal_stat, exposure_df = PanCan_normal_exposure, grouping_type = 'BH_group', sig_type = 'Pan-Cancer', sort = TRUE)
PanCan_exposure_plot_BH <- plot_exposure(stat_df = PanCan_normal_stat, exposure_df = PanCan_normal_exposure, grouping_type = 'BH_group', sig_type = 'Pan-Cancer', sort = TRUE)
PanCan_exposure_plot_BH <- PanCan_exposure_plot_BH %>% insert_bottom(PanCan_color_bar_BH, height = 0.03)
ggsave(filename = 'Pan-Cancer/Figures/exposure_plot_BH_normal.pdf', plot = PanCan_exposure_plot_BH, dpi = 600, width = 30, height = 7, units = 'in')


## 3. PanCan signatures enrichment
#### top 1 enriched signature
PanCan_info_df <- sig_stat_enrich(stat_df = PanCan_stat, enrich_num = 1, sig_type = 'PanCan')
PanCan_info_df$enrich_PanCan_1 <- factor(PanCan_info_df$enrich_PanCan_1, levels = paste0('CX',c(1:17)))
# plot the graph
PanCan_enrich_plot_1 <- ggplot(data = PanCan_info_df, aes(x=BH, y=percentage, fill=enrich_PanCan_1)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='Pan-Cancer signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group_time))
PanCan_enrich_plot_1 <- PanCan_enrich_plot_1 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Pan-Cancer/Figures/PanCan_sig_BH_all.pdf', plot = PanCan_enrich_plot_1, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
PanCan_info_df <- as.data.frame(PanCan_info_df)
write.xlsx(PanCan_info_df, file = 'Pan-Cancer/Pan-Cancer_group.xlsx', sheetName = 'PanCan_1', row.names = FALSE, append = TRUE)

#### top 2 enrich signature
PanCan_info_df <- sig_stat_enrich(stat_df = PanCan_stat, enrich_num = 2, sig_type = 'PanCan')
PanCan_info_df$enrich_PanCan_2 <- factor(PanCan_info_df$enrich_PanCan_2, levels = paste0('CX',c(1:17)))
# plot the graph
PanCan_enrich_plot_2 <- ggplot(data = PanCan_info_df, aes(x=BH, y=percentage, fill=enrich_PanCan_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='Pan-Cancer signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group_time))
PanCan_enrich_plot_2 <- PanCan_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Pan-Cancer/Figures/PanCan_sig_BH_all_2.pdf', plot = PanCan_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
PanCan_info_df <- as.data.frame(PanCan_info_df)
write.xlsx(PanCan_info_df, file = 'Pan-Cancer/Pan-Cancer_group.xlsx', sheetName = 'PanCan_2', row.names = FALSE, append = TRUE)


## 4. VS samples time point exposure
VS_sheet <- read.xlsx(file = 'MaNiLa_All_samplesheet_240311_groups.xlsx', sheetName = 'VS_samples')
VS_stat <- PanCan_stat[which(PanCan_stat$group=='ArchivalVS' | PanCan_stat$group=='MaNiLaVS'),]
for (sampleID in VS_stat$sample) {
  VS_stat[which(VS_stat$sample==sampleID), 'Diagnostic_time'] <- VS_sheet[which(VS_sheet$Library==sampleID), 'Diagnostic_time']
}

### extract the similarity matrix
indices <- rownames(VS_stat)
VS_mat <- PanCan_mat[indices,]
VS_exposure <- sig_exposure(SSmatrix = VS_mat,sig_type = 'PanCan', threshold = 0.5)
VS_color_bar_BH <- color_bar_prep(stat_df = VS_stat, exposure_df = VS_exposure, grouping_type = 'BH_group', sig_type = 'Pan-Cancer')
VS_exposure_plot_BH <- plot_exposure(stat_df = VS_stat, exposure_df = VS_exposure, grouping_type = 'BH_group', sig_type = 'Pan-Cancer')
VS_exposure_plot_BH <- VS_exposure_plot_BH %>% insert_bottom(VS_color_bar_BH, height = 0.03)
ggsave('Pan-Cancer/cosine_similarity/Figures/VS_exposure.pdf', plot = VS_exposure_plot_BH, dpi = 600, width = 10, height = 7, units = 'in')

### for Benign VS samples
Benign_VS_stat <- VS_stat[which(VS_stat$BH=='Benign'),]
unique(Benign_VS_stat$Diagnostic_time) # only "diagnostic"
indices <- rownames(Benign_VS_stat)
Benign_VS_mat <- VS_mat[indices,]
#### exposure
Benign_VS_exposure <- sig_exposure(SSmatrix = Benign_VS_mat,sig_type = 'PanCan',threshold = 0.5)
#### sort by s1
Benign_VS_exposure <- Benign_VS_exposure[order(-Benign_VS_exposure$CX1),]
Benign_VS_exposure <- select(Benign_VS_exposure,!SS_sum)
#### plotting
Benign_VS_plot <- VS_exposure_plot(VS_exposure_df = Benign_VS_exposure, sig_type = 'Pan-Cancer', VS_stat_df = Benign_VS_stat) + labs(title = 'diagnostic') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('Pan-Cancer/cosine_similarity/Figures/Benign_VS_exposure.pdf', plot = Benign_VS_plot, dpi = 600, width = 7, height = 5, units = 'in')

### for RRSO samples
RRSO_VS_stat <- VS_stat[which(VS_stat$BH=='RRSO'),]
unique(RRSO_VS_stat$Diagnostic_time) # "diagnostic"
indices <- rownames(RRSO_VS_stat)
RRSO_VS_mat <- VS_mat[indices,]
#### exposure
RRSO_VS_exposure <- sig_exposure(SSmatrix = RRSO_VS_mat,sig_type = 'PanCan',threshold = 0.5)
#### sort by CX1
RRSO_VS_exposure <- RRSO_VS_exposure[order(-RRSO_VS_exposure$CX1),]
RRSO_VS_exposure <- select(RRSO_VS_exposure,!SS_sum)
#### plotting
RRSO_VS_plot <- VS_exposure_plot(VS_exposure_df = RRSO_VS_exposure, sig_type = 'Pan-Cancer', VS_stat_df = RRSO_VS_stat) + labs(title = 'diagnostic') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('Pan-Cancer/cosine_similarity/Figures/RRSO_VS_exposure.pdf', plot = RRSO_VS_plot, dpi = 600, width = 7, height = 5, units = 'in')

### for HGSC samples
HGSC_VS_stat <- VS_stat[which(VS_stat$BH=='HGSC'),]
unique(HGSC_VS_stat$Diagnostic_time) # "prediagnostic_0-6", "prediagnostic_6+", "diagnostic"
indices <- rownames(HGSC_VS_stat)
HGSC_VS_mat <- VS_mat[indices,]
#### exposure
HGSC_VS_exposure_in <- sig_exposure(SSmatrix = HGSC_VS_mat,sig_type = 'PanCan',threshold = 0.5)
#### sort by CX1
HGSC_VS_exposure <- as.data.frame(matrix(nrow = 0, ncol = 19))
diag_group_num <- c()
for (diagnostic_group in c("prediagnostic_6+", "prediagnostic_0-6", "diagnostic")) {
  indices_val <- which(HGSC_VS_stat$Diagnostic_time==diagnostic_group)
  sort_df <- HGSC_VS_exposure_in[indices_val,]
  sort_df <- sort_df[order(-sort_df$CX1),]
  HGSC_VS_exposure <- rbind(HGSC_VS_exposure, sort_df)
  diag_group_num <- c(diag_group_num, nrow(sort_df))
}
HGSC_VS_exposure <- select(HGSC_VS_exposure,!SS_sum)
# diag_group_num == c(31,24,9)
pre_6_HGSC_exposure <- HGSC_VS_exposure[1:31,]
pre_0_HGSC_exposure <- HGSC_VS_exposure[32:55,]
diag_HGSC_exposure <- HGSC_VS_exposure[56:64,]
#### plotting
#1
pre_6_HGSC_exposure_plot <- VS_exposure_plot(VS_exposure_df = pre_6_HGSC_exposure, sig_type = 'Pan-Cancer', VS_stat_df = HGSC_VS_stat) + labs(title = 'prediagnostic > 6m') + theme(plot.title = element_text(size = 10, hjust = 0.5))
#2
pre_0_HGSC_exposure_plot <- VS_exposure_plot(VS_exposure_df = pre_0_HGSC_exposure, sig_type = 'Pan-Cancer', VS_stat_df = HGSC_VS_stat)  + labs(title = 'prediagnostic 0-6m') + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(size = 10, hjust = 0.5))
#3
diag_HGSC_exposure_plot <- VS_exposure_plot(VS_exposure_df = diag_HGSC_exposure, sig_type = 'Pan-Cancer', VS_stat_df = HGSC_VS_stat)  + labs(title = 'diagnostic') + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(size = 10, hjust = 0.5))
#combine
HGSC_VS_plot <- ggarrange(pre_6_HGSC_exposure_plot, pre_0_HGSC_exposure_plot, diag_HGSC_exposure_plot, ncol = 3, common.legend = 1, legend.grob = get_legend(pre_0_HGSC_exposure_plot), legend = 'right')
ggsave('Pan-Cancer/cosine_similarity/Figures/HGSC_VS_exposure.pdf', plot = HGSC_VS_plot, dpi = 600, width = 10, height = 5, units = 'in')


## 5. reference exposure plots
#### HGSC ffTumor 
HGSC_tumor_indices <- which(PanCan_stat$group=='ffTumor' & PanCan_stat$BH=='HGSC')
HGSC_tumor_mat <- PanCan_mat[HGSC_tumor_indices,]
HGSC_tumor_stat <- PanCan_stat[HGSC_tumor_indices,]
HGSC_tumor_exposure <- sig_exposure(SSmatrix = HGSC_tumor_mat, sig_type = 'PanCan',threshold = 0.5)
#### sort by CX1
HGSC_tumor_exposure <- HGSC_tumor_exposure[order(-HGSC_tumor_exposure$CX1),]
HGSC_tumor_exposure <- select(HGSC_tumor_exposure,!SS_sum)
#### plotting
HGSC_tumor_exposure_plot <- ref_exposure_plot(ref_exposure_df = HGSC_tumor_exposure, sig_type = 'Pan-Cancer', ref_stat_df = HGSC_tumor_stat) + labs(title = 'HGSC ffTumor') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('Pan-Cancer/cosine_similarity/Figures/ref_HGSC_exposure.pdf', plot = HGSC_tumor_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### Benign ffTissue
Benign_tissue_indices <- which(PanCan_stat$group=='ffTissue' & PanCan_stat$BH=='Benign')
Benign_tissue_mat <- PanCan_mat[Benign_tissue_indices,]
Benign_tissue_stat <- PanCan_stat[Benign_tissue_indices,]
Benign_tissue_exposure <- sig_exposure(SSmatrix = Benign_tissue_mat, sig_type = 'PanCan',threshold = 0.5)
#### sort by CX1
Benign_tissue_exposure <- Benign_tissue_exposure[order(-Benign_tissue_exposure$CX1),]
Benign_tissue_exposure <- select(Benign_tissue_exposure, !SS_sum)
#### plotting
Benign_tissue_exposure_plot <- ref_exposure_plot(ref_exposure_df = Benign_tissue_exposure, sig_type = 'Pan-Cancer', ref_stat_df = Benign_tissue_stat) + labs(title = 'Benign ffTissue') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('Pan-Cancer/cosine_similarity/Figures/ref_Benign_exposure.pdf', plot = Benign_tissue_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### RRSO ffTissue+FFPE
RRSO_ref_indices <- which(PanCan_stat$BH=='RRSO' & (PanCan_stat$group=='ffTissue' | PanCan_stat$group=='FFPE'))
RRSO_ref_mat <- PanCan_mat[RRSO_ref_indices,]
RRSO_ref_stat <- PanCan_stat[RRSO_ref_indices,]
RRSO_ref_exposure <- sig_exposure(SSmatrix = RRSO_ref_mat, sig_type = 'PanCan',threshold = 0.5)
#### sort by CX1
RRSO_ref_exposure <- RRSO_ref_exposure[order(-RRSO_ref_exposure$CX1),]
RRSO_ref_exposure <- select(RRSO_ref_exposure, !SS_sum)
#### plotting
RRSO_ref_exposure_plot <- ref_exposure_plot(ref_exposure_df = RRSO_ref_exposure, sig_type = 'Pan-Cancer', ref_stat_df = RRSO_ref_stat) + labs(title = 'RRSO ffTissue&FFPE') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('Pan-Cancer/cosine_similarity/Figures/ref_RRSO_exposure.pdf', plot = RRSO_ref_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### combined
RRSO_ref_exposure_plot_1 <- RRSO_ref_exposure_plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
HGSC_tumor_exposure_plot_1 <- HGSC_tumor_exposure_plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
PanCan_ref_exposure_plot <- ggarrange(Benign_tissue_exposure_plot, RRSO_ref_exposure_plot_1, HGSC_tumor_exposure_plot_1, ncol = 3, common.legend = 1, legend.grob = get_legend(Benign_tissue_exposure_plot), legend = 'right')
ggsave('Pan-Cancer/cosine_similarity/Figures/ref_exposure_all.pdf', plot = PanCan_ref_exposure_plot, dpi = 600, width = 10, height = 5, units = 'in')



###### 3. panConusig signatures ######
# load the dataframe
panConusig_mat <- readRDS('panConusig/panConusig.SSmatrix.rds')
## exclude the samples
for (sampleID in panConusig_mat$sample) {
  if (!(sampleID %in% sample_df$Sample)) {
    panConusig_mat <- filter(panConusig_mat, sample!=sampleID)
  }
}

# calculate the exposures
panConusig_exposure <- sig_exposure(SSmatrix = panConusig_mat,sig_type = 'panConusig',threshold = 0.5)
saveRDS(panConusig_exposure, 'panConusig/panConusig_exposure.rds', compress = FALSE)
# output dataframe
panConusig_df <- select_top_sig(SSmatrix = panConusig_mat, number_of_top = 2, sig_type = 'panConusig')
# calculate the similarity differences between the top 2 signatures
panConusig_out <- delta_val_cal(input_df = panConusig_df, sig_type = 'panConusig', number_of_top = 2)
panConusig_out <- add_df_info(SS_df = panConusig_out, sample_df = sample_df)
panConusig_out <- panConusig_out %>% select('sample', 'patient':'BH_type', 'enrich_panConusig_1':'delta_val')
# output the dataframe
write.xlsx(panConusig_out, file = 'All_sample_signature.xlsx', sheetName = 'panConusig', append = TRUE, row.names = FALSE)

# panConusig signature stats
panConusig_stat <- read.xlsx('All_sample_signature.xlsx',sheetName = 'panConusig')

## 1. distribution of the delta-value
panConusig_dist_delta_plot <- dis_plot_delta(panConusig_stat)
ggsave(filename = 'panConusig/Figures/delta_val_dist.pdf', plot = panConusig_dist_delta_plot, dpi = 600, width = 10, height = 7, units = 'in')
quantile(panConusig_stat$delta_val, 0.1) # 0.013 top 10%
# delta value cut-off at 0.25
panConusig_stat_25 <- panConusig_stat[which(panConusig_stat$delta_val<0.25 & panConusig_stat$panConusig_similarity_2>=0.5),]
#### top 2 enrich signature
panConusig_info_df_25 <- sig_stat_enrich(stat_df = panConusig_stat_25, enrich_num = 2, sig_type = 'panConusig')
panConusig_info_df_25$enrich_panConusig_2 <- factor(panConusig_info_df_25$enrich_panConusig_2, levels = paste0('CN', c(1:25)))
# save the info dataframe
panConusig_info_df_25 <- as.data.frame(panConusig_info_df_25)
write.xlsx(panConusig_info_df_25, file = 'panConusig/panConusig_group.xlsx', sheetName = 'panConusig_025_2',append = TRUE, row.names = FALSE)

# delta value cut-off at 0.1
panConusig_stat_01 <- panConusig_stat[which(panConusig_stat$delta_val<0.1 & panConusig_stat$panConusig_similarity_2>=0.5),]
#### top 2 enrich signature
panConusig_info_df_01 <- sig_stat_enrich(stat_df = panConusig_stat_01, enrich_num = 2, sig_type = 'panConusig')
panConusig_info_df_01$enrich_panConusig_2 <- factor(panConusig_info_df_01$enrich_panConusig_2, levels = paste0('CN', c(1:25)))
# save the info dataframe
panConusig_info_df_01 <- as.data.frame(panConusig_info_df_01)
write.xlsx(panConusig_info_df_01, file = 'panConusig/panConusig_group.xlsx', sheetName = 'panConusig_01_2',append = TRUE, row.names = FALSE)

# after adjusting the delta-value cutoff
## 0.25 for Endome
## 0.1 for all other groups
########## top 1 ##########
panConusig_stat_01 <- panConusig_stat[which(panConusig_stat$delta_val>=0.1),]
panConusig_stat_01 <- panConusig_stat_01 %>% filter(group!='Endome')
panConusig_stat_25 <- panConusig_stat[which(panConusig_stat$delta_val>=0.25),]
panConusig_stat_25 <- panConusig_stat_25 %>% filter(group=='Endome')
panConusig_stat_cut <- rbind(panConusig_stat_01, panConusig_stat_25)
#### top 1 enrich signature
panConusig_info_df_cut <- sig_stat_enrich(stat_df = panConusig_stat_cut, enrich_num = 1, sig_type = 'panConusig')
panConusig_info_df_cut$enrich_panConusig_1 <- factor(panConusig_info_df_cut$enrich_panConusig_1, levels = paste0('CN',c(1:25)))
# plot the graph
panConusig_enrich_plot_1 <- ggplot(data = panConusig_info_df_cut, aes(x=BH, y=percentage, fill=enrich_panConusig_1)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='panConusig signatures') +
  scale_fill_ucscgb() +
  theme_classic() + facet_grid(cols = vars(group_time))
panConusig_enrich_plot_1 <- panConusig_enrich_plot_1 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('panConusig/Figures/panConusig_sig_BH_cut_1.pdf', plot = panConusig_enrich_plot_1, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
panConusig_info_df_cut <- as.data.frame(panConusig_info_df_cut)
write.xlsx(panConusig_info_df_cut, file = 'panConusig/panConusig_group.xlsx', sheetName = 'panConusig_cut_1',append = TRUE, row.names = FALSE)

########## top 2 ##########
panConusig_stat_01 <- panConusig_stat[which(panConusig_stat$delta_val<0.1),]
panConusig_stat_01 <- panConusig_stat_01 %>% filter(group!='Endome')
panConusig_stat_25 <- panConusig_stat[which(panConusig_stat$delta_val<0.25),]
panConusig_stat_25 <- panConusig_stat_25 %>% filter(group=='Endome')
panConusig_stat_cut <- rbind(panConusig_stat_01, panConusig_stat_25)
#### top 2 enrich signature
panConusig_info_df_cut <- sig_stat_enrich(stat_df = panConusig_stat_cut, enrich_num = 2, sig_type = 'panConusig')
panConusig_info_df_cut$enrich_panConusig_2 <- factor(panConusig_info_df_cut$enrich_panConusig_2, levels = paste0('CN', c(1:25)))
# plot the graph
panConusig_enrich_plot_2 <- ggplot(data = panConusig_info_df_cut, aes(x=BH, y=percentage, fill=enrich_panConusig_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='panConusig signatures') +
  scale_fill_ucscgb() +
  theme_classic() + facet_grid(cols = vars(group_time))
panConusig_enrich_plot_2 <- panConusig_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('panConusig/Figures/panConusig_sig_BH_cut_2.pdf', plot = panConusig_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
panConusig_info_df_cut <- as.data.frame(panConusig_info_df_cut)
write.xlsx(panConusig_info_df_cut, file = 'panConusig/panConusig_group.xlsx', sheetName = 'panConusig_cut_2',append = TRUE, row.names = FALSE)


## 2. exposure distribution
### rank by group
panConusig_color_bar_group <- color_bar_prep(stat_df = panConusig_stat, exposure_df = panConusig_exposure, grouping_type = 'group_VS', sig_type = 'panConusig', sort = TRUE) + labs(fill = 'Types')
panConusig_exposure_plot_group <- plot_exposure(stat_df = panConusig_stat, exposure_df = panConusig_exposure, grouping_type = 'group_VS', sig_type = 'panConusig', sort = TRUE)
panConusig_exposure_plot_group <- panConusig_exposure_plot_group %>% insert_bottom(panConusig_color_bar_group, height = 0.03)
ggsave(filename = 'panConusig/Figures/exposure_plot_groups_sort_v1.pdf', plot = panConusig_exposure_plot_group, dpi = 600, width = 40, height = 7, units = 'in')
### rank by BH
panConusig_normal_stat <- panConusig_stat %>% filter(group_VS == 'Normal tissues')
normal_sample <- panConusig_normal_stat$sample
panConusig_normal_exposure <- panConusig_exposure %>% filter(sample %in% normal_sample)
panConusig_color_bar_BH <- color_bar_prep(stat_df = panConusig_normal_stat, exposure_df = panConusig_normal_exposure, grouping_type = 'BH_group', sig_type = 'panConusig', sort = TRUE)
panConusig_exposure_plot_BH <- plot_exposure(stat_df = panConusig_normal_stat, exposure_df = panConusig_normal_exposure, grouping_type = 'BH_group', sig_type = 'panConusig', sort = TRUE)
panConusig_exposure_plot_BH <- panConusig_exposure_plot_BH %>% insert_bottom(panConusig_color_bar_BH, height = 0.03)
ggsave(filename = 'panConusig/Figures/exposure_plot_BH_normal.pdf', plot = panConusig_exposure_plot_BH, dpi = 600, width = 30, height = 7, units = 'in')


## 3. panConusig signatures enrichment
### in this part, combine ffTumor and ffTissue into ffTissue
#### top 1 enriched signature
panConusig_info_df <- sig_stat_enrich(stat_df = panConusig_stat, enrich_num = 1, sig_type = 'panConusig')
panConusig_info_df$enrich_panConusig_1 <- factor(panConusig_info_df$enrich_panConusig_1, levels = paste0('CN',c(1:25)))
# plot the graph
panConusig_enrich_plot_1 <- ggplot(data = panConusig_info_df, aes(x=BH, y=percentage, fill=enrich_panConusig_1)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='panConusig signatures') +
  scale_fill_ucscgb() +
  theme_classic() + facet_grid(cols = vars(group_time))
panConusig_enrich_plot_1 <- panConusig_enrich_plot_1 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('panConusig/Figures/panConusig_sig_BH_all.pdf', plot = panConusig_enrich_plot_1, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
panConusig_info_df <- as.data.frame(panConusig_info_df)
write.xlsx(panConusig_info_df, file = 'panConusig/panConusig_group.xlsx', sheetName = 'panConusig_1',append = TRUE, row.names = FALSE)

#### top 2 enrich signature
panConusig_info_df <- sig_stat_enrich(stat_df = panConusig_stat, enrich_num = 2, sig_type = 'panConusig')
panConusig_info_df$enrich_panConusig_2 <- factor(panConusig_info_df$enrich_panConusig_2, levels = paste0('CN',c(1:25)))
# plot the graph
panConusig_enrich_plot_2 <- ggplot(data = panConusig_info_df, aes(x=BH, y=percentage, fill=enrich_panConusig_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='panConusig signatures') +
  scale_fill_ucscgb() +
  theme_classic() + facet_grid(cols = vars(group_time))
panConusig_enrich_plot_2 <- panConusig_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('panConusig/Figures/panConusig_sig_BH_all_2.pdf', plot = panConusig_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
panConusig_info_df <- as.data.frame(panConusig_info_df)
write.xlsx(panConusig_info_df, file = 'panConusig/panConusig_group.xlsx', sheetName = 'panConusig_2',append = TRUE, row.names = FALSE)


## 4. VS samples time point exposure
VS_sheet <- read.xlsx(file = 'MaNiLa_All_samplesheet_240311_groups.xlsx', sheetName = 'VS_samples')
VS_stat <- panConusig_stat[which(panConusig_stat$group=='ArchivalVS' | panConusig_stat$group=='MaNiLaVS'),]
for (sampleID in VS_stat$sample) {
  VS_stat[which(VS_stat$sample==sampleID), 'Diagnostic_time'] <- VS_sheet[which(VS_sheet$Library==sampleID), 'Diagnostic_time']
}

### extract the similarity matrix
indices <- rownames(VS_stat)
VS_mat <- panConusig_mat[indices,]
VS_exposure <- sig_exposure(SSmatrix = VS_mat,sig_type = 'panConusig',threshold = 0.5)
VS_color_bar_BH <- color_bar_prep(stat_df = VS_stat, exposure_df = VS_exposure, grouping_type = 'BH_group', sig_type = 'panConusig')
VS_exposure_plot_BH <- plot_exposure(stat_df = VS_stat, exposure_df = VS_exposure, grouping_type = 'BH_group', sig_type = 'panConusig')
VS_exposure_plot_BH <- VS_exposure_plot_BH %>% insert_bottom(VS_color_bar_BH, height = 0.03)
ggsave('panConusig/cosine_similarity/Figures/VS_exposure.pdf', plot = VS_exposure_plot_BH, dpi = 600, width = 10, height = 7, units = 'in')

### for Benign VS samples
Benign_VS_stat <- VS_stat[which(VS_stat$BH=='Benign'),]
unique(Benign_VS_stat$Diagnostic_time) # only "diagnostic"
indices <- rownames(Benign_VS_stat)
Benign_VS_mat <- VS_mat[indices,]
#### exposure
Benign_VS_exposure <- sig_exposure(SSmatrix = Benign_VS_mat,sig_type = 'panConusig',threshold = 0.5)
#### sort by CN1
Benign_VS_exposure <- Benign_VS_exposure[order(-Benign_VS_exposure$CN1),]
Benign_VS_exposure <- select(Benign_VS_exposure,!SS_sum)
#### plotting
Benign_VS_plot <- VS_exposure_plot(VS_exposure_df = Benign_VS_exposure, sig_type = 'panConusig', VS_stat_df = Benign_VS_stat) + labs(title = 'diagnostic') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('panConusig/cosine_similarity/Figures/Benign_VS_exposure.pdf', plot = Benign_VS_plot, dpi = 600, width = 7, height = 5, units = 'in')

### for RRSO samples
RRSO_VS_stat <- VS_stat[which(VS_stat$BH=='RRSO'),]
unique(RRSO_VS_stat$Diagnostic_time) # "prediagnostic_6+", "diagnostic"
indices <- rownames(RRSO_VS_stat)
RRSO_VS_mat <- VS_mat[indices,]
#### exposure
RRSO_VS_exposure_in <- sig_exposure(SSmatrix = RRSO_VS_mat,sig_type = 'panConusig',threshold = 0.5)
#### sort by CN1
RRSO_VS_exposure <- as.data.frame(matrix(nrow = 0, ncol = 27))
diag_group_num <- c()
for (diagnostic_group in c("prediagnostic_6+", "diagnostic")) {
  indices_val <- which(RRSO_VS_stat$Diagnostic_time==diagnostic_group)
  sort_df <- RRSO_VS_exposure_in[indices_val,]
  sort_df <- sort_df[order(-sort_df$CN1),]
  RRSO_VS_exposure <- rbind(RRSO_VS_exposure, sort_df)
  diag_group_num <- c(diag_group_num, nrow(sort_df))
}
RRSO_VS_exposure <- select(RRSO_VS_exposure,!SS_sum)
# diag_group_num == c(12, 8)
pre_6_RRSO_exposure <- RRSO_VS_exposure[1:12,]
diag_RRSO_exposure <- RRSO_VS_exposure[13:20,]
#### plotting
#1
pre_6_RRSO_exposure_plot <- VS_exposure_plot(VS_exposure_df = pre_6_RRSO_exposure, sig_type = 'panConusig', VS_stat_df = RRSO_VS_stat) + labs(title = 'prediagnostic > 6m') + theme(plot.title = element_text(size = 10, hjust = 0.5))
#2
diag_RRSO_exposure_plot <- VS_exposure_plot(VS_exposure_df = diag_RRSO_exposure, sig_type = 'panConusig', VS_stat_df = RRSO_VS_stat)  + labs(title = 'diagnostic') + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(size = 10, hjust = 0.5))
#combine
RRSO_VS_plot <- ggarrange(pre_6_RRSO_exposure_plot, diag_RRSO_exposure_plot, ncol = 2, common.legend = 1, legend.grob = get_legend(pre_6_RRSO_exposure_plot), legend = 'right')
ggsave('panConusig/cosine_similarity/Figures/RRSO_VS_exposure.pdf', plot = RRSO_VS_plot, dpi = 600, width = 10, height = 5, units = 'in')

### for HGSC samples
HGSC_VS_stat <- VS_stat[which(VS_stat$BH=='HGSC'),]
unique(HGSC_VS_stat$Diagnostic_time) # "prediagnostic_0-6", "prediagnostic_6+", "diagnostic"
indices <- rownames(HGSC_VS_stat)
HGSC_VS_mat <- VS_mat[indices,]
#### exposure
HGSC_VS_exposure_in <- sig_exposure(SSmatrix = HGSC_VS_mat,sig_type = 'panConusig',threshold = 0.5)
#### sort by CN1
HGSC_VS_exposure <- as.data.frame(matrix(nrow = 0, ncol = 27))
diag_group_num <- c()
for (diagnostic_group in c("prediagnostic_6+", "prediagnostic_0-6", "diagnostic")) {
  indices_val <- which(HGSC_VS_stat$Diagnostic_time==diagnostic_group)
  sort_df <- HGSC_VS_exposure_in[indices_val,]
  sort_df <- sort_df[order(-sort_df$CN1),]
  HGSC_VS_exposure <- rbind(HGSC_VS_exposure, sort_df)
  diag_group_num <- c(diag_group_num, nrow(sort_df))
}
HGSC_VS_exposure <- select(HGSC_VS_exposure,!SS_sum)
# diag_group_num == c(38,7,22)
pre_6_HGSC_exposure <- HGSC_VS_exposure[1:38,]
pre_0_HGSC_exposure <- HGSC_VS_exposure[39:45,]
diag_HGSC_exposure <- HGSC_VS_exposure[46:67,]
#### plotting
#1
pre_6_HGSC_exposure_plot <- VS_exposure_plot(VS_exposure_df = pre_6_HGSC_exposure, sig_type = 'panConusig', VS_stat_df = HGSC_VS_stat) + labs(title = 'prediagnostic > 6m') + theme(plot.title = element_text(size = 10, hjust = 0.5))
#2
pre_0_HGSC_exposure_plot <- VS_exposure_plot(VS_exposure_df = pre_0_HGSC_exposure, sig_type = 'panConusig', VS_stat_df = HGSC_VS_stat)  + labs(title = 'prediagnostic 0-6m') + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(size = 10, hjust = 0.5))
#3
diag_HGSC_exposure_plot <- VS_exposure_plot(VS_exposure_df = diag_HGSC_exposure, sig_type = 'panConusig', VS_stat_df = HGSC_VS_stat)  + labs(title = 'diagnostic') + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(size = 10, hjust = 0.5))
#combine
HGSC_VS_plot <- ggarrange(pre_6_HGSC_exposure_plot, pre_0_HGSC_exposure_plot, diag_HGSC_exposure_plot, ncol = 3, common.legend = 1, legend.grob = get_legend(pre_6_HGSC_exposure_plot), legend = 'right')
ggsave('panConusig/cosine_similarity/Figures/HGSC_VS_exposure.pdf', plot = HGSC_VS_plot, dpi = 600, width = 10, height = 5, units = 'in')

## 5. reference exposure plots
#### HGSC ffTumor 
HGSC_tumor_indices <- which(panConusig_stat$group=='ffTumor' & panConusig_stat$BH=='HGSC')
HGSC_tumor_mat <- panConusig_mat[HGSC_tumor_indices,]
HGSC_tumor_stat <- panConusig_stat[HGSC_tumor_indices,]
HGSC_tumor_exposure <- sig_exposure(SSmatrix = HGSC_tumor_mat, sig_type = 'panConusig',threshold = 0.5)
#### sort by CN1
HGSC_tumor_exposure <- HGSC_tumor_exposure[order(-HGSC_tumor_exposure$CN1),]
HGSC_tumor_exposure <- select(HGSC_tumor_exposure,!SS_sum)
#### plotting
HGSC_tumor_exposure_plot <- ref_exposure_plot(ref_exposure_df = HGSC_tumor_exposure, sig_type = 'panConusig', ref_stat_df = HGSC_tumor_stat) + labs(title = 'HGSC ffTumor') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('panConusig/cosine_similarity/Figures/ref_HGSC_exposure.pdf', plot = HGSC_tumor_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### Benign ffTissue
Benign_tissue_indices <- which(panConusig_stat$group=='ffTissue' & panConusig_stat$BH=='Benign')
Benign_tissue_mat <- panConusig_mat[Benign_tissue_indices,]
Benign_tissue_stat <- panConusig_stat[Benign_tissue_indices,]
Benign_tissue_exposure <- sig_exposure(SSmatrix = Benign_tissue_mat, sig_type = 'panConusig',threshold = 0.5)
#### sort by CN1
Benign_tissue_exposure <- Benign_tissue_exposure[order(-Benign_tissue_exposure$CN1),]
Benign_tissue_exposure <- select(Benign_tissue_exposure, !SS_sum)
#### plotting
Benign_tissue_exposure_plot <- ref_exposure_plot(ref_exposure_df = Benign_tissue_exposure, sig_type = 'panConusig', ref_stat_df = Benign_tissue_stat) + labs(title = 'Benign ffTissue') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('panConusig/cosine_similarity/Figures/ref_Benign_exposure.pdf', plot = Benign_tissue_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### RRSO ffTissue+FFPE
RRSO_ref_indices <- which(panConusig_stat$BH=='RRSO' & (panConusig_stat$group=='ffTissue' | panConusig_stat$group=='FFPE'))
RRSO_ref_mat <- panConusig_mat[RRSO_ref_indices,]
RRSO_ref_stat <- panConusig_stat[RRSO_ref_indices,]
RRSO_ref_exposure <- sig_exposure(SSmatrix = RRSO_ref_mat, sig_type = 'panConusig',threshold = 0.5)
#### sort by CN1
RRSO_ref_exposure <- RRSO_ref_exposure[order(-RRSO_ref_exposure$CN1),]
RRSO_ref_exposure <- select(RRSO_ref_exposure, !SS_sum)
#### plotting
RRSO_ref_exposure_plot <- ref_exposure_plot(ref_exposure_df = RRSO_ref_exposure, sig_type = 'panConusig', ref_stat_df = RRSO_ref_stat) + labs(title = 'RRSO ffTissue&FFPE') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('panConusig/cosine_similarity/Figures/ref_RRSO_exposure.pdf', plot = RRSO_ref_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### combined
RRSO_ref_exposure_plot_1 <- RRSO_ref_exposure_plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
HGSC_tumor_exposure_plot_1 <- HGSC_tumor_exposure_plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
panConusig_ref_exposure_plot <- ggarrange(Benign_tissue_exposure_plot, RRSO_ref_exposure_plot_1, HGSC_tumor_exposure_plot_1, ncol = 3, common.legend = 1, legend.grob = get_legend(Benign_tissue_exposure_plot), legend = 'right')
ggsave('panConusig/cosine_similarity/Figures/ref_exposure_all.pdf', plot = panConusig_ref_exposure_plot, dpi = 600, width = 10, height = 5, units = 'in')



# 4. CerCNsig
# load the dataframe
CerCN_mat <- readRDS('CerCNsig/CerCNsig.All_Cervical.SSmatrix.rds')
## exclude the samples
for (sampleID in CerCN_mat$sample) {
  if (!(sampleID %in% sample_df$Sample)) {
    CerCN_mat <- filter(CerCN_mat, sample!=sampleID)
  }
}

# calculate the exposures
CerCN_exposure <- sig_exposure(SSmatrix = CerCN_mat,sig_type = 'CerCN',threshold=0.5)
saveRDS(CerCN_exposure, file = 'CerCNsig/CerCNsig_exposure.rds',compress = FALSE)
# output dataframe
CerCN_df <- select_top_sig(SSmatrix = CerCN_mat, number_of_top = 2, sig_type = 'CerCN')
# calculate the similarity differences between the top 2 signatures
CerCN_out <- delta_val_cal(input_df = CerCN_df, sig_type = 'CerCN', number_of_top = 2)
# add the sample information (patient, type, group, BH, BH_type)
CerCN_out <- add_df_info(SS_df = CerCN_out, sample_df = sample_df)
CerCN_out <- CerCN_out %>% select('sample', 'patient':'BH_type', 'enrich_CerCN_1':'delta_val')
# output the dataframe
write.xlsx(CerCN_out, file = 'All_sample_signature.xlsx', sheetName = 'CerCN', append = TRUE, row.names = FALSE)

# CN signature info
CerCN_stat <- read.xlsx('All_sample_signature.xlsx',sheetName = 'CerCN')

# VS samples time point exposure
for (sampleID in CerCN_stat$sample) {
  CerCN_stat[which(CerCN_stat$sample==sampleID), 'Diagnostic_time'] <- sample_df[which(sample_df$Sample==sampleID), 'Group_time_detail']
}

### for Benign VS samples
Benign_VS_stat <- CerCN_stat[which(CerCN_stat$BH=='Benign'),]
unique(Benign_VS_stat$Diagnostic_time) # only "diagnostic"
indices <- rownames(Benign_VS_stat)
Benign_VS_mat <- CerCN_mat[indices,]
#### exposure
Benign_VS_exposure <- sig_exposure(SSmatrix = Benign_VS_mat,sig_type = 'CerCN',threshold = 0.5)
#### sort by CerCN1
Benign_VS_exposure <- Benign_VS_exposure[order(-Benign_VS_exposure$CerCN1),]
Benign_VS_exposure <- select(Benign_VS_exposure,!SS_sum)
#### plotting
Benign_VS_plot <- VS_exposure_plot(VS_exposure_df = Benign_VS_exposure, sig_type = 'CerCN', VS_stat_df = Benign_VS_stat) + labs(title = 'diagnostic') + theme(plot.title = element_text(size = 15, hjust = 0.5))
ggsave('CerCNsig/Figures/Benign_VS_exposure.pdf', plot = Benign_VS_plot, dpi = 600, width = 7, height = 5, units = 'in')

### for RRSO samples
RRSO_VS_stat <- CerCN_stat[which(CerCN_stat$BH=='RRSO'),]
unique(RRSO_VS_stat$Diagnostic_time) # "prediagnostic_6+", "diagnostic"
indices <- rownames(RRSO_VS_stat)
RRSO_VS_mat <- CerCN_mat[indices,]
#### exposure
RRSO_VS_exposure_in <- sig_exposure(SSmatrix = RRSO_VS_mat,sig_type = 'CerCN',threshold = 0.5)
#### sort by CerCN1
RRSO_VS_exposure <- as.data.frame(matrix(nrow = 0, ncol = 9))
diag_group_num <- c()
for (diagnostic_group in c("prediagnostic_6+", "diagnostic")) {
  indices_val <- which(RRSO_VS_stat$Diagnostic_time==diagnostic_group)
  sort_df <- RRSO_VS_exposure_in[indices_val,]
  sort_df <- sort_df[order(-sort_df$CerCN1),]
  RRSO_VS_exposure <- rbind(RRSO_VS_exposure, sort_df)
  diag_group_num <- c(diag_group_num, nrow(sort_df))
}
RRSO_VS_exposure <- select(RRSO_VS_exposure,!SS_sum)
# diag_group_num == c(13, 17)
pre_6_RRSO_exposure <- RRSO_VS_exposure[1:13,]
diag_RRSO_exposure <- RRSO_VS_exposure[14:30,]
#### plotting
pre_6_RRSO_exposure_plot <- VS_exposure_plot(VS_exposure_df = pre_6_RRSO_exposure, sig_type = 'CerCN', VS_stat_df = RRSO_VS_stat) + labs(title = 'prediagnostic > 6m') + theme(plot.title = element_text(size = 15, hjust = 0.5))
diag_RRSO_exposure_plot <- VS_exposure_plot(VS_exposure_df = diag_RRSO_exposure, sig_type = 'CerCN', VS_stat_df = RRSO_VS_stat)  + labs(title = 'diagnostic') + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(size = 15, hjust = 0.5))
RRSO_VS_plot <- ggarrange(pre_6_RRSO_exposure_plot, diag_RRSO_exposure_plot, ncol = 2, common.legend = 1, legend.grob = get_legend(pre_6_RRSO_exposure_plot), legend = 'right')
ggsave('CerCNsig/Figures/RRSO_VS_exposure.pdf', plot = RRSO_VS_plot, dpi = 600, width = 7, height = 5, units = 'in')

### for HGSC samples
HGSC_VS_stat <- CerCN_stat[which(CerCN_stat$BH=='HGSC'),]
unique(HGSC_VS_stat$Diagnostic_time) # "prediagnostic_0-6", "prediagnostic_6+", "diagnostic"
indices <- rownames(HGSC_VS_stat)
HGSC_VS_mat <- CerCN_mat[indices,]
#### exposure
HGSC_VS_exposure_in <- sig_exposure(SSmatrix = HGSC_VS_mat,sig_type = 'CerCN',threshold = 0.5)
#### sort by CerCN1
HGSC_VS_exposure <- as.data.frame(matrix(nrow = 0, ncol = 9))
diag_group_num <- c()
for (diagnostic_group in c("prediagnostic_6+", "prediagnostic_0-6", "diagnostic")) {
  indices_val <- which(HGSC_VS_stat$Diagnostic_time==diagnostic_group)
  sort_df <- HGSC_VS_exposure_in[indices_val,]
  sort_df <- sort_df[order(-sort_df$CerCN1),]
  HGSC_VS_exposure <- rbind(HGSC_VS_exposure, sort_df)
  diag_group_num <- c(diag_group_num, nrow(sort_df))
}
HGSC_VS_exposure <- select(HGSC_VS_exposure,!SS_sum)
# diag_group_num == c(62,30,36)
pre_6_HGSC_exposure <- HGSC_VS_exposure[1:62,]
pre_0_HGSC_exposure <- HGSC_VS_exposure[63:92,]
diag_HGSC_exposure <- HGSC_VS_exposure[93:128,]
#### plotting
#1
pre_6_HGSC_exposure_plot <- VS_exposure_plot(VS_exposure_df = pre_6_HGSC_exposure, sig_type = 'CerCN', VS_stat_df = HGSC_VS_stat) + labs(title = 'prediagnostic > 6m') + theme(plot.title = element_text(size = 15, hjust = 0.5))
#2
pre_0_HGSC_exposure_plot <- VS_exposure_plot(VS_exposure_df = pre_0_HGSC_exposure, sig_type = 'CerCN', VS_stat_df = HGSC_VS_stat)  + labs(title = 'prediagnostic 0-6m') + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(size = 15, hjust = 0.5))
#3
diag_HGSC_exposure_plot <- VS_exposure_plot(VS_exposure_df = diag_HGSC_exposure, sig_type = 'CerCN', VS_stat_df = HGSC_VS_stat)  + labs(title = 'diagnostic') + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(size = 15, hjust = 0.5))
#combine
HGSC_VS_plot <- ggarrange(pre_6_HGSC_exposure_plot, pre_0_HGSC_exposure_plot, diag_HGSC_exposure_plot, ncol = 3, common.legend = 1, legend.grob = get_legend(pre_0_HGSC_exposure_plot), legend = 'right')
ggsave('CerCNsig/Figures/HGSC_VS_exposure.pdf', plot = HGSC_VS_plot, dpi = 600, width = 10, height = 5, units = 'in')
