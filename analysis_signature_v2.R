# Title: analysis_signature.R
# Author: Guyuan TANG
# Date: 2024/1/5 - 2024/3/20

# Description: the script is used for analyzing the sample-by-component matrices from the workflow.

library(tidyverse)
library(ggpubr)
library(xlsx)
library(aplot)
library(ggsci)

# set the working directory
setwd("E:/1Lund_Lectures/BINP52_MasterProject/Workflow/Draft/signatures")

# the required groups
groups <- c('ArchivalVS','MaNiLaVS','ffTumor','ffTissue','Endome','ffpe','Blood','ffPlasma')

# Generate the All_sample_signature.vcsv (only contains CN_sig and PanCan_sig)
## for panConusig, only ffTumor samples ran the analysis
# load the sample table
sample_df <- read.table('E:/1Lund_Lectures/BINP52_MasterProject/Workflow/Draft/results/solutions/solution_sample.tsv', sep = '\t', header = 1) %>% select('Sample', 'Patient', 'Type', 'Group')
sample_sheet <- readxl::read_excel('MaNiLa_All_samplesheet_240311_groups.xlsx')
## combine the benign and HGSC information with the samples
sample_df$BH <- NA
sample_df$BH_type <- NA
for (i in 1:nrow(sample_df)) {
  sampleID <- sample_df[i, 'Sample']
  sample_df[i,'BH_type'] <- sample_sheet[which(sample_sheet$Library == sampleID), 'Group']
  if (sample_df[i,'BH_type'] %in% c('Benign', 'Benign_TP53_ddPCR')) {
    sample_df[i,'BH'] <- 'Benign'
  } else if (sample_df[i,'BH_type'] %in% c('HGSC', 'HGSC_gBRCA', 'HGSC_sBRCA')) {
    sample_df[i,'BH'] <- 'HGSC'
  } else if (sample_df[i,'BH_type'] %in% c('BRCA_RRSO', 'RRSO')) {
    sample_df[i,'BH'] <- 'RRSO'
  }
    else {
    sample_df[i,'BH'] <- 'exclude'
  }
}
## exclude the 'corpus' and 'Restaging op' sample (CS3_63 & CS3_72)
sample_df <- sample_df %>% filter(BH != 'exclude')

# BH
BH_group <- c('Benign', 'RRSO','HGSC')


# define a function to calculate the similarity differences between the top 2 signatures
delta_val_cal <- function(input_df, sig_type, number_of_top) {
  if (number_of_top<2){
    print('Number of top values should be at least 2.')
  }
  if (sig_type == 'CN') {
    output_df <- input_df %>% mutate(delta_val = CN_similarity_1 - CN_similarity_2)
  } else if (sig_type == 'PanCan') {
    output_df <- input_df %>% mutate(delta_val = PanCan_similarity_1 - PanCan_similarity_2)
  } else if (sig_type == 'panConusig') {
    output_df <- input_df %>% mutate(delta_val = panConusig_similarity_1 - panConusig_similarity_2)
  }
    
  return(output_df)
}



# define a function to extract the top two most similar signature
select_top_sig <- function(SSmatrix, number_of_top=1, sig_type) {
  SSmatrix <- subset(SSmatrix, select = -enrich)
  output_df <- as.data.frame(matrix(ncol = number_of_top*2+1, nrow = nrow(SSmatrix)))
  colnames(output_df)[1] <- 'sample'
  output_df$sample <- SSmatrix$sample
  for (i in 1:number_of_top) {
    enrich_col <- paste0('enrich_',sig_type,'_',i)
    sim_col <- paste0(sig_type,'_similarity_',i)
    colnames(output_df)[i*2] <- enrich_col
    colnames(output_df)[i*2+1] <- sim_col
  }
  # find the top similar signatures
  for (sampleID in output_df$sample) {
    sig_names <- colnames(SSmatrix)[-1]
    line_content <- as.vector.data.frame(SSmatrix[which(SSmatrix$sample==sampleID),-1])
    line_content_vec <- as.vector(unlist(line_content))
    sort_line <- sort(line_content_vec, decreasing = T)
    for (n in 1:number_of_top) {
      top_val <- sort_line[n]
      for (sig_n in sig_names) {
        if (line_content[sig_n] == top_val) {
          enrich_col <- paste0('enrich_',sig_type,'_',n)
          sim_col <- paste0(sig_type,'_similarity_',n)
          output_df[which(output_df$sample==sampleID),enrich_col] <- sig_n
          output_df[which(output_df$sample==sampleID),sim_col] <- top_val
        }
      }
    }
  }
  return(output_df)
}



# define a function to calculate the similarity exposure (percentage) matrix
sig_exposure <- function(SSmatrix, sig_type) {
  # remove the column 'enrich'
  SSmatrix <- SSmatrix[,-2]
  if (sig_type == 'CN') {
    n = 7+1
  } else if (sig_type == 'PanCan') {
    n = 17+1
  } else if (sig_type == 'panConusig') {
    n = 25+1
  }
  # calculate the exposures (normalization)
  out_df <- SSmatrix %>% mutate(SS_sum = rowSums(.[2:n]))
  for (sampleID in out_df$sample) {
    for (i in 2:n) {
      out_df[which(out_df$sample==sampleID),i] <- out_df[which(out_df$sample==sampleID),i] / out_df[which(out_df$sample==sampleID),'SS_sum']
    }
  }
  return(out_df)
}

# define a function to add information to the output signature dataframe
add_df_info <- function(SS_df, sample_df) {
  for (sampleID in sample_df$Sample) {
    if (sampleID %in% SS_df$sample) {
      SS_df[which(SS_df$sample==sampleID),'patient'] <- sample_df[which(sample_df$Sample==sampleID), 'Patient']
      SS_df[which(SS_df$sample==sampleID),'type'] <- sample_df[which(sample_df$Sample==sampleID), 'Type']
      SS_df[which(SS_df$sample==sampleID),'group'] <- sample_df[which(sample_df$Sample==sampleID), 'Group']
      SS_df[which(SS_df$sample==sampleID),'BH'] <- sample_df[which(sample_df$Sample==sampleID), 'BH']
      SS_df[which(SS_df$sample==sampleID),'BH_type'] <- sample_df[which(sample_df$Sample==sampleID), 'BH_type']
    }
  }
  # remove the samples that do not have a signatures
  SS_df <- filter(SS_df, !is.na(patient))
}


# define a function to draw the distribution plot for delta-value
dis_plot_delta <- function(stat_df) {
  p <- ggplot(stat_df, aes(x=delta_val)) + 
      geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                    binwidth=.5,
                    colour="black", fill="white") +
      scale_x_continuous(limits = c(0,1)) +
      geom_density(alpha=.2, fill="#FF6666")
  return(p)
}

# define a function to prepare the color bars
color_bar_prep <- function(stat_df, exposure_df, grouping_type, sig_type) {
  col_num <- ncol(exposure_df)
  ### create a dataframe to store the information, and rank by either groups or BH_groups
  exposure_out <- as.data.frame(matrix(nrow = 0,ncol = col_num))
  if (grouping_type=='groups') {
    groups <- c('ArchivalVS','MaNiLaVS','ffTumor','ffTissue','Endome','ffpe','Blood','ffPlasma')
    ### rank by groups
    x_lab_name <- 'Sample types'
    for (groupID in groups) {
      indices <- which(stat_df$group==groupID)
      for (i in indices) {
        exposure_out <- rbind(exposure_out, exposure_df[i,])
      }
    }
  } else if (grouping_type == 'BH_group') {
    BH_group <- c('Benign', 'RRSO','HGSC')
    ### rank by BH_group
    x_lab_name <- 'Groups'
    for (BH_ID in BH_group) {
      indices <- which(stat_df$BH==BH_ID)
      for (i in indices) {
        exposure_out <- rbind(exposure_out, exposure_df[i,])
      }
    }
  }
  
  exposure_out <- select(exposure_out, !SS_sum)
  rownames(exposure_out) <- 1:nrow(exposure_out)
  # match the grouping information in order to add the coloring bars
  if (grouping_type=='groups') {
    for (sampleID in exposure_out$sample) {
      exposure_out[which(exposure_out$sample==sampleID),'col_group'] <- stat_df[which(stat_df$sample==sampleID),'group']
    }
    exposure_out$col_group <- factor(exposure_out$col_group, levels = groups)
  } else if (grouping_type=='BH_group') {
    for (sampleID in exposure_out$sample) {
      exposure_out[which(exposure_out$sample==sampleID),'col_group'] <- stat_df[which(stat_df$sample==sampleID),'BH']
    }
    exposure_out$col_group <- factor(exposure_out$col_group, levels = BH_group)
  }
  
  # draw the color plot
  
  exposure_out$sample <- factor(exposure_out$sample, levels = exposure_out$sample)
  exposure_out <- exposure_out %>% rename(c('Groups' = 'col_group'))
  col_plot <- ggplot(exposure_out,aes(x=sample,y=1))+
    geom_tile(aes(fill=Groups))+
    scale_y_continuous(expand = c(0,0))+
    theme(panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10)) +
    scale_fill_brewer(palette = 'Set2')
  return(col_plot)
}

# define a function to plot the exposure of signatures
plot_exposure <- function(stat_df, exposure_df, grouping_type, sig_type) {
  # col25 <- colorRampPalette(RColorBrewer::brewer.pal(8,'Dark2'))(25)
  
  col_num <- ncol(exposure_df)
  ### create a dataframe to store the information, and rank by either groups or BH_groups
  exposure_out <- as.data.frame(matrix(nrow = 0,ncol = col_num))
  if (grouping_type=='groups') {
    groups <- c('ArchivalVS','MaNiLaVS','ffTumor','ffTissue','Endome','ffpe','Blood','ffPlasma')
    ### rank by groups
    x_lab_name <- 'Sample types'
    for (groupID in groups) {
      indices <- which(stat_df$group==groupID)
      for (i in indices) {
        exposure_out <- rbind(exposure_out, exposure_df[i,])
      }
    }
  } else if (grouping_type == 'BH_group') {
    BH_group <- c('Benign', 'RRSO','HGSC')
    ### rank by BH_group
    x_lab_name <- 'Groups'
    for (BH_ID in BH_group) {
      indices <- which(stat_df$BH==BH_ID)
      for (i in indices) {
        exposure_out <- rbind(exposure_out, exposure_df[i,])
      }
    }
  }
  
  exposure_out <- select(exposure_out, !SS_sum)
  rownames(exposure_out) <- 1:nrow(exposure_out)
  # prepare the sample ranking
  sample_rank <- exposure_out$sample
  # prepare the dataframe
  exposure_out <- pivot_longer(exposure_out,!sample, names_to = 'signatures',values_to = 'exposure')
  exposure_out$sample <- factor(exposure_out$sample, levels = sample_rank)
  # set the ranking of the signatures
  if (sig_type == 'Pan-Cancer') {
    exposure_out$signatures <- factor(exposure_out$signatures, levels = c('CX1', 'CX2', 'CX3', 'CX4', 'CX5', 'CX6','CX7','CX8','CX9','CX10','CX11','CX12','CX13','CX14','CX15','CX16','CX17'))
  } else if (sig_type == 'CN') {
    exposure_out$signatures <- factor(exposure_out$signatures, levels = c('s1','s2','s3','s4','s5','s6','s7'))
  } else if (sig_type == 'panConusig') {
    exposure_out$signatures <- factor(exposure_out$signatures, levels = c('CN1','CN2','CN3','CN4','CN5','CN6','CN7','CN8','CN9','CN10','CN11','CN12','CN13','CN14','CN15','CN16','CN17','CN18','CN19','CN20','CN21','CN22','CN23','CN24','CN25'))
  }
  
  exposure_plot <- ggplot(data = exposure_out, aes(x=sample, y=exposure, fill=signatures)) +
    geom_bar(stat = 'identity', width = 0.7, position = position_stack(reverse = T)) +
    labs(fill = paste0(sig_type, ' signatures'), y='Exposure') +
    scale_y_continuous(expand = c(0.01,0)) +
    theme(panel.background = element_blank(),
          axis.title.y = element_text(size = 12, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "right") +
    scale_fill_ucscgb()
  return(exposure_plot)
}

# define a function to prepare the signature enrichment statistical analysis
sig_stat_enrich <- function(stat_df, enrich_num, sig_type) {
  ### in this part, combine ffTumor and ffTissue into ffTissue
  groups2 <- c('ArchivalVS','MaNiLaVS','ffTissue','Endome','ffpe','Blood','ffPlasma')
  BH_group <- c('Benign','RRSO','HGSC')
  for (sampleID in stat_df$sample) {
    group_ID <- stat_df[which(stat_df$sample==sampleID),'group']
    if (group_ID=='ffTumor') {
      stat_df[which(stat_df$sample==sampleID),'group2'] <- 'ffTissue'
    } else {
      stat_df[which(stat_df$sample==sampleID),'group2'] <- stat_df[which(stat_df$sample==sampleID),'group']
    }
  }
  if (enrich_num==1 & sig_type=='CN') {
    info_df <- stat_df %>% group_by(group2, BH, enrich_CN_1) %>%
      summarise(mean_CN_similarity = mean(CN_similarity_1),
                sd_CN_similarity = sd(CN_similarity_1),
                sample_size = n())
  } else if (enrich_num==2 & sig_type=='CN') {
    info_df <- stat_df %>% group_by(group2, BH, enrich_CN_2) %>%
      summarise(mean_CN_similarity = mean(CN_similarity_2),
                sd_CN_similarity = sd(CN_similarity_2),
                sample_size = n())
  } else if (enrich_num==1 & sig_type=='PanCan') {
    info_df <- stat_df %>% group_by(group2, BH, enrich_PanCan_1) %>%
      summarise(mean_PanCan_similarity = mean(PanCan_similarity_1),
                sd_PanCan_similarity = sd(PanCan_similarity_1),
                sample_size = n())
  } else if (enrich_num==2 & sig_type=='PanCan') {
    info_df <- stat_df %>% group_by(group2, BH, enrich_PanCan_2) %>%
      summarise(mean_PanCan_similarity = mean(PanCan_similarity_2),
                sd_PanCan_similarity = sd(PanCan_similarity_2),
                sample_size = n())
  } else if (enrich_num==1 & sig_type=='panConusig') {
    info_df <- stat_df %>% group_by(group2, BH, enrich_panConusig_1) %>%
      summarise(mean_panConusig_similarity = mean(panConusig_similarity_1),
                sd_panConusig_similarity = sd(panConusig_similarity_1),
                sample_size = n())
  } else if (enrich_num==2 & sig_type=='panConusig') {
    info_df <- stat_df %>% group_by(group2, BH, enrich_panConusig_2) %>%
      summarise(mean_panConusig_similarity = mean(panConusig_similarity_2),
                sd_panConusig_similarity = sd(panConusig_similarity_2),
                sample_size = n())
  }
  
  # adjust the sd value
  sd_col_name <- paste0('sd_',sig_type,'_similarity')
  for (i in 1:nrow(info_df)) {
    if (is.na(info_df[i,sd_col_name])) {
      info_df[i,sd_col_name] <- 0
    }
  }
  info_df$group2 <- factor(info_df$group2, levels = groups2)
  info_df$BH <- factor(info_df$BH, levels = BH_group)
  for (groupID in groups2) {
    for (BH_ID in BH_group) {
      info_df[which(info_df$group2==groupID & info_df$BH==BH_ID),"sample_sum"] <- sum(info_df[which(info_df$group2==groupID & info_df$BH==BH_ID),"sample_size"])
    }
  }
  info_df <- info_df %>% mutate(percentage = 100 * sample_size / sample_sum)
  return(info_df)
}




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
CN_exposure <- sig_exposure(SSmatrix = CN_mat,sig_type = 'CN')

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

## 2. exposure distribution
### rank by group
CN_color_bar_group <- color_bar_prep(stat_df = CN_stat, exposure_df = CN_exposure, grouping_type = 'groups', sig_type = 'CN')
CN_exposure_plot_group <- plot_exposure(stat_df = CN_stat, exposure_df = CN_exposure, grouping_type = 'groups', sig_type = 'CN')
CN_exposure_plot_group <- CN_exposure_plot_group %>% insert_bottom(CN_color_bar_group, height = 0.03)
ggsave(filename = 'Brenton/Figures/exposure_plot_groups.pdf', plot = CN_exposure_plot_group, dpi = 600, width = 30, height = 7, units = 'in')
### rank by BH
CN_color_bar_BH <- color_bar_prep(stat_df = CN_stat, exposure_df = CN_exposure, grouping_type = 'BH_group', sig_type = 'CN')
CN_exposure_plot_BH <- plot_exposure(stat_df = CN_stat, exposure_df = CN_exposure, grouping_type = 'BH_group', sig_type = 'CN')
CN_exposure_plot_BH <- CN_exposure_plot_BH %>% insert_bottom(CN_color_bar_BH, height = 0.03)
ggsave(filename = 'Brenton/Figures/exposure_plot_BH.pdf', plot = CN_exposure_plot_BH, dpi = 600, width = 30, height = 7, units = 'in')

## 3. CN signatures enrichment
### in this part, combine ffTumor and ffTissue into ffTissue
groups2 <- c('ArchivalVS','MaNiLaVS','ffTissue','Endome','ffpe','Blood','ffPlasma')
#### top 1 enrich signature
CN_info_df <- sig_stat_enrich(stat_df = CN_stat, enrich_num = 1, sig_type = 'CN')

# plot the graph
CN_enrich_plot_1 <- ggplot(data = CN_info_df, aes(x=BH, y=percentage, fill=enrich_CN_1)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='CN signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group2))
CN_enrich_plot_1 <- CN_enrich_plot_1 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Brenton/Figures/CN_sig_BH_all.pdf', plot = CN_enrich_plot_1, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
# write.xlsx(CN_info_df, file = 'All_sample_signature.xlsx', sheetName = 'CN_group_1', append = TRUE, row.names = FALSE)
write_xlsx(CN_info_df, path = 'CN_group_1.xlsx')


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
PanCan_exposure <- sig_exposure(SSmatrix = PanCan_mat,sig_type = 'PanCan')

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

## 2. exposure distribution
### rank by group
PanCan_color_bar_group <- color_bar_prep(stat_df = PanCan_stat, exposure_df = PanCan_exposure, grouping_type = 'groups', sig_type = 'Pan-Cancer')
PanCan_exposure_plot_group <- plot_exposure(stat_df = PanCan_stat, exposure_df = PanCan_exposure, grouping_type = 'groups', sig_type = 'Pan-Cancer')
PanCan_exposure_plot_group <- PanCan_exposure_plot_group %>% insert_bottom(PanCan_color_bar_group, height = 0.03)
ggsave(filename = 'Pan-Cancer/Figures/exposure_plot_groups.pdf', plot = PanCan_exposure_plot_group, dpi = 600, width = 30, height = 7, units = 'in')
### rank by BH
PanCan_color_bar_BH <- color_bar_prep(stat_df = PanCan_stat, exposure_df = PanCan_exposure, grouping_type = 'BH_group', sig_type = 'Pan-Cancer')
PanCan_exposure_plot_BH <- plot_exposure(stat_df = PanCan_stat, exposure_df = PanCan_exposure, grouping_type = 'BH_group', sig_type = 'Pan-Cancer')
PanCan_exposure_plot_BH <- PanCan_exposure_plot_BH %>% insert_bottom(PanCan_color_bar_BH, height = 0.03)
ggsave(filename = 'Pan-Cancer/Figures/exposure_plot_BH.pdf', plot = PanCan_exposure_plot_BH, dpi = 600, width = 30, height = 7, units = 'in')

## 3. CN signatures enrichment
### in this part, combine ffTumor and ffTissue into ffTissue
groups2 <- c('ArchivalVS','MaNiLaVS','ffTissue','Endome','ffpe','Blood','ffPlasma')
#### top 1 enriched signature
PanCan_info_df <- sig_stat_enrich(stat_df = PanCan_stat, enrich_num = 1, sig_type = 'PanCan')
PanCan_info_df$enrich_PanCan_1 <- factor(PanCan_info_df$enrich_PanCan_1, levels = c('CX1', 'CX2', 'CX3', 'CX5', 'CX15', 'CX17'))
# plot the graph
PanCan_enrich_plot_1 <- ggplot(data = PanCan_info_df, aes(x=BH, y=percentage, fill=enrich_PanCan_1)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='Pan-Cancer signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group2))
PanCan_enrich_plot_1 <- PanCan_enrich_plot_1 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Pan-Cancer/Figures/PanCan_sig_BH_all.pdf', plot = PanCan_enrich_plot_1, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
# write.xlsx(CN_info_df, file = 'All_sample_signature.xlsx', sheetName = 'CN_group_1', append = TRUE, row.names = FALSE)
write_xlsx(PanCan_info_df, path = 'PanCan_group_1.xlsx')
