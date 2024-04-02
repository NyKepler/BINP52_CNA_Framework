# Title: analysis_signature.R
# Author: Guyuan TANG
# Date: 2024/1/5 - 2024/3/20

# Description: the script is used for analyzing the sample-by-component matrices from the workflow.

library(tidyverse)
library(ggpubr)
library(xlsx)
library(aplot)
library(ggsci)
library(writexl)

# set the working directory
setwd("E:/1Lund_Lectures/BINP52_MasterProject/Workflow/Draft/signatures")

# the required groups
groups <- c('ArchivalVS','MaNiLaVS','ffTumor','ffTissue','Endome','ffpe','Blood','ffPlasma')

# Generate the All_sample_signature.vcsv (only contains CN_sig and PanCan_sig)
## for panConusig, only ffTumor samples ran the analysis
# load the sample table
sample_df <- read.table('E:/1Lund_Lectures/BINP52_MasterProject/Workflow/Draft/results/solutions/solution_sample.tsv', sep = '\t', header = 1) %>% select('Sample', 'Patient', 'Type', 'Group')
# sample_df <- read.csv('E:/1Lund_Lectures/BINP52_MasterProject/Workflow/Draft/Snakemake2/config/solution_sample.tsv', sep = '\t') %>% select('Sample','Patient','Type','Group','Ploidy','Cellularity')
sample_sheet <- readxl::read_excel('MaNiLa_All_samplesheet_240311_groups.xlsx',sheet = 'Samples')
## combine the benign and HGSC information with the samples
sample_df$BH <- NA
sample_df$BH_type <- NA
## exclude the 'corpus' and 'Restaging op' sample (CS3_63 & CS3_72)
sample_df <- sample_df %>% filter(Sample != 'CS3_63' & Sample != 'CS3_72')
for (i in 1:nrow(sample_df)) {
  sampleID <- sample_df[i, 'Sample']
  sample_df[i,'Patient'] <- sample_sheet[which(sample_sheet$Library == sampleID), 'Patient']
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
  } else if (sig_type == 'PanCan' | sig_type=='Pan-Cancer') {
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
color_bar_prep <- function(stat_df, exposure_df, grouping_type, sig_type, sort=TRUE) {
  col_num <- ncol(exposure_df)
  ### create a dataframe to store the information, and rank by either groups or BH_groups
  exposure_out <- as.data.frame(matrix(nrow = 0,ncol = col_num))
  if (sort==FALSE | sort==F) {
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
  } else if (sort==TRUE | sort==T) {
    if (grouping_type=='groups') {
      groups <- c('ArchivalVS','MaNiLaVS','ffTumor','ffTissue','Endome','ffpe','Blood','ffPlasma')
      ### rank by groups
      x_lab_name <- 'Sample types'
      for (groupID in groups) {
        indices <- which(stat_df$group==groupID)
        sort_df <- exposure_df[indices,]
        # sort the samples by the first signatures
        if (sig_type == 'CN') {
          sort_df <- sort_df[order(-sort_df$s1),]
        } else if (sig_type=='Pan-Cancer' | sig_type=='PanCan') {
          sort_df <- sort_df[order(-sort_df$CX1),]
        } else if (sig_type=='panConusig') {
          sort_df <- sort_df[order(-sort_df$CN1),]
        }
        exposure_out <- rbind(exposure_out, sort_df)
        }
      } else if (grouping_type == 'BH_group') {
        BH_group <- c('Benign', 'RRSO','HGSC')
        ### rank by BH_group
        x_lab_name <- 'Groups'
        for (BH_ID in BH_group) {
          indices <- which(stat_df$BH==BH_ID)
          sort_df <- exposure_df[indices,]
          # sort the samples by the first signatures
          if (sig_type == 'CN') {
            sort_df <- sort_df[order(-sort_df$s1),]
          } else if (sig_type=='Pan-Cancer' | sig_type=='PanCan') {
            sort_df <- sort_df[order(-sort_df$CX1),]
          } else if (sig_type=='panConusig') {
            sort_df <- sort_df[order(-sort_df$CN1),]
          }
          exposure_out <- rbind(exposure_out, sort_df)
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
plot_exposure <- function(stat_df, exposure_df, grouping_type, sig_type, sort=TRUE) {
  col_num <- ncol(exposure_df)
  ### create a dataframe to store the information, and rank by either groups or BH_groups
  exposure_out <- as.data.frame(matrix(nrow = 0,ncol = col_num))
  if (sort==FALSE | sort==F) {
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
  } else if (sort==TRUE | sort==T) {
    if (grouping_type=='groups') {
      groups <- c('ArchivalVS','MaNiLaVS','ffTumor','ffTissue','Endome','ffpe','Blood','ffPlasma')
      ### rank by groups
      x_lab_name <- 'Sample types'
      for (groupID in groups) {
        indices <- which(stat_df$group==groupID)
        sort_df <- exposure_df[indices,]
        # sort the samples by the first signatures
        if (sig_type == 'CN') {
          sort_df <- sort_df[order(-sort_df$s1),]
        } else if (sig_type=='Pan-Cancer' | sig_type=='PanCan') {
          sort_df <- sort_df[order(-sort_df$CX1),]
        } else if (sig_type=='panConusig') {
          sort_df <- sort_df[order(-sort_df$CN1),]
        }
        exposure_out <- rbind(exposure_out, sort_df)
      }
    } else if (grouping_type == 'BH_group') {
      BH_group <- c('Benign', 'RRSO','HGSC')
      ### rank by BH_group
      x_lab_name <- 'Groups'
      for (BH_ID in BH_group) {
        indices <- which(stat_df$BH==BH_ID)
        sort_df <- exposure_df[indices,]
        # sort the samples by the first signatures
        if (sig_type == 'CN') {
          sort_df <- sort_df[order(-sort_df$s1),]
        } else if (sig_type=='Pan-Cancer' | sig_type=='PanCan') {
          sort_df <- sort_df[order(-sort_df$CX1),]
        } else if (sig_type=='panConusig') {
          sort_df <- sort_df[order(-sort_df$CN1),]
        }
        exposure_out <- rbind(exposure_out, sort_df)
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


# define a function to plot the VS sample exposure
VS_exposure_plot <- function(VS_exposure_df, sig_type, VS_stat_df) {
  rownames(VS_exposure_df) <- 1:nrow(VS_exposure_df)
  # prepare the sample ranking
  sample_rank <- VS_exposure_df$sample
  # prepare the dataframe
  exposure_out <- pivot_longer(VS_exposure_df,!sample, names_to = 'signatures',values_to = 'exposure')
  exposure_out$sample <- factor(exposure_out$sample, levels = sample_rank)
  # set the ranking of the signatures
  if (sig_type == 'Pan-Cancer') {
    exposure_out$signatures <- factor(exposure_out$signatures, levels = c('CX1', 'CX2', 'CX3', 'CX4', 'CX5', 'CX6','CX7','CX8','CX9','CX10','CX11','CX12','CX13','CX14','CX15','CX16','CX17'))
  } else if (sig_type == 'CN') {
    exposure_out$signatures <- factor(exposure_out$signatures, levels = c('s1','s2','s3','s4','s5','s6','s7'))
  } else if (sig_type == 'panConusig') {
    exposure_out$signatures <- factor(exposure_out$signatures, levels = c('CN1','CN2','CN3','CN4','CN5','CN6','CN7','CN8','CN9','CN10','CN11','CN12','CN13','CN14','CN15','CN16','CN17','CN18','CN19','CN20','CN21','CN22','CN23','CN24','CN25'))
  }
  # add the diagnostic timepoints
  exposure_out$diagnostic_time <- NA
  for (sampleID in sample_rank) {
    exposure_out[which(exposure_out$sample==sampleID),'diagnostic_time'] <- VS_stat_df[which(VS_stat_df$sample==sampleID),'Diagnostic_time']
  }
  exposure_out$diagnostic_time <- factor(exposure_out$diagnostic_time, levels = c('prediagnostic_6+', 'prediagnostic_0-6', 'diagnostic'))
  exposure_plot <- ggplot(data = exposure_out, aes(x=sample, y=exposure, fill=signatures)) +
    geom_bar(stat = 'identity', width = 0.8, position = position_stack(reverse = T)) +
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


# define a function to draw reference exposure plots
ref_exposure_plot <- function(ref_exposure_df, sig_type, ref_stat_df) {
  rownames(ref_exposure_df) <- 1:nrow(ref_exposure_df)
  # prepare the sample ranking
  sample_rank <- ref_exposure_df$sample
  # prepare the dataframe
  exposure_out <- pivot_longer(ref_exposure_df,!sample, names_to = 'signatures',values_to = 'exposure')
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
    geom_bar(stat = 'identity', width = 0.8, position = position_stack(reverse = T)) +
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

## 1. differences in ploidy and cellularity
for (sampleID in sample_df$Sample) {
  if (!(sampleID %in% CN_stat$sample)) {
    sample_df[which(sample_df$Sample==sampleID),'enrich_CN_1'] <- NA
  } else if (sampleID %in% CN_stat$sample) {
    sample_df[which(sample_df$Sample==sampleID),'enrich_CN_1'] <- CN_stat[which(CN_stat$sample==sampleID),'enrich_CN_1']
  }
}
sample_df$Group <- factor(sample_df$Group, levels = groups)
### plotting
ploidy_plot <-  ggboxplot(sample_df, x='Group', y='Ploidy', 
                          color = 'black', fill = 'enrich_CN_1', palette = 'Paired', 
                          xlab = 'Group', ylab = 'Ploidy') +
  theme_classic()
ploidy_plot <- ploidy_plot + labs(fill = 'CN signatures')
ggsave('Brenton/Figures/CN_ploidy_all.pdf', plot = ploidy_plot, dpi = 600, width = 10, height = 7, units = 'in')
cellularity_plot <-  ggboxplot(sample_df, x='Group', y='Cellularity', 
                               color = 'black', fill = 'enrich_CN_1', palette = 'Paired', 
                               xlab = 'Group', ylab = 'Cellularity') +
  theme_classic()
cellularity_plot <- cellularity_plot + labs(fill = 'CN signatures')
ggsave('Brenton/Figures/CN_cellularity_all.pdf', plot = cellularity_plot, dpi = 600, width = 10, height = 7, units = 'in') 
### stats analysis
sub_df <- sample_df %>% filter(Group=='ffTumor' & !is.na(enrich_PanCan_1))
# ploidy
kruskal.test(sub_df$Ploidy~sub_df$enrich_CN_1, data = sub_df) # p = 0.26
# cellularity
kruskal.test(sub_df$Cellularity~sub_df$enrich_CN_1, data = sub_df) # p = 0.09

## 2. distribution of the delta-value
CN_dist_delta_plot <- dis_plot_delta(CN_stat)
ggsave(filename = 'Brenton/Figures/delta_val_dist.pdf', plot = CN_dist_delta_plot, dpi = 600, width = 10, height = 7, units = 'in')
quantile(CN_stat$delta_val, 0.1) # 0.079 top 10%
# delta value cut-off at 0.25
for (groupID in groups) {
  print(groupID)
  group_df <- CN_stat[which(CN_stat$group==groupID),]
  sample_sum <- nrow(group_df)
  print(sample_sum)
  group_df <- group_df[which(group_df$delta_val<0.25),]
  sample_size <- nrow(group_df)
  print(sample_size)
  print(100*sample_size/sample_sum)
}
CN_stat_25 <- CN_stat[which(CN_stat$delta_val<0.25),]
groups2 <- c('ArchivalVS','MaNiLaVS','ffTissue','Endome','ffpe','Blood','ffPlasma')
#### top 2 enrich signature
CN_info_df_25 <- sig_stat_enrich(stat_df = CN_stat_25, enrich_num = 2, sig_type = 'CN')
# plot the graph
CN_enrich_plot_2 <- ggplot(data = CN_info_df_25, aes(x=BH, y=percentage, fill=enrich_CN_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='CN signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group2))
CN_enrich_plot_2 <- CN_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Brenton/Figures/CN_sig_BH_025_2.pdf', plot = CN_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
write_xlsx(CN_info_df_25, path = 'CN_group_025_2.xlsx')

# delta value cut-off at 0.1
for (groupID in groups) {
  print(groupID)
  group_df <- CN_stat[which(CN_stat$group==groupID),]
  sample_sum <- nrow(group_df)
  print(sample_sum)
  group_df <- group_df[which(group_df$delta_val<0.1),]
  sample_size <- nrow(group_df)
  print(sample_size)
  print(100*sample_size/sample_sum)
}
CN_stat_01 <- CN_stat[which(CN_stat$delta_val<0.1),]
groups2 <- c('ArchivalVS','MaNiLaVS','ffTissue','Endome','ffpe','Blood','ffPlasma')
#### top 2 enrich signature
CN_info_df_01 <- sig_stat_enrich(stat_df = CN_stat_01, enrich_num = 2, sig_type = 'CN')
# plot the graph
CN_enrich_plot_2 <- ggplot(data = CN_info_df_01, aes(x=BH, y=percentage, fill=enrich_CN_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='CN signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group2))
CN_enrich_plot_2 <- CN_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Brenton/Figures/CN_sig_BH_01_2.pdf', plot = CN_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
write_xlsx(CN_info_df_01, path = 'CN_group_01_2.xlsx')

# after cut-off on delta values
## 0.1 for ffTumor and ffpe
## 0.25 for all other groups
CN_stat_01 <- CN_stat_01 %>% filter(group=='ffTumor' | group=='ffpe')
CN_stat_25 <- CN_stat_25 %>% filter(group!='ffTumor' & group!='ffpe')
CN_stat_cut <- rbind(CN_stat_01, CN_stat_25)
#### top 2 enrich signature
CN_info_df_cut <- sig_stat_enrich(stat_df = CN_stat_cut, enrich_num = 2, sig_type = 'CN')
# plot the graph
CN_enrich_plot_2 <- ggplot(data = CN_info_df_cut, aes(x=BH, y=percentage, fill=enrich_CN_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='CN signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group2))
CN_enrich_plot_2 <- CN_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Brenton/Figures/CN_sig_BH_cut_2.pdf', plot = CN_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
write_xlsx(CN_info_df_cut, path = 'CN_group_cut_2.xlsx')


## 3. exposure distribution
### rank by group
CN_color_bar_group <- color_bar_prep(stat_df = CN_stat, exposure_df = CN_exposure, grouping_type = 'groups', sig_type = 'CN', sort = TRUE)
CN_exposure_plot_group <- plot_exposure(stat_df = CN_stat, exposure_df = CN_exposure, grouping_type = 'groups', sig_type = 'CN', sort = TRUE)
CN_exposure_plot_group <- CN_exposure_plot_group %>% insert_bottom(CN_color_bar_group, height = 0.03)
ggsave(filename = 'Brenton/Figures/exposure_plot_groups_sort.pdf', plot = CN_exposure_plot_group, dpi = 600, width = 30, height = 7, units = 'in')
### rank by BH
CN_color_bar_BH <- color_bar_prep(stat_df = CN_stat, exposure_df = CN_exposure, grouping_type = 'BH_group', sig_type = 'CN', sort = TRUE)
CN_exposure_plot_BH <- plot_exposure(stat_df = CN_stat, exposure_df = CN_exposure, grouping_type = 'BH_group', sig_type = 'CN', sort = TRUE)
CN_exposure_plot_BH <- CN_exposure_plot_BH %>% insert_bottom(CN_color_bar_BH, height = 0.03)
ggsave(filename = 'Brenton/Figures/exposure_plot_BH_sort.pdf', plot = CN_exposure_plot_BH, dpi = 600, width = 30, height = 7, units = 'in')
### details for the ArchivalVS group
CN_archival_stat <- CN_stat %>% filter(group=='ArchivalVS')
archival_sample <- CN_archival_stat$sample
CN_archival_mat <- CN_mat %>% filter(sample %in% archival_sample)
CN_archival_exposure <- sig_exposure(SSmatrix = CN_archival_mat,sig_type = 'CN')
CN_color_bar_BH_AVS <- color_bar_prep(stat_df = CN_archival_stat, exposure_df = CN_archival_exposure, grouping_type = 'BH_group', sig_type = 'CN')
CN_exposure_plot_BH_AVS <- plot_exposure(stat_df = CN_archival_stat, exposure_df = CN_archival_exposure, grouping_type = 'BH_group', sig_type = 'CN')
CN_exposure_plot_BH_AVS <- CN_exposure_plot_BH_AVS %>% insert_bottom(CN_color_bar_BH_AVS, height = 0.03)
ggsave(filename = 'Brenton/Figures/exposure_plot_BH_AVS.pdf', plot = CN_exposure_plot_BH, dpi = 600, width = 30, height = 7, units = 'in')

## 4. CN signatures enrichment
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

#### top 2 enrich signature
CN_info_df <- sig_stat_enrich(stat_df = CN_stat, enrich_num = 2, sig_type = 'CN')
# plot the graph
CN_enrich_plot_2 <- ggplot(data = CN_info_df, aes(x=BH, y=percentage, fill=enrich_CN_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='CN signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group2))
CN_enrich_plot_2 <- CN_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Brenton/Figures/CN_sig_BH_all_2.pdf', plot = CN_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
# write.xlsx(CN_info_df, file = 'All_sample_signature.xlsx', sheetName = 'CN_group_2', append = TRUE, row.names = FALSE)
write_xlsx(CN_info_df, path = 'CN_group_2.xlsx')


## 5. VS samples time point exposure
VS_sheet <- read.xlsx(file = 'MaNiLa_All_samplesheet_240311_groups.xlsx', sheetName = 'VS_samples')
VS_stat <- CN_stat[which(CN_stat$group=='ArchivalVS' | CN_stat$group=='MaNiLaVS'),]
for (sampleID in VS_stat$sample) {
  VS_stat[which(VS_stat$sample==sampleID), 'Diagnostic_time'] <- VS_sheet[which(VS_sheet$Library==sampleID), 'Diagnostic_time']
}

### extract the similarity matrix
indices <- rownames(VS_stat)
VS_mat <- CN_mat[indices,]
VS_exposure <- sig_exposure(SSmatrix = VS_mat,sig_type = 'CN')
VS_color_bar_BH <- color_bar_prep(stat_df = VS_stat, exposure_df = VS_exposure, grouping_type = 'BH_group', sig_type = 'CN')
VS_exposure_plot_BH <- plot_exposure(stat_df = VS_stat, exposure_df = VS_exposure, grouping_type = 'BH_group', sig_type = 'CN')
VS_exposure_plot_BH <- VS_exposure_plot_BH %>% insert_bottom(VS_color_bar_BH, height = 0.03)
ggsave('Brenton/Figures/VS_exposure.pdf', plot = VS_exposure_plot_BH, dpi = 600, width = 10, height = 7, units = 'in')

### for Benign VS samples
Benign_VS_stat <- VS_stat[which(VS_stat$BH=='Benign'),]
unique(Benign_VS_stat$Diagnostic_time) # only "diagnostic"
indices <- rownames(Benign_VS_stat)
Benign_VS_mat <- VS_mat[indices,]
#### exposure
Benign_VS_exposure <- sig_exposure(SSmatrix = Benign_VS_mat,sig_type = 'CN')
#### sort by s1
Benign_VS_exposure <- Benign_VS_exposure[order(-Benign_VS_exposure$s1),]
Benign_VS_exposure <- select(Benign_VS_exposure,!SS_sum)
#### plotting
Benign_VS_plot <- VS_exposure_plot(VS_exposure_df = Benign_VS_exposure, sig_type = 'CN', VS_stat_df = Benign_VS_stat) + labs(title = 'diagnostic') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('Brenton/Figures/Benign_VS_exposure.pdf', plot = Benign_VS_plot, dpi = 600, width = 7, height = 5, units = 'in')

### for RRSO samples
RRSO_VS_stat <- VS_stat[which(VS_stat$BH=='RRSO'),]
unique(RRSO_VS_stat$Diagnostic_time) # "prediagnostic_6+", "diagnostic"
indices <- rownames(RRSO_VS_stat)
RRSO_VS_mat <- VS_mat[indices,]
#### exposure
RRSO_VS_exposure_in <- sig_exposure(SSmatrix = RRSO_VS_mat,sig_type = 'CN')
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
pre_6_RRSO_exposure_plot <- VS_exposure_plot(VS_exposure_df = pre_6_RRSO_exposure, sig_type = 'CN', VS_stat_df = RRSO_VS_stat) + labs(title = 'prediagnostic > 6m') + theme(plot.title = element_text(size = 10, hjust = 0.5))
diag_RRSO_exposure_plot <- VS_exposure_plot(VS_exposure_df = diag_RRSO_exposure, sig_type = 'CN', VS_stat_df = RRSO_VS_stat)  + labs(title = 'diagnostic') + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(size = 10, hjust = 0.5))
RRSO_VS_plot <- ggarrange(pre_6_RRSO_exposure_plot, diag_RRSO_exposure_plot, ncol = 2, common.legend = 1, legend.grob = get_legend(pre_6_RRSO_exposure_plot), legend = 'right')
ggsave('Brenton/Figures/RRSO_VS_exposure.pdf', plot = RRSO_VS_plot, dpi = 600, width = 7, height = 5, units = 'in')

### for HGSC samples
HGSC_VS_stat <- VS_stat[which(VS_stat$BH=='HGSC'),]
unique(HGSC_VS_stat$Diagnostic_time) # "prediagnostic_0-6", "prediagnostic_6+", "diagnostic"
indices <- rownames(HGSC_VS_stat)
HGSC_VS_mat <- VS_mat[indices,]
#### exposure
HGSC_VS_exposure_in <- sig_exposure(SSmatrix = HGSC_VS_mat,sig_type = 'CN')
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
pre_6_HGSC_exposure_plot <- VS_exposure_plot(VS_exposure_df = pre_6_HGSC_exposure, sig_type = 'CN', VS_stat_df = HGSC_VS_stat) + labs(title = 'prediagnostic > 6m') + theme(plot.title = element_text(size = 10, hjust = 0.5))
#2
pre_0_HGSC_exposure_plot <- VS_exposure_plot(VS_exposure_df = pre_0_HGSC_exposure, sig_type = 'CN', VS_stat_df = HGSC_VS_stat)  + labs(title = 'prediagnostic 0-6m') + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(size = 10, hjust = 0.5))
#3
diag_HGSC_exposure_plot <- VS_exposure_plot(VS_exposure_df = diag_HGSC_exposure, sig_type = 'CN', VS_stat_df = HGSC_VS_stat)  + labs(title = 'diagnostic') + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(size = 10, hjust = 0.5))
#combine
HGSC_VS_plot <- ggarrange(pre_6_HGSC_exposure_plot, pre_0_HGSC_exposure_plot, diag_HGSC_exposure_plot, ncol = 3, common.legend = 1, legend.grob = get_legend(pre_0_HGSC_exposure_plot), legend = 'right')
ggsave('Brenton/Figures/HGSC_VS_exposure.pdf', plot = HGSC_VS_plot, dpi = 600, width = 10, height = 5, units = 'in')


## 6. reference exposure plots
#### HGSC ffTumor 
HGSC_tumor_indices <- which(CN_stat$group=='ffTumor' & CN_stat$BH=='HGSC')
HGSC_tumor_mat <- CN_mat[HGSC_tumor_indices,]
HGSC_tumor_stat <- CN_stat[HGSC_tumor_indices,]
HGSC_tumor_exposure <- sig_exposure(SSmatrix = HGSC_tumor_mat, sig_type = 'CN')
#### sort by s1
HGSC_tumor_exposure <- HGSC_tumor_exposure[order(-HGSC_tumor_exposure$s1),]
HGSC_tumor_exposure <- select(HGSC_tumor_exposure,!SS_sum)
#### plotting
HGSC_tumor_exposure_plot <- ref_exposure_plot(ref_exposure_df = HGSC_tumor_exposure, sig_type = 'CN', ref_stat_df = HGSC_tumor_stat) + labs(title = 'HGSC ffTumor') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('Brenton/Figures/ref_HGSC_exposure.pdf', plot = HGSC_tumor_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### Benign ffTissue
Benign_tissue_indices <- which(CN_stat$group=='ffTissue' & CN_stat$BH=='Benign')
Benign_tissue_mat <- CN_mat[Benign_tissue_indices,]
Benign_tissue_stat <- CN_stat[Benign_tissue_indices,]
Benign_tissue_exposure <- sig_exposure(SSmatrix = Benign_tissue_mat, sig_type = 'CN')
#### sort by s1
Benign_tissue_exposure <- Benign_tissue_exposure[order(-Benign_tissue_exposure$s1),]
Benign_tissue_exposure <- select(Benign_tissue_exposure, !SS_sum)
#### plotting
Benign_tissue_exposure_plot <- ref_exposure_plot(ref_exposure_df = Benign_tissue_exposure, sig_type = 'CN', ref_stat_df = Benign_tissue_stat) + labs(title = 'Benign ffTissue') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('Brenton/Figures/ref_Benign_exposure.pdf', plot = Benign_tissue_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### RRSO ffTissue+ffpe
RRSO_ref_indices <- which(CN_stat$BH=='RRSO' & (CN_stat$group=='ffTissue' | CN_stat$group=='ffpe'))
RRSO_ref_mat <- CN_mat[RRSO_ref_indices,]
RRSO_ref_stat <- CN_stat[RRSO_ref_indices,]
RRSO_ref_exposure <- sig_exposure(SSmatrix = RRSO_ref_mat, sig_type = 'CN')
#### sort by s1
RRSO_ref_exposure <- RRSO_ref_exposure[order(-RRSO_ref_exposure$s1),]
RRSO_ref_exposure <- select(RRSO_ref_exposure, !SS_sum)
#### plotting
RRSO_ref_exposure_plot <- ref_exposure_plot(ref_exposure_df = RRSO_ref_exposure, sig_type = 'CN', ref_stat_df = RRSO_ref_stat) + labs(title = 'RRSO ffTissue&ffpe') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('Brenton/Figures/ref_RRSO_exposure.pdf', plot = RRSO_ref_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### combined
RRSO_ref_exposure_plot_1 <- RRSO_ref_exposure_plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
HGSC_tumor_exposure_plot_1 <- HGSC_tumor_exposure_plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
CN_ref_exposure_plot <- ggarrange(Benign_tissue_exposure_plot, RRSO_ref_exposure_plot_1, HGSC_tumor_exposure_plot_1, ncol = 3, common.legend = 1, legend.grob = get_legend(Benign_tissue_exposure_plot), legend = 'right')
ggsave('Brenton/Figures/ref_exposure_all.pdf', plot = CN_ref_exposure_plot, dpi = 600, width = 10, height = 5, units = 'in')



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

## 1. differences in ploidy and cellularity
for (sampleID in sample_df$Sample) {
  if (!(sampleID %in% PanCan_stat$sample)) {
    sample_df[which(sample_df$Sample==sampleID),'enrich_PanCan_1'] <- NA
  } else if (sampleID %in% PanCan_stat$sample) {
    sample_df[which(sample_df$Sample==sampleID),'enrich_PanCan_1'] <- PanCan_stat[which(PanCan_stat$sample==sampleID),'enrich_PanCan_1']
  }
}
PanCan_sample_df <- sample_df %>% filter(!is.na(enrich_PanCan_1))
PanCan_sample_df$Group <- factor(PanCan_sample_df$Group, levels = groups)
PanCan_sample_df$enrich_PanCan_1 <- factor(PanCan_sample_df$enrich_PanCan_1, levels = c('CX1','CX2','CX3','CX5','CX15','CX17'))
### plotting
ploidy_plot <-  ggboxplot(PanCan_sample_df, x='Group', y='Ploidy', 
                          color = 'black', fill = 'enrich_PanCan_1', palette = 'Paired', 
                          xlab = 'Group', ylab = 'Ploidy') +
  theme_classic()
ploidy_plot <- ploidy_plot + labs(fill = 'pan-cancer signatures')
ggsave('Pan-Cancer/Figures/PanCan_ploidy_all.pdf', plot = ploidy_plot, dpi = 600, width = 10, height = 7, units = 'in')
cellularity_plot <-  ggboxplot(PanCan_sample_df, x='Group', y='Cellularity', 
                               color = 'black', fill = 'enrich_PanCan_1', palette = 'Paired', 
                               xlab = 'Group', ylab = 'Cellularity') +
  theme_classic()
cellularity_plot <- cellularity_plot + labs(fill = 'pan-cancer signatures')
ggsave('Pan-Cancer/Figures/PanCan_cellularity_all.pdf', plot = cellularity_plot, dpi = 600, width = 10, height = 7, units = 'in') 
### stats analysis
sub_df <- PanCan_sample_df %>% filter(Group=='ffTumor')
# ploidy
kruskal.test(sub_df$Ploidy~sub_df$enrich_PanCan_1, data = sub_df) # p = 0.84
# cellularity
kruskal.test(sub_df$Cellularity~sub_df$enrich_PanCan_1, data = sub_df) # p = 0.08

## 2. distribution of the delta-value
PanCan_dist_delta_plot <- dis_plot_delta(PanCan_stat)
ggsave(filename = 'Pan-Cancer/Figures/delta_val_dist.pdf', plot = PanCan_dist_delta_plot, dpi = 600, width = 10, height = 7, units = 'in')
quantile(PanCan_stat$delta_val, 0.1) # 0.016 top 10%
# delta value cut-off at 0.25
for (groupID in groups) {
  print(groupID)
  group_df <- PanCan_stat[which(PanCan_stat$group==groupID),]
  sample_sum <- nrow(group_df)
  print(sample_sum)
  group_df <- group_df[which(group_df$delta_val<0.25),]
  sample_size <- nrow(group_df)
  print(sample_size)
  print(100*sample_size/sample_sum)
}
PanCan_stat_25 <- PanCan_stat[which(PanCan_stat$delta_val<0.25),]
groups2 <- c('ArchivalVS','MaNiLaVS','ffTissue','Endome','ffpe','Blood','ffPlasma')
#### top 2 enrich signature
PanCan_info_df_25 <- sig_stat_enrich(stat_df = PanCan_stat_25, enrich_num = 2, sig_type = 'PanCan')
PanCan_info_df_25$enrich_PanCan_2 <- factor(PanCan_info_df_25$enrich_PanCan_2, levels = c('CX1','CX2','CX3','CX5','CX10','CX15','CX17'))
# plot the graph
PanCan_enrich_plot_2 <- ggplot(data = PanCan_info_df_25, aes(x=BH, y=percentage, fill=enrich_PanCan_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='Pan-Cancer signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group2))
PanCan_enrich_plot_2 <- PanCan_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Pan-Cancer/Figures/PanCan_sig_BH_025_2.pdf', plot = PanCan_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
write_xlsx(PanCan_info_df_25, path = 'PanCan_group_025_2.xlsx')

# delta value cut-off at 0.1
for (groupID in groups) {
  print(groupID)
  group_df <- PanCan_stat[which(PanCan_stat$group==groupID),]
  sample_sum <- nrow(group_df)
  print(sample_sum)
  group_df <- group_df[which(group_df$delta_val<0.1),]
  sample_size <- nrow(group_df)
  print(sample_size)
  print(100*sample_size/sample_sum)
}
PanCan_stat_01 <- PanCan_stat[which(PanCan_stat$delta_val<0.1),]
groups2 <- c('ArchivalVS','MaNiLaVS','ffTissue','Endome','ffpe','Blood','ffPlasma')
#### top 2 enrich signature
PanCan_info_df_01 <- sig_stat_enrich(stat_df = PanCan_stat_01, enrich_num = 2, sig_type = 'PanCan')
PanCan_info_df_01$enrich_PanCan_2 <- factor(PanCan_info_df_01$enrich_PanCan_2, levels = c('CX1','CX2','CX3','CX5','CX10','CX15','CX17'))
# plot the graph
PanCan_enrich_plot_2 <- ggplot(data = PanCan_info_df_01, aes(x=BH, y=percentage, fill=enrich_PanCan_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='Pan-Cancer signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group2))
PanCan_enrich_plot_2 <- PanCan_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Pan-Cancer/Figures/PanCan_sig_BH_01_2.pdf', plot = PanCan_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
write_xlsx(PanCan_info_df_01, path = 'PanCan_group_01_2.xlsx')

# after adjusting the delta-value cutoff
## 0.25 for MaNiLaVS, ffTissue, Endome, and Blood
## 0.1 for all other groups
PanCan_stat_01 <- PanCan_stat_01 %>% filter(group!='MaNiLaVS' & group!='ffTissue' & group!='Endome' & group!='Blood')
PanCan_stat_25 <- PanCan_stat_25 %>% filter(group=='MaNiLaVS' | group=='ffTissue' | group=='Endome' | group=='Blood')
PanCan_stat_cut <- rbind(PanCan_stat_01, PanCan_stat_25)
#### top 2 enrich signature
PanCan_info_df_cut <- sig_stat_enrich(stat_df = PanCan_stat_cut, enrich_num = 2, sig_type = 'PanCan')
PanCan_info_df_cut$enrich_PanCan_2 <- factor(PanCan_info_df_cut$enrich_PanCan_2, levels = c('CX1','CX2','CX3','CX5','CX10','CX15','CX17'))
# plot the graph
PanCan_enrich_plot_2 <- ggplot(data = PanCan_info_df_cut, aes(x=BH, y=percentage, fill=enrich_PanCan_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='Pan-Cancer signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group2))
PanCan_enrich_plot_2 <- PanCan_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Pan-Cancer/Figures/PanCan_sig_BH_cut_2.pdf', plot = PanCan_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
write_xlsx(PanCan_info_df_cut, path = 'PanCan_group_cut_2.xlsx')


## 3. exposure distribution
### rank by group
PanCan_color_bar_group <- color_bar_prep(stat_df = PanCan_stat, exposure_df = PanCan_exposure, grouping_type = 'groups', sig_type = 'Pan-Cancer', sort = TRUE)
PanCan_exposure_plot_group <- plot_exposure(stat_df = PanCan_stat, exposure_df = PanCan_exposure, grouping_type = 'groups', sig_type = 'Pan-Cancer', sort = TRUE)
PanCan_exposure_plot_group <- PanCan_exposure_plot_group %>% insert_bottom(PanCan_color_bar_group, height = 0.03)
ggsave(filename = 'Pan-Cancer/Figures/exposure_plot_groups_sort.pdf', plot = PanCan_exposure_plot_group, dpi = 600, width = 30, height = 7, units = 'in')
### rank by BH
PanCan_color_bar_BH <- color_bar_prep(stat_df = PanCan_stat, exposure_df = PanCan_exposure, grouping_type = 'BH_group', sig_type = 'Pan-Cancer', sort = TRUE)
PanCan_exposure_plot_BH <- plot_exposure(stat_df = PanCan_stat, exposure_df = PanCan_exposure, grouping_type = 'BH_group', sig_type = 'Pan-Cancer', sort = TRUE)
PanCan_exposure_plot_BH <- PanCan_exposure_plot_BH %>% insert_bottom(PanCan_color_bar_BH, height = 0.03)
ggsave(filename = 'Pan-Cancer/Figures/exposure_plot_BH_sort.pdf', plot = PanCan_exposure_plot_BH, dpi = 600, width = 30, height = 7, units = 'in')
### details for the ArchivalVS group
PanCan_archival_stat <- PanCan_stat %>% filter(group=='ArchivalVS')
archival_sample <- PanCan_archival_stat$sample
PanCan_archival_mat <- PanCan_mat %>% filter(sample %in% archival_sample)
PanCan_archival_exposure <- sig_exposure(SSmatrix = PanCan_archival_mat,sig_type = 'PanCan')
PanCan_color_bar_BH_AVS <- color_bar_prep(stat_df = PanCan_archival_stat, exposure_df = PanCan_archival_exposure, grouping_type = 'BH_group', sig_type = 'Pan-Cancer')
PanCan_exposure_plot_BH_AVS <- plot_exposure(stat_df = PanCan_archival_stat, exposure_df = PanCan_archival_exposure, grouping_type = 'BH_group', sig_type = 'Pan-Cancer')
PanCan_exposure_plot_BH_AVS <- PanCan_exposure_plot_BH_AVS %>% insert_bottom(PanCan_color_bar_BH_AVS, height = 0.03)
ggsave(filename = 'Pan-Cancer/Figures/exposure_plot_BH_AVS.pdf', plot = PanCan_exposure_plot_BH_AVS, dpi = 600, width = 30, height = 7, units = 'in')

## 4. PanCan signatures enrichment
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
write_xlsx(PanCan_info_df, path = 'PanCan_group_1.xlsx')

#### top 2 enrich signature
PanCan_info_df <- sig_stat_enrich(stat_df = PanCan_stat, enrich_num = 2, sig_type = 'PanCan')
PanCan_info_df$enrich_PanCan_2 <- factor(PanCan_info_df$enrich_PanCan_2, levels = c('CX1', 'CX2', 'CX3', 'CX5', 'CX10', 'CX15', 'CX17'))
# plot the graph
PanCan_enrich_plot_2 <- ggplot(data = PanCan_info_df, aes(x=BH, y=percentage, fill=enrich_PanCan_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='Pan-Cancer signatures') +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() + facet_grid(cols = vars(group2))
PanCan_enrich_plot_2 <- PanCan_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('Pan-Cancer/Figures/PanCan_sig_BH_all_2.pdf', plot = PanCan_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
write_xlsx(PanCan_info_df, path = 'PanCan_group_2.xlsx')


## 5. VS samples time point exposure
VS_sheet <- read.xlsx(file = 'MaNiLa_All_samplesheet_240311_groups.xlsx', sheetName = 'VS_samples')
VS_stat <- PanCan_stat[which(PanCan_stat$group=='ArchivalVS' | PanCan_stat$group=='MaNiLaVS'),]
for (sampleID in VS_stat$sample) {
  VS_stat[which(VS_stat$sample==sampleID), 'Diagnostic_time'] <- VS_sheet[which(VS_sheet$Library==sampleID), 'Diagnostic_time']
}

### extract the similarity matrix
indices <- rownames(VS_stat)
VS_mat <- PanCan_mat[indices,]
VS_exposure <- sig_exposure(SSmatrix = VS_mat,sig_type = 'PanCan')
VS_color_bar_BH <- color_bar_prep(stat_df = VS_stat, exposure_df = VS_exposure, grouping_type = 'BH_group', sig_type = 'Pan-Cancer')
VS_exposure_plot_BH <- plot_exposure(stat_df = VS_stat, exposure_df = VS_exposure, grouping_type = 'BH_group', sig_type = 'Pan-Cancer')
VS_exposure_plot_BH <- VS_exposure_plot_BH %>% insert_bottom(VS_color_bar_BH, height = 0.03)
ggsave('Pan-Cancer/Figures/VS_exposure.pdf', plot = VS_exposure_plot_BH, dpi = 600, width = 10, height = 7, units = 'in')

### for Benign VS samples
Benign_VS_stat <- VS_stat[which(VS_stat$BH=='Benign'),]
unique(Benign_VS_stat$Diagnostic_time) # only "diagnostic"
indices <- rownames(Benign_VS_stat)
Benign_VS_mat <- VS_mat[indices,]
#### exposure
Benign_VS_exposure <- sig_exposure(SSmatrix = Benign_VS_mat,sig_type = 'PanCan')
#### sort by s1
Benign_VS_exposure <- Benign_VS_exposure[order(-Benign_VS_exposure$CX1),]
Benign_VS_exposure <- select(Benign_VS_exposure,!SS_sum)
#### plotting
Benign_VS_plot <- VS_exposure_plot(VS_exposure_df = Benign_VS_exposure, sig_type = 'Pan-Cancer', VS_stat_df = Benign_VS_stat) + labs(title = 'diagnostic') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('Pan-Cancer/Figures/Benign_VS_exposure.pdf', plot = Benign_VS_plot, dpi = 600, width = 7, height = 5, units = 'in')

### for RRSO samples
RRSO_VS_stat <- VS_stat[which(VS_stat$BH=='RRSO'),]
unique(RRSO_VS_stat$Diagnostic_time) # "diagnostic"
indices <- rownames(RRSO_VS_stat)
RRSO_VS_mat <- VS_mat[indices,]
#### exposure
RRSO_VS_exposure <- sig_exposure(SSmatrix = RRSO_VS_mat,sig_type = 'PanCan')
#### sort by CX1
RRSO_VS_exposure <- RRSO_VS_exposure[order(-RRSO_VS_exposure$CX1),]
RRSO_VS_exposure <- select(RRSO_VS_exposure,!SS_sum)
#### plotting
RRSO_VS_plot <- VS_exposure_plot(VS_exposure_df = RRSO_VS_exposure, sig_type = 'Pan-Cancer', VS_stat_df = RRSO_VS_stat) + labs(title = 'diagnostic') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('Pan-Cancer/Figures/RRSO_VS_exposure.pdf', plot = RRSO_VS_plot, dpi = 600, width = 7, height = 5, units = 'in')

### for HGSC samples
HGSC_VS_stat <- VS_stat[which(VS_stat$BH=='HGSC'),]
unique(HGSC_VS_stat$Diagnostic_time) # "prediagnostic_0-6", "prediagnostic_6+", "diagnostic"
indices <- rownames(HGSC_VS_stat)
HGSC_VS_mat <- VS_mat[indices,]
#### exposure
HGSC_VS_exposure_in <- sig_exposure(SSmatrix = HGSC_VS_mat,sig_type = 'PanCan')
#### sort by s1
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
ggsave('Pan-Cancer/Figures/HGSC_VS_exposure.pdf', plot = HGSC_VS_plot, dpi = 600, width = 10, height = 5, units = 'in')

## 6. reference exposure plots
#### HGSC ffTumor 
HGSC_tumor_indices <- which(PanCan_stat$group=='ffTumor' & PanCan_stat$BH=='HGSC')
HGSC_tumor_mat <- PanCan_mat[HGSC_tumor_indices,]
HGSC_tumor_stat <- PanCan_stat[HGSC_tumor_indices,]
HGSC_tumor_exposure <- sig_exposure(SSmatrix = HGSC_tumor_mat, sig_type = 'PanCan')
#### sort by CX1
HGSC_tumor_exposure <- HGSC_tumor_exposure[order(-HGSC_tumor_exposure$CX1),]
HGSC_tumor_exposure <- select(HGSC_tumor_exposure,!SS_sum)
#### plotting
HGSC_tumor_exposure_plot <- ref_exposure_plot(ref_exposure_df = HGSC_tumor_exposure, sig_type = 'Pan-Cancer', ref_stat_df = HGSC_tumor_stat) + labs(title = 'HGSC ffTumor') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('Pan-Cancer/Figures/ref_HGSC_exposure.pdf', plot = HGSC_tumor_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### Benign ffTissue
Benign_tissue_indices <- which(PanCan_stat$group=='ffTissue' & PanCan_stat$BH=='Benign')
Benign_tissue_mat <- PanCan_mat[Benign_tissue_indices,]
Benign_tissue_stat <- PanCan_stat[Benign_tissue_indices,]
Benign_tissue_exposure <- sig_exposure(SSmatrix = Benign_tissue_mat, sig_type = 'PanCan')
#### sort by CX1
Benign_tissue_exposure <- Benign_tissue_exposure[order(-Benign_tissue_exposure$CX1),]
Benign_tissue_exposure <- select(Benign_tissue_exposure, !SS_sum)
#### plotting
Benign_tissue_exposure_plot <- ref_exposure_plot(ref_exposure_df = Benign_tissue_exposure, sig_type = 'Pan-Cancer', ref_stat_df = Benign_tissue_stat) + labs(title = 'Benign ffTissue') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('Pan-Cancer/Figures/ref_Benign_exposure.pdf', plot = Benign_tissue_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### RRSO ffTissue+ffpe
RRSO_ref_indices <- which(PanCan_stat$BH=='RRSO' & (PanCan_stat$group=='ffTissue' | PanCan_stat$group=='ffpe'))
RRSO_ref_mat <- PanCan_mat[RRSO_ref_indices,]
RRSO_ref_stat <- PanCan_stat[RRSO_ref_indices,]
RRSO_ref_exposure <- sig_exposure(SSmatrix = RRSO_ref_mat, sig_type = 'PanCan')
#### sort by CX1
RRSO_ref_exposure <- RRSO_ref_exposure[order(-RRSO_ref_exposure$CX1),]
RRSO_ref_exposure <- select(RRSO_ref_exposure, !SS_sum)
#### plotting
RRSO_ref_exposure_plot <- ref_exposure_plot(ref_exposure_df = RRSO_ref_exposure, sig_type = 'Pan-Cancer', ref_stat_df = RRSO_ref_stat) + labs(title = 'RRSO ffTissue&ffpe') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('Pan-Cancer/Figures/ref_RRSO_exposure.pdf', plot = RRSO_ref_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### combined
RRSO_ref_exposure_plot_1 <- RRSO_ref_exposure_plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
HGSC_tumor_exposure_plot_1 <- HGSC_tumor_exposure_plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
PanCan_ref_exposure_plot <- ggarrange(Benign_tissue_exposure_plot, RRSO_ref_exposure_plot_1, HGSC_tumor_exposure_plot_1, ncol = 3, common.legend = 1, legend.grob = get_legend(Benign_tissue_exposure_plot), legend = 'right')
ggsave('Pan-Cancer/Figures/ref_exposure_all.pdf', plot = PanCan_ref_exposure_plot, dpi = 600, width = 10, height = 5, units = 'in')





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
panConusig_exposure <- sig_exposure(SSmatrix = panConusig_mat,sig_type = 'panConusig')

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

## 1. differences in ploidy and cellularity
for (sampleID in sample_df$Sample) {
  if (!(sampleID %in% panConusig_stat$sample)) {
    sample_df[which(sample_df$Sample==sampleID),'enrich_panConusig_1'] <- NA
  } else if (sampleID %in% panConusig_stat$sample) {
    sample_df[which(sample_df$Sample==sampleID),'enrich_panConusig_1'] <- panConusig_stat[which(panConusig_stat$sample==sampleID),'enrich_panConusig_1']
  }
}
panConusig_sample_df <- sample_df %>% filter(!is.na(enrich_panConusig_1))
panConusig_sample_df$Group <- factor(panConusig_sample_df$Group, levels = groups)
panConusig_sample_df$enrich_panConusig_1 <- factor(panConusig_sample_df$enrich_panConusig_1, levels = c('CN1','CN2','CN7','CN9','CN19','CN20','CN21','CN23'))
### plotting
ploidy_plot <-  ggboxplot(panConusig_sample_df, x='Group', y='Ploidy', 
                          color = 'black', fill = 'enrich_panConusig_1', palette = 'Paired', 
                          xlab = 'Group', ylab = 'Ploidy') +
  theme_classic()
ploidy_plot <- ploidy_plot + labs(fill = 'panConusig signatures')
ggsave('panConusig/Figures/panConusig_ploidy_all.pdf', plot = ploidy_plot, dpi = 600, width = 10, height = 7, units = 'in')
cellularity_plot <-  ggboxplot(panConusig_sample_df, x='Group', y='Cellularity', 
                               color = 'black', fill = 'enrich_panConusig_1', palette = 'Paired', 
                               xlab = 'Group', ylab = 'Cellularity') +
  theme_classic()
cellularity_plot <- cellularity_plot + labs(fill = 'panConusig signatures')
ggsave('panConusig/Figures/panConusig_cellularity_all.pdf', plot = cellularity_plot, dpi = 600, width = 10, height = 7, units = 'in') 
### stats analysis
sub_df <- panConusig_sample_df %>% filter(Group=='ffTumor')
# ploidy
kruskal.test(sub_df$Ploidy~sub_df$enrich_panConusig_1, data = sub_df) # p = 0.13
# cellularity
kruskal.test(sub_df$Cellularity~sub_df$enrich_panConusig_1, data = sub_df) # p = 0.17

## 1. distribution of the delta-value
panConusig_dist_delta_plot <- dis_plot_delta(panConusig_stat)
ggsave(filename = 'panConusig/Figures/delta_val_dist.pdf', plot = panConusig_dist_delta_plot, dpi = 600, width = 10, height = 7, units = 'in')
quantile(panConusig_stat$delta_val, 0.1) # 0.013 top 10%
# delta value cut-off at 0.25
for (groupID in groups) {
  print(groupID)
  group_df <- panConusig_stat[which(panConusig_stat$group==groupID),]
  sample_sum <- nrow(group_df)
  print(sample_sum)
  group_df <- group_df[which(group_df$delta_val<0.25 & group_df$panConusig_similarity_2>=0.5),]
  sample_size <- nrow(group_df)
  print(sample_size)
  print(100*sample_size/sample_sum)
}
panConusig_stat_25 <- panConusig_stat[which(panConusig_stat$delta_val<0.25 & panConusig_stat$panConusig_similarity_2>=0.5),]
groups2 <- c('ArchivalVS','MaNiLaVS','ffTissue','Endome','ffpe','Blood','ffPlasma')
#### top 2 enrich signature
panConusig_info_df_25 <- sig_stat_enrich(stat_df = panConusig_stat_25, enrich_num = 2, sig_type = 'panConusig')
panConusig_info_df_25$enrich_panConusig_2 <- factor(panConusig_info_df_25$enrich_panConusig_2, levels = c('CN1','CN4','CN6','CN7','CN9','CN17','CN18','CN19','CN20','CN21','CN22','CN23','CN24','CN25'))
# plot the graph
panConusig_enrich_plot_2 <- ggplot(data = panConusig_info_df_25, aes(x=BH, y=percentage, fill=enrich_panConusig_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='panConusig signatures') +
  scale_fill_ucscgb() +
  theme_classic() + facet_grid(cols = vars(group2))
panConusig_enrich_plot_2 <- panConusig_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('panConusig/Figures/panConusig_sig_BH_025_2.pdf', plot = panConusig_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
write_xlsx(panConusig_info_df_25, path = 'panConusig_group_025_2.xlsx')

# delta value cut-off at 0.1
for (groupID in groups) {
  print(groupID)
  group_df <- panConusig_stat[which(panConusig_stat$group==groupID),]
  sample_sum <- nrow(group_df)
  print(sample_sum)
  group_df <- group_df[which(group_df$delta_val<0.1 & group_df$panConusig_similarity_2>=0.5),]
  sample_size <- nrow(group_df)
  print(sample_size)
  print(100*sample_size/sample_sum)
}
panConusig_stat_01 <- panConusig_stat[which(panConusig_stat$delta_val<0.1 & panConusig_stat$panConusig_similarity_2>=0.5),]
groups2 <- c('ArchivalVS','MaNiLaVS','ffTissue','Endome','ffpe','Blood','ffPlasma')
#### top 2 enrich signature
panConusig_info_df_01 <- sig_stat_enrich(stat_df = panConusig_stat_01, enrich_num = 2, sig_type = 'panConusig')
panConusig_info_df_01$enrich_panConusig_2 <- factor(panConusig_info_df_01$enrich_panConusig_2, levels = c('CN1','CN6','CN9','CN17','CN18','CN19','CN20','CN21','CN23','CN24','CN25'))
# plot the graph
panConusig_enrich_plot_2 <- ggplot(data = panConusig_info_df_01, aes(x=BH, y=percentage, fill=enrich_panConusig_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='panConusig signatures') +
  scale_fill_ucscgb() +
  theme_classic() + facet_grid(cols = vars(group2))
panConusig_enrich_plot_2 <- panConusig_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('panConusig/Figures/panConusig_sig_BH_01_2.pdf', plot = panConusig_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
write_xlsx(panConusig_info_df_01, path = 'panConusig_group_01_2.xlsx')

# after adjusting the delta-value cutoff
## 0.25 for MaNiLaVS and Endome
## 0.1 for all other groups
panConusig_stat_01 <- panConusig_stat_01 %>% filter(group!='MaNiLaVS' & group!='Endome')
panConusig_stat_25 <- panConusig_stat_25 %>% filter(group=='MaNiLaVS' | group=='Endome')
panConusig_stat_cut <- rbind(panConusig_stat_01, panConusig_stat_25)
#### top 2 enrich signature
panConusig_info_df_cut <- sig_stat_enrich(stat_df = panConusig_stat_cut, enrich_num = 2, sig_type = 'panConusig')
panConusig_info_df_cut$enrich_panConusig_2 <- factor(panConusig_info_df_cut$enrich_panConusig_2, levels = c('CN1','CN6','CN9','CN17','CN18','CN19','CN20','CN21','CN23','CN24','CN25'))
# plot the graph
panConusig_enrich_plot_2 <- ggplot(data = panConusig_info_df_cut, aes(x=BH, y=percentage, fill=enrich_panConusig_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='panConusig signatures') +
  scale_fill_ucscgb() +
  theme_classic() + facet_grid(cols = vars(group2))
panConusig_enrich_plot_2 <- panConusig_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('panConusig/Figures/panConusig_sig_BH_cut_2.pdf', plot = panConusig_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
write_xlsx(panConusig_info_df_cut, path = 'panConusig_group_cut_2.xlsx')




## 2. exposure distribution
### rank by group
panConusig_color_bar_group <- color_bar_prep(stat_df = panConusig_stat, exposure_df = panConusig_exposure, grouping_type = 'groups', sig_type = 'panConusig', sort = TRUE)
panConusig_exposure_plot_group <- plot_exposure(stat_df = panConusig_stat, exposure_df = panConusig_exposure, grouping_type = 'groups', sig_type = 'panConusig', sort = TRUE)
panConusig_exposure_plot_group <- panConusig_exposure_plot_group %>% insert_bottom(panConusig_color_bar_group, height = 0.03)
ggsave(filename = 'panConusig/Figures/exposure_plot_groups_sort.pdf', plot = panConusig_exposure_plot_group, dpi = 600, width = 30, height = 7, units = 'in')
### rank by BH
panConusig_color_bar_BH <- color_bar_prep(stat_df = panConusig_stat, exposure_df = panConusig_exposure, grouping_type = 'BH_group', sig_type = 'panConusig', sort = TRUE)
panConusig_exposure_plot_BH <- plot_exposure(stat_df = panConusig_stat, exposure_df = panConusig_exposure, grouping_type = 'BH_group', sig_type = 'panConusig', sort = TRUE)
panConusig_exposure_plot_BH <- panConusig_exposure_plot_BH %>% insert_bottom(panConusig_color_bar_BH, height = 0.03)
ggsave(filename = 'panConusig/Figures/exposure_plot_BH_sort.pdf', plot = panConusig_exposure_plot_BH, dpi = 600, width = 30, height = 7, units = 'in')
### details for the ArchivalVS group
panConusig_archival_stat <- panConusig_stat %>% filter(group=='ArchivalVS')
archival_sample <- panConusig_archival_stat$sample
panConusig_archival_mat <- panConusig_mat %>% filter(sample %in% archival_sample)
panConusig_archival_exposure <- sig_exposure(SSmatrix = panConusig_archival_mat,sig_type = 'panConusig')
panConusig_color_bar_BH_AVS <- color_bar_prep(stat_df = panConusig_archival_stat, exposure_df = panConusig_archival_exposure, grouping_type = 'BH_group', sig_type = 'panConusig')
panConusig_exposure_plot_BH_AVS <- plot_exposure(stat_df = panConusig_archival_stat, exposure_df = panConusig_archival_exposure, grouping_type = 'BH_group', sig_type = 'panConusig')
panConusig_exposure_plot_BH_AVS <- panConusig_exposure_plot_BH_AVS %>% insert_bottom(panConusig_color_bar_BH_AVS, height = 0.03)
ggsave(filename = 'panConusig/Figures/exposure_plot_BH_AVS.pdf', plot = panConusig_exposure_plot_BH_AVS, dpi = 600, width = 30, height = 7, units = 'in')

## 3. panConusig signatures enrichment
### in this part, combine ffTumor and ffTissue into ffTissue
groups2 <- c('ArchivalVS','MaNiLaVS','ffTissue','Endome','ffpe','ffPlasma')

#### top 1 enriched signature
panConusig_info_df <- sig_stat_enrich(stat_df = panConusig_stat, enrich_num = 1, sig_type = 'panConusig')
panConusig_info_df$enrich_panConusig_1 <- factor(panConusig_info_df$enrich_panConusig_1, levels = c('CN1','CN2','CN7','CN9','CN19','CN20','CN21','CN23'))
# plot the graph
panConusig_enrich_plot_1 <- ggplot(data = panConusig_info_df, aes(x=BH, y=percentage, fill=enrich_panConusig_1)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='panConusig signatures') +
  scale_fill_ucscgb() +
  theme_classic() + facet_grid(cols = vars(group2))
panConusig_enrich_plot_1 <- panConusig_enrich_plot_1 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('panConusig/Figures/panConusig_sig_BH_all.pdf', plot = panConusig_enrich_plot_1, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
write_xlsx(panConusig_info_df, path = 'panConusig_group_1.xlsx')

#### top 2 enrich signature
panConusig_info_df <- sig_stat_enrich(stat_df = panConusig_stat, enrich_num = 2, sig_type = 'panConusig')
panConusig_info_df$enrich_panConusig_2 <- factor(panConusig_info_df$enrich_panConusig_2, levels = c('CN1','CN4','CN6','CN7','CN9','CN17','CN18','CN19','CN20','CN21','CN22','CN23','CN24','CN25'))
# plot the graph
panConusig_enrich_plot_2 <- ggplot(data = panConusig_info_df, aes(x=BH, y=percentage, fill=enrich_panConusig_2)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_stack(reverse = T)) +
  labs(x='Groups', y='Percentage (%)', fill='panConusig signatures') +
  scale_fill_ucscgb() +
  theme_classic() + facet_grid(cols = vars(group2))
panConusig_enrich_plot_2 <- panConusig_enrich_plot_2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave('panConusig/Figures/panConusig_sig_BH_all_2.pdf', plot = panConusig_enrich_plot_2, dpi = 600, width = 10, height = 7, units = 'in')
# save the info dataframe
write_xlsx(panConusig_info_df, path = 'panConusig_group_2.xlsx')


## 4. VS samples time point exposure
VS_sheet <- read.xlsx(file = 'MaNiLa_All_samplesheet_240311_groups.xlsx', sheetName = 'VS_samples')
VS_stat <- panConusig_stat[which(panConusig_stat$group=='ArchivalVS' | panConusig_stat$group=='MaNiLaVS'),]
for (sampleID in VS_stat$sample) {
  VS_stat[which(VS_stat$sample==sampleID), 'Diagnostic_time'] <- VS_sheet[which(VS_sheet$Library==sampleID), 'Diagnostic_time']
}

### extract the similarity matrix
indices <- rownames(VS_stat)
VS_mat <- panConusig_mat[indices,]
VS_exposure <- sig_exposure(SSmatrix = VS_mat,sig_type = 'panConusig')
VS_color_bar_BH <- color_bar_prep(stat_df = VS_stat, exposure_df = VS_exposure, grouping_type = 'BH_group', sig_type = 'panConusig')
VS_exposure_plot_BH <- plot_exposure(stat_df = VS_stat, exposure_df = VS_exposure, grouping_type = 'BH_group', sig_type = 'panConusig')
VS_exposure_plot_BH <- VS_exposure_plot_BH %>% insert_bottom(VS_color_bar_BH, height = 0.03)
ggsave('panConusig/Figures/VS_exposure.pdf', plot = VS_exposure_plot_BH, dpi = 600, width = 10, height = 7, units = 'in')

### for Benign VS samples
Benign_VS_stat <- VS_stat[which(VS_stat$BH=='Benign'),]
unique(Benign_VS_stat$Diagnostic_time) # only "diagnostic"
indices <- rownames(Benign_VS_stat)
Benign_VS_mat <- VS_mat[indices,]
#### exposure
Benign_VS_exposure <- sig_exposure(SSmatrix = Benign_VS_mat,sig_type = 'panConusig')
#### sort by CN1
Benign_VS_exposure <- Benign_VS_exposure[order(-Benign_VS_exposure$CN1),]
Benign_VS_exposure <- select(Benign_VS_exposure,!SS_sum)
#### plotting
Benign_VS_plot <- VS_exposure_plot(VS_exposure_df = Benign_VS_exposure, sig_type = 'panConusig', VS_stat_df = Benign_VS_stat) + labs(title = 'diagnostic') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('panConusig/Figures/Benign_VS_exposure.pdf', plot = Benign_VS_plot, dpi = 600, width = 7, height = 5, units = 'in')

### for RRSO samples
RRSO_VS_stat <- VS_stat[which(VS_stat$BH=='RRSO'),]
unique(RRSO_VS_stat$Diagnostic_time) # "prediagnostic_6+", "diagnostic"
indices <- rownames(RRSO_VS_stat)
RRSO_VS_mat <- VS_mat[indices,]
#### exposure
RRSO_VS_exposure_in <- sig_exposure(SSmatrix = RRSO_VS_mat,sig_type = 'panConusig')
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
ggsave('panConusig/Figures/RRSO_VS_exposure.pdf', plot = RRSO_VS_plot, dpi = 600, width = 10, height = 5, units = 'in')

### for HGSC samples
HGSC_VS_stat <- VS_stat[which(VS_stat$BH=='HGSC'),]
unique(HGSC_VS_stat$Diagnostic_time) # "prediagnostic_0-6", "prediagnostic_6+", "diagnostic"
indices <- rownames(HGSC_VS_stat)
HGSC_VS_mat <- VS_mat[indices,]
#### exposure
HGSC_VS_exposure_in <- sig_exposure(SSmatrix = HGSC_VS_mat,sig_type = 'panConusig')
#### sort by s1
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
ggsave('panConusig/Figures/HGSC_VS_exposure.pdf', plot = HGSC_VS_plot, dpi = 600, width = 10, height = 5, units = 'in')

## 5. reference exposure plots
#### HGSC ffTumor 
HGSC_tumor_indices <- which(panConusig_stat$group=='ffTumor' & panConusig_stat$BH=='HGSC')
HGSC_tumor_mat <- panConusig_mat[HGSC_tumor_indices,]
HGSC_tumor_stat <- panConusig_stat[HGSC_tumor_indices,]
HGSC_tumor_exposure <- sig_exposure(SSmatrix = HGSC_tumor_mat, sig_type = 'panConusig')
#### sort by CX1
HGSC_tumor_exposure <- HGSC_tumor_exposure[order(-HGSC_tumor_exposure$CN1),]
HGSC_tumor_exposure <- select(HGSC_tumor_exposure,!SS_sum)
#### plotting
HGSC_tumor_exposure_plot <- ref_exposure_plot(ref_exposure_df = HGSC_tumor_exposure, sig_type = 'panConusig', ref_stat_df = HGSC_tumor_stat) + labs(title = 'HGSC ffTumor') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('panConusig/Figures/ref_HGSC_exposure.pdf', plot = HGSC_tumor_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### Benign ffTissue
Benign_tissue_indices <- which(panConusig_stat$group=='ffTissue' & panConusig_stat$BH=='Benign')
Benign_tissue_mat <- panConusig_mat[Benign_tissue_indices,]
Benign_tissue_stat <- panConusig_stat[Benign_tissue_indices,]
Benign_tissue_exposure <- sig_exposure(SSmatrix = Benign_tissue_mat, sig_type = 'panConusig')
#### sort by CX1
Benign_tissue_exposure <- Benign_tissue_exposure[order(-Benign_tissue_exposure$CN1),]
Benign_tissue_exposure <- select(Benign_tissue_exposure, !SS_sum)
#### plotting
Benign_tissue_exposure_plot <- ref_exposure_plot(ref_exposure_df = Benign_tissue_exposure, sig_type = 'panConusig', ref_stat_df = Benign_tissue_stat) + labs(title = 'Benign ffTissue') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('panConusig/Figures/ref_Benign_exposure.pdf', plot = Benign_tissue_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### RRSO ffTissue+ffpe
RRSO_ref_indices <- which(panConusig_stat$BH=='RRSO' & (panConusig_stat$group=='ffTissue' | panConusig_stat$group=='ffpe'))
RRSO_ref_mat <- panConusig_mat[RRSO_ref_indices,]
RRSO_ref_stat <- panConusig_stat[RRSO_ref_indices,]
RRSO_ref_exposure <- sig_exposure(SSmatrix = RRSO_ref_mat, sig_type = 'panConusig')
#### sort by CX1
RRSO_ref_exposure <- RRSO_ref_exposure[order(-RRSO_ref_exposure$CN1),]
RRSO_ref_exposure <- select(RRSO_ref_exposure, !SS_sum)
#### plotting
RRSO_ref_exposure_plot <- ref_exposure_plot(ref_exposure_df = RRSO_ref_exposure, sig_type = 'panConusig', ref_stat_df = RRSO_ref_stat) + labs(title = 'RRSO ffTissue&ffpe') + theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave('panConusig/Figures/ref_RRSO_exposure.pdf', plot = RRSO_ref_exposure_plot, dpi = 600, width = 7, height = 5, units = 'in')

#### combined
RRSO_ref_exposure_plot_1 <- RRSO_ref_exposure_plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
HGSC_tumor_exposure_plot_1 <- HGSC_tumor_exposure_plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
panConusig_ref_exposure_plot <- ggarrange(Benign_tissue_exposure_plot, RRSO_ref_exposure_plot_1, HGSC_tumor_exposure_plot_1, ncol = 3, common.legend = 1, legend.grob = get_legend(Benign_tissue_exposure_plot), legend = 'right')
ggsave('panConusig/Figures/ref_exposure_all.pdf', plot = panConusig_ref_exposure_plot, dpi = 600, width = 10, height = 5, units = 'in')

