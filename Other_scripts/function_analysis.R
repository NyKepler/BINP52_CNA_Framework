# Title: function_analysis.R
# Author: Guyuan TANG
# Date: 2024/4/20

# Description: This script contains the functions required for plotting and statistics in signature analysis.


library(tidyverse)
library(ggpubr)
library(xlsx)
library(aplot)
library(ggsci)
library(writexl)

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
  } else if (sig_type=='CerCN') {
    output_df <- input_df %>% mutate(delta_val = CerCN_similarity_1 - CerCN_similarity_2)
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
sig_exposure <- function(SSmatrix, sig_type, threshold=0) {
  # remove the column 'enrich'
  SSmatrix <- SSmatrix[,-2]
  # turn the similarity values that are below threshold into 0
  threshold = as.numeric(threshold)
  for ( i in 1:nrow(SSmatrix)) {
    for (j in 2:ncol(SSmatrix)) {
      if (SSmatrix[i,j]<threshold) {
        SSmatrix[i,j] <- 0
      }
    }
  }
  if (sig_type == 'CN') {
    n = 7+1
  } else if (sig_type == 'PanCan' | sig_type=='Pan-Cancer') {
    n = 17+1
  } else if (sig_type == 'panConusig') {
    n = 25+1
  } else if (sig_type == 'CerCN') {
    n = 6+1
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
      SS_df[which(SS_df$sample==sampleID),'group_VS'] <- sample_df[which(sample_df$Sample==sampleID), 'Group_VS']
      SS_df[which(SS_df$sample==sampleID),'group_time'] <- sample_df[which(sample_df$Sample==sampleID), 'Group_time']
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
    if (grouping_type=='group_VS') {
      groups <- c('Blood','Plasma','Normal tissues','Tumor','Cervical')
      ### rank by groups
      x_lab_name <- 'Sample types'
      for (groupID in groups) {
        indices <- which(stat_df$group_VS==groupID)
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
    if (grouping_type=='group_VS') {
      groups <- c('Blood','Plasma','Normal tissues','Tumor','Cervical')
      ### rank by groups
      x_lab_name <- 'Sample types'
      for (groupID in groups) {
        indices <- which(stat_df$group_VS==groupID)
        sort_df <- exposure_df[indices,]
        # sort the samples by the first signatures
        if (sig_type == 'CN') {
          sort_df <- sort_df[order(-sort_df$s1),]
        } else if (sig_type=='Pan-Cancer' | sig_type=='PanCan') {
          sort_df <- sort_df[order(-sort_df$CX1),]
        } else if (sig_type=='panConusig') {
          sort_df <- sort_df[order(-sort_df$CN1),]
        } else if (sig_type=='CerCN') {
          sort_df <- sort_df[order(-sort_df$CerCN1),]
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
        } else if (sig_type=='CerCN') {
          sort_df <- sort_df[order(-sort_df$CerCN1)]
        }
        exposure_out <- rbind(exposure_out, sort_df)
      }
    }
  }
  exposure_out <- select(exposure_out, !SS_sum)
  rownames(exposure_out) <- 1:nrow(exposure_out)
  # match the grouping information in order to add the coloring bars
  if (grouping_type=='group_VS') {
    for (sampleID in exposure_out$sample) {
      exposure_out[which(exposure_out$sample==sampleID),'col_group'] <- stat_df[which(stat_df$sample==sampleID),'group_VS']
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
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)) +
    scale_fill_brewer(palette = 'Set2')
  return(col_plot)
}

# define a function to plot the exposure of signatures
plot_exposure <- function(stat_df, exposure_df, grouping_type, sig_type, sort=TRUE) {
  col_num <- ncol(exposure_df)
  ### create a dataframe to store the information, and rank by either groups or BH_groups
  exposure_out <- as.data.frame(matrix(nrow = 0,ncol = col_num))
  if (sort==FALSE | sort==F) {
    if (grouping_type=='group_VS') {
      groups <- c('Blood','Plasma','Normal tissues','Tumor','Cervical')
      ### rank by groups
      x_lab_name <- 'Sample types'
      for (groupID in groups) {
        indices <- which(stat_df$group_VS==groupID)
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
    if (grouping_type=='group_VS') {
      groups <- c('Blood','Plasma','Normal tissues','Tumor','Cervical')
      ### rank by groups
      x_lab_name <- 'Sample types'
      for (groupID in groups) {
        indices <- which(stat_df$group_VS==groupID)
        sort_df <- exposure_df[indices,]
        # sort the samples by the first signatures
        if (sig_type == 'CN') {
          sort_df <- sort_df[order(-sort_df$s1),]
        } else if (sig_type=='Pan-Cancer' | sig_type=='PanCan') {
          sort_df <- sort_df[order(-sort_df$CX1),]
        } else if (sig_type=='panConusig') {
          sort_df <- sort_df[order(-sort_df$CN1),]
        } else if (sig_type=='CerCN') {
          sort_df <- sort_df[order(-sort_df$CerCN1),]
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
        } else if (sig_type=='CerCN') {
          sort_df <- sort_df[order(-sort_df$CerCN1),]
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
    exposure_out$signatures <- factor(exposure_out$signatures, levels = paste0('CX',c(1:17)))
  } else if (sig_type == 'CN') {
    exposure_out$signatures <- factor(exposure_out$signatures, levels = paste0('s', c(1:7)))
  } else if (sig_type == 'panConusig') {
    exposure_out$signatures <- factor(exposure_out$signatures, levels = paste0('CN',c(1:25)))
  } else if (sig_type == 'CerCN') {
    exposure_out$signatures <- factor(exposure_out$signatures, levels = paste0('CerCN', c(1:6)))
  }
  
  exposure_plot <- ggplot(data = exposure_out, aes(x=sample, y=exposure, fill=signatures)) +
    geom_bar(stat = 'identity', width = 0.7, position = position_stack(reverse = T)) +
    labs(fill = paste0(sig_type, ' signatures'), y='Exposure') +
    scale_y_continuous(expand = c(0.01,0)) +
    theme(panel.background = element_blank(),
          axis.title.y = element_text(size = 15, hjust = 0.5),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)) +
    scale_fill_ucscgb()
  return(exposure_plot)
}

# define a function to prepare the signature enrichment statistical analysis
sig_stat_enrich <- function(stat_df, enrich_num, sig_type) {
  ### in this part, combine ffTumor and ffTissue into ffTissue
  groups2 <- c('prediagnostic VS','diagnostic VS','ffTissue','Endome','FFPE','Blood','ffPlasma')
  BH_group <- c('Benign','RRSO','HGSC')
  if (enrich_num==1 & sig_type=='CN') {
    info_df <- stat_df %>% group_by(group_time, BH, enrich_CN_1) %>%
      summarise(mean_CN_similarity = mean(CN_similarity_1),
                sd_CN_similarity = sd(CN_similarity_1),
                sample_size = n())
  } else if (enrich_num==2 & sig_type=='CN') {
    info_df <- stat_df %>% group_by(group_time, BH, enrich_CN_2) %>%
      summarise(mean_CN_similarity = mean(CN_similarity_2),
                sd_CN_similarity = sd(CN_similarity_2),
                sample_size = n())
  } else if (enrich_num==1 & sig_type=='PanCan') {
    info_df <- stat_df %>% group_by(group_time, BH, enrich_PanCan_1) %>%
      summarise(mean_PanCan_similarity = mean(PanCan_similarity_1),
                sd_PanCan_similarity = sd(PanCan_similarity_1),
                sample_size = n())
  } else if (enrich_num==2 & sig_type=='PanCan') {
    info_df <- stat_df %>% group_by(group_time, BH, enrich_PanCan_2) %>%
      summarise(mean_PanCan_similarity = mean(PanCan_similarity_2),
                sd_PanCan_similarity = sd(PanCan_similarity_2),
                sample_size = n())
  } else if (enrich_num==1 & sig_type=='panConusig') {
    info_df <- stat_df %>% group_by(group_time, BH, enrich_panConusig_1) %>%
      summarise(mean_panConusig_similarity = mean(panConusig_similarity_1),
                sd_panConusig_similarity = sd(panConusig_similarity_1),
                sample_size = n())
  } else if (enrich_num==2 & sig_type=='panConusig') {
    info_df <- stat_df %>% group_by(group_time, BH, enrich_panConusig_2) %>%
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
  info_df$group_time <- factor(info_df$group_time, levels = groups2)
  info_df$BH <- factor(info_df$BH, levels = BH_group)
  for (groupID in groups2) {
    for (BH_ID in BH_group) {
      info_df[which(info_df$group_time==groupID & info_df$BH==BH_ID),"sample_sum"] <- sum(info_df[which(info_df$group_time==groupID & info_df$BH==BH_ID),"sample_size"])
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
  } else if (sig_type == 'CerCN') {
    exposure_out$signatures <- factor(exposure_out$signatures, levels = paste0('CerCN', c(1:6)))
  }
  # add the diagnostic timepoints
  exposure_out$diagnostic_time <- NA
  for (sampleID in sample_rank) {
    exposure_out[which(exposure_out$sample==sampleID),'diagnostic_time'] <- VS_stat_df[which(VS_stat_df$sample==sampleID),'Diagnostic_time']
  }
  exposure_out$diagnostic_time <- factor(exposure_out$diagnostic_time, levels = c('prediagnostic_6+', 'prediagnostic_0-6', 'diagnostic'))
  
  #use different color palettes for different signatures
  if (sig_type=='CerCN') {
    exposure_plot <- ggplot(data = exposure_out, aes(x=sample, y=exposure, fill=signatures)) +
      geom_bar(stat = 'identity', width = 0.8, position = position_stack(reverse = T)) +
      labs(fill = paste0(sig_type, ' signatures'), y='Exposure') +
      scale_y_continuous(expand = c(0.01,0)) +
      theme(panel.background = element_blank(),
            axis.title.y = element_text(size = 15, hjust = 0.5),
            axis.text.y = element_text(size = 15),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "right",
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15)) +
      scale_fill_brewer(palette = 'Set2')
  } else {
    exposure_plot <- ggplot(data = exposure_out, aes(x=sample, y=exposure, fill=signatures)) +
      geom_bar(stat = 'identity', width = 0.8, position = position_stack(reverse = T)) +
      labs(fill = paste0(sig_type, ' signatures'), y='Exposure') +
      scale_y_continuous(expand = c(0.01,0)) +
      theme(panel.background = element_blank(),
            axis.title.y = element_text(size = 15, hjust = 0.5),
            axis.text.y = element_text(size = 15),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "right",
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15)) +
      scale_fill_ucscgb()
  }
  
  
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
