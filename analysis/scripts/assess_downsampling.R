###########################################################################################
#########                                                                         #########
#########                             Daniel Kunz                                 #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             06/03/2017                                  #########
#########                              EMBL-EBI                                   #########
#########                              Punt-Seq                                   #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                            Microbiome of the Cam                            #########
###########################################################################################
##### Assess the effect of downsampling the number of reads in the different samples   ####
##### for the identification of species.                                               ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

library(data.table)
library(tidyr)
library(plyr)
library(readr)
library(ggplot2)
library(ggthemes)


source("analysis/scripts/rank.R")

set.seed(53043)

#### FUNCTIONS. ####

## Calculate the mean and the standard deviation for each group.
#  From http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
#
# data: a dataframe
# varname: the name of a column containing the variable to be summarized
# groupnames: vector of column names to be used as grouping variables

data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE), sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum <- ddply(data, groupnames, .fun=summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}



#### 1. Input arguments. ####

input_folder <- '/Users/dem44/Desktop/Punt_seq/downsampling/raw_data/' # Folder containing all the input files (centrifuge results)
output_folder <- '/Users/dem44/Desktop/Punt_seq/downsampling/output_data/' # Folder where the table with the results and the plots will be generated

read_steps <- seq(0.1, 1, by=0.1) # % reads that will be included in each step
n_rep <- 10 # Number of times that the downsampling will be performed in each step


#### 2. Run the analysis. ####

files_centrifuge = list.files("analysis/data/", patter="tab$", recursive=T, full.names=T)
final_df <- data.frame(matrix(NA, nrow=length(read_steps)*length(files_centrifuge)*n_rep, ncol=8))
colnames(final_df) <- c('SampleID', 'Time', 'Place', 'Reads_fraction',
                        'Species', 'Genus', 'Family', 'Phylum');

i <- 1
for(f_path in files_centrifuge){ # For each sample

  print(paste0('Processing sample ', strsplit(f_path, '/')[[1]][length(strsplit(f_path, '/')[[1]])], ' ...'))
  out_centrifuge = read_tsv(f_path, col_types="cciiiiii")

  sample_id <- gsub('.tab', '', strsplit(f_path, 'centrifuge_classification_')[[1]][2])
  time <- strsplit(sample_id, '_')[[1]][1]; place <- strsplit(sample_id, '_')[[1]][2]

  for(read_fraction in read_steps){ # For each read fraction

    readIDs = unique(out_centrifuge$readID)
    n_reads_sample = round(length(readIDs)*read_fraction)

      for(r in 1:n_rep){ # Bootstrapping

        print("Start Bootstrapping")

        final_df[i, 1:4] <- c(sample_id, time, place, read_fraction)

        # subsample
        readIDs_sample = sample(readIDs, n_reads_sample)

        subsample = out_centrifuge %>%
                      filter(readID %in% readIDs_sample)

        print("Start Ranking")
        rank_counts = get_rank_counts(subsample, taxID_info)
        final_df[i,c("Species", "Genus", "Family", "Phylum")] = rank_counts[c("n_species", "n_genus", "n_family", "n_phylum")]
        i <- i + 1
      }
  }
}


write_csv(final_df, "analysis/data/downsampling_df.csv")




#### 3. Plot the results. ####


final_df = read_csv("analysis/data/downsampling_df.csv")

summary_species <- data_summary(final_df, varname="Species", groupnames=c('SampleID', 'Time', 'Place', 'Reads_fraction'));

for(month in unique(final_df$Time)){

  summary_species_month <- summary_species[summary_species$Time==month,]

  ## Species plots

  ggplot(data=summary_species_month) +
    geom_pointrange(aes(x=as.numeric(Reads_fraction), y=Species, ymin=Species-sd, ymax=Species+sd, col=SampleID)) +
    geom_line(aes(x=as.numeric(Reads_fraction), y=Species, col=SampleID)) +
    facet_wrap(~SampleID) + xlim(c(0,1)) +
    theme_bw() +
    scale_color_stata(guide=F) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    ylab("Number of species") + xlab("Fraction of reads in the subsample")


  ggsave(paste0('species_plot_', month, '.pdf'),
         height=3*length(unique(summary_species_month$SampleID))/4, width=8);

}


summary_genus <- data_summary(final_df, varname="Genus", groupnames=c('SampleID', 'Time', 'Place', 'Reads_fraction'));

for(month in unique(final_df$Time)){

  summary_genus_month <- summary_genus[summary_genus$Time==month,]

  ## Species plots

  ggplot(data=summary_genus_month) +
    geom_pointrange(aes(x=as.numeric(Reads_fraction), y=Genus, ymin=Genus-sd, ymax=Genus+sd, col=SampleID)) +
    geom_line(aes(x=as.numeric(Reads_fraction), y=Genus, col=SampleID)) +
    facet_wrap(~SampleID) + xlim(c(0,1)) +
    theme_bw() +
    scale_color_stata(guide=F) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    ylab("Number of geni") + xlab("Fraction of reads in the subsample")


  ggsave(paste0('genus_plot_', month, '.pdf'),
         height=3*length(unique(summary_genus_month$SampleID))/4, width=8);

}


#### End of the script ####
