#!/usr/bin/env Rscript
#Input Layout: Rscript calc_fstat.R -f <path_to_f_stat.csv> -m <path_to_master_file> -p <path_to_pop_iids.csv> -t1 <trait1> -t2 <trait2>

library(rio)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(optparse)

#Parse options similar to argparser in Python
option_list = list(
  make_option(c("-f", "--f_stat"), type="character", default=NULL, 
              help="csv for f_stat"),
  make_option(c("-m", "--master"), type="character", default=NULL, 
              help="master table"),
  make_option(c("-p", "--pop_iids"), type="character", default=NULL, 
              help="csv with pop_idds, generated from downsample.py"),
  make_option(c("-t1", "--trait1"), type="character", default=NULL, 
              help="csv of scores for individuals for low f-stat trait"),
  make_option(c("-t2", "--trait2"), type="character", default=NULL, 
              help="csv of scores for individuals for high f-stat trait")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Loading options into variables
f_file <- opt$f_stat
master_file <- opt$master 
trait1_file <- opt$trait1
trait2_file <- opt$trait2
pop_file <- opt$pop_iids

#Ref to local machine for testing 
# f_file <- '~/Desktop/Gini-PGS/f_stat.csv'
# master_file <-'~/Desktop/Gini-PGS/table_PLR_100k_sum_100_zp.txt' 
# trait1_file <- '~/Desktop/Gini-PGS/final/250.1-1KG_PLR.csv'
# trait2_file <- '~/Desktop/Gini-PGS/final/darker_skin0-1KG_PLR.csv'
# pop_file <- '~/Desktop/Gini-PGS/UKB/randomized_pop_iids.csv'

#Loading CSV containing f_stat for each trait
f_df <- import(f_file)

#Load in master dataframe, containing every metric 
master_df <- import(master_file)
master_df <- master_df %>% select(c(prive_code, trait, description))
f_df <- merge(f_df, master_df, by ='prive_code')
rm(master_df)

#Load CSV for scores of each individual for trait1 
trait1_df <- import(trait1_file)
trait1_df <- trait1_df[,2:ncol(trait1_df)] #drop index column
trait1_df$sum <- rowSums(trait1_df[,3:ncol(trait1_df)])
trait1_df <- trait1_df %>% select(c(FID, IID, sum))

#Load CSV for scores of each individual for trait2 
trait2_df <- import(trait2_file)
trait2_df <- trait2_df[,2:ncol(trait2_df)] #drop index column
trait2_df$sum <- rowSums(trait2_df[,3:ncol(trait2_df)])
trait2_df <- trait2_df %>% select(c(FID, IID, sum))

#Load in the CSV with population iids 
pop_iids_df <- import(pop_file)
pop_iids_df <- pop_iids_df[,1:2]

#Merging dataframes to have figure out which ID belongs to which population
trait1_df <- merge(trait1_df, pop_iids_df, by = 'IID')
trait2_df <- merge(trait2_df, pop_iids_df, by = 'IID')

# #Calculate Mean and rank - This is optional - Simply to check the rank - put this elsewhere
# mean_df <- data.frame(ancestry = character(), means = double())
# for (i in unique(final_df$ancestry)) {
#   df_subset <- final_df[final_df$ancestry==i,]
#   amean <- mean(as.numeric(df_subset[,3]))
#   mean_df[nrow(mean_df)+1,] <- c(i, amean)
# }
# mean_df['rank'] <- rank(as.numeric(mean_df[,2]))

#anova to calculate f-stat
one.way <- aov(as.numeric(sum) ~ ancestry, data = trait1_df)
f_val <- round(summary(one.way)[[1]]["F value"][[1]][[1]], digits = 3)
f_val <- round(log(f_val), digits = 3)
#p_val <- summary(one.way)[[1]]["Pr(>F)"][[1]][[1]]

one.way2 <- aov(as.numeric(sum) ~ ancestry, data = trait2_df)
f_val2 <- round(summary(one.way2)[[1]]["F value"][[1]][[1]], digits = 3)
f_val2 <- round(log(f_val2), digits = 3)
#p_val <- summary(one.way)[[1]]["Pr(>F)"][[1]][[1]]

# #Filter to UK only - remove this
# final_df <- subset(final_df, ancestry == 'United Kingdom')

title1 <- paste("log(F) =", f_val, sep=" ")
histogram1 <- ggplot(trait1_df, aes(x=as.numeric(sum), fill=ancestry)) + geom_density(color='#e9ecef', alpha=0.6, position='identity') + theme_bw() + labs(x='Polygenic Score per UKBB Individual', y='Density') + ggtitle(paste('Type 1 Diabetes', title1, sep='\n'))
histogram1 <- histogram1 + theme(legend.title=element_blank())
print(histogram1)

title2 <- paste("log(F) =", f_val2, sep=" ")
histogram2 <- ggplot(trait2_df, aes(x=as.numeric(sum), fill=ancestry)) + geom_density(color='#e9ecef', alpha=0.6, position='identity') + theme_bw() + labs(x='Polygenic Score per UKBB Individual', y='Density') + ggtitle(paste('Skin Colour', title2, sep='\n'))
histogram2 <- histogram2 + theme(legend.title=element_blank()) 

print(histogram2)

ggarrange(histogram1, histogram2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")



