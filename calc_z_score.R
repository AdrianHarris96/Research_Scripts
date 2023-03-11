#Format and visualize PGS for individuals and populations 
#Standardize and visualize D for all pops per trait - most and least divergent
rm(list = ls())
dev.off()

library(pacman)
p_unload(all)
library(rio)
library(ggplot2)
library(tidyverse)

#Input files 
traits_file <- "~/Desktop/Gini-PGS/traits_table.txt"
dosage_dir <- "~/Desktop/Gini-PGS/PGS_traits/"
randomized_pops <- "~/Desktop/Gini-PGS/randomized_pop_iids.csv" #downsample, equal number of individauls for each population
len <- 5L #Top and bottom number of traits to pull

#Load in traits_table 
traits_table <- import(traits_file)
traits_table <- traits_table[order(traits_table$f_stat),]

most_divergent <- tail(traits_table, len)
least_divergent <- head(traits_table, len)

divergent <- rbind(most_divergent, least_divergent)
divergent_traits <- divergent$prive_code

calculate_z_score <- function(traitname) {
  trait_file <- paste(traitname, ".csv", sep="")
  #Load in dataframe of individuals and dosages (either a 0,1,2 at each SNP for 12000+ individuals)
  trait_doses_df <- import(paste(dosage_dir, trait_file, sep=""))
  
  #Drop first two columns 
  trait_doses_df <- trait_doses_df[, 3:ncol(trait_doses_df)]
  
  #Distribute the effect size across the remainder of the row
  skip_cols = c("chrom", "rsID", "A2", "prive_code", "effect_weight", "rsID_allele")
  
  #Creation of final dataframe
  final_df <- data.frame(matrix(ncol=ncol(trait_doses_df), nrow=1))
  colnames(final_df) <- colnames(trait_doses_df)
  final_df <- final_df %>% select(-skip_cols)
  
  #Calculate score per individual for the trait
  for (col in colnames(trait_doses_df)) {
    if (col %in% skip_cols) {
      next
    } else {
      trait_doses_df[col] <- trait_doses_df$effect_weight * trait_doses_df[col]
      individual_sum <- sum(trait_doses_df[col])
      final_df[col] <- individual_sum
    }
  }
  
  rm(trait_doses_df, col, individual_sum)
  
  #Transpose final_df 
  final_df <- as.data.frame(t(final_df))
  final_df <- final_df %>% rename("score" = "V1")
  final_df$IID <- rownames(final_df)
  
  #Mean across scores 
  mean_score <- mean(final_df$score)
  sd_score <- sd(final_df$score)
  final_df$z_score <- ((final_df$score - mean_score) / sd_score)
  
  #Load in group information 
  pop_info_df <- import(randomized_pops)
  pop_info_df <- pop_info_df[, 1:(ncol(pop_info_df)-2)]
  
  #Merge the two dataframes
  final_df <- merge(final_df, pop_info_df, by = "IID")
  
  #Start calculating distributions,means
  pop_means <- final_df %>% group_by(ancestry) %>% summarise_at(vars(z_score), list(pop_mean = mean))
  
  pop_means$prive_code <- traitname
  print(paste(trait, "completed!", sep=" "))
  return(pop_means)
}

means_df <- data.frame(ancestry = character(), pop_mean = double(), prive_code = character())

#Iterate through the most divergent traits 
for (trait in divergent_traits) {
  means_df <- rbind(means_df, calculate_z_score(trait))
}

#Load in proportions for ANC and DER - merge with means_df 
proportion_df <- import("~/Desktop/Gini-PGS/trait_proportions.csv")
means_df <- merge(means_df, proportion_df, by = "prive_code")
means_df$proportion_anc <- round(means_df$proportion_anc, 3)
means_df$proportion_der <- round(means_df$proportion_der, 3)
temp_df <- traits_table %>% select(c(prive_code, description, f_stat))
means_df <- merge(means_df, temp_df, by = "prive_code")

#Minor edits to the dataframe to include information on ancestral proportion for each trait (separated by : from trait name)
means_df$trait_anc <- paste(means_df$description, means_df$proportion_anc, sep=": ")
means_df$trait_der <- paste(means_df$description, means_df$proportion_der, sep=": ")
means_df <- means_df[order(means_df$f_stat),]

#Plot for a group of traits 
anc_plot <- ggplot(data = means_df, mapping = aes(x = pop_mean, y = trait_anc, color = ancestry)) + geom_point(size = 3.0) + theme_light() + labs(x='Z-score', y="Trait_ANC")
der_plot <- ggplot(data = means_df, mapping = aes(x = pop_mean, y = trait_der, color = ancestry)) + geom_point(size = 3.0) + theme_light() + labs(x='Z-score', y="Trait_DER")

ggsave(file = "~/Desktop/proportion_anc_and_z_score.jpeg", units = c("in"), width=12, height=6, dpi=300, anc_plot)
ggsave(file = "~/Desktop/proportion_der_and_z_score.jpeg", units = c("in"), width=12, height=6, dpi=300, der_plot)

