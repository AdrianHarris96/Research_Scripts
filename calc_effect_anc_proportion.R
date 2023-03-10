#Calculate Ancestral-Derived Proportions and Append to Traits Table
rm(list = ls())
dev.off()

library(rio)
library(ggplot2)
library(tidyverse)

#These could be provided as arguments as well
SNPs_file <- '~/Desktop/Gini-PGS/top100bin_SNPs_edited.csv'
proportions_file <- "~/Desktop/trait_proportions.csv"
traits_file <- '~/Desktop/Gini-PGS/traits_table.csv'

#Import table with ancestral and derived alleles
alleles_table <- import(SNPs_file)
#Note: A2 is the effect allele 

#Order the table first by the prive code 
alleles_table <- alleles_table[order(alleles_table$prive_code),]

#Determining the effect_increasing allele (via looping row by row)
alleles_table$effect_increasing_allele <- 'NA'
for (row in 1:nrow(alleles_table)) {
  if (alleles_table[row, 'effect_weight'] > 0) {
    alleles_table[row, 'effect_increasing_allele'] <- alleles_table[row, 'A2']
  } else {
    alleles_table[row, 'effect_increasing_allele'] <- alleles_table[row, 'A1']
  }
  rm(row)
}
#colnames(alleles_table)

#Trim the number of columns (removing redundant info.)
alleles_table <- alleles_table %>% select(c(prive_code, rsID, chrom, chr_position, ANC, REF, ALT, effect_increasing_allele))

#Generate list of unique traits 
traits <- unique(alleles_table$prive_code)

#Empty dataframe to fill 
proportion_df <- data.frame(prive_code = character(), proportion_anc = double(), proportion_der = double())

#Subset allele table and calculate the porportion of ANCvsDER for every prive code/trait
#Conditionals deal with the triallelic scenarios as well
for (trait in traits) {
  subset_df <- alleles_table[(alleles_table$prive_code == trait),]
  total_snps <- 0
  #print(nrow(subset_df))
  effect_ancestral <- 0 #If the effect allele is equal to the ancestral allele, add to this tally
  effect_derived <- 0
  for (row in 1:nrow(subset_df)) {
    if (subset_df[row, 'ANC'] != subset_df[row, 'REF'] & subset_df[row, 'ANC'] != subset_df[row, 'ALT'] & subset_df[row, 'ANC'] != subset_df[row, 'effect_increasing_allele']) {
      next #Skipping the row if it is tri-allelic (condition 1 and 2) AND there is not a match between ANC and effect_increasing_allele (condition 3)
    } else if (subset_df[row, 'effect_increasing_allele'] == subset_df[row, 'ANC']) {
      effect_ancestral <- effect_ancestral + 1
      total_snps <- total_snps + 1
    } else {
      effect_derived <- effect_derived + 1
      total_snps <- total_snps + 1
    }
    rm(row)
  }
  effect_anc <- effect_ancestral / total_snps
  effect_der <- effect_derived / total_snps
  proportion_df[nrow(proportion_df)+1,] <- c(trait, effect_anc, effect_der)
  rm(subset_df, total_snps, effect_ancestral, effect_derived, trait, effect_anc, effect_der)
}

#Write to CSV 
write.csv(proportion_df, file = proportions_file, row.names = FALSE)

#Merge proportion_df with traits_table 
traits_table <- import(traits_file)
traits_table$D_stat <- log10(traits_table$f_stat)
traits_table <- traits_table %>% select(c(prive_code, description, D_stat, group_consolidated))
merged_df <- merge(proportion_df, traits_table, by = 'prive_code')
merged_df$proportion_anc <- as.numeric(merged_df$proportion_anc)
merged_df$proportion_der <- as.numeric(merged_df$proportion_der)

#Plot proportion vs. D-stat 
# plot1 <- ggplot(data=merged_df, mapping = aes(x = D_stat, y = proportion_anc, color = group_consolidated)) + geom_smooth(method = "lm", se = FALSE) + geom_point(size = 1.5) + theme_light() + labs(x='D', y="Proportion in which effect increasing allele = ancestral")
# print(plot1)
plot1 <- ggplot(data=merged_df, mapping = aes(x = D_stat, y = proportion_anc, color = group_consolidated)) + geom_point(size=1) + geom_text(aes(label = description), size=2, nudge_y = 0.005) + theme_light() + labs(x='D', y="Proportion in which effect increasing allele = ancestral")
print(plot1)

plot2 <- ggplot(data=merged_df, mapping = aes(x = D_stat, y = proportion_der, color = group_consolidated)) + geom_point(size=1) + geom_text(aes(label = description), size=2, nudge_y = 0.005) + theme_light() + labs(x='D', y="Proportion in which effect increasing allele = derived")
print(plot2)

library(ggpubr)
g <- ggarrange(plot1, plot2, ncol = 2, nrow=1, legend = FALSE)
print(g)

#Save to JPEG
ggsave(file = "~/Desktop/proportions_anc_der.jpeg", units = c("in"), width=12, height=6, dpi=300, g)



