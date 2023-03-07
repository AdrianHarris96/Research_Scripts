#Generate Portability Example
rm(list = ls())
dev.off()

#Unload prior to loading to avoid conflicts
library(pacman)
p_unload(all)

#Load relavent packages
library(rio)
library(ggplot2)
library(tidyverse)
library(glue)

#Consider accepting these as arguments 
master_file <- '~/Desktop/Gini-PGS/traits_table.csv'
pop_centers_file <- '~/Desktop/Gini-PGS/UKBB-PGS/pop_centers.csv'

#Load in master dataframe
master_df <- read.csv(master_file, sep=',')

#Portability/slope dataframe created
m_df <- data.frame(prive_code = character(), description = character(), ancestry = character(), pcor_ratio = double())

#Calculation of partial correlation ratio
start <- grep('pcor_Caribbean', colnames(master_df))
stop <- grep('pcor_United', colnames(master_df))
for (row in 1:nrow(master_df)) {
  code <- master_df[row, 'prive_code']
  trait_name <- master_df[row, 'description']
  source_pop <- master_df[row, 'pcor_United']
  for (col in start:stop) {
    target_pop <- master_df[row, col]
    col <- colnames(master_df)[col]
    population <- substring(col, 6)
    print(population)
    ratio <- target_pop / source_pop
    m_df[(nrow(m_df) + 1),] <- c(code, trait_name, population, ratio)
  }
  rm(row, col, source_pop, target_pop, ratio, trait_name, code, population)
}

#Quick change from United to UK
m_df["ancestry"][m_df["ancestry"] == 'United'] <- 'United Kingdom' 

#Distances for each ancestry
pop_centers <- read_csv(pop_centers_file) #Population Centers Supplemental File
PC <- pop_centers %>% select(PC1:PC16)
dist_to_UK <- as.matrix(dist(PC))[,1]
distances <- tibble(
  population = pop_centers$Ancestry,
  dist_to_UK = dist_to_UK
)
#print(distances)

#Change colnames of distances 
colnames(distances) <- c('ancestry', 'pc_dist')

#Merge two dataframes
m_df <- merge(m_df, distances, by = 'ancestry')
m_df["ancestry"][m_df["ancestry"] == 'United Kingdom'] <- 'UK British' 

#Load for labeling the populations
library(ggrepel)

#Shallow
i <- "log_bilirubin"
df_subset <- m_df[m_df$prive_code==i,] #Subset
df_subset <- df_subset[!is.na(df_subset$pcor_ratio),] #drop NAs
print(df_subset)
title <- unique(df_subset$description)
df_subset$pc_dist <- as.numeric(df_subset$pc_dist) #x-coordinates
df_subset$pcor_ratio <- as.numeric(df_subset$pcor_ratio) #y-coordinates
fit <- lm((I(pcor_ratio - 1) ~ pc_dist + 0), data = df_subset) #application of linear model
slope <- round(summary(fit)$coefficients[1], digits = 5) #extract slope coefficient
title <- paste(title, slope, sep ='\nm = ') 
predictions <- predict(fit, interval="confidence", level=0.95) #confidence interval for lm
df_subset <- cbind(df_subset, predictions)
plot1 <- ggplot(data = df_subset, aes(x = pc_dist, y = pcor_ratio, label = ancestry)) + geom_point(color='black', size=1.5) + geom_ribbon(aes(ymin = lwr+1, ymax = upr+1), fill = "grey", alpha = 0.6) + geom_line(aes(y=fit+1), color='dodgerblue1') + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=0.25, fill=NA)) + xlim(0,430) + geom_text_repel(size=4) + labs(x="PC Distance from UK", y="PGS Accuracy") + geom_label(label=title, x=100, y = 0.2, colour="black", label.size= NA, size=4)
plot1 <- plot1 + scale_y_continuous(name="PGS Accuracy relative to UK British", breaks=c(0,0.2,0.4, 0.6, 0.8, 1), limits=c(0,1.1))
print(plot1)

#Steep 
i <- "haemoglobin"
df_subset <- m_df[m_df$prive_code==i,]
df_subset <- df_subset[!is.na(df_subset$pcor_ratio),]
print(df_subset)
title <- unique(df_subset$description)
df_subset$pc_dist <- as.numeric(df_subset$pc_dist)
df_subset$pcor_ratio <- as.numeric(df_subset$pcor_ratio)
fit <- lm((I(pcor_ratio - 1) ~ pc_dist + 0), data = df_subset)
slope <- round(summary(fit)$coefficients[1], digits = 5)
title <- paste(title, slope, sep ='\nm = ')
predictions <- predict(fit, interval="confidence", level=0.95)
df_subset <- cbind(df_subset, predictions)
plot2 <- ggplot(data = df_subset, aes(x = pc_dist, y = pcor_ratio, label = ancestry)) + geom_point(color='black', size=1.5) + geom_ribbon(aes(ymin = lwr+1, ymax = upr+1), fill = "grey", alpha = 0.6) + geom_line(aes(y=fit+1), color='dodgerblue1') + xlim(0,430) + geom_text_repel(size=4) + labs(x="PC Distance from UK", y="PGS Accuracy")
plot2 <- plot2 + scale_y_continuous(name="PGS Accuracy relative to UK British", breaks=c(0,0.2,0.4, 0.6, 0.8, 1), limits=c(0,1.1)) + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=0.25, fill=NA)) + geom_label(label=title, x=100, y = 0.2, colour="black", label.size= NA, size=4)
print(plot2)

#Load to arrange in grind
library(ggpubr)
g1 <- ggarrange(plot1, plot2, ncol=1, nrow=2, common.legend = FALSE)
print(g1)

#ggsave(file = "~/Desktop/Figure2.jpeg", units = c("in"), width=5, height=10, dpi=300, g1)
