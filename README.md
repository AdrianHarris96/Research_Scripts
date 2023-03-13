# Research_Scripts README
Handful of scripts, showcasing my work for the Lachance lab

Descriptions 
- `generate_slopes.R`: Calculate partial correlations for each ancestry (in relation to UK pcor) and plot against PC distance from the UK ancestry. Apply a linear model and extract coefficients to generate plots illustrating the "portability" of traits via the shallowness of the slopes. 

Calculate ANC/DER 
- `append_anc_der.py`: Append ancestral and derived state to SNPs.txt. The SNPs.txt is split by chromosome and chromosomes are processed concurrently using multiprocessing for increased efficiency. The output is appended together in SNPs_edited.csv.
- `calc_effect_anc_proportion.R`: Using the SNPs_edited.csv, determine the proportion of ancestral and derived for each trait. Compare these proportions to divergence (D stat) to identify the correlation between higher (or lower) ancestral proportion and divergence between populations for a trait.
- `calc_z_score.R`: Identify most/least divergent traits, calculate the PGS score (based on dosage) for each individual, calculate the mean of these scores, calculate z-score, and finally, group by population via (importing randomized_pops.csv containing the info on grouping)/calculate means of these groups to plot a colored point. 

Generate Example Divergence Plots 
- `divergence_plots.R`: Generates divergence plots for two example traits given the CSV of individual scores for these traits and the population ID CSV (Correlates each IID to Ancestry)