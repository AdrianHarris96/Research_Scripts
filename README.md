# Research_Scripts README
Handful of scripts, showcasing my work for the Lachance lab

Descriptions 
- `generate_slopes.R`: Calculate partial correlations for each ancestry (in relation to UK pcor) and plot against PC distance from the UK ancestry. Apply a linear model and extract coefficients to generate plots illustrating the "portability" of traits via the shallowness of the slopes. 

Calculate ANC/DER 
- `append_anc_der.py`: Append ancestral and derived state to SNPs.txt. The SNPs.txt is split by chromosome and chromosomes are processed concurrently using multiprocessing for increased efficiency. The output is appended together in SNPs_edited.csv.
- `calc_effect_anc_proportion.py`: Using the SNPs_edited.csv, determine the proportion of ancestral and derived for each trait. Compare these proportions to divergence (D stat) to identify the correlation between higher (or lower) ancestral proportion and divergence between populations for a trait.
