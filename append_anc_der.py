#!/usr/bin/env python3

#Example CLI Command: ./append_anc_der.py -t 6 -s ~/Desktop/Gini-PGS/top100bin_SNPs.txt -o ~/Desktop/ -d /Volumes/projects/1000_Genomes_Project/Derived_allele_frequencies/1KG_regions_derived_allele_frequency_population/
#Append the Ancestral and Derived Allele to the top100bin_SNPs.txt
#Be sure you do not have a "tmp" directory in your working directory

import os 
import argparse as ap
import pandas as pd
import multiprocessing
import time

parser = ap.ArgumentParser()
parser = ap.ArgumentParser(description="append ancestral and derived alleles to top100bin_SNPs.txt")
parser.add_argument("-t", "--threads", help="number of threads", metavar='', type=int, default=1, required=False)
parser.add_argument("-s", "--snps", help="input top100bin_SNPs file", required=True)
parser.add_argument("-d", "--derived_freqs", help="directory containing the derived allele frequencies for 1KG Genomes", required=True)
parser.add_argument("-o", "--output_directory", help="output directory", required=True)
args = parser.parse_args()

#Start time
start_time = time.time()

#Path to frequency file
freqs_path = args.derived_freqs
#Grab the working directory to generate temporary directory
wdir = os.getcwd()

#Make directory tmp in the working directory 
if os.path.isdir(wdir + '/tmp/'):
	pass
else:
	os.makedirs('tmp')

#Order the top100bin_SNPs.txt by chromosome (currently, ordered by traits)
SNP_file = pd.read_csv(args.snps, delimiter = "\t")
SNP_file.sort_values(by=['chrom'], inplace=True)

def append_anc_der(chrom_num):
	#First subset the dataframe 
	chrom_df = SNP_file[SNP_file['chrom'] == chrom_num]

	#Access the appropriate file in the SMB drive given the chrom_num and load into dataframe
	freqs = '{0}chr{1}_phase3_derived_frequency_pop.txt'.format(freqs_path, str(chrom_num))

	#Merge the file to chrom_df using the rsID
	print(f'Loading in frequncy dataframe for chromosome {chrom_num}')
	freqs_df = pd.read_csv(freqs, delimiter = '\t', usecols = ['CHROM', 'POS_HG19', 'rsID', 'ANC', 'REF', 'ALT'])

	print(f'Merging for between two dataframes for chromosome {chrom_num}')
	merged_df = chrom_df.merge(freqs_df, how='inner', on='rsID')

	#Free up memory by removing old dataframes 
	del(freqs_df)
	del(chrom_df)

	print(f'Writing merged dataframe to CSV for chromosome {chrom_num}')
	#Generate a file in the tmp directory and write to the file the merged dataframe
	merged_df.to_csv('{0}/tmp/chr{1}tmp.csv'.format(wdir, chrom_num), index = False)

	return None

if __name__ == '__main__':
	chromosomes = list(range(1,23))
	num_process = args.threads
	with multiprocessing.Pool(processes=num_process) as pool:
		results = pool.map(append_anc_der, chromosomes)
	print(results)
	pool.close()
	pool.join()

	print("Total Runtime to Append:", time.time() - start_time)
	final_cols = ['chrom', 'rsID', 'chr_position', 'A1', 'A2', 'prive_code', 'ANC', 'REF', 'ALT']
	final_df = pd.DataFrame(columns=final_cols)

	#Loop through the tmp directory and append to the final_df
	tmp_directory = '{0}/tmp/'.format(wdir)
	file_list = []
	for filename in os.listdir(tmp_directory):
		file = os.path.join(tmp_directory, filename)
		file_list.append(file)

	#Append each file to a master empty dataframe
	for tmp_file in file_list:
		tmp_df = pd.read_csv(tmp_file, delimiter = ',')
		tmp_df = tmp_df[final_cols] #Subsetting tmp_df to final columns
		final_df = final_df.append(tmp_df, ignore_index=True) #Append to the final dataframe
		os.remove(tmp_file) #Remove file with each iteration

	print(final_df)

	#Write the master dataframe to the correct output directory
	output = args.output_directory + 'top100bin_SNPs_edited.csv'
	final_df.to_csv(output, index=False)

	#Remove the tmp directory 
	os.rmdir(tmp_directory)

