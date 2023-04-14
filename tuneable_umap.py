#!/usr/bin/env python3

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import os
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from umap import UMAP
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

#Load CSV with heritabilities segmented into 100 kb chunks
df = pd.read_csv("/Users/adrianharris/Desktop/Gini-PGS/bin_heritability.csv")
df = df.loc[df['population'] == 'United'] #Filter the df to heritabilities in the UK population
df = df.drop('population', 1)
labels = df['prive_code'].tolist() #Labels for the traits
#print(labels)

X = df.drop('prive_code', 1) #Drop this prior to fitting to UMAP
#print(data)

def run_umap(dataframe, neighbors, dist): 
    fit = UMAP(
        n_neighbors=neighbors,
        min_dist=dist
    ) 
    #n_neighbors parameter balances the local vs the global structure. Lower neighbors prioritizes local clusters. 
    #min_dist parameter controls how densely points are packed together. 
    umap_array = fit.fit_transform(X)
    umap_array = umap_array.tolist()

    #Creating a output filename
    filename = '{0}x{1}.csv'.format(str(neighbors), str(dist))
    output = open(filename, 'w')
    for x in range(0, len(umap_array)):
        item = umap_array[x]
        item.insert(0, labels[x]) #Adding the labels to final output
        item = [str(i) for i in item] #converting elements in array to string
        item = ','.join(item)
        item += '\n'
        output.write(item)
    return f"completed {0}".format(filename)

out_dir = '/Users/adrianharris/Desktop/umap_folder'

#Check if the directory exists 
if os.path.exists(out_dir):
    pass
else:
    os.mkdir(out_dir)

os.chdir(out_dir)

neighbors_list = [10, 20, 50, 100, 200]
dist_list = [0.1, 0.25, 0.5, 0.8, 0.99]

for n in neighbors_list:
    for d in dist_list:
        run_umap(X, n, d)
