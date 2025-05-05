import ast
import csv
import os
import collections
from collections import Counter
from datetime import datetime

import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import pytz

from ChemicalCase import ChemicalCase

# Pathway + Timescale Extraction - This code develops a CSV file containing the shortest path between 
# x source species and y target species for each txt file present in a directory. The output is the file name, 
# chemical pathway, and timescale.

df = pd.read_csv("gckpp_EdgeList.csv")
B = nx.DiGraph()
B.add_edges_from([(df["from"][i],df["to"][i]) for i in range(0,len(df["to"]))])

directory = '/Users/lucasvalliere/Desktop/Atmospheric Chemistry Earth Data Science Work/Amazon + Beijing Emitted Species Data/Beijing_L1_2018/'
file_list = sorted(os.listdir(directory))

with open('Beijing2018ISOP_BCTEST.csv', mode='w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['filename', 'chemical_pathway', 'timescale', 'betweenness_centrality'])

    for filename in file_list:
        if filename.endswith(".txt"):
            filepath = os.path.join(directory, filename)
            amazon_surface = ChemicalCase(filepath)

            for i in range(1, 914):
                for u, v in B.in_edges("R" + str(i)):
                    B[u][v]['timescale'] = 0.5 / (amazon_surface.reaction_rates[i - 1] + 1e-20)
                for u, v in B.out_edges("R" + str(i)):
                    B[u][v]['timescale'] = 0.5 / (amazon_surface.reaction_rates[i - 1] + 1e-20)

            if nx.has_path(B, source="ISOP", target="O3"):
                o3_path = nx.shortest_path(B, source="ISOP", target="O3", weight='timescale')
                o3_timescale = nx.path_weight(B, o3_path, weight='timescale')

                subgraph = B.subgraph(o3_path)
                betweenness = nx.betweenness_centrality(subgraph, endpoints=False)
              
                path_betweenness = [betweenness.get(node, 0) for node in o3_path]

                writer.writerow([filename, o3_path, o3_timescale,  path_betweenness])
                print(f"\n{filename}:")
                print("\nChemical Pathway:", o3_path)
                print("\nTimescale value:", o3_timescale)
                print("\nBetweenness Centrality:", path_betweenness) 


# The following code block creates a dictionary to parse through and sum values of respective betweenness centrality scores

file_path= '/Users/lucasvalliere/Desktop/Atmospheric Chemistry Earth Data Science Work/Amazon + Beijing Emitted Species Data/Beijing2018ISOP_BC.csv'

species_bc_sum = defaultdict(float)
    
df = pd.read_csv(file_path)
    
for _, row in df.iterrows():
    pathway = ast.literal_eval(row["chemical_pathway"])
    centralities = ast.literal_eval(row["betweenness_centrality"])
        
    for species, centrality in zip(pathway, centralities):
        if centrality > 0:
            species_bc_sum[species] += centrality
plt.figure(figsize=(16, 9))
sns.barplot(x=list(species_bc_sum.keys()), y=list(species_bc_sum.values()), color='blue')
plt.xlabel("Species")
plt.ylabel("Sum of Betweenness Centrality Scores")
plt.title("Sum of Betweenness Centrality Scores by Species")
plt.xticks(rotation=45)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.savefig("SumBCTest.png")
plt.show()


# The following code extracts the top 5 most influential species on the basis of betweenness centrality.
# This code produces a time series visualization, showing what scores the top 5 most influential species have
# and when they attain them. 


beijing_file_path = '/Users/lucasvalliere/Desktop/Atmospheric Chemistry Earth Data Science Work/Amazon + Beijing Emitted Species Data/Beijing2018ISOP_BCTEST.csv'
beijing_data = pd.read_csv(beijing_file_path)

beijing_data['time'] = beijing_data['filename'].str.extract(r'_(\d{4})\.txt')[0]
beijing_data['date'] = beijing_data['filename'].str.extract(r'_(\d{8})')[0]

beijing_data['datetime_utc'] = pd.to_datetime(beijing_data['date'] + beijing_data['time'], format='%Y%m%d%H%M')

beijing_tz = pytz.timezone('Asia/Shanghai')
beijing_data['datetime_local'] = beijing_data['datetime_utc'].dt.tz_localize('UTC').dt.tz_convert(beijing_tz)

beijing_data['time_in_minutes'] = beijing_data['datetime_local'].dt.hour * 60 + beijing_data['datetime_local'].dt.minute

beijing_data['chemical_pathway'] = beijing_data['chemical_pathway'].apply(eval)

species_bc_sum = defaultdict(float)
for _, row in beijing_data.iterrows():
    pathway = row['chemical_pathway']
    bc_scores = eval(row['betweenness_centrality'])
    for species, bc in zip(pathway, bc_scores):
        species_bc_sum[species] += bc

most_influential_species = sorted(species_bc_sum, key=species_bc_sum.get, reverse=True)[:5]

def get_bc_score(bc_list, pathway_list, species):
    try:
        index = pathway_list.index(species)
        return bc_list[index]
    except ValueError:
        return None

for species in most_influential_species:
    beijing_data[f'BC_{species}'] = beijing_data.apply(lambda row: get_bc_score(eval(row['betweenness_centrality']), row['chemical_pathway'], species), axis=1)

sns.set(style='whitegrid')
plt.figure(figsize=(10, 6))
for species in most_influential_species:
    sns.lineplot(data=beijing_data, x='time_in_minutes', y=f'BC_{species}', marker='o', label=species, errorbar=None)

plt.xlabel('Minutes from Midnight (Local Beijing Time)')
plt.ylabel('Betweenness Centrality Score')
plt.title('Top 5 Most Influential Species Betweenness Centrality, ISOP to O3, Beijing')
plt.legend()
plt.tight_layout()
plt.savefig('Beijing_Top5_BC.png')
plt.show()

