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

with open('Beijing2018ISOP.csv', mode='w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['filename', 'chemical_pathway', 'timescale'])

    for filename in file_list:
        if filename.endswith(".txt"):
            filepath = os.path.join(directory, filename)
            amazon_surface = ChemicalCase(filepath)

            for i in range(1, 914):
                # For all edges coming in and out of a reaction node, set to half the reaction timescale
                for u, v in B.in_edges("R" + str(i)):
                    B[u][v]['timescale'] = 0.5 / (amazon_surface.reaction_rates[i - 1] + 1e-20)
                for u, v in B.out_edges("R" + str(i)):
                    B[u][v]['timescale'] = 0.5 / (amazon_surface.reaction_rates[i - 1] + 1e-20)

            if nx.has_path(B, source="ISOP", target="O3"):
                o3_path = nx.shortest_path(B, source="ISOP", target="O3", weight='timescale')
                o3_timescale = nx.path_weight(B, o3_path, weight='timescale')

                writer.writerow([filename, o3_path, o3_timescale])
                print(f"\n{filename}:")
                print("\nChemical Pathway:", o3_path)
                print("\nTimescale value:", o3_timescale)


# The following code produces several time series visualizations depiciting the average timescale per pathway across
# several files. The visualizations are sectioned by seasonality.

amazon_file_path = '/Users/lucasvalliere/Desktop/Atmospheric Chemistry Earth Data Science Work/Amazon + Beijing Emitted Species Data/Amazon2018ISOP.csv'
amazon_data = pd.read_csv(amazon_file_path)

amazon_data['time'] = amazon_data['filename'].str.extract(r'_(\d{4})\.txt')[0]
amazon_data['date'] = amazon_data['filename'].str.extract(r'_(\d{8})')[0]

amazon_data['datetime_utc'] = pd.to_datetime(amazon_data['date'] + amazon_data['time'], format='%Y%m%d%H%M')

amazon_tz = pytz.timezone('America/Manaus')
amazon_data['datetime_local'] = amazon_data['datetime_utc'].dt.tz_localize('UTC').dt.tz_convert(amazon_tz)

amazon_data['time_in_minutes'] = amazon_data['datetime_local'].dt.hour * 60 + amazon_data['datetime_local'].dt.minute

amazon_data['month'] = amazon_data['date'].str[4:6]

def categorize_season(month):
    month = int(month)
    if month in [1]:  
        return 'Winter'
    elif month in [4]:  
        return 'Spring'
    elif month in [7]:  
        return 'Summer'
    else:  
        return 'Fall'

amazon_data['season'] = amazon_data['month'].apply(categorize_season)

shortest_paths = amazon_data.pivot_table(index='chemical_pathway', columns='season', values='timescale', aggfunc=lambda x: 1 if x.any() else 0).fillna(0)

for season in seasons:
    present_reactions = shortest_paths[season][shortest_paths[season] == 1].index
    season_data = amazon_data[(amazon_data['season'] == season) & (amazon_data['chemical_pathway'].isin(present_reactions))]
    
    season_data['hour'] = season_data['datetime_local'].dt.hour

    hourly_data = season_data.groupby(['hour', 'chemical_pathway'], as_index=False)['timescale'].mean()

    plt.figure(figsize=(20, 12))
    sns.lineplot(data=hourly_data, x='hour', y='timescale', hue='chemical_pathway', marker='o')
    plt.yscale('log')
    plt.xlabel('Hour of Day (Local Amazon Time)')
    plt.ylabel('Timescale [seconds/molec/cc]')
    plt.title(f'Shortest Path Chemical Pathways, Amazon, ISOP to O3 - {season}')
    plt.legend(title="Reaction Pathway")
    plt.tight_layout()
    plt.savefig(f'Seasonal_Timescales_{season}.png', bbox_inches='tight')
    plt.xticks(ticks=range(0, 24, 1))
    plt.show()


# The following code creates a time series for all shortest paths by season for a given location. Code is
# reproducable for different locations and species (as long as information is created within the CSV code block)


amazon_file_path = '/Users/lucasvalliere/Desktop/Atmospheric Chemistry Earth Data Science Work/Amazon + Beijing Emitted Species Data/Amazon2018ISOP.csv'
amazon_data = pd.read_csv(amazon_file_path)

amazon_data['time'] = amazon_data['filename'].str.extract(r'_(\d{4})\.txt')[0]
amazon_data['date'] = amazon_data['filename'].str.extract(r'_(\d{8})')[0]

amazon_data['datetime_utc'] = pd.to_datetime(amazon_data['date'] + amazon_data['time'], format='%Y%m%d%H%M')

amazon_tz = pytz.timezone('America/Manaus')
amazon_data['datetime_local'] = amazon_data['datetime_utc'].dt.tz_localize('UTC').dt.tz_convert(amazon_tz)

amazon_data['time_in_minutes'] = amazon_data['datetime_local'].dt.hour * 60 + amazon_data['datetime_local'].dt.minute

amazon_data['month'] = amazon_data['date'].str[4:6]

def categorize_season(month):
    month = int(month)
    if month in [1]:  
        return 'Winter'
    elif month in [4]:  
        return 'Spring'
    elif month in [7]:  
        return 'Summer'
    else:  
        return 'Fall'

amazon_data['season'] = amazon_data['month'].apply(categorize_season)

sns.set(style='whitegrid')
sns.set_context("poster")
plt.figure(figsize=(10, 6))
sns.lineplot(data=amazon_data, x='time_in_minutes', y='timescale', hue='season', style='season', markers=True, dashes=False, palette='Set1')
plt.yscale('log')
plt.xlabel('Minutes from Midnight (Local Amazon Time)')
plt.ylabel('Timescale [seconds/molec/cc]')
plt.title('Time Scale of ISOP to O3')
plt.legend(title="Season")
plt.tight_layout()
plt.savefig('Amazon2018ISOP.png')
plt.show()
