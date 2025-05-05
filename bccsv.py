import pandas as pd
import networkx as nx
import os

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
