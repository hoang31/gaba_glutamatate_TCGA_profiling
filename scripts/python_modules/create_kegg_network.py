

"""
    Create the kegg metabolism pathway network
"""


#####################################################
## IMPORT THE LIBRARIES
#####################################################


import pandas as pd
import os
import networkx as nx
import matplotlib.pyplot as plt


#####################################################
## READ THE DATA
#####################################################


## read the node information
node_information = pd.read_csv(
    snakemake.input['KEGG_interaction_network_nodes'],
    sep = ','
)

## read the edge information
edge_information = pd.read_csv(
    snakemake.input['KEGG_interaction_network_edges'],
    sep = ','
)

## remove replicates for node and edge information
# edge_information.drop_duplicates(inplace=True)
# node_information.drop_duplicates(inplace=True)

# node_information = node_information.loc[node_information['type'] != 'compound']

#####################################################
## CREATE THE METABOLIC PATHWAY NETWORK
#####################################################


## initialize the graph
kegg_graph = nx.MultiGraph()


#############
## add nodes into the graph
#############

# node_list = []

# for row_index in range(0, len(node_information.index), 1) :

#     ## extract the node corresponding to the row_index-th row
#     node_to_add = node_information.iloc[[row_index]].values.tolist()
    

#     ## if the node doest not exist yet, add into the node_list
#     if node_to_add[0][0] not in node_list :
#         print(node_to_add)
#         node_list.append(node_to_add[0][0])

#         ## add the nodes into the kegg_graph
#         kegg_graph.add_node(
#             node_to_add[0][0]
#         )


print('\n\n\n')
#############
## add edges into the graph
#############

## for each row of the edge information, add the edge in the graph
for row_index in range(0, len(edge_information.index), 1) :

    ## extract edge information associated with the i-th row
    edge_to_add = edge_information.iloc[[row_index],].values.tolist()
    
    # if ((edge_to_add[0][0] in node_list) & (edge_to_add[0][1] in node_list)) :
    print(edge_to_add)

    ## add the edge into the kegg_graph between node1 and node 2
    kegg_graph.add_edge(
        edge_to_add[0][0], # node 2
        edge_to_add[0][1],  # node 2
        nature = edge_to_add[0][2] # add the nature association 
    )


# exit()





#####################################################
## VISUALIZE THE KEGG NETWORK
#####################################################


print("Nodes of graph: ")
print(kegg_graph.number_of_nodes())
print("Edges of graph: ")
print(kegg_graph.number_of_edges())
nx.draw(
    kegg_graph, 
    with_labels=True,
    font_weight='bold'
)
plt.show()

exit()