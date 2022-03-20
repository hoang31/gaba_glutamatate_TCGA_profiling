
# -*- coding: utf-8 -*-

##########################################################

"""
This script will permit to extract the genes associated with a specific KEGG pathways
"""

##########################################################

##### Import library
import urllib.request
import re
import sys
import pandas as pd
import os

##########################################################

## pathways id
GABA_synapse = 'hsa04727' # GABA synapses
GLUTA_synapse = 'hsa04724' # GLUTA synapses
GLUTA_metabolism = 'hsa00250' # Alanine, aspartate and glutamate metabolism
CALCIUM_signaling = 'hsa04020' # Calcium signaling pathway
CALCIUM_endocrine = 'hsa04961' # Endocrine and other factor-regulated calcium reabsorption

## put all the pathways into a same list
kegg_pathways_id_list = [
    GABA_synapse,
    GLUTA_synapse,
    GLUTA_metabolism,
    CALCIUM_signaling,
    CALCIUM_endocrine,

    'hsa04721' # Synaptic vesicle cycle
    # 'hsa00020', # TCA cycles
    # 'hsa01230', # Biosynthesis of amino acids
    # 'hsa01200', # Carbon metabolism
    # 'hsa05230' # Central carbon metabolism in cancer

]
##########################################################


################################################################################
##### function for extrating the KEGG KGML files and informations to data frame. This function will output two outputs : df containing the node information and a df containing the relation ship information between nodes
################################################################################

def get_node_edge_data_from_KGMLfile(kegg_pathways_id) :

    ## url associated with the KEGG pathway
    url = "http://rest.kegg.jp/get/" + kegg_pathways_id + '/kgml'

    ## initialise the the node and edges directories
    kegg_network_nodes = {}
    kegg_network_edges = []

    ## extract the genes
    with urllib.request.urlopen(url) as f:
        lines = f.read().decode('utf-8').splitlines()

        ## for each line        
        for line in lines:

            ## split the line by ' '
            fields = line.split('\" ')
            
            ## reformate the fiekd removing specific character and add underscore
            fields = [field.replace('"', '') for field in fields]
            fields = [field.replace('<', '') for field in fields]
            fields = [field.replace(' ', '_') for field in fields]
            fields = [re.sub('_{2,}', '', field) for field in fields]
            fields = [re.sub('TITLE:', '', field) for field in fields]


            ########################################
            ##### For the relations between entities
            ########################################

            ## look if the line contain the entry id of a entity
            match_object = re.search(
                pattern = 'entry_id',
                string = fields[0]
            )

            ## if the line contain the entry_id, add the information into a dictionnary
            if (match_object):

                ## get the kegg_id of the entity and initialize a dic from this id
                kegg_id = str(fields[0].split('=')[1])
                kegg_network_nodes[kegg_id] = {}

                ## add the id into the dic
                kegg_network_nodes[kegg_id]['kegg_id'] = kegg_id 

                ## for each subfield contained by the field (we pass the first one corresponding to the kegg_id)
                for fields_index in range(1, len(fields), 1):

                    ## split each information of the fields
                    subfield = fields[fields_index].split('=')
                    
                    ## extract the variable name and the variable value
                    variable_type = subfield[0]
                    variable_value = str(subfield[1])

                    ## put the variable name and the variable values into the dictionary
                    kegg_network_nodes[kegg_id][variable_type] = variable_value 

            ## look if the line contain the entry id of a entity
            match_object = re.search(
                pattern = 'graphics_name',
                string = fields[0]
            )

            ## if the line contain the entry_id, add the information into a dictionnary
            if (match_object):
                ## extract the symbole of the entity by splitting by '='
                symbole_name = (fields[0].split('='))[1]
                
                ## extract the first gene name
                symbole_name = (symbole_name.split(','))[0]

                ## reformate the entity_symboles
                symbole_name = symbole_name.replace('_', '')
                symbole_name = symbole_name.replace('...', '')

                ## add the symbole
                kegg_network_nodes[kegg_id]['symbole_name'] = symbole_name 

                # ## create a unique kegg id associated with the name
                # kegg_id_name = kegg_id + '_' + symbole_name
                kegg_network_nodes[kegg_id]['kegg_id'] = symbole_name
 

            ########################################
            ##### For the relations between entities
            ########################################

            ## look the line that contain a relationship
            match_object = re.search(
                pattern = 'relation_entry1',
                string = fields[0]
            )
            
            if match_object:
                ## extract the two nodes involved in the relationship
                node1 = fields[0].split('=')[1]
                node2 = fields[1].split('=')[1]

                # ## add the name into the node name
                node1 = kegg_network_nodes[str(node1)]['kegg_id']
                node2 = kegg_network_nodes[str(node2)]['kegg_id']

            ## look the line that contain the nature of the relationship
            match_object = re.search(
                pattern = 'subtype_name',
                string = fields[0]
            )
            
            if match_object:
                # print('----- RELATION NATURE -----')

                ## extract the nature of the relation ship
                relation_nature = fields[0].split('=')[1]

                ## add the relationship into the edge list
                kegg_network_edges.append([node1, node2, relation_nature])



    

    ########################################
    ##### return data
    ########################################

    ## convert the dictionnary to panda frame
    kegg_network_nodes_df = pd.DataFrame.from_dict(
        kegg_network_nodes,
        orient = 'index',
        dtype = 'string'
    )

    ## renitialize index
    kegg_network_nodes_df = kegg_network_nodes_df.reset_index(drop=True)

    ## transform the list to data frame
    kegg_network_edges_df = pd.DataFrame(kegg_network_edges)
    
    ## rename columns
    kegg_network_edges_df = kegg_network_edges_df.rename(
        columns={
            0: "node1",
            1: "node2",
            2: "nature"
        }
    )
    kegg_network_edges_df.drop_duplicates(inplace=True)
    kegg_network_nodes_df.drop_duplicates(inplace=True)

    return([kegg_network_nodes_df, kegg_network_edges_df])


################################################################################
##### function that merge the node and edge informations generated by the previous function. This function will output two df of merged node information and merge edges informations
################################################################################

def merge_node_edges_information(list1, list2) :
    
    ## merge node data together
    merged_node_information = pd.concat(
        objs = [list1[0], list2[0]]
    )

    ## merge node data together
    merged_edge_information = pd.concat(
        objs = [list1[1], list2[1]]
    )

    ## sort the node and edge data by the 'kegg_id' column and the 'node1' column, respectively
    # merged_node_information = merged_node_information.sort_values(by = ['kegg_id'])
    merged_edge_information = merged_edge_information.sort_values(by = ['node1'])

    ## re initialize the index for the rows
    merged_node_information = merged_node_information.reset_index(drop=True)
    merged_edge_information = merged_edge_information.reset_index(drop=True)

    ## return the two merged
    return([merged_node_information, merged_edge_information])


################################################################################
##### Function for parse KGML files and merged all the data form the kegg id list
################################################################################

def get_node_edge_data_from_KGMLfile_fromList(kegg_pathways_list):

    ## extract node and edge information from the two first id pathways 
    information0 = get_node_edge_data_from_KGMLfile(kegg_pathways_list[0])

    if len(kegg_pathways_list) == 1:
        ## sort the data 
        information0[0] = information0[0].sort_values(by = ['kegg_id'])
        information0[1] = information0[1].sort_values(by = ['node1'])
        return(information0)

    if len(kegg_pathways_list) >= 2 :
        information1 = get_node_edge_data_from_KGMLfile(kegg_pathways_list[1])

        ## merge the two data from the two first kegg pathway id 
        merged_information = merge_node_edges_information(information0, information1)

    if len(kegg_pathways_list) == 2 :
        ## sort the data 
        merged_information[0] = merged_information[0].sort_values(by = ['kegg_id'])
        merged_information[1] = merged_information[1].sort_values(by = ['node1'])
        return(merged_information)

    print(kegg_pathways_list[0])        
    print(kegg_pathways_list[1])        


    ## merge with the other information from the other pathways
    for i in range(2, len(kegg_pathways_list)):

        print(kegg_pathways_list[i])        
        # print(merged_information[0].shape, merged_information[1].shape)

        ## get the node and edge information from the kgnl file form the i-th pathway id
        information3 = get_node_edge_data_from_KGMLfile(kegg_pathways_list[i])
        
        ## merge the information with the merge_information data frame
        merged_information = merge_node_edges_information(merged_information, information3)
    
    ## sort the data 
    merged_information[0] = merged_information[0].sort_values(by = ['kegg_id'])
    merged_information[1] = merged_information[1].sort_values(by = ['node1'])


    print(merged_information[0].shape, merged_information[1].shape)
    return(merged_information)


####################################################
## EXTRACTION OF DATA
####################################################

## extract all the nodes and edges information from the KGML files of each kegg pathway
kegg_network_information = get_node_edge_data_from_KGMLfile_fromList(kegg_pathways_id_list)
## attribute into two independant variable the node and the edge information
kegg_node_information = kegg_network_information[0]
kegg_edge_information = kegg_network_information[1]


####################################################
## write the snakemake output
####################################################

## write the node and edge information into csv files
kegg_node_information.to_csv(
    snakemake.output['KEGG_interaction_network_nodes'],
    index = False
)

kegg_edge_information.to_csv(
    snakemake.output['KEGG_interaction_network_edges'],
    index = False

)


