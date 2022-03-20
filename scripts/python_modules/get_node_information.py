
###########################################
## Clear the interaction data
###########################################

###########################################
## Load the libraries
###########################################


import os
import re
import numpy as np
import pandas as pd
from pandarallel import pandarallel


###########################################
## Load the data
###########################################


## load data from string database
information_edge = pd.read_csv(
    snakemake.input["information_edge"],
    sep = ",",
)

## load genes related to the immune response term ('GO:0002376')
immune_gene_dt = pd.read_csv(
    snakemake.input['immune_response_genes'],
    sep = "\t"
)

## rename the column names
immune_gene_dt.columns = ['source', "gene_name"]

## load the gene that were considered as high expressed
high_expressed_genes = pd.read_csv(
    snakemake.input["genes_highExpressed"],
    sep = ",",
)

###########################################
## Formate edge data
###########################################


## extract the node information from the edge information
information_node = [
    information_edge["node1"].to_list(),
    information_edge["node2"].to_list()
]

## merge the sublist into the same list
information_node = [full_list for i_sublist in information_node for full_list in i_sublist]

## extract only the unique data
information_node = list(set(information_node))

## transform the information node list to data frame
information_node = pd.DataFrame(information_node, columns = ["node_name"]).sort_values(by = ["node_name"])

## reset the index of the information node data frame
information_node = information_node.reset_index(drop=True)


###########################################
## Add attrribute information for the nodes
###########################################


##########
## Add the kegg information
##########

## retreive the genes related to the KEGG metabolic pathways that we have extracted
kegg_file_list = os.listdir(snakemake.params['kegg_pathway_dir'])

## keep only file that ends with .txt
kegg_file_list = [file for file in kegg_file_list if re.search(string = file, pattern = r'.txt')]

## initialize the data frame that will contain the genes for each kegg pathway
df_kegg = pd.DataFrame(columns=['gene_name', 'pathway'])

## load all the kegg file
for i in range(0, len(kegg_file_list)):
    
    ## extract file name
    file_name = kegg_file_list[i]

    ## identify the kegg metabolic pathway name
    kegg_name = file_name.split('.')[0]

    ## create the path of the file
    file_path = snakemake.params['kegg_pathway_dir'] + "/" + file_name

    ## load the file
    file = open(file_path, "r")
    gene_list = file.read().split('\n')

    ## remove the last item from the list that is a empty value
    del gene_list[-1]

    ## transform the genes list to data frame
    df_genes = pd.DataFrame(
        gene_list,
        columns=['gene_name']
    )

    ## add a column that corresponds to the kegg name
    df_genes["pathway"] = kegg_name

    ## concat the df_genes with the df_kegg data frames
    df_kegg = pd.concat(
        [
            df_kegg,
            df_genes
        ],
        axis = 0
    )



#exit()
## merge the kegg metabolic pathway informatiobn with the node data
information_node = pd.merge(
    information_node,
    df_kegg,
    how="left",
    left_on = 'node_name',
    right_on = 'gene_name'
)

## rempve redundant column
information_node = information_node.drop(["gene_name"], axis = 1)

## replace NaN by other
information_node = information_node.replace(np.nan, 'other', regex=True)


##########
## Add information related to the immune response
##########


## extract the genes related to the immune response and put them into a list
immune_gene_list = immune_gene_dt.loc[:, "gene_name"].to_list()

## put the immune information into the information node data frame
information_node.at[information_node['node_name'].isin(immune_gene_list), 'gene_ontology'] = "immune_response"

## replace the nan value to empty string
information_node = information_node.replace(np.nan, '', regex=True)

## cat the kegg pathway and the gene ontology information
information_node[['pathway']] = information_node["pathway"] + "_" + information_node["gene_ontology"]

## remove the _ at the end for every string
information_node = information_node.replace('_$','', regex=True)

## replace the other_immune_responser to immune response
information_node = information_node.replace('other_immune_response', 'immune_response', regex=True)

## remove the gene ontology colum
information_node = information_node.drop(["gene_ontology"], axis = 1)


##########
## Add the information on the genes that were filtered based on their expression
##########

## extract the genes that were highly expressed and put them into a list
high_expressed_genes_list = high_expressed_genes["gene_name"].to_list()

## put into the node information the expression
information_node.at[information_node['node_name'].isin(high_expressed_genes_list), 'pathway'] = "interest"

## replace the na value by "unknow"
information_node = information_node.replace(np.nan, 'unknown', regex=True)


###########################################
## Write the data
###########################################

## write the node information
information_node.to_csv(
    path_or_buf = snakemake.output["information_node"],
    index=False,
    sep = ","
)

















exit()

