
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

## create the dicitonnary
kegg_dictionnary = {
    "hsa04727" : {
        "id" :  "hsa04727",
        "kegg_pathway" : "GABA_synapse"
    },

    "hsa04724" : {
        "id" :  "hsa04724",
        "kegg_pathway" : "GLUTAMATE_synapse"
    },

    "hsa00250" : {
        "id" :  "hsa00250",
        "kegg_pathway" : "GLUTAMATE_metabolism"
    },

    "hsa04020" : {
        "id" :  "hsa04020",
        "kegg_pathway" : "CALCIUM_signaling"
    },

    "hsa04961" : {
        "id" :  "hsa04961",
        "kegg_pathway" : "CALCIUM_endocrine"
    },
}



##########################################################

##### function for extract the genes associated with a KEGG pathways from a url
def get_genes_from_KEGGpathways(kegg_pathways) :

    ## url associated with the KEGG pathway
    url = "http://rest.kegg.jp/get/" + kegg_pathways

    ## initialise the gene list associated with the KEGG pathway
    genes_list = []

    ## extract the genes
    with urllib.request.urlopen(url) as f:
        lines = f.read().decode('utf-8').splitlines()
        want = 0
        for line in lines:
            fields = line.split()

            # print("##############")
            # print(line)

            ## The list of genes starts here
            if fields[0] == 'GENE':
                want = 1
                ## The line with GENE is different
                # print(fields[2].rstrip(';'))
                genes_list.append(fields[2].rstrip(';'))

            ## We reached the next section of the file
            elif want == 1 and re.match('^\S', line):
                break

            ## We're still in the list of genes but not in the first line
            if want == 1 and len(fields)>1 and fields[0] != 'GENE':
                # print(fields[1].rstrip(';'))
                genes_list.append(fields[1].rstrip(';'))

        df = pd.DataFrame(genes_list)
        return(df)

##########################################################


## function for merging all the genes associated with each KEGG pathways
def merge_genes(KEGG_pathways_id_list) :

    ## initialize the panda frame with the first pathways
    kegg_genes_list =  [get_genes_from_KEGGpathways(KEGG_pathways_id_list[0])]

    ## extract the genes for each pathway and put them into the list
    for i in range(1, len(KEGG_pathways_id_list), 1) :
        # print(len(kegg_genes_list))
        kegg_genes_list.append(get_genes_from_KEGGpathways(KEGG_pathways_id_list[i]))
    
    ## merge all the genes into a same data frame
    dt_kegg_genes = pd.concat(kegg_genes_list)

    ## remove duplicated genes
    dt_kegg_genes.drop_duplicates(
        inplace=True
    )

    return(dt_kegg_genes)



##########################################################


## function for extracting all the KEGG gene pathways and save them in several files
def save_KEGG_genes(KEGG_pathways_id_list) :

    ## create the directory will contain the KEGG genes
    os.mkdir(snakemake.output["kegg_pathway_genes_directory"])

    ## extract the directory path
    dir_path = snakemake.output["kegg_pathway_genes_directory"]

    ## extract the genes for each pathway

    for i in range(0, len(KEGG_pathways_id_list), 1) :
        
        ## put the name of the kegg pathway to the directory path
        file_path = dir_path + "/" + kegg_dictionnary[KEGG_pathways_id_list[i]]["kegg_pathway"] + ".txt"

        ## extract the genes for each KEGG pathway
        kegg_genes = get_genes_from_KEGGpathways(KEGG_pathways_id_list[i])

        ## rename the column names
        kegg_genes.columns = ["genes"]

        kegg_genes.to_csv(
            file_path,
            sep = ",",
            header = False,
            index = None
        )
    



##########################################################

## put all the pathways into a same list
kegg_pathways_list = [
    GABA_synapse,
    GLUTA_synapse,
    GLUTA_metabolism,
    CALCIUM_signaling,
    CALCIUM_endocrine
]

## save all the genes independly
save_KEGG_genes(kegg_pathways_list)

## merge akk the kegg pathways genes
dt_kegg_genes_output = merge_genes(kegg_pathways_list)

## return the output 
dt_kegg_genes_output.to_csv(
    path_or_buf = snakemake.output["kegg_pathway_genes"],
    sep = ",",
    index = None,
    header = False
)