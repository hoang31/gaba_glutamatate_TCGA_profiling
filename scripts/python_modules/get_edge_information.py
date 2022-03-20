
###########################################
## Clear the interaction data
###########################################

###########################################
## Load the libraries
###########################################

import pandas as pd
from pandarallel import pandarallel

###########################################
## Load the data
###########################################

## load data from string database
stringdb_protein_protein = pd.read_csv(
    snakemake.input["stringdb_protein_protein"],
    sep = " "
)

## load the protein information for the string database
stringdb_protein_information = pd.read_csv(
    snakemake.input["stringdb_protein_information"],
    sep = "\t"
)

## set the column names for the data
col_names = [
    "RNAInter ID",
    "Interactor1",
    "ID1",
    "Category1",
    "Species1",
    "Interactor2",
    "ID2",
    "Category2",
    "Species2",
    "Score"
]

## load data from the RNAinter database
rnainter_rna_rna = pd.read_csv(
    snakemake.input["rnainter_rna_rna"],
    sep = "\t",
    names = col_names
)

rnainter_rna_protein = pd.read_csv(
    snakemake.input["rnainter_rna_protein_filtered"],
    sep = "\t",
    names = col_names
)

rnainter_rna_compound = pd.read_csv(
    snakemake.input["rnainter_rna_compound"],
    sep = "\t",
    names = col_names
)


###########################################
## Extract the interaction data
###########################################


## set the cutoff for the combined score for the stringdb
stringdb_cutoff_score = 800

## set the cutoff for the RNAinter confidence_score
rnainter_cutoff_score = 0.5

##########
## STRINGdb
##########

## filter the interaction based on the combined score
stringdb_protein_protein = stringdb_protein_protein[stringdb_protein_protein['combined_score'] > stringdb_cutoff_score].rename(columns={"protein1" : "node1", "protein2": "node2", "combined_score" : "score"})

## extract names of the proteins and rename columns
edge_information_stringdb_prot_prot = stringdb_protein_protein[["node1", "node2"]].copy()

## add a colum that correspond to the name of the database of origin
edge_information_stringdb_prot_prot["database"] = ["stringdb"] * len(edge_information_stringdb_prot_prot.index)

## add a column that correspond the interaction type
edge_information_stringdb_prot_prot["interaction_type"] = ["protein_protein"] * len(edge_information_stringdb_prot_prot.index)

## formate the node names
edge_information_stringdb_prot_prot.loc[:,"node1"] = edge_information_stringdb_prot_prot.loc[:,"node1"].apply(lambda x: str(x).split(".")[1])
edge_information_stringdb_prot_prot.loc[:,"node2"] = edge_information_stringdb_prot_prot.loc[:,"node2"].apply(lambda x: str(x).split(".")[1])

## formate the protein id
stringdb_protein_information['protein_external_id'] = stringdb_protein_information['protein_external_id'].apply(lambda x: str(x).split(".")[1])

#edge_information_stringdb_prot_prot = edge_information_stringdb_prot_prot.iloc[1:100,]

## for the parallelization
pandarallel.initialize()

## change the protein ids
edge_information_stringdb_prot_prot["node1"] = edge_information_stringdb_prot_prot["node1"].parallel_apply(lambda x: stringdb_protein_information[stringdb_protein_information["protein_external_id"] == x].loc[:,'preferred_name'].to_string().split(' ')[-1])

edge_information_stringdb_prot_prot["node2"] = edge_information_stringdb_prot_prot["node2"].parallel_apply(lambda x: stringdb_protein_information[stringdb_protein_information["protein_external_id"] == x].loc[:,'preferred_name'].to_string().split(' ')[-1])


##########
## STRINGdb
##########

def get_edge_information_from_rnainter(
    rnainter_pandaframe, # panda frame that contain the interaction data,
    biotype, # list that contain the type of the interactor (if mRNA, miRNA, etc...)
    database_name, # string that correspond to the data base
    interaction_type, # string that correspond to the interaction type
    cutoff # cutoff used for the interaction filtering
):

    ## extract only human interaction
    rnainter_pandaframe = rnainter_pandaframe.loc[rnainter_pandaframe['Species1'] == 'Homo sapiens']
    rnainter_pandaframe = rnainter_pandaframe.loc[rnainter_pandaframe['Species2'] == 'Homo sapiens']

    ## extract only interaction of the interaction type of interest
    rnainter_pandaframe = rnainter_pandaframe[rnainter_pandaframe["Category1"].isin(biotype)]
    rnainter_pandaframe = rnainter_pandaframe[rnainter_pandaframe["Category2"].isin(biotype)]

    ## transform the score column values to numeric
    rnainter_pandaframe["Score"] = pd.to_numeric(rnainter_pandaframe["Score"])

    ## filter by the score
    rnainter_pandaframe = rnainter_pandaframe[rnainter_pandaframe["Score"] > cutoff]

    ## extract the column of interest
    rnainter_pandaframe = rnainter_pandaframe.loc[:,["Interactor1", "Interactor2"]].rename(columns={"Interactor1" : "node1", "Interactor2": "node2"})

    ## add the column related to the database
    rnainter_pandaframe["database"] = [database_name] * len(rnainter_pandaframe.index)

    ## add the column related to the interaction type
    rnainter_pandaframe["interaction_type"] = [interaction_type] * len(rnainter_pandaframe.index)

    ## reset index
    #rnainter_pandaframe = rnainter_pandaframe.reset_index(inplace=True)

    print(rnainter_pandaframe)
    ## return the table that contain the edge informations
    return(rnainter_pandaframe)


## call th function
edge_information_rnainter_rna_rna = get_edge_information_from_rnainter(
    rnainter_pandaframe = rnainter_rna_rna,
    biotype = ['mRNA', 'miRNA'],
    database_name = "RNAinter",
    interaction_type = "rna_rna",
    cutoff = rnainter_cutoff_score
)

edge_information_rnainter_rna_protein = get_edge_information_from_rnainter(
    rnainter_pandaframe = rnainter_rna_protein,
    biotype = ['mRNA', 'miRNA',"protein"],
    database_name = "RNAinter",
    interaction_type = "rna_protein",
    cutoff = rnainter_cutoff_score
)

edge_information_rnainter_rna_compound = get_edge_information_from_rnainter(
    rnainter_pandaframe = rnainter_rna_compound,
    biotype = ['mRNA', 'compound', 'protein'],
    database_name = "RNAinter",
    interaction_type = "rna_compound",
    cutoff = rnainter_cutoff_score
)


###########################################
## Write the data
###########################################

## merging edge information into a same data frame
merged_edge_information = pd.concat(
    [
        edge_information_stringdb_prot_prot,
        edge_information_rnainter_rna_rna,
        edge_information_rnainter_rna_protein,
        edge_information_rnainter_rna_compound
    ]
)

### get the biotype of node1 and node2
#merged_edge_information['biotype_node1'] = merged_edge_information.loc[:,"interaction_type"].parallel_apply(lambda x: str(x).split("_")[0])
#merged_edge_information['biotype_node2'] = merged_edge_information.loc[:,"interaction_type"].parallel_apply(lambda x: str(x).split("_")[1])

### concat the biotype of the node biotype and the node name
#merged_edge_information["node1"] = merged_edge_information["node1"] + "_" + merged_edge_information["biotype_node1"] 
#merged_edge_information["node2"] = merged_edge_information["node2"] + "_" + merged_edge_information["biotype_node2"] 

### remove irrelavant column
#merged_edge_information = merged_edge_information.drop(["biotype_node1", "biotype_node2"], axis = 1)

## write the data
merged_edge_information.to_csv(
    path_or_buf = snakemake.output["information_edge"],
    index=False,
    sep = ","
)




