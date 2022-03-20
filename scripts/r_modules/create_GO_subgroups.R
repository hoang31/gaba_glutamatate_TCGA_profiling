
########################################
## Generate subgroup of GO terms together
########################################

########################################
## Load the libraries
########################################

library(data.table)
library(tidyverse)
library(GO.db)

########################################
## Load the data
########################################

## load enrichment files that are in the correlatred gene enrichment directory 
enrichment_files <- list.files(
    snakemake@input[['dir_correlated_gene_enrichment_results']],
    full.names = T
)

## load some utils fonction
source(snakemake@params[['utils_R']])


########################################
## ANLYSIS
########################################

#####
## Load the information of the gene ontology
#####

## extract go information
go_names <- (Term(GOTERM))
go_ids <- names(go_names)

## merge the information into a same data table
go_information <- data.table(
    go_name = go_names,
    go_id = go_ids
)

#####
## process the data of the enrichment directory
#####

## extract data associated with the BP go terms
GO_BP_all_paths <- grep(
    enrichment_files,
    pattern = 'GO_BP',
    value = T
)

## read the data and put them into a list
GO_BP_all_list <- lapply(
    GO_BP_all_paths,
    function(x) fread(
        x,
        sep = ','
    )
)

## rbind all the data together
GO_BP_all_dt <- as.data.table(Reduce(
    rbind,
    GO_BP_all_list
))

## extract the go term
go_terms <- unique(unlist(GO_BP_all_dt[,ID]))


########################################
## find the levels of each gene ontology terms
########################################

### extract the go term that possess GO term id
goterm_of_interest <- as.list(GOBPANCESTOR)
goterm_of_interest <- names(goterm_of_interest)

## extract go information associated with the go term that have go term id
go_information_subset <- go_information[go_id %in% goterm_of_interest,]

## for each go term, calculate the go term level
for (i_go_id in seq(1, nrow(go_information_subset), 1)) {

    #print('--------------------')
    cat(i_go_id, "/", nrow(go_information_subset), "\n")
    go_term <- (go_information_subset[i_go_id, go_id])
    go_terms_ancestor <- as.list(GOBPANCESTOR[go_term])
    go_terms_ancestor <- Reduce(intersect, go_terms_ancestor)
    
    #go_terms_ancestor <- map_id(
    #    go_terms_ancestor,
    #    go_information,
    #    'go_id',
    #    'go_name'
    #)

    #print(go_terms_ancestor)

    ## add the level of the go term in the data table
    set(
        go_information_subset,
        i = i_go_id,
        j = 'level',
        value = length(go_terms_ancestor)
    )


}

## order the go information with the level
go_information_subset <- go_information_subset[order(level),]

# write the data
fwrite(
    x = go_information_subset,
    file = 'data/gene_ontology/gene_ontology_levels.csv',
    sep = ','
)

gene_ontology_levels <- fread(
    file = 'data/gene_ontology/gene_ontology_levels.csv',
    sep = ','
)

go_term_group <- gene_ontology_levels[level == 3, go_name]

## initialize the data table that will contain the counting for each go term groups
go_term_group_counting <- data.table(
    go_group = go_term_group,
    counting = rep(0, length(go_term_group))
)


#print(sort(go_term_group))

#test <- c(
#    'biological regulation',
#    'cellular process'
#)

#test <- c(
#    'growth',
#    #'interspecies interaction between organisms',
#    'behavior',
#    'metabolic process',
#    'biological adhesion',
#    'locomotion',
#    'signaling',
#    'developmental process',
#    'immune system process',
#    'response to stimulus',
#    'localization'
#)



########################################
## Function for creating GO term subgroup
########################################

## function that return a group a go term that are associeted with the keyword
retreive_go_term <- function(
    go_term_vector, # vector of go term
    key_word # keyword, such as "immune" that will be used for go term grouping
) {

    go_term_of_interest_vector <- c()
    go_terms_to_test <- c()

    for (i_go_terms in seq(1, length(go_term_vector), 1 )) {
    
        print('---------------------------------------------------')
        print(key_word)
        print(i_go_terms)

        ## take the go term from the go terms immune vector
        go_terms_to_test <- c(
            go_term_of_interest_vector,
            go_term_vector[i_go_terms]
        )

        ## check the go term ancestor of the go term to test
        go_terms_ancestor <- as.list(GOBPANCESTOR[go_terms_to_test])
        go_terms_ancestor <- Reduce(intersect, go_terms_ancestor)
        go_terms_ancestor <- map_id(
            go_terms_ancestor,
            go_information,
            'go_id',
            'go_name'
        )


        ## check if there is a go term ancerstor
        if (length(grep(x = go_terms_ancestor, pattern = key_word)) >= 1) {
            go_term_of_interest_vector <- c(
                go_term_of_interest_vector,
                go_term_vector[i_go_terms]
            )

            #print(go_term_of_interest_vector)
        }
    }

    return(go_term_of_interest_vector)
}

## initialize data table that contain all the go group


## initialize the data table that will contain the relation between a go term and a go group
go_child_parent <- data.table(
    go_child = character(),
    go_group = character()
)

## for each go term of interest
for (i_go_terms in seq(1, length(go_terms), 1)) {
    
    print('-------------------')
    print(i_go_terms)
    print(go_terms[i_go_terms])

    ## extract the go term ancestor
    go_terms_ancestor <- as.list(GOBPANCESTOR[go_terms[i_go_terms]])
    go_terms_ancestor <- Reduce(intersect, go_terms_ancestor)
    go_terms_ancestor <- map_id(
        go_terms_ancestor,
        go_information,
        'go_id',
        'go_name'
    )

    ## extract the go term associated with the go group
    go_of_interest <- go_terms_ancestor[go_terms_ancestor %in% go_term_group]

    ## create a data table that contain the relation of the go term child with the go_of_interest
    sub_dt <- data.table(
        go_child = rep(go_terms[i_go_terms], length(go_of_interest)),
        go_group = go_of_interest
    )

    ## merge with the go_child_parent data table
    go_child_parent <- rbind(
        go_child_parent,
        sub_dt
    )

    ## actualize the counting of the go term group counting
    go_term_group_counting[go_group %in% go_of_interest, counting := counting + 1]
}


## merge the data of the go parent with the enrichment statistical information
go_child_parent <- merge(
    go_child_parent,
    GO_BP_all_dt,
    by.x = 'go_child',
    by.y = 'ID',
    all.x = T,
    sort = F
)

## for the go group that ocntain only just one go term, put 'others'
go_term_others <- unlist(as.data.table((table(go_child_parent[,go_group])))[N < 2, V1])
go_child_parent[go_group %in% go_term_others, go_group := 'others']


########################################
## Create group of go terms
########################################

### initialize the list that will contain all the go term associated with the go term list
#go_term_group_list <- list()
##go_terms <- go_terms[1:50]

### for each go term group
#for (i_go_term_group in seq(1, length(go_term_group), 1)) {

#    ## extract all the go term associated with the go term group
#    go_terms_of_group <- retreive_go_term(
#        go_term_vector = go_terms,
#        key_word = go_term_group[i_go_term_group]
#    )

#    ## add the go terms into the list
#    go_term_group_list[go_term_group[i_go_term_group]] <- list(go_terms_of_group)

#}

### re order the list by the length
#go_term_group_list <- go_term_group_list[order(sapply(go_term_group_list, length))]

#print(go_term_group_list)

##print(GO_BP_all_dt)
#for (i_go_term_group_list in seq(1, length(go_term_group_list), 1)) {

#    ## extract the go term group
#    go_group_name <- names(go_term_group_list[i_go_term_group_list]) 

#    ## extract the go term associated with the go term group
#    go_terms_of_group <- unlist(go_term_group_list[i_go_term_group_list])

#    ## add the go group name for each go term of the group into the data table
#    GO_BP_all_dt[ID %in% go_terms_of_group, type := go_group_name]

#}

#GO_BP_all_dt[is.na(GO_BP_all_dt[, type]), type := 'others']  


#print(table(GO_BP_all_dt[, type]))

















########################################
## Write the data table into a file
########################################


#fwrite(
#    GO_BP_all_dt,
#    snakemake@output[['go_groups']],
#    sep = ','
#)

fwrite(
    go_child_parent,
    snakemake@output[['go_groups']],
    sep = ','
)