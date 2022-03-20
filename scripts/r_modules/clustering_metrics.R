
###################################################
###################################################


##### The scripts calcules some metrics for the clustering validation


###################################################
###################################################


##### load library
library(data.table)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)


###################################################
###################################################


# Case_ID

##### ENTROPY


## calclute the entropy of a clustering.
calulate_entropy <- function(
    ## inputs
    HEATMAP_INPUT, # heatmap generated with the function "create_heatmap2"
    EXPRESSION_DATA_INPUT, # expression data used for the clustering
    VARIABLE_DATA_INPUT, # table containing all the variables data used for the heatmap generation
    VARIABLE_INPUT, # variables used in the heatmap
    ID_SAMPLE_VARIABLE_NAMES_INPUT # name of the column of the clinical data that contains the sample id
    ) {

    ## create a data table containg all the variable for each sample
    cluster_data <- VARIABLE_DATA_INPUT[
        ,
        append(
            ID_SAMPLE_VARIABLE_NAMES_INPUT,
            VARIABLE_INPUT
        ),
        with = F
    ]
    
    ## rename the NA value by unknown
    cluster_data[, names(cluster_data) := lapply(.SD, function(x) {x[is.na(x)] <- "unknown" ; x})]

    ## extract sample id
    sample_data <- colnames(EXPRESSION_DATA_INPUT)[!(colnames(EXPRESSION_DATA_INPUT) == "Gene_Name")]

    ## extract the column cluster from the heatmap
    column_clusters <- column_order(HEATMAP_INPUT)

    if (length(column_clusters) == length(sample_data)) {
        cluster_data[
            ,
            cluster := 1
        ]
    }
    else {
        ## add the cluster number to othe cluster data
        for (cluster_n in seq(1, length(column_clusters), 1)) {
            cluster_data[
                unlist(cluster_data[, ID_SAMPLE_VARIABLE_NAMES_INPUT, with = F]) %in% sample_data[column_clusters[[cluster_n]]],
                cluster := cluster_n
            ]
        }

        # for (cluster_n in seq(1, length(column_clusters), 1)) {
        #     cluster_data[
        #         Case_ID %in% sample_data[column_clusters[[cluster_n]]],
        #         cluster := cluster_n
        #     ]
        # }
    }


    ## initialize the entropy list
    entropy_list <- list()

    ## calculate the entropy for each cluster
    for (cluster_n in seq(1, length(table(cluster_data[,cluster])), 1)) {
        # cat("########## cluster ", cluster_n, "##########", "\n")

        ## initialize the list of entropy values associated with a cluster
        entropy_cluster_list <- c()

        ## extract the data associated with the cluster
        cluster_samples <- cluster_data[cluster == cluster_n, !c(ID_SAMPLE_VARIABLE_NAMES_INPUT,"cluster"), with = F]

        # cat(nrow(cluster_samples), "samples in this cluster\n")

        ## for each variable, calculate the entropy
        for (variable in colnames(cluster_samples)) {
            # cat("---", variable, "---", "\n")
            # print((cluster_samples[, variable, with = F]))
            # print((table(cluster_samples[, variable, with = F])))

            ## apply the furmula of the entropy
            variable_entropy <- -sum(
                sapply(
                    table(cluster_samples[, variable, with = F])/nrow(cluster_samples),
                    function(x) x * log2(x)
                )
            )

            ## add the entropy value associated with the variable in the list
            entropy_cluster_list <- append(
                entropy_cluster_list,
                variable_entropy
            )
        }
        
        ## add the entropy cluster value (corresponding to the sum of all the variable entropy) to the list
        entropy_list[[paste("cluster_", cluster_n, sep = "")]]["entropy"] = sum(entropy_cluster_list)
        entropy_list[[paste("cluster_", cluster_n, sep = "")]]["weight"] = nrow(cluster_samples)/length(sample_data)

    }

    ## caculate the entropy for all the clustering
    cluster_entropy <- sum(
        sapply(
            entropy_list,
            function(x) (x[1] * x[2])
        )
    )
    # cat("----", cluster_entropy, "----\n")
    return(cluster_entropy)
}


## calculate entropy of clustering varying the final number of clusters. Return a ggplot
calculate_entropy_varying_clusterNumber <- function(
    CLUSTER_NUMBER_MIN_INPUT = 2, # number of clusters chosen for beggining the entropy calculation
    CLUSTER_NUMBER_MAX_INPUT, # last number of clusters that we will use for the entropy calculation
    EXPRESSION_DATA_INPUT, # expression input used for the entropy calculation
    VARIABLE_DATA_INPUT, # table containing all the variables data used for the heatmap
    VARIABLE_INPUT, # variables used in the heatmap
    ID_SAMPLE_VARIABLE_NAMES_INPUT, # name of the column of the clinical data that contains the sample id
    PATH_INPUT  # if there is a path, save the plot of the entropy
) {

    ## if there is no path, we do not generate the figures
    if (missing(PATH_INPUT)) {
        generate_figure <- FALSE
    }
    else {
        generate_figure <- TRUE
    }

    ## initialize the entropy data table
    entropy_data_output <- data.table(
    col_split = numeric(),
    entropy = numeric()
    )

    ## add the 0 value for the clustering, corresponding to no colulum spliting
    cluster_number_list <- append(
        0,
        seq(CLUSTER_NUMBER_MIN_INPUT, CLUSTER_NUMBER_MAX_INPUT, 1)
    )

    ## for each number of cluster, compute the entropy
    for (i in cluster_number_list) {
        ## create the heatmap
        ht <- create_heatmap2(
        expression_data = EXPRESSION_DATA_INPUT,
        clinical = VARIABLE_DATA_INPUT,
        variables = VARIABLE_INPUT,
        genes = unlist(EXPRESSION_DATA_INPUT[, Gene_Name]),
        pseudo_count = 0.1,
        # path_file = snakemake@output[["heatmap_IDHwt"]],
        log_normalization = F,
        col_split = i,
        clustering_distance_rows = "minkowski",
        clustering_method_rows = "ward.D2",
        clustering_distance_columns = "spearman",
        clustering_method_columns = "ward.D2"
        )
        
        ## calculate the entropy for the generated heatmap
        entropy_output <- calulate_entropy(
        HEATMAP_INPUT = ht,
        EXPRESSION_DATA_INPUT = EXPRESSION_DATA_INPUT,
        VARIABLE_DATA_INPUT = VARIABLE_DATA_INPUT,
        VARIABLE_INPUT = VARIABLE_INPUT,
        ID_SAMPLE_VARIABLE_NAMES_INPUT = ID_SAMPLE_VARIABLE_NAMES_INPUT
        )

        ## add all the entropy values to the data table
        if (i == 0) {
            entropy_data_output <- rbind(
                entropy_data_output,
                data.table(
                col_split = 1,
                entropy = entropy_output
                )
            )
        }
        else {
            entropy_data_output <- rbind(
                entropy_data_output,
                data.table(
                col_split = i,
                entropy = entropy_output
                )
            )
        }
        
        
    }

    ## extract the last cluster number that there is a big decrease of entropy
    count <- 0
    for (i in seq(1, nrow(entropy_data_output), 1)) {
        if (i != nrow(entropy_data_output)) {
            i0 <- unlist(entropy_data_output[i,"entropy"])
            i1 <- unlist(entropy_data_output[i + 1,"entropy"])

            entropy_diff <- i0 - i1
            cat(i, "->", i+1, "entropy diff =", entropy_diff, "\n")

            ## count the number of time where the difference of entropy are lower than 0.1, showing that the curve reach a plateau
            if (entropy_diff < 0.1) {
                count <- count + 1
            }

            if (entropy_diff > 0.1 & count < 3) {
                optimal_cluster_nb <- i+1
            }
        }
    }


    ## create the figure
    p <- ggplot(data = entropy_data_output, aes(x = col_split, y = entropy)) +
        geom_line(color = "black", size = 3, alpha = 0.8) +
        geom_point(color = "blue3", size = 6) +
        labs(
            title = "Optimal number of clusters",
            subtitle = "Entropy method",
            x = "Number of Clusters",
            y = "Entropy Average"
        ) +
        scale_x_continuous(breaks = seq(1, CLUSTER_NUMBER_MAX_INPUT, 1)) +
        theme_minimal(base_size = 30)
    print(p)

    ## if there is a PATH, create a plot of the entropy depending to the cluster number
    if (generate_figure) {
        ## save the figures
        ggsave(
            filename = PATH_INPUT,
            plot = p,
            device = "svg",
            width = 12,
            height = 8
        )
    }

    ## return the optimal cluster number
    return(optimal_cluster_nb)
}


###################################################
###################################################


##### CLUSTERING METHOD for the following function.

clustering_method <- function(
    x,
    k
) {
    clusters <- hcut(
    x = x,
    k = k,
    method = "ward.D2",
    hc_func = "hclust",
    hc_metric = "spearman"
    )
}


###################################################
###################################################


##### ELBOW METHOD

## calcule within-cluster sum of squares for the elbow methods. Return a ggplot
calculate_elbow <- function(
    EXPRESSION_DATA_INPUT, # expression input used for the entropy calculation, containing for the first column : "Gene_Name"
    CLUSTER_NUMBER_MAX_INPUT, # last number of clusters that we will use for the entropy calculation
    PATH_INPUT # if there is a path, save the plot
) {
    ## if there is a path in the input function
    if (missing(PATH_INPUT)) {
        generate_figure <- FALSE
    }
    else {
        generate_figure <- TRUE
    }

    p <- fviz_nbclust(
        x = t(EXPRESSION_DATA_INPUT[, !c("Gene_Name"), with = FALSE]),
        FUNcluster = clustering_method,
        method = "wss",
        k.max = CLUSTER_NUMBER_MAX_INPUT,
        linecolor = "blue3"
    ) +
        # geom_vline(xintercept = 3, linetype = 2) +
        labs(subtitle = "Elbow method") +
        theme_minimal()
    p

    if (generate_figure) {
        ggsave(
            filename = PATH_INPUT,
            plot = p,
            device = "jpeg"
        )
    }

    return(p)
    
}

###################################################
###################################################


##### SILHOUETTE METHOD

## calcule SILHOUETTE of the clustering. Return a ggplot
calculate_silhouette <- function(
    EXPRESSION_DATA_INPUT, # expression input used for the entropy calculation, containing for the first column : "Gene_Name"
    CLUSTER_NUMBER_MAX_INPUT, # last number of clusters that we will use for the entropy calculation
    PATH_INPUT  # if there is a path, save the plot
) {
    ## if there is a path in the input function
    if (missing(PATH_INPUT)) {
        generate_figure <- FALSE
    }
    else {
        generate_figure <- TRUE
    }

    p <- fviz_nbclust(
        x = t(EXPRESSION_DATA_INPUT[, !c("Gene_Name"), with = FALSE]),
        FUNcluster = clustering_method,
        method = "silhouette",
        k.max = CLUSTER_NUMBER_MAX_INPUT,
        linecolor = "blue3"
    ) +
        labs(subtitle = "Silhouette method") +
        theme_minimal()
    p

    if (generate_figure) {
        ggsave(
            filename = PATH_INPUT,
            plot = p,
            device = "jpeg"
        )
    }

    return(p)
}

###################################################
###################################################

##### GAB STATISTIC METHOD

## calcule GAP STATISTIC of the clustering. Return a ggplot
calculate_gapstatistic <- function(
    EXPRESSION_DATA_INPUT, # expression input used for the entropy calculation, containing for the first column : "Gene_Name"
    CLUSTER_NUMBER_MAX_INPUT, # last number of clusters that we will use for the entropy calculation
    PATH_INPUT  # if there is a path, save the plot
) {
    ## if there is a path in the input function
    if (missing(PATH_INPUT)) {
        generate_figure <- FALSE
    }
    else {
        generate_figure <- TRUE
    }

    ## copy the input
    EXPRESSION_DATA_INPUT = copy(EXPRESSION_DATA_INPUT)

    p <- fviz_nbclust(
        x = t(EXPRESSION_DATA_INPUT[, !c("Gene_Name"), with = FALSE]),
        FUNcluster = clustering_method,
        method = "gap_stat",
        k.max = CLUSTER_NUMBER_MAX_INPUT,
        nboot = 500,
        linecolor = "blue3"
    ) +
        labs(subtitle = "Gap statistic method") +
        theme_minimal()
    p

    if (generate_figure) {
        ggsave(
            filename = PATH_INPUT,
            plot = p,
            device = "jpeg"
        )
    }

    return(p)
}

###################################################
###################################################

##### NBCLUST for finding the optimal number of clustering

calculate_nbclust <- function(
    EXPRESSION_DATA_INPUT, # expression input used for the entropy calculation, containing for the first column : "Gene_Name"
    CLUSTER_NUMBER_MAX_INPUT, # last number of clusters that we will use for the entropy calculation
    INDEX_INPUT, # indice to compute
    WHERE_INPUT, # choose between columns or rows
    DISTANCE_INPUT, # choose the method for the correlation matrix
    PATH_INPUT  # if there is a path, save the plot
) {
    ## if there is a path in the input function
    if (missing(PATH_INPUT)) {
        generate_figure <- FALSE
    }
    else {
        generate_figure <- TRUE
    }

    ## copy the input
    EXPRESSION_DATA_INPUT = copy(EXPRESSION_DATA_INPUT)

    ## keep the gene name into a vector
    gene_list <- unlist(EXPRESSION_DATA_INPUT[,"Gene_Name"])
    names(gene_list) <- NULL # remive "gene"

    
    ## transform the data table to data frame and put the gene name as row names in the expresson data
    setDF(EXPRESSION_DATA_INPUT, rownames = gene_list)

    ## remove the "Gene_Name" column from the expression data
    EXPRESSION_DATA_INPUT <- subset(EXPRESSION_DATA_INPUT, select = -Gene_Name)

    ## depending to the WHERE_INPUT
    if (WHERE_INPUT == "column") {
        EXPRESSION_DATA_INPUT = EXPRESSION_DATA_INPUT
    }
    else if (WHERE_INPUT == "row") {
        EXPRESSION_DATA_INPUT = t(EXPRESSION_DATA_INPUT)
    }
    else {
        stop("Choose between row or column !")
    }
    
    ## if for we want to look the columns
    if (WHERE_INPUT == "column") {
        corr_matrix <- as.dist(
            1 - cor(
                (EXPRESSION_DATA_INPUT),
                method = DISTANCE_INPUT
            ),
        )
        # print(WHERE_INPUT)
    }
    
    ## if for we want to look the rows
    if (WHERE_INPUT == "row") {
        corr_matrix <- dist(
            t(EXPRESSION_DATA_INPUT),
            method = DISTANCE_INPUT
        )
    }

    
    ## compute the indice
    p <- NbClust(
        data = t(EXPRESSION_DATA_INPUT),
        diss = corr_matrix,
        distance = NULL,
        method = "ward.D2",
        min.nc = 2,
        max.nc = CLUSTER_NUMBER_MAX_INPUT,
        index = INDEX_INPUT
    )
    return(p)
}

## function for calculate the indexes from a vector. Return a ggplot.
calculate_nbclust_for_indexes <- function(
    EXPRESSION_DATA_INPUT, # expression input used for the entropy calculation, containing for the first column : "Gene_Name"
    INDEX_VECTOR_INPUT, # vector containing the index that we want to compute
    SUBTITLE_INPUT, # string corresponding to the sub title name of the ggplot
    WHERE_INPUT, # choose between columns or rows
    DISTANCE_INPUT, #choose the method for the distance matrix
    PATH_INPUT  # if there is a path, save the plot
) {
    ## if there is a path
    if (missing(PATH_INPUT)) {
        generate_figure <- FALSE
    }
    else {
        generate_figure <- TRUE
    }

    ## initialize the results output which contains all the results of each index
    results_output <- list()

    ## compute for all the index in the index vector input
    for (index in seq(1, length(INDEX_VECTOR_INPUT), 1)) {
        
        cat("------- INDEX :", INDEX_VECTOR_INPUT[index], "-------\n")

        ## compute for the index index-th
        result <- try(calculate_nbclust(
            EXPRESSION_DATA_INPUT = EXPRESSION_DATA_INPUT,
            CLUSTER_NUMBER_MAX_INPUT = 20,
            WHERE_INPUT = WHERE_INPUT,
            INDEX_INPUT = INDEX_VECTOR_INPUT[index],
            DISTANCE_INPUT = DISTANCE_INPUT
        ), silent = TRUE)

        if(inherits(result, "try-error"))
        {
            warning("calculation of the ",INDEX_VECTOR_INPUT[index], "index is not possible")
            next
        }

        ## put the results in the result list
        results_output[[INDEX_VECTOR_INPUT[index]]] <- result
    }
  

    ## extract the best number of cluster for each indice
    results_output <- as.data.table(
        t(sapply(results_output, function(x) {
            return(x$Best.nc)
        })),
        keep.rownames = T
    )[, setnames(.SD, "rn", "index")]

    ## remove the -inf cluster keeping (i do not know why it exists)
    results_output <- results_output[Number_clusters > 0, ]

    ## count the frequency of each indices
    indice_frequencies <- as.data.table(table(results_output[, Number_clusters]))[, setnames(.SD, "V1", "Number_clusters")]

    indice_frequencies[N == max(N), n_clusters := "optimal"][N != max(N), n_clusters := "no_optimal"]
    # print(indice_frequencies)

    ## create the plot
    p <- ggplot(
        data = indice_frequencies,
        aes(
            x = Number_clusters,
            y = N,
            fill = n_clusters
        )
    ) +
        geom_bar(
            stat = "identity",
            color = "black",
            size = 1.5
    ) +
        labs(
            title = "Optimal Number of Clusters",
            subtitle = SUBTITLE_INPUT,
            x = "Number of clusters",
            y = "Frequency among all indices"
    ) +
        scale_fill_manual("Legend", values = c("no_optimal" = "blue3", "optimal" = "red3")) +
        theme_minimal(
            base_size = 18
        ) + 
        theme(legend.position="none")

    print(p)

    if (generate_figure) {
        ggsave(
            filename = PATH_INPUT,
            plot = p,
            device = "svg"
        )
    }


    ## return the results ouput
    return(as.numeric(indice_frequencies[N == max(N), Number_clusters]))
}