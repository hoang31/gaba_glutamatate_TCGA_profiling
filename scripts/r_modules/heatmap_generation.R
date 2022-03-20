
###################################################
###################################################


###### Create heatmap from the expressiom data


###################################################
###################################################


##### input for the function :
## expression data, with in row corresponding to X genes ; in column corresponding to Y samples +1 for the column "Gene_Name"
## clinical data in row corresponding to each Y samples and N variables in columns
## variable vector corresponding to the variables shown in the heatmap
## gene vector corresponding to the genes used for the clustering for in the heatmap
## pseudo count used

#### input (optionnal)
## path of the files
## col_split, row_split corresponding to the number of groups in columns or in rows, respectively


###################################################
###################################################


## load packages
library(data.table)
library(ComplexHeatmap)
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))

###################################################
###################################################


##### function for correlation
correlation_dist <- function(mat) {
    output <- as.dist(1 - cor(t(mat), method = "spearman"))
    return(output)
}

###################################################

##### function for minkowski distance
minkowski_dist <- function(mat) {
    output <- dist(mat, method = "minkowski")
    return(output)
}

###################################################

##### Function for chosing the heatmap annotation colors for each variables (for be more clear)
create_annotation <- function(variable, clinicalData,
horizontal_legend = F) {
    ## create the colors
    all_colors <- c(
        brewer.pal(9, "Set1"),
        brewer.pal(8, "Dark2"),
        brewer.pal(8, "Set2"),
        brewer.pal(12, "Paired")
    )
    # print(all_colors)

    ## extract data table with the specific variables
    clinicalData <- clinicalData[, variable, with = F]

    ## initialise the color list
    colors_for_variables <- list()

    ## for each data associated with each variables, do a counting of occurence and extract colors
    for (i in (1:length(variable)))
    {

        ## extract the number of catagories in the variable i-th
        categories <- names(table(clinicalData[, variable[i], with = F]))

        ## take the i-th first colors from the "all_colors" variable
        colors_for_categories <- all_colors[1:length(categories)]

        ## remove the colors used from the initial color vector
        all_colors <- all_colors[!(all_colors %in% colors_for_categories)]
 

        ## name each color by the category
        names(colors_for_categories) <- categories

        ## append into the list the colors categories
        colors_for_variables[[i]] <- colors_for_categories
    }

    ## give variable names for each vector in the list
    names(colors_for_variables) <- variable

    #return(list(clinicalData, colors_for_variables))

    if (horizontal_legend == FALSE) {
        ## create the annotations
        annotation <- HeatmapAnnotation(
            df = clinicalData,
            col = colors_for_variables
        )
    }

    if (horizontal_legend == TRUE) {
    ## create the annotations
        annotation <- HeatmapAnnotation(
            df = clinicalData,
            col = colors_for_variables,
            annotation_legend_param = list(
                nrow = 1
            )
        )
    }


    return(annotation)
}

###################################################

create_annotation_row <- function(variable, clinicalData) {
    ## create the colors
    all_colors <- c(
        brewer.pal(9, "Set1"),
        brewer.pal(12, "Paired"),
        brewer.pal(8, "Set2"),
        brewer.pal(8, "Dark2")
    )
    # print(all_colors)

    ## extract data table with the specific variables
    clinicalData <- clinicalData[, variable, with = F]

    ## initialise the color list
    colors_for_variables <- list()

    ## for each data associated with each variables, do a counting of occurence and extract colors
    for (i in (1:length(variable)))
    {

        ## extract the number of catagories in the variable i-th
        categories <- names(table(clinicalData[, variable[i], with = F]))

        ## take the i-th first colors from the "all_colors" variable
        colors_for_categories <- all_colors[1:length(categories)]

        ## remove the colors used from the initial color vector
        all_colors <- all_colors[!(all_colors %in% colors_for_categories)]
 

        ## name each color by the category
        names(colors_for_categories) <- categories

        ## append into the list the colors categories
        colors_for_variables[[i]] <- colors_for_categories
    }

    ## give variable names for each vector in the list
    names(colors_for_variables) <- variable

    ## create the annotations
    annotation <- rowAnnotation(
        df = clinicalData,
        col = colors_for_variables
    )
    return(annotation)
}

###################################################

##### function for creating a heatmap
create_heatmap2 <- function(expression_data,
                            clinical,
                            variables,
                            genes,
                            pseudo_count,
                            path_file,
                            col_split,
                            row_split,
                            return_matrix,
                            log_normalization = T,
                            clustering_distance_rows = "euclidean",
                            clustering_method_rows = "complete",
                            clustering_distance_columns = "euclidean",
                            clustering_method_columns = "complete") {

    ## generate the heatmap file
    generation_heatmap <- FALSE

    ## return the count matrix ?
    if (missing(return_matrix)) {
        return_matrix <- FALSE
    }

    if (!(missing(path_file))) {
        generation_heatmap <- TRUE
    }

    ## boolean for splitting or not the rows and columns with the clustering
    split_col <- FALSE
    split_row <- FALSE

    ## if the argument col_split or row_split are present
    if (!(missing(col_split))) {
        split_col <- TRUE
    }

    if (!(missing(row_split))) {
        split_row <- TRUE
    }

    ## if there is a column spliting
    if (col_split == 0) {
        split_col <- FALSE
    }

    if (split_row == 0) {
        split_row <- FALSE
    }

    ## extract the expression related to the specific genes
    expression_filtered <- expression_data[Gene_Name %in% genes, ]

    ## transform the data table to matrix removing the "Gene_Name" column
    expression_filtered_matrix <- as.matrix(
        x = expression_filtered[, !c("Gene_Name")],
        rownames = expression_filtered[, Gene_Name]
    )

    ## if we return the filtered matrix
    if (return_matrix == TRUE) {
        return(expression_filtered_matrix)
    }

    if (log_normalization == TRUE) {
        # print("log normalization...")
        ## transform the matrix with a log + pseudocount
        expression_filtered_matrix_log <- apply(expression_filtered_matrix, c(1, 2), function(x) log2(x + pseudo_count))
    }

    if (log_normalization == FALSE) {
        # print("no log normalization...")
        ## transform the matrix with a log + pseudocount
        expression_filtered_matrix_log <- copy(expression_filtered_matrix)
    }

    ## extract and make the annotations matrix
    annotation <- create_annotation(variables, clinical)

    ## write the file
    if (generation_heatmap) { svg(path_file, width = 10, height = 8) } ## write the file

    ##### params

    ## color
    # cat(min(expression_filtered_matrix),"- ", median(expression_filtered_matrix),"- ", max(expression_filtered_matrix), "\n")
    ht_colors <- colorRamp2(
        c(min(expression_filtered_matrix), mean(expression_filtered_matrix), max(expression_filtered_matrix)),
        c("#ffffff", "#eeff00", "#ff0000")
    )
    
    ht_colors <- colorRamp2(
        c(4, 15),
        c("#000000", "#ff0000")
    )

    ht_colors <- colorRamp2(
        c(6, 9, 14),
        c("blue3", "white", "red3")
    )


    column_title <- "SAMPLES"
    row_title <- "GENES"
    column_title_gp <- gpar(fontsize = 22, fontface = "bold")
    row_title_gp <- gpar(fontsize = 22, fontface = "bold")

    ## make the heatmap with the row and column splitting
    if ((split_col == T) & (split_row == T)) {
        # print("split rows and columns")
        ht <- (Heatmap(expression_filtered_matrix_log,
            col = ht_colors,
            column_title = column_title,
            row_title = row_title,
            column_title_gp = column_title_gp,
            row_title_gp = row_title_gp,

            column_gap = unit(5, "mm"),
            row_gap = unit(5, "mm"),

            border = TRUE,
            show_column_names = FALSE,
            show_row_names = FALSE,
            column_split = col_split,
            row_split = row_split,

            top_annotation = annotation,
            heatmap_legend_param = list(
                title = "EXPRESSION LEVELS",
                title_position = "leftcenter-rot",
                at = c(0, 5, 10, 15),
                grid_width = unit(0.75, "cm"),
                legend_height = unit(6, "cm"),
                border = "black"
            ),

            clustering_distance_rows = clustering_distance_rows,
            clustering_method_rows = clustering_method_rows,
            clustering_distance_columns = clustering_distance_columns,
            clustering_method_columns = clustering_method_columns,

            row_dend_reorder = FALSE,
            column_dend_reorder = TRUE
        )
        )
    }

    ## make the heatmap with the column splitting
    if ((split_col == T) & (split_row == F)) {
        # print("split columns")
        ht <- (Heatmap((expression_filtered_matrix_log),
            col = ht_colors,
            column_title = column_title,
            row_title = row_title,
            column_title_gp = column_title_gp,
            row_title_gp = row_title_gp,


            #rect_gp = gpar(type = "none"),



            border = TRUE,
            show_column_names = FALSE,
            show_row_names = FALSE,
            column_split = col_split,
            top_annotation = annotation,
            heatmap_legend_param = list(
                title = "EXPRESSION LEVELS",
                title_position = "leftcenter-rot",
                at = c(0, 5, 10, 15),
                grid_width = unit(0.75, "cm"),
                legend_height = unit(6, "cm"),
                border = "black"
            ),

            column_gap = unit(5, "mm"),

            clustering_distance_rows = clustering_distance_rows,
            clustering_method_rows = clustering_method_rows,
            clustering_distance_columns = clustering_distance_columns,
            clustering_method_columns = clustering_method_columns,

            row_dend_reorder = FALSE,
            column_dend_reorder = TRUE
        )
        )
    }

    ## make the heatmap with the row splitting
    if ((split_col == F) & (split_row == T)) {
        # print("split rows")
        ht <- (Heatmap((expression_filtered_matrix_log),
            border = TRUE,
            col = ht_colors,
            column_title = column_title,
            row_title = row_title,
            column_title_gp = column_title_gp,
            row_title_gp = row_title_gp,

            show_column_names = FALSE,
            show_row_names = FALSE,

            row_split = row_split,
            row_gap = unit(5, "mm"),
            top_annotation = annotation,
            heatmap_legend_param = list(
                title = "EXPRESSION LEVELS",
                title_position = "leftcenter-rot",
                at = c(0, 5, 10, 15),
                grid_width = unit(0.75, "cm"),
                legend_height = unit(6, "cm"),
                border = "black"
            ),

            clustering_distance_rows = clustering_distance_rows,
            clustering_method_rows = clustering_method_rows,
            clustering_distance_columns = clustering_distance_columns,
            clustering_method_columns = clustering_method_columns,

            row_dend_reorder = FALSE,
            column_dend_reorder = TRUE
        )
        )
    }

    ## make the heatmap without the splitting
    if ((split_col == F) & (split_row == F)) {
        # print("no split rows and columns")
        ht <- (Heatmap((expression_filtered_matrix_log),
            border = TRUE,
            col = ht_colors,
            column_title = column_title,
            row_title = row_title,
            column_title_gp = column_title_gp,
            row_title_gp = row_title_gp,

            show_column_names = FALSE,
            show_row_names = FALSE,
            

            ## legends for the heatmap
            top_annotation = annotation,
            heatmap_legend_param = list(
                title = "EXPRESSION LEVELS",
                title_position = "leftcenter-rot",
                at = c(0, 5, 10, 15),
                grid_width = unit(0.75, "cm"),
                legend_height = unit(6, "cm"),
                border = "black"
            ),

            clustering_distance_rows = clustering_distance_rows,
            clustering_method_rows = clustering_method_rows,
            clustering_distance_columns = clustering_distance_columns,
            clustering_method_columns = clustering_method_columns,

            row_dend_reorder = F,
            column_dend_reorder = F
        )
        )
    }


    if (generation_heatmap) {
        print(ht)
        dev.off()
    }

    return(ht)

    print("###############################################################")
}

###################################################
###################################################




create_heatmap3 <- function(expression_data,
                            clinical,
                            variables,
                            genes,
                            pseudo_count,
                            path_file,
                            col_split,
                            row_split,
                            return_matrix,
                            log_normalization = T,
                            clustering_distance_rows = "euclidean",
                            clustering_method_rows = "complete",
                            clustering_distance_columns = "euclidean",
                            clustering_method_columns = "complete",
                            cluster_name
) {

    ## generate the heatmap file
    generation_heatmap <- FALSE

    ## return the count matrix ?
    if (missing(return_matrix)) {
        return_matrix <- FALSE
    }

    if (!(missing(path_file))) {
        generation_heatmap <- TRUE
    }

    ## boolean for splitting or not the rows and columns with the clustering
    split_col <- FALSE
    split_row <- FALSE

    ## if the argument col_split or row_split are present
    if (!(missing(col_split))) {
        split_col <- TRUE
    }

    if (!(missing(row_split))) {
        split_row <- TRUE
    }

    ## if there is a column spliting
    if (col_split == 0) {
        split_col <- FALSE
    }

    if (split_row == 0) {
        split_row <- FALSE
    }

    ## extract the expression related to the specific genes
    expression_filtered <- expression_data[Gene_Name %in% genes, ]

    ## transform the data table to matrix removing the "Gene_Name" column
    expression_filtered_matrix <- as.matrix(
        x = expression_filtered[, !c("Gene_Name")],
        rownames = expression_filtered[, Gene_Name]
    )

    ## if we return the filtered matrix
    if (return_matrix == TRUE) {
        return(expression_filtered_matrix)
    }

    if (log_normalization == TRUE) {
        # print("log normalization...")
        ## transform the matrix with a log + pseudocount
        expression_filtered_matrix_log <- apply(expression_filtered_matrix, c(1, 2), function(x) log2(x + pseudo_count))
    }

    if (log_normalization == FALSE) {
        # print("no log normalization...")
        ## transform the matrix with a log + pseudocount
        expression_filtered_matrix_log <- copy(expression_filtered_matrix)
    }

    # expression_filtered_matrix_log <- t(expression_filtered_matrix_log)

    ## extract and make the annotations matrix
    annotation <- create_annotation(variables, clinical)

    ## save the heatmap
    # path_file = paste("/home/hoang/Desktop/", name_file, sep = "") ## path of the file
    # if (generation_heatmap) { jpeg(path_file, width = 1250, height = 1050) }
    ## write the file
    if (generation_heatmap) { svg(path_file, width = 12, height = 10) } ## write the file



    # print(names(hc1))
    
    # print(hc1$order)
    ##### params

    ## color
    # cat(min(expression_filtered_matrix),"- ", median(expression_filtered_matrix),"- ", max(expression_filtered_matrix), "\n")
    ht_colors <- colorRamp2(
        c(min(expression_filtered_matrix), mean(expression_filtered_matrix), max(expression_filtered_matrix)),
        c("#000000", "#eeff00", "#ff0000")
    )
    
    ht_colors <- colorRamp2(
        c(0, 7, 15),
        c("#ffffff", "#eeff00", "#ff0000")
    )

    ht_colors <- colorRamp2(
        c(0, 9, 15),
        c("blue3", "white", "red3")
    )



    column_title <- "SAMPLES"
    row_title <- "GENES"
    column_title_gp <- gpar(fontsize = 18, fontface = "bold")
    row_title_gp <- gpar(fontsize = 18, fontface = "bold")



    ## make the heatmap without the splitting
 
    # print("no split rows and columns")
    ht <- (Heatmap((expression_filtered_matrix_log),
        border = TRUE,
        col = ht_colors,
        # column_title = column_title,
        # row_title = row_title,
        column_title_gp = column_title_gp,
        row_title_gp = row_title_gp,

        show_column_names = FALSE,
        show_row_names = FALSE,
        
        ## legends for the heatmap
        top_annotation = annotation,
        heatmap_legend_param = list(
            title = "EXPRESSION LEVELS",
            title_position = "leftcenter-rot",
            at = c(0, 5, 10, 15),
            grid_width = unit(0.75, "cm"),
            legend_height = unit(6, "cm"),
            border = "black"
        ),

        clustering_distance_rows = clustering_distance_rows,
        clustering_method_rows = clustering_method_rows,
        clustering_distance_columns = clustering_distance_columns,
        clustering_method_columns = clustering_method_columns,

        column_gap = unit(5, "mm"),
        column_split = cluster_name,
        cluster_column_slices = FALSE,
        column_title = c(
            "IDH_MUT\ncodel",
            "IDH_MUT\nnoncodel",
            # "IDH_WT\nMIXTED",
            "IDH_WT",
            "MIXTED"
        )



    )
    )


    if (generation_heatmap) {
        draw(
            ht, 
            row_title = "GENES",
            column_title = "SAMPLES",
            column_title_gp = gpar(fontsize = 28, fontface = "bold"),
            row_title_gp = gpar(fontsize = 28, fontface = "bold")
        )
        dev.off()
    }

    return(ht)

    print("###############################################################")
}






######################################################################

##### function for extract clusters generated with the generation of the heatmap

extract_clusters <- function(
    HEATMAP_OBJECT_INPUT, # object generated with the the function "heatmap2"
    EXPRESSION_DATA_INPUT # expression data used for the generation of the heatmap
) {

    ## extract the clustering groups of the heatmap
    clustering_groups <- column_order(HEATMAP_OBJECT_INPUT)

    ## extract the samples from the expression data
    samples_groups <- as.data.table(
        list(colnames(EXPRESSION_DATA_INPUT)[colnames(EXPRESSION_DATA_INPUT) != "Gene_Name"])
    )[
        ,
        setnames(
            .SD,
            "V1",
            "sample_id"
        )
    ]

    ## for each clustering groups, add the group number to the "samples_groups" data frames
    for (i in seq(1, length(clustering_groups), 1)) {
        # clustering_groups[clustering_groups[[i]], cluster := i]

        set(
            samples_groups,
            i = clustering_groups[[i]],
            j = "cluster",
            value = toString(i)
        )
    }

    return(samples_groups)

}



##### function for creating a heatmap with sample and gene annotations
create_heatmap_annotation <- function(expression_data,
                            clinical,
                            variables,
                            row_variables = 'missing',
                            genes,
                            pseudo_count,
                            path_file,
                            col_split,
                            row_split,
                            return_matrix,
                            log_normalization = T,
                            clustering_distance_rows = "euclidean",
                            clustering_method_rows = "complete",
                            clustering_distance_columns = "euclidean",
                            clustering_method_columns = "complete",
                            horizontal_legend = T
) {

    ## generate the heatmap file
    generation_heatmap <- FALSE

    ## return the count matrix ?
    if (missing(return_matrix)) {
        return_matrix <- FALSE
    }

    if (!(missing(path_file))) {
        generation_heatmap <- TRUE
    }

    ## boolean for splitting or not the rows and columns with the clustering
    split_col <- FALSE
    split_row <- FALSE

    ## if the argument col_split or row_split are present
    if (!(missing(col_split))) {
        split_col <- TRUE
    }

    if (!(missing(row_split))) {
        split_row <- TRUE
    }

    ## if there is a column spliting
    if (col_split == 0) {
        split_col <- FALSE
    }

    if (split_row == 0) {
        split_row <- FALSE
    }

    ## extract the expression related to the specific genes
    expression_filtered <- expression_data[Gene_Name %in% genes, ]

    ## transform the data table to matrix removing the "Gene_Name" column
    expression_filtered_matrix <- as.matrix(
        x = expression_filtered[, !c("Gene_Name")],
        rownames = expression_filtered[, Gene_Name]
    )

    ## if we return the filtered matrix
    if (return_matrix == TRUE) {
        return(expression_filtered_matrix)
    }

    if (log_normalization == TRUE) {
        # print("log normalization...")
        ## transform the matrix with a log + pseudocount
        expression_filtered_matrix_log <- apply(expression_filtered_matrix, c(1, 2), function(x) log2(x + pseudo_count))
    }

    if (log_normalization == FALSE) {
        # print("no log normalization...")
        ## transform the matrix with a log + pseudocount
        expression_filtered_matrix_log <- copy(expression_filtered_matrix)
    }


    ## extract and make the annotations matrix
    annotation <- create_annotation(variables, clinical, horizontal_legend = horizontal_legend)

    if (row_variables != 'missing') {
        row_annotation <- create_annotation_row('status', row_variables)
    }
    

    ## save the heatmap
    # path_file = paste("/home/hoang/Desktop/", name_file, sep = "") ## path of the file
    # if (generation_heatmap) { jpeg(path_file, width = 1050, height = 1050) }
    ## write the file
    if (generation_heatmap) { svg(path_file, width = 10, height = 8) } ## write the file

    ##### params

    ## color
    # cat(min(expression_filtered_matrix),"- ", median(expression_filtered_matrix),"- ", max(expression_filtered_matrix), "\n")
    ht_colors <- colorRamp2(
        c(min(expression_filtered_matrix), mean(expression_filtered_matrix), max(expression_filtered_matrix)),
        c("#ffffff", "#eeff00", "#ff0000")
    )
    
    ht_colors <- colorRamp2(
        c(4, 15),
        c("#000000", "#ff0000")
    )

    ht_colors <- colorRamp2(
        c(6, 9, 14),
        c("blue3", "white", "red3")
    )

    #ht_colors <- colorRamp2(
    #    c(0, 9, 20),
    #    c("blue3", "white", "red3")
    #)



    column_title <- "SAMPLES"
    row_title <- "GENES"
    column_title_gp <- gpar(fontsize = 22, fontface = "bold")
    row_title_gp <- gpar(fontsize = 22, fontface = "bold")

    ## make the heatmap with the row and column splitting
    if ((split_col == T) & (split_row == T) & (row_variables != 'missing')) {
        # print("split rows and columns")
        ht <- (Heatmap(expression_filtered_matrix_log,
            col = ht_colors,
            column_title = column_title,
            row_title = row_title,
            column_title_gp = column_title_gp,
            row_title_gp = row_title_gp,

            column_gap = unit(5, "mm"),
            row_gap = unit(5, "mm"),

            border = TRUE,
            show_column_names = FALSE,
            show_row_names = FALSE,
            column_split = col_split,
            row_split = row_split,

            top_annotation = annotation,
            left_annotation = row_annotation,

            heatmap_legend_param = list(
                title = "EXPRESSION LEVELS",
                title_position = "leftcenter-rot",
                at = c(0, 5, 10, 15),
                grid_width = unit(0.75, "cm"),
                legend_height = unit(6, "cm"),
                border = "black"
            ),

            clustering_distance_rows = clustering_distance_rows,
            clustering_method_rows = clustering_method_rows,
            clustering_distance_columns = clustering_distance_columns,
            clustering_method_columns = clustering_method_columns,

            row_dend_reorder = FALSE,
            column_dend_reorder = TRUE
        )
        )
    }

    if ((split_col == T) & (split_row == T) & (row_variables == 'missing')) {
        # print("split rows and columns")
        ht <- (Heatmap(expression_filtered_matrix_log,
            col = ht_colors,
            column_title = column_title,
            row_title = row_title,
            column_title_gp = column_title_gp,
            row_title_gp = row_title_gp,

            column_gap = unit(5, "mm"),
            row_gap = unit(5, "mm"),

            rect_gp = gpar(type = "none"),

            border = TRUE,
            show_column_names = FALSE,
            show_row_names = FALSE,
            column_split = col_split,
            row_split = row_split,

            top_annotation = annotation,
            #left_annotation = row_annotation,

            heatmap_legend_param = list(
                title = "EXPRESSION LEVELS",
                title_position = "leftcenter-rot",
                at = c(0, 5, 10, 15),
                grid_width = unit(0.75, "cm"),
                legend_height = unit(6, "cm"),
                border = "black"

            ),

            clustering_distance_rows = clustering_distance_rows,
            clustering_method_rows = clustering_method_rows,
            clustering_distance_columns = clustering_distance_columns,
            clustering_method_columns = clustering_method_columns,

            row_dend_reorder = FALSE,
            column_dend_reorder = TRUE
        )
        )
    }

    ## make the heatmap with the column splitting
    if ((split_col == T) & (split_row == F)) {
        # print("split columns")
        ht <- (Heatmap((expression_filtered_matrix_log),
            col = ht_colors,
            column_title = column_title,
            row_title = row_title,
            column_title_gp = column_title_gp,
            row_title_gp = row_title_gp,


            #rect_gp = gpar(type = "none"),



            border = TRUE,
            show_column_names = FALSE,
            show_row_names = FALSE,
            column_split = col_split,
            top_annotation = annotation,
            heatmap_legend_param = list(
                title = "EXPRESSION LEVELS",
                title_position = "leftcenter-rot",
                at = c(0, 5, 10, 15),
                grid_width = unit(0.75, "cm"),
                legend_height = unit(6, "cm"),
                border = "black"
            ),

            column_gap = unit(5, "mm"),

            clustering_distance_rows = clustering_distance_rows,
            clustering_method_rows = clustering_method_rows,
            clustering_distance_columns = clustering_distance_columns,
            clustering_method_columns = clustering_method_columns,

            row_dend_reorder = FALSE,
            column_dend_reorder = TRUE
        )
        )
    }

    ## make the heatmap with the row splitting
    if ((split_col == F) & (split_row == T)) {
        # print("split rows")
        ht <- (Heatmap((expression_filtered_matrix_log),
            border = TRUE,
            col = ht_colors,
            column_title = column_title,
            row_title = row_title,
            column_title_gp = column_title_gp,
            row_title_gp = row_title_gp,

            show_column_names = FALSE,
            show_row_names = FALSE,

            row_split = row_split,
            row_gap = unit(5, "mm"),
            top_annotation = annotation,
            heatmap_legend_param = list(
                title = "EXPRESSION LEVELS",
                title_position = "leftcenter-rot",
                at = c(0, 5, 10, 15),
                grid_width = unit(0.75, "cm"),
                legend_height = unit(6, "cm"),
                border = "black"
            ),

            clustering_distance_rows = clustering_distance_rows,
            clustering_method_rows = clustering_method_rows,
            clustering_distance_columns = clustering_distance_columns,
            clustering_method_columns = clustering_method_columns,

            row_dend_reorder = FALSE,
            column_dend_reorder = TRUE
        )
        )
    }

    ## make the heatmap without the splitting
    if ((split_col == F) & (split_row == F)) {
        # print("no split rows and columns")
        ht <- (Heatmap((expression_filtered_matrix_log),
            border = TRUE,
            col = ht_colors,
            column_title = column_title,
            row_title = row_title,
            column_title_gp = column_title_gp,
            row_title_gp = row_title_gp,

            show_column_names = FALSE,
            show_row_names = FALSE,
            

            ## legends for the heatmap
            top_annotation = annotation,
            heatmap_legend_param = list(
                title = "EXPRESSION LEVELS",
                title_position = "leftcenter-rot",
                at = c(0, 5, 10, 15),
                grid_width = unit(0.75, "cm"),
                legend_height = unit(6, "cm"),
                border = "black"
            ),

            clustering_distance_rows = clustering_distance_rows,
            clustering_method_rows = clustering_method_rows,
            clustering_distance_columns = clustering_distance_columns,
            clustering_method_columns = clustering_method_columns,

            row_dend_reorder = F,
            column_dend_reorder = F
        )
        )
    }


    if (generation_heatmap) {
        print(ht)
        dev.off()
    }

    return(ht)

    print("###############################################################")
}


