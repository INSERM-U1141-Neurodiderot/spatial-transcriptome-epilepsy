.libPaths(c("usr/local/lib/R/site-library", .libPaths()))
library(magrittr)
library(SCPA)
library(Seurat)
library(SingleCellExperiment)
library(SpatialExperiment)
library(ggspavis)
library(msigdbr)
library(withr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(backports)
library(ggpubr)
library(memoise, attach.required=T)
library(ggtree, attach.required=T)
library(gprofiler2)
library(grr)
library(orthogene)
library(RCurl)
library(ComplexHeatmap)
library(tzdb)
library(readr)
library(forcats)
library(lubridate)
library(tidyverse)
library(multicross)
library(circlize)

##############
# Parameters #
##############
work_dir <- "/home/ronan.jouanard/NeuroDev_ADD/Share_ronan.jouanard/EpiReg_suite"
old_work_dir <- "/home/ronan.jouanard/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"

input_model_organism <- "rat"
output_model_organism <- "human"

# For human
output_dir_human <- sprintf("/home/ronan.jouanard/NeuroDev_ADD/Share_ronan.jouanard/EpiReg_suite/1-Pathways_Analysis/output/%s", output_model_organism)

# For rat
output_dir_rat <- sprintf("/home/ronan.jouanard/NeuroDev_ADD/Share_ronan.jouanard/EpiReg_suite/1-Pathways_Analysis/output/%s", input_model_organism)

# Parameters
seurat_custom_cluster.v <- c("0.0","0.1","0.2","0.3","1","2","3.0","3.1","3.2","3.3","4")
time_points.v <- c("5","10","20","40")
pathway_collections <- c("HALLMARK", "REACTOME", "GOBP", "KEGG")

# Filtering and ranking parameters
sorting_methods <- c("mean", "max", "median")
filtering_methods <- c("q-val", "top-n")
thresholds_q_val <- 1:4 %>% as.character()
thresholds_top_n <- seq(10, 50, by = 10) %>% as.character()

# Heatmap design parameters
breaks_gradient <- c(0, 5, 10)
colors_gradient <- c("white", "#85c6ad", "#304d36")
color_gradient_mapping <- colorRamp2(breaks_gradient, colors_gradient)

heatmap_structure <- c("No_ordered_by_time_no_clusters_hierachisation", "ordered_by_time_no_clusters_hierachisation", "ordered_by_time_clusters_hierachisation")

##############################
# Functions for the analysis #
##############################

main_pathways <- function(SCPA_output, collection, colnames = NULL, 
                          sorting_method, 
                          filtering_method,
                          top_n_pathways = 50, 
                          threshold = 2){
  
  # Columns extraction
    if (length(colnames) == 0) {
    data_extracted <- SCPA_output
        } else {data_extracted <- SCPA_output[,colnames]}
    
    # Remove name of collection pathways in rownames
    rownames(data_extracted) <- gsub(paste(collection, "_", sep=""),"", rownames(data_extracted))
    
    # Compute mean and median values of q-value for each pathway in data_extracted
    mean_q_val_data_extracted <- apply(data_extracted, 1, function(x){mean(x)})
    median_q_val_data_extracted <- apply(data_extracted, 1, function(x){median(x)})
    max_q_val_data_extracted <- apply(data_extracted, 1, function(x){max(x)})

    # Store them in two last columns of stats_data_extracted
    stats_data_extracted <- data.frame(data_extracted,
                                       mean_q_val_data_extracted,
                                       median_q_val_data_extracted,
                                       max_q_val_data_extracted)


 if ((sorting_method == "mean") & (filtering_method == "q-val")) {
     
     stats_data_extracted_sorted <- stats_data_extracted[order(mean_q_val_data_extracted,
                                                               decreasing = TRUE), ]
     selection_data_extracted <- stats_data_extracted_sorted %>% filter(mean_q_val_data_extracted > threshold)

 } else if ((sorting_method == "median") & (filtering_method == "q-val")) {
     
     stats_data_extracted_sorted <- stats_data_extracted[order(median_q_val_data_extracted,
                                                               decreasing = TRUE), ]
     selection_data_extracted <- stats_data_extracted_sorted %>% filter(median_q_val_data_extracted > threshold)
     
 } else if ((sorting_method == "max") & (filtering_method == "q-val")) {
     
     stats_data_extracted_sorted <- stats_data_extracted[order(max_q_val_data_extracted,
                                                               decreasing = TRUE), ]
     
     selection_data_extracted <- stats_data_extracted_sorted %>% filter(max_q_val_data_extracted > threshold)
     
 } else if ((sorting_method == "mean") & (filtering_method == "top-n")) {
     
     stats_data_extracted_sorted <- stats_data_extracted[order(mean_q_val_data_extracted,
                                                               decreasing = TRUE), ]
     selection_data_extracted <- stats_data_extracted_sorted[1:top_n_pathways,]

 } else if ((sorting_method == "median") & (filtering_method == "top-n")) {
     
     stats_data_extracted_sorted <- stats_data_extracted[order(median_q_val_data_extracted,
                                                               decreasing = TRUE), ]
     selection_data_extracted <- stats_data_extracted_sorted[1:top_n_pathways,]

 } else {
     
     stats_data_extracted_sorted <- stats_data_extracted[order(max_q_val_data_extracted,
                                                               decreasing = TRUE), ]
     selection_data_extracted <- stats_data_extracted_sorted[1:top_n_pathways,]            
}

selection_data_extracted_final <- selection_data_extracted[, 1:44]
return(selection_data_extracted_final) 
}

Heatmap_structure <- function(SCPA_output, cluster_hierarchisation = FALSE, time_ordered = FALSE) {
    # Removing X prefix in colnames
    colnames(SCPA_output) <- gsub("X", "", colnames(SCPA_output))
    
    # Clusters and time levels extraction
    
    cluster <- colnames(SCPA_output) %>% strsplit("_") %>% sapply(function(x) x[1]) %>% as.numeric() %>% as.factor()
    names(cluster) <- colnames(SCPA_output)                                                             
    
    time <- colnames(SCPA_output) %>% strsplit("_") %>% sapply(function(x) x[length(x)]) %>% as.numeric() %>% as.factor() 
    names(time) <- colnames(SCPA_output)
    
    # Conditional structuration of heatmap columns
    if ((cluster_hierarchisation == FALSE) & (time_ordered == FALSE)) {
        
        clusters_order <- cluster
        times_order <- time
        SCPA_output_order <- SCPA_output
        
    } else if ((cluster_hierarchisation == TRUE) & (time_ordered == FALSE)) {
        
        clusters_order <- cluster[order(cluster)]
        times_order <- time[order(cluster)]
        SCPA_output_order <- SCPA_output[,order(cluster)]
        
    } else if ((cluster_hierarchisation == TRUE) & (time_ordered == TRUE)) {
        
        clusters_order <- cluster[order(cluster, time)]
        times_order <- time[order(cluster, time)]
        SCPA_output_order <- SCPA_output[,order(cluster, time)]
        
    } else {
        
        clusters_order <- cluster[order(time)]
        times_order <- time[order(time)]
        SCPA_output_order <- SCPA_output[,order(time)]
    }
        # Get back to character vectors
        clusters_order <- as.character(clusters_order)
        times_order <- as.character(times_order)
                                                                 
        # Get back to initial cluster IDs
        clusters_order[clusters_order == "0"] <- "0.0"
        clusters_order[clusters_order == "3"] <- "3.0"
    
                                                               
  return(list(clusters_order, times_order, SCPA_output_order))
}


# Heatmap - Most representative genesets

SCPA_human_output_dir <- sprintf("%s/SCPA_analysis", output_dir_human)

for (collection in pathway_collections) {
    SCPA_output <- paste("SCPA", collection, sep = "_")
    assign(SCPA_output, read.csv(file = sprintf("%s/%s.csv", SCPA_human_output_dir, collection)))
    
    # Convert the first column to row names
    df <- get(SCPA_output) %>% column_to_rownames(var = "X")
    
    # Reassign the modified dataframe back to the variable
    assign(SCPA_output, df)
}

## Filtering and ranking SCPA_output


list_filtered_scpa <- list()

for (collection in pathway_collections) {
    for (sorting_method in sorting_methods) {
        for (filtering_method in filtering_methods) {
            if (filtering_method == "q-val") {
                for (threshold_q_val in thresholds_q_val) {
                    result_qval_name <- paste("SCPA", collection, sorting_method, filtering_method, threshold_q_val, sep = "_")
                    result_qval <- main_pathways(
                        get(paste0("SCPA_", collection)), 
                        collection = collection, 
                        sorting_method = sorting_method, 
                        filtering_method = filtering_method, 
                        threshold = threshold_q_val
                    )
                    list_filtered_scpa[[result_qval_name]] <- result_qval
                }
            } else {
                for (threshold_top_n in thresholds_top_n) {
                    result_topn_name <- paste("SCPA", collection, sorting_method, filtering_method, threshold_top_n, sep = "_")
                    result_topn <- main_pathways(
                        get(paste0("SCPA_", collection)), 
                        collection = collection, 
                        sorting_method = sorting_method, 
                        filtering_method = filtering_method, 
                        top_n_pathways = threshold_top_n
                    )
                    list_filtered_scpa[[result_topn_name]] <- result_topn
                }
            }
        }
    }
}

## Removing empty dataframes from SCPA filtered list

empty_dataframes <- sapply(list_filtered_scpa, function(x) nrow(x) == 0)
names_of_empty_dataframes <- names(empty_dataframes[empty_dataframes])
                           
if (length(names_of_empty_dataframes) > 0) {
  cat("Filtering of following SCPA outputs leads to the removal of all pathways:\n")
  print(names_of_empty_dataframes)
}

list_filtered_scpa <- list_filtered_scpa[!empty_dataframes]

## Heatmaps generation

ht_opt$message <- FALSE
                           
time_legend <- Legend(at = c("5", "10", "20", "40"), 
                      legend_gp = gpar(fill = c("#feedde", "#fdbe85", "#fd8d3c", "#d94701")),
                      title = "Time")

SCPA_heatmap_human_output_dir <- sprintf("%s/visualisation/heatmaps", SCPA_human_output_dir)

for (structure in heatmap_structure) {
    structure_output_dir <- sprintf("%s/%s", SCPA_heatmap_human_output_dir, structure)
    
    if (!dir.exists(structure_output_dir)) {
        dir.create(structure_output_dir, recursive = TRUE, mode = "0775")
    }
    
    if (structure == heatmap_structure[1]) {
        structure_param <- c(FALSE, FALSE)
    } else if (structure == heatmap_structure[2]) {
        structure_param <- c(TRUE, FALSE)
    } else {
        structure_param <- c(TRUE, TRUE)
    }
    
    for (sorting_method in sorting_methods) {
        sorting_method_output_dir <- sprintf("%s/%s", structure_output_dir, sorting_method)

        if (!dir.exists(sorting_method_output_dir)) {
            dir.create(sorting_method_output_dir, recursive = TRUE, mode = "0775")
        }
        
        for (filtering_method in filtering_methods) {
            filtering_method_output_dir <- sprintf("%s/%s", sorting_method_output_dir, filtering_method)

            if (!dir.exists(filtering_method_output_dir)) {
                dir.create(filtering_method_output_dir, recursive = TRUE, mode = "0775")
            }
            
            if (filtering_method == "q-val") {
                for (threshold_q_val in thresholds_q_val) {
                    threshold_q_val_output_dir <- sprintf("%s/%s", filtering_method_output_dir, threshold_q_val)

                    if (!dir.exists(threshold_q_val_output_dir)) {
                        dir.create(threshold_q_val_output_dir, recursive = TRUE, mode = "0775")
                    }
                    
                    for (collection in pathway_collections) {
                        SCPA_name_qval <- sprintf("SCPA_%s_%s_%s_%s", collection, sorting_method, filtering_method, threshold_q_val)
                        SCPA_output_qval <- list_filtered_scpa[[SCPA_name_qval]]
                        
                        if (!is.null(SCPA_output_qval)) {
                            pdf_path <- file.path(threshold_q_val_output_dir, paste0("heatmap_", SCPA_name_qval, ".pdf"))

                            pdf(pdf_path, width = 8, height = 10)
                            
                            Heatmap_structure_list <- Heatmap_structure(SCPA_output_qval, time_ordered = structure_param[1], cluster_hierarchisation = structure_param[2])
                            clusters_order <- Heatmap_structure_list[[1]]
                            times_order <- Heatmap_structure_list[[2]]
                            SCPA_output_order <- Heatmap_structure_list[[3]]
                            
                            column_annotation <- HeatmapAnnotation(
                                Cluster = clusters_order,
                                Time = times_order,
                                col = list(Cluster = c("1" = "#cccccc", "2" = "#f768a1", "4" = "#bae1ff", "0.0" = "#238b45", "0.1" = "#74c476", "0.2" = "#bae4b3", "0.3" = "#edf9e9", "3.0" = "#fdf498", "3.1" = "#807dba", "3.2" = "#4a1486", "3.3" = "#ac8fff"),
                                           Time = c("5" = "#feedde", "10" = "#fdbe85", "20" = "#fd8d3c", "40" = "#d94701")),
                                gp = gpar(col = "white", lwd = 0.05),
                                annotation_name_gp = gpar(fontsize = 9),
                                simple_anno_size = unit(3, "mm"),
                                show_legend = c(TRUE, FALSE)
                            )

                            show_row_names.bool <- isTRUE(nrow(SCPA_output_order) <= 1000)
                            
                            ht <- Heatmap(data.matrix(SCPA_output_order), 
                                          name = "Qval",
                                          col = color_gradient_mapping,
                                          column_title = SCPA_name_qval,
                                          show_row_names = show_row_names.bool, 
                                          top_annotation = column_annotation,
                                          border = TRUE,
                                          show_row_dend = FALSE,
                                          show_column_dend = TRUE,
                                          show_column_names = TRUE,
                                          row_names_gp = gpar(fontsize = 6),
                                          row_dend_side = "right",
                                          cluster_columns = FALSE)

                            draw(ht, annotation_legend_list = list(time_legend))
                            dev.off()
                        }
                    }
                }
            }
            
            if (filtering_method == "top-n") {
                for (threshold_top_n in thresholds_top_n) {
                    threshold_top_n_output_dir <- sprintf("%s/%s", filtering_method_output_dir, threshold_top_n)

                    if (!dir.exists(threshold_top_n_output_dir)) {
                        dir.create(threshold_top_n_output_dir, recursive = TRUE, mode = "0775")
                    }
                    
                    for (collection in pathway_collections) {
                        SCPA_name_topn <- sprintf("SCPA_%s_%s_%s_%s", collection, sorting_method, filtering_method, threshold_top_n)
                        SCPA_output_topn <- list_filtered_scpa[[SCPA_name_topn]]
                        
                        if (!is.null(SCPA_output_topn)) {
                            pdf_path <- file.path(threshold_top_n_output_dir, paste0("heatmap_", SCPA_name_topn, ".pdf"))

                            pdf(pdf_path, width = 8, height = 10)
                            
                            Heatmap_structure_list <- Heatmap_structure(SCPA_output_topn, time_ordered = structure_param[1], cluster_hierarchisation = structure_param[2])
                            clusters_order <- Heatmap_structure_list[[1]]
                            times_order <- Heatmap_structure_list[[2]]
                            SCPA_output_order <- Heatmap_structure_list[[3]]
                            
                            column_annotation <- HeatmapAnnotation(
                                Cluster = clusters_order,
                                Time = times_order,
                                col = list(Cluster = c("1" = "#cccccc", "2" = "#f768a1", "4" = "#bae1ff", "0.0" = "#238b45", "0.1" = "#74c476", "0.2" = "#bae4b3", "0.3" = "#edf9e9", "3.0" = "#fdf498", "3.1" = "#807dba", "3.2" = "#4a1486", "3.3" = "#ac8fff"),
                                           Time = c("5" = "#feedde", "10" = "#fdbe85", "20" = "#fd8d3c", "40" = "#d94701")),
                                gp = gpar(col = "white", lwd = 0.05),
                                annotation_name_gp = gpar(fontsize = 9),
                                simple_anno_size = unit(3, "mm"),
                                show_legend = c(TRUE, FALSE)
                            )

                            show_row_names.bool <- isTRUE(nrow(SCPA_output_order) <= 1000)
                            
                            ht <- Heatmap(data.matrix(SCPA_output_order), 
                                          name = "Qval",
                                          col = color_gradient_mapping,
                                          column_title = SCPA_name_topn,
                                          show_row_names = show_row_names.bool, 
                                          top_annotation = column_annotation,
                                          border = TRUE,
                                          show_row_dend = FALSE,
                                          show_column_dend = TRUE,
                                          show_column_names = TRUE,
                                          row_names_gp = gpar(fontsize = 6),
                                          row_dend_side = "right",
                                          cluster_columns = FALSE)

                            draw(ht, annotation_legend_list = list(time_legend))
                            dev.off()
                        }
                    }
                }
            }
        }
    }
}