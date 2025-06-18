.libPaths(c("usr/local/lib/R/site-library", .libPaths()))
library(magrittr)
library(SCPA)
library(Seurat)
library(msigdbr)
library(withr)
library(ggplot2)
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
filtering_methods <- c("top_n", "q_val")
thresholds_q_val <- 1:4
thresholds_top_n <- seq(10, 50, by = 10)

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
    assign(SCPA_output, read.table(file = sprintf("%s/%s.csv", SCPA_output_dir_human, collection)))
    }

## Filtering and ranking SCPA_output

for (collection in pathway_collections) {
    for (sorting_method in sorting_methods) {
        for (filtering_method in filtering_methods) {
            if (filtering_method == "q-val") {
                for (threshold_q_val in thresholds_q_val) {
                    result_name <- paste("SCPA", collection, sorting_method, filtering_method, threshold_q_val, sep = "_")
                    result <- main_pathways(get(paste0("SCPA_", collection)), collection = collection, sorting_method = sorting_method, filtering_method = filtering_method, threshold = threshold_q_val)
                    list_filtered_scpa[[result_name]] <- result
                }
            } else {
                for (threshold_top_n in thresholds_top_n) {
                    result_name <- paste("SCPA", collection, sorting_method, filtering_method, threshold_top_n, sep = "_")
                    result <- main_pathways(get(paste0("SCPA_", collection)), collection = collection, sorting_method = sorting_method, filtering_method = filtering_method, top_n_pathways = threshold_top_n)
                    list_filtered_scpa[[result_name]] <- result
                }
            }
        }
    }
}

