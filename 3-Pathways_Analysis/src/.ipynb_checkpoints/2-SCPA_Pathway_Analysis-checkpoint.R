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
output_model_organisms <- c("human", "mouse")

# For human
output_dir_human <- sprintf("/home/ronan.jouanard/NeuroDev_ADD/Share_ronan.jouanard/EpiReg_suite/1-Pathways_Analysis/output/%s", output_model_organisms[1])

# For mouse
output_dir_mouse <- sprintf("/home/ronan.jouanard/NeuroDev_ADD/Share_ronan.jouanard/EpiReg_suite/1-Pathways_Analysis/output/%s", output_model_organisms[2])

# Parameters
seurat_custom_cluster.v <- c("0.0","0.1","0.2","0.3","1","2","3.0","3.1","3.2","3.3","4")
time_points.v <- c("5","10","20","40") 

##############################
# Functions for the analysis #
##############################

get_qvals <- function(scpa_out, name) {
  
  df <- list()
  for (i in names(scpa_out)) {
    df[[i]] <- scpa_out[[i]] %>%
      select(Pathway, qval)
  }
  
  col_names <- names(df)
  for (i in 1:length(df)) {
    df[[i]] <- set_colnames(df[[i]], c("pathway", paste(name, col_names[[i]], sep = "_")))
  }
  
  return(df)
  
}

scpa_subclusters <- function(split_seurat, subclusters, time_points, organism, pathway_collection, specific_pathways = NULL, downsample = 500,
                             min_genes = 15, max_genes = 500, parallel = FALSE, cores = NULL){
  
  # Initialize a SCPA outputs list for a given cluster 
  scpa_pways <- list()

  # Extracting pathways from a given collection
  pathways_list <- msigdbr(organism) %>%
    filter(grepl(pathway_collection, x = gs_name, ignore.case = TRUE)) %>%
    format_pathways()

  # Filter pathways if specific_pathways are provided
  if (!is.null(specific_pathways)) {
    pathways_list <- lapply(pathways_list, function(tbl) {
      tbl %>%
        filter(any(sapply(specific_pathways, function(y) grepl(y, Pathway, ignore.case = TRUE))))
    })
    # Remove empty tibbles from the list
    pathways_list <- pathways_list[sapply(pathways_list, nrow) > 0]
  }

  # Loop over each subcluster
  for (cluster_id in subclusters) {
    index <- which(cluster_id == subclusters)
    se_vs_ctrl <- list()
    
    # Loop over each time point
    for (time in time_points) {
      # Extract data for SE condition
      se <- seurat_extract(split_seurat$SE, 
                           meta1 = "time", value_meta1 = time,
                           meta2 = "seurat_custom_clusters", value_meta2 = cluster_id)
  
      # Extract data for CTRL condition
      ctrl <- seurat_extract(split_seurat$CTRL, 
                             meta1 = "time", value_meta1 = time,
                             meta2 = "seurat_custom_clusters", value_meta2 = cluster_id)
      
      # Compare pathways
      cat("Comparing cluster", cluster_id, "at time", time, "\n")
      se_vs_ctrl[[time]] <- suppressWarnings(compare_pathways(samples = list(se, ctrl), pathways = pathways_list, downsample = downsample,
                                             min_genes = min_genes, max_genes = max_genes, parallel = parallel, cores = cores))
    }
    
    # Combine results for the current subcluster
    scpa_pways[[index]] <- Reduce(full_join, c(get_qvals(se_vs_ctrl, "SE"))) %>%
      column_to_rownames("pathway") %>%
      set_colnames(paste(cluster_id, colnames(.), sep = "_")) %>%
      rownames_to_column("pathway")
  }
  
  # Combine all subcluster results into a single dataframe
  all_data <- Reduce(function(x, y) full_join(x, y, by = "pathway"), scpa_pways) %>%
    column_to_rownames("pathway")

  # Remove pathways collection name from rownames in all_data
  rownames(all_data) <- gsub(paste(pathway_collection, "_"), "", rownames(all_data))
                     
  return(all_data)
}

# load data
## Seurat objects
### Human Genes Annotations

output_dir_seurat_human <- sprintf("%s/seurat_obj", output_dir_human)
seurat_orth_human_file <- file.path(output_dir_seurat_human, "seurat_orth_human.RData")

load(seurat_orth_human_file)

split_condition_seurat_orth_human <- SplitObject(seurat_orth_human, split.by = "condition")

### Mouse Genes Annotations
output_dir_seurat_mouse <- sprintf("%s/seurat_obj", output_dir_mouse)
seurat_orth_mouse_file <- file.path(output_dir_seurat_mouse, "seurat_orth_mouse.RData")

load(seurat_orth_mouse_file)

split_condition_seurat_orth_mouse <- SplitObject(seurat_orth_mouse, split.by = "condition")

# Pathways Analysis
## Human
### Hallmark
SCPA_human_HALLMARK <- scpa_subclusters(split_seurat = split_condition_seurat_orth_human, subclusters = seurat_custom_cluster.v , time_points = time_points.v, organism = "Homo sapiens", pathway_collection = "HALLMARK", parallel = TRUE)

### Kegg
SCPA_human_KEGG <- scpa_subclusters(split_seurat = split_condition_seurat_orth_human, subclusters = seurat_custom_cluster.v , time_points = time_points.v, organism = "Homo sapiens", pathway_collection = "KEGG", parallel = TRUE)

### Reactome
SCPA_human_REACTOME <- scpa_subclusters(split_seurat = split_condition_seurat_orth_human, subclusters = seurat_custom_cluster.v , time_points = time_points.v, organism = "Homo sapiens", pathway_collection = "REACTOME", parallel = TRUE)

### GOBP
SCPA_human_GOBP <- scpa_subclusters(split_seurat = split_condition_seurat_orth_human, subclusters = seurat_custom_cluster.v , time_points = time_points.v, organism = "Homo sapiens", pathway_collection = "GOBP", parallel = TRUE)


# Save seurat objects
## Human
output_dir_SCPA_human <- sprintf("%s/SCPA_analysis", output_dir_human)
if (! dir.exists(output_dir_SCPA_human)) {
    dir.create(output_dir_SCPA_human, recursive=TRUE, mode="0775")
}

### HALLMARK
                     
write.csv(SCPA_human_HALLMARK, file = sprintf("%s/HALLMARK.csv", output_dir_SCPA_human), row.names = TRUE)

### KEGG

write.csv(SCPA_human_KEGG, file = sprintf("%s/KEGG.csv", output_dir_SCPA_human), row.names = TRUE)

### REACTOME
write.csv(x = SCPA_human_REACTOME, file = sprintf("%s/REACTOME.csv", output_dir_SCPA_human), row.names = TRUE)
