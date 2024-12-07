.libPaths(c("usr/local/lib/R/site-library", .libPaths()))
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(SingleR)
library(scater)
library(scran)
library(scRNAseq)
library(SingleCellExperiment)
library(SpatialExperiment)
library(cowplot)


##############
# Parameters #
##############
work_dir <- "/home/ronan.jouanard/NeuroDev_ADD/Share_ronan.jouanard/EpiReg_suite"
old_work_dir <- "/home/ronan.jouanard/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"

input_model_organism <- "rat"
output_model_organism <- "human"

# For human
output_dir_human <- sprintf("/home/ronan.jouanard/NeuroDev_ADD/Share_ronan.jouanard/EpiReg_suite/2-Cell_type_mapping/output/%s", output_model_organism)

if (!dir.exists(output_dir_human)) {
        dir.create(output_dir_human, recursive = TRUE, mode = "0775")
    }

### To update
output_dir_human_seurat <- sprintf("/home/ronan.jouanard/NeuroDev_ADD/Share_ronan.jouanard/EpiReg_suite/1-Pathways_Analysis/output/%s", output_model_organism)

# For rat
output_dir_rat <- sprintf("/home/ronan.jouanard/NeuroDev_ADD/Share_ronan.jouanard/EpiReg_suite/2-Cell_type_mapping/output/%s", input_model_organism)

if (!dir.exists(output_dir_rat)) {
        dir.create(output_dir_rat, recursive = TRUE, mode = "0775")
    }

# Parameters
seurat_custom_cluster.v <- c("0.0","0.1","0.2","0.3","1","2","3.0","3.1","3.2","3.3","4")
time_points.v <- c("5","10","20","40")

# Color per cell type - Darmanis et al., 2015
col_cell_type <- c('neurons' = '#FF0000',
                   'oligodendrocytes' = '#0366FF',
                   'astrocytes' = '#BFDD00',
                   'hybrid' = '#03FFFF',
                   'endothelial' = '#FF04F0',
                   'fetal_quiescent' = '#AB04FF',
                   'OPC' = '#DB00FE',
                   'fetal_replicating' = '#00AAFF',
                   'microglia' = '#FFA600',
                   'unselected' = 'gray50')

## Replicates samples - CTRL
sample_ids_1st_repli_CTRL <- c("B_L1_S2", "B_L2_S6", "D_L1_S4", "D_L2_S8")
sample_ids_2nd_repli_CTRL <- c("B_L3_S10", "B_L4_S14", "D_L4_S16", "D_L3_S12")

## Replicates samples - SE
sample_ids_1st_repli_SE <- c("A_L1_S1", "A_L2_S5", "C_L1_S3", "C_L2_S7")
sample_ids_2nd_repli_SE <- c("A_L3_S9", "A_L4_S13", "C_L4_S15", "C_L3_S11")

# Seurat
## Human Gene Annotations
human_seurat_file <- sprintf("%s/seurat_obj/seurat_orth_human.RData", output_dir_human_seurat)
load(human_seurat_file)

# SingleR - Spot Annotation prediction
## Human Gene Annotations
### Darmanis et al., 2015
#### Retrieve dataset
sce_darmanis <- fetchDataset("darmanis-brain-2015", "2023-12-21")

# Remove unlabelled cells or cells without a clear label and normalise counts
sce_darmanis <- sce_darmanis[,!is.na(sce_darmanis$cell.type)] 

sce_darmanis_norm <- logNormCounts(sce_darmanis)

#To transfer cell types annotations between the dataset of reference (Darmanis et al., 2015) and our study spatial transcriptomics dataset (seurat_object_whole_clustering) with SingleR() we have to do it on SingleCellExperiment object. I will choose as reference and test datasets :

#ref : sce_darmanis (already a SingleCellExperiment object)
#test : seurat_object_whole_clustering (needs to be converted into a SingleCellExperiment object)

# seurat_object_whole_clustering conversion into SingleCellExperiment on processed data
sce2annotate <- as.SingleCellExperiment(seurat_orth_human, assay = "SCT")

# Cell-type labels prediction on our spatial dataset
pred.whole_clustering <- SingleR(test=sce2annotate, ref=sce_darmanis_norm, labels=sce_darmanis_norm$cell.type, de.method="wilcox", de.n = 50)

# Add cell type labels as metadata
pred.whole_clustering_df <- pred.whole_clustering %>% as.data.frame()
pred.whole_clustering_labels <- pred.whole_clustering_df[, "labels"]

seurat_orth_human <- seurat_orth_human %>% AddMetaData(metadata = pred.whole_clustering_labels, col.name = "labels")

# Darmanis et al., 2015 - Cell Type deconvolution and visualisation
cell_type_vis_output_dir <- sprintf("%s/cell_type_visualisation", output_dir_human)
    
    if (!dir.exists(cell_type_vis_output_dir)) {
        dir.create(cell_type_vis_output_dir, recursive = TRUE, mode = "0775")
    }

# Custom clustering for the visualisation
# Hippocampus clusters (3.0, 3.1, 3.2, 3.3)

cell_type_hippocampus_vis_output_dir <- sprintf("%s/cell_type_visualisation/hippocampus", output_dir_human)
    
    if (!dir.exists(cell_type_hippocampus_vis_output_dir)) {
        dir.create(cell_type_hippocampus_vis_output_dir, recursive = TRUE, mode = "0775")
    }

hippocampus_clusters <- c("3.0", "3.1", "3.2", "3.3")

for (cluster in hippocampus_clusters) {
    
    cell_type_vis_hippocampus_cluster_dir <- sprintf("%s/%s", cell_type_hippocampus_vis_output_dir, cluster)
    if (!dir.exists(cell_type_vis_hippocampus_cluster_dir)) {
        dir.create(cell_type_vis_hippocampus_cluster_dir, recursive = TRUE, mode = "0775")
    }
}
# White matter cluster (1)
cell_type_white_matter_vis_output_dir <- sprintf("%s/cell_type_visualisation/white_matter", output_dir_human)
    
    if (!dir.exists(cell_type_white_matter_vis_output_dir)) {
        dir.create(cell_type_white_matter_vis_output_dir, recursive = TRUE, mode = "0775")
    }
# All clusters
cell_type_all_clusters_vis_output_dir <- sprintf("%s/cell_type_visualisation/all_clusters", output_dir_human)
    
    if (!dir.exists(cell_type_all_clusters_vis_output_dir)) {
        dir.create(cell_type_all_clusters_vis_output_dir, recursive = TRUE, mode = "0775")
    }

Idents(seurat_orth_human) <- "seurat_custom_clusters"

## Hippocampus - Cell Type visualisation
###3.0 - 1st replicates - CTRL
list_ct_vis <- list()
seurat_orth_human@meta.data$labels <- NULL

# Get spots corresponding to cluster 3.0
selected_cluster_3.0_spots <- CellsByIdentities(seurat_orth_human, idents = "3.0") %>% unlist()
all_spots <- colnames(seurat_orth_human)
unselected_spots <- setdiff(all_spots, selected_cluster_3.0_spots)

# Update labels in pred.whole_clustering_df
pred.whole_clustering_df$labels[!rownames(pred.whole_clustering_df) %in% selected_cluster_3.0_spots] <- "unselected"

# Add the updated labels to the Seurat object metadata
pred.whole_clustering_labels <- pred.whole_clustering_df[, "labels"]
seurat_orth_human <- seurat_orth_human %>% AddMetaData(metadata = pred.whole_clustering_labels, col.name = "labels")

# Create spatial plots for each image in sample_ids_1st_repli_CTRL
for (image in sample_ids_1st_repli_CTRL) {
    plot <- SpatialDimPlot(seurat_orth_human, group.by = "labels", cols = col_cell_type, images = image, pt.size.factor = 25, image.alpha = 0) +
    ggtitle(image) +
    theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", face = "bold", size = 15), legend.position = "bottom")
    
    plot_name <- sprintf("CTRL_%s_cluster_3.0", image)
    list_ct_vis[[plot_name]] <- plot
}

# Extract the legend from the second plot in the list
legend_combined_plots <-  get_plot_component(list_ct_vis[[2]], "guide-box-bottom", return_all = TRUE)

# Remove legends from the individual plots
for (i in 1:4) {
    list_ct_vis[[i]] <- list_ct_vis[[i]] + theme(legend.position = "none")
}

# Combine the plots into a single row
combined_plot <- plot_grid(list_ct_vis[[1]], list_ct_vis[[2]], list_ct_vis[[3]], list_ct_vis[[4]], nrow = 1)

# Combine the plots with the legend, placing the legend below
combined_plot <- plot_grid(combined_plot, legend_combined_plots, ncol = 1, rel_heights = c(1, 0.1)) +
ggtitle("Cluster 3.0 - 1st replicates - Control") +
theme(plot.title = element_text(hjust = 0.5))

# Display the final plot
pdf_path <- file.path(cell_type_hippocampus_vis_output_dir, paste0("3.0/", "1st_replicates_ctrl", ".pdf"))
pdf(pdf_path, width = 8, height = 10)
print(combined_plot)
dev.off()

###3.0 - 1st replicates - SE
list_ct_vis <- list()
seurat_orth_human@meta.data$labels <- NULL

# Get spots corresponding to cluster 3.0
selected_cluster_3.0_spots <- CellsByIdentities(seurat_orth_human, idents = "3.0") %>% unlist()
all_spots <- colnames(seurat_orth_human)
unselected_spots <- setdiff(all_spots, selected_cluster_3.0_spots)

# Update labels in pred.whole_clustering_df
pred.whole_clustering_df$labels[!rownames(pred.whole_clustering_df) %in% selected_cluster_3.0_spots] <- "unselected"

# Add the updated labels to the Seurat object metadata
pred.whole_clustering_labels <- pred.whole_clustering_df[, "labels"]
seurat_orth_human <- seurat_orth_human %>% AddMetaData(metadata = pred.whole_clustering_labels, col.name = "labels")

# Create spatial plots for each image in sample_ids_1st_repli_CTRL
for (image in sample_ids_1st_repli_SE) {
    plot <- SpatialDimPlot(seurat_orth_human, group.by = "labels", cols = col_cell_type, images = image, pt.size.factor = 25, image.alpha = 0) +
    ggtitle(image) +
    theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", face = "bold", size = 15), legend.position = "bottom")
    
    plot_name <- sprintf("SE_%s_cluster_3.0", image)
    list_ct_vis[[plot_name]] <- plot
}

# Extract the legend from the second plot in the list
legend_combined_plots <-  get_plot_component(list_ct_vis[[2]], "guide-box-bottom", return_all = TRUE)

# Remove legends from the individual plots
for (i in 1:4) {
    list_ct_vis[[i]] <- list_ct_vis[[i]] + theme(legend.position = "none")
}

# Combine the plots into a single row
combined_plot <- plot_grid(list_ct_vis[[1]], list_ct_vis[[2]], list_ct_vis[[3]], list_ct_vis[[4]], nrow = 1)

# Combine the plots with the legend, placing the legend below
combined_plot <- plot_grid(combined_plot, legend_combined_plots, ncol = 1, rel_heights = c(1, 0.1)) +
ggtitle("Cluster 3.0 - 1st replicates - Status Epilepticus") +
theme(plot.title = element_text(hjust = 0.5))

# Display the final plot
pdf_path <- file.path(cell_type_hippocampus_vis_output_dir, paste0("3.0/", "1st_replicates_se", ".pdf"))
pdf(pdf_path, width = 8, height = 10)
print(combined_plot)
dev.off()

## White matter - Cell Type visualisation
###1st replicates - CTRL
list_ct_vis <- list()
seurat_orth_human@meta.data$labels <- NULL

# Get spots corresponding to cluster 1
selected_cluster_1_spots <- CellsByIdentities(seurat_orth_human, idents = "1") %>% unlist()
all_spots <- colnames(seurat_orth_human)
unselected_spots <- setdiff(all_spots, selected_cluster_1_spots)

# Update labels in pred.whole_clustering_df
pred.whole_clustering_df <- pred.whole_clustering %>% as.data.frame()
pred.whole_clustering_df$labels[!rownames(pred.whole_clustering_df) %in% selected_cluster_1_spots] <- "unselected"

# Add the updated labels to the Seurat object metadata
pred.whole_clustering_labels <- pred.whole_clustering_df[, "labels"]
seurat_orth_human <- seurat_orth_human %>% AddMetaData(metadata = pred.whole_clustering_labels, col.name = "labels")

# Create spatial plots for each image in sample_ids_1st_repli_CTRL
for (image in sample_ids_1st_repli_CTRL) {
    plot <- SpatialDimPlot(seurat_orth_human, group.by = "labels", cols = col_cell_type, images = image, pt.size.factor = 25, image.alpha = 0) +
    ggtitle(image) +
    theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", face = "bold", size = 15), legend.position = "bottom")
    
    plot_name <- sprintf("CTRL_%s_cluster_1", image)
    list_ct_vis[[plot_name]] <- plot
}

# Extract the legend from the second plot in the list
legend_combined_plots <-  get_plot_component(list_ct_vis[[2]], "guide-box-bottom", return_all = TRUE)

# Remove legends from the individual plots
for (i in 1:4) {
    list_ct_vis[[i]] <- list_ct_vis[[i]] + theme(legend.position = "none")
}

# Combine the plots into a single row
combined_plot <- plot_grid(list_ct_vis[[1]], list_ct_vis[[2]], list_ct_vis[[3]], list_ct_vis[[4]], nrow = 1)

# Combine the plots with the legend, placing the legend below
combined_plot <- plot_grid(combined_plot, legend_combined_plots, ncol = 1, rel_heights = c(1, 0.1)) +
ggtitle("Cluster 1 - 1st replicates - Control") +
theme(plot.title = element_text(hjust = 0.5))

# Display the final plot
pdf_path <- file.path(cell_type_white_matter_vis_output_dir, paste0("1st_replicates_ctrl", ".pdf"))
pdf(pdf_path, width = 8, height = 10)
print(combined_plot)
dev.off()


###1st replicates - SE
list_ct_vis <- list()
seurat_orth_human@meta.data$labels <- NULL

# Get spots corresponding to cluster 1
selected_cluster_1_spots <- CellsByIdentities(seurat_orth_human, idents = "1") %>% unlist()
all_spots <- colnames(seurat_orth_human)
unselected_spots <- setdiff(all_spots, selected_cluster_1_spots)

# Update labels in pred.whole_clustering_df
pred.whole_clustering_df <- pred.whole_clustering %>% as.data.frame()
pred.whole_clustering_df$labels[!rownames(pred.whole_clustering_df) %in% selected_cluster_1_spots] <- "unselected"

# Add the updated labels to the Seurat object metadata
pred.whole_clustering_labels <- pred.whole_clustering_df[, "labels"]
seurat_orth_human <- seurat_orth_human %>% AddMetaData(metadata = pred.whole_clustering_labels, col.name = "labels")

# Create spatial plots for each image in sample_ids_1st_repli_SE
for (image in sample_ids_1st_repli_SE) {
    plot <- SpatialDimPlot(seurat_orth_human, group.by = "labels", cols = col_cell_type, images = image, pt.size.factor = 25, image.alpha = 0) +
    ggtitle(image) +
    theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", face = "bold", size = 15), legend.position = "bottom")
    
    plot_name <- sprintf("CTRL_%s_cluster_1", image)
    list_ct_vis[[plot_name]] <- plot
}

# Extract the legend from the second plot in the list
legend_combined_plots <-  get_plot_component(list_ct_vis[[2]], "guide-box-bottom", return_all = TRUE)

# Remove legends from the individual plots
for (i in 1:4) {
    list_ct_vis[[i]] <- list_ct_vis[[i]] + theme(legend.position = "none")
}

# Combine the plots into a single row
combined_plot <- plot_grid(list_ct_vis[[1]], list_ct_vis[[2]], list_ct_vis[[3]], list_ct_vis[[4]], nrow = 1)

# Combine the plots with the legend, placing the legend below
combined_plot <- plot_grid(combined_plot, legend_combined_plots, ncol = 1, rel_heights = c(1, 0.1)) +
ggtitle("Cluster 1 - 1st replicates - Status Epilepticus") +
theme(plot.title = element_text(hjust = 0.5))

# Display the final plot
pdf_path <- file.path(cell_type_white_matter_vis_output_dir, paste0("1st_replicates_se", ".pdf"))
pdf(pdf_path, width = 8, height = 10)
print(combined_plot)
dev.off()

## All clusters - Cell Type visualisation
###1st replicates - CTRL
list_ct_vis <- list()
seurat_orth_human@meta.data$labels <- NULL

# Update labels in pred.whole_clustering_df
pred.whole_clustering_df <- pred.whole_clustering %>% as.data.frame()

# Add the updated labels to the Seurat object metadata
pred.whole_clustering_labels <- pred.whole_clustering_df[, "labels"]
seurat_orth_human <- seurat_orth_human %>% AddMetaData(metadata = pred.whole_clustering_labels, col.name = "labels")

# Create spatial plots for each image in sample_ids_1st_repli_CTRL
for (image in sample_ids_1st_repli_CTRL) {
    plot <- SpatialDimPlot(seurat_orth_human, group.by = "labels", cols = col_cell_type, images = image, pt.size.factor = 25, image.alpha = 0) +
    ggtitle(image) +
    theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", face = "bold", size = 15), legend.position = "bottom")
    
    plot_name <- sprintf("CTRL_%s_all_clusters", image)
    list_ct_vis[[plot_name]] <- plot
}

# Extract the legend from the second plot in the list
legend_combined_plots <-  get_plot_component(list_ct_vis[[2]], "guide-box-bottom", return_all = TRUE)

# Remove legends from the individual plots
for (i in 1:4) {
    list_ct_vis[[i]] <- list_ct_vis[[i]] + theme(legend.position = "none")
}

# Combine the plots into a single row
combined_plot <- plot_grid(list_ct_vis[[1]], list_ct_vis[[2]], list_ct_vis[[3]], list_ct_vis[[4]], nrow = 1)

# Combine the plots with the legend, placing the legend below
combined_plot <- plot_grid(combined_plot, legend_combined_plots, ncol = 1, rel_heights = c(1, 0.1)) +
ggtitle("All clusters - 1st replicates - Control") +
theme(plot.title = element_text(hjust = 0.5))

# Display the final plot
pdf_path <- file.path(cell_type_all_clusters_vis_output_dir, paste0("1st_replicates_ctrl", ".pdf"))
pdf(pdf_path, width = 8, height = 10)
print(combined_plot)
dev.off()

#1st replicates - SE
list_ct_vis <- list()
seurat_orth_human@meta.data$labels <- NULL

# Update labels in pred.whole_clustering_df
pred.whole_clustering_df <- pred.whole_clustering %>% as.data.frame()

# Add the updated labels to the Seurat object metadata
pred.whole_clustering_labels <- pred.whole_clustering_df[, "labels"]
seurat_orth_human <- seurat_orth_human %>% AddMetaData(metadata = pred.whole_clustering_labels, col.name = "labels")

# Create spatial plots for each image in sample_ids_1st_repli_SE
for (image in sample_ids_1st_repli_SE) {
    plot <- SpatialDimPlot(seurat_orth_human, group.by = "labels", cols = col_cell_type, images = image, pt.size.factor = 25, image.alpha = 0) +
    ggtitle(image) +
    theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", face = "bold", size = 15), legend.position = "bottom")
    
    plot_name <- sprintf("SE_%s_all_clusters", image)
    list_ct_vis[[plot_name]] <- plot
}

# Extract the legend from the second plot in the list
legend_combined_plots <-  get_plot_component(list_ct_vis[[2]], "guide-box-bottom", return_all = TRUE)

# Remove legends from the individual plots
for (i in 1:4) {
    list_ct_vis[[i]] <- list_ct_vis[[i]] + theme(legend.position = "none")
}

# Combine the plots into a single row
combined_plot <- plot_grid(list_ct_vis[[1]], list_ct_vis[[2]], list_ct_vis[[3]], list_ct_vis[[4]], nrow = 1)

# Combine the plots with the legend, placing the legend below
combined_plot <- plot_grid(combined_plot, legend_combined_plots, ncol = 1, rel_heights = c(1, 0.1)) +
ggtitle("All clusters - 1st replicates - Status Epilepticus") +
theme(plot.title = element_text(hjust = 0.5))

# Display the final plot
pdf_path <- file.path(cell_type_all_clusters_vis_output_dir, paste0("1st_replicates_se", ".pdf"))
pdf(pdf_path, width = 8, height = 10)
print(combined_plot)
dev.off()

# Stacked bar blot - Cell type proportion
