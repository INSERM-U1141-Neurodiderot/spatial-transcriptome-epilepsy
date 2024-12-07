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
hippocampus_clusters <- c("3.0", "3.1", "3.2", "3.3")
epil_conditions <- c("CTRL", "SE")

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

# Convert Seurat object and extract metadata
metadata_sce2annotate <- as.SingleCellExperiment(seurat_orth_human, assay = "RNA") %>%
                         colData() %>% as.data.frame()

### - Hippocampus - iterate over spatial clusters, time points, and epil conditions
cell.type_fraction <- list()
for (spatial_cluster in hippocampus_clusters) {
    for (time_point in time_points.v) {
        for (epil_condition in epil_conditions) {
            
            # Filter barcodes based on condition, cluster, and time
            barcodes_cluster_condition_time <- metadata_sce2annotate %>%
                filter(condition == !!epil_condition, seurat_custom_clusters == !!spatial_cluster, time == !!time_point) %>%
                rownames()
                                
            # Get cell type barcodes
            cell.type_barcodes <- pred.whole_clustering_df %>%
                filter(rownames(pred.whole_clustering_df) %in% barcodes_cluster_condition_time) %>%
                rownames_to_column("row_names") %>%
                select(row_names, cell_type = labels) %>%
                mutate(cell_type = factor(cell_type))
            
            # Calculate proportions
            cluster_condition_time_cell.type.prop <- prop.table(table(cell.type_barcodes$cell_type))
            
            fract_category <- sprintf("%s_%s_%s", spatial_cluster, epil_condition, time_point)
            
            cell.type_fraction[[fract_category]] <- as.data.frame(cluster_condition_time_cell.type.prop)
        }
    }
}

# Combine the list of data frames into one data frame
combined_df <- bind_rows(cell.type_fraction, .id = "cluster_condition_time")

# Change column names
colnames(combined_df) <- c("category", "cell_type", "fraction")

# Reorder categories for barplots
combined_df$category <- factor(combined_df$category, levels=names(cell.type_fraction))

levels_time_condition <- levels(combined_df$category) %>% sub("^[^_]*_", "", .) %>% unique()

combined_df$cluster <- sub("_.*", "", combined_df$category)

combined_df$time_condition <- factor(sub("^[^_]*_", "", combined_df$category), levels = levels_time_condition)

# Collect all unique cell types across clusters for the comprehensive legend
all_cell_types <- unique(combined_df$cell_type)

# Separate the data by clusters and create plots
plots <- list()
unique_clusters <- unique(combined_df$cluster)
for (i in seq_along(unique_clusters)) {
    i_cluster <- unique_clusters[i]
    plot_data <- combined_df %>% filter(cluster == !!i_cluster)
    
    p <- ggplot(plot_data, aes(x = time_condition, y = fraction, fill = cell_type)) +
        geom_bar(stat = "identity", position = "stack") +
        labs(x = "", y = i_cluster) +
        scale_fill_manual(values = col_cell_type) +
        theme_minimal() +
        theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm"))

    # Apply x axis theme conditionally
    if (i < length(unique_clusters)) {
        p <- p + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank())
    }

    plots[[i_cluster]] <- p
}

# Combine the plots vertically
pcol <- plot_grid(plotlist = plots, align = "v", ncol = 1)


# Create a dummy plot to generate a comprehensive legend
dummy_data <- data.frame(cell_type = factor(all_cell_types, levels = all_cell_types), fraction = 1)
dummy_plot <- ggplot(dummy_data, aes(x = cell_type, y = fraction, fill = cell_type)) +
              geom_bar(stat = "identity") +
              scale_fill_manual(values = col_cell_type) +
              theme_void() +
              theme(legend.position = "bottom")
legend <- get_plot_component(dummy_plot, "guide-box-bottom", return_all = TRUE)
hippocampus_combined_plot <- plot_grid(pcol, legend, ncol = 1, rel_heights = c(1, .1))

pdf_path <- file.path(output_dir_human, paste0("barplot_cell_type_proportion/st_barplot_cell_type_proportion_darmanis_annot", ".pdf"))
pdf(pdf_path, width = 8, height = 10)
print(hippocampus_combined_plot)
dev.off()

### - Hippocampus - iterate over spatial clusters, time points, and epil conditions
cell.type_fraction <- list()
for (spatial_cluster in hippocampus_clusters) {
    for (time_point in time_points.v) {
        for (epil_condition in epil_conditions) {
            
            # Filter barcodes based on condition, cluster, and time
            barcodes_cluster_condition_time <- metadata_sce2annotate %>%
                filter(condition == !!epil_condition, seurat_custom_clusters == !!spatial_cluster, time == !!time_point) %>%
                rownames()
                                
            # Get cell type barcodes
            cell.type_barcodes <- pred.whole_clustering_df %>%
                filter(rownames(pred.whole_clustering_df) %in% barcodes_cluster_condition_time) %>%
                rownames_to_column("row_names") %>%
                select(row_names, cell_type = labels) %>%
                mutate(cell_type = factor(cell_type))
            
            # Calculate proportions
            cluster_condition_time_cell.type.prop <- prop.table(table(cell.type_barcodes$cell_type))
            
            fract_category <- sprintf("%s_%s_%s", spatial_cluster, epil_condition, time_point)
            
            cell.type_fraction[[fract_category]] <- as.data.frame(cluster_condition_time_cell.type.prop)
        }
    }
}

# Combine the list of data frames into one data frame
combined_df <- bind_rows(cell.type_fraction, .id = "cluster_condition_time")

# Change column names
colnames(combined_df) <- c("category", "cell_type", "fraction")

# Reorder categories for barplots
combined_df$category <- factor(combined_df$category, levels=names(cell.type_fraction))

levels_time_condition <- levels(combined_df$category) %>% sub("^[^_]*_", "", .) %>% unique()

combined_df$cluster <- sub("_.*", "", combined_df$category)

combined_df$time_condition <- factor(sub("^[^_]*_", "", combined_df$category), levels = levels_time_condition)

# Collect all unique cell types across clusters for the comprehensive legend
all_cell_types <- unique(combined_df$cell_type)

# Separate the data by clusters and create plots
plots <- list()
unique_clusters <- unique(combined_df$cluster)
for (i in seq_along(unique_clusters)) {
    i_cluster <- unique_clusters[i]
    plot_data <- combined_df %>% filter(cluster == !!i_cluster)
    
    p <- ggplot(plot_data, aes(x = time_condition, y = fraction, fill = cell_type)) +
        geom_bar(stat = "identity", position = "stack") +
        labs(x = "", y = i_cluster) +
        scale_fill_manual(values = col_cell_type) +
        theme_minimal() +
        theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm"))

    # Apply x axis theme conditionally
    if (i < length(unique_clusters)) {
        p <- p + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank())
    }

    plots[[i_cluster]] <- p
}

# Combine the plots vertically
pcol <- plot_grid(plotlist = plots, align = "v", ncol = 1)


# Create a dummy plot to generate a comprehensive legend
dummy_data <- data.frame(cell_type = factor(all_cell_types, levels = all_cell_types), fraction = 1)
dummy_plot <- ggplot(dummy_data, aes(x = cell_type, y = fraction, fill = cell_type)) +
              geom_bar(stat = "identity") +
              scale_fill_manual(values = col_cell_type) +
              theme_void() +
              theme(legend.position = "bottom")
legend <- get_plot_component(dummy_plot, "guide-box-bottom", return_all = TRUE)
hippocampus_combined_plot <- plot_grid(pcol, legend, ncol = 1, rel_heights = c(1, .1))

pdf_path <- file.path(output_dir_human, paste0("barplot_cell_type_proportion/st_barplot_hippocampus_cell_type_proportion_darmanis_annot", ".pdf"))
pdf(pdf_path, width = 8, height = 10)
print(hippocampus_combined_plot)
dev.off()

### - White matter - iterate over spatial clusters, time points, and epil conditions
cell.type_fraction <- list()
for (time_point in time_points.v) {
    for (epil_condition in epil_conditions) {
        
        # Filter barcodes based on condition, cluster, and time
        barcodes_cluster_condition_time <- metadata_sce2annotate %>%
            filter(condition == !!epil_condition, seurat_custom_clusters == "1", time == !!time_point) %>%
            rownames()
                            
        # Get cell type barcodes
        cell.type_barcodes <- pred.whole_clustering_df %>%
            filter(rownames(pred.whole_clustering_df) %in% barcodes_cluster_condition_time) %>%
            rownames_to_column("row_names") %>%
            select(row_names, cell_type = labels) %>%
            mutate(cell_type = factor(cell_type))
        
        # Calculate proportions
        cluster_condition_time_cell.type.prop <- prop.table(table(cell.type_barcodes$cell_type))
        
        fract_category <- sprintf("%s_%s_%s", spatial_cluster, epil_condition, time_point)
        
        cell.type_fraction[[fract_category]] <- as.data.frame(cluster_condition_time_cell.type.prop)
    }
}


# Combine the list of data frames into one data frame
combined_df <- bind_rows(cell.type_fraction, .id = "cluster_condition_time")

# Change column names
colnames(combined_df) <- c("category", "cell_type", "fraction")

# Reorder categories for barplots
combined_df$category <- factor(combined_df$category, levels=names(cell.type_fraction))

levels_time_condition <- levels(combined_df$category) %>% sub("^[^_]*_", "", .) %>% unique()

combined_df$cluster <- sub("_.*", "", combined_df$category)

combined_df$time_condition <- factor(sub("^[^_]*_", "", combined_df$category), levels = levels_time_condition)

# Collect all unique cell types across clusters for the comprehensive legend
all_cell_types <- unique(combined_df$cell_type)

# Separate the data by clusters and create plots
plots <- list()
unique_clusters <- unique(combined_df$cluster)
for (i in seq_along(unique_clusters)) {
    i_cluster <- unique_clusters[i]
    plot_data <- combined_df %>% filter(cluster == !!i_cluster)
    
    p <- ggplot(plot_data, aes(x = time_condition, y = fraction, fill = cell_type)) +
        geom_bar(stat = "identity", position = "stack") +
        labs(x = "", y = i_cluster) +
        scale_fill_manual(values = col_cell_type) +
        theme_minimal() +
        theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm"))

    # Apply x axis theme conditionally
    if (i < length(unique_clusters)) {
        p <- p + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank())
    }

    plots[[i_cluster]] <- p
}

# Combine the plots vertically
pcol <- plot_grid(plotlist = plots, align = "v", ncol = 1)


# Create a dummy plot to generate a comprehensive legend
dummy_data <- data.frame(cell_type = factor(all_cell_types, levels = all_cell_types), fraction = 1)
dummy_plot <- ggplot(dummy_data, aes(x = cell_type, y = fraction, fill = cell_type)) +
              geom_bar(stat = "identity") +
              scale_fill_manual(values = col_cell_type) +
              theme_void() +
              theme(legend.position = "bottom")
legend <- get_plot_component(dummy_plot, "guide-box-bottom", return_all = TRUE)
white_matter_combined_plot <- plot_grid(pcol, legend, ncol = 1, rel_heights = c(1, .1))

pdf_path <- file.path(output_dir_human, paste0("barplot_cell_type_proportion/st_barplot_white_matter_cell_type_proportion_darmanis_annot", ".pdf"))
pdf(pdf_path, width = 8, height = 10)
print(white_matter_combined_plot)
dev.off()

### Statistical summary between samples

sce.ident <- "seurat_custom_clusters"

# Convert Seurat object to SingleCellExperiment
sce_orth_human <- as.SingleCellExperiment(seurat_orth_human, assay = "RNA")

# Convert colData to data frame
metadata_df <- colData(sce_orth_human) %>% as.data.frame()

# Initialise list for results
cell.type_spots_counts <- list()

# Loop over clusters, cell types, conditions, and time points
for (cluster in unique(metadata_df[[sce.ident]])) {
  for (cell_type in unique(metadata_df$labels)) {
    for (epil_condition in unique(metadata_df$condition)) {
      for (time_point in unique(metadata_df$time)) {
        # Filter and count spots
        barcodes <- metadata_df %>%
                    filter(!!sym(sce.ident) == cluster, time == time_point, 
                           condition == epil_condition, labels == cell_type) %>% rownames()
        cluster_size <- nrow(metadata_df %>% 
                             filter(!!sym(sce.ident) == cluster, condition == epil_condition, time == time_point))
        name_condition <- sprintf("%s_%s_%s_%s_%d", epil_condition, time_point, cluster, cell_type, cluster_size)
        cell.type_spots_counts[[name_condition]] <- barcodes
      }
    }
  }
}

# Convert the list to a dataframe
spots_stat_df <- do.call(rbind, lapply(names(cell.type_spots_counts), function(name_condition) {
  condition_info <- strsplit(name_condition, "_")[[1]]
  epil_condition <- condition_info[1]
  time_point <- condition_info[2]
  cluster <- condition_info[3]
  cell_type <- condition_info[4]
  total_size <- as.numeric(condition_info[length(condition_info)])
  
  spot_count <- length(cell.type_spots_counts[[name_condition]])
  
  data.frame(
    epil_condition = epil_condition,
    time_point = time_point,
    cluster = as.numeric(cluster),  # Ensure cluster is treated as numeric for sorting
    cell_type = cell_type,
    spot_count = spot_count,
    total_size = total_size,
    proportion = spot_count / total_size
  )
}))

# Sort the dataframe by the cluster column in increasing order
spots_stat_df <- spots_stat_df %>% arrange(cluster)

# Replace "-" by "_" in cell_type
spots_stat_df$cell_type <- gsub("-", "_", spots_stat_df$cell_type)

# Remove rownames
rownames(spots_stat_df) <- NULL

### z-test on cell-type proportions per cluster

# Initialize a list to store the results
z_test_results <- list()

# Loop through clusters, time points, and cell types
for (given.cluster in unique(spots_stat_df$cluster)) {
  for (time in unique(spots_stat_df$time_point)) {
    for (cell.type in unique(spots_stat_df$cell_type)) {
      
      # Filter data for the specified cluster, time_point, and cell_type
      filtered_data <- spots_stat_df %>%
                       filter(cluster == given.cluster, time_point == time, cell_type == cell.type)
      
      # Extract data for control and SE conditions
      ctrl_spot_stat <- filtered_data %>% filter(epil_condition == "CTRL") %>% select(spot_count, total_size) %>% unlist()
      se_spot_stat <- filtered_data %>% filter(epil_condition == "SE") %>% select(spot_count, total_size) %>% unlist()

      # Perform a proportion test, suppressing warnings
      z_test <- suppressWarnings(prop.test(x = c(se_spot_stat[1], ctrl_spot_stat[1]), 
                                           n = c(se_spot_stat[2], ctrl_spot_stat[2])))

      # Extract p-value and test name
      p_value <- ifelse(is.nan(z_test[["p.value"]]), NA, z_test[["p.value"]])
      test_name <- z_test[["method"]]
      
      # Store the result in a data frame and add to the results list
      z_test_result_df <- data.frame(
        cluster = given.cluster,
        time_point = time,
        cell_type = cell.type,
        test = test_name,
        p_value = p_value
      )
      z_test_results <- append(z_test_results, list(z_test_result_df))
    }
  }
}

# Combine all results into a single dataframe
z_test_results_df <- do.call(rbind, z_test_results)

# Remove NA values from p_values column
z_test_results_df <- z_test_results_df %>% filter(!is.na(p_value))

# Keep only significative p_value
z_test_results_signif_df <- z_test_results_df %>% filter(p_value <= 0.05)

# Save table
csv_path <- file.path(output_dir_human, "barplot_cell_type_proportion")
write.csv(z_test_results_signif_df, file = sprintf("%s/cell_types_per_cluster_proportion_test_signif_results.csv", csv_path), row.names = FALSE)