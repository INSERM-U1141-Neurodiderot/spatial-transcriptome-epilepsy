.libPaths(c(.libPaths(), "/home/christophe.lepriol/NeuroDev_ADD/R/r_4.1.0"))
library(ggplot2)
library(ggpubr)
library(cowplot)
library(patchwork)
library(RColorBrewer)
library(Seurat)



##############
# Parameters #
##############
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"
src_dir <- sprintf("%s/20-Count_analysis/10-EpiReg/src", work_dir)
genome_name <- "RefSeq108"

epireg_visium_metadata_df <- data.frame(sample=c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16"),
                                       condition=factor(c(rep("SE", 4), rep("CTRL", 4), rep("SE", 4), rep("CTRL", 4))), 
                                       time=factor(c(rep(c(5, 10), 4), rep(c(20, 40, 40, 20), 2))))

# sections
sections <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")
timepoint <- NA
#timepoint <- 5

# parameters
## count matrix: Space Ranger or ST Pipeline
input_matrix <- "Space Ranger" # allowed values: 'ST Pipeline' or 'Space Ranger'
## QC
spot_sum_of_counts_threshold <- 1000
spot_detected_genes_threshold <- 500
spot_mito_gene_count_pct_threshold <- 100
spot_hb_gene_count_pct_threshold <- 100
variables_to_regress_in_path <- "mito"
scale_data_split_by_in_path <- "orig"
integration_method <- "Seurat"
integration_group_by_in_path <- "orig"

prep_SCT_FindMarkers <- FALSE
prep_SCT_FindMarkers_subclustering <- FALSE

one_timepoint <- "5"
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))


#############
# Functions #
#############
#source(sprintf("%s/6X-integration_functions.R", src_dir))

marker_figure_plot <- function(seurat_obj, gene, palette, axis_text_x_size) {
    #### set color palette based on gene max expression
    gene_max_exp <- max(seurat_obj@assays$SCT@data[gene,])
    sc <- scale_fill_gradientn(colours=palette(100), limits=c(0, gene_max_exp))
    #### violin plot
    violin_plot <- VlnPlot(object=seurat_obj, features=gene)
    violin_plot <- violin_plot + labs(x="Cluster") + theme(legend.position="none") + theme(axis.text.x=element_text(angle=0, hjust=0.55))
    if (! is.na(axis_text_x_size)) {
        violin_plot <- violin_plot + theme(axis.text.x=element_text(size=axis_text_x_size))
    }
    #### expression on tissue section
    spatial_plot <- SpatialFeaturePlot(object=seurat_obj, features=gene, alpha=c(0.1,1), ncol=2)
    spatial_plot <- spatial_plot & sc
    ##### get legend
    spatial_plot_legend <- get_legend(spatial_plot)
    ##### spatial plot without legend
    spatial_plot <- spatial_plot & theme(plot.title=element_text(size=10), legend.position="none")
    ##### combine spatial plot and legend
    spatial_plot_figure <- ggdraw() +
        draw_plot(spatial_plot, 0, 0.15, 1, 0.85) +
        draw_plot(spatial_plot_legend, 0, 0, 1, 0.15)
    marker_figure <- violin_plot / spatial_plot_figure + plot_layout(heights=c(1,3))
    #svg(sprintf("%s/whole_clustering_figure_markers_cortex.svg", one_fig_output_dir, category_colname))
    #print(marker_figure)
    #dev.off()
    return(marker_figure)
}

############
# Analysis #
############

# process parameters
## sections
if (! is.na(timepoint)) {
    sections <- epireg_visium_metadata_df[which(epireg_visium_metadata_df$time==timepoint), "sample"]
}
nb_sections <- length(sections)

## input matrix
if (input_matrix == "ST Pipeline") {
    input_dirname <- "00-ST_P"
    visium_params_in_path <- TRUE
} else {
    if (input_matrix == "Space Ranger") {
        input_dirname <- "10-SR"
        visium_params_in_path <- FALSE
    }
}

output_dirname <- ifelse(is.na(timepoint), "all", sprintf("D%d", timepoint))

# load whole clustering Seurat object and list of subclustering Seurat objects
## build directory path
output_dir <- sprintf("%s/20-Count_analysis/10-EpiReg/output/%s/TSO_polyA_R1trim1_ov5_n2_min20/%s", work_dir, input_dirname, genome_name)
if (visium_params_in_path) {
    output_dir <- sprintf("%s/Visium_params", output_dir) 
}
#output_dir <- sprintf("%s/20-Integration/%s/sub_markers_20230527", output_dir, output_dirname)
integration_output_dir <- sprintf("%s/20-Integration/%s/sub_markers/Norm_%s-scale_%s-int_%s/QC-sum%d_det%s_mito%d_hb%d/10-QC_filtering/%s", output_dir, output_dirname, variables_to_regress_in_path, scale_data_split_by_in_path, integration_group_by_in_path, spot_sum_of_counts_threshold, spot_detected_genes_threshold, spot_mito_gene_count_pct_threshold, spot_hb_gene_count_pct_threshold, integration_method)
## load data: run PrepSCTFindMarkers if needed
#load(sprintf("%s/whole_clustering_Seurat_object_after_PrepSCTFindMarkers.RData", integration_output_dir))
if (prep_SCT_FindMarkers) {
    ## load Seurat objects: PrepSCTFindMarkers() function was not run
    load(sprintf("%s/whole_clustering_subclustering_Seurat_objects.RData", integration_output_dir))
    ## run seurat_object_whole_clustering()
    seurat_object_whole_clustering <- PrepSCTFindMarkers(seurat_object_whole_clustering)
    ## save whole clustering Seurat object after PrepSCTFindMarkers() call
    save(seurat_object_whole_clustering, file=sprintf("%s/whole_clustering_Seurat_object_after_PrepSCTFindMarkers.RData", integration_output_dir))
} else {
    ## load Seurat objects: PrepSCTFindMarkers() function was run
    load(sprintf("%s/whole_clustering_Seurat_object_after_PrepSCTFindMarkers.RData", integration_output_dir))
}


fig_output_dir <- sprintf("%s/20-Integration/%s/sub_markers/Norm_%s-scale_%s-int_%s/QC-sum%d_det%s_mito%d_hb%d/10-QC_filtering/%s", output_dir, output_dirname, variables_to_regress_in_path, scale_data_split_by_in_path, integration_group_by_in_path, spot_sum_of_counts_threshold, spot_detected_genes_threshold, spot_mito_gene_count_pct_threshold, spot_hb_gene_count_pct_threshold, integration_method)


# whole clustering
one_fig_output_dir <- sprintf("%s/whole/fig", fig_output_dir)
#pdf(sprintf("%s/whole_clustering_figure.pdf", one_fig_output_dir))
if (! dir.exists(one_fig_output_dir)) {
    dir.create(one_fig_output_dir, recursive=TRUE, mode="0775")
}

## UMAP plots
dimplot_list <- list()
for (category_colname in c("seurat_clusters", "condition", "time", "orig.ident")) {
    p <- DimPlot(seurat_object_whole_clustering, reduction="umap", group.by=category_colname, label=FALSE)
    if (length(levels(seurat_object_whole_clustering@meta.data[[category_colname]])) > 1) {
        p <- p + ggtitle(sprintf("Mean silhouette: %.2f", mean(seurat_object_whole_clustering@meta.data[[sprintf("silhouette_%s", category_colname)]]))) + theme(plot.title=element_text(size=12)) + theme(legend.position="bottom")
    }
    if (category_colname == "orig.ident") {
       p <- p + theme(legend.text=element_text(size=5)) +
        theme(legend.key.size=unit(1, 'mm'))
    }
    dimplot_list[[category_colname]] <- p
    svg(sprintf("%s/whole_clustering_figure_umap_%s.svg", one_fig_output_dir, category_colname))
    print(p)
    dev.off()
    pdf(sprintf("%s/whole_clustering_figure_umap_%s.pdf", one_fig_output_dir, category_colname))
    print(p)
    dev.off()
}
umap_figure <- (dimplot_list[["seurat_clusters"]] + dimplot_list[["condition"]]) / (dimplot_list[["time"]] + dimplot_list[["orig.ident"]])
svg(sprintf("%s/whole_clustering_figure_umap.svg", one_fig_output_dir))
print(umap_figure)
dev.off()
pdf(sprintf("%s/whole_clustering_figure_umap.pdf", one_fig_output_dir))
print(umap_figure)
dev.off()

## clustering on tissue sections
one_timepoint_metadata_df <- seurat_object_whole_clustering@meta.data[which(seurat_object_whole_clustering@meta.data$time==one_timepoint),]
one_timepoint_metadata_df <- droplevels(one_timepoint_metadata_df)
sections <- levels(one_timepoint_metadata_df$orig.ident)
### plot to get legend
clustering_splatial_plot <- SpatialDimPlot(seurat_object_whole_clustering, images=sections[1], label=TRUE, label.size=2)
clustering_splatial_plot <- clustering_splatial_plot + labs(fill="Cluster") + theme(legend.position="bottom")
clustering_splatial_plot_legend <- get_legend(clustering_splatial_plot)
### plot with all timepoint sections, without legend
clustering_splatial_plot <- SpatialDimPlot(seurat_object_whole_clustering, images=sections, ncol=2, label=TRUE, label.size=2) & labs(title="") & theme(legend.position="none")
### combine spatial plot and legend
clustering_splatial_plot_figure <- ggdraw() +
    draw_plot(clustering_splatial_plot, 0, 0.1, 1, 0.9) +
    draw_plot(clustering_splatial_plot_legend, 0, 0, 1, 0.1)
svg(sprintf("%s/whole_clustering_figure_spatial.svg", one_fig_output_dir))
print(clustering_splatial_plot_figure)
dev.off()
pdf(sprintf("%s/whole_clustering_figure_spatial.pdf", one_fig_output_dir))
print(clustering_splatial_plot_figure)
dev.off()

## marker figures
### subset Seurat object to two control samples
ctrl_sample_one_timepoint <- epireg_visium_metadata_df[which(epireg_visium_metadata_df$time==one_timepoint & epireg_visium_metadata_df$condition=="CTRL"), "sample"]
barcodes2use <- rownames(seurat_object_whole_clustering@meta.data[which(seurat_object_whole_clustering@meta.data$orig.ident %in% ctrl_sample_one_timepoint),])
seurat_object_whole_clustering_subset <- subset(seurat_object_whole_clustering, cells=barcodes2use)
#### after subsetting, only keep image of the samples in the subset
seurat_object_whole_clustering_subset@meta.data <- droplevels(seurat_object_whole_clustering_subset@meta.data)
seurat_object_whole_clustering_subset@images <- seurat_object_whole_clustering_subset@images[ctrl_sample_one_timepoint]

DefaultAssay(seurat_object_whole_clustering_subset) <- "SCT"

### cortex
one_gene <- "Tbr1"
marker_figure <- marker_figure_plot(seurat_object_whole_clustering_subset, one_gene, myPalette, NA)
svg(sprintf("%s/whole_clustering_figure_marker_cortex.svg", one_fig_output_dir))
print(marker_figure)
dev.off()
pdf(sprintf("%s/whole_clustering_figure_marker_cortex.pdf", one_fig_output_dir))
print(marker_figure)
dev.off()

### hippocampus
one_gene <- "Cabp7"
marker_figure <- marker_figure_plot(seurat_object_whole_clustering_subset, one_gene, myPalette, NA)
svg(sprintf("%s/whole_clustering_figure_marker_hippocampus.svg", one_fig_output_dir))
print(marker_figure)
dev.off()
pdf(sprintf("%s/whole_clustering_figure_marker_hippocampus.pdf", one_fig_output_dir))
print(marker_figure)
dev.off()

# subclusterings
one_fig_output_dir <- sprintf("%s/sub", fig_output_dir)
Idents(seurat_object_whole_clustering_subset) <- "seurat_subclusters"
subcluster_markers <- list("0"=c("Rprm", "Fezf2"), "1"=c("Ermn", "Mfrp"), "2"=c("Calb2", "Calb2"), "3"=c("Rgs14", "Kctd4", "Rspo2"), "4"=c("Adora2a", "Meis2"))

## load data: run PrepSCTFindMarkers if needed
seurat_object_subclustering_list2use <- list()
if (prep_SCT_FindMarkers_subclustering) {
    ## load Seurat objects: PrepSCTFindMarkers() function was not run
    load(sprintf("%s/whole_clustering_subclustering_Seurat_objects.RData", integration_output_dir))
    for (one_cluster in names(seurat_object_subclustering_list)) {
        seurat_object_subclustering_list2use[[one_cluster]] <- PrepSCTFindMarkers(seurat_object_subclustering_list[[one_cluster]])
    }
    ## save subclustering Seurat objects after PrepSCTFindMarkers() call
    save(seurat_object_subclustering_list2use, file=sprintf("%s/subclustering_Seurat_objects_after_PrepSCTFindMarkers.RData", integration_output_dir))
} else {
    ## load Seurat objects: PrepSCTFindMarkers() function was run
    load(sprintf("%s/subclustering_Seurat_objects_after_PrepSCTFindMarkers.RData", integration_output_dir))
}

for (one_cluster in names(seurat_object_subclustering_list2use)) {
#for (one_cluster in c("0", "3")) {
    one_cluster_fig_output_dir <- sprintf("%s/c%s/fig", one_fig_output_dir, one_cluster)
    if (! dir.exists(one_cluster_fig_output_dir)) {
        dir.create(one_cluster_fig_output_dir, recursive=TRUE, mode="0775")
    }
    ## UMAP plots
    dimplot_list <- list()
    for (category_colname in c("seurat_clusters", "condition", "time", "orig.ident")) {
        p <- DimPlot(seurat_object_subclustering_list2use[[one_cluster]], reduction="umap", group.by=category_colname, label=FALSE)
        if (length(levels(seurat_object_subclustering_list2use[[one_cluster]]@meta.data[[category_colname]])) > 1) {
            p <- p + ggtitle(sprintf("Mean silhouette: %.2f", mean(seurat_object_subclustering_list2use[[one_cluster]]@meta.data[[sprintf("silhouette_%s", category_colname)]]))) + theme(plot.title=element_text(size=12)) + theme(legend.position="bottom")
        }
        if (category_colname == "orig.ident") {
           p <- p + theme(legend.text=element_text(size=5)) +
            theme(legend.key.size=unit(1, 'mm'))
        }
        dimplot_list[[category_colname]] <- p
        svg(sprintf("%s/c%s_subclustering_figure_umap_%s.svg", one_cluster_fig_output_dir, one_cluster, category_colname))
        print(p)
        dev.off()
        pdf(sprintf("%s/c%s_subclustering_figure_umap_%s.pdf", one_cluster_fig_output_dir, one_cluster, category_colname))
        print(p)
        dev.off()
    }
    umap_figure <- (dimplot_list[["seurat_clusters"]] + dimplot_list[["condition"]]) / (dimplot_list[["time"]] + dimplot_list[["orig.ident"]])
    svg(sprintf("%s/c%s_subclustering_figure_umap.svg", one_cluster_fig_output_dir, one_cluster))
    print(umap_figure)
    dev.off()
    pdf(sprintf("%s/c%s_subclustering_figure_umap.pdf", one_cluster_fig_output_dir, one_cluster))
    print(umap_figure)
    dev.off()
    
    ## clustering on tissue sections
    one_timepoint_metadata_df <- seurat_object_subclustering_list2use[[one_cluster]]@meta.data[which(seurat_object_subclustering_list2use[[one_cluster]]@meta.data$time==one_timepoint),]
    one_timepoint_metadata_df <- droplevels(one_timepoint_metadata_df)
    sections <- levels(one_timepoint_metadata_df$orig.ident)
    ### plot to get legend
    clustering_splatial_plot <- SpatialDimPlot(seurat_object_subclustering_list2use[[one_cluster]], images=sections[1], label=TRUE, label.size=2)
    clustering_splatial_plot <- clustering_splatial_plot + labs(fill="Cluster") + theme(legend.position="bottom")
    clustering_splatial_plot_legend <- get_legend(clustering_splatial_plot)
    ### plot with all timepoint sections, without legend
    clustering_splatial_plot <- SpatialDimPlot(seurat_object_subclustering_list2use[[one_cluster]], images=sections, ncol=2, label=TRUE, label.size=2) & theme(plot.title=element_text(size=10), legend.position="none")
    ### combine spatial plot and legend
    clustering_splatial_plot_figure <- ggdraw() +
        draw_plot(clustering_splatial_plot, 0, 0.1, 1, 0.9) +
    draw_plot(clustering_splatial_plot_legend, 0, 0, 1, 0.1)
    svg(sprintf("%s/c%s_subclustering_figure_spatial.svg", one_cluster_fig_output_dir, one_cluster))
    print(clustering_splatial_plot_figure)
    dev.off()
    pdf(sprintf("%s/c%s_subclustering_figure_spatial.pdf", one_cluster_fig_output_dir, one_cluster))
    print(clustering_splatial_plot_figure)
    dev.off()
    
    ## markers
    one_cluster_markers <- subcluster_markers[[one_cluster]]
    for (one_marker in one_cluster_markers) {
        marker_figure <- marker_figure_plot(seurat_object_whole_clustering_subset, one_marker, myPalette, 7)
        svg(sprintf("%s/c%s_subclustering_figure_marker_%s.svg", one_cluster_fig_output_dir, one_cluster, one_marker))
        print(marker_figure)
        dev.off()
        pdf(sprintf("%s/c%s_subclustering_figure_marker_%s.pdf", one_cluster_fig_output_dir, one_cluster, one_marker))
        print(marker_figure)
        dev.off()
    }
}




# custom clustering, DEG and GO term counts
one_fig_output_dir <- sprintf("%s/custom", fig_output_dir)
Idents(seurat_object_whole_clustering) <- "seurat_custom_clusters"
## clustering on tissue sections
one_timepoint_metadata_df <- seurat_object_whole_clustering@meta.data[which(seurat_object_whole_clustering@meta.data$time==one_timepoint),]
one_timepoint_metadata_df <- droplevels(one_timepoint_metadata_df)
sections <- levels(one_timepoint_metadata_df$orig.ident)
### plot to get legend
clustering_splatial_plot <- SpatialDimPlot(seurat_object_whole_clustering, images=sections[1], label=TRUE, label.size=2)
clustering_splatial_plot <- clustering_splatial_plot + labs(fill="Cluster")
clustering_splatial_plot_legend <- get_legend(clustering_splatial_plot)
### plot with all timepoint sections, without legend
clustering_splatial_plot <- SpatialDimPlot(seurat_object_whole_clustering, images=sections, ncol=2, label=TRUE, label.size=2) & theme(plot.title=element_text(size=10), legend.position="none")
### combine spatial plot and legend
clustering_splatial_plot_figure <- ggdraw() +
    draw_plot(clustering_splatial_plot, 0, 0, 0.9, 1) +
    draw_plot(clustering_splatial_plot_legend, 0.9, 0, 0.1, 1)
svg(sprintf("%s/whole_clustering_figure_spatial_custom.svg", one_fig_output_dir))
print(clustering_splatial_plot_figure)
dev.off()
pdf(sprintf("%s/whole_clustering_figure_spatial_custom.pdf", one_fig_output_dir))
print(clustering_splatial_plot_figure)
dev.off()

## DEG counts: upset plots
### custom clustering
library(UpSetR)
library(dplyr)
one_fig_output_dir <- sprintf("%s/custom", fig_output_dir)
DEG_output_dir <- sprintf("%s/10-DEGs/negbinom/origident", one_fig_output_dir)
all_times_DE_results2compare_df <- read.csv(sprintf("%s/all_times_all_clusters_DE_analysis.csv", DEG_output_dir))
all_times_DE_results2compare_only_DEGs_df <- all_times_DE_results2compare_df[which(all_times_DE_results2compare_df$p_val_adj < 0.05),]

upset_figure_list <- list()
for (one_time in unique(all_times_DE_results2compare_df$time)) {
    print(sprintf("time: %d", one_time))
    all_times_DE_results2compare_only_DEGs_one_time_df <- all_times_DE_results2compare_only_DEGs_df[which(all_times_DE_results2compare_only_DEGs_df$time == one_time),]
    for (one_dysregulation in c("up", "down")) {
        print(sprintf("dysregulation: %s", one_dysregulation))
        if (one_dysregulation == "up") {
            all_times_DE_results2compare_only_DEGs_one_time_dys_df <- all_times_DE_results2compare_only_DEGs_one_time_df[which(all_times_DE_results2compare_only_DEGs_one_time_df$avg_log2FC > 0),]
        } else {
            if (one_dysregulation == "down") {
                all_times_DE_results2compare_only_DEGs_one_time_dys_df <- all_times_DE_results2compare_only_DEGs_one_time_df[which(all_times_DE_results2compare_only_DEGs_one_time_df$avg_log2FC < 0),]
            }
        }
        nb_clusters_with_DEGs <- length(unique(all_times_DE_results2compare_only_DEGs_one_time_dys_df$cluster))
        
        #### prepare data frame for upset plot
        upset_plot_df <- all_times_DE_results2compare_only_DEGs_one_time_dys_df %>%
            group_by(gene) %>%
            summarise(cluster = list(cluster))
        upset_plot_df2 <- upset_plot_df %>%
            group_by(gene) %>%
            summarise(cluster = unlist(lapply(cluster, function(x) {paste(x, collapse="&")})))
        expression_upset_input <- as.numeric(table(upset_plot_df2$cluster))
        names(expression_upset_input) <- names(table(upset_plot_df2$cluster))
        upset_plot <- upset(fromExpression(expression_upset_input), nsets=nb_clusters_with_DEGs, order.by = "freq", point.size=4, line.size=1.5)
        upset_figure_list[[sprintf("%s_%s", one_time, one_dysregulation)]] <- upset_plot
    }
}

library(grid)
pdf(sprintf("%s/all_times_all_clusters_DEG_upsets.pdf", DEG_output_dir), width=9)
for (one_time_one_dysregulation in names(upset_figure_list)) {
  print(upset_figure_list[[one_time_one_dysregulation]])
  grid.text(sprintf("%s", one_time_one_dysregulation), x = 0.65, y=0.95, gp=gpar(fontsize=14))
}
dev.off()


### whole clustering
one_fig_output_dir <- sprintf("%s/whole", fig_output_dir)
DEG_output_dir <- sprintf("%s/01-SCTransform_v2/qc_mito_percent_RNA/featnb2000/Seurat/10_dims/res0_03/10-DEGs/negbinom/origident/comp/", one_fig_output_dir)
all_times_DE_results2compare_df <- read.csv(sprintf("%s/SCTransform_v2_featnb2000_Seurat_10dims_res0_03_DEGs_all_times_all_clusters.csv", DEG_output_dir))
all_times_DE_results2compare_only_DEGs_df <- all_times_DE_results2compare_df[which(all_times_DE_results2compare_df$p_val_adj < 0.05),]

upset_figure_list <- list()
for (one_time in unique(all_times_DE_results2compare_df$time)) {
    print(sprintf("time: %d", one_time))
    all_times_DE_results2compare_only_DEGs_one_time_df <- all_times_DE_results2compare_only_DEGs_df[which(all_times_DE_results2compare_only_DEGs_df$time == one_time),]
    for (one_dysregulation in c("up", "down")) {
        print(sprintf("dysregulation: %s", one_dysregulation))
        if (one_dysregulation == "up") {
            all_times_DE_results2compare_only_DEGs_one_time_dys_df <- all_times_DE_results2compare_only_DEGs_one_time_df[which(all_times_DE_results2compare_only_DEGs_one_time_df$avg_log2FC > 0),]
        } else {
            if (one_dysregulation == "down") {
                all_times_DE_results2compare_only_DEGs_one_time_dys_df <- all_times_DE_results2compare_only_DEGs_one_time_df[which(all_times_DE_results2compare_only_DEGs_one_time_df$avg_log2FC < 0),]
            }
        }
        nb_clusters_with_DEGs <- length(unique(all_times_DE_results2compare_only_DEGs_one_time_dys_df$cluster))
        
        #### prepare data frame for upset plot
        upset_plot_df <- all_times_DE_results2compare_only_DEGs_one_time_dys_df %>%
            group_by(gene) %>%
            summarise(cluster = list(cluster))
        upset_plot_df2 <- upset_plot_df %>%
            group_by(gene) %>%
            summarise(cluster = unlist(lapply(cluster, function(x) {paste(x, collapse="&")})))
        expression_upset_input <- as.numeric(table(upset_plot_df2$cluster))
        names(expression_upset_input) <- names(table(upset_plot_df2$cluster))
        upset_plot <- upset(fromExpression(expression_upset_input), nsets=nb_clusters_with_DEGs, order.by = "freq", point.size=4, line.size=1.5)
        upset_figure_list[[sprintf("%s_%s", one_time, one_dysregulation)]] <- upset_plot
    }
}

pdf(sprintf("%s/all_times_all_clusters_DEG_upsets.pdf", DEG_output_dir), width=9)
for (one_time_one_dysregulation in names(upset_figure_list)) {
  print(upset_figure_list[[one_time_one_dysregulation]])
  grid.text(sprintf("%s", one_time_one_dysregulation), x = 0.65, y=0.95, gp=gpar(fontsize=14))
}
dev.off()



