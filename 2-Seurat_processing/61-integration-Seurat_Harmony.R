.libPaths(c(.libPaths(), "/home/christophe.lepriol/NeuroDev_ADD/R/r_4.1.0"))
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(cluster) # silhouette() function
library(cowplot)
#library(kBET)
library(harmony)

##############
# Parameters #
##############
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"
genome_name <- "NCBIRefSeq108_NCBIRefSeq108GTF"
load_image <- TRUE

epireg_visium_metadata_df <- data.frame(sample=c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16"),
                                       condition=factor(c(rep("SE", 4), rep("CTRL", 4), rep("SE", 4), rep("CTRL", 4))), 
                                       time=factor(c(rep(c(5, 10), 4), rep(c(20, 40, 40, 20), 2))))

# sections
sections <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")
#sections <- c("C_L1_S3", "C_L4_S15", "D_L1_S4", "D_L4_S16")
#sections <- c("A_L3_S9", "B_L3_S10")
#timepoint <- NA
timepoint <- 20
if (! is.na(timepoint)) {
    sections <- epireg_visium_metadata_df[which(epireg_visium_metadata_df$time==timepoint), "sample"]
}
nb_sections <- length(sections)

# output directory name
if (! is.na(timepoint)) {
    output_dirname <- sprintf("%dd_samples", timepoint)
} else {
    #if (nb_sections <= 4) {
    #    output_dirname <- paste(sections, collapse="_")
    #}
    output_dirname <- "all_samples"
    #output_dirname <- "SE_samples"
    #output_dirname <- "Ctrl_samples"
    #output_dirname <- "ABC_D_L1_L2_L3_samples_RAM_no_images_2000"
}
output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/20-Integration/%s", work_dir, genome_name, output_dirname)
if (! dir.exists(output_dir)) {
    dir.create(output_dir, recursive=TRUE, mode="0775")
}

# parameters
#nb_integration_features <- 2000
#feature_residual_variance_threshold <- 1.1
#feature_residual_variance_threshold <- NA
#feature_residual_variance_threshold_vector <- c(1, 1.1, 1.3, 1.5, 2, 2.5)
#features_param_vector <- c(1.1, 1.3, 1.5, 2)
#features_param_vector <- c(500, 1000, 1500, 2000, 2500)
#features_param_vector <- c(250, 500, 1000, 1500, 2000)
#features_param_vector <- c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000)
#features_param_vector <- c(500, 1000)
features_param_vector <- c(500)
#dimensions_nb_param_vector <- c(10, 15, 20)
#dimensions_nb_param_vector <- c(10, 15)
dimensions_nb_param_vector <- c(20)
#clustering_resolution_vector <- c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2)
#clustering_resolution_vector <- c(0.4, 0.6, 0.8, 1)
#clustering_resolution_vector <- c(0.4, 1)
#clustering_resolution_vector <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)
#clustering_resolution_vector <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)
#clustering_resolution_vector <- c(0.1, 0.2, 0.4, 0.5, 0.8, 1)
clustering_resolution_vector <- c(0.4)
#features_param <- 500
#dimensions_nb_param <- 20
#clustering_resolution <- 0.1
cluster_markers <- FALSE
DE_analysis <- FALSE
conserved_markers <- FALSE

# number of samples per page when using SpatialFeaturePlot() function
nb_SpatialFeaturePlot_per_page <- 2
# number of samples per page when using SpatialFeaturePlot() function
nb_cluster_markers_to_plot <- 4

nb_DEGs_to_plot <- 30
nb_DEGs_per_page <- 3

metadata_categories <- c("lame", "zone", "condition", "orig.ident")

#############
# Functions #
#############
load_data <- function(sample, st_dir, sr_dir, add_image) {
    # expression matrices
    ## ST Pipeline
    st_pipeline_matrix_file <- sprintf("%s/%s/10-Pipeline/%s_stdata.tsv", st_dir, sample, sample)
    st_pipeline_matrix <- read.table(st_pipeline_matrix_file, sep="\t", header=TRUE, quote="", row.names=1)
    # Space Ranger output directory for tissue positions list and image
    space_ranger_output_dir <- sprintf("%s/%s/10-Pipeline/outs", sr_dir, sample)

    # build SpatialExperiment object
    spatial_coordinates_file <- file.path(space_ranger_output_dir, "spatial", "tissue_positions_list.csv")
    spatial_coordinates_df <- read.csv(spatial_coordinates_file, header=FALSE, quote="")
    colnames(spatial_coordinates_df) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    rownames(spatial_coordinates_df) <- spatial_coordinates_df$barcode
    nb_spots <- dim(spatial_coordinates_df)[1]
    ## add lame, zone and condition
    sample_strsplit <- unlist(strsplit(sample, "_"))
    zone <- sample_strsplit[1]
    lame <- sample_strsplit[2]
    condition <- ifelse(zone=="A" | zone=="C", "SE", "CTRL")
    spatial_coordinates_df <- cbind(spatial_coordinates_df, lame=rep(lame, nb_spots), zone=rep(zone, nb_spots), condition=rep(condition, nb_spots))
    spatial_coordinates_df$zone <- as.factor(spatial_coordinates_df$zone)
    spatial_coordinates_df$lame <- as.factor(spatial_coordinates_df$lame)
    spatial_coordinates_df$condition <- as.factor(spatial_coordinates_df$condition)
    ## matrix: convert coordinates to barcodes
    barcodes <- unlist(lapply(rownames(st_pipeline_matrix), function(x, coord2barcode=spatial_coordinates_df) {
        coordinates <- unlist(strsplit(x, "x"))
        return(coord2barcode[which(coord2barcode$array_col==as.integer(coordinates[1])-1 & coord2barcode$array_row==as.integer(coordinates[2])-1), "barcode"])
    }))
    rownames(st_pipeline_matrix) <- barcodes
    
    # create Seurat object
    seurat_obj <- CreateSeuratObject(counts=t(st_pipeline_matrix), project=sample, meta.data=spatial_coordinates_df)
    ## add image
    if (add_image) {
        img <- Read10X_Image(image.dir=file.path(space_ranger_output_dir, "spatial"))
        img@assay <- Assays(seurat_obj) # set assay identical to Seurat object
        Key(img) <- sprintf("%s_", gsub("_", "", sample)) # add key for SpatialFeaturePlot() and SpatialPlot() functions: Keys should be one or more alphanumeric characters followed by an underscore
        img_list <- list(img) # the image must be stored in a list
        names(img_list) <- sample # change image name in the list: used as title by SpatialFeaturePlot() function for example
        seurat_obj@images <- img_list
    }
    ## keep spot over tissue
    seurat_obj <- seurat_obj[, seurat_obj$in_tissue==1]
    
    return(seurat_obj)
}

compute_silhouette <- function(seurat_obj, reduc, dim_nb, metadata_colname) {
    distance_matrix <- dist(Embeddings(seurat_obj, reduction=reduc)[, 1:dim_nb])
    categories <- as.factor(seurat_obj@meta.data[[metadata_colname]])
    silhouette_obj <- silhouette(as.numeric(categories), dist=distance_matrix)
    seurat_obj@meta.data[[sprintf("silhouette_%s", metadata_colname)]] <- silhouette_obj[,3]
    return(seurat_obj)
}

clustering <- function(seurat_obj_int, rez, reduc, dim_nb) {
    # compute clustering
    seurat_obj_int <- FindClusters(seurat_obj_int, resolution=rez, verbose=FALSE)
    seurat_obj_int <- RunUMAP(seurat_obj_int, reduction=reduc, dims=1:dim_nb)
    # comupute silhouette scores
    seurat_obj_int <- compute_silhouette(seurat_obj_int, reduc, dim_nb, "seurat_clusters")
    
    # kBET
    ### clusters
    #kBET_values <- kBET(Embeddings(seurat_objects_merge$pca), seurat_objects_merge@meta.data$seurat_clusters)
    ### lame
    #### 50 PCs
    #kBET_values_lame <- kBET(Embeddings(seurat_objects_merge$pca), seurat_objects_merge@meta.data$lame)
    #### 20 PCs
    #kBET_values_inetgrated_lame_20PCs <- kBET(Embeddings(seurat_obj_int$pca)[,1:20], seurat_obj_int@meta.data$lame)

    return(seurat_obj_int)
}

clustering_plots <- function(seurat_obj_int, load_image, plots_per_page, plot_title, out_dir, out_name) {
    # plot clusters onto UMAP or onto the tissue section
    pdf(sprintf("%s/%s.pdf", out_dir, out_name))
    dimplot_list <- list()
    for (category_colname in c("seurat_clusters", "lame", "zone", "condition", "orig.ident")) {
        p <- DimPlot(seurat_obj_int, reduction="umap", group.by=category_colname, label=TRUE)
        if (length(levels(seurat_obj_int@meta.data[[category_colname]])) > 1) {
            p <- p + ggtitle(sprintf("Mean silhouette for %s: %f", category_colname, mean(seurat_obj_int@meta.data[[sprintf("silhouette_%s", category_colname)]]))) + theme(plot.title=element_text(size=14))
        }
        dimplot_list[[category_colname]] <- p
    }
    print(dimplot_list[["seurat_clusters"]] / dimplot_list[["lame"]] + plot_annotation(title=sprintf("%s", plot_title), theme=theme(plot.title=element_text(size=16)))) 
    print(dimplot_list[["seurat_clusters"]] / dimplot_list[["zone"]] + plot_annotation(title=sprintf("%s", plot_title), theme=theme(plot.title=element_text(size=16))))
    print(dimplot_list[["seurat_clusters"]] / dimplot_list[["condition"]] + plot_annotation(title=sprintf("%s", plot_title), theme=theme(plot.title=element_text(size=16))))
    print(dimplot_list[["seurat_clusters"]] / dimplot_list[["orig.ident"]] + plot_annotation(title=sprintf("%s", plot_title), theme=theme(plot.title=element_text(size=16))))
    
    # silhouette: https://romanhaa.github.io/projects/scrnaseq_workflow/#silhouette-plot
    mean_silhouette_score <- mean(seurat_obj_int@meta.data[["silhouette_seurat_clusters"]])
    df2ggplot <- seurat_obj_int@meta.data %>%
      mutate(barcode = rownames(.)) %>%
      arrange(seurat_clusters,-silhouette_seurat_clusters) %>%
      mutate(barcode = factor(barcode, levels = barcode))
    write.csv(df2ggplot, file=sprintf("%s/%s_silhouette.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
    
    silhouette_score_plot <- ggplot(df2ggplot) +
      geom_col(aes(barcode, silhouette_seurat_clusters, fill = seurat_clusters), show.legend = FALSE) +
      geom_hline(yintercept = mean_silhouette_score, color = 'red', linetype = 'dashed') +
      labs(title=sprintf("Mean silhouette score: %f", mean_silhouette_score)) +
      scale_x_discrete(name = 'Cells') +
      scale_y_continuous(name = 'Silhouette score') +
      #scale_fill_manual(values = custom_colors$discrete) +
      theme_bw() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    print(silhouette_score_plot)
    
    if (load_image) {
        nb_sections <- length(unique(seurat_obj_int@meta.data$orig.ident))
        for (i in 1:(nb_sections/plots_per_page)) {
            print(SpatialDimPlot(seurat_obj_int, images=names(seurat_obj_int@images)[((i-1)*plots_per_page+1):(i*plots_per_page)], label=TRUE, label.size=2) + plot_annotation(title=sprintf("%s", plot_title), theme=theme(plot.title=element_text(size=18))))
        }
        for (i in 1:nb_sections) {
            print(SpatialDimPlot(seurat_obj_int, images=names(seurat_obj_int@images)[i], label=TRUE, label.size=3) + plot_annotation(title=sprintf("%s: %s", plot_title, names(seurat_obj_int@images)[i]), theme=theme(plot.title=element_text(size=16))))
        }
    }
    dev.off()
}

identify_cluster_markers <- function(seurat_obj_int, assay2use, load_image, nb_plots, plots_per_page, out_dir, out_name) {
    print("Identify cluster markers")
    # set default assay: must be either 'RNA' or 'SCT'
    DefaultAssay(seurat_obj_int) <- assay2use
    markers <- FindAllMarkers(seurat_obj_int, logfc.threshold=0.25, min.pct=0.25, only.pos=TRUE)
    write.csv(markers, file=sprintf("%s/%s.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
    pdf(sprintf("%s/%s.pdf", out_dir, out_name))
    ## expression heatmap
    top10_markers <- markers %>%
        group_by(cluster) %>%
        top_n(n=10, wt=avg_log2FC)
    print(DoHeatmap(seurat_obj_int, features=top10_markers$gene) + NoLegend())
    ## get top markers for each cluster
    if (load_image) {
        top_markers <- markers %>%
            group_by(cluster) %>%
            slice_max(n=nb_plots, order_by=avg_log2FC)
        nb_sections <- length(unique(seurat_obj_int@meta.data$orig.ident))
        for (one_cluster in levels(markers$cluster)) {
            print(sprintf("cluster: %s", one_cluster))
            one_cluster_markers <- top_markers[which(top_markers$cluster==one_cluster), "gene", drop=TRUE]
            for (i in 1:(nb_plots/2)) {
                for (j in 1:(nb_sections/plots_per_page)) {
                    print(SpatialFeaturePlot(object=seurat_obj_int, features=one_cluster_markers[((i-1)*2+1):(i*2)], images=names(seurat_obj_int@images)[((j-1)*plots_per_page+1):(j*plots_per_page)], alpha=c(0.1,1), ncol=3) + plot_annotation(title=sprintf("Cluster: %s", one_cluster), theme=theme(plot.title=element_text(size=16))))
                }
            }
        }
    }
    dev.off()
}

identify_DEGs <- function(seurat_obj_int, assay2use, cluster_col, condition_col, out_dir, out_name) {
    print("Identify DEGs")
    # set default assay: must be either 'RNA' or 'SCT'
    DefaultAssay(seurat_obj_int) <- assay2use
    # average expression plots
    print("average expression plots")
    theme_set(theme_cowplot())
    ## set identity classes to clusters
    Idents(seurat_obj_int) <- seurat_obj_int[[cluster_col]]
    seurat_obj_int.all_clusters.avg <- data.frame()
    pdf(sprintf("%s/%s_all_clusters_average_expression.pdf", out_dir, out_name))
    for (one_cluster in levels(seurat_obj_int[[cluster_col]])) {
        print(sprintf("cluster: %s", one_cluster))
        seurat_obj_int.one_cluster <- subset(seurat_obj_int, idents=one_cluster)
        ## set identity classes to conditions
        Idents(seurat_obj_int.one_cluster) <- condition_col
        seurat_obj_int.one_cluster.avg <- as.data.frame(log1p(AverageExpression(seurat_obj_int.one_cluster, verbose=FALSE)[[assay2use]]))
        seurat_obj_int.one_cluster.avg$gene <- rownames(seurat_obj_int.one_cluster.avg)
        seurat_obj_int.one_cluster.avg$diff <- seurat_obj_int.one_cluster.avg$SE - seurat_obj_int.one_cluster.avg$CTRL
        seurat_obj_int.all_clusters.avg <- rbind(seurat_obj_int.all_clusters.avg, cbind(seurat_obj_int.one_cluster.avg, cluster=rep(one_cluster, dim(seurat_obj_int.one_cluster.avg)[1])))
        genes.to.label <- c(rownames(seurat_obj_int.one_cluster.avg[order(seurat_obj_int.one_cluster.avg$diff, decreasing=TRUE),])[1:10], tail(rownames(seurat_obj_int.one_cluster.avg[order(seurat_obj_int.one_cluster.avg$diff, decreasing=TRUE),]), n=10))
        p <- ggplot(seurat_obj_int.one_cluster.avg, aes(CTRL, SE)) + geom_point() + ggtitle(sprintf("Cluster %s", one_cluster))
        p <- LabelPoints(p, points=genes.to.label, repel=TRUE)
        print(p)
    }
    dev.off()
    write.csv(seurat_obj_int.all_clusters.avg, file=sprintf("%s/%s_all_clusters_average_expression.csv", out_dir), quote=FALSE, row.names=TRUE)

    # differentially expressed genes in different conditions for cells of the same type
    print("differentially expressed genes between conditions")
    ## column in metadata table to hold both the cell type and stimulation information and switch the current ident to that column
    seurat_obj_int[[sprintf("%s.%s", cluster_col, condition_col)]] <- paste(Idents(seurat_obj_int), seurat_obj_int[[condition_col]], sep="_")
    ## set identity classes to clusters and conditions
    Idents(seurat_obj_int) <- seurat_obj_int[[sprintf("%s.%s", cluster_col, condition_col)]]
    ## differentially expressed genes between SE and CTRL spots
    seurat_obj_int.DE.all_clusters <- data.frame()
    for (one_cluster in levels(seurat_obj_int[[cluster_col]])) {
        print(sprintf("cluster: %s", one_cluster))
        seurat_obj_int.DE.one_cluster <- FindMarkers(seurat_obj_int, ident.1=sprintf("%s_SE", one_cluster), ident.2=sprintf("%s_CTRL", one_cluster), verbose=FALSE)
        write.csv(seurat_obj_int.DE.one_cluster, file=sprintf("%s/%s_cluster_%s.csv", out_dir, out_name, one_cluster), quote=FALSE, row.names=TRUE)
        seurat_obj_int.DE.all_clusters <- rbind(seurat_obj_int.DE.all_clusters, cbind(seurat_obj_int.DE.one_cluster, cluster=rep(one_cluster, dim(seurat_obj_int.DE.one_cluster)[1])))
        #### visualizations
        pdf(sprintf("%s/%s_cluster_%s.pdf", out_dir, out_name, one_cluster))
        for(i in seq(1, nb_DEGs_to_plot, nb_DEGs_per_page)) {
            print(FeaturePlot(seurat_obj_int, features=rownames(seurat_obj_int.DE.one_cluster)[i:(i+nb_DEGs_per_page-1)], split.by=condition_col, max.cutoff=3, cols=c("grey", "red")))
            plots <- VlnPlot(seurat_obj_int, features=rownames(seurat_obj_int.DE.one_cluster)[i:(i+nb_DEGs_per_page-1)], split.by=condition_col, group.by=cluster_col, pt.size=0, combine=FALSE)
            print(wrap_plots(plots=plots, ncol=1))
        }
        dev.off()
    }
    write.csv(seurat_obj_int.DE.all_clusters, file=sprintf("%s/%s_all_clusters.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
}

identify_conserved_markers <- function(seurat_obj_int, assay2use, cluster_col, condition_col, out_dir, out_name) {
    print("Identify conserved markers")
    # set default assay: must be either 'RNA' or 'SCT'
    DefaultAssay(seurat_obj_int) <- assay2use
    # set identity classes to clusters
    Idents(seurat_obj_int) <- seurat_obj_int[[cluster_col]]
    pdf(sprintf("%s/%s_umap_dotplots.pdf", out_dir, out_name))
    for (one_cluster in levels(Idents(seurat_obj_int))) {
        print(sprintf("cluster: %s", one_cluster))
        cluster_conserved_markers <- FindConservedMarkers(seurat_obj_int, ident.1=one_cluster, grouping.var=condition_col, verbose=FALSE)
        write.csv(cluster_conserved_markers, file=sprintf("%s/%s_cluster_%s.csv", out_dir, out_name, one_cluster), quote=FALSE, row.names=TRUE)
        print(FeaturePlot(seurat_obj_int, features=rownames(cluster_conserved_markers)[1:9], min.cutoff="q9") + plot_annotation(title=sprintf("Cluster %s", one_cluster), theme=theme(plot.title=element_text(size=16))))
        print(DotPlot(seurat_obj_int, features=rownames(cluster_conserved_markers)[1:30], cols=c("blue", "red"), dot.scale=8, split.by=condition_col) + RotatedAxis() + plot_annotation(title=sprintf("Cluster %s", one_cluster), theme=theme(plot.title=element_text(size=16))))
    }
    dev.off()
}


############
# Analysis #
############

integration_batch_mean_silhouette_df <- data.frame()
integration_cluster_mean_silhouette_df <- data.frame()

# hippocampus manual annotations
all_sections_hippocampus_barcodes_Seurat_merge <- c()
all_sections_hippocampus_barcodes_Seurat_integrated <- c()
sample_number <- 1
for (sample in sections) {
    space_ranger_output_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs", work_dir, genome_name, sample)
    manual_annotations_file <- sprintf("%s/manual_annotations.csv", space_ranger_output_dir)
    manual_annotations_df <- read.csv(file=manual_annotations_file)
    hippocampus_barcodes <- manual_annotations_df[which(manual_annotations_df$manual.annotations=="hippocampus"), "Barcode"]
    all_sections_hippocampus_barcodes_Seurat_merge <- c(all_sections_hippocampus_barcodes_Seurat_merge, sprintf("%s_%s", sample, hippocampus_barcodes))
    all_sections_hippocampus_barcodes_Seurat_integrated <- c(all_sections_hippocampus_barcodes_Seurat_integrated, sprintf("%s_%d", hippocampus_barcodes, sample_number))
    sample_number <- sample_number + 1
}

# load data
print("load data")
st_pipeline_dir <- sprintf("%s/10-ST_analysis/00-ST_Pipeline/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples", work_dir, genome_name)
space_ranger_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples", work_dir, genome_name)
seurat_objects <- lapply(sections, load_data, st_dir=st_pipeline_dir, sr_dir=space_ranger_dir, add_image=load_image)
names(seurat_objects) <- sections

# quality control
## merge datasets
print("merge datasets")
seurat_objects_merge <- merge(seurat_objects[[1]], y=unlist(seurat_objects)[2:nb_sections], add.cell.ids=sections)
if (load_image) {
    print("quality control")
    ## statistics: number of counts, number of features
    pdf(sprintf("%s/QC.pdf", output_dir))
    VlnPlot(seurat_objects_merge, features=c("nCount_RNA", "nFeature_RNA"), pt.size=0.1, ncol=2) + NoLegend()
    ### SpatialFeaturePlot: impossible to change legend parameters due to combine=TRUE
    #### SpatialFeaturePlot(seurat_objects_merge, features=c("nCount_RNA", "nFeature_RNA")) + theme(legend.title=element_text(size=8), legend.text=element_text(size=6)): the parameters are only applied to the last of plot
    #### only solution to make the legend readable: plot only 2 samples per page
    for (i in 1:(nb_sections/nb_SpatialFeaturePlot_per_page)) {
        print(SpatialFeaturePlot(seurat_objects_merge, features=c("nCount_RNA", "nFeature_RNA"), images=names(seurat_objects_merge@images)[((i-1)*nb_SpatialFeaturePlot_per_page+1):(i*nb_SpatialFeaturePlot_per_page)]))
    }
    ## filter
    ### select all spots with less than 25% mitochondrial reads, less than 20% hb-reads and more than 500 detected genes
    seurat_objects_merge <- seurat_objects_merge[, seurat_objects_merge$nFeature_RNA > 500]
    for (i in 1:(nb_sections/nb_SpatialFeaturePlot_per_page)) {
        print(SpatialFeaturePlot(seurat_objects_merge, features=c("nCount_RNA", "nFeature_RNA"), images=names(seurat_objects_merge@images)[((i-1)*nb_SpatialFeaturePlot_per_page+1):(i*nb_SpatialFeaturePlot_per_page)]))
    }
    dev.off()
}

# non-integrated datasets
integration_method <- "no integration"
for (features_param in features_param_vector) {
    ## normalization
    if (features_param > 10) {
        ### features_param is a number of features
        print(sprintf("features param: %d", features_param))
        normalization_output_dir <- sprintf("%s/00-SCTransform/featurenumber%d", output_dir, features_param)
        if (! dir.exists(normalization_output_dir)) {
            dir.create(normalization_output_dir, recursive=TRUE, mode="0775")
        }
        normalization_output_name <- sprintf("SCTransform_featurenumber%d", features_param)
        seurat_objects_merge <- SCTransform(seurat_objects_merge, assay="RNA", verbose=TRUE, method="poisson", variable.features.n=features_param, return.only.var.genes=FALSE)
    } else {
        ### features_param is a residual variance threshold
        print(sprintf("features param: %.1f", features_param))
        normalization_output_dir <- sprintf("%s/00-SCTransform/residualvariance%s", output_dir, sub("[.]", "_", features_param))
        if (! dir.exists(normalization_output_dir)) {
            dir.create(normalization_output_dir, recursive=TRUE, mode="0775")
        }
        normalization_output_name <- sprintf("SCTransform_residualvariance%s", sub("[.]", "_", features_param))
        seurat_objects_merge <- SCTransform(seurat_objects_merge, assay="RNA", verbose=TRUE, method="poisson", variable.features.n=NULL, variable.features.rv.th=features_param, return.only.var.genes=FALSE)
    }
    
    ## dimensionality reduction
    seurat_objects_merge <- RunPCA(seurat_objects_merge, verbose=FALSE)
    reduction2use <- "pca"
    for (dimensions_nb_param in dimensions_nb_param_vector) {
        print(sprintf("dimension number: %d", dimensions_nb_param))
        dimensions_output_dir <- sprintf("%s/00-No_integration/%d_PCs", normalization_output_dir, dimensions_nb_param)
        if (! dir.exists(dimensions_output_dir)) {
            dir.create(dimensions_output_dir, recursive=TRUE, mode="0775")
        }
        dimensions_output_name <- sprintf("%s_noIntegration_%dPCs", normalization_output_name, dimensions_nb_param)
        
        seurat_objects_merge <- FindNeighbors(seurat_objects_merge, reduction=reduction2use, dims=1:dimensions_nb_param)
        
        ### compute silhouette scores for lames, zones, conditions and samples
        mean_silhouette_vector <- c()
        for (metadata_category_colname in metadata_categories) {
            seurat_objects_merge <- compute_silhouette(seurat_objects_merge, reduction2use, dimensions_nb_param, metadata_category_colname)
            mean_silhouette <- mean(seurat_objects_merge@meta.data[[sprintf("silhouette_%s", metadata_category_colname)]])
            mean_silhouette_vector <- c(mean_silhouette_vector, mean_silhouette)
        }
        mean_silhouette_df <- as.data.frame(t(mean_silhouette_vector))
        colnames(mean_silhouette_df) <- sprintf("mean_silhouette_%s", metadata_categories)
        mean_silhouette_df <- cbind(normalization="SCTransform", features=features_param, integration=integration_method, dimensions=dimensions_nb_param, mean_silhouette_df)
        integration_batch_mean_silhouette_df <- rbind(integration_batch_mean_silhouette_df, mean_silhouette_df)
        
        ## clustering
        for (clustering_resolution in clustering_resolution_vector) {
            print(sprintf("cluster resolution: %.1f", clustering_resolution))
            resolution_output_dir <- sprintf("%s/resolution%s", dimensions_output_dir, sub("[.]", "_", clustering_resolution))
            if (! dir.exists(resolution_output_dir)) {
                dir.create(resolution_output_dir, recursive=TRUE, mode="0775")
            }
            resolution_output_name <- sprintf("%s_resolution%s", dimensions_output_name, sub("[.]", "_", clustering_resolution))
            
            seurat_objects_merge <- clustering(seurat_objects_merge, clustering_resolution, reduction2use, dimensions_nb_param)
            write.csv(seurat_objects_merge@meta.data, file=sprintf("%s/%s_metadata.csv", resolution_output_dir, resolution_output_name), quote=FALSE, row.names=TRUE)
            cluster_mean_silhouette <- mean(seurat_objects_merge$silhouette_seurat_clusters)
            ### plot clusters onto UMAP or onto the tissue section
            clustering_plots(seurat_objects_merge, load_image, nb_SpatialFeaturePlot_per_page, "Non-integrated datasets, all spots", resolution_output_dir, sprintf("%s_clustering", resolution_output_name))
            
            ### hippocampus spots
            seurat_objects_merge_hippocampus <- subset(seurat_objects_merge, cells=all_sections_hippocampus_barcodes_Seurat_merge)
            hippocampus_cluster_mean_silhouette <- mean(seurat_objects_merge_hippocampus$silhouette_seurat_clusters)
            #### plot clusters onto UMAP or onto the tissue section
            clustering_plots(seurat_objects_merge_hippocampus, load_image, nb_SpatialFeaturePlot_per_page, "Non-integrated datasets, hippocampus spots", resolution_output_dir, sprintf("%s_clustering_hippocampus", resolution_output_name))
            
            ### get mean silhouette scores: all spots and hippocampus spots
            mean_silhouette_df <- data.frame(normalization="SCTransform", features=features_param, integration=integration_method, dimensions=dimensions_nb_param, resolution=clustering_resolution, mean_silhouette_seurat_clusters=cluster_mean_silhouette, mean_silhouette_seurat_clusters_hippocampus=hippocampus_cluster_mean_silhouette)
            integration_cluster_mean_silhouette_df <- rbind(integration_cluster_mean_silhouette_df, mean_silhouette_df)
        }
    }
}
rm(seurat_objects_merge)
rm(seurat_objects_merge_hippocampus)
gc()


# integrated datasets
## Seurat
### need to set maxSize to avoid error when identifying conserved markers: Error in getGlobalsAndPackages(expr, envir = envir, globals = globals): The total size of the 3 globals exported for future expression (‘FUN()’) is 1.38 GiB.. This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize').
#options(future.globals.maxSize=2000 * 1024^2) # set allowed size to 2K MiB
integration_method <- "Seurat"
for (features_param in features_param_vector) {
    # normalization
    if (features_param > 10) {
        ## features_param is a number of features
        print(sprintf("features param: %d", features_param))
        normalization_output_dir <- sprintf("%s/00-SCTransform/featurenumber%d", output_dir, features_param)
        if (! dir.exists(normalization_output_dir)) {
            dir.create(normalization_output_dir, recursive=TRUE, mode="0775")
        }
        normalization_output_name <- sprintf("SCTransform_featurenumber%d", features_param)
        seurat_objects <- lapply(seurat_objects, FUN=SCTransform, variable.features.n=features_param, return.only.var.genes=FALSE)
        integration_features_nb <- features_param
    } else {
        ## features_param is a residual variance threshold
        print(sprintf("features param: %.1f", features_param))
        normalization_output_dir <- sprintf("%s/00-SCTransform/residualvariance%s", output_dir, sub("[.]", "_", features_param))
        if (! dir.exists(normalization_output_dir)) {
            dir.create(normalization_output_dir, recursive=TRUE, mode="0775")
        }
        normalization_output_name <- sprintf("SCTransform_residualvariance%s", sub("[.]", "_", features_param))
        seurat_objects <- lapply(seurat_objects, FUN=SCTransform, variable.features.n=NULL, variable.features.rv.th=features_param, return.only.var.genes=FALSE)
        ## get minimum number of variable features
        min_variable_features <- length(seurat_objects[[sections[1]]]@assays$SCT@var.features)
        for (sample in sections[2:length(sections)]) {
            if (length(seurat_objects[[sample]]@assays$SCT@var.features) < min_variable_features) {
                min_variable_features <- length(seurat_objects[[sample]]@assays$SCT@var.features)
            }
        }
        integration_features_nb <- min_variable_features
    }

    ## integration
    integration_output_dir <- sprintf("%s/10-Seurat", normalization_output_dir)
    if (! dir.exists(integration_output_dir)) {
        dir.create(integration_output_dir, recursive=TRUE, mode="0775")
    }
    integration_output_name <- sprintf("%s_%s", normalization_output_name, integration_method)
    features <- SelectIntegrationFeatures(seurat_objects, nfeatures=integration_features_nb)
    seurat_objects <- PrepSCTIntegration(object.list=seurat_objects, anchor.features=features)
    int.anchors <- FindIntegrationAnchors(object.list=seurat_objects, normalization.method="SCT", anchor.features=features)
    seurat_objects_integrated <- IntegrateData(anchorset=int.anchors, normalization.method="SCT")
    rm(int.anchors)
    gc()
    
    ## dimensionality reduction
    seurat_objects_integrated <- RunPCA(seurat_objects_integrated, verbose=FALSE)
    reduction2use <- "pca"
    for (dimensions_nb_param in dimensions_nb_param_vector) {
        print(sprintf("dimension number: %d", dimensions_nb_param))
        dimensions_output_dir <- sprintf("%s/%d_PCs", integration_output_dir, dimensions_nb_param)
        if (! dir.exists(dimensions_output_dir)) {
            dir.create(dimensions_output_dir, recursive=TRUE, mode="0775")
        }
        dimensions_output_name <- sprintf("%s_%dPCs", integration_output_name, dimensions_nb_param)
                
        seurat_objects_integrated <- FindNeighbors(seurat_objects_integrated, reduction=reduction2use, dims=1:dimensions_nb_param)
        ### compute silhouette scores for lames, zones, conditions and samples
        mean_silhouette_vector <- c()
        for (metadata_category_colname in metadata_categories) {
            seurat_objects_integrated <- compute_silhouette(seurat_objects_integrated, reduction2use, dimensions_nb_param, metadata_category_colname)
            mean_silhouette <- mean(seurat_objects_integrated@meta.data[[sprintf("silhouette_%s", metadata_category_colname)]])
            mean_silhouette_vector <- c(mean_silhouette_vector, mean_silhouette)
        }
        mean_silhouette_df <- as.data.frame(t(mean_silhouette_vector))
        colnames(mean_silhouette_df) <- sprintf("mean_silhouette_%s", metadata_categories)
        mean_silhouette_df <- cbind(normalization="SCTransform", features=features_param, integration=integration_method, dimensions=dimensions_nb_param, mean_silhouette_df)
        integration_batch_mean_silhouette_df <- rbind(integration_batch_mean_silhouette_df, mean_silhouette_df)
        
        ## clustering
        for (clustering_resolution in clustering_resolution_vector) {
            print(sprintf("cluster resolution: %.1f", clustering_resolution))
            resolution_output_dir <- sprintf("%s/resolution%s", dimensions_output_dir, sub("[.]", "_", clustering_resolution))
            if (! dir.exists(resolution_output_dir)) {
                dir.create(resolution_output_dir, recursive=TRUE, mode="0775")
            }
            resolution_output_name <- sprintf("%s_resolution%s", dimensions_output_name, sub("[.]", "_", clustering_resolution))

            ### set default assay to integrated before clustering
            if (DefaultAssay(seurat_objects_integrated) != "integrated") {
                DefaultAssay(seurat_objects_integrated) <- "integrated"
            }
            
            seurat_objects_integrated <- clustering(seurat_objects_integrated, clustering_resolution, reduction2use, dimensions_nb_param)
            write.csv(seurat_objects_integrated@meta.data, file=sprintf("%s/%s_metadata.csv", resolution_output_dir, resolution_output_name), quote=FALSE, row.names=TRUE)
            cluster_mean_silhouette <- mean(seurat_objects_integrated$silhouette_seurat_clusters)
            ### plot clusters onto UMAP or onto the tissue section
            clustering_plots(seurat_objects_integrated, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integrated datasets: %s, all spots", integration_method), resolution_output_dir, sprintf("%s_clustering", resolution_output_name))
            
            ### hippocampus spots
            seurat_objects_integrated_hippocampus <- subset(seurat_objects_integrated, cells=all_sections_hippocampus_barcodes_Seurat_integrated)
            hippocampus_cluster_mean_silhouette <- mean(seurat_objects_integrated_hippocampus$silhouette_seurat_clusters)
            #### plot clusters onto UMAP or onto the tissue section
            clustering_plots(seurat_objects_integrated_hippocampus, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integrated datasets: %s, hippocampus spots", integration_method), resolution_output_dir, sprintf("%s_clustering_hippocampus", resolution_output_name))
            
            ### get mean silhouette scores: all spots and hippocampus spots
            mean_silhouette_df <- data.frame(normalization="SCTransform", features=features_param, integration=integration_method, dimensions=dimensions_nb_param, resolution=clustering_resolution, mean_silhouette_seurat_clusters=cluster_mean_silhouette, mean_silhouette_seurat_clusters_hippocampus=hippocampus_cluster_mean_silhouette)
            integration_cluster_mean_silhouette_df <- rbind(integration_cluster_mean_silhouette_df, mean_silhouette_df)
            
            ### cluster markers
            if (cluster_markers) {
                # find markers for every cluster compared to all remaining spots, report only the positive ones
                markers_output_dir <- sprintf("%s/00-Cluster_markers", resolution_output_dir)
                if (! dir.exists(markers_output_dir)) {
                    dir.create(markers_output_dir, recursive=TRUE, mode="0775")
                }
                identify_cluster_markers(seurat_objects_integrated, "SCT", load_image, nb_cluster_markers_to_plot, nb_SpatialFeaturePlot_per_page, markers_output_dir, sprintf("%s_cluster_marker_genes", resolution_output_name))
            }

            ### differentially expressed genes across conditions
            if (DE_analysis) {
                DE_output_dir <- sprintf("%s/10-DEGs", resolution_output_dir)
                if (! dir.exists(DE_output_dir)) {
                    dir.create(DE_output_dir, recursive=TRUE, mode="0775")
                }
                identify_DEGs(seurat_objects_integrated, "SCT", "seurat_clusters", "condition", DE_output_dir, sprintf("%s_DEGs", resolution_output_name))
            }

            ### conserved cell type markers
            if (conserved_markers) {
                conserved_markers_output_dir <- sprintf("%s/20-Conserved_markers", resolution_output_dir)
                if (! dir.exists(conserved_markers_output_dir)) {
                    dir.create(conserved_markers_output_dir, recursive=TRUE, mode="0775")
                }
                identify_conserved_markers(seurat_objects_integrated, "SCT", "seurat_clusters", "condition", conserved_markers_output_dir, sprintf("%s_conserved_markers", resolution_output_name))
            }
        }
    }
    rm(seurat_objects_integrated)
    rm(seurat_objects_integrated_hippocampus)
    gc()
}

## Harmony
integration_method <- "Harmony"
### Harmony integration function requires a merged Seurat object
seurat_objects_merge <- merge(seurat_objects[[1]], y=unlist(seurat_objects)[2:nb_sections], add.cell.ids=sections)
#### orig.ident, lame, zone and condition are not factors in the merged objeas.factor(spatial_coordinates_df$zone)
seurat_objects_merge@meta.data$orig.ident <- as.factor(seurat_objects_merge@meta.data$orig.ident)
seurat_objects_merge@meta.data$lame <- as.factor(seurat_objects_merge@meta.data$lame)
seurat_objects_merge@meta.data$zone <- as.factor(seurat_objects_merge@meta.data$zone)
seurat_objects_merge@meta.data$condition <- as.factor(seurat_objects_merge@meta.data$condition)
rm(seurat_objects)
gc()

for (features_param in features_param_vector) {
    ### normalization
    if (features_param > 10) {
        #### features_param is a number of features
        print(sprintf("features param: %d", features_param))
        normalization_output_dir <- sprintf("%s/00-SCTransform/featurenumber%d", output_dir, features_param)
        if (! dir.exists(normalization_output_dir)) {
            dir.create(normalization_output_dir, recursive=TRUE, mode="0775")
        }
        normalization_output_name <- sprintf("SCTransform_featurenumber%d", features_param)
        seurat_objects_merge <- SCTransform(seurat_objects_merge, variable.features.n=features_param, return.only.var.genes=FALSE)
    } else {
        #### features_param is a residual variance threshold
        print(sprintf("features param: %.1f", features_param))
        normalization_output_dir <- sprintf("%s/00-SCTransform/residualvariance%s", output_dir, sub("[.]", "_", features_param))
        if (! dir.exists(normalization_output_dir)) {
            dir.create(normalization_output_dir, recursive=TRUE, mode="0775")
        }
        normalization_output_name <- sprintf("SCTransform_residualvariance%s", sub("[.]", "_", features_param))
        seurat_objects_merge <- SCTransform(seurat_objects_merge, variable.features.n=NULL, variable.features.rv.th=features_param, return.only.var.genes=FALSE)
        ## TODO: ## get minimum number of variable features
    }

    ### compute PCs 
    seurat_objects_merge <- RunPCA(seurat_objects_merge, verbose=FALSE)
    for (dimensions_nb_param in dimensions_nb_param_vector) {
        print(sprintf("dimension number: %d", dimensions_nb_param))
        dimensions_output_dir <- sprintf("%s/11-Harmony/%d_PCs", normalization_output_dir, dimensions_nb_param)
        if (! dir.exists(dimensions_output_dir)) {
            dir.create(dimensions_output_dir, recursive=TRUE, mode="0775")
        }
        dimensions_output_name <- sprintf("%s_%s_%dPCs", normalization_output_name, integration_method, dimensions_nb_param)

        ### integration and dimensionality reduction
        seurat_objects_integrated <- RunHarmony(seurat_objects_merge, group.by.vars="orig.ident", dims=1:dimensions_nb_param)
        reduction2use <- "harmony"
        seurat_objects_integrated <- FindNeighbors(seurat_objects_integrated, reduction=reduction2use, dims=1:dimensions_nb_param)
        ### compute silhouette scores for lames, zones, conditions and samples
        mean_silhouette_vector <- c()
        for (metadata_category_colname in metadata_categories) {
            seurat_objects_integrated <- compute_silhouette(seurat_objects_integrated, reduction2use, dimensions_nb_param, metadata_category_colname)
            mean_silhouette <- mean(seurat_objects_integrated@meta.data[[sprintf("silhouette_%s", metadata_category_colname)]])
            mean_silhouette_vector <- c(mean_silhouette_vector, mean_silhouette)
        }
        mean_silhouette_df <- as.data.frame(t(mean_silhouette_vector))
        colnames(mean_silhouette_df) <- sprintf("mean_silhouette_%s", metadata_categories)
        mean_silhouette_df <- cbind(normalization="SCTransform", features=features_param, integration=integration_method_covariate, dimensions=dimensions_nb_param, mean_silhouette_df)
        integration_batch_mean_silhouette_df <- rbind(integration_batch_mean_silhouette_df, mean_silhouette_df)

        ### clustering
        for (clustering_resolution in clustering_resolution_vector) {
            print(sprintf("cluster resolution: %.1f", clustering_resolution))
            resolution_output_dir <- sprintf("%s/resolution%s", dimensions_output_dir, sub("[.]", "_", clustering_resolution))
            if (! dir.exists(resolution_output_dir)) {
                dir.create(resolution_output_dir, recursive=TRUE, mode="0775")
            }
            resolution_output_name <- sprintf("%s_resolution%s", dimensions_output_name, sub("[.]", "_", clustering_resolution))

            seurat_objects_integrated <- clustering(seurat_objects_integrated, clustering_resolution, reduction2use, dimensions_nb_param)
            write.csv(seurat_objects_integrated@meta.data, file=sprintf("%s/%s_metadata.csv", resolution_output_dir, resolution_output_name), quote=FALSE, row.names=TRUE)
            cluster_mean_silhouette <- mean(seurat_objects_integrated$silhouette_seurat_clusters)
            ### plot clusters onto UMAP or onto the tissue section
            clustering_plots(seurat_objects_integrated, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integrated datasets: %s, all spots", integration_method_covariate), resolution_output_dir, sprintf("%s_clustering", resolution_output_name))

            ### hippocampus spots
            seurat_objects_integrated_hippocampus <- subset(seurat_objects_integrated, cells=all_sections_hippocampus_barcodes_Seurat_merge)
            hippocampus_cluster_mean_silhouette <- mean(seurat_objects_integrated_hippocampus$silhouette_seurat_clusters)
            #### plot clusters onto UMAP or onto the tissue section
            clustering_plots(seurat_objects_integrated_hippocampus, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integrated datasets: %s, hippocampus spots", integration_method_covariate), resolution_output_dir, sprintf("%s_clustering_hippocampus", resolution_output_name))

            ### get mean silhouette scores: all spots and hippocampus spots
            mean_silhouette_df <- data.frame(normalization="SCTransform", features=features_param, integration=integration_method_covariate, dimensions=dimensions_nb_param, resolution=clustering_resolution, mean_silhouette_seurat_clusters=cluster_mean_silhouette, mean_silhouette_seurat_clusters_hippocampus=hippocampus_cluster_mean_silhouette)
            integration_cluster_mean_silhouette_df <- rbind(integration_cluster_mean_silhouette_df, mean_silhouette_df)

            ### cluster markers
            if (cluster_markers) {
                # find markers for every cluster compared to all remaining spots, report only the positive ones
                markers_output_dir <- sprintf("%s/00-Cluster_markers", resolution_output_dir)
                if (! dir.exists(markers_output_dir)) {
                    dir.create(markers_output_dir, recursive=TRUE, mode="0775")
                }
                identify_cluster_markers(seurat_objects_integrated, "SCT", load_image, nb_cluster_markers_to_plot, nb_SpatialFeaturePlot_per_page, markers_output_dir, sprintf("%s_cluster_marker_genes", resolution_output_name))
            }

            ### differentially expressed genes across conditions
            if (DE_analysis) {
                DE_output_dir <- sprintf("%s/10-DEGs", resolution_output_dir)
                if (! dir.exists(DE_output_dir)) {
                    dir.create(DE_output_dir, recursive=TRUE, mode="0775")
                }
                identify_DEGs(seurat_objects_integrated, "SCT", "seurat_clusters", "condition", DE_output_dir, sprintf("%s_DEGs", resolution_output_name))
            }

            ### conserved cell type markers
            if (conserved_markers) {
                conserved_markers_output_dir <- sprintf("%s/20-Conserved_markers", resolution_output_dir)
                if (! dir.exists(conserved_markers_output_dir)) {
                    dir.create(conserved_markers_output_dir, recursive=TRUE, mode="0775")
                }
                identify_conserved_markers(seurat_objects_integrated, "SCT", "seurat_clusters", "condition", conserved_markers_output_dir, sprintf("%s_conserved_markers", resolution_output_name))
            }
        }
        rm(seurat_objects_integrated)
        gc()
    }
}
rm(seurat_objects_merge)
gc()
write.csv(integration_cluster_mean_silhouette_df, file=sprintf("%s/00-SCTransform/SCTransform_featuresnumber_integration_PCs_resolution_cluster_mean_silhouette.csv", output_dir), quote=FALSE, row.names=FALSE)
integration_cluster_mean_silhouette_df$integration <- factor(integration_cluster_mean_silhouette_df$integration, levels=c("no integration", "Seurat", "Harmony"))
integration_cluster_mean_silhouette_df$features <- factor(integration_cluster_mean_silhouette_df$features, levels=features_param_vector)
integration_cluster_mean_silhouette_df$dimensions <- as.factor(integration_cluster_mean_silhouette_df$dimensions)
integration_cluster_mean_silhouette_df$resolution <- as.factor(integration_cluster_mean_silhouette_df$resolution)
pdf(sprintf("%s/00-SCTransform/SCTransform_featuresnumber_integration_PCs_resolution_cluster_mean_silhouette.pdf", output_dir))
for (one_category in c("seurat_clusters", "seurat_clusters_hippocampus")) {
    p <- ggplot(integration_cluster_mean_silhouette_df, aes(x=resolution, y=.data[[sprintf("mean_silhouette_%s", one_category)]], fill=integration)) +
    #p <- ggplot(integration_cluster_mean_silhouette_df, aes(x=resolution, y=!!sym(sprintf("mean_silhouette_%s", one_category)), fill=integration)) +
        geom_bar(stat="identity", position="dodge") +
        facet_grid(dimensions~features) +
        labs(title=sprintf("Mean silhouette: %s", one_category), x="Resolution", y="Mean silhouette", fill="Integration method") +
        theme_bw() +
        theme(legend.position="bottom") +
        theme(axis.text.x=element_text(angle=60, hjust=1)) +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
}
dev.off()

write.csv(integration_batch_mean_silhouette_df, file=sprintf("%s/00-SCTransform/SCTransform_featuresnumber_integration_PCs_resolution_batch_mean_silhouette.csv", output_dir), quote=FALSE, row.names=FALSE)
integration_batch_mean_silhouette_df$integration <- factor(integration_batch_mean_silhouette_df$integration, levels=c("no integration", "Seurat", "Harmony"))
integration_batch_mean_silhouette_df$features <- factor(integration_batch_mean_silhouette_df$features, levels=features_param_vector)
integration_batch_mean_silhouette_df$dimensions <- as.factor(integration_batch_mean_silhouette_df$dimensions)
pdf(sprintf("%s/00-SCTransform/SCTransform_featuresnumber_integration_PCs_resolution_batch_mean_silhouette.pdf", output_dir))
for (one_category in c("lame", "zone", "condition", "orig.ident")) {
    p <- ggplot(integration_batch_mean_silhouette_df, aes(x=dimensions, y=.data[[sprintf("mean_silhouette_%s", one_category)]], fill=integration)) +
        geom_bar(stat="identity", position="dodge") +
        facet_grid(~features) +
        labs(title=sprintf("Mean silhouette: %s", one_category), x="Number of dimensions", y="Mean silhouette", fill="Integration method") +
        theme_bw() +
        theme(legend.position="bottom") +
        theme(axis.text.x=element_text(angle=60, hjust=1)) +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
}

for (one_category in c("lame", "zone", "condition", "orig.ident")) {
    p <- ggplot(integration_batch_mean_silhouette_df, aes(x=integration, y=.data[[sprintf("mean_silhouette_%s", one_category)]], fill=integration)) +
        geom_bar(stat="identity", position="dodge") +
        facet_grid(dimensions~features) +
        labs(title=sprintf("Mean silhouette: %s", one_category), x="Integration method", y="Mean silhouette", fill="Integration method") +
        theme_bw() +
        theme(legend.position="bottom") +
        theme(axis.text.x=element_text(angle=60, hjust=1)) +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
}
dev.off()











for (one_integration_method in levels(integration_cluster_mean_silhouette_df$integration)) {
    p <- ggplot(integration_cluster_mean_silhouette_df[which(integration_cluster_mean_silhouette_df$integration==one_integration_method),], aes(x=resolution, y=mean_silhouette_seurat_clusters, fill=resolution)) +
        geom_bar(stat="identity", position="dodge") +
        facet_grid(dimensions~features) +
        labs(title=sprintf("Cluster mean silhouette: %s", one_integration_method), x="Number of features", y="Number of dimensions", fill="Resolution") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=60, hjust=1)) +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
}


for (one_feature_number in levels(integration_cluster_mean_silhouette_df$features)) {
    p <- ggplot(integration_cluster_mean_silhouette_df[which(integration_cluster_mean_silhouette_df$features==one_feature_number),], aes(x=resolution, y=mean_silhouette_seurat_clusters, fill=resolution)) +
        geom_bar(stat="identity", position="dodge") +
        facet_grid(integration~dimensions) +
        labs(title=sprintf("Cluster mean silhouette: %s", one_feature_number), x="Number of features", y="Number of dimensions", fill="Resolution") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=60, hjust=1)) +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
}





p <- ggplot(integration_cluster_mean_silhouette_df, aes(x=features_param, y=silhouette_clusters, fill=features_param)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(integration~resolution) +
    labs(title="Cluster mean silhouette", x="Resolution", y="Integration method", fill="Number of\nfeatures") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=60, hjust=1)) +
    theme(panel.border=element_rect(color="grey50"))
print(p)
dev.off()








#one_feature <- "Cd74"
#SpatialFeaturePlot(seurat_objects_integrated, features=one_feature, images=names(seurat_objects_integrated@images)[c(3,7)], alpha=c(0.1,1), ncol=2) + plot_annotation(title=sprintf("Cluster: %s", one_cluster), theme=theme(plot.title=element_text(size=16)))
#SpatialFeaturePlot(seurat_objects_integrated, features=one_feature, images=names(seurat_objects_integrated@images)[c(4,8)], alpha=c(0.1,1), ncol=2) + plot_annotation(title=sprintf("Cluster: %s", one_cluster), theme=theme(plot.title=element_text(size=16)))





print(SpatialFeaturePlot(seurat_objects_integrated, features="Serpina11", images=names(seurat_objects_integrated@images), alpha=c(0.1,1), ncol=2) + plot_annotation(title=sprintf("Cluster: %s", one_cluster), theme=theme(plot.title=element_text(size=16))))


pct_threshold <- 0.7
seurat_objects_integrated.DE.all_clusters.ordered.pct <- seurat_objects_integrated.DE.all_clusters.ordered[which(seurat_objects_integrated.DE.all_clusters.ordered$pct.1 > pct_threshold & seurat_objects_integrated.DE.all_clusters.ordered$pct.2 > pct_threshold),]


one_feature <- "Tmem106a"
SpatialFeaturePlot(seurat_objects_integrated, features=one_feature, images=names(seurat_objects_integrated@images)[c(1,2)], alpha=c(0.1,1), ncol=2) + plot_annotation(title=sprintf("Cluster: %s", one_cluster), theme=theme(plot.title=element_text(size=16)))
SpatialFeaturePlot(seurat_objects_integrated, features=one_feature, images=names(seurat_objects_integrated@images)[c(3,4)], alpha=c(0.1,1), ncol=2) + plot_annotation(title=sprintf("Cluster: %s", one_cluster), theme=theme(plot.title=element_text(size=16)))








# identification of spatially variable genes
## Finding differentially expressed features (cluster biomarkers): https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
print("cluster marker genes")
### find markers for every cluster compared to all remaining spots, report only the positive ones
markers <- FindAllMarkers(seurat_objects_integrated, logfc.threshold=0.25, min.pct=0.25, only.pos=TRUE)
write.csv(markers, file=sprintf("%s/30-Cluster_marker_genes.csv", output_dir), quote=FALSE, row.names=TRUE)
pdf(sprintf("%s/30-Cluster_marker_genes.pdf", output_dir))
#### expression heatmap
top10_markers <- markers %>%
    group_by(cluster) %>%
    top_n(n=10, wt=avg_log2FC)
DoHeatmap(seurat_objects_integrated, features=top10_markers$gene) + NoLegend()
#### get top markers for each cluster
if (load_image) {
    top_markers <- markers %>%
        group_by(cluster) %>%
        slice_max(n=nb_cluster_markers_to_plot, order_by=avg_log2FC)
    for (one_cluster in levels(markers$cluster)) {
        one_cluster_markers <- top_markers[which(top_markers$cluster==one_cluster), "gene", drop=TRUE]
        for (i in 1:(nb_cluster_markers_to_plot/2)) {
            for (j in 1:(nb_sections/nb_SpatialFeaturePlot_per_page)) {
                print(SpatialFeaturePlot(object=seurat_objects_integrated, features=one_cluster_markers[((i-1)*2+1):(i*2)], images=names(seurat_objects_integrated@images)[((j-1)*nb_SpatialFeaturePlot_per_page+1):(j*nb_SpatialFeaturePlot_per_page)], alpha=c(0.1,1), ncol=3) + plot_annotation(title=sprintf("Cluster: %s", one_cluster), theme=theme(plot.title=element_text(size=16))))
            }
        }
    }
}
dev.off()

## cluster marker genes
#print("cluster marker genes")
## differential expression between clusters 10 and 17 vs all others respectively
#write.csv(markers, file=sprintf("%s/30-Cluster_marker_genes.csv", output_dir), quote=FALSE, row.names=TRUE)
#pdf(sprintf("%s/30-Cluster_marker_genes.pdf", output_dir))
#for (one_cluster in c(10, 17)) {
#    print(sprintf("cluster: %d"))
#    print("find markers")
#    de_markers <- FindMarkers(merge_seurat_object.integrated, ident.1=one_cluster, ident.2=NULL)
#    print("plots")
#    for (i in 1:(nb_cluster_markers_to_plot/2)) {
#        for (j in 1:(nb_sections/nb_SpatialFeaturePlot_per_page)) {
#            print(SpatialFeaturePlot(object=merge_seurat_object.integrated, features=rownames(de_markers)[((i-1)*2+1):(i*2)], images=names(merge_seurat_object@images)[((j-1)*nb_SpatialFeaturePlot_per_page+1):(j*nb_SpatialFeaturePlot_per_page)], alpha=c(0.1,1), ncol=3) + plot_annotation(title=sprintf("Cluster: %d", one_cluster), theme=theme(plot.title=element_text(size=16))))
#        }
#    }
#}
#dev.off()

## genes with spatial patterning
### find spatially variable genes
seurat_objects_integrated <- FindSpatiallyVariableFeatures(seurat_objects_integrated, assay="integrated", features=VariableFeatures(seurat_objects_integrated)[1:1000], selection.method="markvariogram")
### visualize the expression of the top6 features
pdf(sprintf("%s/31-Spatially_variable_genes.pdf", output_dir))
top.features <- head(SpatiallyVariableFeatures(seurat_objects_integrated, selection.method="markvariogram"), 3)
for (i in 1:(nb_sections/nb_SpatialFeaturePlot_per_page)) {
    print(SpatialFeaturePlot(seurat_objects_integrated, features=top.features, images=names(seurat_objects_integrated@images)[((i-1)*nb_SpatialFeaturePlot_per_page+1):(i*nb_SpatialFeaturePlot_per_page)], alpha=c(0.1,1), ncol=3))
}
dev.off()




