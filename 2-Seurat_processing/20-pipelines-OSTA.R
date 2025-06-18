.libPaths(c("/home/christophe.lepriol/NeuroDev_ADD/R/r_4.1.0", .libPaths()))
library(SpatialExperiment)
library(ggspavis)
library(scater)
library(patchwork) # for wrap_plots() function
library(scran) # for buildSNNGraph() function
library(RColorBrewer)
library(pheatmap)


##############
# Parameters #
##############
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"
genome_name <- "NCBIRefSeq108_NCBIRefSeq108GTF"
samples <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")

color_palette <- c(brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 8, name = "Dark2"), brewer.pal(n = 8, name = "Accent"))


#############
# Functions #
#############
modif_point_size_rm_legend_title_modif_title <- function(one_plot, point_size, one_title) {
    if (! is.na(point_size)) {
        one_plot$layers[[1]]$aes_params$size <- point_size
    }
    one_plot <- one_plot + labs(title=one_title) + theme(legend.title=element_blank())
}

modif_point_size_rm_title <- function(one_plot, point_size) {
    one_plot$layers[[1]]$aes_params$size <- point_size
    one_plot <- one_plot + labs(title="")
}


############
# Analysis #
############

dataset_plot_list <- list()
for (sample in samples) {
    print(sprintf("sample: %s", sample))
    output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s/20-pipelines/00-OSTA", work_dir, genome_name, sample)
    if (! dir.exists(output_dir)) {
        dir.create(output_dir, recursive=TRUE, mode="0775")
    }
    sample_plot_list <- list()
    
    #############
    # Load data #
    #############

    # expression matrices
    ## Space Ranger
    #h5_file <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs/raw_feature_bc_matrix.h5", work_dir, genome_name, sample)
    #h5_data <- Read10X_h5(h5_file)
    ## ST Pipeline
    st_pipeline_matrix_file <- sprintf("%s/10-ST_analysis/00-ST_Pipeline/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s/10-Pipeline/%s_stdata.tsv", work_dir, genome_name, sample, sample)
    st_pipeline_matrix <- read.table(st_pipeline_matrix_file, sep="\t", header=TRUE, quote="", row.names=1)
    # Space Ranger output directory for tissue positions list and image
    space_ranger_output_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs", work_dir, genome_name, sample)

    # build SpatialExperiment object
    spatial_coordinates_file <- file.path(space_ranger_output_dir, "spatial", "tissue_positions_list.csv")
    spatial_coordinates_df <- read.csv(spatial_coordinates_file, header=FALSE, quote="")
    colnames(spatial_coordinates_df) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    rownames(spatial_coordinates_df) <- spatial_coordinates_df$barcode
    ## matrix: convert coordinates to barcodes
    barcodes <- unlist(lapply(rownames(st_pipeline_matrix), function(x, coord2barcode=spatial_coordinates_df) {
        coordinates <- unlist(strsplit(x, "x"))
        return(coord2barcode[which(coord2barcode$array_col==as.integer(coordinates[1])-1 & coord2barcode$array_row==as.integer(coordinates[2])-1), "barcode"])
    }))
    rownames(st_pipeline_matrix) <- barcodes
    gene_data <- data.frame(gene_name=colnames(st_pipeline_matrix))
    ## image data
    img <- readImgData(imageSources=file.path(space_ranger_output_dir, "spatial", "tissue_lowres_image.png"), scaleFactors=file.path(space_ranger_output_dir, "spatial", "scalefactors_json.json"), sample_id=sample)

    spe <- SpatialExperiment(assay=list(counts=t(st_pipeline_matrix)), rowData=gene_data, colData=spatial_coordinates_df[rownames(st_pipeline_matrix),], imgData=img, spatialDataNames=c("barcode", "in_tissue", "array_row", "array_col"), spatialCoordsNames=c("pxl_col_in_fullres", "pxl_row_in_fullres"), sample_id=sample)

    # keep spot over tissue
    spe <- spe[, colData(spe)$in_tissue==1]

    ######
    # QC #
    ######

    # identify mitochondrial genes
    is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
    table(is_mito)

    # calculate per-spot QC metrics and store in colData
    assay_type <- "counts"
    spe <- addPerCellQC(spe, assay.type=assay_type, subsets=list(mito=is_mito))
    ## rename columns
    for (one_stat in c("sum", "detected", "total")) {
        colnames(colData(spe))[colnames(colData(spe)) == one_stat] <- sprintf("%s.%s", assay_type, one_stat)
    }
    head(colData(spe))

    # plot sum of counts and number of detected genes per spot
    pdf(sprintf("%s/00-QC.pdf", output_dir))
    ## library size
    hist(colData(spe)[, sprintf("%s.sum", assay_type)], breaks=20, main="Library size per spot", xlab="Count")
    plot_spots <- plotSpots(spe, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate=sprintf("%s.sum", assay_type), palette=c("white", "black"), size=2)
    plot_spots <- plot_spots + labs(title="Library size per spot") + theme(legend.title=element_blank())
    print(plot_spots)
    plot_visium <- plotVisium(spe, fill=sprintf("%s.sum", assay_type), x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    plot_visium <- plot_visium + labs(title="Library size per spot") + theme(legend.title=element_blank())
    print(plot_visium)
    ## number of expressed genes
    hist(colData(spe)[, sprintf("%s.sum", assay_type)], breaks=20, main="Number of expressed genes per spot", xlab="Count")
    plot_spots <- plotSpots(spe, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate=sprintf("%s.detected", assay_type), palette=c("white", "black"), size=2)
    plot_spots <- plot_spots + labs(title="Number of expressed genes per spot") + theme(legend.title=element_blank())
    print(plot_spots)
    plot_visium <- plotVisium(spe, fill=sprintf("%s.detected", assay_type), x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    plot_visium <- plot_visium + labs(title="Number of expressed genes per spot") + theme(legend.title=element_blank())
    print(plot_visium)

    # thresholds
    ## library size
    threshold_lib_size <- 2000
    qc_lib_size <- colData(spe)[, sprintf("%s.sum", assay_type)] < threshold_lib_size
    table(qc_lib_size)
    colData(spe)$qc_lib_size <- qc_lib_size
    colnames(colData(spe))[which(colnames(colData(spe)) == "qc_lib_size")] <- sprintf("qc_lib_size_%d", threshold_lib_size)
    plot_qc <- plotQC(spe, type="spots", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", discard=sprintf("qc_lib_size_%d", threshold_lib_size))
    plot_qc <- modif_point_size_rm_legend_title_modif_title(plot_qc, 2, sprintf("QC spots: library size < %d", threshold_lib_size))
    print(plot_qc)
    ## number of expressed features
    threshold_exp_genes <- 1000
    qc_detected <- colData(spe)[, sprintf("%s.detected", assay_type)] < threshold_exp_genes
    table(qc_detected)
    colData(spe)$qc_detected <- qc_detected
    colnames(colData(spe))[which(colnames(colData(spe)) == "qc_detected")] <- sprintf("qc_detected_%d", threshold_exp_genes)
    plot_qc <- plotQC(spe, type="spots", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", discard=sprintf("qc_detected_%d", threshold_exp_genes))
    plot_qc <- modif_point_size_rm_legend_title_modif_title(plot_qc, 2, sprintf("QC spots: number of expressed genes < %d", threshold_exp_genes))
    print(plot_qc)

    # remove low-quality spots
    ## number of discarded spots for each metric
    apply(cbind(colData(spe)$qc_lib_size_2000, colData(spe)$qc_detected_1000), 2, sum)
    ## combined set of discarded spots
    discard <- colData(spe)$qc_lib_size_2000 | colData(spe)$qc_detected_1000
    table(discard)
    ## store in object
    colData(spe)$discard <- discard
    plot_qc <- plotQC(spe, type="scatter", metric_x=sprintf("%s.sum", assay_type), metric_y=sprintf("%s.detected", assay_type), threshold_x=threshold_lib_size, threshold_y=threshold_exp_genes)
    plot_qc <- plot_qc + labs(x="Library size", y="Number of expressed genes")
    print(plot_qc)
    plot_qc <- plotQC(spe, type="spots", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", discard="discard")
    sample_plot_list[["qc_discard"]] <- plot_qc
    plot_qc <- modif_point_size_rm_legend_title_modif_title(plot_qc, 2, "QC spots: discarded spots")
    print(plot_qc)
    ## remove combined set of low-quality spots
    spe <- spe[, !colData(spe)$discard]
    dim(spe)
    plot_spots <- plotSpots(spe, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate=sprintf("%s.sum", assay_type), palette=c("white", "black"), size=2)
    plot_spots <- plot_spots + labs(title="Library size per spot after low-quality spot filtering") + theme(legend.title=element_blank())
    print(plot_spots)
    plot_spots <- plotSpots(spe, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate=sprintf("%s.detected", assay_type), palette=c("white", "black"), size=2)
    plot_spots <- plot_spots + labs(title="Number of expressed genes per spot after low-quality spot filtering") + theme(legend.title=element_blank())
    print(plot_spots)
    dev.off()
    write.csv(colData(spe), file=sprintf("%s/00-QC.csv", output_dir), quote=FALSE, row.names=TRUE)

    #################
    # Normalization #
    #################
    # quick clustering for pool-based size factors
    set.seed(123)
    qclus <- quickCluster(spe)
    table(qclus)
    # calculate size factors and store in object
    spe <- computeSumFactors(spe, cluster=qclus)
    summary(sizeFactors(spe))
    # calculate logcounts (log-transformed normalized counts) and store in object
    spe <- logNormCounts(spe)
    # check
    assayNames(spe)
    dim(counts(spe))
    dim(logcounts(spe))

    # calculate per-spot QC metrics and store in colData
    ## normalized counts
    assay_type <- "logcounts"
    spe <- addPerCellQC(spe, assay.type=assay_type)
    for (one_stat in c("sum", "detected", "total")) {
        colnames(colData(spe))[colnames(colData(spe)) == one_stat] <- sprintf("%s.%s", assay_type, one_stat)
    }
    head(colData(spe))
    write.csv(colData(spe), file=sprintf("%s/10-Normalization.csv", output_dir), quote=FALSE, row.names=TRUE)

    # plots
    pdf(sprintf("%s/10-Normalization.pdf", output_dir))
    hist(sizeFactors(spe), breaks=20, main="Normalization size factors", xlab="Value")
    ## total counts and number of detected genes
    hist(colData(spe)$counts.sum, breaks=20, main="Library size per spot", xlab="Count")
    hist(colData(spe)$logcounts.sum, breaks=20, main="Library size per spot after normalization", xlab="Count")
    hist(colData(spe)$counts.detected, breaks=20, main="Number of expressed genes per spot", xlab="Count")
    hist(colData(spe)$logcounts.detected, breaks=20, main="Number of expressed genes per spot after normalization", xlab="Count")
    ## splot plots
    ### raw counts
    raw_lib_size_plot_spots <- plotSpots(spe, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="counts.sum", palette=c("white", "black"), size=2)
    raw_lib_size_plot_visium <- plotVisium(spe, fill="counts.sum", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    raw_exp_genes_plot_spots <- plotSpots(spe, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="counts.detected", palette=c("white", "black"), size=2)
    raw_exp_genes_plot_visium <- plotVisium(spe, fill="counts.detected", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    ### normalized counts
    norm_lib_size_plot_spots <- plotSpots(spe, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="logcounts.sum", palette=c("white", "black"), size=2)
    sample_plot_list[["norm_lib_size"]] <- norm_lib_size_plot_spots
    norm_lib_size_plot_visium <- plotVisium(spe, fill="logcounts.sum", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    norm_exp_genes_plot_spots <- plotSpots(spe, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="logcounts.detected", palette=c("white", "black"), size=2)
    sample_plot_list[["norm_exp_genes"]] <- norm_exp_genes_plot_spots
    norm_exp_genes_plot_visium <- plotVisium(spe, fill="logcounts.detected", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    ### raw and normalized counts
    #### library size
    raw_lib_size_plot_spots <- modif_point_size_rm_legend_title_modif_title(raw_lib_size_plot_spots, 0.5, "Raw counts")
    raw_lib_size_plot_visium <- modif_point_size_rm_legend_title_modif_title(raw_lib_size_plot_visium, NA, "Raw counts")
    norm_lib_size_plot_spots <- modif_point_size_rm_legend_title_modif_title(norm_lib_size_plot_spots, 0.5, "Log counts")
    norm_lib_size_plot_visium <- modif_point_size_rm_legend_title_modif_title(norm_lib_size_plot_visium, NA, "Log counts")
    print(wrap_plots(raw_lib_size_plot_spots, norm_lib_size_plot_spots, raw_lib_size_plot_visium, norm_lib_size_plot_visium) + plot_annotation(title="Library size per spot"))
    #### number of expressed genes
    raw_exp_genes_plot_spots <- modif_point_size_rm_legend_title_modif_title(raw_exp_genes_plot_spots, 0.5, "Raw counts")
    raw_exp_genes_plot_visium <- modif_point_size_rm_legend_title_modif_title(raw_exp_genes_plot_visium, NA, "Raw counts")
    norm_exp_genes_plot_spots <- modif_point_size_rm_legend_title_modif_title(norm_exp_genes_plot_spots, 0.5, "Log counts")
    norm_exp_genes_plot_visium <- modif_point_size_rm_legend_title_modif_title(norm_exp_genes_plot_visium, NA, "Log counts")
    print(wrap_plots(raw_exp_genes_plot_spots, norm_exp_genes_plot_spots, raw_exp_genes_plot_visium, norm_exp_genes_plot_visium) + plot_annotation(title="Number of expressed genes per spot"))
    dev.off()

    ##############
    # Clustering #
    ##############

    # feature selection
    ## HVGs: highly variable genes
    ### fit mean-variance relationship
    dec <- modelGeneVar(spe)
    ### visualize mean-variance relationship
    fit <- metadata(dec)
    pdf(sprintf("%s/20-Clustering.pdf", output_dir))
    plot(fit$mean, fit$var, main="Fit mean-variance relationship", xlab="Mean of log-expression", ylab="Variance of log-expression")
    curve(fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)
    ## select top HVGs: top 10%
    top_hvgs <- getTopHVGs(dec, prop=0.1)
    length(top_hvgs)

    # dimensionality reduction
    ## apply of PCA to HVGs and retain the top 50 PCs
    set.seed(123)
    spe <- runPCA(spe, subset_row=top_hvgs)
    reducedDimNames(spe)
    dim(reducedDim(spe, "PCA"))

    ## run UMAP on the set of top 50 PCs and retain the top 2 UMAP components
    set.seed(123)
    spe <- runUMAP(spe, dimred="PCA")
    reducedDimNames(spe)
    dim(reducedDim(spe, "UMAP"))
    ### update column names for easier plotting
    colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

    ## visualizations
    ### plot top 2 PCA dimensions
    plotDimRed(spe, type="PCA")
    ### plot top 2 UMAP dimensions
    plotDimRed(spe, type="UMAP")

    # clustering
    ## graph-based clustering using the Walktrap method applied to the top 50 PCs calculated on the set of top HVGs
    set.seed(123)
    k <- 10
    g <- buildSNNGraph(spe, k=k, use.dimred="PCA")
    g_walk <- igraph::cluster_walktrap(g)
    clus <- g_walk$membership
    table(clus)
    ### store cluster labels in column 'label' in colData
    colLabels(spe) <- factor(clus)
    write.csv(colData(spe), file=sprintf("%s/20-Clustering.csv", output_dir), quote=FALSE, row.names=TRUE)

    ## visualizations
    ### plot clusters in PCA reduced dimensions
    clusters_pca <- plotDimRed(spe, type="PCA", annotate="label", palette=color_palette, size=0.5)
    clusters_pca <- clusters_pca + labs(title="Clusters in PCA reduced dimensions", color="Cluster")
    print(clusters_pca)
    ### plot clusters in UMAP reduced dimensions
    clusters_umap <- plotDimRed(spe, type="UMAP", annotate="label", palette=color_palette, size=0.5)
    clusters_umap <- clusters_umap + labs(title="Clusters in UMAP reduced dimensions", color="Cluster")
    sample_plot_list[["clusters_umap"]] <- clusters_umap
    print(clusters_umap)
    ### plot clusters in spatial x-y coordinates
    clusters_spatial <- plotSpots(spe, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="label", palette=color_palette, size=2)
    clusters_spatial <- clusters_spatial + labs(title="Clusters in spatial coordinates", color="Cluster")
    sample_plot_list[["clusters_spatial"]] <- clusters_spatial
    print(clusters_spatial)
    #### modify point size after plot creation
    clusters_pca <- modif_point_size_rm_title(clusters_pca, 0.3)
    clusters_umap <- modif_point_size_rm_title(clusters_umap, 0.3)
    clusters_spatial <- modif_point_size_rm_title(clusters_spatial, 0.5)
    print(wrap_plots(clusters_pca, clusters_spatial) + plot_annotation(title="Clusters in PCA reduced dimensions and spatial coordinates"))
    print(wrap_plots(clusters_umap, clusters_spatial) + plot_annotation(title="Clusters in UMAP reduced dimensions and spatial coordinates"))
    dev.off()

    ################
    # Marker genes #
    ################

    # identify marker genes by testing for differential expression between clusters
    ## test for marker genes
    markers <- findMarkers(spe, test="binom", direction="up")

    ## plot log-fold changes for one cluster over all other clusters
    pdf(sprintf("%s/30-Marker_genes.pdf", output_dir))
    pheatmap_plots <- list()
    ### first, get heatmap for all clusters (trick to avoid empty pdf file)
    for (i in 1:length(names(table(clus)))) {
        one_cluster_markers <- markers[[i]]
        write.csv(one_cluster_markers, file=sprintf("%s/30-Marker_genes_cluster_%d.csv", output_dir, i), quote=FALSE, row.names=TRUE)
        best_set <- one_cluster_markers[one_cluster_markers$Top <= 5, ]
        logFCs <- getMarkerEffects(best_set)
        pheatmap_plot <- pheatmap(logFCs, breaks=seq(-5, 5, length.out=101), silent=TRUE)
        pheatmap_plots[[i]] <- pheatmap_plot[[4]]
    }
    ### plots
    pdf(sprintf("%s/30-Marker_genes.pdf", output_dir))
    for (i in 1:length(names(table(clus)))) {
        #### heatmap
        grid::grid.newpage()
        grid::grid.draw(pheatmap_plots[[i]])
        one_cluster_markers <- markers[[i]]
        #### plot log-transformed normalized expression of top genes for one cluster
        top_genes <- head(rownames(one_cluster_markers))
        print(plotExpression(spe, x="label", features=top_genes))
    }
    dev.off()
    dataset_plot_list[[sample]] <- sample_plot_list
}

# dataset plots
output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/10-Dataset/20-pipelines/00-OSTA", work_dir, genome_name)
if (! dir.exists(output_dir)) {
    dir.create(output_dir, recursive=TRUE, mode="0775")
}

pdf(sprintf("%s/00-QC.pdf", output_dir))
qc_plots_dataset_list <- list()
for (sample in names(dataset_plot_list)) {
    sample_plot <- dataset_plot_list[[sample]][["qc_discard"]]
    sample_plot <- modif_point_size_rm_legend_title_modif_title(sample_plot, 0.01, sample)
    sample_plot <- sample_plot + theme(plot.title=element_text(size=10)) + theme(legend.text=element_text(size=8))
    qc_plots_dataset_list[[sample]] <- sample_plot
}
print(wrap_plots(qc_plots_dataset_list) + plot_annotation(title="Discarded spots"))
dev.off()

pdf(sprintf("%s/10-Normalization.pdf", output_dir))
norm_plots_dataset_list <- list()
for (sample in names(dataset_plot_list)) {
    sample_plot <- dataset_plot_list[[sample]][["norm_lib_size"]]
    sample_plot <- modif_point_size_rm_legend_title_modif_title(sample_plot, 0.01, sample)
    sample_plot <- sample_plot + theme(plot.title=element_text(size=10)) + theme(legend.text=element_text(size=8))
    norm_plots_dataset_list[[sample]] <- sample_plot
}
print(wrap_plots(norm_plots_dataset_list) + plot_annotation(title="Library size per spot after normalization"))
norm_plots_dataset_list <- list()
for (sample in names(dataset_plot_list)) {
    sample_plot <- dataset_plot_list[[sample]][["norm_exp_genes"]]
    sample_plot <- modif_point_size_rm_legend_title_modif_title(sample_plot, 0.01, sample)
    sample_plot <- sample_plot + theme(plot.title=element_text(size=10)) + theme(legend.text=element_text(size=8))
    norm_plots_dataset_list[[sample]] <- sample_plot
}
print(wrap_plots(qc_plots_dataset_list) + plot_annotation(title="Number of expressed genes per spot after normalization"))
dev.off()

pdf(sprintf("%s/20-Clustering.pdf", output_dir))
clusters_plots_dataset_list <- list()
for (sample in names(dataset_plot_list)) {
    sample_plot <- dataset_plot_list[[sample]][["clusters_umap"]]
    sample_plot <- modif_point_size_rm_legend_title_modif_title(sample_plot, 0.01, sample)
    sample_plot <- sample_plot + theme(plot.title=element_text(size=10)) + theme(axis.title=element_text(size=9)) + theme(legend.text=element_text(size=8))
    clusters_plots_dataset_list[[sample]] <- sample_plot
}
print(wrap_plots(clusters_plots_dataset_list) + plot_annotation(title="Clusters in UMAP reduced dimensions"))
clusters_plots_dataset_list <- list()
for (sample in names(dataset_plot_list)) {
    sample_plot <- dataset_plot_list[[sample]][["clusters_spatial"]]
    sample_plot <- modif_point_size_rm_legend_title_modif_title(sample_plot, 0.01, sample)
    sample_plot <- sample_plot + theme(plot.title=element_text(size=10)) + theme(axis.title=element_text(size=9)) + theme(legend.text=element_text(size=8))
    clusters_plots_dataset_list[[sample]] <- sample_plot
}
print(wrap_plots(clusters_plots_dataset_list) + plot_annotation(title="Clusters in spatial coordinates"))
dev.off()
