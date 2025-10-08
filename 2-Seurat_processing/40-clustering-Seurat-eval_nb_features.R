.libPaths(c("/home/christophe.lepriol/NeuroDev_ADD/R/r_4.1.0", .libPaths()))
library(Seurat)
library(scater)
library(patchwork)
library(scran) # for buildSNNGraph() function
library(intrinsicDimension) # for maxLikGlobalDimEst() function
library(bluster) # for bootstrapStability() function
library(tidyr) # for pivot_longer() function
library(tibble) # for rownames_to_column() function
library(cluster) # silhouette() function
library(dplyr) # for mutate() function

#library(RColorBrewer)
#library(pheatmap)


##############
# Parameters #
##############
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"
genome_name <- "NCBIRefSeq108_NCBIRefSeq108GTF"
#samples <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")
sample <- "A_L1_S1"

#color_palette <- c(brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 8, name = "Dark2"), brewer.pal(n = 8, name = "Accent"))


#############
# Functions #
#############
modif_point_size_rm_legend_title_modif_title <- function(one_plot, point_size, one_title) {
    if (! is.na(point_size)) {
        one_plot$layers[[1]]$aes_params$size <- point_size
    }
    one_plot <- one_plot + labs(title=one_title) + theme(legend.title=element_blank())
}


############
# Analysis #
############

library(SpatialExperiment)
library(ggspavis)

print(sprintf("sample: %s", sample))
output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s/20-pipelines/10-Seurat", work_dir, genome_name, sample)
if (! dir.exists(output_dir)) {
    dir.create(output_dir, recursive=TRUE, mode="0775")
}

# Load data

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

# QC

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
discarded_barcodes <- rownames(colData(spe)[colData(spe)$discard,])
plot_qc <- plotQC(spe, type="scatter", metric_x=sprintf("%s.sum", assay_type), metric_y=sprintf("%s.detected", assay_type), threshold_x=threshold_lib_size, threshold_y=threshold_exp_genes)
plot_qc <- plot_qc + labs(x="Library size", y="Number of expressed genes")
print(plot_qc)
plot_qc <- plotQC(spe, type="spots", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", discard="discard")
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

# create Seurat object
spe_seurat <- CreateSeuratObject(counts=counts(spe))
## add image
img <- Read10X_Image(image.dir=file.path(space_ranger_output_dir, "spatial"))
Key(img) <- "image_" # add key for SpatialFeaturePlot() and SpatialPlot() functions
spe_seurat@images <- list(sample=img)
### remove discarded barcodes from image coordinates
spe_seurat@images$sample@coordinates <- spe_seurat@images$sample@coordinates[! rownames(spe_seurat@images$sample@coordinates) %in% discarded_barcodes,]
#spe_seurat@images$sample@assay <- "RNA"
#test <- Load10X_Spatial(data.dir=space_ranger_output_dir)

# unload SpatialExperiment package and namespace to solve conflict between SpatialImage class name from SpatialExperiment and Seurat (makes an error when using SpatialFeaturePlot() function)
detach("package:ggspavis", unload=TRUE) # also unload ggspavis package which imports ‘SpatialExperiment’ namespace
detach("package:SpatialExperiment", unload=TRUE)

# Normalization
normalization_output_dir <- sprintf("%s/10-Normalization", output_dir)
if (! dir.exists(normalization_output_dir)) {
    dir.create(normalization_output_dir, recursive=TRUE, mode="0775")
}

nb_features_list <- list()
for (nb_features in c(100, 250, 500, 1000, 1500, 2000, 3000)) {
    print(sprintf("nb features: %d", nb_features))
    normalization_output_sudir <- sprintf("%s/nb_features_evaluation/%s_features", normalization_output_dir, nb_features)
    if (! dir.exists(normalization_output_sudir)) {
        dir.create(normalization_output_sudir, recursive=TRUE, mode="0775")
    }
    
    # normalization
    spe_seurat <- SCTransform(spe_seurat, assay = "RNA", variable.features.n=nb_features)
    # counts per spot
    pdf(sprintf("%s/10-Normalization.pdf", normalization_output_sudir))
    ## raw counts
    nCount_Spatial_plot1 <- VlnPlot(spe_seurat, features = "nCount_RNA", pt.size = 0.1) + NoLegend()
    nCount_Spatial_plot2 <- SpatialFeaturePlot(spe_seurat, features = "nCount_RNA") + theme(legend.position = "right")
    print(wrap_plots(nCount_Spatial_plot1, nCount_Spatial_plot2) + plot_annotation(title="Raw counts"))
    ## normalized counts
    nCount_SCT_plot1 <- VlnPlot(spe_seurat, features = "nCount_SCT", pt.size = 0.1) + NoLegend()
    nCount_SCT_plot2 <- SpatialFeaturePlot(spe_seurat, features = "nCount_SCT") + theme(legend.position = "right")
    print(wrap_plots(nCount_SCT_plot1, nCount_SCT_plot2) + plot_annotation(title="Normalized counts"))
    ## raw and normalized nCount plots in the same page
    print(wrap_plots(nCount_Spatial_plot2, nCount_SCT_plot2) + plot_annotation(title="Raw and normalized counts"))

    # comparison with log-normalization: https://satijalab.org/seurat/articles/spatial_vignette.html
    ## rerun normalization to store sctranform residuals for all genes
    spe_seurat_norm_comp <- SCTransform(spe_seurat, assay = "RNA", variable.features.n=nb_features, return.only.var.genes=FALSE, verbose=FALSE)
    ## also run standard log normalization for comparison
    spe_seurat_norm_comp <- NormalizeData(spe_seurat_norm_comp, assay="RNA", verbose=FALSE)
    ## compute the correlation between the log normalized data and sctransform residuals with the number of UMIs
    spe_seurat_norm_comp <- GroupCorrelation(spe_seurat_norm_comp, group.assay="RNA", assay="RNA", slot="data", do.plot=FALSE)
    spe_seurat_norm_comp <- GroupCorrelation(spe_seurat_norm_comp, group.assay="RNA", assay="SCT", slot="scale.data", do.plot=FALSE)
    spe_seurat <- GroupCorrelation(spe_seurat, group.assay="RNA", assay="SCT", slot="scale.data", do.plot=FALSE)
    p1 <- GroupCorrelationPlot(spe_seurat_norm_comp, assay="RNA", cor="nCount_RNA_cor") + ggtitle("Log normalization") + theme(plot.title=element_text(hjust=0.5))
    p2 <- GroupCorrelationPlot(spe_seurat_norm_comp, assay="SCT", cor="nCount_RNA_cor") + ggtitle("SCTransform Normalization\nscale data for all genes") + theme(plot.title=element_text(hjust=0.5))
    p3 <- GroupCorrelationPlot(spe_seurat, assay="SCT", cor="nCount_RNA_cor") + ggtitle("SCTransform Normalization\n scale data only for the\nvariable genes") + theme(plot.title=element_text(hjust=0.5))
    print(p1 + p2)
    print(p1 + p3)
    print(p2 + p3)
    dev.off()
    
    # Clustering
    clustering_output_sudir <- sprintf("%s/20-Clustering/nb_features_evaluation/%s_features", output_dir, nb_features)
    if (! dir.exists(clustering_output_sudir)) {
        dir.create(clustering_output_sudir, recursive=TRUE, mode="0775")
    }
    
    # dimensionality reduction and clustering
    total_pcs <- 50
    seurat_spatial_vignette_pcs <- 30 # https://satijalab.org/seurat/articles/spatial_vignette.html
    spe_seurat <- RunPCA(spe_seurat, assay="SCT", npcs=total_pcs, verbose=FALSE)
    ## estimate intrinsic dimension
    dim_estimate <- maxLikGlobalDimEst(spe_seurat@reductions$pca@cell.embeddings, k=10)
    nb_pcs <- round(dim_estimate$dim.est) + 5
    df2ggplot <- data.frame(PC=1:total_pcs, stdev=spe_seurat@reductions$pca@stdev)
    write.csv(df2ggplot, file=sprintf("%s/PCs_stdev.csv", clustering_output_sudir), quote=FALSE, row.names=TRUE)
    pdf(sprintf("%s/PCs_stdev.pdf", clustering_output_sudir))
    p <- ggplot(df2ggplot, aes(PC, stdev)) +
        geom_point() +
        geom_point() +
        geom_vline(xintercept = round(dim_estimate$dim.est), color = "blue") +
        geom_vline(xintercept = nb_pcs, color = "red") +
        geom_vline(xintercept = seurat_spatial_vignette_pcs, color = "green") +
        theme_bw() +
        labs(x = "Principal components", y = "Standard deviation")
    print(p)
    dev.off()
    
    spe_seurat <- FindNeighbors(spe_seurat, reduction="pca", dims=1:nb_pcs)
    spe_seurat <- FindClusters(spe_seurat, verbose=FALSE)
    spe_seurat <- RunUMAP(spe_seurat, reduction="pca", dims=1:nb_pcs)
    SCT_cluster_plot1 <- DimPlot(spe_seurat, reduction="umap", label=TRUE)
    SCT_cluster_plot2 <- SpatialDimPlot(spe_seurat, label=TRUE, label.size=3)
    pdf(sprintf("%s/clustering_dimension_spatial.pdf", clustering_output_sudir))
    print(SCT_cluster_plot1 + SCT_cluster_plot2 + plot_annotation(title="Normalized counts"))
    print(SCT_cluster_plot1)
    print(SCT_cluster_plot2)
    dev.off()
    one_nb_features_list <- list()
    one_nb_features_list[["spatilal_dimplot"]] <-  SCT_cluster_plot2
    SCT_cluster_plot2_no_legend <- SCT_cluster_plot2 + NoLegend()
    one_nb_features_list[["spatilal_dimplot_no_legend"]] <-  SCT_cluster_plot2_no_legend

    # cluster stability: https://romanhaa.github.io/projects/scrnaseq_workflow/#cluster-stability
    sce <- as.SingleCellExperiment(spe_seurat)

    #reducedDim(sce, 'PCA_sub') <- reducedDim(sce, 'PCA')[,1:15, drop = FALSE]

    #ass_prob <- bootstrapStability(sce, FUN = function(x) {
    #    g <- buildSNNGraph(x, use.dimred = "PCA_sub")
    #    igraph::cluster_walktrap(g)$membership
    #}, clusters = sce$seurat_clusters)
    ass_prob <- bootstrapStability(sce, FUN = function(x) {
        g <- buildSNNGraph(x, use.dimred = "PCA")
        igraph::cluster_walktrap(g)$membership
    }, clusters = sce$seurat_clusters, adjusted=FALSE)

    df2ggplot <- ass_prob %>%
      as_tibble() %>%
      rownames_to_column(var = 'cluster_1') %>%
      pivot_longer(
        cols = 2:ncol(.),
        names_to = 'cluster_2',
        values_to = 'probability'
      ) %>%
      mutate(
        cluster_1 = as.character(as.numeric(cluster_1) - 1),
        cluster_1 = factor(cluster_1, levels = rev(unique(cluster_1))),
        cluster_2 = factor(cluster_2, levels = unique(cluster_2))
      )
    write.csv(df2ggplot, file=sprintf("%s/clustering_stability.csv", clustering_output_sudir), quote=FALSE, row.names=TRUE)

    p <- ggplot(df2ggplot, aes(cluster_2, cluster_1, fill = probability)) +
      geom_tile(color = 'white') +
      geom_text(aes(label = round(probability, digits = 2)), size = 2.5) +
      scale_x_discrete(name = 'Cluster', position = 'top') +
      scale_y_discrete(name = 'Cluster') +
      scale_fill_gradient(
        name = 'Probability', low = 'white', high = '#c0392b', na.value = '#bdc3c7',
        limits = c(0,1),
        guide = guide_colorbar(
          frame.colour = 'black', ticks.colour = 'black', title.position = 'left',
          title.theme = element_text(hjust = 1, angle = 90),
          barwidth = 0.75, barheight = 10
        )
      ) +
      coord_fixed() +
      theme_bw() +
      theme(
        legend.position = 'right',
        panel.grid.major = element_blank()
      )
    pdf(sprintf("%s/clustering_stability.pdf", clustering_output_sudir))
    print(p)
    dev.off()

    # silhouette: https://romanhaa.github.io/projects/scrnaseq_workflow/#silhouette-plot
    distance_matrix <- dist(Embeddings(spe_seurat[['pca']]))
    clusters <- spe_seurat@meta.data$seurat_clusters
    silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
    spe_seurat@meta.data$silhouette_score <- silhouette[,3]
    mean_silhouette_score <- mean(spe_seurat@meta.data$silhouette_score)
    print(sprintf("mean silhouette score: %f", mean_silhouette_score))
    one_nb_features_list[["mean_silhouette_score"]] <- mean_silhouette_score

    df2ggplot <- spe_seurat@meta.data %>%
      mutate(barcode = rownames(.)) %>%
      arrange(seurat_clusters,-silhouette_score) %>%
      mutate(barcode = factor(barcode, levels = barcode))
    write.csv(df2ggplot, file=sprintf("%s/clustering_silhouette.csv", clustering_output_sudir), quote=FALSE, row.names=TRUE)

    p <- ggplot(df2ggplot) +
      geom_col(aes(barcode, silhouette_score, fill = seurat_clusters), show.legend = FALSE) +
      geom_hline(yintercept = mean_silhouette_score, color = 'red', linetype = 'dashed') +
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
    pdf(sprintf("%s/clustering_silhouette.pdf", clustering_output_sudir))
    print(p)
    dev.off()

    # cluster similarity: https://romanhaa.github.io/projects/scrnaseq_workflow/#cluster-similarity
    sce <- as.SingleCellExperiment(spe_seurat)
    g <- buildSNNGraph(sce, use.dimred = "PCA")
    ratio <- pairwiseModularity(g, spe_seurat@meta.data$seurat_clusters, as.ratio = TRUE)
    ratio_to_plot <- log10(ratio+1)

    df2ggplot <- ratio_to_plot %>%
      as_tibble() %>%
      rownames_to_column(var = 'cluster_1') %>%
      pivot_longer(
        cols = 2:ncol(.),
        names_to = 'cluster_2',
        values_to = 'probability'
      ) %>%
      mutate(
        cluster_1 = as.character(as.numeric(cluster_1) - 1),
        cluster_1 = factor(cluster_1, levels = rev(unique(cluster_1))),
        cluster_2 = factor(cluster_2, levels = unique(cluster_2))
      )
    write.csv(df2ggplot, file=sprintf("%s/clustering_similarity.csv", clustering_output_sudir), quote=FALSE, row.names=TRUE)

    p <- ggplot(df2ggplot, aes(cluster_2, cluster_1, fill = probability)) +
      geom_tile(color = 'white') +
      geom_text(aes(label = round(probability, digits = 2)), size = 2.5) +
      scale_x_discrete(name = 'Cluster', position = 'top') +
      scale_y_discrete(name = 'Cluster') +
      scale_fill_gradient(
        name = 'log10(ratio)', low = 'white', high = '#c0392b', na.value = '#bdc3c7',
        guide = guide_colorbar(
          frame.colour = 'black', ticks.colour = 'black', title.position = 'left',
          title.theme = element_text(hjust = 1, angle = 90),
          barwidth = 0.75, barheight = 10
        )
      ) +
      coord_fixed() +
      theme_bw() +
      theme(
        legend.position = 'right',
        panel.grid.major = element_blank()
      )
    pdf(sprintf("%s/clustering_similarity.pdf", clustering_output_sudir))
    print(p)
    dev.off()

    nb_features_list[[sprintf("%d_features", nb_features)]] <- one_nb_features_list
}

spatial_dimplot <- list()
mean_silhouette_score_vector <- c()
for (nb_features in c(100, 250, 500, 1000, 1500, 2000, 3000)) {
    # add title 
    spatial_dimplot[[sprintf("%d_features", nb_features)]] <- nb_features_list[[sprintf("%d_features", nb_features)]][["spatilal_dimplot_no_legend"]] + labs(title=sprintf("%d features: %f", nb_features, nb_features_list[[sprintf("%d_features", nb_features)]][["mean_silhouette_score"]]))
    #print(spatial_dimplot[[sprintf("%d_PCs", dim2use)]])
    mean_silhouette_score_vector <- c(mean_silhouette_score_vector, nb_features_list[[sprintf("%d_features", nb_features)]][["mean_silhouette_score"]])
}
clustering_output_dir <- sprintf("%s/20-Clustering/nb_features_evaluation", output_dir)
pdf(sprintf("%s/nb_features_clustering_dimension_spatial.pdf", clustering_output_dir))
print(wrap_plots(spatial_dimplot) + plot_annotation(title="Number of features comparison"))
dev.off()


df2ggplot <- data.frame(features=c(100, 250, 500, 1000, 1500, 2000, 3000), silhouette=mean_silhouette_score_vector)
write.csv(df2ggplot, file=sprintf("%s/nb_features_mean_silhouette_score.csv", clustering_output_dir), quote=FALSE, row.names=TRUE)
pdf(sprintf("%s/nb_features_mean_silhouette_score.pdf", clustering_output_dir))
p <- ggplot(df2ggplot, aes(x=features, y=silhouette)) +
    geom_point() +
    geom_line() +
    labs(x="Number of PCs", y="mean silhouette score") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50"))
print(p)
dev.off()
