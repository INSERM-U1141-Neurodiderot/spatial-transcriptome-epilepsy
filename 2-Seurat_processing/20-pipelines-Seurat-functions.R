load_ST_pipeline_data <- function(sample, st_dir, sr_dir, metadata_df) {
    # input data
    ## ST Pipeline count matrix
    st_pipeline_matrix_file <- sprintf("%s/%s/10-Pipeline/%s_stdata.tsv", st_dir, sample, sample)
    st_pipeline_matrix <- read.table(st_pipeline_matrix_file, sep="\t", header=TRUE, quote="", row.names=1)
    ## Space Ranger output directory for tissue positions list and image
    space_ranger_out_dir <- sprintf("%s/%s/10-Pipeline/outs", sr_dir, sample)

    # build SpatialExperiment object
    spatial_coordinates_file <- file.path(space_ranger_out_dir, "spatial", "tissue_positions_list.csv")
    spatial_coordinates_df <- read.csv(spatial_coordinates_file, header=FALSE, quote="")
    colnames(spatial_coordinates_df) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    rownames(spatial_coordinates_df) <- spatial_coordinates_df$barcode
    ## add lame, zone and condition
    nb_spots <- dim(spatial_coordinates_df)[1]
    sample_strsplit <- unlist(strsplit(sample, "_"))
    zone <- sample_strsplit[1]
    lame <- sample_strsplit[2]
    condition <- metadata_df[which(metadata_df$sample==sample), "condition"]
    timepoint <- metadata_df[which(metadata_df$sample==sample), "time"]
    spatial_coordinates_df <- cbind(spatial_coordinates_df, lame=rep(lame, nb_spots), zone=rep(zone, nb_spots), condition=rep(condition, nb_spots), time=rep(timepoint, nb_spots))
    spatial_coordinates_df$zone <- as.factor(spatial_coordinates_df$zone)
    spatial_coordinates_df$lame <- as.factor(spatial_coordinates_df$lame)
    spatial_coordinates_df$condition <- as.factor(spatial_coordinates_df$condition)
    spatial_coordinates_df$time <- as.factor(spatial_coordinates_df$time)
    spatial_coordinates_df <- droplevels(spatial_coordinates_df)
    ## matrix: convert coordinates to barcodes
    barcodes <- unlist(lapply(rownames(st_pipeline_matrix), function(x, coord2barcode=spatial_coordinates_df) {
        coordinates <- unlist(strsplit(x, "x"))
        return(coord2barcode[which(coord2barcode$array_col==as.integer(coordinates[1])-1 & coord2barcode$array_row==as.integer(coordinates[2])-1), "barcode"])
    }))
    rownames(st_pipeline_matrix) <- barcodes
    gene_data <- data.frame(gene_name=colnames(st_pipeline_matrix))
    ## image data
    img <- readImgData(imageSources=file.path(space_ranger_out_dir, "spatial", "tissue_lowres_image.png"), scaleFactors=file.path(space_ranger_out_dir, "spatial", "scalefactors_json.json"), sample_id=sample)
    
    # create SpatialExperiment object
    spe <- SpatialExperiment(assay=list(counts=t(st_pipeline_matrix)), rowData=gene_data, colData=spatial_coordinates_df[rownames(st_pipeline_matrix),], imgData=img, spatialDataNames=c("barcode", "in_tissue", "array_row", "array_col"), spatialCoordsNames=c("pxl_col_in_fullres", "pxl_row_in_fullres"), sample_id=sample)

    # keep spot over tissue
    spe <- spe[, colData(spe)$in_tissue==1]
    
    # remove non-expressed genes: genes with 0 count for all spots
    non_expressed <- apply(spe@assays@data$counts, 1, sum) == 0
    spe <- spe[!non_expressed, ]
    
    return(spe)
}

######
# QC #
######
density_plot <- function(data2plot, plot_title, plot_x_lab) {
    df2ggplot <- data.frame(value=data2plot)
    percentiles <- quantile(data2plot, probs=seq(0,1,0.05))
    percentiles2plot <- c("5%", "10%", "25%", "50%", "75%", "90%", "95%")
    percentiles_df2ggplot <- data.frame(percentile=percentiles2plot, exp=percentiles[percentiles2plot])
    density_plot <- ggplot(df2ggplot, aes(x=value)) +
        geom_line(stat="density") +
        geom_vline(data=percentiles_df2ggplot[which(percentiles_df2ggplot$percentile == "50%"),], aes(xintercept=exp), linetype="longdash", alpha=0.7) +
        geom_vline(data=percentiles_df2ggplot[which(percentiles_df2ggplot$percentile %in% c("25%", "75%")),], aes(xintercept=exp), linetype="dashed", alpha=0.7) +
        geom_vline(data=percentiles_df2ggplot[which(percentiles_df2ggplot$percentile %in% c("10%", "90%")),], aes(xintercept=exp), linetype="dotdash", alpha=0.7) +
        geom_vline(data=percentiles_df2ggplot[which(percentiles_df2ggplot$percentile %in% c("5%", "95%")),], aes(xintercept=exp), linetype="dotted", alpha=0.7) +
        labs(title=plot_title, x=plot_x_lab, y="Density") +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50"))
    print(density_plot)
    
    return(percentiles)
}

qc_threshold <- function(sp_exp, colname, thres, gt_bool, qc_colname) {
    if (gt_bool) {
        spot_qc <- colData(sp_exp)[[colname]] > thres
    } else {
        spot_qc <- colData(sp_exp)[[colname]] < thres
    }
    colData(sp_exp)$spot_qc <- spot_qc
    colnames(colData(sp_exp))[which(colnames(colData(sp_exp)) == "spot_qc")] <- qc_colname
    return(sp_exp)
}

qc_threshold_plot <- function(sp_exp, colname) {
    plot_qc <- plotQC(sp_exp, type="spots", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", discard=colname)
    plot_qc$layers[[1]]$aes_params$size <- 2
    print(plot_qc)
    return(plot_qc)
}



modif_point_size_rm_legend_title_modif_title <- function(one_plot, point_size, one_title) {
    if (! is.na(point_size)) {
        one_plot$layers[[1]]$aes_params$size <- point_size
    }
    one_plot <- one_plot + labs(title=one_title) + theme(legend.title=element_blank())
}

QC <- function(spe, lib_size, nb_exp_genes, space_ranger_out_dir, out_dir, out_basename) {
    library(SpatialExperiment)
    library(ggspavis)
    
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
    pdf(sprintf("%s/%s.pdf", out_dir, out_basename))
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
    qc_lib_size <- colData(spe)[, sprintf("%s.sum", assay_type)] < lib_size
    table(qc_lib_size)
    colData(spe)$qc_lib_size <- qc_lib_size
    colnames(colData(spe))[which(colnames(colData(spe)) == "qc_lib_size")] <- sprintf("qc_lib_size_%d", lib_size)
    plot_qc <- plotQC(spe, type="spots", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", discard=sprintf("qc_lib_size_%d", lib_size))
    plot_qc <- modif_point_size_rm_legend_title_modif_title(plot_qc, 2, sprintf("QC spots: library size < %d", lib_size))
    print(plot_qc)
    ## number of expressed features
    nb_exp_genes <- 1000
    qc_detected <- colData(spe)[, sprintf("%s.detected", assay_type)] < nb_exp_genes
    table(qc_detected)
    colData(spe)$qc_detected <- qc_detected
    colnames(colData(spe))[which(colnames(colData(spe)) == "qc_detected")] <- sprintf("qc_detected_%d", nb_exp_genes)
    plot_qc <- plotQC(spe, type="spots", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", discard=sprintf("qc_detected_%d", nb_exp_genes))
    plot_qc <- modif_point_size_rm_legend_title_modif_title(plot_qc, 2, sprintf("QC spots: number of expressed genes < %d", nb_exp_genes))
    print(plot_qc)

    # remove low-quality spots
    ## number of discarded spots for each metric
    apply(cbind(colData(spe)[[sprintf("qc_lib_size_%d", lib_size)]], colData(spe)[[sprintf("qc_detected_%d", nb_exp_genes)]]), 2, sum)
    ## combined set of discarded spots
    discard <- colData(spe) [[sprintf("qc_lib_size_%d", lib_size)]] | colData(spe)[[sprintf("qc_detected_%d", nb_exp_genes)]]
    table(discard)
    ## store in object
    colData(spe)$discard <- discard
    discarded_barcodes <- rownames(colData(spe)[colData(spe)$discard,])
    plot_qc <- plotQC(spe, type="scatter", metric_x=sprintf("%s.sum", assay_type), metric_y=sprintf("%s.detected", assay_type), threshold_x=lib_size, threshold_y=nb_exp_genes)
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
    write.csv(colData(spe), file=sprintf("%s/%s.csv", out_dir, out_basename), quote=FALSE, row.names=TRUE)

    # create Seurat object
    seurat_obj <- CreateSeuratObject(counts=counts(spe))
    ## add image
    img <- Read10X_Image(image.dir=file.path(space_ranger_out_dir, "spatial"))
    Key(img) <- "image_" # add key for SpatialFeaturePlot() and SpatialPlot() functions
    seurat_obj@images <- list(sample=img)
    ### remove discarded barcodes from image coordinates
    seurat_obj@images$sample@coordinates <- seurat_obj@images$sample@coordinates[! rownames(seurat_obj@images$sample@coordinates) %in% discarded_barcodes,]
    #spe_seurat@images$sample@assay <- "RNA"
    seurat_obj@images$sample@assay <- seurat_obj@active.assay
    #test <- Load10X_Spatial(data.dir=space_ranger_out_dir)

    # unload SpatialExperiment package and namespace to solve conflict between SpatialImage class name from SpatialExperiment and Seurat (makes an error when using SpatialFeaturePlot() function)
    detach("package:ggspavis", unload=TRUE) # also unload ggspavis package which imports ‘SpatialExperiment’ namespace
    detach("package:SpatialExperiment", unload=TRUE)
    
    return(seurat_obj)
}

normalization_plots <- function(seurat_obj, feature_name, suffix) {
    plots2return <- list()
    violin_plot <- VlnPlot(seurat_obj, features = sprintf("%s_%s", feature_name, suffix), pt.size = 0.1) + NoLegend()
    spatial_feature_plot <- SpatialFeaturePlot(seurat_obj, features = sprintf("%s_%s", feature_name, suffix)) + theme(legend.position = "right")
    plots2return[["violin"]] <- violin_plot
    plots2return[["spatial"]] <- spatial_feature_plot
    return(plots2return)
}

get_normalization_subdir_basename <- function(norm_method) {
    list2return <- list()
    norm_out_name <- sprintf("Normalization-%s", norm_method)
    if (norm_method == "SCTransform") {
        norm_method_dir_index <- "00"
    } else {
        if (norm_method == "LogNormalize") {
            norm_method_dir_index <- "10"
        } else {
            if (norm_method == "CLR") {
                norm_method_dir_index <- "20"
            } else {
                if (norm_method == "RC") {
                    norm_method_dir_index <- "30"
                }
            }
        }
    }
    norm_subdir <- sprintf("%s-%s", norm_method_dir_index, norm_out_name)
    list2return[["subdir"]] <- norm_subdir
    list2return[["basename"]] <- norm_out_name
    return(list2return)
}

get_features_subdir_basename <- function(features_method, features_nb, features_resvar) {
    list2return <- list()
    if (features_method == "HVGs") {
        features_method_dir_index <- "00"
        features_method_out_basename <- sprintf("Features-%s", features_method)
        if (!is.na(features_nb) & is.na(features_resvar)) {
            features_param_out_basename <- sprintf("nb%d", features_nb)
        } else {
            if (!is.na(features_resvar) & is.na(features_nb)) {
                features_resvar_in_dir <- sub("[.]", "_", features_resvar)
                features_param_out_basename <- sprintf("residualvariance%s", features_resvar_in_dir)
            }
        }
        features_out_basename <- sprintf("%s_%s", features_method_out_basename, features_param_out_basename)
    } else {
        if (features_method == "SVGs") {
            features_method_dir_index <- "10"
            features_method_out_basename <- sprintf("00-Features-%s", features_method)
            # TODO
        }
    }
    features_method_subdir <- sprintf("%s-%s", features_method_dir_index, features_method_out_basename)
    
    list2return[["method_subdir"]] <- features_method_subdir
    list2return[["param_subdir"]] <- features_param_out_basename
    list2return[["method_basename"]] <- features_method_out_basename
    list2return[["param_basename"]] <- features_param_out_basename
    return(list2return)
}

normalization_feature_selection <- function(seurat_obj, norm_method, features_method, features_nb, features_resvar, norm_out_dir, norm_out_name, features_out_dir, features_out_name) {
    
    if (norm_method == "SCTransform") {
        if (features_method == "HVGs") {
            if (! is.na(features_nb) & is.na(features_resvar)) {
                seurat_obj <- SCTransform(seurat_obj, assay="RNA", variable.features.n=features_nb)
            } else {
                if (! is.na(features_resvar) & is.na(features_nb)) {
                    seurat_obj <- SCTransform(seurat_obj, assay="RNA", variable.features.n=NULL, variable.features.rv.th=features_resvar, return.only.var.genes=FALSE)
                }
            }
        }
        features_suffix2plot <- "SCT"
        selection_method <- "sct"
        assay_name <- "SCT"
    } else {
        # NormalizeData() function
        seurat_obj <- NormalizeData(seurat_obj, normalization.method=norm_method)
        # compute sum of counts and number of expressed genes after normalization
        counts <- as.matrix(seurat_obj@assays$RNA@data)
        nCount <- colSums(counts)
        seurat_obj$nCount_norm <- nCount
        nFeature <- apply(counts, 2, function(x) { length(x[x>0]) })
        seurat_obj$nFeature_norm <- nFeature
        features_suffix2plot <- "norm"
        selection_method <- "vst"
        assay_name <- "RNA"
        if (features_method == "HVGs") {
            if (! is.na(features_nb) & is.na(features_resvar)) {
                seurat_obj <- FindVariableFeatures(seurat_obj, selection.method=selection_method, nfeatures=features_nb)
            }
        }
    }
    
    # plots
    ## normalization
    pdf(sprintf("%s/%s.pdf", norm_out_dir, norm_out_name))
    ## counts per spot
    ### raw counts
    raw_plots <- normalization_plots(seurat_obj, "nCount", "RNA")
    print((wrap_plots(raw_plots[["violin"]], raw_plots[["spatial"]])) + plot_annotation(title="Raw counts"))
    ### normalized counts
    norm_plots <- normalization_plots(seurat_obj, "nCount", features_suffix2plot)
    print((wrap_plots(norm_plots[["violin"]], norm_plots[["spatial"]])) + plot_annotation(title="Normalized counts"))
    ## number of expressed genes per spot
    ### raw counts
    raw_plots <- normalization_plots(seurat_obj, "nFeature", "RNA")
    print((wrap_plots(raw_plots[["violin"]], raw_plots[["spatial"]])) + plot_annotation(title="Raw counts"))
    ### normalized counts
    norm_plots <- normalization_plots(seurat_obj, "nFeature", features_suffix2plot)
    print((wrap_plots(norm_plots[["violin"]], norm_plots[["spatial"]])) + plot_annotation(title="Normalized counts"))

    if (norm_method == "SCTransform") {
        # add comparison with log-normalization: https://satijalab.org/seurat/articles/spatial_vignette.html
        ## rerun normalization to store sctranform residuals for all genes
        seurat_obj_norm_comp <- SCTransform(seurat_obj, assay = "RNA", return.only.var.genes=FALSE, verbose=FALSE)
        ## also run standard log normalization for comparison
        seurat_objnorm_comp <- NormalizeData(seurat_obj_norm_comp, assay="RNA", verbose=FALSE)
        ## compute the correlation between the log normalized data and sctransform residuals with the number of UMIs
        seurat_obj_norm_comp <- GroupCorrelation(seurat_obj_norm_comp, group.assay="RNA", assay="RNA", slot="data", do.plot=FALSE)
        seurat_obj_norm_comp <- GroupCorrelation(seurat_obj_norm_comp, group.assay="RNA", assay="SCT", slot="scale.data", do.plot=FALSE)
        seurat_obj <- GroupCorrelation(seurat_obj, group.assay="RNA", assay="SCT", slot="scale.data", do.plot=FALSE)
        p1 <- GroupCorrelationPlot(seurat_obj_norm_comp, assay="RNA", cor="nCount_RNA_cor") + ggtitle("Log normalization") + theme(plot.title=element_text(hjust=0.5))
        p2 <- GroupCorrelationPlot(seurat_obj_norm_comp, assay="SCT", cor="nCount_RNA_cor") + ggtitle("SCTransform Normalization\nscale data for all genes") + theme(plot.title=element_text(hjust=0.5))
        p3 <- GroupCorrelationPlot(seurat_obj, assay="SCT", cor="nCount_RNA_cor") + ggtitle("SCTransform Normalization\n scale data only for the\nvariable genes") + theme(plot.title=element_text(hjust=0.5))
        print(p1 + p2)
        print(p1 + p3)
        print(p2 + p3)
    }
    dev.off()

    ## feature selection
    ### Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(seurat_obj), 10)
    ### plot variable features with and without labels
    pdf(sprintf("%s/%s.pdf", features_out_dir, features_out_name))
    plot1 <- VariableFeaturePlot(seurat_obj, selection.method=selection_method, assay=assay_name)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    print(plot2)
    dev.off()

    # scale data
    if (norm_method != "SCTransform" & features_method == "HVGs") {
        # linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA
        all.genes <- rownames(seurat_obj)
        seurat_obj <- ScaleData(seurat_obj, features=all.genes)
    }
    
    return(seurat_obj)
}

bootstrap_myclusters <- function(x, FUN, clusters=NULL, transposed=FALSE, n.cells=5000, 
                                 iterations=50, ...) {
  if (is.null(clusters)) {
    clusters <- FUN(x, ...)
  }
  cluster.ids <- as.character(sort(unique(clusters)))
  output <- matrix(0, length(cluster.ids), length(cluster.ids))
  output[lower.tri(output)] <- NA_real_
  dimnames(output) <- list(cluster.ids, cluster.ids)
  
  for (i in seq_len(iterations)) {
    if (transposed) {
      chosen <- sample(nrow(x), n.cells)
      resampled <- x[chosen,]
    } else {
      chosen <- sample(ncol(x), n.cells)
      resampled <- x[,chosen]
    }
    reclusters <- FUN(resampled, ...)
    tab <- table(clusters[chosen], reclusters)
    for (j1 in seq_along(cluster.ids)) {
      spread1 <- tab[cluster.ids[j1],]
      spread1 <- spread1/sum(spread1)
      for (j2 in seq_len(j1)) {
        spread2 <- tab[cluster.ids[j2],]
        spread2 <- spread2/sum(spread2)
        output[j2,j1] <- output[j2,j1] + sum(spread1 * spread2)/iterations
      }
    }
  }
  
  output
}

seurat_clustering <- function(x, dims, knn, rez) {
    x <- FindNeighbors(x, reduction='pca', dims=1:dims, k.param=knn) 
    x <- FindClusters(x, resolution=rez)
    return(as.numeric(x$seurat_clusters))
}

clustering_analysis <- function(seurat_obj, dims, knn, rez, nb_spots, iter, out_dir, out_name) {
    sce <- as.SingleCellExperiment(seurat_obj)
    # cluster stability: https://romanhaa.github.io/projects/scrnaseq_workflow/#cluster-stability
    #ass_prob <- bootstrapStability(sce, FUN = function(x) {
    #    g <- buildSNNGraph(x, use.dimred = "PCA")
    #    igraph::cluster_walktrap(g)$membership
    #}, clusters = sce$seurat_clusters, adjusted=FALSE)

    #df2ggplot <- ass_prob %>%
    #  as_tibble() %>%
    #  rownames_to_column(var = 'cluster_1') %>%
    #  pivot_longer(
    #    cols = 2:ncol(.),
    #    names_to = 'cluster_2',
    #    values_to = 'probability'
    #  ) %>%
    #  mutate(
    #    cluster_1 = as.character(as.numeric(cluster_1) - 1),
    #    cluster_1 = factor(cluster_1, levels = rev(unique(cluster_1))),
    #    cluster_2 = factor(cluster_2, levels = unique(cluster_2))
    #  )
    #write.csv(df2ggplot, file=sprintf("%s/%s_stability.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)

    #p <- ggplot(df2ggplot, aes(cluster_2, cluster_1, fill = probability)) +
    #  geom_tile(color = 'white') +
    #  geom_text(aes(label = round(probability, digits = 2)), size = 2.5) +
    #  labs(title="Cluster stability") +
    #  scale_x_discrete(name = 'Cluster', position = 'top') +
    #  scale_y_discrete(name = 'Cluster') +
    #  scale_fill_gradient(
    #    name = 'Probability', low = 'white', high = '#c0392b', na.value = '#bdc3c7',
    #    limits = c(0,1),
    #    guide = guide_colorbar(
    #      frame.colour = 'black', ticks.colour = 'black', title.position = 'left',
    #      title.theme = element_text(hjust = 1, angle = 90),
    #      barwidth = 0.75, barheight = 10
    #    )
    #  ) +
    #  coord_fixed() +
    #  theme_bw() +
    #  theme(
    #    legend.position = 'right',
    #    panel.grid.major = element_blank()
    #  )
    #pdf(sprintf("%s/%s_plots.pdf", out_dir, out_name))
    #print(p)
    
    # cluster stability: https://github.com/PrashINRA/BootStrap_SingleCell
    originals<- seurat_obj$seurat_clusters #This is the cluster or CellType information, you already have stored in Seurat object
    coassign <- bootstrap_myclusters(seurat_obj, clusters = originals, FUN = seurat_clustering, dims=dims, knn=knn, rez=rez, n.cells=nb_spots, iterations=iter)
    heatmap_plot <- pheatmap(coassign, cluster_row=F, cluster_col=F, main= "Coassignment probabilities", angle_col = 45,
             color=rev(viridis::magma(100)), breaks=seq(0, 1, length.out=101))
    
    # silhouette: https://romanhaa.github.io/projects/scrnaseq_workflow/#silhouette-plot
    distance_matrix <- dist(Embeddings(seurat_obj[['pca']]))
    clusters <- seurat_obj@meta.data$seurat_clusters
    silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
    seurat_obj@meta.data$silhouette_score <- silhouette[,3]
    mean_silhouette_score <- mean(seurat_obj@meta.data$silhouette_score)
    print(sprintf("mean silhouette score: %f", mean_silhouette_score))

    df2ggplot <- seurat_obj@meta.data %>%
      mutate(barcode = rownames(.)) %>%
      arrange(seurat_clusters,-silhouette_score) %>%
      mutate(barcode = factor(barcode, levels = barcode))
    write.csv(df2ggplot, file=sprintf("%s/%s_silhouette.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
    
    silhouette_score_plot <- ggplot(df2ggplot) +
      geom_col(aes(barcode, silhouette_score, fill = seurat_clusters), show.legend = FALSE) +
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
    
    ## silhouette score on UMAP plot
    cluster_plot <- DimPlot(seurat_obj, reduction="umap", label=TRUE) + labs(title="Clustering labels")
    silhouette_plot <- FeaturePlot(seurat_obj, features="silhouette_score", label=TRUE) + scale_colour_gradientn(colours=rev(brewer.pal(n=11, name="RdYlBu")), limits=c(-1,1)) + labs(title="Silhouette score")
    silhouette_umap_plot <- cluster_plot / silhouette_plot

    # cluster similarity: https://romanhaa.github.io/projects/scrnaseq_workflow/#cluster-similarity
    g <- buildSNNGraph(sce, use.dimred = "PCA")
    ratio <- pairwiseModularity(g, seurat_obj@meta.data$seurat_clusters, as.ratio = TRUE)
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
    write.csv(df2ggplot, file=sprintf("%s/%s_similarity.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)

    similarity_plot <- ggplot(df2ggplot, aes(cluster_2, cluster_1, fill = probability)) +
      geom_tile(color = 'white') +
      geom_text(aes(label = round(probability, digits = 2)), size = 2.5) +
      labs(title="Cluster similarity") +
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
    
    pdf(sprintf("%s/%s_plots.pdf", out_dir, out_name))
    print(heatmap_plot)
    print(silhouette_score_plot)
    print(silhouette_umap_plot)
    print(similarity_plot)
    dev.off()
}

get_clustering_subdir_basename <- function(dim_param, knn, rez, clust_method, out_name) {
    list2return <- list()
    if (!is.na(dim_param) & !is.na(knn)) {
        if (dim_param == "estimation") {
            dim_subdir <- "00-maxLikGlobalDimEst_p5"
            dim_nb_in_dir <- "estimp5"
        } else {
            dim_nb_in_dir <- dim_param
            dim_subdir <- sprintf("10-%d_PCs", dim_param)
        }
        knn_subdir <- sprintf("k%d", knn)
        clust_out_basename <- sprintf("%s-Clustering-%sPCs_%s", out_name, dim_nb_in_dir, knn_subdir) # here out_name is set to sample name
    } else {
        dim_subdir <- dim_nb_in_dir <- knn_subdir <- NA
        clust_out_basename <- out_name # here out_name is set to sample name, dim_param and knn
    }
    
    if (!is.na(rez) & !is.na(clust_method)) {
        rez_in_dir <- sprintf("resolution%s", sub("[.]", "_", rez))

        if (clust_method == 1) {
            clust_method_dir_index <- "00"
            clust_method_name <- "Louvain"
        } else {
           if (clust_method == 2) {
               clust_method_dir_index <- "10"
               clust_method_name <- "Louvain_multilevel"
           } else {
               if (clust_method == 3) {
                   clust_method_dir_index <- "20"
                   clust_method_name <- "SLM"
               } else {
                   if (clust_method == 4) {
                       clust_method_dir_index <- "30"
                       clust_method_name <- "Leiden"
                   }
               }
           } 
        }
        clust_method_subdir <- sprintf("%s-%s", clust_method_dir_index, clust_method_name)
        clust_out_basename <- sprintf("%s_%s_%s", clust_out_basename, rez_in_dir, clust_method_name)
    } else {
        rez_in_dir <- clust_method_subdir <- NA
    }
    
    list2return[["dimensions_subdir"]] <- dim_subdir
    list2return[["knn_subdir"]] <- knn_subdir
    list2return[["resolution_subdir"]] <- rez_in_dir
    list2return[["clustering_subdir"]] <- clust_method_subdir
    list2return[["basename"]] <- clust_out_basename
    
    return(list2return)
}

clustering <- function(seurat_obj, dim_param, knn, rez, method, analysis, nb_spots, iter, out_dir, out_name) {
    #clust_out <- get_clustering_subdir_basename(dim_param, knn, rez, method)
    #clust_out_basename <- clust_out[["basename"]]
    #clust_out_dir <- sprintf("%s/%s/%s/%s/%s", out_dir, clust_out[["dimensions_subdir"]], clust_out[["knn_subdir"]], clust_out[["resolution_subdir"]], clust_out[["clustering_subdir"]])
    #clust_out_name <- sprintf("%s-%s", out_name, clust_out_basename)
    
    if (!is.na(dim_param) & !is.na(knn)) {
        clust_out <- get_clustering_subdir_basename(dim_param, knn, NA, NA, out_name)
        
        if ("SCT" %in% names(seurat_obj@assays)) {
            seurat_obj_assay="SCT"
        } else {
            # NormalizeData() function
            seurat_obj_assay="RNA"
        }
        seurat_obj <- RunPCA(seurat_obj, assay=seurat_obj_assay, verbose=FALSE)

        if (dim_param == "estimation") {
            # estimate intrinsic dimension
            dim_estimate <- maxLikGlobalDimEst(seurat_obj@reductions$pca@cell.embeddings, k=10)
            dim_nb <- round(dim_estimate$dim.est)+5
            df2ggplot <- data.frame(PC=1:length(seurat_obj@reductions$pca), stdev=seurat_obj@reductions$pca@stdev)
            dim_out_dir <- sprintf("%s/%s", out_dir, clust_out[["dimensions_subdir"]])
            if (! dir.exists(dim_out_dir)) {
                dir.create(dim_out_dir, recursive=TRUE, mode="0775")
            }
            write.csv(df2ggplot, file=sprintf("%s/%s-PCs_stdev.csv", dim_out_dir, out_name), quote=FALSE, row.names=TRUE)
            pdf(sprintf("%s/%s-PCs_stdev.pdf", dim_out_dir, out_name))
            p <- ggplot(df2ggplot, aes(PC, stdev)) +
                geom_point() +
                geom_point() +
                geom_vline(xintercept = round(dim_estimate$dim.est), color = "blue") +
                geom_vline(xintercept = dim_nb, color = "red") +
                theme_bw() +
                labs(x = "Principal components", y = "Standard deviation")
            print(p)
            dev.off()
        } else {
            dim_nb <- dim_param
        }
        
        seurat_obj <- FindNeighbors(seurat_obj, reduction="pca", dims=1:dim_nb, k.param=knn)
    } else {
        dim_knn_out_dir <- out_dir
    }
    
    if (!is.na(rez) & !is.na(method)) {
        # get numbers of dimensions and nearest neighbors
        if ("SCT" %in% names (seurat_obj@assays)) {
            dim_nb <- seurat_obj@commands$FindNeighbors.SCT.pca@params$dims[length(seurat_obj@commands$FindNeighbors.SCT.pca@params$dims)]
            knn <- seurat_obj@commands$FindNeighbors.SCT.pca@params$k.param
        } else {
            dim_nb <- seurat_obj@commands$FindNeighbors.RNA.pca@params$dims[length(seurat_obj@commands$FindNeighbors.RNA.pca@params$dims)]
            knn <- seurat_obj@commands$FindNeighbors.RNA.pca@params$k.param
        }
        
        # clustering
        seurat_obj <- FindClusters(seurat_obj, resolution=rez, algorithm=method, verbose=FALSE)
        seurat_obj <- RunUMAP(seurat_obj, reduction="pca", dims=1:dim_nb)
        SCT_cluster_plot1 <- DimPlot(seurat_obj, reduction="umap", label=TRUE)
        SCT_cluster_plot2 <- SpatialDimPlot(seurat_obj, label=TRUE, label.size=3)
        pdf(sprintf("%s/%s.pdf", out_dir, out_name))
        print(SCT_cluster_plot1 + SCT_cluster_plot2 + plot_annotation(title="Normalized counts"))
        print(SCT_cluster_plot1)
        print(SCT_cluster_plot2)
        dev.off()
        
        if (analysis) {
            clustering_analysis(seurat_obj, dim_nb, knn, rez, nb_spots, iter, out_dir, sprintf("%s_analysis", out_name))
        }
    }
    
    return(seurat_obj)
}

marker_genes <- function(seurat_obj, test, logFC, pct, out_dir, out_name) {
    library(dplyr)
    # Finding differentially expressed features (cluster biomarkers): https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
    ## find markers for every cluster compared to all remaining spots, report only the positive ones
    markers <- FindAllMarkers(seurat_obj, test.use=test, logfc.threshold=logFC, min.pct=pct, only.pos=TRUE)
    write.csv(markers, file=sprintf("%s/%s.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
    ### visualizations
    pdf(sprintf("%s/%s.pdf", out_dir, out_name))
    #### expression heatmap
    markers <- droplevels(markers) # removes clusters without markers to avoid errors
    top10_markers <- markers %>%
        group_by(cluster) %>%
        top_n(n=10, wt=avg_log2FC)
    DoHeatmap(seurat_obj, features=top10_markers$gene) + NoLegend()
    #### get top6 markers for each cluster
    top6_markers <- markers %>%
        group_by(cluster) %>%
        slice_max(n=6, order_by=avg_log2FC)
    for (one_cluster in levels(markers$cluster)) {
        one_cluster_markers <- top6_markers[which(top6_markers$cluster==one_cluster), "gene", drop=TRUE]
        #### violin plot
        print(VlnPlot(seurat_obj, features=one_cluster_markers) + plot_annotation(title=sprintf("Cluster: %s", one_cluster), theme=theme(plot.title=element_text(size=16))))
        #### UMAP reduced dimensions
        print(FeaturePlot(seurat_obj, features=one_cluster_markers) + plot_annotation(title=sprintf("Cluster: %s", one_cluster), theme=theme(plot.title=element_text(size=16))))
        #### spatial coordinates
        print(SpatialFeaturePlot(seurat_obj, features=one_cluster_markers, alpha=c(0.1,1)) + plot_annotation(title=sprintf("Cluster: %s", one_cluster), theme=theme(plot.title=element_text(size=14))))
    }
    dev.off()
    
    return(seurat_obj)
}

Seurat_pipeline <- function(spe_seurat, libsize, exp_genes, space_ranger_out_dir, norm_method, features_method, features_nb, features_resvar, dim_param, knn, clust_resolution, clust_method, clust_analysis, nb_spots, iter, test, logFC, pct, out_dir, out_name) {
    # parameters
    ## spe_seurat: either a SpatialExperiment object for QC() function or a Seurat for all the other steps
    ## if libsize and exp_genes are set to NA, then do not perform QC
    ## if norm_method, features_method and (features_nb and features_resvar) are set to NA, then do not perform normalization and feature selection
    ## if dim_param, knn, clust_method and clust_analysis parameters are set to NA, then do not perform clustering
    ## if test, logFC and pct are set to NA, then do not perform marker genes
    
    # QC
    if (!is.na(libsize) & !is.na(exp_genes)) {
        qc_out_basename <- sprintf("QC-libsize%d_exp%d", libsize, exp_genes)
        qc_out_dir <- sprintf("%s/00-%s", out_dir, qc_out_basename)
        if (! dir.exists(qc_out_dir)) {
            dir.create(qc_out_dir, recursive=TRUE, mode="0775")
        }
        qc_out_basename <- sprintf("%s-%s", out_name, qc_out_basename)
        spe_seurat <- QC(spe_seurat, libsize, exp_genes, space_ranger_out_dir, qc_out_dir, qc_out_basename)
    } else {
        qc_out_dir <- out_dir
    }
    
    # normalization and feature selection
    if (!is.na(norm_method) & !is.na(features_method) & !(is.na(features_nb) & is.na(features_resvar))) {
        ## get output subdirectories and basenames
        ### normalization
        norm_out <- get_normalization_subdir_basename(normalization_method)
        norm_out_subdir <- norm_out[["subdir"]]
        norm_out_basename <- norm_out[["basename"]]
        norm_out_dir <- sprintf("%s/%s", qc_out_dir, norm_out_subdir)
        if (! dir.exists(norm_out_dir)) {
            dir.create(norm_out_dir, recursive=TRUE, mode="0775")
        }
        norm_out_name <- sprintf("%s-%s", out_name, norm_out_basename)
        ### feature selection
        features_out <- get_features_subdir_basename(features_method, features_nb, features_resvar)
        features_out_subdir <- sprintf("%s/%s", features_out[["method_subdir"]], features_out[["param_subdir"]])
        features_out_basename <- sprintf("%s_%s", features_out[["method_basename"]], features_out[["param_basename"]])
        features_out_dir <- sprintf("%s/%s", norm_out_dir, features_out_subdir)
        if (! dir.exists(features_out_dir)) {
            dir.create(features_out_dir, recursive=TRUE, mode="0775")
        }
        features_out_name <- sprintf("%s-%s", out_name, features_out_basename)

        spe_seurat <- normalization_feature_selection(spe_seurat, norm_method, features_method, features_nb, features_resvar, norm_out_dir, norm_out_name, features_out_dir, features_out_name)
    } else {
        features_out_dir <- out_dir
    }
    
    # clustering
    if ((!is.na(dim_param) & !is.na(knn)) | (!is.na(clust_resolution) & !is.na(clust_method) & !is.na(clust_analysis))) {
        if (!is.na(dim_param) & !is.na(knn)) {
            clust_out_dir <- features_out_dir
            clust_out <- get_clustering_subdir_basename(dim_param, knn, NA, NA, out_name)
            clust_out_basename <- clust_out[["basename"]]
            clust_out_name <- sprintf("%s-%s", out_name, clust_out_basename)
        } else {
            if (!is.na(clust_resolution) & !is.na(clust_method) & !is.na(clust_analysis)) {
                clust_out <- get_clustering_subdir_basename(NA, NA, clust_resolution, clust_method, out_name)
                clust_out_dir <- sprintf("%s/%s/%s", out_dir, clust_out[["resolution_subdir"]], clust_out[["clustering_subdir"]])
                if (! dir.exists(clust_out_dir)) {
                    dir.create(clust_out_dir, recursive=TRUE, mode="0775")
                }
                clust_out_name <- clust_out[["basename"]]
            }
        }
        spe_seurat <- clustering(spe_seurat, dim_param, knn, clust_resolution, clust_method, clust_analysis, nb_spots, iter, clust_out_dir, clust_out_name)
    } else {
        clust_out_dir <- out_dir
    }
    
    # marker genes
    if (!is.na(test) & !is.na(logFC) & !is.na(pct)) {
        if (test == "wilcox") {
            mark_test_dir_index <- "00"

        } else {
            if (mark_test == "bimod") {
                mark_test_dir_index <- "10"
            }
        }
        mark_test_out_basename <- sprintf("Markers-%s", test)
        mark_param_out_basename <- sprintf("logFC%s_minpct%s", sub("[.]", "_", logFC), sub("[.]", "_", pct))
        mark_out_basename <- sprintf("%s_%s", mark_test_out_basename, mark_param_out_basename)
        mark_out_dir <- sprintf("%s/%s-%s/%s", clust_out_dir, mark_test_dir_index, mark_test_out_basename, mark_param_out_basename)
        if (! dir.exists(mark_out_dir)) {
            dir.create(mark_out_dir, recursive=TRUE, mode="0775")
        }
        mark_out_basename <- sprintf("%s-%s", out_name, mark_out_basename)

        spe_seurat <- marker_genes(spe_seurat, test, logFC, pct, mark_out_dir, mark_out_basename)
    }
    
    return(spe_seurat)
}

top_counts <- function(df, col, rank_col, n) {
    parameter_values <- levels(df[[col]])
    count_vector <- c()
    for (one_parameter_value in parameter_values) {
        nb_top <- sum(df[which(df[[col]]==one_parameter_value), rank_col] <= n)
        count_vector <- c(count_vector, nb_top)
    }
    df2return <- data.frame(value=parameter_values, count=count_vector)
    return(df2return)
}

top_counts_barplot <- function(df, x_label, fill_label, n) {
    p <- ggplot(df, aes(x=value, y=count, fill=value)) +
        geom_bar(stat="identity") +
        labs(x=x_label, y=sprintf("#top%d", n), fill=fill_label) +
        theme_bw() +
        theme(legend.position="top") +
        theme(axis.title.x=element_text(size=10), axis.title.y=element_text(size=10), legend.title=element_text(size=8)) +
        theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), legend.text=element_text(size=6))
        theme(panel.border=element_rect(color="grey50"))
    return(p)
}

mean_silhouette_plots <- function(df, col_silhouette, col_rank, out_dir, out_name) {
    
    pdf(sprintf("%s/%s.pdf", out_dir, out_name))
    p <- ggplot(df, aes(x=param_setting, y=.data[[col_silhouette]])) +
        geom_line(group=1) +
        geom_point() +
        labs(x="Parameter setting", y="Mean silhouette score") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=90, hjust=1, colour=ifelse(df$rank < 10, "red", "black"), size=6)) +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
    
    p <- ggplot(df, aes(x=normalization, y=.data[[col_silhouette]], fill=features_resvar)) +
        geom_boxplot() +
        facet_wrap(~dimensions) +
        labs(x="Normalization method", y="Mean silhouette score", fill="Residual variance threshold\nfor feature selection") +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50"))
    print(p) 
    
    for (normalization_method in normalization_method_vector) {
        p <- ggplot(df[which(df$normalization==normalization_method),], aes(x=features_resvar, y=.data[[col_silhouette]], fill=dimensions)) +
            geom_boxplot() +
            facet_wrap(~knn) +
            labs(title=sprintf("Normalization: %s", normalization_method), x="Residual variance threshold for feature selection", y="Mean silhouette score", fill="Number of\ndimensions") +
            theme_bw() +
            theme(panel.border=element_rect(color="grey50"))
        print(p)

        for (feature_residual_variance_threshold in feature_residual_variance_threshold_vector) {
            p <- ggplot(df[which(df$normalization=="SCTransform" & df$features_resvar==feature_residual_variance_threshold),], aes(x=dimensions, y=.data[[col_silhouette]], fill=knn)) +
                geom_boxplot() +
                facet_wrap(~resolution) +
                labs(title=sprintf("Normalization: %s\nresidual variance threshold for feature selection: %.1f", normalization_method, feature_residual_variance_threshold), x="Number of dimensions", y="Mean silhouette score", fill="Number of nearest\nneighbors") +
                theme_bw() +
                theme(panel.border=element_rect(color="grey50"))
            print(p)

            for (dimensions_nb_param in dimensions_nb_param_vector) {
                if (dimensions_nb_param == 0) {
                    dimensions_nb_param <- "estimation"
                }
                p <- ggplot(df[which(df$normalization=="SCTransform" & df$features_resvar==feature_residual_variance_threshold & df$dimensions==dimensions_nb_param),], aes(x=knn, y=.data[[col_silhouette]], fill=resolution)) +
                    geom_bar(stat="identity", position="dodge") +
                    facet_wrap(~clustering_method) +
                    labs(title=sprintf("Normalization: %s\nresidual variance threshold for feature selection: %.1f\nnumber of dimensions: %s", normalization_method, feature_residual_variance_threshold, dimensions_nb_param), x="Number of nearest\nneighbors", y="Mean silhouette score", fill="Resolution") +
                    theme_bw() +
                    theme(panel.border=element_rect(color="grey50"))
                print(p)
            }
        }
    }
    dev.off()
}

mean_silhouette_plots_top <- function(df, col_rank, out_dir, out_name) {
    
    pdf(sprintf("%s/%s.pdf", out_dir, out_name))
    # number of times in the top parameter settings
    top_df <- data.frame()
    for (top in c(3, 10, 20, 50)) {
        df2ggplot <- top_counts(df, "normalization", col_rank, top)
        normalization_plot <- top_counts_barplot(df2ggplot, "Normalization method", "Normalization\nmethod", top)
        top_df <- rbind(top_df, cbind(parameter=rep("normalization method", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        df2ggplot <- top_counts(df, "features_resvar", col_rank, top)
        features_resvar_plot <- top_counts_barplot(df2ggplot, "Feature residual variance threshold", "Feature residual\nvariance threshold", top)
        top_df <- rbind(top_df, cbind(parameter=rep("feature residual variance threshold", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        df2ggplot <- top_counts(df, "dimensions", col_rank, top)
        dimensions_plot <- top_counts_barplot(df2ggplot, "Number of dimensions", "Number of\ndimensions", top)
        top_df <- rbind(top_df, cbind(parameter=rep("number of dimensions", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        df2ggplot <- top_counts(df, "knn", col_rank, top)
        knn_plot <- top_counts_barplot(df2ggplot, "Number of nearest neighbors", "Number of nearest\nneighbors", top)
        top_df <- rbind(top_df, cbind(parameter=rep("number of nearest neighbors", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        df2ggplot <- top_counts(df, "resolution", col_rank, top)
        resolution_plot <- top_counts_barplot(df2ggplot, "Clustering resolution", "Clustering\nresolution", top)
        top_df <- rbind(top_df, cbind(parameter=rep("resolution", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        df2ggplot <- top_counts(df, "clustering_method", col_rank, top)
        clustering_method_plot <- top_counts_barplot(df2ggplot, "Clustering method", "Clustering\nmethod", top)
        top_df <- rbind(top_df, cbind(parameter=rep("clustering method", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        multiplot <- ggdraw() +
          draw_plot(normalization_plot, 0, 0.5, 0.33, 0.5) +
          draw_plot(features_resvar_plot, 0.33, 0.5, 0.33, 0.5) +
          draw_plot(dimensions_plot, 0.66, 0.5, 0.33, 0.5) +
          draw_plot(knn_plot, 0, 0, 0.33, 0.5) +
          draw_plot(resolution_plot, 0.33, 0, 0.33, 0.5) +
          draw_plot(clustering_method_plot, 0.66, 0, 0.33, 0.5)
        print(multiplot)
    }
    dev.off()
    write.csv(top_df, file=sprintf("%s/%s.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
}
