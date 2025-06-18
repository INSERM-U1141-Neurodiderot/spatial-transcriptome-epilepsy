#####################
# load count matrix #
#####################

ST_Pipeline2SpatialExperiment <- function(sample, st_dir, sr_dir, metadata_df) {
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

Space_Ranger2SpatialExperiment <- function(sample, sr_dir, metadata_df) {
    # Space Ranger output directory for tissue positions list and image
    space_ranger_out_dir <- sprintf("%s/%s/10-Pipeline/outs", sr_dir, sample)
    
    # Space Ranger filtered feature count matrix: contains only tissue-associated barcodes
    library(Seurat)
    space_ranger_matrix_file <- sprintf("%s/filtered_feature_bc_matrix.h5", space_ranger_out_dir)
    space_ranger_matrix <- Read10X_h5(space_ranger_matrix_file)
    ## detach Seurat and SeuratObject packages to avoid conflict with SpatialExperiment package: Found more than one class "SpatialImage" in cache
    detach("package:Seurat", unload=TRUE)
    detach("package:SeuratObject", unload=TRUE)
    
    # meta data
    ## spatial coordinates
    spatial_coordinates_file <- file.path(space_ranger_out_dir, "spatial", "tissue_positions_list.csv")
    spatial_coordinates_df <- read.csv(spatial_coordinates_file, header=FALSE, quote="")
    colnames(spatial_coordinates_df) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    rownames(spatial_coordinates_df) <- spatial_coordinates_df$barcode
    ## add lame, zone, condition and time
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
    ## gene names
    gene_data <- data.frame(gene_name=rownames(space_ranger_matrix))
    ## image data
    img <- readImgData(imageSources=file.path(space_ranger_out_dir, "spatial", "tissue_lowres_image.png"), scaleFactors=file.path(space_ranger_out_dir, "spatial", "scalefactors_json.json"), sample_id=sample)
    
    # create SpatialExperiment object
    spe <- SpatialExperiment(assay=list(counts=space_ranger_matrix), rowData=gene_data, colData=spatial_coordinates_df[colnames(space_ranger_matrix),], imgData=img, spatialDataNames=c("barcode", "in_tissue", "array_row", "array_col"), spatialCoordsNames=c("pxl_col_in_fullres", "pxl_row_in_fullres"), sample_id=sample)
    
    # remove non-expressed genes: genes with 0 count for all spots
    non_expressed <- apply(spe@assays@data$counts, 1, sum) == 0
    spe <- spe[!non_expressed, ]
    
    return(spe)
}

#load_data <- function(sample, st_dir, sr_dir, add_image) {
#    # expression matrices
#    ## ST Pipeline
#    st_pipeline_matrix_file <- sprintf("%s/%s/10-Pipeline/%s_stdata.tsv", st_dir, sample, sample)
#    st_pipeline_matrix <- read.table(st_pipeline_matrix_file, sep="\t", header=TRUE, quote="", row.names=1)
#    # Space Ranger output directory for tissue positions list and image
#    space_ranger_output_dir <- sprintf("%s/%s/10-Pipeline/outs", sr_dir, sample)
#
#    # build SpatialExperiment object
#    spatial_coordinates_file <- file.path(space_ranger_output_dir, "spatial", "tissue_positions_list.csv")
#    spatial_coordinates_df <- read.csv(spatial_coordinates_file, header=FALSE, quote="")
#    colnames(spatial_coordinates_df) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
#    rownames(spatial_coordinates_df) <- spatial_coordinates_df$barcode
#    nb_spots <- dim(spatial_coordinates_df)[1]
#    ## add lame, zone and condition
#    sample_strsplit <- unlist(strsplit(sample, "_"))
#    zone <- sample_strsplit[1]
#    lame <- sample_strsplit[2]
#    condition <- ifelse(zone=="A" | zone=="C", "SE", "CTRL")
#    spatial_coordinates_df <- cbind(spatial_coordinates_df, lame=rep(lame, nb_spots), zone=rep(zone, nb_spots), condition=rep(condition, nb_spots))
#    spatial_coordinates_df$zone <- as.factor(spatial_coordinates_df$zone)
#    spatial_coordinates_df$lame <- as.factor(spatial_coordinates_df$lame)
#    spatial_coordinates_df$condition <- as.factor(spatial_coordinates_df$condition)
#    ## matrix: convert coordinates to barcodes
#    barcodes <- unlist(lapply(rownames(st_pipeline_matrix), function(x, coord2barcode=spatial_coordinates_df) {
#        coordinates <- unlist(strsplit(x, "x"))
#        return(coord2barcode[which(coord2barcode$array_col==as.integer(coordinates[1])-1 & coord2barcode$array_row==as.integer(coordinates[2])-1), "barcode"])
#    }))
#    rownames(st_pipeline_matrix) <- barcodes
#    
#    # create Seurat object
#    seurat_obj <- CreateSeuratObject(counts=t(st_pipeline_matrix), project=sample, meta.data=spatial_coordinates_df)
#    ## add image
#    if (add_image) {
#        img <- Read10X_Image(image.dir=file.path(space_ranger_output_dir, "spatial"))
#        img@assay <- Assays(seurat_obj) # set assay identical to Seurat object
#        Key(img) <- sprintf("%s_", gsub("_", "", sample)) # add key for SpatialFeaturePlot() and SpatialPlot() functions: Keys should be one or more alphanumeric characters followed by an underscore
#        img_list <- list(img) # the image must be stored in a list
#        names(img_list) <- sample # change image name in the list: used as title by SpatialFeaturePlot() function for example
#        seurat_obj@images <- img_list
#    }
#    ## keep spot over tissue
#    seurat_obj <- seurat_obj[, seurat_obj$in_tissue==1]
#    
#    return(seurat_obj)
#}

gtf_gene_regex <- function(GRanges_obj, colname, regex, out_dir) {
    # convert GRanges object to DFrame object
    dframe_obj <- mcols(GRanges_obj)
    dframe_obj_genes <- dframe_obj[which(grepl(regex, dframe_obj[[colname]])),]
    regex_sub <- gsub("[(]|[)]", "", regex)
    regex_sub <- gsub("[|]| ", "_", regex_sub)
    out_name <- sprintf("%s_genes", regex_sub)
    write.table(dframe_obj_genes, file=sprintf("%s/%s.tsv", out_dir, out_name), sep="\t", quote=FALSE, row.names=FALSE)
    return(dframe_obj_genes$gene_id)
}

######
# QC #
######
add_qc_stats <- function(sp_exp, mito_genes, hb_genes, col_prefix, out_dir, out_name) {
    # identify mitochondrial genes
    #is_mito <- grepl("(^MT(-|.))|(^mt(-|.))", rowData(sp_exp)$gene_name)
    is_mito <- rowData(sp_exp)$gene_name %in% mito_genes
    # identify hemoglobin genes
    #is_hb <- grepl("^Hb(a|b|q|s|z)", rowData(sp_exp)$gene_name)
    is_hb <- rowData(sp_exp)$gene_name %in% hb_genes
    # calculate per-spot QC metrics and store in colData
    sp_exp <- addPerCellQC(sp_exp, subsets=list(mito=is_mito, hb=is_hb))
    col2rename <- c("sum", "detected", "subsets_mito_sum", "subsets_mito_detected", "subsets_mito_percent", "subsets_hb_sum", "subsets_hb_detected", "subsets_hb_percent")
    colnames(colData(sp_exp))[colnames(colData(sp_exp)) %in% col2rename] <- sprintf("%s_%s", col_prefix, col2rename)
    write.csv(colData(sp_exp), file=sprintf("%s/%s_%s.csv", out_dir, sp_exp@int_metadata$imgData$sample_id, out_name), quote=FALSE, row.names=TRUE)
    return(sp_exp)
}

plot_density <- function(data2plot, plot_title, plot_x_lab) {
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

qc_stat_plots <- function(sp_exp, col_prefix, out_dir, out_name) {
    pdf(sprintf("%s/%s_%s.pdf", out_dir, sp_exp@int_metadata$imgData$sample_id, out_name))
    # sum of counts
    colname <- sprintf("%s_sum", col_prefix)
    print(plotQC(sp_exp, type="scatter", metric_x=sprintf("%s_sum", col_prefix), metric_y=sprintf("%s_detected", col_prefix)))
    ## distribution
    spot_data <- colData(sp_exp)[[colname]]
    percentiles <- plot_density(spot_data, "Detected genes per spot", "Sum of counts")    
    ## spatial plots
    print(plotSpots(sp_exp, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate=colname, palette=c("white", "black"), size=2))
    print(plotVisium(sp_exp, fill=colname, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    
    # number of expressed features
    colname <- sprintf("%s_detected", col_prefix)
    ## distribution
    spot_data <- colData(sp_exp)[[colname]]
    percentiles <- plot_density(spot_data, "Sum of counts per spot", "Detected genes")    
    ## spatial plots
    print(plotSpots(sp_exp, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate=colname, palette=c("white", "black"), size=2))
    print(plotVisium(sp_exp, fill=colname, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    
    # percentage of mitochondrial gene counts
    colname <- sprintf("%s_subsets_mito_percent", col_prefix)
    print(plotQC(sp_exp, type="scatter", metric_x=sprintf("%s_sum", col_prefix), metric_y=colname))
    ## distribution
    spot_data <- colData(sp_exp)[[colname]]
    percentiles <- plot_density(spot_data, "Percentage of mitochondrial gene counts per spot", "Percentage of mitochondrial gene counts")    
    ## spatial plots
    print(plotSpots(sp_exp, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate=colname, palette=c("white", "black"), size=2))
    print(plotVisium(sp_exp, fill=colname, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    
    # percentage of hemoglobin gene counts
    colname <- sprintf("%s_subsets_hb_percent", col_prefix)
    print(plotQC(sp_exp, type="scatter", metric_x=sprintf("%s_sum", col_prefix), metric_y=colname))
    ## distribution
    spot_data <- colData(sp_exp)[[colname]]
    percentiles <- plot_density(spot_data, "Percentage of hemoglobin gene counts per spot", "Percentage of hemoglobin gene counts")    
    ## spatial plots
    print(plotSpots(sp_exp, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate=colname, palette=c("white", "black"), size=2))
    print(plotVisium(sp_exp, fill=colname, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    dev.off()
}

gene_spatial_mean_exp_plots <- function(seurat_obj, gene, assay2plot, plots_per_page) {
    nb_images <- length(names(seurat_obj@images))
    # spatial feature plot
    plot_title <- sprintf("%s\nmean: %.2f, expressed: %.2f%%, rank: %d", gene, seurat_obj@assays[[assay2plot]]@meta.features[gene, "mean"], seurat_obj@assays[[assay2plot]]@meta.features[gene, "detected"], floor(seurat_obj@assays[[assay2plot]]@meta.features[gene, "mean_rank"]))
    if (nb_images == 1) {
        p <- SpatialFeaturePlot(seurat_obj, features=gene, alpha=c(0.5,1)) + plot_annotation(title=plot_title, theme=theme(plot.title=element_text(size=14)))
        print(p)
    } else {
        if (nb_images > 4) {
            # one page per timepoint
            timepoints <- levels(seurat_obj@meta.data$time)
            for (one_timepoint in timepoints) {
                one_timepoint_metadata_df <- seurat_obj@meta.data[which(seurat_obj@meta.data$time==one_timepoint),]
                one_timepoint_metadata_df <- droplevels(one_timepoint_metadata_df)
                sections <- levels(one_timepoint_metadata_df$orig.ident)
                p <- SpatialFeaturePlot(seurat_obj, features=gene, alpha=c(0.5,1), images=sections, ncol=2) + plot_annotation(title=sprintf("%s\ntimepoint: %s", plot_title, one_timepoint), theme=theme(plot.title=element_text(size=14)))
                print(p)
            }
        } else {
            for (i in 1:(nb_images/plots_per_page)) {
                p <- SpatialFeaturePlot(seurat_obj, features=gene, alpha=c(0.5,1), images=names(seurat_obj@images)[((i-1)*plots_per_page+1):(i*plots_per_page)], ncol=2) + plot_annotation(title=plot_title, theme=theme(plot.title=element_text(size=14)))
                print(p)
            }
        }
    }

    # mean expression density plot
    gene_log_exp <- log(seurat_obj@assays[[assay2plot]]@counts[gene,], 2)
    gene_percentiles <- plot_density(gene_log_exp, sprintf("%s expression", gene), "log2(Expression)")

    return(gene_percentiles)
}

qc_spatial_mean_exp_plots <- function(seurat_obj, gene_vec, mean_order, plots_per_page, out_dir, out_name) {
    sample <- seurat_obj@project.name
    active_assay <- seurat_obj@active.assay
    gene_vec <- gene_vec[gene_vec %in% rownames(seurat_obj@assays[[active_assay]]@meta.features)]
    gene_stats_df <- seurat_obj@assays[[active_assay]]@meta.features[gene_vec,]
    if (mean_order) {
        gene_stats_df <- gene_stats_df[order(gene_stats_df$mean, decreasing=TRUE),]
        if (length(gene_vec) > 20) {
            # get top 20 most expressed genes
            genes2plot <- rownames(gene_stats_df[1:20,])
        } else {
            genes2plot <- rownames(gene_stats_df)
        }
    } else {
        genes2plot <- gene_vec
    }
    genes2plot_df <- data.frame()
    pdf(sprintf("%s/%s_%s.pdf", out_dir, sample, out_name))
    for (one_gene in genes2plot) {
        one_gene_mean_exp_percentiles <- gene_spatial_mean_exp_plots(seurat_obj, one_gene, active_assay, plots_per_page)
        genes2plot_df <- rbind(genes2plot_df, c(one_gene, one_gene_mean_exp_percentiles))
    }
    dev.off()
    colnames(genes2plot_df) <- c("gene", names(one_gene_mean_exp_percentiles))
    write.csv(genes2plot_df, file=sprintf("%s/%s_%s.csv", out_dir, sample, out_name), quote=FALSE, row.names=FALSE)
    
    return(cbind(sample=rep(sample, dim(genes2plot_df)[1]), genes2plot_df))
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

qc_filter <- function(sp_exp, raw_col_prefix, qc_col_prefix, sum_thres, dectected_thres, mito_pct_thres, hb_pct_thres, out_dir, out_name) {
    pdf(sprintf("%s/%s_%s.pdf", out_dir, sp_exp@int_metadata$imgData$sample_id, out_name))
    # identify spots to filter according different statistics
    ## sum of counts
    colname <- "sum"
    raw_colname <- sprintf("%s_%s", raw_col_prefix, colname)
    threshold <- sum_thres
    sum_qc_colname <- sprintf("%s_%s_%d", qc_col_prefix, colname, threshold)
    sp_exp <- qc_threshold(sp_exp, raw_colname, threshold, FALSE, sum_qc_colname)
    qc_threshold_plot(sp_exp, sum_qc_colname)
    ## number of expressed features
    colname <- "detected"
    raw_colname <- sprintf("%s_%s", raw_col_prefix, colname)
    threshold <- dectected_thres
    detected_qc_colname <- sprintf("%s_%s_%d", qc_col_prefix, colname, threshold)
    sp_exp <- qc_threshold(sp_exp, raw_colname, threshold, FALSE, detected_qc_colname)
    qc_threshold_plot(sp_exp, detected_qc_colname)
    ## percentage of mitochondrial gene counts
    colname <- "subsets_mito_percent"
    raw_colname <- sprintf("%s_%s", raw_col_prefix, colname)
    threshold <- mito_pct_thres
    mito_qc_colname <- sprintf("%s_%s_%d", qc_col_prefix, colname, threshold)
    sp_exp <- qc_threshold(sp_exp, raw_colname, threshold, TRUE, mito_qc_colname)
    qc_threshold_plot(sp_exp, mito_qc_colname)
    ## percentage of hemoglobin gene counts
    colname <- "subsets_hb_percent"
    raw_colname <- sprintf("%s_%s", raw_col_prefix, colname)
    threshold <- hb_pct_thres
    hb_qc_colname <- sprintf("%s_%s_%d", qc_col_prefix, colname, threshold)
    sp_exp <- qc_threshold(sp_exp, raw_colname, threshold, TRUE, hb_qc_colname)
    qc_threshold_plot(sp_exp, hb_qc_colname)
    ## all statistics
    discard <- colData(sp_exp)[[sum_qc_colname]] | colData(sp_exp)[[detected_qc_colname]] | colData(sp_exp)[[mito_qc_colname]] | colData(sp_exp)[[hb_qc_colname]]
    colData(sp_exp)$discard <- discard
    qc_threshold_plot(sp_exp, "discard")
    dev.off()
    
    # remove spots to filter
    sp_exp <- sp_exp[,! colData(sp_exp)$discard]
    return(sp_exp)
}

addFeatureStats2SpatialExperiment <- function(sp_exp) {
    # if addFeatureStats2SpatialExperiment() has already been performed, then remove 'mean', 'detected' and 'mean_rank' columns from rowData(sp_exp)
    rowData(sp_exp) <- rowData(sp_exp)[, ! colnames(rowData(sp_exp)) %in% c("mean", "detected", "mean_rank")]
    # here rowData(sp_exp) has only one column which name is 'value': change it to 'gene_name'
    colnames(rowData(sp_exp)) <- "gene_name"
    sp_exp <- addPerFeatureQC(sp_exp)
    rowData(sp_exp)$mean_rank <- rank(-rowData(sp_exp)$mean)
    return(sp_exp)
}

addFeatureStats2Seurat <- function(seurat_obj, assay2use, slot2use) {
    # if addFeatureStats2Seurat() has already been performed, then rename 'mean', 'detected' and 'mean_rank' columns from seurat_obj@assays[[assay2use]]@meta.features by adding '.old' suffix
    for (one_colname in c("mean", "detected", "mean_rank")) {
        if (one_colname %in% colnames(seurat_obj@assays[[assay2use]]@meta.features)) {
            colnames(seurat_obj@assays[[assay2use]]@meta.features)[which(colnames(seurat_obj@assays[[assay2use]]@meta.features) == one_colname)] <- sprintf("%s.old", one_colname)
        }
    }
    mean_exp <- apply(slot(seurat_obj@assays[[assay2use]], slot2use), 1, mean)
    nb_spots <- dim(slot(seurat_obj@assays[[assay2use]], slot2use))[2]
    detected_spots <- apply(slot(seurat_obj@assays[[assay2use]], slot2use), 1, function(x, spots=nb_spots) { length(x[x > 0])/ spots * 100 })
    if (dim(seurat_obj@assays[[assay2use]]@meta.features)[2] > 0) {
        stats_df <- data.frame(gene_name=rownames(slot(seurat_obj@assays[[assay2use]], slot2use)), mean=mean_exp, detected=detected_spots)
        gene_stats_df <- merge(seurat_obj@assays[[assay2use]]@meta.features, stats_df, by="gene_name")
        rownames(gene_stats_df) <- gene_stats_df$gene_name
        gene_stats_df <- gene_stats_df[rownames(slot(seurat_obj@assays[[assay2use]], slot2use)),]
    } else {
        gene_stats_df <- data.frame(gene_name=rownames(slot(seurat_obj@assays[[assay2use]], slot2use)), mean=mean_exp, detected=detected_spots)
    }
    gene_stats_df$mean_rank <- rank(-gene_stats_df$mean)
    seurat_obj@assays[[assay2use]]@meta.features <- gene_stats_df
    return(seurat_obj)
}

SpatialExperiment_Seurat_conversion <- function(sp_exp, sr_dir, add_image) {
    # Space Ranger output directory for tissue positions list and image
    sample <- sp_exp@int_metadata$imgData$sample_id
    space_ranger_out_dir <- sprintf("%s/%s/10-Pipeline/outs", sr_dir, sample)
    
    # meta data
    ## colnames(colData(sp_exp))[grepl("_percent", colnames(colData(sp_exp)))]
    metadata_df <- merge(as.data.frame(colData(sp_exp))[, c("barcode", "in_tissue", "array_row", "array_col", "lame", "zone", "condition", "time", colnames(colData(sp_exp))[grepl("_percent", colnames(colData(sp_exp)))])], as.data.frame(sp_exp@int_colData$spatialCoords), by.x="barcode", by.y="row.names")
    rownames(metadata_df) <- metadata_df$barcode
    
    # create Seurat object
    #library(Seurat)
    seurat_obj <- CreateSeuratObject(counts=counts(sp_exp), project=sample, meta.data=metadata_df)
    ## SpatialExperiment package is required for Seurat object creation, if it is loaded before creating Seurat object, then it is automatically loaded during Seurat object creation
    #detach("package:SpatialExperiment", unload=TRUE)
    ## add image
    if (add_image) {
        img <- Read10X_Image(image.dir=file.path(space_ranger_out_dir, "spatial"))
        img@assay <- Assays(seurat_obj) # set assay identical to Seurat object
        Key(img) <- sprintf("%s_", gsub("_", "", sample)) # add key for SpatialFeaturePlot() and SpatialPlot() functions: Keys should be one or more alphanumeric characters followed by an underscore
        img_list <- list(img) # the image must be stored in a list
        names(img_list) <- sample # change image name in the list: used as title by SpatialFeaturePlot() function for example
        seurat_obj@images <- img_list
    }
    ## remove discarded barcodes from image coordinates
    discarded_barcodes <- rownames(seurat_obj@images[[sample]]@coordinates)[! rownames(seurat_obj@images[[sample]]@coordinates) %in% rownames(colData(sp_exp))]
    seurat_obj@images[[sample]]@coordinates <- seurat_obj@images[[sample]]@coordinates[! rownames(seurat_obj@images[[sample]]@coordinates) %in% discarded_barcodes,]
    #spe_seurat@images$sample@assay <- "RNA"
    active_assay <- seurat_obj@active.assay
    seurat_obj@images[[sample]]@assay <- active_assay
    ## add feature data
    ### replace '_' by '-' in feature names: “Feature names cannot have underscores ('_'), replacing with dashes ('-')”
    feature_data_df <- as.data.frame(rowData(sp_exp))
    feature_data_df$gene_name <- gsub("_", "-", feature_data_df$gene_name)
    rownames(feature_data_df) <- feature_data_df$gene_name
    seurat_obj@assays[[active_assay]]@meta.features <- feature_data_df
    #test <- Load10X_Spatial(data.dir=space_ranger_out_dir)
    
    return(seurat_obj)
}

addPercentageFeatureSet2Seurat <- function(seurat_obj, mito_genes, hb_genes, assay2use, col_prefix) {
    seurat_obj <- PercentageFeatureSet(seurat_obj, pattern=NULL, features=mito_genes[mito_genes %in% rownames(seurat_obj@assays[[assay2use]]@counts)], col.name=sprintf("%s_mito_percent_%s", col_prefix, assay2use), assay2use)
    seurat_obj <- PercentageFeatureSet(seurat_obj, pattern=NULL, features=hb_genes[hb_genes %in% rownames(seurat_obj@assays[[assay2use]]@counts)], col.name=sprintf("%s_hb_percent_%s", col_prefix, assay2use), assay2use)
    return(seurat_obj)
}

seurat_feature_stat_plots <- function(seurat_obj, feature_name, suffix, plots_per_page, plot_title, plot_legend_title) {
    nb_images <- length(names(seurat_obj@images))
    plots2return <- list()
    violin_plot <- VlnPlot(seurat_obj, features=sprintf("%s_%s", feature_name, suffix), group.by="orig.ident", pt.size=0.1) + NoLegend()
    
    if (nb_images == 1) {
        spatial_feature_plot_list <- list(plot=SpatialFeaturePlot(seurat_obj, features=sprintf("%s_%s", feature_name, suffix)) + theme(legend.position="top"))
    } else {
        if (nb_images > 4) {
            spatial_feature_plot_list <- list()
            # one page per timepoint
            timepoints <- levels(seurat_obj@meta.data$time)
            for (one_timepoint in timepoints) {
                one_timepoint_metadata_df <- seurat_obj@meta.data[which(seurat_obj@meta.data$time==one_timepoint),]
                one_timepoint_metadata_df <- droplevels(one_timepoint_metadata_df)
                sections <- levels(one_timepoint_metadata_df$orig.ident)
                spatial_feature_plot_list[[one_timepoint]] <- SpatialFeaturePlot(seurat_obj, features=sprintf("%s_%s", feature_name, suffix), images=sections, ncol=2) + plot_annotation(title=sprintf("%s, timepoint: %s", plot_title, one_timepoint), theme=theme(plot.title=element_text(size=16))) + theme(legend.position="top")
            }
        } else {
            for (i in 1:(nb_sections/plots_per_page)) {
                spatial_feature_plot_list <- list(plot=SpatialFeaturePlot(seurat_obj, features=sprintf("%s_%s", feature_name, suffix), images=names(seurat_obj@images)[((i-1)*plots_per_page+1):(i*plots_per_page)], ncol=2) + plot_annotation(title=plot_title, theme=theme(plot.title=element_text(size=16))) + theme(legend.position="top"))
            }
        }
    }
    
    #if (nb_images == 1) {
    #    spatial_feature_plot <- SpatialFeaturePlot(seurat_obj, features = sprintf("%s_%s", feature_name, suffix)) + theme(legend.position = "top")
    #} else {
    #    spatial_feature_plot <- SpatialFeaturePlot(seurat_obj, features = sprintf("%s_%s", feature_name, suffix), ncol=2) + theme(legend.position = "top")
    #}
    
    if (!is.na(plot_legend_title)) {
        violin_plot <- violin_plot + labs(title=plot_legend_title) + theme(plot.title=element_text(size=12))
        spatial_feature_plot_list <- lapply(spatial_feature_plot_list, function(x) {
            x <- x +
                scale_fill_gradientn(colors=rev(x=brewer.pal(n=11, name="Spectral")), name=plot_legend_title) + # need to recreate color palette to change legend title (see https://github.com/satijalab/seurat/blob/HEAD/R/visualization.R) 
                theme(legend.text=element_text(size=8))
            return(x)
        })
        
        #spatial_feature_plot <- spatial_feature_plot + 
        #    scale_fill_gradientn(colors=rev(x=brewer.pal(n=11, name="Spectral")), name=plot_legend_title) + # need to recreate color palette to change legend title (see https://github.com/satijalab/seurat/blob/HEAD/R/visualization.R) 
        #    theme(legend.text=element_text(size=8))
    }
    plots2return[["violin"]] <- violin_plot
    plots2return[["spatial"]] <- spatial_feature_plot_list
    return(plots2return)
}

seurat_all_feature_stat_plots <- function(seurat_obj, assay2use, col_prefix, plots_per_page, plot_title, out_dir, out_name) {
    sample <- seurat_obj@project.name
    nb_images <- length(names(seurat_obj@images))
    pdf(sprintf("%s/%s_%s_Seurat_features.pdf", out_dir, sample, out_name))
    # counts per spot
    plot_legend_title <- "count sum"
    raw_plots <- seurat_feature_stat_plots(seurat_obj, "nCount", assay2use, plots_per_page, plot_title, plot_legend_title)
    if (nb_images == 1) {
        plot_name <- names(raw_plots[["spatial"]])
        print(wrap_plots(raw_plots[["violin"]], raw_plots[["spatial"]][[plot_name]]))
    } else {
        print(raw_plots[["violin"]] + plot_annotation(title=plot_title))
        for (plot_name in names(raw_plots[["spatial"]])) {
            print(raw_plots[["spatial"]][[plot_name]])
        }
    }
    # number of expressed genes per spot
    plot_legend_title <- "detected genes"
    raw_plots <- seurat_feature_stat_plots(seurat_obj, "nFeature", assay2use, plots_per_page, plot_title, plot_legend_title)
    if (nb_images == 1) {
        plot_name <- names(raw_plots[["spatial"]])
        print(wrap_plots(raw_plots[["violin"]], raw_plots[["spatial"]][[plot_name]]))
    } else {
        print(raw_plots[["violin"]] + plot_annotation(title=plot_title))
        for (plot_name in names(raw_plots[["spatial"]])) {
            print(raw_plots[["spatial"]][[plot_name]])
        }
    }
    # mitochondrial gene count percentage per spot
    plot_legend_title <- "mitochondrial gene\ncount percentage"
    raw_plots <- seurat_feature_stat_plots(seurat_obj, sprintf("%s_mito_percent", col_prefix), assay2use, plots_per_page, plot_title, plot_legend_title)
    if (nb_images == 1) {
        plot_name <- names(raw_plots[["spatial"]])
        print(wrap_plots(raw_plots[["violin"]], raw_plots[["spatial"]][[plot_name]]))
    } else {
        print(raw_plots[["violin"]] + plot_annotation(title=plot_title))
        for (plot_name in names(raw_plots[["spatial"]])) {
            print(raw_plots[["spatial"]][[plot_name]])
        }
    }
    # hemoglobin gene count percentage per spot
    plot_legend_title <- "hemoglobin gene\ncount percentage"
    raw_plots <- seurat_feature_stat_plots(seurat_obj, sprintf("%s_hb_percent", col_prefix), assay2use, plots_per_page, plot_title, plot_legend_title)
    if (nb_images == 1) {
        plot_name <- names(raw_plots[["spatial"]])
        print(wrap_plots(raw_plots[["violin"]], raw_plots[["spatial"]][[plot_name]]))
    } else {
        print(raw_plots[["violin"]] + plot_annotation(title=plot_title))
        for (plot_name in names(raw_plots[["spatial"]])) {
            print(raw_plots[["spatial"]][[plot_name]])
        }
    }
    dev.off()
}

#################
# Normalization #
#################

seurat_list_normalization <- function(seurat_obj_list, method, vars2regress, feat_param, only_var, mito_genes, hb_genes, out_dir) {
    if (method == "SCTransform" | method == "SCTransform_v2") {
        vst_flavor <- switch(as.integer(method == "SCTransform_v2")+1, NULL, "v2")
        vars_to_regress <- switch(as.integer(is.na(vars2regress))+1, vars2regress, NULL)
        variable_features_n <- switch(as.integer(feat_param > 10)+1, NULL, feat_param) # if feat_param > 10, then it is a feature number, otherwise it is residual variance threshold
        variable_features_rv_th <- switch(as.integer(feat_param > 10)+1, feat_param, NULL) # if feat_param > 10, then it is a feature number, otherwise it is residual variance threshold
        print(sprintf("VST flavor: %s, variables to regress: %s, variable features n: %s, variable features rv th: %s", ifelse(is.null(vst_flavor), "NULL", vst_flavor), ifelse(is.null(vars_to_regress), "NULL", vars_to_regress), ifelse(is.null(variable_features_n), "NULL", variable_features_n), ifelse(is.null(variable_features_rv_th), "NULL", variable_features_rv_th)))
        seurat_objects_norm <- lapply(seurat_obj_list, FUN=SCTransform, vst.flavor=vst_flavor, vars.to.regress=vars_to_regress, variable.features.n=variable_features_n, variable.features.rv.th=variable_features_rv_th, return.only.var.genes=only_var, verbose=TRUE)
        
        method_dir_index <- ifelse(method == "SCTransform", "00", "01")
        vars2regress_in_dirname <- ifelse(is.na(vars2regress), "no_vars_to_regress", paste(gsub("[.]", "", vars2regress), collapse="-"))
        if (feat_param > 10) {
            # feat_param is a number of features
            integration_features_nb <- feat_param
            normalization_out_dir <- sprintf("%s/%s-%s/%s/featnb%d", out_dir, method_dir_index, method, vars2regress_in_dirname, feat_param)
            normalization_out_name <- sprintf("%s_featnb%d", method, feat_param)
        } else {
            # feat_param is a residual variance threshold
            ## get minimum number of variable features
            min_variable_features <- length(seurat_objects_norm[[sections[1]]]@assays$SCT@var.features)
            for (sample in sections[2:length(sections)]) {
                if (length(seurat_objects_norm[[sample]]@assays$SCT@var.features) < min_variable_features) {
                    min_variable_features <- length(seurat_objects_norm[[sample]]@assays$SCT@var.features)
                }
            }
            integration_features_nb <- min_variable_features
            normalization_out_dir <- sprintf("%s/%s-%s/%s/resvar%s", out_dir, method_dir_index, method, vars2regress_in_dirname, sub("[.]", "_", feat_param))
            normalization_out_name <- sprintf("%s_resvar%s", method, sub("[.]", "_", feat_param))
        }
        integration_normalization_method <- default_assay <- "SCT"
        selection_method <- "sct"
        feature_stat_slot <- "counts"
    } else {
        if (method == "LogNormalize") {
            selection_method <- "vst"
            # feat_param is a number of features
            seurat_objects_norm <- lapply(seurat_objects, function(x, nb_features=feat_param, selec=selection_method) {
                x <- NormalizeData(x)
                x <- FindVariableFeatures(x, selection.method=selec, nfeatures=nb_features)
            })
            integration_features_nb <- feat_param
            integration_normalization_method <- method
            default_assay <- "RNA"
            feature_stat_slot <- "data"
            normalization_out_dir <- sprintf("%s/10-%s/featnb%d", out_dir, method, feat_param)
            normalization_out_name <- sprintf("%s_featnb%d", method, feat_param)
        }
    }
    
    # add feature statistics to Seurat object after normalization
    seurat_objects_norm <- lapply(seurat_objects_norm, addFeatureStats2Seurat, assay2use=default_assay, slot2use=feature_stat_slot)
    
    # QC plots after normalization
    if (! dir.exists(normalization_out_dir)) {
        dir.create(normalization_out_dir, recursive=TRUE, mode="0775")
    }
    ## feature stat plots
    norm_column_prefix <- "norm"
    seurat_objects_norm <- lapply(seurat_objects_norm, addPercentageFeatureSet2Seurat, mito_genes=mito_genes, hb_genes=hb_genes, assay2use=default_assay, col_prefix=norm_column_prefix)
    library(RColorBrewer)
    lapply(seurat_objects_norm, seurat_all_feature_stat_plots, assay2use=default_assay, col_prefix=norm_column_prefix, plots_per_page=1, plot_title="Normalized counts", out_dir=normalization_out_dir, out_name=sprintf("%s_normalized_QC_stats", normalization_out_name)) 
    detach("package:RColorBrewer")
    ## mitochondrial and hemoglobin gene spatial plots after filerting
    ### hemoglobin genes
    df_list <- lapply(seurat_objects_norm, qc_spatial_mean_exp_plots, gene_vec=hb_genes, mean_order=TRUE, plots_per_page=NA, out_dir=normalization_out_dir, out_name=sprintf("%s_hemoglobin_genes", normalization_out_name))
    ### mitochondrial genes
    df_list <- lapply(seurat_objects_norm, qc_spatial_mean_exp_plots, gene_vec=mito_genes, mean_order=TRUE, plots_per_page=NA, out_dir=normalization_out_dir, out_name=sprintf("%s_mitochondrial_genes", normalization_out_name))
    
    return(list(seurat_list=seurat_objects_norm, features_nb=integration_features_nb, norm_method=integration_normalization_method, assay=default_assay, selec=selection_method, out_dir=normalization_out_dir, out_name=normalization_out_name))
}

seurat_merge_normalization <- function(seurat_merge, method, vars2regress, feat_param, only_var, mito_genes, hb_genes, plots_per_page, out_dir) {
    if (method == "SCTransform" | method == "SCTransform_v2") {
        vst_flavor <- switch(as.integer(method == "SCTransform_v2")+1, NULL, "v2")
        vars_to_regress <- switch(as.integer(is.na(vars2regress))+1, vars2regress, NULL)
        variable_features_n <- switch(as.integer(feat_param > 10)+1, NULL, feat_param) # if feat_param > 10, then it is a feature number, otherwise it is residual variance threshold
        variable_features_rv_th <- switch(as.integer(feat_param > 10)+1, feat_param, NULL) # if feat_param > 10, then it is a feature number, otherwise it is residual variance threshold
        print(sprintf("VST flavor: %s, variables to regress: %s, variable features n: %s, variable features rv th: %s", ifelse(is.null(vst_flavor), "NULL", vst_flavor), ifelse(is.null(vars_to_regress), "NULL", vars_to_regress), ifelse(is.null(variable_features_n), "NULL", variable_features_n), ifelse(is.null(variable_features_rv_th), "NULL", variable_features_rv_th)))
        seurat_merge_norm <- SCTransform(seurat_merge, vst.flavor=vst_flavor, vars.to.regress=vars_to_regress, variable.features.n=variable_features_n, variable.features.rv.th=variable_features_rv_th, return.only.var.genes=only_var, verbose=TRUE)
        
        method_dir_index <- ifelse(method == "SCTransform", "00", "01")
        vars2regress_in_dirname <- ifelse(is.na(vars2regress), "no_vars_to_regress", paste(gsub("[.]", "", vars2regress), collapse="-"))
        if (feat_param > 10) {
            # feat_param is a number of features
            normalization_out_dir <- sprintf("%s/%s-%s/%s/featnb%d", out_dir, method_dir_index, method, vars2regress_in_dirname, feat_param)
            normalization_out_name <- sprintf("%s_featnb%d", method, feat_param)
            integration_features_nb <- feat_param
        } else {
            # feat_param is a residual variance threshold
            normalization_out_dir <- sprintf("%s/%s-%s/%s/resvar%s", out_dir, method_dir_index, method, vars2regress_in_dirname, sub("[.]", "_", feat_param))
            normalization_out_name <- sprintf("%s_resvar%s", method, sub("[.]", "_", feat_param))
            ### get minimum number of variable features
            #### TODO
            #integration_features_nb <- min_variable_features
        }
        integration_normalization_method <- default_assay <- "SCT"
        selection_method <- "sct"
        feature_stat_slot <- "counts"
    } else {
        if (method == "LogNormalize") {
            selection_method <- "vst"
            ### feat_param is a number of features
            print(sprintf("features param: %d", feat_param))
            normalization_out_dir <- sprintf("%s/10-%s/featnb%d", out_dir, method, feat_param)
            if (! dir.exists(normalization_out_dir)) {
                dir.create(normalization_out_dir, recursive=TRUE, mode="0775")
            }
            normalization_out_name <- sprintf("%s_featnb%d", method, feat_param)

            seurat_merge_norm <- NormalizeData(seurat_merge)
            seurat_merge_norm <- FindVariableFeatures(seurat_merge_norm, selection.method=selection_method, nfeatures=feat_param)
            
            integration_features_nb <- feat_param
            integration_normalization_method <- method
            default_assay <- "RNA"
            feature_stat_slot <- "data"
        }
    }
    
    # add feature statistics to Seurat object after SCTransform
    seurat_merge_norm <- addFeatureStats2Seurat(seurat_merge_norm, default_assay, feature_stat_slot)
    
    # QC plots after normalization
    if (! dir.exists(normalization_out_dir)) {
        dir.create(normalization_out_dir, recursive=TRUE, mode="0775")
    }
    ## feature stat plots
    norm_column_prefix <- "norm"
    seurat_merge_norm <- addPercentageFeatureSet2Seurat(seurat_merge_norm, mito_genes, hb_genes, default_assay, norm_column_prefix)
    library(RColorBrewer)
    seurat_all_feature_stat_plots(seurat_merge_norm, default_assay, norm_column_prefix, plots_per_page, "Normalized counts", normalization_out_dir, "normalized_QC_stats") 
    detach("package:RColorBrewer")
    ## mitochondrial and hemoglobin gene spatial plots after filerting
    ### hemoglobin genes
    df_list <- qc_spatial_mean_exp_plots(seurat_merge_norm, hb_genes, TRUE, plots_per_page, normalization_out_dir, "normalized_hemoglobin_genes")
    ### mitochondrial genes
    df_list <- qc_spatial_mean_exp_plots(seurat_merge_norm, mito_genes, TRUE, plots_per_page, normalization_out_dir, "normalized_mitochondrial_genes")
    
    return(list(seurat_merge=seurat_merge_norm, features_nb=integration_features_nb, norm_method=integration_normalization_method, assay=default_assay, selec=selection_method, out_dir=normalization_out_dir, out_name=normalization_out_name))
}

hvg_plots <- function(seurat_obj, assay2use, selec, nb_hvgs, plots_per_page, out_dir, out_name) {
    sample <- seurat_obj@project.name
    
    # identify the 20 most highly variable genes
    top20 <- head(VariableFeatures(seurat_obj), 20)
    ## plot variable features with and without labels
    pdf(sprintf("%s/%s_%s_residual_mean_variance.pdf", out_dir, sample, out_name))
    plot1 <- VariableFeaturePlot(seurat_obj, selection.method=selec, assay=assay2use)
    plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
    print(plot2)
    dev.off()
    
    # get feature residual variances
    if (assay2use == "SCT") {
        feature_attributes_df <- seurat_obj@assays[[assay2use]]@SCTModel.list$model1@feature.attributes
        ## order according to decreasing residual variance
        feature_attributes_df <- feature_attributes_df[order(feature_attributes_df$residual_variance, decreasing=TRUE),]
    } else {
        feature_attributes_df <- seurat_obj@assays[[assay2use]]@meta.features
        ## order according to decreasing residual variance
        feature_attributes_df <- feature_attributes_df[order(feature_attributes_df$vst.variance.standardized, decreasing=TRUE),]
    }
    write.csv(feature_attributes_df, file=sprintf("%s/%s_%s_residual_mean_variance.csv", out_dir, sample, out_name), quote=FALSE, row.names=TRUE)
    
    # get HVGs
    hvgs <- rownames(feature_attributes_df)[1:nb_hvgs]
    
    # HVG spatial plots
    ## HVGs to plot: top 20, middle gene, last 20
    hvgs2plot <- c(hvgs[1:20], hvgs[floor(nb_hvgs/2)], hvgs[(nb_hvgs-20+1):nb_hvgs])
    exp_df <- qc_spatial_mean_exp_plots(seurat_obj, hvgs2plot, FALSE, plots_per_page, out_dir, out_name)
    
    return(exp_df)
}

#get_HVGs_merged_Seurat_SCT <- function(seurat_obj_merge, n_hvgs) {
#    hvgs <- c()
#    for (one_model in names(seurat_obj_merge@assays$SCT@SCTModel.list)) {
#        seurat_obj_merge_one_model_feature_attributes <- seurat_obj_merge@assays$SCT@SCTModel.list[[one_model]]@feature.attributes
#        hvgs <- c(hvgs, rownames(seurat_obj_merge_one_model_feature_attributes[order(seurat_obj_merge_one_model_feature_attributes$residual_variance, decreasing=TRUE),])[1:n_hvgs])
#    }
#    hvgs <- unique(hvgs)
#    # some genes are not present in all SCTModel: only keep those present in all models
#    for (one_model in names(seurat_obj_merge@assays$SCT@SCTModel.list)) {
#        hvgs <- hvgs[hvgs %in% rownames(seurat_obj_merge@assays$SCT@SCTModel.list[[one_model]]@feature.attributes)]
#    }
#    return(hvgs)
#}

#get_HVGs_merged_Seurat_LogNormalize <- function(seurat_obj_merge, seurat_obj_norm_list, n_hvgs) {
#    hvgs <- c()
#    for (one_sample in names(seurat_obj_merge@images)) {
#        hvgs <- c(hvgs, VariableFeatures(seurat_obj_norm_list[[one_sample]]))
#    }
#    hvgs <- unique(hvgs)
#    #### some genes are not expressed in all samples: only keep those expressed in all samples
#    for (one_sample in names(seurat_obj_merge@images)) {
#        hvgs <- hvgs[hvgs %in% rownames(seurat_obj_norm_list[[one_sample]]@assays$RNA@counts)]
#    }
#    return(hvgs)
#}

merge_stats_plots_after_norm <- function(seurat_obj_list, project_name, norm_method, feat_nb, mito_genes, hb_genes, plots_per_page, out_dir, out_name) {
    samples <- names(seurat_obj_list)
    nb_samples <- length(samples)
    seurat_object_merge <- merge(seurat_obj_list[[1]], y=unlist(seurat_obj_list)[2:nb_samples], add.cell.ids=samples)
    # set project name
    seurat_object_merge@project.name <- project_name
    levels(seurat_object_merge@assays$SCT) <- samples
    # orig.ident, lame, zone, condition and time are not factors in the merged object meta data
    seurat_object_merge@meta.data$orig.ident <- as.factor(seurat_object_merge@meta.data$orig.ident)
    seurat_object_merge@meta.data$lame <- as.factor(seurat_object_merge@meta.data$lame)
    seurat_object_merge@meta.data$zone <- as.factor(seurat_object_merge@meta.data$zone)
    seurat_object_merge@meta.data$condition <- as.factor(seurat_object_merge@meta.data$condition)
    seurat_object_merge@meta.data$time <- as.factor(seurat_object_merge@meta.data$time)
    # set variable features: after merging datasets, variable features are not set in the merged Seurat object, which makes RunPCA() fails: Error in PrepDR(object = object, features = features, verbose = verbose): Variable features haven't been set. Run FindVariableFeatures() or provide a vector of feature names.
    ## use SelectIntegrationFeatures() function
    seurat_object_merge_hvgs <- SelectIntegrationFeatures(seurat_obj_list, nfeatures=feat_nb)
    VariableFeatures(seurat_object_merge) <- seurat_object_merge_hvgs
    
    # HVGs after merge
    seurat_object_merge_hvgs <- VariableFeatures(seurat_object_merge)
    for (sample in names(seurat_object_merge@assays$SCT@SCTModel.list)) {
        sample_hvgs_after_merge_df <- seurat_object_merge@assays$SCT@SCTModel.list[[sample]]@feature.attributes[seurat_object_merge_hvgs,]
        write.csv(sample_hvgs_after_merge_df, file=sprintf("%s/%s_%s_HVGs_residual_mean_variance_after_merge.csv", out_dir, sample, out_name), quote=FALSE, row.names=TRUE)
    }
    
    # QC plots after normalization
    ## feature stat plots
    default_assay <- DefaultAssay(seurat_object_merge)
    ## add feature statistics to merged Seurat object
    feature_stat_slot <- ifelse(norm_method == "LogNormalize", "data", "counts")
    seurat_object_merge <- addFeatureStats2Seurat(seurat_object_merge, default_assay, feature_stat_slot)
    norm_column_prefix <- "norm"
    seurat_object_merge <- addPercentageFeatureSet2Seurat(seurat_object_merge, mito_genes, hb_genes, default_assay, norm_column_prefix)
    library(RColorBrewer)
    seurat_all_feature_stat_plots(seurat_object_merge, default_assay, norm_column_prefix, plots_per_page, "Normalized counts", out_dir, "QC_stats_normalization") 
    detach("package:RColorBrewer")
    ## mitochondrial and hemoglobin gene spatial plots after filerting
    ### hemoglobin genes
    df_list <- qc_spatial_mean_exp_plots(seurat_object_merge, hb_genes, TRUE, plots_per_page, out_dir, "normalized_hemoglobin_genes")
    ### mitochondrial genes
    df_list <- qc_spatial_mean_exp_plots(seurat_object_merge, mito_genes, TRUE, plots_per_page, out_dir, "normalized_mitochondrial_genes")
    
    return(seurat_object_merge)
}

set_SVGs_variable_features <- function(seurat_obj, assay2use) {
    seurat_obj_svgs <- FindSpatiallyVariableFeatures(seurat_obj, assay=assay2use, selection.method="markvariogram")
    # considering that seurat_obj_svgs is a character vector storing SVG names
    VariableFeatures(seurat_obj) <- seurat_obj_svgs
    return(seurat_obj)
}

#######################
# Dimension reduction #
#######################

dimension_reduction_outputs <- function(seurat_obj_reduc, dimred, out_dir, out_name) {
    # cell embeddings, feature loadings and dimension standard deviation
    write.csv(seurat_obj_reduc[[dimred]]@cell.embeddings, file=sprintf("%s/%s_cell_embeddings.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
    write.csv(seurat_obj_reduc[[dimred]]@feature.loadings, file=sprintf("%s/%s_feature_loadings.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
    pc_stdev_out_name <- sprintf("%s_stdev", out_name)
    df2ggplot <- data.frame(dimension=colnames(seurat_obj_reduc[[dimred]]@cell.embeddings), stdev=seurat_obj_reduc[[dimred]]@stdev)
    df2ggplot$dimension <- factor(df2ggplot$dimension, levels=df2ggplot$dimension)
    write.csv(df2ggplot, file=sprintf("%s/%s.csv", out_dir, pc_stdev_out_name), quote=FALSE, row.names=FALSE)
    pdf(sprintf("%s/%s.pdf", out_dir, pc_stdev_out_name))
    p <- ggplot(df2ggplot, aes(x=dimension, y=stdev)) +
        geom_bar(stat="identity", fill="steelblue1", color="steelblue") +
        labs(title=sprintf("Standard deviation of each dimension"), x="Dimension", y="Standard deviation percentage") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=60, hjust=1)) +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
    dev.off()
}

##############
# Clustering #
##############

compute_silhouette <- function(seurat_obj, reduc, dim_nb, metadata_colname) {
    categories <- as.factor(seurat_obj@meta.data[[metadata_colname]])
    if (length(levels(categories)) > 1) {
        distance_matrix <- dist(Embeddings(seurat_obj, reduction=reduc)[, 1:dim_nb])
        silhouette_obj <- silhouette(as.numeric(categories), dist=distance_matrix)
        seurat_obj@meta.data[[sprintf("silhouette_%s", metadata_colname)]] <- silhouette_obj[,3]
    } else {
        # only a single level: can not compute silhouette
        seurat_obj@meta.data[[sprintf("silhouette_%s", metadata_colname)]] <- rep(NA, length(dim(seurat_obj@meta.data)[1]))
    }
    return(seurat_obj)
}

compute_kBET <- function(embeddings, batch, pca, dims) {
    sample_size <- length(batch)
    df2return <- data.frame()
    for (one_pct in c(1, seq(5, 25, 5))) {
        k0 <- floor(sample_size * one_pct / 100)
        print(sprintf("pct sample size: %d, k0: %d", one_pct, k0))
        kBET_values <- kBET(embeddings, batch, k0=k0, do.pca=pca, dim.pca=dims, heuristic=FALSE, verbose=FALSE, plot=FALSE)
        if (length(kBET_values) > 1) {
            df2return <- rbind(df2return, data.frame(pct_sample_size=rep(one_pct, length(kBET_values$stats$kBET.observed)), kBET.rejection=kBET_values$stats$kBET.observed))
        }
    }
    df2return$kBET.acceptance <- 1 - df2return$kBET.rejection
    df2return$pct_sample_size <- as.factor(df2return$pct_sample_size)
    return(df2return)
}

kBET_plots <- function(kBET_df, out_dir, out_name) {
    # plot acceptance rate per sample size percentage
    write.csv(kBET_df, file=sprintf("%s/%s_kBET.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
    pdf(sprintf("%s/%s.pdf", out_dir, out_name))
    p <- ggplot(kBET_df, aes(x=pct_sample_size, y=kBET.acceptance, fill=pct_sample_size)) +
        geom_boxplot() +
        labs(title="kBET analysis", x="% sample size", y="Acceptance rate") +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
    
    # mean and median accpetance rate per sample size percentage
    mean_kBET_df <- aggregate(kBET_df$kBET.acceptance, by=list(kBET_df$pct_sample_size), FUN=mean)
    colnames(mean_kBET_df) <- c("pct_sample_size", "mean_kBET.acceptance")
    p <- ggplot(mean_kBET_df, aes(x=pct_sample_size, y=mean_kBET.acceptance, group=1)) +
        geom_point() +
        geom_line() +
        labs(title="kBET analysis", x="% sample size", y="Mean acceptance rate") +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
    median_kBET_df <- aggregate(kBET_df$kBET.acceptance, by=list(kBET_df$pct_sample_size), FUN=median)
    colnames(median_kBET_df) <- c("pct_sample_size", "median_kBET.acceptance")
    p <- ggplot(median_kBET_df, aes(x=pct_sample_size, y=median_kBET.acceptance, group=1)) +
        geom_point() +
        geom_line() +
        labs(title="kBET analysis", x="% sample size", y="Median acceptance rate") +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
    dev.off()
    
    df2return <- merge(mean_kBET_df, median_kBET_df, by="pct_sample_size")
    return(df2return)
}

clustering <- function(seurat_obj_int, rez, reduc, dim_nb) {
    # compute clustering
    seurat_obj_int <- FindClusters(seurat_obj_int, resolution=rez, verbose=FALSE)
    seurat_obj_int <- RunUMAP(seurat_obj_int, reduction=reduc, dims=1:dim_nb)
    # comupute silhouette scores
    seurat_obj_int <- compute_silhouette(seurat_obj_int, reduc, dim_nb, "seurat_clusters")
    
    # kBET
    ## batch
    #kBET_values <- kBET(Embeddings(seurat_obj_int$pca)[, 1:dim_nb], seurat_obj_int@meta.data$orig.ident, do.pca=FALSE, verbose=TRUE)
    
    #kBET_values <- kBET(as.matrix(t(seurat_obj_int@assays$SCT@counts)), seurat_obj_int@meta.data$orig.ident, do.pca=TRUE, dim.pca=dim_nb, verbose=TRUE)
    
    ### lame
    #### 50 PCs
    #kBET_values_lame <- kBET(Embeddings(seurat_objects_merge$pca), seurat_objects_merge@meta.data$lame)
    #### 20 PCs
    #kBET_values_inetgrated_lame_20PCs <- kBET(Embeddings(seurat_obj_int$pca)[,1:20], seurat_obj_int@meta.data$lame)

    return(seurat_obj_int)
}

clustering_plots <- function(seurat_obj_int, metadata_cat_vec, load_image, plots_per_page, plot_title, out_dir, out_name) {
    # plot clusters onto UMAP or onto the tissue section
    pdf(sprintf("%s/%s.pdf", out_dir, out_name))
    dimplot_list <- list()
    for (category_colname in c("seurat_clusters", metadata_cat_vec)) {
        p <- DimPlot(seurat_obj_int, reduction="umap", group.by=category_colname, label=TRUE)
        if (length(levels(seurat_obj_int@meta.data[[category_colname]])) > 1) {
            p <- p + ggtitle(sprintf("Mean silhouette for %s: %f", category_colname, mean(seurat_obj_int@meta.data[[sprintf("silhouette_%s", category_colname)]]))) + theme(plot.title=element_text(size=12))
        }
        dimplot_list[[category_colname]] <- p
    }
    for (one_metadata_cat in metadata_cat_vec) {
        print(dimplot_list[["seurat_clusters"]] / dimplot_list[[one_metadata_cat]] + plot_annotation(title=sprintf("%s", plot_title), theme=theme(plot.title=element_text(size=16)))) 
    }
    
    # silhouette: https://romanhaa.github.io/projects/scrnaseq_workflow/#silhouette-plot
    if (length(levels(seurat_obj_int@meta.data$seurat_clusters)) > 1) {
        ## silhouette only more than 1 cluster
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
    }
    
    if (load_image) {
        nb_sections <- length(unique(seurat_obj_int@meta.data$orig.ident))
        if (nb_sections > 4) {
            # one page per timepoint
            timepoints <- levels(seurat_obj_int@meta.data$time)
            for (one_timepoint in timepoints) {
                one_timepoint_metadata_df <- seurat_obj_int@meta.data[which(seurat_obj_int@meta.data$time==one_timepoint),]
                one_timepoint_metadata_df <- droplevels(one_timepoint_metadata_df)
                sections <- levels(one_timepoint_metadata_df$orig.ident)
                print(SpatialDimPlot(seurat_obj_int, images=sections, ncol=2, label=TRUE, label.size=2) + plot_annotation(title=sprintf("%s\ntimepoint: %s", plot_title, one_timepoint), theme=theme(plot.title=element_text(size=18))))
            }
        } else {
            for (i in 1:(nb_sections/plots_per_page)) {
                print(SpatialDimPlot(seurat_obj_int, images=names(seurat_obj_int@images)[((i-1)*plots_per_page+1):(i*plots_per_page)], ncol=2, label=TRUE, label.size=2) + plot_annotation(title=sprintf("%s", plot_title), theme=theme(plot.title=element_text(size=18))))
            }
        }
        for (i in 1:nb_sections) {
            print(SpatialDimPlot(seurat_obj_int, images=names(seurat_obj_int@images)[i], label=TRUE, label.size=3) + plot_annotation(title=sprintf("%s: %s", plot_title, names(seurat_obj_int@images)[i]), theme=theme(plot.title=element_text(size=16))))
            print(SpatialDimPlot(seurat_obj_int, images=names(seurat_obj_int@images)[i], label=FALSE, label.size=3, alpha=0) + plot_annotation(title=sprintf("%s: %s", plot_title, names(seurat_obj_int@images)[i]), theme=theme(plot.title=element_text(size=16))))
            print(SpatialDimPlot(seurat_obj_int, images=names(seurat_obj_int@images)[i], label=TRUE, label.size=3, alpha=0.1) + plot_annotation(title=sprintf("%s: %s", plot_title, names(seurat_obj_int@images)[i]), theme=theme(plot.title=element_text(size=16))))
            print(SpatialDimPlot(seurat_obj_int, images=names(seurat_obj_int@images)[i], label=TRUE, label.size=3, alpha=0.2) + plot_annotation(title=sprintf("%s: %s", plot_title, names(seurat_obj_int@images)[i]), theme=theme(plot.title=element_text(size=16))))
        }
    }
    
    # number of spots per cluster
    spot_count_df <- aggregate(seurat_obj_int@meta.data$seurat_clusters, by=list(seurat_obj_int@meta.data$orig.ident, seurat_obj_int@meta.data$seurat_clusters), FUN=length)
    colnames(spot_count_df) <- c("sample", "cluster", "count")
    spot_count_pct_df <- data.frame()
    for (one_sample in levels(spot_count_df$sample)) {
        total_count <- sum(spot_count_df[which(spot_count_df$sample == one_sample), "count"])
        spot_count_pct_df <- rbind(spot_count_pct_df, cbind(spot_count_df[which(spot_count_df$sample == one_sample),], pct=spot_count_df[which(spot_count_df$sample == one_sample), "count"] / total_count * 100))
    }
    write.csv(spot_count_pct_df, file=sprintf("%s/%s_spots.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
    
    p <- ggplot(spot_count_pct_df, aes(x=cluster, y=count, fill=sample)) +
        geom_bar(stat="identity", position="dodge") +
        labs(title="Number of spots per cluster", x="Cluster", y="Number of spots", fill="Sample") +
        theme_bw() +
        theme(legend.position="bottom") +
        theme(panel.border=element_rect(color="grey50"))
    print(p)

    p <- ggplot(spot_count_pct_df, aes(x=cluster, y=pct, fill=sample)) +
        geom_bar(stat="identity", position="dodge") +
        labs(title="Percentage of spots per cluster", x="Cluster", y="Percentage of spots", fill="Sample") +
        theme_bw() +
        theme(legend.position="bottom") +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
    
    p <- ggplot(spot_count_pct_df, aes(x=sample, y=pct, fill=cluster)) +
        geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
        labs(title="Percentage of spots per cluster", x="Cluster", y="Percentage of spots", fill="Sample") +
        theme_bw() +
        theme(legend.position="bottom") +
        theme(axis.text.x=element_text(angle=30, hjust=1)) +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
    
    p <- ggplot(spot_count_pct_df, aes(x=sample, y=count, fill=cluster)) +
        geom_bar(stat="identity", position=position_fill(reverse=TRUE)) +
        labs(title="Percentage of spots per cluster", x="Cluster", y="Percentage of spots", fill="Sample") +
        theme_bw() +
        theme(legend.position="bottom") +
        theme(axis.text.x=element_text(angle=30, hjust=1)) +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
    
    dev.off()
}


clustering_kBET <- function(seurat_obj, out_dir, out_name) {
    
    out_name <- sprintf("")
    cluster_kBET_df2ggplot <- data.frame()
    for (one_cluster in levels(seurat_obj@meta.data$seurat_clusters)) {
        print(sprintf("cluster: %s", one_cluster))
        barcodes2use <- rownames(seurat_obj@meta.data[which(seurat_obj@meta.data$seurat_clusters == one_cluster),])
        one_cluster_embeddings <- Embeddings(seurat_obj@reductions[[reduction2use]])[barcodes2use,]
        one_cluster_batch <- seurat_obj@meta.data[barcodes2use,"orig.ident"]
        kBET_df <- compute_kBET(one_cluster_embeddings, one_cluster_batch, kBET_pca_param, kBET_dims_param)
        cluster_kBET_df2ggplot <- rbind(cluster_kBET_df2ggplot, cbind(cluster=rep(one_cluster, dim(kBET_df)[1]), kBET_df))
    }
    write.csv(cluster_kBET_df2ggplot, file=sprintf("%s/%s_kBET_cluster.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
    
    pdf(sprintf("%s/%s_kBET_cluster.pdf", out_dir, out_name))
    cluster_kBET_df2ggplot$cluster <- as.factor(cluster_kBET_df2ggplot$cluster)
    for (one_pct in levels(cluster_kBET_df2ggplot$pct_sample_size)) {
        cluster_kBET_one_pct_df2ggplot <- cluster_kBET_df2ggplot[which(cluster_kBET_df2ggplot$pct_sample_size == one_pct),]
        mean_kBET <- mean(cluster_kBET_one_pct_df2ggplot$kBET.acceptance)
        median_kBET <- median(cluster_kBET_one_pct_df2ggplot$kBET.acceptance)
        p <- ggplot(cluster_kBET_one_pct_df2ggplot, aes(x=cluster, y=kBET.acceptance, fill=cluster)) +
            geom_boxplot() +
            geom_hline(yintercept = mean_kBET, color="red", linetype="dashed") +
            geom_hline(yintercept = median_kBET, color="blue", linetype="dashed") +
            labs(title=sprintf("kBET analysis per cluster: sample size=%s%%", one_pct), x="Cluster", y="Acceptance rate") +
            theme_bw() +
            theme(panel.border=element_rect(color="grey50"))
        print(p)
    }
    for (one_cluster in levels(cluster_kBET_df2ggplot$cluster)) {
        cluster_kBET_one_cluster_df2ggplot <- cluster_kBET_df2ggplot[which(cluster_kBET_df2ggplot$cluster == one_cluster),]
        mean_kBET <- mean(cluster_kBET_one_cluster_df2ggplot$kBET.acceptance)
        median_kBET <- median(cluster_kBET_one_cluster_df2ggplot$kBET.acceptance)
        p <- ggplot(cluster_kBET_one_cluster_df2ggplot, aes(x=pct_sample_size, y=kBET.acceptance, fill=pct_sample_size)) +
            geom_boxplot() +
            geom_hline(yintercept = mean_kBET, color="red", linetype="dashed") +
            geom_hline(yintercept = median_kBET, color="blue", linetype="dashed") +
            labs(title=sprintf("kBET analysis per cluster: cluster %s", one_cluster), x="% sample size", y="Acceptance rate") +
            theme_bw() +
            theme(panel.border=element_rect(color="grey50"))
        print(p)
    }

    mean_cluster_kBET_df2ggplot <- aggregate(cluster_kBET_df2ggplot$kBET.acceptance, by=list(cluster_kBET_df2ggplot$pct_sample_size, cluster_kBET_df2ggplot$cluster), FUN=mean)
    colnames(mean_cluster_kBET_df2ggplot) <- c("pct_sample_size", "cluster", "mean_kBET.acceptance")
    p <- ggplot(mean_cluster_kBET_df2ggplot, aes(x=pct_sample_size, y=mean_kBET.acceptance, group=cluster, col=cluster)) +
        geom_point() +
        geom_line() +
        labs(title="kBET analysis per cluster", x="Sample size", y="Mean acceptance rate", col="Cluster") +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
    median_cluster_kBET_df2ggplot <- aggregate(cluster_kBET_df2ggplot$kBET.acceptance, by=list(cluster_kBET_df2ggplot$pct_sample_size, cluster_kBET_df2ggplot$cluster), FUN=median)
    colnames(median_cluster_kBET_df2ggplot) <- c("pct_sample_size", "cluster", "median_kBET.acceptance")
    p <- ggplot(median_cluster_kBET_df2ggplot, aes(x=pct_sample_size, y=median_kBET.acceptance, group=cluster, col=cluster)) +
        geom_point() +
        geom_line() +
        labs(title="kBET analysis per cluster", x="Sample size", y="Median acceptance rate", col="Cluster") +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
    dev.off()
    
    df2return <- merge(mean_cluster_kBET_df2ggplot, median_cluster_kBET_df2ggplot, by=c("pct_sample_size", "cluster"))
    return(df2return)
}


###########
# Markers #
###########

identify_cluster_markers <- function(seurat_obj_int, assay2use, condition_col, condition_subset, cluster_col, scale_split, scale_regress, test, latent_vars, load_image, nb_plots, out_dir, out_name) {
    # subset Seurat object to samples of condition
    Idents(seurat_obj_int) <- seurat_obj_int@meta.data[[condition_col]]
    seurat_obj_int_subset <- subset(seurat_obj_int, idents=condition_subset)
    ## after subsetting, only keep image of the samples in the subset
    seurat_obj_int_subset@meta.data <- droplevels(seurat_obj_int_subset@meta.data)
    samples <- levels(seurat_obj_int_subset@meta.data$orig.ident)
    seurat_obj_int_subset@images <- seurat_obj_int_subset@images[samples]
    
    # set identity classes to clusters
    Idents(seurat_obj_int_subset) <- seurat_obj_int_subset@meta.data[[cluster_col]]
    latent_vars_param <- ifelse(length(latent_vars) == 1, switch(as.integer(is.na(latent_vars))+1, latent_vars, NULL), latent_vars)
    print("Identify cluster markers")
    if (DefaultAssay(seurat_obj_int_subset) == "integrated" & assay2use == "RNA") {
        # scale.data slot of RNA assay is empty, which makes DoHeatmap() function fails: No requested features found in the scale.data slot for the RNA assay.
        if (! is.na(scale_split)) {
            if (!is.na(scale_regress)) {
                seurat_obj_int_subset_to_heatmap <- ScaleData(seurat_obj_int_subset, split.by=scale_split, vars.to.regress=scale_regress)
            } else {
                seurat_obj_int_subset_to_heatmap <- ScaleData(seurat_obj_int_subset, split.by=scale_split)
            }
        } else {
            seurat_obj_int_subset_to_heatmap <- ScaleData(seurat_obj_int_subset)
        }
    } else {
        seurat_obj_int_subset_to_heatmap <- seurat_obj_int_subset
    }
    
    # set default assay: must be either 'RNA' or 'SCT'
    DefaultAssay(seurat_obj_int_subset) <- assay2use
    if (test == "negbinom") {
        mean_function <- function(x, pseudocount.use=1, base=2) {
            return(log(x = rowMeans(x = x) + pseudocount.use, base = base))
        }
        markers <- FindAllMarkers(seurat_obj_int_subset, logfc.threshold=0, test.use=test, latent.vars=latent_vars_param, mean.fxn=mean_function, only.pos=TRUE, recorrect_umi=FALSE)
    } else {
        markers <- FindAllMarkers(seurat_obj_int_subset, logfc.threshold=0, test.use=test, latent.vars=latent_vars_param, only.pos=TRUE, recorrect_umi=FALSE)
    }
    write.csv(markers, file=sprintf("%s/%s.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
    pdf(sprintf("%s/%s.pdf", out_dir, out_name))
    ## expression heatmap
    top10_markers <- markers %>%
        group_by(cluster) %>%
        top_n(n=10, wt=avg_log2FC)
    print(DoHeatmap(seurat_obj_int_subset_to_heatmap, features=top10_markers$gene) + NoLegend())
    dev.off()
    ## get top markers for each cluster
    if (load_image) {
        top_markers <- markers %>%
            group_by(cluster) %>%
            slice_max(n=nb_plots, order_by=avg_log2FC)
        top_markers_pct2 <- markers[which(markers$pct.2 < 0.3),] %>%
            group_by(cluster) %>%
            slice_max(n=nb_plots, order_by=avg_log2FC)
        nb_sections <- length(unique(seurat_obj_int_subset@meta.data$orig.ident))
        for (one_cluster in levels(markers$cluster)) {
            print(sprintf("cluster: %s", one_cluster))
            # overall top cluster markers
            one_cluster_markers <- top_markers[which(top_markers$cluster==one_cluster), "gene", drop=TRUE]
            pdf(sprintf("%s/%s_c%s_overall.pdf", out_dir, out_name, sub("[.]", "_", one_cluster)))
            for (i in 1:nb_plots) {
                violin_plot <- VlnPlot(object=seurat_obj_int_subset, features=one_cluster_markers[i], log=TRUE)
                #spatial_plot <- SpatialFeaturePlot(object=seurat_obj_int_subset, features=one_cluster_markers[i], images=names(seurat_obj_int_subset@images)[((j-1)*plots_per_page+1):(j*plots_per_page)], alpha=c(0.1,1), ncol=plots_per_page) + plot_annotation(title=sprintf("Cluster: %s", one_cluster), theme=theme(plot.title=element_text(size=16)))
                spatial_plot <- SpatialFeaturePlot(object=seurat_obj_int_subset, features=one_cluster_markers[i], alpha=c(0.1,1), ncol=4) + plot_annotation(title=sprintf("Cluster: %s", one_cluster), theme=theme(plot.title=element_text(size=16)))
                print(violin_plot / spatial_plot + plot_layout(heights=c(1,3)))
            }
            dev.off()
            # top cluster markers expressed in less than 30% of the spots of all other clusters
            one_cluster_markers <- top_markers_pct2[which(top_markers_pct2$cluster==one_cluster), "gene", drop=TRUE]
            pdf(sprintf("%s/%s_c%s_pct2_30.pdf", out_dir, out_name, sub("[.]", "_", one_cluster)))
            for (i in 1:nb_plots) {
                violin_plot <- VlnPlot(object=seurat_obj_int_subset, features=one_cluster_markers[i], log=TRUE)
                #spatial_plot <- SpatialFeaturePlot(object=seurat_obj_int_subset, features=one_cluster_markers[i], images=names(seurat_obj_int_subset@images)[((j-1)*plots_per_page+1):(j*plots_per_page)], alpha=c(0.1,1), ncol=plots_per_page) + plot_annotation(title=sprintf("Cluster: %s", one_cluster), theme=theme(plot.title=element_text(size=16)))
                spatial_plot <- SpatialFeaturePlot(object=seurat_obj_int_subset, features=one_cluster_markers[i], alpha=c(0.1,1), ncol=4) + plot_annotation(title=sprintf("Cluster: %s", one_cluster), theme=theme(plot.title=element_text(size=16)))
                print(violin_plot / spatial_plot + plot_layout(heights=c(1,3)))
            }
            dev.off()
        }
    }
}

GO_analysis_plots <- function(go_results, gsea, sim_method, sem_data, out_dir, out_name) {
    pdf(sprintf("%s/%s.pdf", out_dir, out_name), width=12)
    signif_go_terms <- dim(go_results@result[which(go_results@result$p.adjust < 0.05),])[1]
    # set max number of GO terms to plot to 200 if more than 200 enriched GO terms (to speed up calculations)
    maxterms2plot <- ifelse(signif_go_terms > 200, 200, signif_go_terms)
    # barplot / dotplot
    if (signif_go_terms > 30) {
        for (i in seq(1, maxterms2plot, 30)) {
            if (i + 30 < maxterms2plot) {
                goterms2plot <- go_results$Description[i:(i+30-1)]
                #print(dotplot(go_results, showCategory=go_results$Description[i:(i+30-1)]))
                #print(barplot(go_results, showCategory=go_results$Description[i:(i+30-1)], font.size=6))
            } else {
                goterms2plot <- go_results$Description[i:maxterms2plot]
                #print(dotplot(go_results, showCategory=go_results$Description[i:maxterms2plot]))
                #print(barplot(go_results, showCategory=go_results$Description[i:(i+30-1)], font.size=6))
            }
            if (gsea) {
                print(dotplot(go_results, showCategory=goterms2plot, font.size=6))
            } else {
                print(barplot(go_results, showCategory=goterms2plot, font.size=6))
            }
        }
    } else {
        if (gsea) {
            print(dotplot(go_results, showCategory=maxterms2plot))
        } else {
            print(barplot(go_results, showCategory=maxterms2plot))
        }
        
    }
    
    # compute GO term similarity, if needed
    if (dim(go_results@termsim)[1] == 0) {
        go_results_termsim <- pairwise_termsim(go_results, method=sim_method, semData=sem_data, showCategory=maxterms2plot)
    } else {
        go_results_termsim <- go_results
    }
    
    # tree plot: similarity between enriched GO terms
    if (signif_go_terms > 4) {
        p <- treeplot(go_results_termsim, offset.params=list(bar_tree=rel(5), tiplab=rel(3), extend=0.3, hexpand=0.1)) # TODO: understant and use label_format_tiplab parameter to make labels readable
        p <- p + labs(title=sprintf("30/%d most enriched GO terms, clustering method: kmeans", signif_go_terms))
        print(p)
        ## print maximum number of enriched GO terms
        for (nb_clusters in c(5, 10, 15)) {
            if (maxterms2plot > nb_clusters) {
                p <- treeplot(go_results_termsim, showCategory=maxterms2plot, cluster.params=list(n=nb_clusters), offset.params=list(bar_tree=rel(5), tiplab=rel(3), extend=0.3, hexpand=0.1))
                p <- p + labs(title=sprintf("%d/%d enriched GO terms, clustering method: kmeans", maxterms2plot, signif_go_terms))
                print(p)
            }
        }
        ## clustering method: average
        p <- treeplot(go_results_termsim, cluster.params=list(method="average"), offset.params=list(bar_tree=rel(5), tiplab=rel(3), extend=0.3, hexpand=0.1))
        p <- p + labs(title=sprintf("30/%d most enriched GO terms, clustering method: average", signif_go_terms))
        print(p)
        ## print maximum number of enriched GO terms
        for (nb_clusters in c(5, 10, 15)) {
            if (maxterms2plot > nb_clusters) {
                p <- treeplot(go_results_termsim, showCategory=maxterms2plot, cluster.params=list(method="average", n=nb_clusters), offset.params=list(bar_tree=rel(5), tiplab=rel(3), extend=0.3, hexpand=0.1))
                p <- p + labs(title=sprintf("%d/%d enriched GO terms, clustering method: average", maxterms2plot, signif_go_terms))
                print(p)
            }
        }
    }
    
    # enrichment map
    p <- emapplot(go_results_termsim, cex.params=list(category_label=0.5, category_node=0.5, line=0.5))
    print(p)
    ## with clusters
    if (maxterms2plot > 2 * 5) {
        p <- emapplot(go_results_termsim, cex.params=list(category_label=0.5, category_node=1, line=0.5), cluster.params=list(cluster=TRUE, legend=TRUE, n=5, label_words_n=3))
        p <- p + labs(title=sprintf("30/%d most enriched GO terms", signif_go_terms))
        print(p)
    }
    ## with clusters, print maximum number of enriched GO terms
    for (nb_clusters in c(5, 10, 15)) {
        if (maxterms2plot > 2 * nb_clusters) {
            p <- emapplot(go_results_termsim, showCategory=maxterms2plot, cex.params=list(category_label=0.5, category_node=1, line=0.5), cluster.params=list(cluster=TRUE, legend=TRUE, n=nb_clusters, label_words_n=3))
            p <- p + labs(title=sprintf("%d/%d most enriched GO terms", maxterms2plot, signif_go_terms))
            print(p)
        }
    }
    
    if (gsea) {
        # ridgeline plot for expression distribution    
        print(ridgeplot(go_results))
        # GSEA plots
        nb_go_terms_to_plot <- ifelse(signif_go_terms < 30, signif_go_terms, 30)
        for (i in 1:nb_go_terms_to_plot) {
            print(gseaplot(go_results, geneSetID=i, title=go_results$Description[i]))
        }
        nb_go_terms_to_plot <- ifelse(signif_go_terms < 10, signif_go_terms, 10)
        print(gseaplot2(go_results, geneSetID=1:nb_go_terms_to_plot))
    }
    dev.off()
}

GO_ORA_analysis <- function(de_entrez_df, universe_entrez_df, sim_method, sem_data, out_dir, out_name) {
    ego <- enrichGO(gene=de_entrez_df$ENTREZID[!is.na(de_entrez_df$ENTREZID)], universe=universe_entrez_df$ENTREZID[!is.na(universe_entrez_df$ENTREZID)], OrgDb="org.Rn.eg.db", ont="BP", pAdjustMethod="BH", readable=TRUE)
    write.table(ego, file=sprintf("%s/%s.tsv", out_dir, out_name), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
    ego_signif <- ego@result[which(ego@result$p.adjust < 0.05),]
    signif_go_terms <- dim(ego_signif)[1]
    if (signif_go_terms > 0) {
        ### visualization
        GO_analysis_plots(ego, FALSE, sim_method, sem_data, out_dir, out_name)
        
        ### simplify
        #### compute similarity for all significant GO terms
        ego_termsim <- pairwise_termsim(ego, method=sim_method, semData=sem_data, showCategory=signif_go_terms)
        write.csv(ego_termsim@termsim, file=sprintf("%s/%s_sim_%s.csv", out_dir, out_name, sim_method), quote=FALSE, row.names=TRUE)
        #### simplify
        ego_simp_list <- list()
        for (simp_cutoff in c(0.5, 0.6, 0.7, 0.8)) {
            ego_termsim_simplify <- simplify(ego_termsim, measure=sim_method, cutoff=simp_cutoff, semData=sem_data)
            simplify_out_name <- sprintf("%s_simp_%s_%s", out_name, sim_method, sub("[.]", "_", simp_cutoff))
            write.table(ego_termsim_simplify, file=sprintf("%s/%s.tsv", out_dir, simplify_out_name), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
            ego_termsim_simplify_signif <- ego_termsim_simplify@result[which(ego_termsim_simplify@result$p.adjust < 0.05),]
            signif_go_terms <- dim(ego_termsim_simplify_signif)[1]
            ### visualization
            if (signif_go_terms > 0) {
                GO_analysis_plots(ego_termsim_simplify, FALSE, NULL, NULL, out_dir, simplify_out_name)
            }
            ego_simp_list[[sprintf("%.1f", simp_cutoff)]] <- ego_termsim_simplify_signif
        }
    } else {
        ego_simp_list <- list()
    }
    
    return(list(whole=ego_signif, simplify=ego_simp_list))
}

GSEA_analysis <- function(de_df, entrez_df, pval_bool, sem_data, sim_method, out_dir, out_name) {
    de_entrez_df <- merge(entrez_df, de_df, by.x="SYMBOL", by.y="row.names")
    adj_pval_no_null <- de_entrez_df$p_val_adj
    adj_pval_no_null[adj_pval_no_null  == 0] <- min(adj_pval_no_null[adj_pval_no_null != 0])
    if (pval_bool) {
        de_geneList <- de_entrez_df$avg_log2FC * (-log(adj_pval_no_null, 10))
    } else {
        de_geneList <- de_entrez_df$avg_log2FC
    }
    de_entrez_df$GSEA_val <- de_geneList
    names(de_geneList) <- de_entrez_df$ENTREZID
    de_geneList <- sort(de_geneList, decreasing=TRUE)
    de_entrez_df2GSEA <- de_entrez_df[match(names(de_geneList), de_entrez_df$ENTREZID),]
    write.csv(de_entrez_df2GSEA, file=sprintf("%s/%s_input.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
    gsea <- gseGO(geneList=de_geneList, OrgDb="org.Rn.eg.db", ont="BP", pvalueCutoff=0.05, verbose=FALSE)
    write.table(gsea, file=sprintf("%s/%s_results.tsv", out_dir, out_name), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
    gsea_signif <- gsea@result[which(gsea@result$p.adjust < 0.05),]
    
    # split results depending on up or down regulation (NES sign)
    gseaResult_up_down_list <- list()
    ## up
    gsea_up <- gsea[(gsea$p.adjust < 0.05 & gsea$NES > 0), asis=TRUE]
    gseaResult_up_down_list$up <- gsea_up
    ## down
    gsea_down <- gsea[(gsea$p.adjust < 0.05 & gsea$NES < 0), asis=TRUE]
    gseaResult_up_down_list$down <- gsea_down
    
    gsea_up_down_list <- list()
    for (one_gsea_up_down in names(gseaResult_up_down_list)) {
        one_gsea_up_down_list <- list()
        up_down_out_name <- sprintf("%s_%s", out_name, one_gsea_up_down)
        gsea_up_down <- gseaResult_up_down_list[[one_gsea_up_down]]
        gsea_up_down_signif <- gsea_up_down@result[which(gsea_up_down@result$p.adjust < 0.05),]
        signif_go_terms <- dim(gsea_up_down_signif)[1]
        one_gsea_up_down_list[["whole"]] <- gsea_up_down_signif
        if (signif_go_terms > 0) {
            ### visualization
            #GSEA_plots(gsea, out_dir, out_name)
            GO_analysis_plots(gsea_up_down, TRUE, sim_method, sem_data, out_dir, up_down_out_name)

            ### simplify
            one_gsea_up_down_simp_list <- list()
            #### compute similarity for all significant GO terms
            gsea_up_down_termsim <- pairwise_termsim(gsea_up_down, method=sim_method, semData=sem_data, showCategory=signif_go_terms)
            write.csv(gsea_up_down_termsim@termsim, file=sprintf("%s/%s_sim_%s.csv", out_dir, up_down_out_name, sim_method), quote=FALSE, row.names=TRUE)
            #### simplify
            gsea_up_down_simp_list <- list()
            for (simp_cutoff in c(0.5, 0.6, 0.7, 0.8)) {
                gsea_up_down_termsim_simplify <- simplify(gsea_up_down_termsim, measure=sim_method, cutoff=simp_cutoff, semData=sem_data)
                simplify_out_name <- sprintf("%s_simp_%s_%s", up_down_out_name, sim_method, sub("[.]", "_", simp_cutoff))
                write.table(gsea_up_down_termsim_simplify, file=sprintf("%s/%s.tsv", out_dir, simplify_out_name), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
                gsea_up_down_termsim_simplify_signif <- gsea_up_down_termsim_simplify@result[which(gsea_up_down_termsim_simplify@result$p.adjust < 0.05),]
                signif_go_terms <- dim(gsea_up_down_termsim_simplify_signif)[1]
                ### visualization
                if (signif_go_terms > 0) {
                    GO_analysis_plots(gsea_up_down_termsim_simplify, TRUE, NULL, NULL, out_dir, simplify_out_name)
                }
                one_gsea_up_down_simp_list[[sprintf("%.1f", simp_cutoff)]] <- gsea_up_down_termsim_simplify_signif
            }
        } else {
            one_gsea_up_down_simp_list <- list()
        }
        one_gsea_up_down_list[["simplify"]] <- one_gsea_up_down_simp_list
        gsea_up_down_list[[one_gsea_up_down]] <- one_gsea_up_down_list
    }
    
    return(gsea_up_down_list)
}

GO_analysis <- function(de_df, fc_threshold, sim_method, sem_data, gsea_bool, out_dir, out_name) {
    # get DE genes and universe gene set
    up_results_df <- de_df[which(de_df$p_val_adj < 0.05 & de_df$avg_log2FC > fc_threshold),]
    write.csv(up_results_df, file=sprintf("%s/%s_up_DE_genes.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
    up_genes <- rownames(up_results_df)
    down_results_df <- de_df[which(de_df$p_val_adj < 0.05 & de_df$avg_log2FC < fc_threshold),]
    write.csv(down_results_df, file=sprintf("%s/%s_down_DE_genes.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
    down_genes <- rownames(down_results_df)
    universe_genes <- rownames(de_df)
    
    # convert RefSeq IDs into Entrez IDs
    up_genes_entrez_df <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Rn.eg.db", drop=FALSE)
    write.csv(up_genes_entrez_df, file=sprintf("%s/%s_up_DE_genes_Entrez.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
    down_genes_entrez_df <- bitr(down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Rn.eg.db", drop=FALSE)
    write.csv(down_genes_entrez_df, file=sprintf("%s/%s_down_DE_genes_Entrez.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
    universe_genes_entrez_df <- bitr(universe_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Rn.eg.db", drop=FALSE)
    write.csv(universe_genes_entrez_df, file=sprintf("%s/%s_universe_genes_Entrez.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
    
    # GO term enrichment analysis
    ## up DE genes
    go_out_name <- sprintf("%s_GO_ORA_up", out_name)
    if (dim(up_genes_entrez_df)[1] > 0) {
        up_genes_ego_list <- GO_ORA_analysis(up_genes_entrez_df, universe_genes_entrez_df, sim_method, sem_data, out_dir, go_out_name)
    } else {
        ### empty results
        up_genes_ego_list <- list(whole=data.frame(), simplify=list())
    }
    ## down DE genes
    go_out_name <- sprintf("%s_GO_ORA_down", out_name)
    if (dim(down_genes_entrez_df)[1] > 0) {
        down_genes_ego_list <- GO_ORA_analysis(down_genes_entrez_df, universe_genes_entrez_df, sim_method, sem_data, out_dir, go_out_name)
    } else {
        ### empty results
        down_genes_ego_list <- list(whole=data.frame(), simplify=list())
    }
    
    # GSEA
    if (gsea_bool) {
        ## all expressed genes
        ### do not use p-value to order genes
        gsea_out_name <- sprintf("%s_GSEA", out_name)
        gsea <- GSEA_analysis(de_df, universe_genes_entrez_df, FALSE, sem_data, sim_method, out_dir, gsea_out_name)
        ### use p-value to order genes
        if (dim(up_genes_entrez_df)[1] > 0 & dim(down_genes_entrez_df)[1] > 0) {
            gsea_out_name <- sprintf("%s_GSEA_pval", out_name)
            gsea_pval <- GSEA_analysis(de_df, universe_genes_entrez_df, TRUE, sem_data, sim_method, out_dir, gsea_out_name)
        } else {
            ### empty results
            gsea_pval <- list(up=list(whole=data.frame(), simplify=list()), down=list(whole=data.frame(), simplify=list()))
        }
    } else {
        # empty results
        gsea <- gsea_pval <- list(up=list(whole=data.frame(), simplify=list()), down=list(whole=data.frame(), simplify=list()))
    }
    
    return(list(go_ora_up=up_genes_ego_list[["whole"]], go_ora_up_simp=up_genes_ego_list[["simplify"]], go_ora_down=down_genes_ego_list[["whole"]], go_ora_down_simp=down_genes_ego_list[["simplify"]], gsea_up=gsea[["up"]][["whole"]], gsea_up_simp=gsea[["up"]][["simplify"]], gsea_down=gsea[["down"]][["whole"]], gsea_down_simp=gsea[["down"]][["simplify"]], gsea_pval_up=gsea_pval[["up"]][["whole"]], gsea_pval_up_simp=gsea_pval[["up"]][["simplify"]], gsea_pval_down=gsea_pval[["down"]][["whole"]], gsea_pval_down_simp=gsea_pval[["down"]][["simplify"]]))
}

go_term_count_time_cluster <- function(df, category_cols, plot_title, out_dir, out_name) {
    aggregate_list <- list()
    for (one_category_col in category_cols) {
        aggregate_list[[length(aggregate_list)+1]] <- df[[one_category_col]]
    }
    count_df <- aggregate(df$ID, by=aggregate_list, length)
    colnames(count_df) <- c(category_cols, "count")
    count_df$time <- factor(count_df$time, levels=sort(as.numeric(unique(count_df$time))))
    count_df$cluster <- as.factor(count_df$cluster)
    count_df$category <- as.factor(count_df$category)
    count_out_name <- sprintf("%s_counts", out_name)
    write.csv(count_df, file=sprintf("%s/%s.csv", out_dir, count_out_name), quote=FALSE, row.names=TRUE)
    
    pdf(sprintf("%s/%s.pdf", out_dir, count_out_name))
    p <- ggplot(count_df, aes(x=time, y=count, fill=cluster)) +
        geom_bar(stat="identity", position="dodge")
    if ("similarity_cutoff" %in% category_cols) {
        p <- p + facet_grid(category~similarity_cutoff)
    } else {
        p <- p + facet_grid(~category)
    }
    p <- p + labs(title=plot_title, x="Time", y="Count", fill="Cluster") +
        theme_bw() +
        theme(legend.position="bottom") +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
    dev.off()
}

volcano_plots <- function(de_df, fc_threshold, plot_title, out_dir, out_name) {
    if (fc_threshold == 0) {
        DE_cat <- ifelse(de_df$p_val_adj < 0.05, ifelse(de_df$avg_log2FC > 0, "up", "down"), "non-DE")
    } else {
        DE_cat <- ifelse(de_df$p_val_adj < 0.05, ifelse(de_df$avg_log2FC > 0.25, "up", ifelse(de_df$avg_log2FC < -0.25, "down", "non-DE")), "non-DE")
    }
    DE_counts <- table(DE_cat)
    n_up <- DE_counts["up"]
    n_down <- DE_counts["down"]
    DE_label <- c(rownames(de_df)[1:30], rep(NA, length(rownames(de_df))-30))
    DE_cols <- c("up"="Red", "down"="Blue", "non-DE"="Black")
    DE_alphas <- c("up"=1, "down"=1, "non-DE"=0.5)
    x_limit <- max(abs(min(de_df$avg_log2FC)), max(de_df$avg_log2FC))
    volcano_df2ggplot <- data.frame(fc=de_df$avg_log2FC, pval=-log(de_df$p_val_adj, 10), category=DE_cat, gene=rownames(de_df), de_label=DE_label)
    write.csv(volcano_df2ggplot, file=sprintf("%s/%s_volcano.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
    p <- ggplot(volcano_df2ggplot, aes(x=fc, y=pval, col=category, alpha=category, label=de_label)) +
        geom_point() +
        geom_text(nudge_y=1, check_overlap=TRUE) +
        geom_hline(yintercept=-log(0.05, 10), linetype='dashed') +
        xlim(-x_limit, x_limit) +
        scale_colour_manual(values=DE_cols) +
        scale_alpha_manual(values=DE_alphas) +
        labs(title=sprintf("%s\n%s up genes, %s down genes", plot_title, n_up, n_down), x="log2(fold change)", y="-log10(p-value)", col="Category") +
        guides(alpha="none") +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50"))
    if (fc_threshold == 0) {
        p <- p + geom_vline(xintercept=0, linetype='dashed')
    } else {
        p <- p + geom_vline(xintercept=c(-fc_threshold, fc_threshold), linetype='dotdash')
    }
    print(p)
}

identify_DE_markers <- function(seurat_obj, cluster_col, condition_col, test, latent_vars, recorrect, nb_plots, plots_per_page, out_dir, out_name) {
    # set identity classes to clusters and conditions
    cluster_condition_col <- sprintf("%s.%s", cluster_col, condition_col)
    Idents(seurat_obj) <- seurat_obj@meta.data[[cluster_condition_col]]
    # differentially expressed genes between SE and CTRL spots
    seurat_obj.DE.all_clusters <- data.frame()
    for (one_cluster in levels(seurat_obj@meta.data[[cluster_col]])) {
        print(sprintf("cluster: %s", one_cluster))
        SE_ident <- sprintf("%s_SE", one_cluster)
        CTRL_ident <- sprintf("%s_CTRL", one_cluster)
        ## analysis only for clusters present in both conditions
        if (SE_ident %in% levels(Idents(seurat_obj)) & CTRL_ident %in% levels(Idents(seurat_obj))) {
            ### analysis only clusters with at least 3 spots for each condition
            SE_spots <- dim(seurat_obj@meta.data[which(seurat_obj@meta.data[[cluster_condition_col]] == SE_ident),])[1]
            CTRL_spots <- dim(seurat_obj@meta.data[which(seurat_obj@meta.data[[cluster_condition_col]] == CTRL_ident),])[1]
            if (SE_spots > 3 & CTRL_spots > 3) {
                cluster_out_dir <- sprintf("%s/c%s", out_dir, one_cluster)
                if (! dir.exists(cluster_out_dir)) {
                    dir.create(cluster_out_dir, recursive=TRUE, mode="0775")
                }
                cluster_out_name <- sprintf("%s_c%s", out_name, one_cluster)
                if (test == "negbinom") {
                    mean_function <- function(x, pseudocount.use=1, base=2) {
                        return(log(x = rowMeans(x = x) + pseudocount.use, base = base))
                    }
                    seurat_obj.DE.one_cluster <- FindMarkers(seurat_obj, ident.1=SE_ident, ident.2=CTRL_ident, logfc.threshold=0, test.use=test, latent.vars=latent_vars, mean.fxn=mean_function, recorrect_umi=recorrect, verbose=FALSE)
                } else {
                    seurat_obj.DE.one_cluster <- FindMarkers(seurat_obj, ident.1=SE_ident, ident.2=CTRL_ident, logfc.threshold=0, test.use=test, latent.vars=latent_vars, recorrect_umi=recorrect, verbose=FALSE)
                }
                write.csv(seurat_obj.DE.one_cluster, file=sprintf("%s/%s.csv", cluster_out_dir, cluster_out_name), quote=FALSE, row.names=TRUE)
                seurat_obj.DE.all_clusters <- rbind(seurat_obj.DE.all_clusters, cbind(cluster=rep(one_cluster, dim(seurat_obj.DE.one_cluster)[1]), gene=rownames(seurat_obj.DE.one_cluster), data.frame(seurat_obj.DE.one_cluster, row.names=NULL)))
                #### visualizations
                pdf(sprintf("%s/%s.pdf", cluster_out_dir, cluster_out_name))
                ##### volcano plots
                one_fc_threshold <- 0
                volcano_plots(seurat_obj.DE.one_cluster, one_fc_threshold, sprintf("Cluster: %s", one_cluster), cluster_out_dir, sprintf("%s_%s", cluster_out_name, one_fc_threshold))
                one_fc_threshold <- 0.25
                volcano_plots(seurat_obj.DE.one_cluster, one_fc_threshold, sprintf("Cluster: %s", one_cluster), cluster_out_dir, sprintf("%s_%s", cluster_out_name, sub("[.]", "_", one_fc_threshold)))              
                ##### violin plots
                for(i in seq(1, nb_plots, plots_per_page)) {
                    print(FeaturePlot(seurat_obj, features=rownames(seurat_obj.DE.one_cluster)[i:(i+plots_per_page-1)], split.by=condition_col, max.cutoff=3, cols=c("grey", "red")))
                    plots <- VlnPlot(seurat_obj, features=rownames(seurat_obj.DE.one_cluster)[i:(i+plots_per_page-1)], split.by=condition_col, group.by=cluster_col, pt.size=0, combine=FALSE)
                    print(wrap_plots(plots=plots, ncol=1))
                    plots <- VlnPlot(seurat_obj, features=rownames(seurat_obj.DE.one_cluster)[i:(i+plots_per_page-1)], split.by="orig.ident", group.by=cluster_col, pt.size=0, combine=FALSE)
                    print(wrap_plots(plots=plots, ncol=1))
                }
                ##### tissue sections
                genes2plot <- rownames(seurat_obj.DE.one_cluster)[1:10]
                nb_sections <- length(levels(seurat_obj@meta.data$orig.ident))
                if (nb_sections > 4) {
                    # TODO
                    # one page per timepoint
                    #timepoints <- levels(seurat_obj@meta.data$time)
                    #for (one_timepoint in timepoints) {
                    #    one_timepoint_metadata_df <- seurat_obj_int@meta.data[which(seurat_obj_int@meta.data$time==one_timepoint),]
                    #    one_timepoint_metadata_df <- droplevels(one_timepoint_metadata_df)
                    #    sections <- levels(one_timepoint_metadata_df$orig.ident)
                    #    print(SpatialDimPlot(seurat_obj_int, images=sections, ncol=2, label=TRUE, label.size=2) + plot_annotation(title=sprintf("%s\ntimepoint: %s", plot_title, one_timepoint), theme=theme(plot.title=element_text(size=18))))
                    #}
                } else {
                    for (one_gene in genes2plot) {
                        one_gene_pval <- seurat_obj.DE.one_cluster[one_gene, "p_val"]
                        one_gene_fc <- seurat_obj.DE.one_cluster[one_gene, "avg_log2FC"]
                        print(SpatialFeaturePlot(object=seurat_obj, features=one_gene, alpha=c(0.1,1), ncol=2) + plot_annotation(title=sprintf("Cluster: %s\np-val: %.2e, FC: %.2f", one_cluster, one_gene_pval, one_gene_fc), theme=theme(plot.title=element_text(size=16))))
                    }
                }
                dev.off()
                #### GO analysis
                #GO_analysis(seurat_obj.DE.one_cluster, 0, cluster_out_dir, cluster_out_name)
            }
        }
    }
    write.csv(seurat_obj.DE.all_clusters, file=sprintf("%s/%s_all_clusters.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
    return(seurat_obj.DE.all_clusters)
}

identify_DEGs <- function(seurat_obj_int, assay2use, cluster_col, condition_col, timepoint_analysis, test, latent_vars, nb_plots, plots_per_page, out_dir, out_name) {
    print("Identify DEGs")
    # set default assay: must be either 'RNA' or 'SCT'
    DefaultAssay(seurat_obj_int) <- assay2use
    # average expression plots
    print("average expression plots")
    theme_set(theme_cowplot())
    ## set identity classes to clusters
    Idents(seurat_obj_int) <- seurat_obj_int@meta.data[[cluster_col]]
    seurat_obj_int.all_clusters.avg <- data.frame()
    pdf(sprintf("%s/%s_avg_exp.pdf", out_dir, out_name))
    for (one_cluster in levels(seurat_obj_int@meta.data[[cluster_col]])) {
        print(sprintf("cluster: %s", one_cluster))
        seurat_obj_int.one_cluster <- subset(seurat_obj_int, idents=one_cluster)
        ## set identity classes to conditions
        Idents(seurat_obj_int.one_cluster) <- condition_col
        if ("SE" %in% levels(Idents(seurat_obj_int.one_cluster)) & "CTRL" %in% levels(Idents(seurat_obj_int.one_cluster))) {
            seurat_obj_int.one_cluster.avg <- as.data.frame(log1p(AverageExpression(seurat_obj_int.one_cluster, verbose=FALSE)[[assay2use]]))
            seurat_obj_int.one_cluster.avg$gene <- rownames(seurat_obj_int.one_cluster.avg)
            seurat_obj_int.one_cluster.avg$diff <- seurat_obj_int.one_cluster.avg$SE - seurat_obj_int.one_cluster.avg$CTRL
            seurat_obj_int.all_clusters.avg <- rbind(seurat_obj_int.all_clusters.avg, cbind(seurat_obj_int.one_cluster.avg, cluster=rep(one_cluster, dim(seurat_obj_int.one_cluster.avg)[1])))
            genes.to.label <- c(rownames(seurat_obj_int.one_cluster.avg[order(seurat_obj_int.one_cluster.avg$diff, decreasing=TRUE),])[1:10], tail(rownames(seurat_obj_int.one_cluster.avg[order(seurat_obj_int.one_cluster.avg$diff, decreasing=TRUE),]), n=10))
            p <- ggplot(seurat_obj_int.one_cluster.avg, aes(CTRL, SE)) + geom_point() + ggtitle(sprintf("Cluster %s", one_cluster))
            p <- LabelPoints(p, points=genes.to.label, repel=TRUE)
            print(p) 
        }
    }
    dev.off()
    write.csv(seurat_obj_int.all_clusters.avg, file=sprintf("%s/%s_avg_exp.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)

    # differentially expressed genes in different conditions for cells of the same type
    print("differentially expressed genes between conditions")
    latent_vars_param <- ifelse(length(latent_vars) == 1, switch(as.integer(is.na(latent_vars))+1, latent_vars, NULL), latent_vars)
    ## column in metadata table to hold both the cell type and stimulation information and switch the current ident to that column
    seurat_obj_int[[sprintf("%s.%s", cluster_col, condition_col)]] <- paste(Idents(seurat_obj_int), seurat_obj_int@meta.data[[condition_col]], sep="_")
    if (timepoint_analysis) {
        # analysis over all timepoints
        #over_timepoints_latent_vars <- c(latent_vars, "time")
        #over_timepoints_out_dir <- sprintf("%s/all", out_dir)
        #if (! dir.exists(over_timepoints_out_dir)) {
        #    dir.create(over_timepoints_out_dir, recursive=TRUE, mode="0775")
        #}
        #all_clusters_DE_results_df <- identify_DE_markers(seurat_obj_int, cluster_col, condition_col, test, over_timepoints_latent_vars, TRUE, nb_plots, plots_per_page, over_timepoints_out_dir, out_name)
        # per-timepoint analysis
        all_times_all_clusters_DE_results_df <- data.frame()
        for (one_time in levels(seurat_obj_int@meta.data$time)) {
            barcodes2use <- rownames(seurat_obj_int@meta.data[which(seurat_obj_int@meta.data$time == one_time),])
            ## subset
            seurat_obj_int_subset <- subset(seurat_obj_int, cells=barcodes2use)
            timepoint_out_dir <- sprintf("%s/D%s", out_dir, one_time)
            if (! dir.exists(timepoint_out_dir)) {
                dir.create(timepoint_out_dir, recursive=TRUE, mode="0775")
            }
            ### after subsetting, only keep image of the samples in the subset
            seurat_obj_int_subset@meta.data <- droplevels(seurat_obj_int_subset@meta.data)
            samples <- levels(seurat_obj_int_subset@meta.data$orig.ident)
            seurat_obj_int_subset@images <- seurat_obj_int_subset@images[samples]
            timepoint_out_name <- sprintf("%s_D%s", out_name, one_time)
            ## identify DE genes
            all_clusters_DE_results_df <- identify_DE_markers(seurat_obj_int_subset, cluster_col, condition_col, test, latent_vars, FALSE, nb_plots, plots_per_page, timepoint_out_dir, timepoint_out_name)
            all_times_all_clusters_DE_results_df <- rbind(all_times_all_clusters_DE_results_df, cbind(time=rep(one_time, dim(all_clusters_DE_results_df)[1]), all_clusters_DE_results_df))
        }
        timepoint_comparison_out_dir <- sprintf("%s/comp", out_dir)
        if (! dir.exists(timepoint_comparison_out_dir)) {
            dir.create(timepoint_comparison_out_dir, recursive=TRUE, mode="0775")
        }
        ## comparisons across timepoints
        comp_out_name <- sprintf("%s_all_times_all_clusters", out_name)
        write.csv(all_times_all_clusters_DE_results_df, file=sprintf("%s/%s.csv", timepoint_comparison_out_dir, comp_out_name), quote=FALSE, row.names=FALSE)
        all_times_all_clusters_DE_results_only_DEGs_df <- all_times_all_clusters_DE_results_df[which(all_times_all_clusters_DE_results_df$p_val_adj < 0.05),]
        ### DE counts per time and cluster
        up_DEG_count_df <- aggregate(all_times_all_clusters_DE_results_only_DEGs_df$avg_log2FC, by=list(all_times_all_clusters_DE_results_only_DEGs_df$time, all_times_all_clusters_DE_results_only_DEGs_df$cluster), function(x) { length(x[x > 0]) })
        colnames(up_DEG_count_df) <- c("time", "cluster", "count")
        up_DEG_count_df <- cbind(category=rep("up", dim(up_DEG_count_df)[1]), up_DEG_count_df)
        up_DEG_count_df$category <- as.factor(up_DEG_count_df$category)
        up_DEG_count_df$time <- factor(up_DEG_count_df$time, levels=sort(as.numeric(unique(up_DEG_count_df$time))))
        up_DEG_count_df$cluster <- as.factor(up_DEG_count_df$cluster)
        down_DEG_count_df <- aggregate(all_times_all_clusters_DE_results_only_DEGs_df$avg_log2FC, by=list(all_times_all_clusters_DE_results_only_DEGs_df$time, all_times_all_clusters_DE_results_only_DEGs_df$cluster), function(x) { length(x[x < 0]) })
        colnames(down_DEG_count_df) <- c("time", "cluster", "count")
        down_DEG_count_df <- cbind(category=rep("down", dim(down_DEG_count_df)[1]), down_DEG_count_df)
        down_DEG_count_df$category <- as.factor(down_DEG_count_df$category)
        down_DEG_count_df$time <- factor(down_DEG_count_df$time, levels=sort(as.numeric(unique(down_DEG_count_df$time))))
        down_DEG_count_df$cluster <- as.factor(down_DEG_count_df$cluster)
        DE_count_df <- rbind(up_DEG_count_df, down_DEG_count_df)
        write.csv(all_times_all_clusters_DE_results_df, file=sprintf("%s/%s_DEG_counts.csv", timepoint_comparison_out_dir, comp_out_name), quote=FALSE, row.names=FALSE)
        pdf(sprintf("%s/%s_DEG_counts.pdf", timepoint_comparison_out_dir, comp_out_name))
        p <- ggplot(DE_count_df, aes(x=time, y=count, fill=cluster)) +
            geom_bar(stat="identity", position="dodge") +
            facet_grid(~category) +
            labs(title="DE genes per time and cluster", x="Time", y="Count", fill="Cluster") +
            theme_bw() +
            theme(legend.position="bottom") +
            theme(panel.border=element_rect(color="grey50"))
        print(p)
        dev.off()
    } else {
        all_clusters_DE_results_df <- identify_DE_markers(seurat_obj_int, cluster_col, condition_col, test, latent_vars, TRUE, nb_plots, plots_per_page, out_dir, out_name)
    }
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

#prep_marker_analysis <- function (seurat_obj_int, mito_genes, hb_genes, plots_per_page, out_dir, out_name) {
#    seurat_obj_int <- PrepSCTFindMarkers(seurat_obj_int)
#    feature_stat_slot <- "counts"
#    seurat_obj_int <- addFeatureStats2Seurat(seurat_obj_int, "SCT", feature_stat_slot)
#    norm_column_prefix <- "prep_sct_norm"
#    seurat_obj_int <- addPercentageFeatureSet2Seurat(seurat_obj_int, mito_genes, hb_genes, "SCT", norm_column_prefix)
#    library(RColorBrewer)
#    seurat_all_feature_stat_plots(seurat_obj_int, "SCT", norm_column_prefix, plots_per_page, "Normalized counts - prep SCT markers", out_dir, sprintf("%s_prep_markers", out_name))
#    detach("package:RColorBrewer")
#    return(seurat_obj_int)
#}

marker_analysis <- function(seurat_obj_int, assay2use, mito_genes, hb_genes, test, latent_vars, clusters, condition_subset, scale_split, scale_regress, load_image, nb_plots, plots_per_page, DEGs, timepoint_analysis, DEG_nb_plots, DEG_plots_per_page, conserved, cluster_col, condition_col, out_dir, out_name) {
    # prepare for marker identification
    if (clusters | DEGs | conserved) {
        if (assay2use == "SCT") {
            if (length(seurat_obj_int@assays$SCT@SCTModel.list) > 1) {
                seurat_obj_int <- PrepSCTFindMarkers(seurat_obj_int)
                feature_stat_slot <- "counts"
                seurat_obj_int <- addFeatureStats2Seurat(seurat_obj_int, assay2use, feature_stat_slot)
                norm_column_prefix <- "prep_sct_norm"
                seurat_obj_int <- addPercentageFeatureSet2Seurat(seurat_obj_int, mito_genes, hb_genes, assay2use, norm_column_prefix)
            }
        }
    }
    
    # cluster markers
    if (clusters) {
        ## find markers for every cluster compared to all remaining spots, report only the positive ones
        latent_vars_in_path <- ifelse(length(latent_vars > 1), paste(gsub("[.]", "", latent_vars), collapse="-"), ifelse(is.na(latent_vars), "no_var", gsub("[.]", "", latent_vars)))
        markers_out_dir <- sprintf("%s/00-Markers/%s/%s/%s", out_dir, test, latent_vars_in_path, sub("seurat_", "", cluster_col))
        if (! dir.exists(markers_out_dir)) {
            dir.create(markers_out_dir, recursive=TRUE, mode="0775")
        }
        library(RColorBrewer)
        seurat_all_feature_stat_plots(seurat_obj_int, assay2use, norm_column_prefix, plots_per_page, "Normalized counts - prep SCT markers", markers_out_dir, sprintf("%s_prep_markers", out_name))
        detach("package:RColorBrewer")
        identify_cluster_markers(seurat_obj_int, assay2use, condition_col, condition_subset, cluster_col, scale_split, scale_regress, test, latent_vars, load_image, nb_plots, markers_out_dir, sprintf("%s_cluster_marker_genes", out_name))
    }

    # differentially expressed genes across conditions
    if (DEGs) {
        latent_vars_in_path <- ifelse(length(latent_vars > 1), paste(gsub("[.]", "", latent_vars), collapse="-"), ifelse(is.na(latent_vars), "no_var", gsub("[.]", "", latent_vars)))
        DEG_out_dir <- sprintf("%s/10-DEGs/%s/%s", out_dir, test, latent_vars_in_path)
        if (! dir.exists(DEG_out_dir)) {
            dir.create(DEG_out_dir, recursive=TRUE, mode="0775")
        }
        library(RColorBrewer)
        seurat_all_feature_stat_plots(seurat_obj_int, assay2use, norm_column_prefix, plots_per_page, "Normalized counts - prep SCT markers", DEG_out_dir, sprintf("%s_prep_markers", out_name))
        detach("package:RColorBrewer")
        identify_DEGs(seurat_obj_int, assay2use, cluster_col, condition_col, timepoint_analysis, test, latent_vars, DEG_nb_plots, DEG_plots_per_page, DEG_out_dir, sprintf("%s_DEGs", out_name))
    }

    # conserved cell type markers
    if (conserved) {
        conserved_markers_out_dir <- sprintf("%s/20-Conserved_markers", out_dir)
        if (! dir.exists(conserved_markers_out_dir)) {
            dir.create(conserved_markers_out_dir, recursive=TRUE, mode="0775")
        }
        library(RColorBrewer)
        seurat_all_feature_stat_plots(seurat_obj_int, assay2use, norm_column_prefix, plots_per_page, "Normalized counts - prep SCT markers", conserved_markers_out_dir, sprintf("%s_prep_markers", out_name))
        detach("package:RColorBrewer")
        identify_conserved_markers(seurat_obj_int, assay2use, cluster_col, condition_col, conserved_markers_out_dir, sprintf("%s_conserved_markers", out_name))
    }
}

#################
# Subclustering #
#################

clustering_top_params <- function(eval_df, eval_col, method) {
    eval_subdf <- eval_df[which(eval_df$normalization=="SCTransform_v2" & eval_df$integration == method),]
    eval_subdf_order <- eval_subdf[order(eval_subdf[[eval_col]], decreasing=TRUE),]
    eval_normalization <- eval_subdf_order[1, "normalization"]
    eval_features <- eval_subdf_order[1, "features"]
    eval_integration <- eval_subdf_order[1, "integration"]
    eval_dimensions <- eval_subdf_order[1, "dimensions"]
    eval_resolution <- eval_subdf_order[1, "resolution"]
    return(list(normalization=eval_normalization, features=eval_features, integration=eval_integration, dimensions=eval_dimensions, resolution=eval_resolution))
}

graph_clustering <- function(seurat_obj_int_sub, dim_nb, red2use, res_vec, metadata_cat_vec, load_image, plots_per_page, plot_title, out_dir, out_name) {
    all_params_mean_silhouette_df <- data.frame()
    print(sprintf("dimension number: %d", dim_nb))
    seurat_object_to_cluster_sub <- FindNeighbors(seurat_obj_int_sub, reduction=red2use, dims=1:dim_nb)
    for (subclustering_resolution in res_vec) {
        print(sprintf("reclustering resolution: %.1f", subclustering_resolution))
        resolution_out_dir <- sprintf("%s/res%s", out_dir, sub("[.]", "_", subclustering_resolution))
        if (! dir.exists(resolution_out_dir)) {
            dir.create(resolution_out_dir, recursive=TRUE, mode="0775")
        }
        resolution_output_name <- sprintf("%s_res%s", out_name, sub("[.]", "_", subclustering_resolution))

        seurat_object_clustering_sub <- clustering(seurat_object_to_cluster_sub, subclustering_resolution, red2use, dim_nb)
        write.csv(seurat_object_clustering_sub@meta.data, file=sprintf("%s/%s_metadata.csv", resolution_out_dir, resolution_output_name), quote=FALSE, row.names=TRUE)
        mean_silhouette <- mean(seurat_object_clustering_sub$silhouette_seurat_clusters)
        nb_clusters <- length(levels(seurat_object_clustering_sub@meta.data$seurat_clusters))
        ### plot clusters onto UMAP or onto the tissue section
        clustering_plots(seurat_object_clustering_sub, metadata_cat_vec, load_image, plots_per_page, plot_title, resolution_out_dir, sprintf("%s_clustering", resolution_output_name))

        ### get mean silhouette scores: all spots and hippocampus spots
        mean_silhouette_df <- data.frame(dimensions=dim_nb, resolution=subclustering_resolution, nb_clusters_sub=nb_clusters, mean_silhouette_seurat_clusters_sub=mean_silhouette)
        all_params_mean_silhouette_df <- rbind(all_params_mean_silhouette_df, mean_silhouette_df)

        rm(seurat_object_clustering_sub)
        gc()
    }
    rm(seurat_object_to_cluster_sub)
    gc()
    return(all_params_mean_silhouette_df)
}

norm_dimred_clustering <- function(seurat_obj_int_sub, norm_method, feat_vec, var2regress, only_var, mito_genes, hb_genes, integration, integration_method, harmony_seed, integration_var, dim_vec, metadata_cat_vec, res_vec, load_image, plots_per_page, plot_title, out_dir) {
    project_name <- seurat_obj_int_sub@project.name
    # split integrated object into a list of objects per sample
    seurat_objects_sub <- SplitObject(seurat_obj_int_sub, split.by=integration_var)
    ## after splitting, all the objects still have all the images: only keep the image of the corresponding sample
    for (one_sample in names(seurat_objects_sub)) {
        seurat_objects_sub[[one_sample]]@images <- list(seurat_objects_sub[[one_sample]]@images[[one_sample]])
        names(seurat_objects_sub[[one_sample]]@images) <- one_sample
    }
    if (norm_method %in% c("SCTransform", "SCTransform_v2")) {
        ## after splitting, all the objects have 4 models in SCTModel.list slot but only one of them has cells: only keep this one
        for (one_sample in names(seurat_objects_sub)) {
            for (one_model in names(seurat_objects_sub[[one_sample]]@assays$SCT@SCTModel.list)) {
                if (dim(seurat_objects_sub[[one_sample]]@assays$SCT@SCTModel.list[[one_model]]@cell.attributes)[1] > 0) {
                    seurat_objects_sub[[one_sample]]@assays$SCT@SCTModel.list <- list(seurat_objects_sub[[one_sample]]@assays$SCT@SCTModel.list[[one_model]])
                    names(seurat_objects_sub[[one_sample]]@assays$SCT@SCTModel.list) <- one_model
                    break
                }
            }
        }
    }

    # set default assay to RNA and project name to sample
    for (one_sample in names(seurat_objects_sub)) {
        DefaultAssay(seurat_objects_sub[[one_sample]]) <- "RNA"
        seurat_objects_sub[[one_sample]]@project.name <- one_sample
    }
    
    all_params_batch_silhouette_df <- all_params_cluster_silhouette_df <- data.frame()
    # normalization and feature selection
    for (feat_nb in feat_vec) {
        returned_list <- seurat_list_normalization(seurat_objects_sub, norm_method, var2regress, feat_nb, only_var, mito_genes, hb_genes, out_dir)
        seurat_objects_norm <- returned_list[["seurat_list"]]
        integration_features_nb <- returned_list[["features_nb"]]
        integration_normalization_method <- returned_list[["norm_method"]]
        default_assay <- returned_list[["assay"]]
        selection_method <- returned_list[["selec"]]
        normalization_output_dir <- returned_list[["out_dir"]]
        normalization_output_name <- returned_list[["out_name"]]
        rm(returned_list)
        gc()
        
        ## HVG plots
        hvg_df_list <- lapply(seurat_objects_norm, hvg_plots, assay2use=default_assay, selec=selection_method, nb_hvgs=feat_nb, plots_per_page=plots_per_page, out_dir=normalization_output_dir, out_name=sprintf("%s_HVGs", normalization_output_name))
        
        ## merge datasets, feature statistics and QC plots after normalization: for Seurat integration, merge is not necessary but still run merge_stats_plots_after_norm() for feature statistics and QC plots after normalization
        seurat_object_merge_norm <- merge_stats_plots_after_norm(seurat_objects_norm, project_name, norm_method, feat_nb, mito_genes, hb_genes, plots_per_page, normalization_output_dir, normalization_output_name)
        
        ## get variable features after subclustering
        if (integration_method == "Seurat") {
            rm(seurat_object_merge_norm)
            gc()
            variable_features <- SelectIntegrationFeatures(seurat_objects_norm, nfeatures=integration_features_nb)
        } else {
            variable_features <- VariableFeatures(seurat_object_merge_norm)
        }
        
        for (dimensions_nb_subclustering in dim_vec) {
            if (integration)  {
                dimensions_output_dir <- sprintf("%s/%s/%d_dims", normalization_output_dir, integration_method, dimensions_nb_subclustering)
                if (! dir.exists(dimensions_output_dir)) {
                    dir.create(dimensions_output_dir, recursive=TRUE, mode="0775")
                }
                dimensions_output_name <- sprintf("%s_%s_%d_dims", normalization_output_name, integration_method, dimensions_nb_subclustering)
                
                if (integration_method == "Seurat") {
                    if (norm_method %in% c("SCTransform", "SCTransform_v2")) {
                        seurat_objects_norm_prep <- PrepSCTIntegration(object.list=seurat_objects_norm, anchor.features=variable_features)
                    }
                    if (dimensions_nb_subclustering == dim_vec[length(dim_vec)]) {
                        # last Seurat integration: Seurat objects after normalization are not needed any more
                        rm(seurat_objects_norm)
                        gc()
                    }
                    int.anchors <- FindIntegrationAnchors(object.list=seurat_objects_norm_prep, anchor.features=variable_features, normalization.method=integration_normalization_method, dims=1:dimensions_nb_subclustering)
                    rm(seurat_objects_norm_prep)
                    gc()
                    
                    ### get minimum number of spots per sample
                    min_spot_nb <- min(table(seurat_obj_int_sub@meta.data$orig.ident))
                    if (min_spot_nb < 100) { # IntegrateData() k.weight parameter default value
                        seurat_obj_to_dimred <- IntegrateData(anchorset=int.anchors, normalization.method=integration_normalization_method, dims=1:dimensions_nb_subclustering, k.weight=min_spot_nb)
                    } else {
                        seurat_obj_to_dimred <- IntegrateData(anchorset=int.anchors, normalization.method=integration_normalization_method, dims=1:dimensions_nb_subclustering)
                    }
                    rm(int.anchors)
                    gc()
                    
                    ## set project name to output dirname
                    seurat_obj_to_dimred@project.name <- project_name
                    ## orig.ident, lame, zone, condition and time are not factors in the integrated object meta data
                    seurat_obj_to_dimred@meta.data$orig.ident <- as.factor(seurat_obj_to_dimred@meta.data$orig.ident)
                    seurat_obj_to_dimred@meta.data$lame <- as.factor(seurat_obj_to_dimred@meta.data$lame)
                    seurat_obj_to_dimred@meta.data$zone <- as.factor(seurat_obj_to_dimred@meta.data$zone)
                    seurat_obj_to_dimred@meta.data$condition <- as.factor(seurat_obj_to_dimred@meta.data$condition)
                    seurat_obj_to_dimred@meta.data$time <- as.factor(seurat_obj_to_dimred@meta.data$time)
                    
                    ## HVGs after integration
                    seurat_obj_to_dimred_hvgs <- VariableFeatures(seurat_obj_to_dimred)
                    seurat_obj_to_dimred_df <- seurat_obj_to_dimred@assays$integrated@SCTModel.list$refmodel@feature.attributes[seurat_obj_to_dimred_hvgs,]
                    write.csv(seurat_obj_to_dimred_df, file=sprintf("%s/%s_%s_HVGs_residual_mean_variance_after_integration.csv", dimensions_output_dir, seurat_obj_to_dimred@project.name, dimensions_output_name), quote=FALSE, row.names=TRUE)
                    
                    reduction2use <- "pca"
                } else {
                    if (integration_method == "Harmony") {
                        seurat_obj_to_dimred <- seurat_object_merge_norm
                    }
                }
            } else {
                dimensions_output_dir <- sprintf("%s/%d_dims", normalization_output_dir, dimensions_nb_subclustering)
                if (! dir.exists(dimensions_output_dir)) {
                    dir.create(dimensions_output_dir, recursive=TRUE, mode="0775")
                }
                dimensions_output_name <- sprintf("%s_%d_dims", normalization_output_name, dimensions_nb_subclustering)
                seurat_obj_to_dimred <- seurat_obj_int_sub
            }
            
            # dimension reduction
            seurat_obj_dimred <- RunPCA(seurat_obj_to_dimred, features=variable_features, verbose=TRUE)
            rm(seurat_obj_to_dimred)
            gc()
            ## cell embeddings, feature loadings and dimension standard deviation
            dimension_reduction_outputs(seurat_obj_dimred@reductions, "pca", dimensions_output_dir, sprintf("%s_PC", dimensions_output_name))
            
            if (integration)  {
                if (integration_method == "Harmony") {
                    ##### dims.use parameter of RunHarmony() function is not used: build a new DimReduc object with the desired number of dimensions
                    ###### https://github.com/immunogenomics/harmony/issues/82
                    ###### https://github.com/immunogenomics/harmony/issues/151
                    seurat_obj_dimred@reductions$pca_dims <- seurat_obj_dimred@reductions$pca
                    seurat_obj_dimred@reductions$pca_dims@cell.embeddings <- seurat_obj_dimred@reductions$pca_dims@cell.embeddings[, 1:dimensions_nb_subclustering]
                    ##### Harmony results are not reproducible: set seed
                    ###### https://github.com/immunogenomics/harmony/issues/13
                    ###### https://github.com/immunogenomics/harmony/issues/56    
                    set.seed(harmony_seed)
                    seurat_obj_dimred <- RunHarmony(seurat_obj_dimred, group.by.vars=integration_var, reduction="pca_dims")
                    reduction2use <- "harmony"
                    ## cell embeddings, feature loadings and dimension standard deviation
                    dimension_reduction_outputs(seurat_obj_dimred@reductions, reduction2use, dimensions_output_dir, sprintf("%s_%s", dimensions_output_name, reduction2use))
                }
                
                ## compute silhouette scores for lames, zones, conditions and samples
                mean_silhouette_vector <- c()
                for (metadata_category_colname in metadata_cat_vec) {
                    seurat_obj_dimred <- compute_silhouette(seurat_obj_dimred, reduction2use, dimensions_nb_subclustering, metadata_category_colname)
                    mean_silhouette <- mean(seurat_obj_dimred@meta.data[[sprintf("silhouette_%s", metadata_category_colname)]])
                    mean_silhouette_vector <- c(mean_silhouette_vector, mean_silhouette)
                }
                mean_silhouette_df <- as.data.frame(t(mean_silhouette_vector))
                colnames(mean_silhouette_df) <- sprintf("mean_silhouette_%s", metadata_cat_vec)
                params_mean_silhouette_df <- cbind(normalization=norm_method, features=feat_nb, integration=integration_method, dimensions=dimensions_nb_subclustering, mean_silhouette_df)
                all_params_batch_silhouette_df <- rbind(all_params_batch_silhouette_df, params_mean_silhouette_df)
            } else {
                # set reduction2use if integration is NULL, i.e. no integration is done and the reduction to use is the one done during the whole clustering
                if (integration_method == "Seurat") {
                    reduction2use <- "pca"
                } else {
                    if (integration_method == "Harmony") {
                        reduction2use <- "harmony"
                    }
                }
            }
            
            # graph and clustering
            silhouette_df <- graph_clustering(seurat_obj_dimred, dimensions_nb_subclustering, reduction2use, res_vec, metadata_cat_vec, load_image, plots_per_page, plot_title, dimensions_output_dir, dimensions_output_name)
            all_params_cluster_silhouette_df <- rbind(all_params_cluster_silhouette_df, cbind(data.frame(normalization=rep(norm_method, dim(silhouette_df)[1]), features=rep(feat_nb, dim(silhouette_df)[1]), integration=rep(integration_method, dim(silhouette_df)[1])), silhouette_df))
            rm(seurat_obj_dimred)
            gc()
        }
        rm(seurat_objects_norm)
        rm(seurat_object_merge_norm)
        gc()   
    }
    
    return(list(batch=all_params_batch_silhouette_df, cluster=all_params_cluster_silhouette_df))
}


###################
# Entire workflow #
###################

norm2markers_workflow <- function(seurat_obj_list, param_list, out_dir) {
    
    # get parameters
    project_name <- param_list$project_name
    normalization <- param_list$normalization_method
    normalization_var2regress <- param_list$normalization_variables_to_regress
    features <- param_list$feature_nb
    mito_genes <- param_list$mitochondrial_genes
    hb_genes <- param_list$hemoglobin_genes
    integration <- param_list$integration_method
    dimensions <- param_list$dimension_nb
    scale_var2regress <- param_list$scale_variables_to_regress
    sample_variable <- param_list$sample_variable
    harmony_seed <- param_list$harmony_seed
    metadata_categories <- param_list$metadata_categories
    resolution <- param_list$resolution
    plots_per_page <- param_list$plots_per_page
    markers_test <- param_list$markers_test
    markers_latent_vars <- param_list$markers_latent_variables
    cluster_markers <- param_list$markers_cluster
    cluster_markers_condition <- param_list$markers_cluster_condition
    cluster_markers_plot_nb <- param_list$markers_cluster_plot_nb
    DE_analysis <- param_list$markers_DE
    timepoint_analysis <- param_list$markers_DE_timepoint_analysis
    nb_DEGs_to_plot <- param_list$markers_DE_plot_nb
    nb_DEGs_per_page <- param_list$markers_DE_plots_per_page
    conserved_markers <- param_list$markers_conserved
    markers_condition <- param_list$markers_condition
    
    # workflow
    print(sprintf("normalization method: %s", normalization))
    return_only_var_genes <- ifelse(normalization %in% c("SCTransform", "SCTransform_v2"), FALSE, NA)

    ## normalization and feature selection
    returned_list <- seurat_list_normalization(seurat_obj_list, normalization, normalization_var2regress, features, return_only_var_genes, mito_genes, hb_genes, out_dir)
    seurat_objects_norm <- returned_list[["seurat_list"]]
    integration_features_nb <- returned_list[["features_nb"]]
    integration_normalization_method <- returned_list[["norm_method"]]
    default_assay <- returned_list[["assay"]]
    selection_method <- returned_list[["selec"]]
    normalization_out_dir <- returned_list[["out_dir"]]
    normalization_out_name <- returned_list[["out_name"]]
    rm(returned_list)
    gc()

    ### HVG plots
    hvg_df_list <- lapply(seurat_objects_norm, hvg_plots, assay2use=default_assay, selec=selection_method, nb_hvgs=features, plots_per_page=NA, out_dir=normalization_out_dir, out_name=sprintf("%s_HVGs", normalization_out_name))

    ### merge datasets, feature statistics and QC plots after normalization
    seurat_object_merge_norm <- merge_stats_plots_after_norm(seurat_objects_norm, project_name, normalization, features, mito_genes, hb_genes, plots_per_page, normalization_out_dir, normalization_out_name)

    ## integration
    integration_out_dir <- sprintf("%s/%s", normalization_out_dir, gsub(" ", "_", integration))
    if (! dir.exists(integration_out_dir)) {
        dir.create(integration_out_dir, recursive=TRUE, mode="0775")
    }
    ### dimensionality reduction
    print(sprintf("dimension number: %d", dimensions))
    dimensions_out_dir <- sprintf("%s/%d_dims", integration_out_dir, dimensions)
    if (! dir.exists(dimensions_out_dir)) {
        dir.create(dimensions_out_dir, recursive=TRUE, mode="0775")
    }
    dimensions_out_name <- sprintf("%s_%s_%ddims", normalization_out_name, gsub(" ", "_", integration), dimensions)

    if (integration == "Seurat") {
        rm(seurat_object_merge_norm)
        gc()

        ### Seurat integration
        features <- SelectIntegrationFeatures(seurat_objects_norm, nfeatures=integration_features_nb)
        if (normalization %in% c("SCTransform", "SCTransform_v2")) {
            seurat_objects_norm <- PrepSCTIntegration(object.list=seurat_objects_norm, anchor.features=features)
        }
        int.anchors <- FindIntegrationAnchors(object.list=seurat_objects_norm, anchor.features=features, normalization.method=integration_normalization_method, dims=1:dimensions)
        
        #### get minimum number of spots per sample
        min_spot_nb <- min(unlist(lapply(seurat_objects_norm, function(x) { dim(x@meta.data)[1] })))
        rm(seurat_objects_norm)
        gc()
        if (min_spot_nb < 100) { # IntegrateData() k.weight parameter default value
            seurat_object_integration <- IntegrateData(anchorset=int.anchors, normalization.method=integration_normalization_method, dims=1:dimensions, k.weight=min_spot_nb)
        } else {
            seurat_object_integration <- IntegrateData(anchorset=int.anchors, normalization.method=integration_normalization_method, dims=1:dimensions)
        }
        rm(int.anchors)
        gc()
        #### set project name to output dirname
        seurat_object_integration@project.name <- project_name
        #### orig.ident, lame, zone, condition and time are not factors in the integrated object meta data
        seurat_object_integration@meta.data$orig.ident <- as.factor(seurat_object_integration@meta.data$orig.ident)
        seurat_object_integration@meta.data$lame <- as.factor(seurat_object_integration@meta.data$lame)
        seurat_object_integration@meta.data$zone <- as.factor(seurat_object_integration@meta.data$zone)
        seurat_object_integration@meta.data$condition <- as.factor(seurat_object_integration@meta.data$condition)
        seurat_object_integration@meta.data$time <- as.factor(seurat_object_integration@meta.data$time)

        #### HVGs after integration
        seurat_object_integration_hvgs <- VariableFeatures(seurat_object_integration)
        seurat_object_integration_df <- seurat_object_integration@assays$integrated@SCTModel.list$refmodel@feature.attributes[seurat_object_integration_hvgs,]
        write.csv(seurat_object_integration_df, file=sprintf("%s/%s_%s_HVGs_residual_mean_variance_after_merge.csv", dimensions_out_dir, seurat_object_integration@project.name, dimensions_out_name), quote=FALSE, row.names=TRUE)
    } else {
        seurat_object_integration <- seurat_object_merge_norm
        rm(seurat_objects_norm)
        rm(seurat_object_merge_norm)
        gc()
    }

    ### scaling
    if (normalization == "LogNormalize") {
        if (! is.na(sample_variable)) {
            if (!is.na(scale_var2regress)) {
                seurat_object_integration <- ScaleData(seurat_object_integration, split.by=sample_variable, vars.to.regress=scale_var2regress)
                dimensions_out_dir <- sprintf("%s/%s_%s", dimensions_out_dir, gsub("[.]", "", sample_variable), paste(gsub("[.]", "", scale_var2regress), collapse="-"))
            } else {
                seurat_object_integration <- ScaleData(seurat_object_integration, split.by=sample_variable)
                dimensions_out_dir <- sprintf("%s/%s", dimensions_out_dir, gsub("[.]", "", sample_variable))
            }
        } else {
            seurat_object_integration <- ScaleData(seurat_object_integration)
        }
    }

    ### dimension reduction
    seurat_object_integration <- RunPCA(seurat_object_integration, verbose=FALSE)
    #### cell embeddings, feature loadings and dimension standard deviation
    dimension_reduction_outputs(seurat_object_integration@reductions, "pca", dimensions_out_dir, sprintf("%s_PC", dimensions_out_name))
    
    #### Harmony integration and dimensionality reduction
    if (integration == "Harmony") {
        ##### dims.use parameter of RunHarmony() function is not used: build a new DimReduc object with the desired number of dimensions
        ###### https://github.com/immunogenomics/harmony/issues/82
        ###### https://github.com/immunogenomics/harmony/issues/151
        seurat_object_integration@reductions$pca_dims <- seurat_object_integration@reductions$pca
        seurat_object_integration@reductions$pca_dims@cell.embeddings <- seurat_object_integration@reductions$pca_dims@cell.embeddings[, 1:dimensions]
        ##### Harmony results are not reproducible: set seed
        ###### https://github.com/immunogenomics/harmony/issues/13
        ###### https://github.com/immunogenomics/harmony/issues/56    
        set.seed(harmony_seed)
        seurat_object_integration <- RunHarmony(seurat_object_integration, group.by.vars=sample_variable, reduction="pca_dims")
        reduction2use <- "harmony"
        #### cell embeddings, feature loadings and dimension standard deviation
        dimension_reduction_outputs(seurat_object_integration@reductions, reduction2use, dimensions_out_dir, sprintf("%s_%s", dimensions_out_name, reduction2use))
    } else {
        reduction2use <- "pca"
    }

    ### compute silhouette scores for lames, zones, conditions and samples
    for (metadata_category_colname in metadata_categories) {
        seurat_object_integration <- compute_silhouette(seurat_object_integration, reduction2use, dimensions, metadata_category_colname)
    }
    
    ## clustering
    seurat_object_to_cluster <- FindNeighbors(seurat_object_integration, reduction=reduction2use, dims=1:dimensions)
    rm(seurat_object_integration)
    gc()
    
    print(sprintf("cluster resolution: %.1f", resolution))
    resolution_out_dir <- sprintf("%s/res%s", dimensions_out_dir, sub("[.]", "_", resolution))
    if (! dir.exists(resolution_out_dir)) {
        dir.create(resolution_out_dir, recursive=TRUE, mode="0775")
    }
    resolution_out_name <- sprintf("%s_res%s", dimensions_out_name, sub("[.]", "_", resolution))

    seurat_object_clustering <- clustering(seurat_object_to_cluster, resolution, reduction2use, dimensions)
    write.csv(seurat_object_clustering@meta.data, file=sprintf("%s/%s_metadata.csv", resolution_out_dir, resolution_out_name), quote=FALSE, row.names=TRUE)
    rm(seurat_object_to_cluster)
    gc()
    ### plot clusters onto UMAP or onto the tissue section
    clustering_plots(seurat_object_clustering, metadata_categories, TRUE, plots_per_page, sprintf("Integrated method: %s", integration), resolution_out_dir, sprintf("%s_clustering", resolution_out_name))
    
    #if (cluster_markers | DE_analysis | conserved_markers) {
    #    if (default_assay == "SCT") {
    #        if (length(seurat_object_clustering@assays$SCT@SCTModel.list) > 1) {
    #            seurat_object_clustering <- prep_marker_analysis(seurat_object_clustering, mito_genes, hb_genes, plots_per_page, resolution_out_dir, resolution_out_name)
    #        }
    #    }
    #}
    
    ### cluster markers, differentially expressed genes across conditions and conserved cell type markers
    marker_analysis(seurat_object_clustering, default_assay, mito_genes, hb_genes, markers_test, markers_latent_vars, cluster_markers, cluster_markers_condition, sample_variable, scale_var2regress, TRUE, cluster_markers_plot_nb, plots_per_page, DE_analysis, timepoint_analysis, nb_DEGs_to_plot, nb_DEGs_per_page, conserved_markers, "seurat_clusters", markers_condition, resolution_out_dir, resolution_out_name)
    
    return(seurat_object_clustering)
}


#################################
# GO analysis: compareCluster() #
#################################

entrezid2GSEA <- function(de_entrez_df, list_var, entrez_col, fc_col, pval_col, pval_bool, out_dir, out_name) {
    entrezid_list2return <- list()
    for (one_level in levels(de_entrez_df[[list_var]])) {
        one_level_de_entrez_df <- de_entrez_df[which(de_entrez_df[[list_var]]==one_level),]
        if (pval_bool) {
            # replace null p-values
            adj_pval_no_null <- one_level_de_entrez_df[[pval_col]]
            adj_pval_no_null[adj_pval_no_null  == 0] <- min(adj_pval_no_null[adj_pval_no_null != 0])
            # compute log2FC*(-log10(p-value))
            one_level_de_entrez_df$FC_pval <- one_level_de_entrez_df[[fc_col]] * (-log(adj_pval_no_null, 10))
            # order
            one_level_de_entrez_df_order <- one_level_de_entrez_df[order(-one_level_de_entrez_df$FC_pval),]
            val2GSEA <- sprintf("%s_%s", fc_col, pval_col)
            colnames(one_level_de_entrez_df_order)[colnames(one_level_de_entrez_df_order) == "FC_pval"] <- val2GSEA
        } else {
            # order
            one_level_de_entrez_df_order <- one_level_de_entrez_df[order(-one_level_de_entrez_df[[fc_col]]),]
            val2GSEA <- fc_col
        }
        # remove genes without Entrez ID
        one_level_de_entrez_df_order <- one_level_de_entrez_df_order[which(! is.na(one_level_de_entrez_df_order[[entrez_col]])),]
        rownames(one_level_de_entrez_df_order) <- one_level_de_entrez_df_order$ENTREZID
        one_level_de_entrez_df_order_gene_list <- one_level_de_entrez_df_order[[val2GSEA]]
        # add one_level to entrezid_list2return if at least one gene with positive value and one gene with negative value
        if (length(one_level_de_entrez_df_order_gene_list[one_level_de_entrez_df_order_gene_list > 0]) > 0 & length(one_level_de_entrez_df_order_gene_list[one_level_de_entrez_df_order_gene_list < 0]) > 0) {
            names(one_level_de_entrez_df_order_gene_list) <- one_level_de_entrez_df_order[[entrez_col]]
            entrezid_list2return[[sprintf("%s", one_level)]] <- one_level_de_entrez_df_order_gene_list
        }
    }
    # write list of Entrez IDs
    entrezid2GSEA_df <- data.frame()
    for (one_level in names(entrezid_list2return)) {
        entrezid_order <- names(entrezid_list2return[[one_level]])
        entrezid2GSEA_df <- rbind(entrezid2GSEA_df, cbind(time=rep(one_level, length(entrezid_list2return[[one_level]])), entrez=entrezid_order, FC=one_level_de_entrez_df_order[entrezid_order, fc_col], pval=one_level_de_entrez_df_order[entrezid_order, pval_col], val2GSEA=unname(entrezid_list2return[[one_level]])))
    }
    colnames(entrezid2GSEA_df)[colnames(entrezid2GSEA_df) == "FC"] <- fc_col
    colnames(entrezid2GSEA_df)[colnames(entrezid2GSEA_df) == "pval"] <- pval_col
    write.csv(entrezid2GSEA_df, file=sprintf("%s/%s.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
    
    return(entrezid_list2return)
}

compareCluster_prepare2plot <- function(compareClusterResult_obj, gsea_bool, comp_col) {
    # replace 'p.adjust' column
    compareClusterResult_obj@compareClusterResult$p.adjust.old <- compareClusterResult_obj@compareClusterResult$p.adjust
    compareClusterResult_obj@compareClusterResult$p.adjust <- -log(compareClusterResult_obj@compareClusterResult$p.adjust, 10)   
    # order terms according to min p-value
    min_pval_df <- aggregate(compareClusterResult_obj@compareClusterResult$p.adjust.old, by=list(compareClusterResult_obj@compareClusterResult$Description), min)
    colnames(min_pval_df) <- c("Description", "min_pval")
    terms2plot_order <- min_pval_df[order(min_pval_df$min_pval), "Description"]
    compareClusterResult_obj@compareClusterResult <- compareClusterResult_obj@compareClusterResult[order(match(compareClusterResult_obj@compareClusterResult$Description, terms2plot_order)),]
    if (gsea_bool) {
        # add DE category
        DE_cat <- ifelse(compareClusterResult_obj@compareClusterResult$NES > 0, "up", "down")
        compareClusterResult_obj@compareClusterResult <- cbind(compareClusterResult_obj@compareClusterResult, category=DE_cat)
        compareClusterResult_obj@compareClusterResult$category <- as.factor(compareClusterResult_obj@compareClusterResult$category)
        # rename 'Cluster' column into comp_col
        colnames(compareClusterResult_obj@compareClusterResult)[which(colnames(compareClusterResult_obj@compareClusterResult) == "Cluster")] <- comp_col
        # new Cluster column: comp_col.category
        compareClusterResult_obj@compareClusterResult <- cbind(Cluster=sprintf("%s.%s", compareClusterResult_obj@compareClusterResult[[comp_col]], compareClusterResult_obj@compareClusterResult$category), compareClusterResult_obj@compareClusterResult)
        compareClusterResult_obj@compareClusterResult$Cluster <- as.factor(compareClusterResult_obj@compareClusterResult$Cluster)
    } else {
        # GO ORA: add comp_col and 'category' columns if needed
        if (! comp_col %in% colnames(compareClusterResult_obj@compareClusterResult)) {
            # create a comp_col column based on the 'Cluster' column which has the "comp_col.category" format
            compareClusterResult_obj@compareClusterResult$new_col <- unlist(lapply(as.character(compareClusterResult_obj@compareClusterResult$Cluster), function(x) { return(unlist(strsplit(x, "[.]"))[1]) }))
            colnames(compareClusterResult_obj@compareClusterResult)[colnames(compareClusterResult_obj@compareClusterResult) == "new_col"] <- comp_col
        }
        if (! "category" %in% colnames(compareClusterResult_obj@compareClusterResult)) {
            # create a 'category' column based on the 'Cluster' column which has the "time.category" format
            compareClusterResult_obj@compareClusterResult$category <- unlist(lapply(as.character(compareClusterResult_obj@compareClusterResult$Cluster), function(x) { return(unlist(strsplit(x, "[.]"))[2]) }))
        }
    }
    # reorder columns: make the "Cluster", comp_col and "category" columns the 3 first columns
    first_columns <- c("Cluster", comp_col, "category")
    compareClusterResult_obj@compareClusterResult <- compareClusterResult_obj@compareClusterResult[, c(first_columns, colnames(compareClusterResult_obj@compareClusterResult)[! colnames(compareClusterResult_obj@compareClusterResult) %in% first_columns])]
    
    # make 'time' and 'category' columns as factors
    compareClusterResult_obj@compareClusterResult[[comp_col]] <- as.factor(compareClusterResult_obj@compareClusterResult[[comp_col]])
    compareClusterResult_obj@compareClusterResult$category <- as.factor(compareClusterResult_obj@compareClusterResult$category)
    
    return(compareClusterResult_obj)
}


compareCluster_emapplot_terms <- function(compareClusterResult_obj, go_terms, plot_title, comp_col_level_col) {
    nb_terms <- length(go_terms)
    # emapplot
    if (length(go_terms) > 1) {
        p <- emapplot(compareClusterResult_obj, showCategory=go_terms, pie.params=list(pie="count"), cex.params=list(category_label=0.5, category_node=1, line=0.2))
    } else {
        if (is.na(go_terms)) {
            p <- emapplot(compareClusterResult_obj, pie.params=list(pie="count"), cex.params=list(category_label=0.5, category_node=1, line=0.2))
        }
    }
    p <- p + labs(title=plot_title) +
        theme(legend.position="bottom")
    if (length(comp_col_level_col) > 1) {
        p <- p + scale_fill_manual(values=comp_col_level_col)
    }
    print(p)
    # emapplot with clusters
    for (nb_clusters in c(5, 10, 15)) {
        if (nb_terms > 2 * nb_clusters) {
            if (length(go_terms) > 1) {
                p <- emapplot(compareClusterResult_obj, showCategory=go_terms, pie.params=list(pie="count"), cex.params=list(category_label=0.5, category_node=1, line=0.2), cluster.params=list(cluster=TRUE, legend=TRUE, n=nb_clusters, label_words_n=3))
            } else {
                if (is.na(go_terms)) {
                    p <- emapplot(compareClusterResult_obj, pie.params=list(pie="count"), cex.params=list(category_label=0.5, category_node=1, line=0.2), cluster.params=list(cluster=TRUE, legend=TRUE, n=nb_clusters, label_words_n=3))
                }
            }
            p <- p + labs(title=sprintf("%s, nb clusters: %d", plot_title, nb_clusters))
            if (length(comp_col_level_col) > 1) {
                p <- p + scale_fill_manual(values=comp_col_level_col)
            }
            print(p)
        }
    }
}

compareCluster_plots <- function(compareClusterResult_obj, comp_col, comp_col_level_col, out_dir, out_name) {
    # get all terms
    terms2plot <- unique(compareClusterResult_obj@compareClusterResult$Description)
    maxterms2plot <- length(terms2plot)
    ## some GO terms may have NA as similarity value with all other GO terms: remove them
    na_sim_values_terms <- rownames(compareClusterResult_obj@termsim)[apply(compareClusterResult_obj@termsim, 1, function(x) { sum(is.na(x)) } ) == dim(compareClusterResult_obj@termsim)[1]]
    compareClusterResult_obj_no_na_sim <- compareClusterResult_obj
    compareClusterResult_obj_no_na_sim@compareClusterResult <- compareClusterResult_obj_no_na_sim@compareClusterResult[which(! compareClusterResult_obj_no_na_sim@compareClusterResult$Description %in% na_sim_values_terms),]
    terms2plot_sim <- unique(compareClusterResult_obj_no_na_sim@compareClusterResult$Description)
    maxterms2plot_sim <- length(terms2plot_sim)
    
    # 30 most significant results dotplot
    pdf(sprintf("%s/%s_dotplot_top30.pdf", out_dir, out_name), height=9)
    goterms2plot <- terms2plot[1:30]
    p <- dotplot(compareClusterResult_obj, showCategory=goterms2plot, font.size=8)
    p <- p + theme(axis.text.y=element_text(size=6))
    print(p)
    dev.off()
    
    # all dotplots
    pdf(sprintf("%s/%s_dotplots.pdf", out_dir, out_name), height=9)
    for (i in seq(1, maxterms2plot, 30)) {
        if (i + 30 < maxterms2plot) {
            goterms2plot <- terms2plot[i:(i+30-1)]
        } else {
            goterms2plot <- terms2plot[i:maxterms2plot]
        }
        p <- dotplot(compareClusterResult_obj, showCategory=goterms2plot, font.size=8)
        p <- p + theme(axis.text.y=element_text(size=6))
        print(p)
    }
    dev.off()
    
    # all plots
    pdf(sprintf("%s/%s.pdf", out_dir, out_name), width=12, height=9)
    # dotplot
    for (i in seq(1, maxterms2plot, 30)) {
        if (i + 30 < maxterms2plot) {
            goterms2plot <- terms2plot[i:(i+30-1)]
        } else {
            goterms2plot <- terms2plot[i:maxterms2plot]
        }
        p <- dotplot(compareClusterResult_obj, showCategory=goterms2plot, font.size=8)
        p <- p + theme(axis.text.y=element_text(size=6))
        print(p)
    }
    # emapplot
    ## showCategory: 30 (default value)
    plot_title <- sprintf("showCategory: 30 (default value), total GO terms: %d", maxterms2plot_sim)
    compareCluster_emapplot_terms(compareClusterResult_obj_no_na_sim, NA, plot_title, comp_col_level_col)
    
    ## GO terms with p-value < 0.01
    terms2plot_sim_pval <- unique(compareClusterResult_obj_no_na_sim@compareClusterResult[which(compareClusterResult_obj_no_na_sim@compareClusterResult$p.adjust.old < 0.01), "Description"])
    maxterms2plot_sim_pval <- length(terms2plot_sim_pval)
    if (maxterms2plot_sim_pval > 30) {
        plot_title <- sprintf("showCategory: GO terms with adj pval < 0.01, nb GO terms: %d", maxterms2plot_sim_pval)
        compareCluster_emapplot_terms(compareClusterResult_obj_no_na_sim, terms2plot_sim_pval, plot_title, comp_col_level_col)
    }
    
    ## all GO terms
    if (maxterms2plot_sim_pval > 30) {
        plot_title <- sprintf("showCategory: all GO terms, nb GO terms: %d", maxterms2plot_sim)
        compareCluster_emapplot_terms(compareClusterResult_obj_no_na_sim, terms2plot_sim, plot_title, comp_col_level_col)
    }
    
    # emappplot per comparison category level and dysregulation
    if (comp_col == "cluster") {
        one_level_prefix_out_name <- "c"
    } else {
        if (comp_col == "time") {
            one_level_prefix_out_name <- "D"
        }
    }
    for (one_comp_col_level in levels(compareClusterResult_obj_no_na_sim@compareClusterResult[[comp_col]])) {
        for (dysregulation in levels(compareClusterResult_obj_no_na_sim@compareClusterResult$category)) {
            ## per timepoint and dysregulation
            print(sprintf("%s: %s, dysregulation: %s", comp_col, one_comp_col_level, dysregulation))
            ### copy compareClusterResult object
            one_comp_col_level_dys_compareClusterResult_obj_no_na_sim <- compareClusterResult_obj_no_na_sim
            one_comp_col_level_dys_compareClusterResult_obj_no_na_sim@compareClusterResult <- one_comp_col_level_dys_compareClusterResult_obj_no_na_sim@compareClusterResult[which(one_comp_col_level_dys_compareClusterResult_obj_no_na_sim@compareClusterResult[[comp_col]] == one_comp_col_level & one_comp_col_level_dys_compareClusterResult_obj_no_na_sim@compareClusterResult$category == dysregulation),]
            one_comp_col_level_dys_compareClusterResult_obj_no_na_sim@compareClusterResult <- droplevels(one_comp_col_level_dys_compareClusterResult_obj_no_na_sim@compareClusterResult)            
            
            terms2plot <- unique(one_comp_col_level_dys_compareClusterResult_obj_no_na_sim@compareClusterResult$Description)
            maxterms2plot <- length(terms2plot)
            if (maxterms2plot > 0) {
                ### duplicate last GO term with a new fake cluster (timepoint.dysregulation) to avoid error when printing emapplot with cluster.params
                #one_comp_col_level_dys_compareClusterResult_obj_no_na_sim@compareClusterResult <- rbind(one_comp_col_level_dys_compareClusterResult_obj_no_na_sim@compareClusterResult, cbind(Cluster="NA", comp="NA", category="NA", one_comp_col_level_dys_compareClusterResult_obj_no_na_sim@compareClusterResult[maxterms2plot, 4:dim(one_comp_col_level_dys_compareClusterResult_obj_no_na_sim@compareClusterResult)[2]]))
                first_columns <- c("NA", "NA", "NA")
                names(first_columns) <- c("Cluster", comp_col, "category")
                one_comp_col_level_dys_compareClusterResult_obj_no_na_sim@compareClusterResult <- rbind(one_comp_col_level_dys_compareClusterResult_obj_no_na_sim@compareClusterResult, cbind(as.data.frame(c(first_columns, one_comp_col_level_dys_compareClusterResult_obj_no_na_sim@compareClusterResult[maxterms2plot, 4:dim(one_comp_col_level_dys_compareClusterResult_obj_no_na_sim@compareClusterResult)[2]]))))
                
                ### showCategory: 30 (default value)
                plot_title <- sprintf("%s%s %s\nshowCategory: 30 (default value), total GO terms: %d", one_level_prefix_out_name, one_comp_col_level, dysregulation, maxterms2plot)
                compareCluster_emapplot_terms(one_comp_col_level_dys_compareClusterResult_obj_no_na_sim, NA, plot_title, NA)
                
                if (maxterms2plot > 30) {
                    ### showCategory: all GO terms
                    plot_title <- sprintf("%s%s %s\nshowCategory: all GO terms, total GO terms: %d", one_level_prefix_out_name, one_comp_col_level, dysregulation, maxterms2plot)
                    compareCluster_emapplot_terms(one_comp_col_level_dys_compareClusterResult_obj_no_na_sim, terms2plot, plot_title, NA)
                }
            }
        }
    }
    dev.off()
}

compareCluster_simplify2plots <- function(compareClusterResult_obj, sim_method, GOSemSimDATA_obj, sim_cutoff, gsea_bool, comp_col, plot_col, out_dir, out_name) {
    # compute similarity for all GO terms
    terms <- unique(compareClusterResult_obj@compareClusterResult$Description)
    compareClusterResult_obj_termsim <- pairwise_termsim(compareClusterResult_obj, method=sim_method, semData=GOSemSimDATA_obj, showCategory=terms)
    write.table(compareClusterResult_obj_termsim@compareClusterResult, file=sprintf("%s/%s_compareClusterResult.tsv", out_dir, out_name), sep="\t", quote=FALSE, row.names=FALSE)
    sim_out_name <- sprintf("%s_termsim_%s", out_name, sim_method)
    write.csv(compareClusterResult_obj_termsim@termsim, file=sprintf("%s/%s.csv", out_dir, sim_out_name), quote=FALSE, row.names=TRUE)
    
    # simplify
    simp_out_name <- sprintf("%s_simp_%s_%s", out_name, sim_method, sub("[.]", "_", sim_cutoff))
    compareClusterResult_obj_termsim_simp <- simplify(compareClusterResult_obj_termsim, measure=sim_method, cutoff=sim_cutoff, semData=GOSemSimDATA_obj)
    # compute similarity for all GO terms
    terms <- unique(compareClusterResult_obj_termsim_simp@compareClusterResult$Description)
    compareClusterResult_obj_termsim_simp_termsim <- pairwise_termsim(compareClusterResult_obj_termsim_simp, method=sim_method, semData=GOSemSimDATA_obj, showCategory=terms)
    
    # prepare compareClusterResult object to plot
    compareClusterResult_obj_termsim_simp_termsim <- compareCluster_prepare2plot(compareClusterResult_obj_termsim_simp_termsim, gsea_bool, comp_col)
    write.table(compareClusterResult_obj_termsim_simp_termsim@compareClusterResult, file=sprintf("%s/%s_compareClusterResult.tsv", out_dir, simp_out_name), sep="\t", quote=FALSE, row.names=FALSE)
    write.csv(compareClusterResult_obj_termsim_simp_termsim@termsim, file=sprintf("%s/%s_termsim.csv", out_dir, simp_out_name), quote=FALSE, row.names=TRUE)
    # plots
    compareCluster_plots(compareClusterResult_obj_termsim_simp_termsim, comp_col, plot_col, out_dir, simp_out_name)
}

compareCluster_launch <- function(de_df, comp_col, group_col, sim_method, GOSemSimDATA_obj, sim_cutoff, comp_dys_colors, gsea_bool, out_dir, out_name) {
    
    if (group_col == "cluster") {
        one_level_prefix_out_name <- "c"
    } else {
        if (group_col == "time") {
            one_level_prefix_out_name <- "D"
        }
    }
    
    de_only_DEGs_df <- de_df[which(de_df$p_val_adj < 0.05),]
    for (one_group_col_level in levels(de_df[[group_col]])) {
        print(sprintf("%s: %s", group_col, one_group_col_level))
        one_group_col_level_out_dir <- sprintf("%s/%s%s", out_dir, one_level_prefix_out_name, sub("[.]", "_", one_group_col_level))
        if (! dir.exists(one_group_col_level_out_dir)) {
            dir.create(one_group_col_level_out_dir, recursive=TRUE, mode="0775")
        }
        one_group_col_level_comp_out_name <- sprintf("%s_%s%s", out_name, one_level_prefix_out_name, sub("[.]", "_", one_group_col_level))
        
        ### compareCluster
        #### ORA: formula interface
        comp_cluster_GO_ORA_form_out_name <- sprintf("%s_compareCluster_GO_ORA_form", one_group_col_level_comp_out_name)
        ##### subcluster Entrez IDs
        ###### get one cluster DE results
        one_group_col_level_de_only_DEGs_df <- de_only_DEGs_df[which(de_only_DEGs_df[[group_col]]==one_group_col_level),]
        one_group_col_level_de_only_DEGs_df <- droplevels(one_group_col_level_de_only_DEGs_df)
        write.csv(one_group_col_level_de_only_DEGs_df, file=sprintf("%s/%s_Entrez.csv", one_group_col_level_out_dir, one_group_col_level_comp_out_name), quote=FALSE, row.names=FALSE)
        ###### add DE category
        DE_cat <- ifelse(one_group_col_level_de_only_DEGs_df$p_val_adj < 0.05, ifelse(one_group_col_level_de_only_DEGs_df$avg_log2FC > 0, "up", "down"), "non-DE")
        one_group_col_level_de_only_DEGs_df <- cbind(one_group_col_level_de_only_DEGs_df, category=DE_cat)
        one_group_col_level_de_only_DEGs_df$category <- as.factor(one_group_col_level_de_only_DEGs_df$category)
        ###### background genes
        universe_genes <- unique(de_df[which(de_df[[group_col]] == one_group_col_level), "ENTREZID"])
        ##### compareCluster call
        comp_cluster_GO_ORA_form <- compareCluster(as.formula(paste("ENTREZID~", comp_col, "+category")), data=one_group_col_level_de_only_DEGs_df, OrgDb="org.Rn.eg.db", ont="BP", universe=universe_genes)
        if (! is.null(comp_cluster_GO_ORA_form)) {
            #### simplify and plots
            compareCluster_simplify2plots(comp_cluster_GO_ORA_form, sim_method, GOSemSimDATA_obj, sim_cutoff, FALSE, comp_col, comp_dys_colors, one_group_col_level_out_dir, comp_cluster_GO_ORA_form_out_name)
        }

        #### GSEA
        if (gsea_bool) {
            one_group_col_level_de_df <- de_df[which(de_df[[group_col]] == one_group_col_level),]
            ##### order genes by log2FC
            comp_cluster_GSEA_FC_out_name <- sprintf("%s_compareCluster_GSEA_FC", one_group_col_level_comp_out_name)
            ###### list of Entrez IDs per time
            entrezid_list <- entrezid2GSEA(one_group_col_level_de_df, comp_col, "ENTREZID", "avg_log2FC", "p_val_adj", FALSE, one_group_col_level_out_dir, sprintf("%s_inputs", comp_cluster_GSEA_FC_out_name))
            ###### compareCluster call
            comp_cluster_GSEA_FC <- compareCluster(entrezid_list, fun="gseGO", OrgDb="org.Rn.eg.db", ont="BP")
            if (! is.null(comp_cluster_GSEA_FC)) {
                ####### simplify and plots
                compareCluster_simplify2plots(comp_cluster_GSEA_FC, sim_method, GOSemSimDATA_obj, sim_cutoff, TRUE, comp_col, comp_dys_colors, one_group_col_level_out_dir, comp_cluster_GSEA_FC_out_name)
            }

            ##### order genes by log2FC*(-log10(p-value))
            comp_cluster_GSEA_FCpval_out_name <- sprintf("%s_compareCluster_GSEA_FCpval", one_group_col_level_comp_out_name)
            ###### list of Entrez IDs per time
            entrezid_list_pval <- entrezid2GSEA(one_group_col_level_de_df, comp_col, "ENTREZID", "avg_log2FC", "p_val_adj", TRUE, one_group_col_level_out_dir, sprintf("%s_inputs", comp_cluster_GSEA_FCpval_out_name))
            ###### compareCluster call
            comp_cluster_GSEA_FCpval <- compareCluster(entrezid_list_pval, fun="gseGO", OrgDb="org.Rn.eg.db", ont="BP")
            if (! is.null(comp_cluster_GSEA_FCpval)) {
                ####### simplify and plots
                compareCluster_simplify2plots(comp_cluster_GSEA_FCpval, sim_method, GOSemSimDATA_obj, sim_cutoff,  TRUE, comp_col, comp_dys_colors, one_group_col_level_out_dir, comp_cluster_GSEA_FCpval_out_name)
            }
        }
    }
}

GO_term_comparison_analysis <- function(de_df, sim_method, GOSemSimDATA_obj, sim_cutoff, gsea_bool, out_dir, out_name) {
    # convert symbol to Entrez ID for all expressed genes
    de_symbol2entrez_df <- bitr(de_df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Rn.eg.db", drop=FALSE)
    de_entrez_df <- merge(de_df, de_symbol2entrez_df, by.x="gene", by.y="SYMBOL")
    de_entrez_df$cluster <- as.factor(de_entrez_df$cluster)
    de_entrez_df$time <- as.factor(de_entrez_df$time)

    ## cluster comparison per time
    nb_clusters <- length(levels(de_entrez_df$cluster))
    nb_clusters_updown <- 2 * nb_clusters
    if (nb_clusters_updown > 12 ) {
        cluster_dys_col <- brewer.pal(12, "Paired")
        for (i in 1:(nb_clusters-6)) {
            cluster_dys_col <- c(cluster_dys_col, brewer.pal(8, "Set2")[i], brewer.pal(8, "Dark2")[i])
        }
    } else {
        cluster_dys_col <- brewer.pal(nb_clusters_updown, "Paired")
    }
    cluster_dys_col_names <- c()
    for (one_cluster in levels(de_entrez_df$cluster)) {
        cluster_dys_col_names <- c(cluster_dys_col_names, c(sprintf("%s.down", one_cluster), sprintf("%s.up", one_cluster)))
    }
    names(cluster_dys_col) <- cluster_dys_col_names
    compareCluster_launch(de_entrez_df, "cluster", "time", sim_method, GOSemSimDATA_obj, sim_cutoff, cluster_dys_col, gsea_bool, out_dir, out_name)

    ## time comparison per cluster
    timedys_col <- brewer.pal(8, "Paired")
    timedys_col_names <- c()
    for (one_time in levels(de_entrez_df$time)) {
        timedys_col_names <- c(timedys_col_names, c(sprintf("%s.up", one_time), sprintf("%s.down", one_time)))
    }
    names(timedys_col) <- timedys_col_names
    compareCluster_launch(de_entrez_df, "time", "cluster", sim_method, GOSemSimDATA_obj, sim_cutoff, timedys_col, gsea_bool, out_dir, out_name)
}

GO_analysis_all_clusters_all_times <- function(metadata_df, nb_clusters, analysis, sim_method, GOSemSimDATA_obj, sim_cutoff, gsea_bool, comp_analysis, out_dir, out_name) {
    timepoint_comparison_out_dir <- sprintf("%s/comp", out_dir)
    timepoint_comparison_out_name <- sprintf("%s_all_times_all_clusters", out_name)
    # GO analysis
    if (analysis) {
        all_times_all_clusters_GO_ORA_df <- all_times_all_clusters_GO_ORA_simp_df <- data.frame()
        if (gsea_bool) {
            all_times_all_clusters_GSEA_df <- all_times_all_clusters_GSEA_simp_df <- data.frame()
            all_times_all_clusters_GSEA_pval_df <- all_times_all_clusters_GSEA_pval_simp_df <- data.frame()
        }
        for (one_time in levels(metadata_df$time)) {
            print(sprintf("time: %s", one_time))
            one_time_out_dir <- sprintf("%s/D%s", out_dir, one_time)
            all_clusters_GO_ORA_df <- all_clusters_GO_ORA_simp_df <- data.frame()
            all_clusters_GSEA_df <- all_clusters_GSEA_simp_df <- data.frame()
            all_clusters_GSEA_pval_df <- all_clusters_GSEA_pval_simp_df <- data.frame()
            for (one_cluster in 0:(nb_clusters-1)) {
                print(sprintf("cluster: %d", one_cluster))
                one_cluster_one_time_out_dir <- sprintf("%s/c%s", one_time_out_dir, one_cluster)
                one_cluster_one_time_out_name <- sprintf("%s_D%s_c%s", out_name, one_time, one_cluster)
                one_cluster_one_time_out_file <- sprintf("%s/%s.csv", one_cluster_one_time_out_dir, one_cluster_one_time_out_name)
                one_cluster_one_time_out_df <- read.csv(one_cluster_one_time_out_file, row.names=1)
                go_analysis_result_list <- GO_analysis(one_cluster_one_time_out_df, 0, sim_method, GOSemSimDATA_obj, gsea_bool, one_cluster_one_time_out_dir, one_cluster_one_time_out_name)

                ## get GO analysis results
                ### GO ORA
                #### whole
                ##### up
                nb_terms <- dim(go_analysis_result_list$go_ora_up)[1]
                if (nb_terms > 0) {
                    all_clusters_GO_ORA_df <- rbind(all_clusters_GO_ORA_df, cbind(cluster=rep(one_cluster, nb_terms), category=rep("up", nb_terms), data.frame(go_analysis_result_list$go_ora_up, row.names=NULL)))
                }
                ##### down
                nb_terms <- dim(go_analysis_result_list$go_ora_down)[1]
                if (nb_terms > 0) {
                    all_clusters_GO_ORA_df <- rbind(all_clusters_GO_ORA_df, cbind(cluster=rep(one_cluster, nb_terms), category=rep("down", nb_terms), data.frame(go_analysis_result_list$go_ora_down, row.names=NULL)))
                }
                #### simplify
                ##### up
                for (one_similarity_cutoff in names(go_analysis_result_list$go_ora_up_simp)) {
                    nb_terms <- dim(go_analysis_result_list$go_ora_up_simp[[one_similarity_cutoff]])[1]
                    if (nb_terms > 0) {
                        all_clusters_GO_ORA_simp_df <- rbind(all_clusters_GO_ORA_simp_df, cbind(cluster=rep(one_cluster, nb_terms), category=rep("up", nb_terms), similarity_cutoff=rep(one_similarity_cutoff, nb_terms), data.frame(go_analysis_result_list$go_ora_up_simp[[one_similarity_cutoff]], row.names=NULL)))
                    }
                }
                ##### down
                for (one_similarity_cutoff in names(go_analysis_result_list$go_ora_down_simp)) {
                    nb_terms <- dim(go_analysis_result_list$go_ora_down_simp[[one_similarity_cutoff]])[1]
                    if (nb_terms > 0) {
                        all_clusters_GO_ORA_simp_df <- rbind(all_clusters_GO_ORA_simp_df, cbind(cluster=rep(one_cluster, nb_terms), category=rep("down", nb_terms), similarity_cutoff=rep(one_similarity_cutoff, nb_terms), data.frame(go_analysis_result_list$go_ora_down_simp[[one_similarity_cutoff]], row.names=NULL)))
                    }
                }
                ### GSEA
                if (gsea_bool) {
                    #### log2FC
                    ##### whole
                    ###### up
                    nb_terms <- dim(go_analysis_result_list$gsea_up)[1]
                    if (nb_terms > 0) {
                        all_clusters_GSEA_df <- rbind(all_clusters_GSEA_df, cbind(cluster=rep(one_cluster, nb_terms), category=rep("up", nb_terms), data.frame(go_analysis_result_list$gsea_up, row.names=NULL)))
                    }
                    ###### down
                    nb_terms <- dim(go_analysis_result_list$gsea_down)[1]
                    if (nb_terms > 0) {
                        all_clusters_GSEA_df <- rbind(all_clusters_GSEA_df, cbind(cluster=rep(one_cluster, nb_terms), category=rep("down", nb_terms), data.frame(go_analysis_result_list$gsea_down, row.names=NULL)))
                    }
                    ##### simplify
                    ###### up
                    for (one_similarity_cutoff in names(go_analysis_result_list$gsea_up_simp)) {
                        nb_terms <- dim(go_analysis_result_list$gsea_up_simp[[one_similarity_cutoff]])[1]
                        if (nb_terms > 0) {
                            all_clusters_GSEA_simp_df <- rbind(all_clusters_GSEA_simp_df, cbind(cluster=rep(one_cluster, nb_terms), category=rep("up", nb_terms), similarity_cutoff=rep(one_similarity_cutoff, nb_terms), data.frame(go_analysis_result_list$gsea_up_simp[[one_similarity_cutoff]], row.names=NULL)))
                        }
                    }
                    ###### down
                    for (one_similarity_cutoff in names(go_analysis_result_list$gsea_down_simp)) {
                        nb_terms <- dim(go_analysis_result_list$gsea_down_simp[[one_similarity_cutoff]])[1]
                        if (nb_terms > 0) {
                            all_clusters_GSEA_simp_df <- rbind(all_clusters_GSEA_simp_df, cbind(cluster=rep(one_cluster, nb_terms), category=rep("down", nb_terms), similarity_cutoff=rep(one_similarity_cutoff, nb_terms), data.frame(go_analysis_result_list$gsea_down_simp[[one_similarity_cutoff]], row.names=NULL)))
                        }
                    }
                    ### GSEA
                    #### log2FClog10pval
                    ##### whole
                    ###### up
                    nb_terms <- dim(go_analysis_result_list$gsea_pval_up)[1]
                    if (nb_terms > 0) {
                        all_clusters_GSEA_pval_df <- rbind(all_clusters_GSEA_pval_df, cbind(cluster=rep(one_cluster, nb_terms), category=rep("up", nb_terms), data.frame(go_analysis_result_list$gsea_pval_up, row.names=NULL)))
                    }
                    ###### down
                    nb_terms <- dim(go_analysis_result_list$gsea_pval_down)[1]
                    if (nb_terms > 0) {
                        all_clusters_GSEA_pval_df <- rbind(all_clusters_GSEA_pval_df, cbind(cluster=rep(one_cluster, nb_terms), category=rep("down", nb_terms), data.frame(go_analysis_result_list$gsea_pval_down, row.names=NULL)))
                    }
                    ##### simplify
                    ###### up
                    for (one_similarity_cutoff in names(go_analysis_result_list$gsea_pval_up_simp)) {
                        nb_terms <- dim(go_analysis_result_list$gsea_pval_up_simp[[one_similarity_cutoff]])[1]
                        if (nb_terms > 0) {
                            all_clusters_GSEA_pval_simp_df <- rbind(all_clusters_GSEA_pval_simp_df, cbind(cluster=rep(one_cluster, nb_terms), category=rep("up", nb_terms), similarity_cutoff=rep(one_similarity_cutoff, nb_terms), data.frame(go_analysis_result_list$gsea_pval_up_simp[[one_similarity_cutoff]], row.names=NULL)))
                        }
                    }
                    ###### down
                    for (one_similarity_cutoff in names(go_analysis_result_list$gsea_pval_down_simp)) {
                        print(one_similarity_cutoff)
                        nb_terms <- dim(go_analysis_result_list$gsea_pval_down_simp[[one_similarity_cutoff]])[1]
                        if (nb_terms > 0) {
                            all_clusters_GSEA_pval_simp_df <- rbind(all_clusters_GSEA_pval_simp_df, cbind(cluster=rep(one_cluster, nb_terms), category=rep("down", nb_terms), similarity_cutoff=rep(one_similarity_cutoff, nb_terms), data.frame(go_analysis_result_list$gsea_pval_down_simp[[one_similarity_cutoff]], row.names=NULL)))
                        }
                    }
                }
            }
            all_times_all_clusters_GO_ORA_df <- rbind(all_times_all_clusters_GO_ORA_df, cbind(time=rep(one_time, dim(all_clusters_GO_ORA_df)[1]), all_clusters_GO_ORA_df))
            all_times_all_clusters_GO_ORA_simp_df <- rbind(all_times_all_clusters_GO_ORA_simp_df, cbind(time=rep(one_time, dim(all_clusters_GO_ORA_simp_df)[1]), all_clusters_GO_ORA_simp_df))
            if (gsea_bool) {
                all_times_all_clusters_GSEA_df <- rbind(all_times_all_clusters_GSEA_df, cbind(time=rep(one_time, dim(all_clusters_GSEA_df)[1]), all_clusters_GSEA_df))
                all_times_all_clusters_GSEA_simp_df <- rbind(all_times_all_clusters_GSEA_simp_df, cbind(time=rep(one_time, dim(all_clusters_GSEA_simp_df)[1]), all_clusters_GSEA_simp_df))
                all_times_all_clusters_GSEA_pval_df <- rbind(all_times_all_clusters_GSEA_pval_df, cbind(time=rep(one_time, dim(all_clusters_GSEA_pval_df)[1]), all_clusters_GSEA_pval_df))
                all_times_all_clusters_GSEA_pval_simp_df <- rbind(all_times_all_clusters_GSEA_pval_simp_df, cbind(time=rep(one_time, dim(all_clusters_GSEA_pval_simp_df)[1]), all_clusters_GSEA_pval_simp_df))
            }
        }

        ## GO term counts per time and cluster
        ### GO ORA
        #### whole GO term set
        if (dim(all_times_all_clusters_GO_ORA_df)[1] > 0) {
            GO_out_name <- sprintf("%s_GO_ORA", timepoint_comparison_out_name)
            write.table(all_times_all_clusters_GO_ORA_df, file=sprintf("%s/%s.tsv", timepoint_comparison_out_dir, GO_out_name), sep="\t", quote=FALSE, row.names=FALSE)
            category_columns <- c("time", "cluster", "category")
            go_term_count_time_cluster(all_times_all_clusters_GO_ORA_df, category_columns, "GO ORA: enriched GO terms per time and cluster", timepoint_comparison_out_dir, GO_out_name)
        }
        #### simplified
        if (dim(all_times_all_clusters_GO_ORA_simp_df)[1] > 0) {
            GO_out_name <- sprintf("%s_GO_ORA_simp", timepoint_comparison_out_name)
            write.table(all_times_all_clusters_GO_ORA_simp_df, file=sprintf("%s/%s.tsv", timepoint_comparison_out_dir, GO_out_name), sep="\t", quote=FALSE, row.names=FALSE)
            category_columns <- c("time", "cluster", "category", "similarity_cutoff")
            go_term_count_time_cluster(all_times_all_clusters_GO_ORA_simp_df, category_columns, "Simplified GO ORA: enriched GO terms per time and cluster", timepoint_comparison_out_dir, GO_out_name)
        }
        ### GSEA
        if (gsea_bool) {
            #### log2FC
            ##### whole GO term set
            if (dim(all_times_all_clusters_GSEA_df)[1] > 0) {
                GO_out_name <- sprintf("%s_GSEA", timepoint_comparison_out_name)
                write.table(all_times_all_clusters_GSEA_df, file=sprintf("%s/%s.tsv", timepoint_comparison_out_dir, GO_out_name), sep="\t", quote=FALSE, row.names=FALSE)
                category_columns <- c("time", "cluster", "category")
                go_term_count_time_cluster(all_times_all_clusters_GSEA_df, category_columns, "GSEA: enriched GO terms per time and cluster", timepoint_comparison_out_dir, GO_out_name)
            }
            #### simplified
            if (dim(all_times_all_clusters_GSEA_simp_df)[1] > 0) {
                GO_out_name <- sprintf("%s_GSEA_simp", timepoint_comparison_out_name)
                write.table(all_times_all_clusters_GSEA_simp_df, file=sprintf("%s/%s.tsv", timepoint_comparison_out_dir, GO_out_name), sep="\t", quote=FALSE, row.names=FALSE)
                category_columns <- c("time", "cluster", "category", "similarity_cutoff")
                go_term_count_time_cluster(all_times_all_clusters_GSEA_simp_df, category_columns, "Simplified GSEA: enriched GO terms per time and cluster", timepoint_comparison_out_dir, GO_out_name)
            }
            ### GSEA
            #### log2FClog10pval
            ##### whole GO term set
            if (dim(all_times_all_clusters_GSEA_pval_df)[1] > 0) {
                GO_out_name <- sprintf("%s_GSEA_pval", timepoint_comparison_out_name)
                write.table(all_times_all_clusters_GSEA_pval_df, file=sprintf("%s/%s.tsv", timepoint_comparison_out_dir, GO_out_name), sep="\t", quote=FALSE, row.names=FALSE)
                category_columns <- c("time", "cluster", "category")
                go_term_count_time_cluster(all_times_all_clusters_GSEA_pval_df, category_columns, "GSEA p-value: enriched GO terms per time and cluster", timepoint_comparison_out_dir, GO_out_name)
            }
            #### simplified
            if (dim(all_times_all_clusters_GSEA_pval_simp_df)[1] > 0) {
                GO_out_name <- sprintf("%s_GSEA_pval_simp", timepoint_comparison_out_name)
                write.table(all_times_all_clusters_GSEA_pval_simp_df, file=sprintf("%s/%s.tsv", timepoint_comparison_out_dir, GO_out_name), sep="\t", quote=FALSE, row.names=FALSE)
                category_columns <- c("time", "cluster", "category", "similarity_cutoff")
                go_term_count_time_cluster(all_times_all_clusters_GSEA_pval_simp_df, category_columns, "Simplified GSEA p-value: enriched GO terms per time and cluster", timepoint_comparison_out_dir, GO_out_name)
            }
        }
    }
    
    # GO term comparison: compareCluster() calls
    if (comp_analysis) {
        all_times_all_clusters_DE_results_df <- read.csv(sprintf("%s/%s.csv", timepoint_comparison_out_dir, timepoint_comparison_out_name))
        GO_term_comparison_analysis(all_times_all_clusters_DE_results_df, sim_method, GOSemSimDATA_obj, sim_cutoff, gsea_bool, timepoint_comparison_out_dir, out_name)
    }
}

spe2speGOpval <- function(spe_obj, go_id, go_df, go_id_num)  {
    print(sprintf("GO ID: %s", go_id))
    go_id_go_df <- go_df[which(go_df$ID == go_id),]
    go_id_go_df <- droplevels(go_id_go_df)
    go_term_pval_df <- data.frame()
    # get p-value per time and cluster
    for (one_time in levels(go_id_go_df$time)) {
        one_time_go_id_go_df <- go_id_go_df[which(go_id_go_df$time == one_time),]
        one_time_go_id_go_df <- droplevels(one_time_go_id_go_df)
        for (one_cluster in levels(one_time_go_id_go_df$cluster)) {
            one_cluster_one_time_go_id_go_df <- one_time_go_id_go_df[which(one_time_go_id_go_df$cluster == one_cluster),]
            if (dim(one_cluster_one_time_go_id_go_df)[1] == 2) {
                # significant p-value for both down and up regulation: take the most significant p-value
                min_pvalue <- min(one_cluster_one_time_go_id_go_df$p.adjust)
                one_pvalue <- one_cluster_one_time_go_id_go_df[which(one_cluster_one_time_go_id_go_df$p.adjust == min_pvalue), "p.adjust"]
                one_category <- one_cluster_one_time_go_id_go_df[which(one_cluster_one_time_go_id_go_df$p.adjust == min_pvalue), "category"]
            } else {
                if (dim(one_cluster_one_time_go_id_go_df)[1] == 1) {
                    # significant p-value only for down or up regualtion
                    one_pvalue <- one_cluster_one_time_go_id_go_df$p.adjust
                    one_category <- one_cluster_one_time_go_id_go_df$category
                }
            }
            # get barcodes corresponding to time and subcluster
            barcodes <- rownames(colData(spe_obj)[which(colData(spe_obj)$time == one_time & colData(spe_obj)$seurat_custom_clusters == one_cluster),])
            go_term_pval_df <- rbind(go_term_pval_df, data.frame(barcode=barcodes, GO_pvalue=rep(one_pvalue, length(barcodes)), category=rep(one_category, length(barcodes))))
        }
    }
    # log10
    go_term_pval_df$GO_pvalue <- -log(go_term_pval_df$GO_pvalue, 10)
    # log10(p-value) to plot: negative for downregulation, positive for upregulation
    go_term_pval_df$GO_pvalue2plot <- go_term_pval_df$GO_pvalue * ifelse(go_term_pval_df$category == "down", -1, 1)
    na_barcodes <- all_barcodes[! all_barcodes %in% go_term_pval_df$barcode]
    na_barcodes_nb <- length(na_barcodes)
    go_term_pval_df <- rbind(go_term_pval_df, data.frame(barcode=na_barcodes, GO_pvalue=rep(NA, na_barcodes_nb), category=rep(NA, na_barcodes_nb), GO_pvalue2plot=rep(NA, na_barcodes_nb)))
    ## add 'GO_pvalue' column to scolData(spe_obj): order go_term_pval_df
    rownames(go_term_pval_df) <- go_term_pval_df$barcode
    go_term_pval_df <- go_term_pval_df[rownames(colData(spe_obj)),]
    colData(spe_obj)$GO_pvalue2plot <- go_term_pval_df$GO_pvalue2plot
    
    # change sample IDs before combining multiple SpatialExperiment objects
    spe_obj@colData@rownames <- sprintf("%s_%s", go_id_num, spe_obj@colData@rownames)
    spe_obj@colData$sample_id <- sprintf("%s_%s", go_id_num, spe_obj@colData$sample_id)
    spe_obj@int_metadata$imgData$sample_id <- sprintf("%s_%s", go_id_num, spe_obj@int_metadata$imgData$sample_id)
    
    return(spe_obj)
}

go_pval_time_spatial_plots_2 <- function(spe_obj, go_df, go_ids, out_dir, out_name) {
    sample_ids <- spe_obj@int_metadata$imgData$sample_id
    sample_id_times <- c(sample_id_times, c("5", "10", "20", "40"))
    names(sample_id_times) <- sample_ids
    
    go_ids_df <- go_df[which(go_df$ID %in% go_ids),]
    go_ids2plot <- unique(go_ids_df$ID)
    pdf(sprintf("%s/%s.pdf", out_dir, out_name), width=12)
    for (one_go_id in go_ids2plot) {
        print(sprintf("GO ID: %s", one_go_id))
        one_go_id_description <- unique(go_df[which(go_df$ID == one_go_id), "Description"])
        
        one_go_id_go_df <- go_df[which(go_df$ID == one_go_id),]
        one_go_id_go_df <- droplevels(one_go_id_go_df)
        go_term_pval_df <- data.frame()
        # get p-value per time and cluster
        for (one_time in levels(one_go_id_go_df$time)) {
            one_time_one_go_id_go_df <- one_go_id_go_df[which(one_go_id_go_df$time == one_time),]
            one_time_one_go_id_go_df <- droplevels(one_time_one_go_id_go_df)
            for (one_cluster in levels(one_time_one_go_id_go_df$cluster)) {
                one_cluster_one_time_one_go_id_go_df <- one_time_one_go_id_go_df[which(one_time_one_go_id_go_df$cluster == one_cluster),]
                if (dim(one_cluster_one_time_one_go_id_go_df)[1] == 2) {
                    # significant p-value for both down and up regulation: take the most significant p-value
                    min_pvalue <- min(one_cluster_one_time_one_go_id_go_df$p.adjust)
                    one_pvalue <- one_cluster_one_time_one_go_id_go_df[which(one_cluster_one_time_one_go_id_go_df$p.adjust == min_pvalue), "p.adjust"]
                    one_category <- one_cluster_one_time_one_go_id_go_df[which(one_cluster_one_time_one_go_id_go_df$p.adjust == min_pvalue), "category"]
                } else {
                    if (dim(one_cluster_one_time_one_go_id_go_df)[1] == 1) {
                        # significant p-value only for down or up regualtion
                        one_pvalue <- one_cluster_one_time_one_go_id_go_df$p.adjust
                        one_category <- one_cluster_one_time_one_go_id_go_df$category
                    }
                }
                # get barcodes corresponding to time and subcluster
                barcodes <- rownames(colData(spe_obj)[which(colData(spe_obj)$time == one_time & colData(spe_obj)$seurat_custom_clusters == one_cluster),])
                go_term_pval_df <- rbind(go_term_pval_df, data.frame(barcode=barcodes, GO_pvalue=rep(one_pvalue, length(barcodes)), category=rep(one_category, length(barcodes))))
            }
        }
        # log10
        go_term_pval_df$GO_pvalue <- -log(go_term_pval_df$GO_pvalue, 10)
        # log10(p-value) to plot: negative for downregulation, positive for upregulation
        go_term_pval_df$GO_pvalue2plot <- go_term_pval_df$GO_pvalue * ifelse(go_term_pval_df$category == "down", -1, 1)
        na_barcodes <- all_barcodes[! all_barcodes %in% go_term_pval_df$barcode]
        na_barcodes_nb <- length(na_barcodes)
        go_term_pval_df <- rbind(go_term_pval_df, data.frame(barcode=na_barcodes, GO_pvalue=rep(NA, na_barcodes_nb), category=rep(NA, na_barcodes_nb), GO_pvalue2plot=rep(NA, na_barcodes_nb)))
        ## add 'GO_pvalue' column to scolData(spe_obj): order go_term_pval_df
        rownames(go_term_pval_df) <- go_term_pval_df$barcode
        go_term_pval_df <- go_term_pval_df[rownames(colData(spe_obj)),]
        colData(spe_obj)$GO_pvalue2plot <- go_term_pval_df$GO_pvalue2plot
        
        min_pval <- min(colData(spe_obj)$GO_pvalue2plot, na.rm=TRUE)
        if (min_pval > 0) {
            min_color <- "grey50"
            min_pval <- 0
            scale_close_to_zero_inf <- 0
        } else {
            min_color <- "blue"
            scale_close_to_zero_inf <- -1
        }
        max_pval <- max(colData(spe_obj)$GO_pvalue2plot, na.rm=TRUE)
        if (max_pval < 0) {
            max_color <- "grey50"
            max_pval <- 0
            scale_close_to_zero_sup <- 0
        } else {
            max_color <- "red"
            scale_close_to_zero_sup <- 1
        }
        plot_title <- sprintf("%s: %s", one_go_id, one_go_id_description)
        p <- plotVisium(spe_obj, fill="GO_pvalue2plot", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
        p <- p + scale_fill_gradientn(colours=c(min_color, "grey50", max_color), na.value="white", values=scales::rescale(c(min_pval, scale_close_to_zero_inf, 0, scale_close_to_zero_sup, max_pval))) +
            labs(title=plot_title) +
            theme(legend.position="bottom", legend.key.size=unit(1, "lines"))
            #theme(legend.position="bottom", legend.key.size=unit(1, "lines"), legend.title=element_text(size=8), legend.text=element_text(size=7))
        ## modify facet number of columns
        p$facet$params$ncol <- 4
        ## modify facet titles
        p$facet$params$labeller <- labeller(sample_id=sample_id_times)
        ## modify point size
        p$layers[[5]]$aes_params$size <- 1
        ## modify legend title
        p$guides$fill$title <- "-1^{0,1}*log10(p-value)"
        print(p)
    }
    dev.off()
    
    return(go_ids_df)
}

