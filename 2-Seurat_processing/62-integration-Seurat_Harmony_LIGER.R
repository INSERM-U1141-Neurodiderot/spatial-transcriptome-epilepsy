.libPaths(c(.libPaths(), "/home/christophe.lepriol/NeuroDev_ADD/R/r_4.1.0"))
library(ggplot2)
library(patchwork)
library(dplyr)
#library(cowplot)
#library(kBET)


##############
# Parameters #
##############
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"
src_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/src", work_dir)
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
timepoint <- 10
if (! is.na(timepoint)) {
    sections <- epireg_visium_metadata_df[which(epireg_visium_metadata_df$time==timepoint), "sample"]
}
nb_sections <- length(sections)

# output directory name
if (! is.na(timepoint)) {
    output_dirname <- sprintf("D%d_samples", timepoint)
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
## QC
spot_sum_of_counts_threshold <- 1000
spot_detected_genes_threshold <- 500
spot_mito_gene_count_pct_threshold <- 100
spot_hb_gene_count_pct_threshold <- 10
## normalization, feature selection, dimension reduction, clustering
normalization_method_vector <- c("SCTransform", "LogNormalize")
#normalization_method_vector <- c("LogNormalize")
#feature_residual_variance_threshold <- 1.1
#feature_residual_variance_threshold <- NA
#feature_residual_variance_threshold_vector <- c(1, 1.1, 1.3, 1.5, 2, 2.5)
#features_param_vector <- c(1.1, 1.3, 1.5, 2)
#features_param_vector <- c(500, 1000, 1500, 2000, 2500)
#features_param_vector <- c(250, 500, 1000, 1500, 2000)
#features_param_vector <- c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000)
features_param_vector <- c(500, 1000, 1500, 2000, 2500, 3000)
#features_param_vector <- c(500)
#dimensions_nb_param_vector <- c(10, 15, 20)
dimensions_nb_param_vector <- c(15, 20, 30)
#dimensions_nb_param_vector <- c(10)
#clustering_resolution_vector <- c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2)
#clustering_resolution_vector <- c(0.4, 0.6, 0.8, 1)
#clustering_resolution_vector <- c(0.4, 1)
#clustering_resolution_vector <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)
#clustering_resolution_vector <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)
clustering_resolution_vector <- c(0.1, 0.3, 0.5, 0.8, 1, 1.2)
#clustering_resolution_vector <- c(0.4)
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
source(sprintf("%s/20-pipelines-Seurat-functions.R", src_dir))

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

add_qc_stats <- function(sp_exp, mito_genes, hb_genes, out_dir, out_name) {
    # identify mitochondrial genes
    #is_mito <- grepl("(^MT(-|.))|(^mt(-|.))", rowData(sp_exp)$gene_name)
    is_mito <- rowData(sp_exp)$gene_name %in% mito_genes
    # identify hemoglobin genes
    #is_hb <- grepl("^Hb(a|b|q|s|z)", rowData(sp_exp)$gene_name)
    is_hb <- rowData(sp_exp)$gene_name %in% hb_genes
    # calculate per-spot QC metrics and store in colData
    sp_exp <- addPerCellQC(sp_exp, subsets=list(mito=is_mito, hb=is_hb))
    write.csv(colData(sp_exp), file=sprintf("%s/%s_%s.csv", out_dir, sp_exp@int_metadata$imgData$sample_id, out_name), quote=FALSE, row.names=TRUE)
    return(sp_exp)
}

qc_stat_plots <- function(sp_exp, out_dir, out_name) {
    pdf(sprintf("%s/%s_%s.pdf", out_dir, sp_exp@int_metadata$imgData$sample_id, out_name))
    # sum of counts
    colname <- "sum"
    print(plotQC(sp_exp, type="scatter", metric_x="sum", metric_y="detected"))
    ## distribution
    spot_data <- colData(sp_exp)[[colname]]
    percentiles <- density_plot(spot_data, "Detected genes per spot", "Sum of counts")    
    ## spatial plots
    print(plotSpots(sp_exp, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate=colname, palette=c("white", "black"), size=2))
    print(plotVisium(sp_exp, fill=colname, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    
    # number of expressed features
    colname <- "detected"
    ## distribution
    spot_data <- colData(sp_exp)[[colname]]
    percentiles <- density_plot(spot_data, "Sum of counts per spot", "Detected genes")    
    ## spatial plots
    print(plotSpots(sp_exp, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate=colname, palette=c("white", "black"), size=2))
    print(plotVisium(sp_exp, fill=colname, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    
    # percentage of mitochondrial gene counts
    colname <- "subsets_mito_percent"
    print(plotQC(sp_exp, type="scatter", metric_x="sum", metric_y=colname))
    ## distribution
    spot_data <- colData(sp_exp)[[colname]]
    percentiles <- density_plot(spot_data, "Percentage of mitochondrial gene counts per spot", "Percentage of mitochondrial gene counts")    
    ## spatial plots
    print(plotSpots(sp_exp, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate=colname, palette=c("white", "black"), size=2))
    print(plotVisium(sp_exp, fill=colname, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    
    # percentage of hemoglobin gene counts
    colname <- "subsets_hb_percent"
    print(plotQC(sp_exp, type="scatter", metric_x="sum", metric_y=colname))
    ## distribution
    spot_data <- colData(sp_exp)[[colname]]
    percentiles <- density_plot(spot_data, "Percentage of hemoglobin gene counts per spot", "Percentage of hemoglobin gene counts")    
    ## spatial plots
    print(plotSpots(sp_exp, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate=colname, palette=c("white", "black"), size=2))
    print(plotVisium(sp_exp, fill=colname, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    dev.off()
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

gene_spatial_mean_exp_plots <- function(seurat_obj, gene, assay2plot) {
    # spatial feature plot
    p <- SpatialFeaturePlot(seurat_obj, features=gene, alpha=c(0.5,1)) + 
        labs(title=sprintf("%s\nmean expression: %.2f, expressed in %.2f%% spots\nmean expression rank: %d", gene, seurat_obj@assays[[assay2plot]]@meta.features[gene, "mean"], seurat_obj@assays[[assay2plot]]@meta.features[gene, "detected"], floor(seurat_obj@assays[[assay2plot]]@meta.features[gene, "mean_rank"]))) +
        theme(legend.position="bottom")
    print(p)

    # mean expression density plot
    gene_log_exp <- log(seurat_obj@assays[[assay2plot]]@counts[gene,], 2)
    gene_percentiles <- plot_density(gene_log_exp, sprintf("%s expression", gene), "log2(Expression)")

    return(gene_percentiles)
}

qc_spatial_mean_exp_plots <- function(seurat_obj, gene_vec, mean_order, out_dir, out_name) {
    sample <- seurat_obj@project.name
    active_assay <- seurat_obj@active.assay
    gene_vec <- gsub("-", ".", gene_vec) # '-' were replaced by '.' during Seurat object creation
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
        one_gene_mean_exp_percentiles <- gene_spatial_mean_exp_plots(seurat_obj, one_gene, active_assay)
        genes2plot_df <- rbind(genes2plot_df, c(one_gene, one_gene_mean_exp_percentiles))
    }
    dev.off()
    colnames(genes2plot_df) <- c("gene", names(one_gene_mean_exp_percentiles))
    write.csv(genes2plot_df, file=sprintf("%s/%s_%s.csv", out_dir, sample, out_name), quote=FALSE, row.names=FALSE)
    
    return(cbind(sample=rep(sample, dim(genes2plot_df)[1]), genes2plot_df))
}

qc_filter <- function(sp_exp, sum_thres, dectected_thres, mito_pct_thres, hb_pct_thres, out_dir, out_name) {
    pdf(sprintf("%s/%s_%s.pdf", out_dir, sp_exp@int_metadata$imgData$sample_id, out_name))
    # identify spots to filter according different statistics
    ## sum of counts
    colname <- "sum"
    threshold <- sum_thres
    sum_qc_colname <- sprintf("qc_%s_%d", colname, threshold)
    sp_exp <- qc_threshold(sp_exp, colname, threshold, FALSE, sum_qc_colname)
    qc_threshold_plot(sp_exp, sum_qc_colname)
    ## number of expressed features
    colname <- "detected"
    threshold <- dectected_thres
    detected_qc_colname <- sprintf("qc_%s_%d", colname, threshold)
    sp_exp <- qc_threshold(sp_exp, colname, threshold, FALSE, detected_qc_colname)
    qc_threshold_plot(sp_exp, detected_qc_colname)
    ## percentage of mitochondrial gene counts
    colname <- "subsets_mito_percent"
    threshold <- mito_pct_thres
    mito_qc_colname <- sprintf("qc_%s_%d", colname, threshold)
    sp_exp <- qc_threshold(sp_exp, colname, threshold, TRUE, mito_qc_colname)
    qc_threshold_plot(sp_exp, mito_qc_colname)
    ## percentage of hemoglobin gene counts
    colname <- "subsets_hb_percent"
    threshold <- hb_pct_thres
    hb_qc_colname <- sprintf("qc_%s_%d", colname, threshold)
    sp_exp <- qc_threshold(sp_exp, colname, threshold, TRUE, hb_qc_colname)
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

addFeatureStats2Seurat <- function(seurat_obj, assay2use) {
    mean_exp <- apply(seurat_obj@assays[[assay2use]]@counts, 1, mean)
    nb_spots <- dim(seurat_obj@assays[[assay2use]]@counts)[2]
    detected_spots <- apply(seurat_obj@assays[[assay2use]]@counts, 1, function(x, spots=nb_spots) { length(x[x > 0])/ spots * 100 })
    gene_stats_df <- data.frame(gene_name=rownames(seurat_obj@assays[[assay2use]]@counts), mean=mean_exp, detected=detected_spots)
    gene_stats_df$mean_rank <- rank(-gene_stats_df$mean)
    seurat_obj@assays[[assay2use]]@meta.features <- gene_stats_df
    return(seurat_obj)
}

SpatialExperiment_Seurat_conversion <- function(sp_exp, sr_dir, add_image) {
    # Space Ranger output directory for tissue positions list and image
    sample <- sp_exp@int_metadata$imgData$sample_id
    space_ranger_out_dir <- sprintf("%s/%s/10-Pipeline/outs", sr_dir, sample)
    
    # meta data
    metadata_df <- merge(as.data.frame(colData(sp_exp))[, c("barcode", "in_tissue", "array_row", "array_col", "lame", "zone", "condition")], as.data.frame(sp_exp@int_colData$spatialCoords), by.x="barcode", by.y="row.names")
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

seurat_list_normalization <- function(seurat_obj_list, method, feat_param, only_var, out_dir) {
    if (method == "SCTransform") {
        if (feat_param > 10) {
            ### feat_param is a number of features
            print(sprintf("features param: %d", feat_param))
            normalization_out_dir <- sprintf("%s/00-SCTransform/featurenumber%d", out_dir, feat_param)
            if (! dir.exists(normalization_out_dir)) {
                dir.create(normalization_out_dir, recursive=TRUE, mode="0775")
            }
            normalization_out_name <- sprintf("%s_featurenumber%d", method, feat_param)
            seurat_objects_norm <- lapply(seurat_obj_list, FUN=SCTransform, variable.features.n=feat_param, return.only.var.genes=only_var)
            integration_features_nb <- feat_param
        } else {
            ### feat_param is a residual variance threshold
            print(sprintf("features param: %.1f", feat_param))
            normalization_out_dir <- sprintf("%s/00-SCTransform/residualvariance%s", out_dir, sub("[.]", "_", feat_param))
            if (! dir.exists(normalization_out_dir)) {
                dir.create(normalization_out_dir, recursive=TRUE, mode="0775")
            }
            normalization_out_name <- sprintf("%s_residualvariance%s", method, sub("[.]", "_", feat_param))
            seurat_objects_norm <- lapply(seurat_obj_list, FUN=SCTransform, variable.features.n=NULL, variable.features.rv.th=feat_param, return.only.var.genes=only_var)
            ### get minimum number of variable features
            min_variable_features <- length(seurat_objects_norm[[sections[1]]]@assays$SCT@var.features)
            for (sample in sections[2:length(sections)]) {
                if (length(seurat_objects_norm[[sample]]@assays$SCT@var.features) < min_variable_features) {
                    min_variable_features <- length(seurat_objects_norm[[sample]]@assays$SCT@var.features)
                }
            }
            integration_features_nb <- min_variable_features
        }
        integration_normalization_method <- default_assay <- "SCT"
        selection_method <- "sct"

        ## add feature statistics to Seurat object after SCTransform
        seurat_objects_norm <- lapply(seurat_objects_norm, addFeatureStats2Seurat, assay2use=default_assay)
    } else {
        if (method == "LogNormalize") {
            selection_method <- "vst"
            ### feat_param is a number of features
            print(sprintf("features param: %d", feat_param))
            normalization_out_dir <- sprintf("%s/10-%s/featurenumber%d", out_dir, method, feat_param)
            if (! dir.exists(normalization_out_dir)) {
                dir.create(normalization_out_dir, recursive=TRUE, mode="0775")
            }
            normalization_out_name <- sprintf("%s_featurenumber%d", method, feat_param)

            seurat_objects_norm <- lapply(seurat_objects, function(x, nb_features=feat_param) {
                x <- NormalizeData(x)
                x <- FindVariableFeatures(x, selection.method=selection_method, nfeatures=nb_features)
            })
            integration_features_nb <- feat_param
            integration_normalization_method <- method
            default_assay <- "RNA"
        }
    }
    
    return(list(seurat_list=seurat_objects_norm, features_nb=integration_features_nb, norm_method=integration_normalization_method, assay=default_assay, selec=selection_method, out_dir=normalization_out_dir, out_name=normalization_out_name))
}

hvg_plots <- function(seurat_obj, assay2use, selec, nb_hvgs, out_dir, out_name) {
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
    exp_df <- qc_spatial_mean_exp_plots(seurat_obj, hvgs2plot, FALSE, out_dir, out_name)
    
    return(exp_df)
}

get_HVGs_merged_Seurat_SCT <- function(seurat_obj_merge, n_hvgs) {
    hvgs <- c()
    for (one_model in names(seurat_obj_merge@assays$SCT@SCTModel.list)) {
        seurat_obj_merge_one_model_feature_attributes <- seurat_obj_merge@assays$SCT@SCTModel.list[[one_model]]@feature.attributes
        hvgs <- c(hvgs, rownames(seurat_obj_merge_one_model_feature_attributes[order(seurat_obj_merge_one_model_feature_attributes$residual_variance, decreasing=TRUE),])[1:n_hvgs])
    }
    hvgs <- unique(hvgs)
    # some genes are not present in all SCTModel: only keep those present in all models
    for (one_model in names(seurat_obj_merge@assays$SCT@SCTModel.list)) {
        hvgs <- hvgs[hvgs %in% rownames(seurat_obj_merge@assays$SCT@SCTModel.list[[one_model]]@feature.attributes)]
    }
    return(hvgs)
}

get_HVGs_merged_Seurat_LogNormalize <- function(seurat_obj_merge, seurat_obj_norm_list, n_hvgs) {
    hvgs <- c()
    for (one_sample in names(seurat_obj_merge@images)) {
        hvgs <- c(hvgs, VariableFeatures(seurat_obj_norm_list[[one_sample]]))
    }
    hvgs <- unique(hvgs)
    #### some genes are not expressed in all samples: only keep those expressed in all samples
    for (one_sample in names(seurat_obj_merge@images)) {
        hvgs <- hvgs[hvgs %in% rownames(seurat_obj_norm_list[[one_sample]]@assays$RNA@counts)]
    }
    return(hvgs)
}

set_SVGs_variable_features <- function(seurat_obj, assay2use) {
    seurat_obj_svgs <- FindSpatiallyVariableFeatures(seurat_obj, assay=assay2use, selection.method="markvariogram")
    # considering that seurat_obj_svgs is a character vector storing SVG names
    VariableFeatures(seurat_obj) <- seurat_obj_svgs
    return(seurat_obj)
}

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

marker_analysis <- function(seurat_obj_int, assay2use, clusters, load_image, nb_plots, plots_per_page, DEGs, conserved, cluster_col, condition_col, out_dir, out_name) {
    # cluster markers
    if (clusters) {
        ## find markers for every cluster compared to all remaining spots, report only the positive ones
        markers_out_dir <- sprintf("%s/00-Cluster_markers", out_dir)
        if (! dir.exists(markers_out_dir)) {
            dir.create(markers_out_dir, recursive=TRUE, mode="0775")
        }
        identify_cluster_markers(seurat_obj_int, assay2use, load_image, nb_plots, plots_per_page, markers_out_dir, sprintf("%s_cluster_marker_genes", out_name))
    }

    # differentially expressed genes across conditions
    if (DEGs) {
        DEG_out_dir <- sprintf("%s/10-DEGs", out_dir)
        if (! dir.exists(DEG_out_dir)) {
            dir.create(DEG_out_dir, recursive=TRUE, mode="0775")
        }
        identify_DEGs(seurat_obj_int, assay2use, cluster_col, condition_col, DEG_out_dir, sprintf("%s_DEGs", out_name))
    }

    # conserved cell type markers
    if (conserved) {
        conserved_markers_out_dir <- sprintf("%s/20-Conserved_markers", out_dir)
        if (! dir.exists(conserved_markers_out_dir)) {
            dir.create(conserved_markers_out_dir, recursive=TRUE, mode="0775")
        }
        identify_conserved_markers(seurat_obj_int, assay2use, cluster_col, condition_col, conserved_markers_out_dir, sprintf("%s_conserved_markers", out_name))
    }
}


############
# Analysis #
############

# load data
print("load data")
library(SpatialExperiment)
library(ggspavis)
st_pipeline_dir <- sprintf("%s/10-ST_analysis/00-ST_Pipeline/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples", work_dir, genome_name)
space_ranger_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples", work_dir, genome_name)
#seurat_objects <- lapply(sections, load_data, st_dir=st_pipeline_dir, sr_dir=space_ranger_dir, add_image=load_image)
#names(seurat_objects) <- sections
spe_objects <- lapply(sections, load_ST_pipeline_data, st_dir=st_pipeline_dir, sr_dir=space_ranger_dir)
names(spe_objects) <- sections

# quality control
library(scater) # Loading required package: scuttle
qc_output_dir <- sprintf("%s/Norm_merge-QC-sum%d_detected%s_mito%d_hb%d", output_dir, spot_sum_of_counts_threshold, spot_detected_genes_threshold, spot_mito_gene_count_pct_threshold, spot_hb_gene_count_pct_threshold)
if (! dir.exists(qc_output_dir)) {
    dir.create(qc_output_dir, recursive=TRUE, mode="0775")
}
## feature statistics
spe_objects <- lapply(spe_objects, addFeatureStats2SpatialExperiment)

## identify mitochondrial and hemoglobin genes
library(rtracklayer)
raw_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/data/ref_genome/Rattus_norvegicus"
gtf_file <- sprintf("%s/mRatBN7.2/annotations/NCBI_RefSeq/release_108/GCF_015227675.2_mRatBN7.2_genomic.gtf.gz", raw_dir)
gtf_GRanges <- import(gtf_file, format="gtf")
gtf_GRanges_gene <- gtf_GRanges[which(gtf_GRanges$type=="gene"),]
## mitochondrial genes
mt_genes <- gtf_gene_regex(gtf_GRanges_gene, "gene_id", "mt-", qc_output_dir)
#mitochondrial_genes <- gtf_gene_regex(gtf_GRanges_gene, "description", "mitochondrial", qc_output_dir)
#mitochondrial_geneset <- c(mt_genes, mitochondrial_genes)
mitochondrial_geneset <- mt_genes
mitochondrial_geneset_spe <- gsub("-", ".", mitochondrial_geneset) # '-' were replaced by '.' during SpatialExperiment object creation
## hemoglobin genes
#hemoglobin_genes <- gtf_gene_regex(gtf_GRanges_gene, "description", "(hemo|neuro|myo|cyto|hapto| )globin", qc_output_dir)
hemoglobin_genes <- gtf_gene_regex(gtf_GRanges_gene, "description", "hemoglobin", qc_output_dir)
### add LOC120093065: https://www.ncbi.nlm.nih.gov/gene/120093065
hemoglobin_geneset <- c(hemoglobin_genes, "LOC120093065")
hemoglobin_geneset_spe <- gsub("-", ".", hemoglobin_geneset) # '-' were replaced by '.' during SpatialExperiment object creation
## globin genes
#globin_genes <- gtf_gene_regex(gtf_GRanges_gene, "description", "globin", qc_output_dir)
rm(gtf_GRanges)
rm(gtf_GRanges_gene)
gc()
detach("package:rtracklayer", unload=TRUE)

## add QC statistics
qc_stats_output_dir <- sprintf("%s/00-QC_stats", qc_output_dir)
if (! dir.exists(qc_stats_output_dir)) {
    dir.create(qc_stats_output_dir, recursive=TRUE, mode="0775")
}
spe_objects <- lapply(spe_objects, add_qc_stats, mito_genes=mitochondrial_geneset_spe, hb_genes=hemoglobin_geneset_spe, out_dir=qc_stats_output_dir, out_name="QC_stats")
lapply(spe_objects, qc_stat_plots, out_dir=qc_stats_output_dir, out_name="QC_stats")

## mitochondrial and hemoglobin gene spatial plots
### SpatialExperiment to Seurat object conversion
#### detach SpatialExperiment package to avoid conflict with Seurat package
detach("package:ggspavis", unload=TRUE) # also unload ggspavis package which imports ‘SpatialExperiment’ namespace
detach("package:SpatialExperiment", unload=TRUE)
library(Seurat) # also attaching SeuratObject package
seurat_objects <- lapply(spe_objects, SpatialExperiment_Seurat_conversion, sr_dir=space_ranger_dir, add_image=load_image) # Loading required package: SpatialExperiment
#### SpatialExperiment package is required for Seurat object creation, if it is not loaded before creating Seurat object, then it is automatically loaded during Seurat object creation
#### detach SpatialExperiment package to avoid conflict with Seurat package
detach("package:SpatialExperiment", unload=TRUE)
### hemoglobin genes
df_list <- lapply(seurat_objects, qc_spatial_mean_exp_plots, gene_vec=hemoglobin_geneset_spe, mean_order=TRUE, out_dir=qc_stats_output_dir, out_name="hemoglobin_genes")
### mitochondrial genes
df_list <- lapply(seurat_objects, qc_spatial_mean_exp_plots, gene_vec=mitochondrial_geneset_spe, mean_order=TRUE, out_dir=qc_stats_output_dir, out_name="mitochondrial_genes")
### detach Seurat and SeuratObject packages to avoid conflict with SpatialExperiment package
detach("package:Seurat", unload=TRUE)
detach("package:SeuratObject", unload=TRUE)

## QC filtering
library(SpatialExperiment)
library(ggspavis)
qc_filtering_output_dir <- sprintf("%s/10-QC_filtering", qc_output_dir)
if (! dir.exists(qc_filtering_output_dir)) {
    dir.create(qc_filtering_output_dir, recursive=TRUE, mode="0775")
}
spe_objects <- lapply(spe_objects, qc_filter, sum_thres=spot_sum_of_counts_threshold, dectected_thres=spot_detected_genes_threshold, mito_pct_thres=spot_mito_gene_count_pct_threshold, hb_pct_thres=spot_hb_gene_count_pct_threshold, out_dir=qc_filtering_output_dir, out_name="filtered_spots")
### QC statistics after filerting
lapply(spe_objects, function(x, out_dir=qc_filtering_output_dir, out_name="filtered_QC_stats") {
    write.csv(colData(x), file=sprintf("%s/%s_%s.csv", out_dir, x@int_metadata$imgData$sample_id, out_name), quote=FALSE, row.names=TRUE)
})
lapply(spe_objects, qc_stat_plots, out_dir=qc_filtering_output_dir, out_name="filtered_QC_stats")
### feature statistics after filtering
spe_objects <- lapply(spe_objects, addFeatureStats2SpatialExperiment)
### mitochondrial and hemoglobin gene spatial plots after filerting
#### SpatialExperiment to Seurat object conversion
##### detach SpatialExperiment package to avoid conflict with Seurat package
detach("package:ggspavis", unload=TRUE) # also unload ggspavis package which imports ‘SpatialExperiment’ namespace
detach("package:SpatialExperiment", unload=TRUE)
library(Seurat) # also attaching SeuratObject package
seurat_objects <- lapply(spe_objects, SpatialExperiment_Seurat_conversion, sr_dir=space_ranger_dir, add_image=load_image) # does not load SpatialExperiment package !
#detach("package:SpatialExperiment", unload=TRUE)
rm(spe_objects)
gc()
#### hemoglobin genes
df_list <- lapply(seurat_objects, qc_spatial_mean_exp_plots, gene_vec=hemoglobin_geneset_spe, mean_order=TRUE, out_dir=qc_filtering_output_dir, out_name="filtered_hemoglobin_genes")
#### mitochondrial genes
df_list <- lapply(seurat_objects, qc_spatial_mean_exp_plots, gene_vec=mitochondrial_geneset_spe, mean_order=TRUE, out_dir=qc_filtering_output_dir, out_name="filtered_mitochondrial_genes")
detach("package:scater", unload=TRUE)
unloadNamespace("DropletUtils") # also unload DropletUtils package which imports ‘scuttle’ namespace
#detach("package:DropletUtils", unload=TRUE) # also unload DropletUtils package which imports ‘scuttle’ namespace
detach("package:scuttle", unload=TRUE)


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

integration_batch_mean_silhouette_df <- data.frame()
integration_cluster_mean_silhouette_df <- data.frame()

library(cluster) # silhouette() function

# non-integrated datasets
integration_method <- "no integration"
integration_output_dir <- sprintf("%s/00-No_integration", qc_filtering_output_dir)
if (! dir.exists(integration_output_dir)) {
    dir.create(integration_output_dir, recursive=TRUE, mode="0775")
}
for (normalization_method in normalization_method_vector) {
    print(sprintf("normalization method: %s", normalization_method))
    for (features_param in features_param_vector) {
        ## normalization
        return_only_var_genes <- ifelse(normalization_method == "SCTransform", FALSE, NA)
        returned_list <- seurat_list_normalization(seurat_objects, normalization_method, features_param, return_only_var_genes, integration_output_dir)
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
        hvg_df_list <- lapply(seurat_objects_norm, hvg_plots, assay2use=default_assay, selec=selection_method, nb_hvgs=features_param, out_dir=normalization_output_dir, out_name= "HVGs")
        
        ## merge datasets
        seurat_objects_merge <- merge(seurat_objects_norm[[1]], y=unlist(seurat_objects_norm)[2:nb_sections], add.cell.ids=sections)
        ### orig.ident, lame, zone and condition are not factors in the merged object meta data
        seurat_objects_merge@meta.data$orig.ident <- as.factor(seurat_objects_merge@meta.data$orig.ident)
        seurat_objects_merge@meta.data$lame <- as.factor(seurat_objects_merge@meta.data$lame)
        seurat_objects_merge@meta.data$zone <- as.factor(seurat_objects_merge@meta.data$zone)
        seurat_objects_merge@meta.data$condition <- as.factor(seurat_objects_merge@meta.data$condition)
        ### set variable features: after merging datasets, variable features are not set in the merged Seurat object, which makes RunPCA() fails: Error in PrepDR(object = object, features = features, verbose = verbose): Variable features haven't been set. Run FindVariableFeatures() or provide a vector of feature names.
        #### get the union of the HVGs of each dataset
        if (normalization_method == "SCTransform") {
            seurat_objects_merge_hvgs <- get_HVGs_merged_Seurat_SCT(seurat_objects_merge, features_param)
        } else {
            seurat_objects_merge_hvgs <- get_HVGs_merged_Seurat_LogNormalize(seurat_objects_merge, seurat_objects_norm, features_param)
        }
        VariableFeatures(seurat_objects_merge) <- seurat_objects_merge_hvgs
        
        rm(seurat_objects_norm)
        gc()
        
        if (normalization_method == "LogNormalize") {
            seurat_objects_merge <- ScaleData(seurat_objects_merge)
        }
        
        ## dimensionality reduction
        seurat_objects_merge <- RunPCA(seurat_objects_merge, verbose=FALSE)
        reduction2use <- "pca"
        for (dimensions_nb_param in dimensions_nb_param_vector) {
            print(sprintf("dimension number: %d", dimensions_nb_param))
            dimensions_output_dir <- sprintf("%s/%d_PCs", normalization_output_dir, dimensions_nb_param)
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
            mean_silhouette_df <- cbind(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_param, mean_silhouette_df)
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
                mean_silhouette_df <- data.frame(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_param, resolution=clustering_resolution, mean_silhouette_seurat_clusters=cluster_mean_silhouette, mean_silhouette_seurat_clusters_hippocampus=hippocampus_cluster_mean_silhouette)
                integration_cluster_mean_silhouette_df <- rbind(integration_cluster_mean_silhouette_df, mean_silhouette_df)
                
                rm(seurat_objects_merge_hippocampus)
                gc()
            }
        }
        rm(seurat_objects_merge)
        gc()
    }
}

# integrated datasets
## Seurat
### need to set maxSize to avoid error when identifying conserved markers: Error in getGlobalsAndPackages(expr, envir = envir, globals = globals): The total size of the 3 globals exported for future expression (‘FUN()’) is 1.38 GiB.. This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize').
#options(future.globals.maxSize=2000 * 1024^2) # set allowed size to 2K MiB
integration_method <- "Seurat"
integration_output_dir <- sprintf("%s/10-Seurat", qc_filtering_output_dir)
if (! dir.exists(integration_output_dir)) {
    dir.create(integration_output_dir, recursive=TRUE, mode="0775")
}
for (normalization_method in normalization_method_vector) {
    print(sprintf("normalization method: %s", normalization_method))
    for (features_param in features_param_vector) {
        ## normalization
        return_only_var_genes <- ifelse(normalization_method == "SCTransform", FALSE, NA)
        returned_list <- seurat_list_normalization(seurat_objects, normalization_method, features_param, return_only_var_genes, integration_output_dir)
        seurat_objects_norm <- returned_list[["seurat_list"]]
        integration_features_nb <- returned_list[["features_nb"]]
        integration_normalization_method <- returned_list[["norm_method"]]
        default_assay <- returned_list[["assay"]]
        selection_method <- returned_list[["selec"]]
        normalization_output_dir <- returned_list[["out_dir"]]
        normalization_output_name <- returned_list[["out_name"]]
        rm(returned_list)
        gc()
        
        # HVG plots
        hvg_df_list <- lapply(seurat_objects_norm, hvg_plots, assay2use=default_assay, selec=selection_method, nb_hvgs=features_param, out_dir=normalization_output_dir, out_name= "HVGs")
        
        ## integration
        integration_output_name <- sprintf("%s_%s", normalization_output_name, integration_method)
        features <- SelectIntegrationFeatures(seurat_objects_norm, nfeatures=integration_features_nb)
        if (normalization_method == "SCTransform") {
            seurat_objects_norm <- PrepSCTIntegration(object.list=seurat_objects_norm, anchor.features=features)
        }
        int.anchors <- FindIntegrationAnchors(object.list=seurat_objects_norm, normalization.method=integration_normalization_method, anchor.features=features)
        rm(seurat_objects_norm)
        gc()
        seurat_objects_integrated <- IntegrateData(anchorset=int.anchors, normalization.method=integration_normalization_method)
        rm(int.anchors)
        gc()
        
        ### orig.ident, lame, zone and condition are not factors in the integrated object meta data
        seurat_objects_integrated@meta.data$orig.ident <- as.factor(seurat_objects_integrated@meta.data$orig.ident)
        seurat_objects_integrated@meta.data$lame <- as.factor(seurat_objects_integrated@meta.data$lame)
        seurat_objects_integrated@meta.data$zone <- as.factor(seurat_objects_integrated@meta.data$zone)
        seurat_objects_integrated@meta.data$condition <- as.factor(seurat_objects_integrated@meta.data$condition)
        
        if (normalization_method == "LogNormalize") {
            seurat_objects_integrated <- ScaleData(seurat_objects_integrated)
        }

        ## dimensionality reduction
        seurat_objects_integrated <- RunPCA(seurat_objects_integrated, verbose=FALSE)
        reduction2use <- "pca"
        for (dimensions_nb_param in dimensions_nb_param_vector) {
            print(sprintf("dimension number: %d", dimensions_nb_param))
            dimensions_output_dir <- sprintf("%s/%d_PCs", normalization_output_dir, dimensions_nb_param)
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
            mean_silhouette_df <- cbind(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_param, mean_silhouette_df)
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
                rm(seurat_objects_integrated_hippocampus)
                gc()

                ### get mean silhouette scores: all spots and hippocampus spots
                mean_silhouette_df <- data.frame(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_param, resolution=clustering_resolution, mean_silhouette_seurat_clusters=cluster_mean_silhouette, mean_silhouette_seurat_clusters_hippocampus=hippocampus_cluster_mean_silhouette)
                integration_cluster_mean_silhouette_df <- rbind(integration_cluster_mean_silhouette_df, mean_silhouette_df)

                ### cluster markers, differentially expressed genes across conditions and conserved cell type markers
                marker_analysis(seurat_objects_integrated, default_assay, cluster_markers, load_image, nb_cluster_markers_to_plot, nb_SpatialFeaturePlot_per_page, DE_analysis, conserved_markers, "seurat_clusters", "condition", resolution_output_dir, resolution_output_name)
            }
        }
        rm(seurat_objects_integrated)
        gc()
    }
}

## Harmony
library(harmony)
integration_method <- "Harmony"
integration_output_dir <- sprintf("%s/11-Harmony", qc_filtering_output_dir)
if (! dir.exists(integration_output_dir)) {
    dir.create(integration_output_dir, recursive=TRUE, mode="0775")
}
for (normalization_method in normalization_method_vector) {
    print(sprintf("normalization method: %s", normalization_method))
    for (features_param in features_param_vector) {
        ## normalization
        return_only_var_genes <- ifelse(normalization_method == "SCTransform", FALSE, NA)
        returned_list <- seurat_list_normalization(seurat_objects, normalization_method, features_param, return_only_var_genes, integration_output_dir)
        seurat_objects_norm <- returned_list[["seurat_list"]]
        integration_features_nb <- returned_list[["features_nb"]]
        integration_normalization_method <- returned_list[["norm_method"]]
        default_assay <- returned_list[["assay"]]
        selection_method <- returned_list[["selec"]]
        normalization_output_dir <- returned_list[["out_dir"]]
        normalization_output_name <- returned_list[["out_name"]]
        rm(returned_list)
        gc()
        
        # HVG plots
        hvg_df_list <- lapply(seurat_objects_norm, hvg_plots, assay2use=default_assay, selec=selection_method, nb_hvgs=features_param, out_dir=normalization_output_dir, out_name= "HVGs")
        
        ## merge datasets
        seurat_objects_merge <- merge(seurat_objects_norm[[1]], y=unlist(seurat_objects_norm)[2:nb_sections], add.cell.ids=sections)
        ### orig.ident, lame, zone and condition are not factors in the merged object meta data
        seurat_objects_merge@meta.data$orig.ident <- as.factor(seurat_objects_merge@meta.data$orig.ident)
        seurat_objects_merge@meta.data$lame <- as.factor(seurat_objects_merge@meta.data$lame)
        seurat_objects_merge@meta.data$zone <- as.factor(seurat_objects_merge@meta.data$zone)
        seurat_objects_merge@meta.data$condition <- as.factor(seurat_objects_merge@meta.data$condition)
        ### set variable features: after merging datasets, variable features are not set in the merged Seurat object, which makes RunPCA() fails: Error in PrepDR(object = object, features = features, verbose = verbose): Variable features haven't been set. Run FindVariableFeatures() or provide a vector of feature names.
        #### get the union of the HVGs of each dataset
        if (normalization_method == "SCTransform") {
            seurat_objects_merge_hvgs <- get_HVGs_merged_Seurat_SCT(seurat_objects_merge, features_param)
        } else {
            seurat_objects_merge_hvgs <- get_HVGs_merged_Seurat_LogNormalize(seurat_objects_merge, seurat_objects_norm, features_param)
        }
        VariableFeatures(seurat_objects_merge) <- seurat_objects_merge_hvgs
        
        rm(seurat_objects_norm)
        gc()
        
        if (normalization_method == "LogNormalize") {
            seurat_objects_merge <- ScaleData(seurat_objects_merge)
        }
        
        ### compute PCs 
        seurat_objects_merge <- RunPCA(seurat_objects_merge, verbose=FALSE)
        for (dimensions_nb_param in dimensions_nb_param_vector) {
            print(sprintf("dimension number: %d", dimensions_nb_param))
            dimensions_output_dir <- sprintf("%s/%d_PCs", normalization_output_dir, dimensions_nb_param)
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
            mean_silhouette_df <- cbind(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_param, mean_silhouette_df)
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
                clustering_plots(seurat_objects_integrated, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integrated datasets: %s, all spots", integration_method), resolution_output_dir, sprintf("%s_clustering", resolution_output_name))

                ### hippocampus spots
                seurat_objects_integrated_hippocampus <- subset(seurat_objects_integrated, cells=all_sections_hippocampus_barcodes_Seurat_merge)
                hippocampus_cluster_mean_silhouette <- mean(seurat_objects_integrated_hippocampus$silhouette_seurat_clusters)
                #### plot clusters onto UMAP or onto the tissue section
                clustering_plots(seurat_objects_integrated_hippocampus, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integrated datasets: %s, hippocampus spots", integration_method), resolution_output_dir, sprintf("%s_clustering_hippocampus", resolution_output_name))
                rm(seurat_objects_integrated_hippocampus)
                gc()

                ### get mean silhouette scores: all spots and hippocampus spots
                mean_silhouette_df <- data.frame(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_param, resolution=clustering_resolution, mean_silhouette_seurat_clusters=cluster_mean_silhouette, mean_silhouette_seurat_clusters_hippocampus=hippocampus_cluster_mean_silhouette)
                integration_cluster_mean_silhouette_df <- rbind(integration_cluster_mean_silhouette_df, mean_silhouette_df)

                ### cluster markers, differentially expressed genes across conditions and conserved cell type markers
                marker_analysis(seurat_objects_integrated, default_assay, cluster_markers, load_image, nb_cluster_markers_to_plot, nb_SpatialFeaturePlot_per_page, DE_analysis, conserved_markers, "seurat_clusters", "condition", resolution_output_dir, resolution_output_name)
            }
            rm(seurat_objects_integrated)
            gc()
        }
        rm(seurat_objects_merge)
        gc()
    }
}
detach("package:harmony", unload=TRUE)

## LIGER
library(rliger)
library(SeuratWrappers)
integration_method <- "LIGER"
integration_output_dir <- sprintf("%s/12-LIGER", qc_filtering_output_dir)
if (! dir.exists(integration_output_dir)) {
    dir.create(integration_output_dir, recursive=TRUE, mode="0775")
}
for (normalization_method in normalization_method_vector) {
    print(sprintf("normalization method: %s", normalization_method))
    for (features_param in features_param_vector) {
        ## normalization
        return_only_var_genes <- ifelse(normalization_method == "SCTransform", TRUE, NA) # SCTransform: only return variable genes, otherwise error with RunOptimizeALS() function: Error in dimnames(x) <- dn: length of 'dimnames' [2] not equal to array extent
        returned_list <- seurat_list_normalization(seurat_objects, normalization_method, features_param, return_only_var_genes, integration_output_dir)
        seurat_objects_norm <- returned_list[["seurat_list"]]
        integration_features_nb <- returned_list[["features_nb"]]
        integration_normalization_method <- returned_list[["norm_method"]]
        default_assay <- returned_list[["assay"]]
        selection_method <- returned_list[["selec"]]
        normalization_output_dir <- returned_list[["out_dir"]]
        normalization_output_name <- returned_list[["out_name"]]
        rm(returned_list)
        gc()
        
        # HVG plots
        hvg_df_list <- lapply(seurat_objects_norm, hvg_plots, assay2use=default_assay, selec=selection_method, nb_hvgs=features_param, out_dir=normalization_output_dir, out_name= "HVGs")
        
        ## merge datasets
        seurat_objects_merge <- merge(seurat_objects_norm[[1]], y=unlist(seurat_objects_norm)[2:nb_sections], add.cell.ids=sections)
        ### orig.ident, lame, zone and condition are not factors in the merged object meta data
        seurat_objects_merge@meta.data$orig.ident <- as.factor(seurat_objects_merge@meta.data$orig.ident)
        seurat_objects_merge@meta.data$lame <- as.factor(seurat_objects_merge@meta.data$lame)
        seurat_objects_merge@meta.data$zone <- as.factor(seurat_objects_merge@meta.data$zone)
        seurat_objects_merge@meta.data$condition <- as.factor(seurat_objects_merge@meta.data$condition)
        ### set variable features: after merging datasets, variable features are not set in the merged Seurat object, which makes RunPCA() fails: Error in PrepDR(object = object, features = features, verbose = verbose): Variable features haven't been set. Run FindVariableFeatures() or provide a vector of feature names.
        #### get the union of the HVGs of each dataset
        if (normalization_method == "SCTransform") {
            seurat_objects_merge_hvgs <- get_HVGs_merged_Seurat_SCT(seurat_objects_merge, features_param)
        } else {
            seurat_objects_merge_hvgs <- get_HVGs_merged_Seurat_LogNormalize(seurat_objects_merge, seurat_objects_norm, features_param)
        }
        VariableFeatures(seurat_objects_merge) <- seurat_objects_merge_hvgs
        
        rm(seurat_objects_norm)
        gc()
        
        if (normalization_method == "LogNormalize") {
            seurat_objects_merge <- ScaleData(seurat_objects_merge, split.by="orig.ident", do.center=FALSE)
        }

        for (dimensions_nb_param in dimensions_nb_param_vector) {
            print(sprintf("dimension number: %d", dimensions_nb_param))
            dimensions_output_dir <- sprintf("%s/%d_PCs", normalization_output_dir, dimensions_nb_param)
            if (! dir.exists(dimensions_output_dir)) {
                dir.create(dimensions_output_dir, recursive=TRUE, mode="0775")
            }
            dimensions_output_name <- sprintf("%s_%s_%dPCs", integration_output_name, integration_method, dimensions_nb_param)
            
            ### integration and dimensionality reduction
            seurat_objects_integrated <- RunOptimizeALS(seurat_objects_merge, k=dimensions_nb_param, lambda=5, split.by="orig.ident")
            seurat_objects_integrated <- RunQuantileNorm(seurat_objects_integrated, split.by="orig.ident")
            reduction2use <- "iNMF"

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
            mean_silhouette_df <- cbind(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_param, mean_silhouette_df)
            integration_batch_mean_silhouette_df <- rbind(integration_batch_mean_silhouette_df, mean_silhouette_df)

            ## clustering
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
                clustering_plots(seurat_objects_integrated, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integrated datasets: %s, all spots", integration_method), resolution_output_dir, sprintf("%s_clustering", resolution_output_name))

                ### hippocampus spots
                seurat_objects_integrated_hippocampus <- subset(seurat_objects_integrated, cells=all_sections_hippocampus_barcodes_Seurat_merge)
                hippocampus_cluster_mean_silhouette <- mean(seurat_objects_integrated_hippocampus$silhouette_seurat_clusters)
                #### plot clusters onto UMAP or onto the tissue section
                clustering_plots(seurat_objects_integrated_hippocampus, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integrated datasets: %s, hippocampus spots", integration_method), resolution_output_dir, sprintf("%s_clustering_hippocampus", resolution_output_name))
                rm(seurat_objects_integrated_hippocampus)
                gc()

                ### get mean silhouette scores: all spots and hippocampus spots
                mean_silhouette_df <- data.frame(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_param, resolution=clustering_resolution, mean_silhouette_seurat_clusters=cluster_mean_silhouette, mean_silhouette_seurat_clusters_hippocampus=hippocampus_cluster_mean_silhouette)
                integration_cluster_mean_silhouette_df <- rbind(integration_cluster_mean_silhouette_df, mean_silhouette_df)
                
                ### cluster markers, differentially expressed genes across conditions and conserved cell type markers
                marker_analysis(seurat_objects_integrated, default_assay, cluster_markers, load_image, nb_cluster_markers_to_plot, nb_SpatialFeaturePlot_per_page, DE_analysis, conserved_markers, "seurat_clusters", "condition", resolution_output_dir, resolution_output_name)
            }
            rm(seurat_objects_integrated)
            gc()
        }
        rm(seurat_objects_merge)
        gc()
    }
}
rm(seurat_objects)
gc()
write.csv(integration_cluster_mean_silhouette_df, file=sprintf("%s/featuresnumber_integration_PCs_resolution_cluster_mean_silhouette.csv", qc_filtering_output_dir), quote=FALSE, row.names=FALSE)
integration_cluster_mean_silhouette_df$integration <- factor(integration_cluster_mean_silhouette_df$integration, levels=c("no integration", "Seurat", "Harmony", "LIGER"))
integration_cluster_mean_silhouette_df$integration_normalization <- sprintf("%s %s", integration_cluster_mean_silhouette_df$integration, integration_cluster_mean_silhouette_df$normalization)
integration_cluster_mean_silhouette_df$integration_normalization <- factor(integration_cluster_mean_silhouette_df$integration_normalization, levels=c("no integration SCTransform", "no integration LogNormalize", "Seurat SCTransform", "Seurat LogNormalize", "Harmony SCTransform", "Harmony LogNormalize", "LIGER SCTransform", "LIGER LogNormalize"))
levels(integration_cluster_mean_silhouette_df$integration_normalization) <- unlist(lapply(levels(integration_cluster_mean_silhouette_df$integration_normalization), function(x) { return(gsub(" SCTransform", "\nSCTransform", x)) }))
levels(integration_cluster_mean_silhouette_df$integration_normalization) <- unlist(lapply(levels(integration_cluster_mean_silhouette_df$integration_normalization), function(x) { return(gsub(" LogNormalize", "\nLogNormalize", x)) }))
integration_cluster_mean_silhouette_df$features <- factor(integration_cluster_mean_silhouette_df$features, levels=features_param_vector)
integration_cluster_mean_silhouette_df$dimensions <- as.factor(integration_cluster_mean_silhouette_df$dimensions)
integration_cluster_mean_silhouette_df$resolution <- as.factor(integration_cluster_mean_silhouette_df$resolution)
pdf(sprintf("%s/featuresnumber_integration_PCs_resolution_cluster_mean_silhouette.pdf", qc_filtering_output_dir))
for (one_category in c("seurat_clusters", "seurat_clusters_hippocampus")) {
    if (one_category == "seurat_clusters") {
        plot_title_category <- "all spots"
    } else {
        if (one_category == "seurat_clusters_hippocampus") {
            plot_title_category <- "hippocampus spots"
        }
    }
    p <- ggplot(integration_cluster_mean_silhouette_df, aes(x=resolution, y=.data[[sprintf("mean_silhouette_%s", one_category)]], fill=integration_normalization)) +
        geom_bar(stat="identity", position="dodge") +
        facet_grid(dimensions~features) +
        labs(title=sprintf("Mean silhouette: %s", plot_title_category), x="Resolution", y="Mean silhouette", fill="Integration\nnormalization\nmethods") +
        theme_bw() +
        theme(legend.position="bottom") +
        theme(axis.text.x=element_text(angle=60, hjust=1)) +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
    
    for (one_dimensions in levels(integration_cluster_mean_silhouette_df$dimensions)) {
        dim_df2ggplot <- integration_cluster_mean_silhouette_df[which(integration_cluster_mean_silhouette_df$dimensions==one_dimensions),]
        p <- ggplot(dim_df2ggplot, aes(x=resolution, y=.data[[sprintf("mean_silhouette_%s", one_category)]], fill=integration_normalization)) +
            geom_bar(stat="identity", position="dodge") +
            facet_grid(~features) +
            labs(title=sprintf("Mean silhouette: %s, dimensions: %s", plot_title_category, one_dimensions), x="Resolution", y="Mean silhouette", fill="Integration\nnormalization\nmethods") +
            theme_bw() +
            theme(legend.position="bottom") +
            theme(axis.text.x=element_text(angle=60, hjust=1)) +
            theme(panel.border=element_rect(color="grey50"))
        print(p)
        
        for (one_features in levels(dim_df2ggplot$features)) {
            df2ggplot <- dim_df2ggplot[which(dim_df2ggplot$features==one_features),]
            p <- ggplot(df2ggplot, aes(x=resolution, y=.data[[sprintf("mean_silhouette_%s", one_category)]], fill=integration_normalization)) +
                geom_bar(stat="identity", position="dodge") +
                labs(title=sprintf("Mean silhouette: %s, dimensions: %s, features: %s", plot_title_category, one_dimensions, one_features), x="Resolution", y="Mean silhouette", fill="Integration\nnormalization\nmethods") +
                theme_bw() +
                theme(legend.position="bottom") +
                theme(axis.text.x=element_text(angle=60, hjust=1)) +
                theme(panel.border=element_rect(color="grey50"))
            print(p)
            
        }
    }
}
dev.off()

write.csv(integration_batch_mean_silhouette_df, file=sprintf("%s/featuresnumber_integration_PCs_resolution_batch_mean_silhouette.csv", qc_filtering_output_dir), quote=FALSE, row.names=FALSE)
integration_batch_mean_silhouette_df$integration <- factor(integration_batch_mean_silhouette_df$integration, levels=c("no integration", "Seurat", "Harmony", "LIGER"))
integration_batch_mean_silhouette_df$integration_normalization <- sprintf("%s %s", integration_batch_mean_silhouette_df$integration, integration_batch_mean_silhouette_df$normalization)
integration_batch_mean_silhouette_df$integration_normalization <- factor(integration_batch_mean_silhouette_df$integration_normalization, levels=c("no integration SCTransform", "no integration LogNormalize", "Seurat SCTransform", "Seurat LogNormalize", "Harmony SCTransform", "Harmony LogNormalize", "LIGER SCTransform", "LIGER LogNormalize"))
levels(integration_batch_mean_silhouette_df$integration_normalization) <- unlist(lapply(levels(integration_batch_mean_silhouette_df$integration_normalization), function(x) { return(gsub(" SCTransform", "\nSCTransform", x)) }))
levels(integration_batch_mean_silhouette_df$integration_normalization) <- unlist(lapply(levels(integration_batch_mean_silhouette_df$integration_normalization), function(x) { return(gsub(" LogNormalize", "\nLogNormalize", x)) }))
integration_batch_mean_silhouette_df$features <- factor(integration_batch_mean_silhouette_df$features, levels=features_param_vector)
integration_batch_mean_silhouette_df$dimensions <- as.factor(integration_batch_mean_silhouette_df$dimensions)
pdf(sprintf("%s/featuresnumber_integration_PCs_resolution_batch_mean_silhouette.pdf", qc_filtering_output_dir))
for (one_category in c("lame", "zone", "condition", "orig.ident")) {
    p <- ggplot(integration_batch_mean_silhouette_df, aes(x=dimensions, y=.data[[sprintf("mean_silhouette_%s", one_category)]], fill=integration_normalization)) +
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
    p <- ggplot(integration_batch_mean_silhouette_df, aes(x=integration, y=.data[[sprintf("mean_silhouette_%s", one_category)]], fill=integration_normalization)) +
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




