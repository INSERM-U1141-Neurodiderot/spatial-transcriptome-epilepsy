.libPaths(c("/home/christophe.lepriol/NeuroDev_ADD/R/r_4.1.0", .libPaths()))
library(SpatialExperiment)
library(ggspavis)
library(scater)
library(ggplot2)
library(cowplot)
library(patchwork) # for wrap_plots() function


##############
# Parameters #
##############
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"
genome_name <- "NCBIRefSeq108_NCBIRefSeq108GTF"
samples <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")
#sample <- "A_L1_S1"
#samples <- c("A_L1_S1", "A_L2_S5")

# vectors and list for dataset plots
## library size and number of expressed genes densities
sample_lib_size_vector <- count_sum_vector <- detected_vector <- c()
## discarded spots barplots
sample_discarded_vector <- metric_vector <- threshold_vector <- count_vector <- c()
## plot QC list
plot_qc_list <- list()

for (sample in samples) {
    print(sprintf("sample: %s", sample))
    
    output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s/00-QC", work_dir, genome_name, sample)
    if (! dir.exists(output_dir)) {
      dir.create(output_dir, recursive=TRUE, mode="0775")
    }
    space_ranger_output_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs", work_dir, genome_name, sample)

    # load data
    ## Visium V1 barcode coordinates
    #visium_v1_coordinates_file <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/data/barcodes/visium-v1_coordinates.txt"
    #visium_v1_coordinates_df <- read.table(visium_v1_coordinates_file, sep="\t", quote="")

    ## expression matrices
    ### Space Ranger
    #h5_file <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs/raw_feature_bc_matrix.h5", work_dir, genome_name, sample)
    #h5_data <- Read10X_h5(h5_file)
    ### ST Pipeline
    st_pipeline_matrix_file <- sprintf("%s/10-ST_analysis/00-ST_Pipeline/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s/10-Pipeline/%s_stdata.tsv", work_dir, genome_name, sample, sample)
    st_pipeline_matrix <- read.table(st_pipeline_matrix_file, sep="\t", header=TRUE, quote="", row.names=1)

    # build SpatialExperiment object
    ## spatial coordinates
    #visium_v1_coordinates_df <- read.table(visium_v1_coordinates_file, sep="\t", quote="")
    #colnames(visium_v1_coordinates_df) <- c("barcode", "x", "y")
    #rownames(visium_v1_coordinates_df) <- visium_v1_coordinates_df$barcode
    spatial_coordinates_file <- file.path(space_ranger_output_dir, "spatial", "tissue_positions_list.csv")
    spatial_coordinates_df <- read.csv(spatial_coordinates_file, header=FALSE, quote="")
    colnames(spatial_coordinates_df) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    rownames(spatial_coordinates_df) <- spatial_coordinates_df$barcode

    ### matrix: convert coordinates to barcodes
    #barcodes <- unlist(lapply(rownames(st_pipeline_matrix), function(x, coord2barcode=visium_v1_coordinates_df) {
    #    coordinates <- unlist(strsplit(x, "x"))
    #    return(coord2barcode[which(coord2barcode$x==as.integer(coordinates[1]) & coord2barcode$y==as.integer(coordinates[2])), "barcode"])
    #}))
    barcodes <- unlist(lapply(rownames(st_pipeline_matrix), function(x, coord2barcode=spatial_coordinates_df) {
        coordinates <- unlist(strsplit(x, "x"))
        return(coord2barcode[which(coord2barcode$array_col==as.integer(coordinates[1])-1 & coord2barcode$array_row==as.integer(coordinates[2])-1), "barcode"])
    }))
    rownames(st_pipeline_matrix) <- barcodes

    gene_data <- data.frame(gene_name=colnames(st_pipeline_matrix))

    ### image data
    img <- readImgData(imageSources=file.path(space_ranger_output_dir, "spatial", "tissue_lowres_image.png"), scaleFactors=file.path(space_ranger_output_dir, "spatial", "scalefactors_json.json"), sample_id=sample)

    #spe <- SpatialExperiment(assay=list(counts=t(st_pipeline_matrix)), rowData=gene_data, colData=spatial_coordinates_df[rownames(st_pipeline_matrix),], imgData=img, spatialCoordsNames=c("pxl_col_in_fullres", "pxl_row_in_fullres"), sample_id=sample)
    spe <- SpatialExperiment(assay=list(counts=t(st_pipeline_matrix)), rowData=gene_data, colData=spatial_coordinates_df[rownames(st_pipeline_matrix),], imgData=img, spatialDataNames=c("barcode", "in_tissue", "array_row", "array_col"), spatialCoordsNames=c("pxl_col_in_fullres", "pxl_row_in_fullres"), sample_id=sample)

    pdf(sprintf("%s/%s_lib_size_expressed_genes_loess_curve.pdf", output_dir, sample))
    plotSpots(spe, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    # keep spot over tissue
    spe_in_tissue <- spe[, colData(spe)$in_tissue==1]

    # identify mitochondrial genes
    is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe_in_tissue)$gene_name)
    table(is_mito)
    #rowData(spe)$gene_name[is_mito]

    # calculate per-spot QC metrics and store in colData
    spe_in_tissue <- addPerCellQC(spe_in_tissue)
    head(colData(spe_in_tissue))

    hist(colData(spe_in_tissue)$sum, breaks=20, main="Sum of counts per spot", xlab="Sum of counts")
    count_sum_vector <- c(count_sum_vector, colData(spe_in_tissue)$sum)
    hist(colData(spe_in_tissue)$detected, breaks=20, main="Detected genes per spot", xlab="Detected genes")
    detected_vector <- c(detected_vector, colData(spe_in_tissue)$detected)
    sample_lib_size_vector <- c(sample_lib_size_vector, rep(sample, length(colData(spe_in_tissue)$detected)))
    ## plot sum of counts and number of detected genes per spot
    print(plotSpots(spe_in_tissue, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="sum", palette=c("white", "black"), size=2))
    print(plotVisium(spe_in_tissue, fill="sum", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    print(plotSpots(spe_in_tissue, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="detected", palette=c("white", "black"), size=2))
    print(plotVisium(spe_in_tissue, fill="detected", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))

    # selecting thresholds
    sample_plot_qc_list <- list()
    ## library size
    ### set threshold and check the number of spots below this threshold
    thresholds <- c(1000, 2000, 3000, 4000, 5000)
    for (threshold in thresholds) {
        qc_lib_size <- colData(spe_in_tissue)$sum < threshold
        table(qc_lib_size)
        count_vector <- c(count_vector, unname(table(qc_lib_size)["TRUE"]))
        colData(spe_in_tissue)$qc_lib_size <- qc_lib_size
        qc_colname <- sprintf("qc_lib_size_%d", threshold)
        colnames(colData(spe_in_tissue))[which(colnames(colData(spe_in_tissue)) == "qc_lib_size")] <- qc_colname
        plot_qc <- plotQC(spe_in_tissue, type="spots", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", discard=qc_colname)
        plot_qc$layers[[1]]$aes_params$size <- 2
        print(plot_qc)
        sample_plot_qc_list[[qc_colname]] <- plot_qc
    }
    threshold_vector <- c(threshold_vector, thresholds)
    metric_vector <- c(metric_vector, rep("library size", length(thresholds)))
    sample_discarded_vector <- c(sample_discarded_vector, rep(sample, length(thresholds)))
    ## number of expressed features
    ### set threshold and check the number of spots below this threshold
    thresholds <- c(500, 1000, 1500, 2000, 2500)
    for (threshold in thresholds) {
        qc_detected <- colData(spe_in_tissue)$detected < threshold
        table(qc_detected)
        count_vector <- c(count_vector, unname(table(qc_detected)["TRUE"]))
        colData(spe_in_tissue)$qc_detected <- qc_detected
        qc_colname <- sprintf("qc_detected_%d", threshold)
        colnames(colData(spe_in_tissue))[which(colnames(colData(spe_in_tissue)) == "qc_detected")] <- qc_colname
        plot_qc <- plotQC(spe_in_tissue, type="spots", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", discard=qc_colname)
        plot_qc$layers[[1]]$aes_params$size <- 2
        print(plot_qc)
        sample_plot_qc_list[[qc_colname]] <- plot_qc
    }
    threshold_vector <- c(threshold_vector, thresholds)
    metric_vector <- c(metric_vector, rep("expressed genes", length(thresholds)))
    sample_discarded_vector <- c(sample_discarded_vector, rep(sample, length(thresholds)))
    ## remove low-quality spots
    ### number of discarded spots for each metric
    apply(cbind(colData(spe_in_tissue)$qc_lib_size_2000, colData(spe_in_tissue)$qc_detected_1000), 2, sum)
    ### combined set of discarded spots
    discard <- colData(spe_in_tissue)$qc_lib_size_2000 | colData(spe_in_tissue)$qc_detected_1000
    table(discard)
    count_vector <- c(count_vector, unname(table(discard)["TRUE"]))
    threshold_vector <- c(threshold_vector, NA)
    metric_vector <- c(metric_vector, "discard")
    sample_discarded_vector <- c(sample_discarded_vector, sample)
    ### store in object
    colData(spe_in_tissue)$discard <- discard
    print(plotQC(spe_in_tissue, type="scatter", metric_x="sum", metric_y="detected", threshold_x=2000, threshold_y=1000))
    plot_qc <- plotQC(spe_in_tissue, type="spots", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", discard="discard")
    plot_qc$layers[[1]]$aes_params$size <- 2
    print(plot_qc)
    sample_plot_qc_list[["discard"]] <- plot_qc
    ### remove combined set of low-quality spots
    spe_in_tissue <- spe_in_tissue[, !colData(spe_in_tissue)$discard]
    dim(spe_in_tissue)
    print(plotSpots(spe_in_tissue, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="sum", palette=c("white", "black"), size=2))
    print(plotSpots(spe_in_tissue, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="detected", palette=c("white", "black"), size=2))
    #dev.off()
    
    
    # filter spots with high values of sum of read counts without an increase of number of expressed genes
    ## reproduce loess regression plotQC() function
    df <- data.frame(barcode=rownames(colData(spe_in_tissue)), sum=colData(spe_in_tissue)$sum, detected=colData(spe_in_tissue)$detected)
    loess_fit <- loess(detected ~ sum, data=df, span=0.5)
    loess_df <- cbind(df, fitted=loess_fit$fitted, sample=rep(sample, length(loess_fit$fitted)))
    loess_df <- loess_df[order(loess_df$sum),]
    loess_curve <- ggplot(loess_df, aes(x=sum)) +
        geom_point(aes(y=detected)) +
        geom_line(aes(y=fitted), color="blue", size=1.5) +
        labs(title="Loess regression", x="Read count per spot", y="Number of expressed genes per spot") +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50"))
    print(loess_curve)
    ## compute loess curve slope between two consecutive points
    ### slope approximation: difference formulas to compute a numerical derivative at any point
    #### irregularly spaced data: slope=diff(y)/diff(x)
    ##### https://blogs.sas.com/content/iml/2018/07/05/derivatives-nonparametric-regression.html
    ##### https://blogs.sas.com/content/iml/2012/05/02/the-dif-function.html
    sum_diff <- diff(loess_df$sum)
    fitted_diff <- diff(loess_df$fitted)
    slope <- fitted_diff / sum_diff
    #df2ggplot_3 <- data.frame(index=index_vector, slope=slope, window=rep("diff", length(df_diff)))
    loess_fit_slope_df <- cbind(loess_df[-1,], slope=slope)
    rownames(loess_fit_slope_df) <- loess_fit_slope_df$barcode
    loess_curve_slope <- ggplot(loess_fit_slope_df, aes(x=sum, y=slope)) +
        geom_point() +
        geom_line() +
        geom_hline(yintercept=c(0, 0.01, 0.02), alpha=0.5, linetype="longdash") +
        labs(title="Expressed genes read counts loess regression slope", x="Read count per spot", y="Slope") +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50"))
    print(loess_curve_slope)
    ### add fitted and slope values to colData(spe_in_tissue)
    colData(spe_in_tissue)$fitted <- loess_fit_slope_df[colData(spe_in_tissue)$barcode, "fitted"]
    colData(spe_in_tissue)$slope <- loess_fit_slope_df[colData(spe_in_tissue)$barcode, "slope"]
    
    thresholds <- c(0, 0.01, 0.02)
    for (threshold in thresholds) {
        ### loess curve slope only
        qc_colname <- sprintf("slope_%s", sub("\\.", "_", threshold))
        loess_fit_slope_df <- cbind(loess_fit_slope_df, slope_threshold=loess_fit_slope_df$slope < threshold)
        count_vector <- c(count_vector, unname(table(loess_fit_slope_df$slope_threshold)["TRUE"]))
        colnames(loess_fit_slope_df)[which(colnames(loess_fit_slope_df) == "slope_threshold")] <- qc_colname
        colData(spe_in_tissue)$slope_threshold <- loess_fit_slope_df[colData(spe_in_tissue)$barcode, qc_colname]
        colnames(colData(spe_in_tissue))[which(colnames(colData(spe_in_tissue)) == "slope_threshold")] <- qc_colname
        plot_qc <- plotQC(spe_in_tissue, type="spots", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", discard=qc_colname)
        plot_qc$layers[[1]]$aes_params$size <- 2
        print(plot_qc)
        sample_plot_qc_list[[qc_colname]] <- plot_qc
        ### filter spots based on loess curve slope and if detected < fitted
        qc_colname <- sprintf("slope_%s_fitted", sub("\\.", "_", threshold))
        loess_fit_slope_df <- cbind(loess_fit_slope_df, slope_threshold_fitted=loess_fit_slope_df$slope < threshold & loess_fit_slope_df$detected < loess_fit_slope_df$fitted)
        count_vector <- c(count_vector, unname(table(loess_fit_slope_df$slope_threshold_fitted)["TRUE"]))
        colnames(loess_fit_slope_df)[which(colnames(loess_fit_slope_df) == "slope_threshold_fitted")] <- qc_colname
        colData(spe_in_tissue)$slope_threshold_fitted <- loess_fit_slope_df[colData(spe_in_tissue)$barcode, qc_colname]
        colnames(colData(spe_in_tissue))[which(colnames(colData(spe_in_tissue)) == "slope_threshold_fitted")] <- qc_colname
        plot_qc <- plotQC(spe_in_tissue, type="spots", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", discard=qc_colname)
        plot_qc$layers[[1]]$aes_params$size <- 2
        print(plot_qc)
        sample_plot_qc_list[[qc_colname]] <- plot_qc
        threshold_vector <- c(threshold_vector, rep(threshold, 2))
    }
    plot_qc_list[[sample]] <- sample_plot_qc_list
    metric_vector <- c(metric_vector, rep(c("loess slope", "loess slope fitted"), length(thresholds)))
    sample_discarded_vector <- c(sample_discarded_vector, rep(sample, 2*length(thresholds)))
    dev.off()
    write.csv(colData(spe_in_tissue), file=sprintf("%s/%s_lib_size_expressed_genes_loess_curve.csv", output_dir, sample), quote=FALSE, row.names=FALSE)
}

# dataset plots
output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/10-Dataset", work_dir, genome_name)
if (! dir.exists(output_dir)) {
    dir.create(output_dir, recursive=TRUE, mode="0775")
}

## library size and number of expressed genes per spot distributions
pdf(sprintf("%s/lib_size_expressed_genes_distributions.pdf", output_dir))
df2ggplot <- data.frame(sample=sample_lib_size_vector, count_sum=count_sum_vector, detected=detected_vector)
write.csv(df2ggplot, file=sprintf("%s/lib_size_expressed_genes_distributions.csv", output_dir), quote=FALSE, row.names=FALSE)
density_plot <- ggplot(df2ggplot, aes(x=count_sum, color=sample)) +
    geom_line(stat="density") +
    geom_vline(xintercept=2000, alpha=0.5, linetype="longdash") +
    labs(title="Count sum per spot", x="Library size", y="Density", color="Sample") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50"))
print(density_plot)
density_plot <- ggplot(df2ggplot, aes(x=detected, color=sample)) +
    geom_line(stat="density") +
    geom_vline(xintercept=1000, alpha=0.5, linetype="longdash") +
    labs(title="Number of expressed genes per spot", x="Number of expressed genes", y="Density", color="Sample") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50"))
print(density_plot)
dev.off()

## discarded spots: library size, number of expressed genes
pdf(sprintf("%s/lib_size_expressed_genes_discarded_spots.pdf", output_dir))
### library size
#### number of discarded spots
df2ggplot <- data.frame(sample=sample_discarded_vector, metric=metric_vector, threshold=factor(threshold_vector), count=count_vector)
write.csv(df2ggplot, file=sprintf("%s/lib_size_expressed_genes_discarded_spots.csv", output_dir), quote=FALSE, row.names=FALSE)
barplot <- ggplot(df2ggplot[which(df2ggplot$metric=="library size"),], aes(x=sample, y=count, fill=threshold)) +
    geom_bar(position="dodge", stat="identity") +
    labs(title="Number of discarded spots genes per spot: library size", x="Sample", y="# spots", fill="Threshold") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.text.x=element_text(angle=30, hjust=1))
print(barplot)
#### all samples QC plots on a single page
##### modify plots to adjust with multiple plots on a single page: decrease point size, remove legend, change title
for (threshold in c(1000, 2000, 3000, 4000, 5000)) {
    qc_plots_dataset_list <- list()
    qc_name <- sprintf("qc_lib_size_%d", threshold)
    for (sample in names(plot_qc_list)) {
        sample_plot <- plot_qc_list[[sample]][[qc_name]]
        sample_plot$layers[[1]]$aes_params$size <- 0.01
        sample_plot <- sample_plot + theme(legend.position="none") + labs(title=sample)
        qc_plots_dataset_list[[sample]] <- sample_plot
    }
    print(wrap_plots(qc_plots_dataset_list) + plot_annotation(title=sprintf("Library size threshold: %d", threshold)))
}

### number of expressed genes
#### number of discarded spots
barplot <- ggplot(df2ggplot[which(df2ggplot$metric=="expressed genes"),], aes(x=sample, y=count, fill=threshold)) +
    geom_bar(position="dodge", stat="identity") +
    labs(title="Number of discarded spots genes per spot: expressed genes", x="Sample", y="# spots", fill="Threshold") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.text.x=element_text(angle=30, hjust=1))
print(barplot)
#### all samples QC plots on a single page
##### modify plots to adjust with multiple plots on a single page: decrease point size, remove legend, change title
for (threshold in c(500, 1000, 1500, 2000, 2500)) {
    qc_plots_dataset_list <- list()
    qc_name <- sprintf("qc_detected_%d", threshold)
    for (sample in names(plot_qc_list)) {
        sample_plot <- plot_qc_list[[sample]][[qc_name]]
        sample_plot$layers[[1]]$aes_params$size <- 0.01
        sample_plot <- sample_plot + theme(legend.position="none") + labs(title=sample)
        qc_plots_dataset_list[[sample]] <- sample_plot
    }
    print(wrap_plots(qc_plots_dataset_list) + plot_annotation(title=sprintf("Number of expressed genes threshold: %d", threshold)))
}

### library size and number of expressed genes
#### number of discarded spots
barplot <- ggplot(df2ggplot[which(df2ggplot$metric=="discard"),], aes(x=sample, y=count)) +
    geom_bar(stat="identity") +
    labs(title="Number of discarded spots genes per spot: expressed genes", x="Sample", y="# spots") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.text.x=element_text(angle=30, hjust=1))
print(barplot)
#### all samples QC plots on a single page
##### modify plots to adjust with multiple plots on a single page: decrease point size, remove legend, change title
qc_plots_dataset_list <- list()
for (sample in names(plot_qc_list)) {
    sample_plot <- plot_qc_list[[sample]][["discard"]]
    sample_plot$layers[[1]]$aes_params$size <- 0.01
    sample_plot <- sample_plot + theme(legend.position="none") + labs(title=sample)
    qc_plots_dataset_list[[sample]] <- sample_plot
}
print(wrap_plots(qc_plots_dataset_list) + plot_annotation(title="Discarded spots"))

### high values of sum of read counts without an increase of number of expressed genes
#### number of discarded spots
##### loess curve slope only
barplot <- ggplot(df2ggplot[which(df2ggplot$metric=="loess slope"),], aes(x=sample, y=count, fill=threshold)) +
    geom_bar(position="dodge", stat="identity") +
    labs(title="Number of discarded spots genes per spot: loess slope", x="Sample", y="# spots", fill="Threshold") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.text.x=element_text(angle=30, hjust=1))
print(barplot)
##### filter spots based on loess curve slope and if detected < fitted
barplot <- ggplot(df2ggplot[which(df2ggplot$metric=="loess slope fitted"),], aes(x=sample, y=count, fill=threshold)) +
    geom_bar(position="dodge", stat="identity") +
    labs(title="Number of discarded spots genes per spot: loess slope and detected < fitted", x="Sample", y="# spots", fill="Threshold") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.text.x=element_text(angle=30, hjust=1))
print(barplot)
#### all samples QC plots on a single page
##### modify plots to adjust with multiple plots on a single page: decrease point size, remove legend, change title
for (threshold in c(0, 0.01, 0.02)) {
    ###### loess curve slope only
    qc_plots_dataset_list <- list()
    qc_name <- sprintf("slope_%s", sub("\\.", "_", threshold))
    for (sample in names(plot_qc_list)) {
        sample_plot <- plot_qc_list[[sample]][[qc_name]]
        sample_plot$layers[[1]]$aes_params$size <- 0.01
        sample_plot <- sample_plot + theme(legend.position="none") + labs(title=sample)
        qc_plots_dataset_list[[sample]] <- sample_plot
    }
    print(wrap_plots(qc_plots_dataset_list) + plot_annotation(title=sprintf("Loess curve slope threshold: %s", threshold)))
    ###### filter spots based on loess curve slope and if detected < fitted
    qc_plots_dataset_list <- list()
    qc_name <- sprintf("slope_%s_fitted", sub("\\.", "_", threshold))
    for (sample in names(plot_qc_list)) {
        sample_plot <- plot_qc_list[[sample]][[qc_name]]
        sample_plot$layers[[1]]$aes_params$size <- 0.01
        sample_plot <- sample_plot + theme(legend.position="none") + labs(title=sample)
        qc_plots_dataset_list[[sample]] <- sample_plot
    }
    print(wrap_plots(qc_plots_dataset_list) + plot_annotation(title=sprintf("Loess curve slope threshold: %s and detected < fitted", threshold)))
}
dev.off()












########################################################################################################################

# filter spots with high values of sum of read counts without an increase of number of expressed genes
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"
genome_name <- "NCBIRefSeq108_NCBIRefSeq108GTF"
samples <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")
#samples <- c("A_L1_S1", "A_L2_S5")
#sample <- "A_L2_S5"
output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/NCBIRefSeq108_NCBIRefSeq108GTF/00-Visium_recommended/10-Dataset/00-QC", work_dir)
if (! dir.exists(output_dir)) {
    dir.create(output_dir, recursive=TRUE, mode="0775")
}
loess_fit_df <- data.frame()
loess_fit_slope_df <- data.frame()
#pdf(sprintf("%s/expressed_genes_read_counts_loess.pdf", output_dir))
for (sample in samples) {
    print(sprintf("sample: %s", sample))
    space_ranger_output_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs", work_dir, genome_name, sample)

    # load data
    ## expression matrices
    ### ST Pipeline
    st_pipeline_matrix_file <- sprintf("%s/10-ST_analysis/00-ST_Pipeline/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s/10-Pipeline/%s_stdata.tsv", work_dir, genome_name, sample, sample)
    st_pipeline_matrix <- read.table(st_pipeline_matrix_file, sep="\t", header=TRUE, quote="", row.names=1)

    # build SpatialExperiment object
    ## spatial coordinates
    spatial_coordinates_file <- file.path(space_ranger_output_dir, "spatial", "tissue_positions_list.csv")
    spatial_coordinates_df <- read.csv(spatial_coordinates_file, header=FALSE, quote="")
    colnames(spatial_coordinates_df) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    rownames(spatial_coordinates_df) <- spatial_coordinates_df$barcode

    ### matrix: convert coordinates to barcodes
    barcodes <- unlist(lapply(rownames(st_pipeline_matrix), function(x, coord2barcode=spatial_coordinates_df) {
        coordinates <- unlist(strsplit(x, "x"))
        return(coord2barcode[which(coord2barcode$array_col==as.integer(coordinates[1])-1 & coord2barcode$array_row==as.integer(coordinates[2])-1), "barcode"])
    }))
    rownames(st_pipeline_matrix) <- barcodes

    gene_data <- data.frame(gene_name=colnames(st_pipeline_matrix))

    ### image data
    img <- readImgData(imageSources=file.path(space_ranger_output_dir, "spatial", "tissue_lowres_image.png"), scaleFactors=file.path(space_ranger_output_dir, "spatial", "scalefactors_json.json"), sample_id=sample)

    spe <- SpatialExperiment(assay=list(counts=t(st_pipeline_matrix)), rowData=gene_data, colData=spatial_coordinates_df[rownames(st_pipeline_matrix),], imgData=img, spatialDataNames=c("barcode", "in_tissue", "array_row", "array_col"), spatialCoordsNames=c("pxl_col_in_fullres", "pxl_row_in_fullres"), sample_id=sample)

    # keep spot over tissue
    spe_in_tissue <- spe[, colData(spe)$in_tissue==1]

    # identify mitochondrial genes
    is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe_in_tissue)$gene_name)
    table(is_mito)
    #rowData(spe)$gene_name[is_mito]

    # calculate per-spot QC metrics and store in colData
    spe_in_tissue <- addPerCellQC(spe_in_tissue)

    # reproduce loess regression plotQC() function
    df <- data.frame(barcode=rownames(colData(spe_in_tissue)), sum=colData(spe_in_tissue)$sum, detected=colData(spe_in_tissue)$detected)
    loess_fit <- loess(detected ~ sum, data=df, span=0.5)
    #df_smooth <- predict(df_loess)
    loess_df <- cbind(df, fitted=loess_fit$fitted, sample=rep(sample, length(loess_fit$fitted)))
    loess_df <- loess_df[order(loess_df$sum),]
    loess_fit_df <- rbind(loess_fit_df, loess_df)
    
    loess_curve <- ggplot(loess_df, aes(x=sum)) +
        geom_point(aes(y=detected)) +
        geom_line(aes(y=fitted), color="blue", size=2) +
        labs(title="Loess regression", x="Read count per spot", y="Number of expressed genes per spot") +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50"))
        
    ## compute loess curve slope between two consecutive points
    ### slope approximation: difference formulas to compute a numerical derivative at any point
    #### irregularly spaced data: slope=diff(y)/diff(x)
    ##### https://blogs.sas.com/content/iml/2018/07/05/derivatives-nonparametric-regression.html
    ##### https://blogs.sas.com/content/iml/2012/05/02/the-dif-function.html
    sum_diff <- diff(loess_df$sum)
    fitted_diff <- diff(loess_df$fitted)
    slope <- fitted_diff / sum_diff
    #df2ggplot_3 <- data.frame(index=index_vector, slope=slope, window=rep("diff", length(df_diff)))
    sample_loess_fit_slope_df <- cbind(loess_df[-1,], slope=slope)
    loess_fit_slope_df <- rbind(loess_fit_slope_df, sample_loess_fit_slope_df)
    
    loess_curve_slope <- ggplot(sample_loess_fit_slope_df, aes(x=sum, y=slope)) +
        geom_point() +
        geom_line() +
        geom_hline(yintercept=c(0, 0.01, 0.02), alpha=0.5, linetype="longdash") +
        labs(title="Expressed genes read counts loess regression slope", x="Read count per spot", y="Slope") +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50"))
    
    multiplot <- ggdraw() +
        draw_plot(loess_curve, 0, 0, 0.50, 0.96) +
        draw_plot(loess_curve_slope, 0.50, 0, 0.50, 0.96) +
        draw_label(sample, x=0.02, y=0.97, hjust=0, vjust=0)
    print(multiplot)

}
# all samples plots
p <- ggplot(loess_fit_df, aes(x=sum)) +
    geom_point(aes(y=detected)) +
    geom_line(aes(y=fitted), color="blue", size=2) +
    facet_wrap(~sample, scales="free") +
    labs(title="Loess regression", x="Read count per spot", y="Number of expressed genes per spot") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50"))
print(p)
p <- ggplot(loess_fit_slope_df, aes(x=sum, y=slope, color=sample)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept=c(0, 0.01, 0.02), alpha=0.5, linetype="longdash") +
    labs(title="Expressed genes read counts loess regression slope", x="Read count per spot", y="Slope", color="Sample") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50"))
print(p)
dev.off()

# number of filtered spots per sample according to threshold values
threshold_vector <- sample_vector <- filtered_spots_vector <- c()
for (threshold in c(0, 0.01, 0.02)) {
    filtered_spots <- table(loess_fit_slope_df[which(loess_fit_slope_df$slope < threshold), "sample"])
    filtered_spots_vector <- c(filtered_spots_vector, as.numeric(unname(filtered_spots)))
    sample_vector <- c(sample_vector, names(filtered_spots))
    threshold_vector <- c(threshold_vector, rep(threshold, length(filtered_spots)))
}
filtered_spots_df2ggplot <- data.frame(threshold=as.factor(threshold_vector), sample=sample_vector, spots=filtered_spots_vector)
p <- ggplot(filtered_spots_df2ggplot, aes(x=sample, y=spots, fill=threshold)) +
    geom_bar(stat="identity", position=position_dodge()) +
    labs(title="Number of filtered spots according to loess curve slope", x="Sample", y="#spots", fill="Threshold") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.text.x=element_text(angle=20, hjust=1))
print(p)


one_sample <- "D_L4_S16"
one_sample_loess_fit_slope_df <- loess_fit_slope_df[which(loess_fit_slope_df$sample==one_sample),]
threshold <- 0.04
one_sample_loess_fit_slope_df <- cbind(one_sample_loess_fit_slope_df, slope_threshold=one_sample_loess_fit_slope_df$slope < threshold)
rownames(one_sample_loess_fit_slope_df) <- one_sample_loess_fit_slope_df$barcode
colData(spe_in_tissue)$slope_threshold <- one_sample_loess_fit_slope_df[colData(spe_in_tissue)$barcode, "slope_threshold"]
qc_colname <- sprintf("slope_%s", sub("\\.", "_", threshold))
colnames(colData(spe_in_tissue))[which(colnames(colData(spe_in_tissue)) == "slope_threshold")] <- qc_colname
plot_qc <- plotQC(spe_in_tissue, type="spots", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", discard=qc_colname)
plot_qc$layers[[1]]$aes_params$size <- 2
print(plot_qc)
                                      
















