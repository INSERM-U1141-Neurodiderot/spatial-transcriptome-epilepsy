.libPaths(c("/home/christophe.lepriol/NeuroDev_ADD/R/r_4.1.0", .libPaths()))
library(scater)
library(ggplot2)
library(cowplot)
library(patchwork) # for wrap_plots() function
#library(Seurat)


##############
# Parameters #
##############
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"
src_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/src", work_dir)
genome_name <- "NCBIRefSeq108_NCBIRefSeq108GTF"
st_pipeline_dir <- sprintf("%s/10-ST_analysis/00-ST_Pipeline/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples", work_dir, genome_name)
space_ranger_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples", work_dir, genome_name)
samples <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")
#sample <- "A_L1_S1"
#samples <- c("A_L1_S1", "A_L2_S5")


#############
# Functions #
#############
source(sprintf("%s/20-pipelines-Seurat-functions.R", src_dir))

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

spot_stat_plots <- function(sp_exp, colname, plot_title, plot_x_lab, thres_vec, gt_bool, plot_list) {
    # distribution
    spot_data <- colData(sp_exp)[[colname]]
    percentiles <- plot_density(spot_data, plot_title, plot_x_lab)    
    # spatial plots
    print(plotSpots(sp_exp, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate=colname, palette=c("white", "black"), size=2))
    print(plotVisium(sp_exp, fill=colname, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    # selecting thresholds and plot discarded spots
    count_vector <- c()
    for (threshold in thres_vec) {
        qc_colname <- sprintf("qc_%s_%d", colname, threshold)
        sp_exp <- qc_threshold(sp_exp, colname, threshold, gt_bool, qc_colname)
        count_vector <- c(count_vector, unname(table(colData(sp_exp)[[qc_colname]])["TRUE"]))
        plot_list[[qc_colname]] <- qc_threshold_plot(sp_exp, qc_colname)
    }
    
    return(list(sp_data=sp_exp, stat=spot_data, counts=count_vector, plots=plot_list))
}

group_density_plots <- function(data_df, data_col, group, v_thres, plot_title, plot_x_lab) {
    density_plot <- ggplot(data_df, aes(x=.data[[data_col]], color=.data[[group]])) +
        geom_line(stat="density") +
        labs(title=plot_title, x=plot_x_lab, y="Density", color="Sample") +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50"))
    if (! is.na(v_thres)) {
        density_plot <- density_plot + geom_vline(xintercept=v_thres, alpha=0.5, linetype="longdash")
    }
    print(density_plot)
}

qc_statistic_thresholds_discarded_spots <- function(data_df, x_col, stat_col, fill_col, stat_name, x_col_name, fill_col_name, thres_vec, plot_list) {
    df2ggplot <- data_df[which(data_df$metric==stat_col),]
    # number of discarded spot barplot
    if (! is.na(fill_col)) {
        barplot <- ggplot(df2ggplot, aes(x=.data[[x_col]], y=count, fill=.data[[fill_col]])) +
            geom_bar(position="dodge", stat="identity") +
            labs(title=sprintf("Number of discarded spots: %s", stat_name), x=x_col_name, y="# spots", fill=fill_col_name)
    } else {
        barplot <- ggplot(df2ggplot, aes(x=.data[[x_col]], y=count)) +
            geom_bar(position="dodge", stat="identity") +
            labs(title=sprintf("Number of discarded spots: %s", stat_name), x=x_col_name, y="# spots")
    }
    barplot <- barplot +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50")) +
        theme(axis.text.x=element_text(angle=30, hjust=1))
    print(barplot)
    # all samples QC plots on a single page
    if (! is.na(thres_vec)) {
        for (threshold in thres_vec) {
            if (threshold - floor(threshold) == 0) {
                qc_name <- sprintf("qc_%s_%d", stat_col, threshold)
                qc_discarded_spots_plot(plot_list, qc_name, sprintf("%s threshold: %d", stat_name, threshold))
            } else {
                qc_name <- sprintf("qc_%s_%s", stat_col, sub("\\.", "_", threshold))
                qc_discarded_spots_plot(plot_list, qc_name, sprintf("%s threshold: %.2f", stat_name, threshold))
            }
        }
    } else {
        qc_name <- sprintf("qc_%s", stat_col)
        qc_discarded_spots_plot(plot_list, qc_name, sprintf("%s", stat_name))
    }
}

qc_discarded_spots_plot <- function(plot_list, qc_name, plot_title) {
    qc_plot_list <- list()
    for (sample in names(plot_list)) {
        sample_plot <- plot_list[[sample]][[qc_name]]
        ## modify plots to adjust with multiple plots on a single page: decrease point size, remove legend, change title
        sample_plot$layers[[1]]$aes_params$size <- 0.01
        sample_plot <- sample_plot + theme(legend.position="none") + labs(title=sample)
        qc_plot_list[[sample]] <- sample_plot
    }
    print(wrap_plots(qc_plot_list) + plot_annotation(title=plot_title))
}


############
# Analysis #
############

dataset_output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/10-Dataset/01-QC", work_dir, genome_name)
if (! dir.exists(dataset_output_dir)) {
    dir.create(dataset_output_dir, recursive=TRUE, mode="0775")
}

# vectors and list for dataset plots
## library size and number of expressed genes densities
spot_stat_vector <- count_sum_vector <- detected_vector <- mito_percent_vector <- hb_percent_vector <- c()
## discarded spots barplots
sample_discarded_vector <- metric_vector <- threshold_vector <- discarded_spot_count_vector <- c()
## plot QC list
plot_qc_list <- list()
## gene expression statistics
gene_stat_vector <- gene_mean_expression_vector <- gene_detected_spot_vector <- c()
mitochondrial_gene_exp_df <- hb_gene_exp_df <- data.frame()

# identify mitochondrial and hemoglobin genes
library(rtracklayer)
raw_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/data/ref_genome/Rattus_norvegicus"
gtf_file <- sprintf("%s/mRatBN7.2/annotations/NCBI_RefSeq/release_108/GCF_015227675.2_mRatBN7.2_genomic.gtf.gz", raw_dir)
gtf_GRanges <- import(gtf_file, format="gtf")
gtf_GRanges_gene <- gtf_GRanges[which(gtf_GRanges$type=="gene"),]
## mitochondrial genes
mt_genes <- gtf_gene_regex(gtf_GRanges_gene, "gene_id", "mt-", dataset_output_dir)
#mitochondrial_genes <- gtf_gene_regex(gtf_GRanges_gene, "description", "mitochondrial", dataset_output_dir)
#mitochondrial_geneset <- c(mt_genes, mitochondrial_genes)
mitochondrial_geneset <- mt_genes
## hemoglobin genes
#hemoglobin_genes <- gtf_gene_regex(gtf_GRanges_gene, "description", "(hemo|neuro|myo|cyto|hapto| )globin", dataset_output_dir)
hemoglobin_genes <- gtf_gene_regex(gtf_GRanges_gene, "description", "hemoglobin", dataset_output_dir)
### add LOC120093065: https://www.ncbi.nlm.nih.gov/gene/120093065
hemoglobin_geneset <- c(hemoglobin_genes, "LOC120093065")
detach("package:rtracklayer", unload=TRUE)

for (sample in samples) {
    sample_plot_qc_list <- list()
    library(SpatialExperiment)
    library(ggspavis)
    print(sprintf("sample: %s", sample))
    
    output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s/01-QC", work_dir, genome_name, sample)
    if (! dir.exists(output_dir)) {
      dir.create(output_dir, recursive=TRUE, mode="0775")
    }
    space_ranger_output_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs", work_dir, genome_name, sample)

    # Load data
    #library(SpatialExperiment)
    spe <- load_ST_pipeline_data(sample, st_pipeline_dir, space_ranger_dir)
    
    # identify mitochondrial genes
    #is_mito <- grepl("(^MT(-|.))|(^mt(-|.))", rowData(spe)$gene_name)
    mitochondrial_geneset_spe <- gsub("-", ".", mitochondrial_geneset) # '-' were replaced by '.' during SpatialExperiment object creation
    is_mito <- rowData(spe)$gene_name %in% mitochondrial_geneset_spe
    # identify hemoglobin genes
    #is_hb <- grepl("^Hb(a|b|q|s|z)", rowData(spe)$gene_name)
    hemoglobin_geneset_spe <- gsub("-", ".", hemoglobin_geneset) # '-' were replaced by '.' during SpatialExperiment object creation
    is_hb <- rowData(spe)$gene_name %in% hemoglobin_geneset_spe
    
    # spot statistics
    ## calculate per-spot QC metrics and store in colData
    spe <- addPerCellQC(spe, subsets=list(mito=is_mito, hb=is_hb))
    
    pdf(sprintf("%s/%s_spot_statistics.pdf", output_dir, sample))
    print(plotQC(spe, type="scatter", metric_x="sum", metric_y="detected"))
    ## sum of counts
    spot_data_colname <- "sum"
    thresholds <- c(1000, 2000, 3000, 4000, 5000)
    return_list <- spot_stat_plots(spe, spot_data_colname, "Sum of counts per spot", "Sum of counts", thresholds, FALSE, sample_plot_qc_list)
    spe <- return_list$sp_data
    count_sum_vector <- c(count_sum_vector, return_list$stat)
    sample_plot_qc_list <- return_list$plots
    spot_stat_vector <- c(spot_stat_vector, rep(sample, length(return_list$stat)))
    discarded_spot_count_vector <- c(discarded_spot_count_vector, return_list$counts)
    threshold_vector <- c(threshold_vector, thresholds)
    metric_vector <- c(metric_vector, rep(spot_data_colname, length(thresholds)))
    sample_discarded_vector <- c(sample_discarded_vector, rep(sample, length(thresholds)))
    
    ## number of expressed features
    spot_data_colname <- "detected"
    thresholds <- c(500, 1000, 1500, 2000, 2500)
    return_list <- spot_stat_plots(spe, spot_data_colname, "Detected genes per spot", "Detected genes", thresholds, FALSE, sample_plot_qc_list)
    spe <- return_list$sp_data
    detected_vector <- c(detected_vector, return_list$stat)
    sample_plot_qc_list <- return_list$plots
    discarded_spot_count_vector <- c(discarded_spot_count_vector, return_list$counts)
    threshold_vector <- c(threshold_vector, thresholds)
    metric_vector <- c(metric_vector, rep(spot_data_colname, length(thresholds)))
    sample_discarded_vector <- c(sample_discarded_vector, rep(sample, length(thresholds)))
    
    ## percentage of mitochondrial gene counts
    spot_data_colname <- "subsets_mito_percent"
    print(plotQC(spe, type="scatter", metric_x="sum", metric_y=spot_data_colname))
    thresholds <- c(5, 10, 15, 20, 25)
    return_list <- spot_stat_plots(spe, spot_data_colname, "Percentage of mitochondrial gene counts per spot", "Percentage of mitochondrial gene counts", thresholds, TRUE, sample_plot_qc_list)
    spe <- return_list$sp_data
    mito_percent_vector <- c(mito_percent_vector, return_list$stat)
    sample_plot_qc_list <- return_list$plots
    discarded_spot_count_vector <- c(discarded_spot_count_vector, return_list$counts)
    threshold_vector <- c(threshold_vector, thresholds)
    metric_vector <- c(metric_vector, rep(spot_data_colname, length(thresholds)))
    sample_discarded_vector <- c(sample_discarded_vector, rep(sample, length(thresholds)))
    
    ## percentage of hemoglobin gene counts
    spot_data_colname <- "subsets_hb_percent"
    print(plotQC(spe, type="scatter", metric_x="sum", metric_y=spot_data_colname))
    thresholds <- c(3, 5, 10, 15)
    return_list <- spot_stat_plots(spe, spot_data_colname, "Percentage of hemoglobin gene counts per spot", "Percentage of hemoglobin gene counts", thresholds, TRUE, sample_plot_qc_list)
    spe <- return_list$sp_data
    hb_percent_vector <- c(hb_percent_vector, return_list$stat)
    sample_plot_qc_list <- return_list$plots
    discarded_spot_count_vector <- c(discarded_spot_count_vector, return_list$counts)
    threshold_vector <- c(threshold_vector, thresholds)
    metric_vector <- c(metric_vector, rep(spot_data_colname, length(thresholds)))
    sample_discarded_vector <- c(sample_discarded_vector, rep(sample, length(thresholds)))
    
    ## remove low-quality spots
    ### number of discarded spots for each metric
    apply(cbind(colData(spe)$qc_sum_2000, colData(spe)$qc_detected_1000), 2, sum)
    ### combined set of discarded spots
    discard <- colData(spe)$qc_sum_2000 | colData(spe)$qc_detected_1000
    table(discard)
    discarded_spot_count_vector <- c(discarded_spot_count_vector, unname(table(discard)["TRUE"]))
    threshold_vector <- c(threshold_vector, NA)
    metric_vector <- c(metric_vector, "discard")
    sample_discarded_vector <- c(sample_discarded_vector, sample)
    ### store in object
    colData(spe)$discard <- discard
    plot_qc <- plotQC(spe, type="spots", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", discard="discard")
    plot_qc$layers[[1]]$aes_params$size <- 2
    print(plot_qc)
    sample_plot_qc_list[["qc_discard"]] <- plot_qc
    ### remove combined set of low-quality spots
    spe <- spe[, !colData(spe)$discard]
    dim(spe)
    print(plotSpots(spe, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="sum", palette=c("white", "black"), size=2))
    print(plotSpots(spe, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="detected", palette=c("white", "black"), size=2))
    
    
    # filter spots with high values of sum of read counts without an increase of number of expressed genes
    ## reproduce loess regression plotQC() function
    df <- data.frame(barcode=rownames(colData(spe)), sum=colData(spe)$sum, detected=colData(spe)$detected)
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
    ### add fitted and slope values to colData(spe)
    colData(spe)$fitted <- loess_fit_slope_df[colData(spe)$barcode, "fitted"]
    colData(spe)$slope <- loess_fit_slope_df[colData(spe)$barcode, "slope"]
    
    thresholds <- c(0, 0.01, 0.02)
    for (threshold in thresholds) {
        ### loess curve slope only
        qc_colname <- sprintf("qc_loess_slope_%s", sub("\\.", "_", threshold))
        loess_fit_slope_df <- cbind(loess_fit_slope_df, slope_threshold=loess_fit_slope_df$slope < threshold)
        discarded_spot_count_vector <- c(discarded_spot_count_vector, unname(table(loess_fit_slope_df$slope_threshold)["TRUE"]))
        colnames(loess_fit_slope_df)[which(colnames(loess_fit_slope_df) == "slope_threshold")] <- qc_colname
        colData(spe)$slope_threshold <- loess_fit_slope_df[colData(spe)$barcode, qc_colname]
        colnames(colData(spe))[which(colnames(colData(spe)) == "slope_threshold")] <- qc_colname
        plot_qc <- plotQC(spe, type="spots", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", discard=qc_colname)
        plot_qc$layers[[1]]$aes_params$size <- 2
        print(plot_qc)
        sample_plot_qc_list[[qc_colname]] <- plot_qc
        ### filter spots based on loess curve slope and if detected < fitted
        qc_colname <- sprintf("qc_loess_slope_fitted_%s", sub("\\.", "_", threshold))
        loess_fit_slope_df <- cbind(loess_fit_slope_df, slope_threshold_fitted=loess_fit_slope_df$slope < threshold & loess_fit_slope_df$detected < loess_fit_slope_df$fitted)
        discarded_spot_count_vector <- c(discarded_spot_count_vector, unname(table(loess_fit_slope_df$slope_threshold_fitted)["TRUE"]))
        colnames(loess_fit_slope_df)[which(colnames(loess_fit_slope_df) == "slope_threshold_fitted")] <- qc_colname
        colData(spe)$slope_threshold_fitted <- loess_fit_slope_df[colData(spe)$barcode, qc_colname]
        colnames(colData(spe))[which(colnames(colData(spe)) == "slope_threshold_fitted")] <- qc_colname
        plot_qc <- plotQC(spe, type="spots", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", discard=qc_colname)
        plot_qc$layers[[1]]$aes_params$size <- 2
        print(plot_qc)
        sample_plot_qc_list[[qc_colname]] <- plot_qc
        threshold_vector <- c(threshold_vector, rep(threshold, 2))
    }
    plot_qc_list[[sample]] <- sample_plot_qc_list
    metric_vector <- c(metric_vector, rep(c("loess_slope", "loess_slope_fitted"), length(thresholds)))
    sample_discarded_vector <- c(sample_discarded_vector, rep(sample, 2*length(thresholds)))
    dev.off()
    write.csv(colData(spe), file=sprintf("%s/%s_spot_statistics.csv", output_dir, sample), quote=FALSE, row.names=FALSE)
    
    # gene statistics and plots
    nb_genes <- dim(spe@assays@data$counts)[1]
    nb_spots <- dim(spe@assays@data$counts)[2]
    ## calculate per-gene metrics and store in rowData
    spe <- addPerFeatureQC(spe)
    rowData(spe)$mean_rank <- rank(-rowData(spe)$mean)
    
    # create Seurat object
    library(Seurat)
    seurat_obj <- CreateSeuratObject(counts=counts(spe))
    ## SpatialExperiment package is required for Seurat object creation, if it is loaded before creating Seurat object, then it is automatically loaded during Seurat object creation
    detach("package:ggspavis", unload=TRUE) # also unload ggspavis package which imports ‘SpatialExperiment’ namespace
    detach("package:SpatialExperiment", unload=TRUE)
    ## add image
    img <- Read10X_Image(image.dir=file.path(space_ranger_output_dir, "spatial"))
    Key(img) <- "image_" # add key for SpatialFeaturePlot() and SpatialPlot() functions
    seurat_obj@images <- list(sample=img)
    ### remove discarded barcodes from image coordinates
    discarded_barcodes <- rownames(seurat_obj@images$sample@coordinates)[! rownames(seurat_obj@images$sample@coordinates) %in% rownames(colData(spe))]
    seurat_obj@images$sample@coordinates <- seurat_obj@images$sample@coordinates[! rownames(seurat_obj@images$sample@coordinates) %in% discarded_barcodes,]
    #spe_seurat@images$sample@assay <- "RNA"
    active_assay <- seurat_obj@active.assay
    seurat_obj@images$sample@assay <- active_assay
    ## add feature data
    ### replace '_' by '-' in feature names: “Feature names cannot have underscores ('_'), replacing with dashes ('-')”
    feature_data_df <- as.data.frame(rowData(spe))
    feature_data_df$gene_name <- gsub("_", "-", feature_data_df$gene_name)
    rownames(feature_data_df) <- feature_data_df$gene_name
    seurat_obj@assays[[active_assay]]@meta.features <- feature_data_df
    #test <- Load10X_Spatial(data.dir=space_ranger_out_dir)
    
    ## mitochondrial genes
    mitochondrial_genes <- mitochondrial_geneset_spe[mitochondrial_geneset_spe %in% rownames(seurat_obj@assays[[active_assay]]@meta.features)]
    #mitochondrial_genes <- rownames(seurat_obj@assays[[active_assay]]@meta.features[is_mito,])
    sample_mitochondrial_gene_exp_df <- data.frame()
    pdf(sprintf("%s/%s_mitochondrial_genes.pdf", output_dir, sample))
    for (one_gene in mitochondrial_genes) {
        one_gene_mean_exp_percentiles <- gene_spatial_mean_exp_plots(seurat_obj, one_gene, active_assay)
        sample_mitochondrial_gene_exp_df <- rbind(sample_mitochondrial_gene_exp_df, c(one_gene, one_gene_mean_exp_percentiles))
    }
    dev.off()
    colnames(sample_mitochondrial_gene_exp_df) <- c("gene", names(one_gene_mean_exp_percentiles))
    write.csv(sample_mitochondrial_gene_exp_df, file=sprintf("%s/%s_mitochondrial_genes.csv", output_dir, sample), quote=FALSE, row.names=FALSE)
    mitochondrial_gene_exp_df <- rbind(mitochondrial_gene_exp_df, cbind(sample=rep(sample, dim(sample_mitochondrial_gene_exp_df)[1]), sample_mitochondrial_gene_exp_df))
    
    ## hemoglobin genes
    hb_genes <- hemoglobin_geneset_spe[hemoglobin_geneset_spe %in% rownames(seurat_obj@assays[[active_assay]]@meta.features)]
    #hb_genes <- rownames(seurat_obj@assays[[active_assay]]@meta.features[is_hb,])
    sample_hb_gene_exp_df <- data.frame()
    pdf(sprintf("%s/%s_hemoglobin_genes.pdf", output_dir, sample))
    for (one_gene in hb_genes) {
        one_gene_mean_exp_percentiles <- gene_spatial_mean_exp_plots(seurat_obj, one_gene, active_assay)
        sample_hb_gene_exp_df <- rbind(sample_hb_gene_exp_df, c(one_gene, one_gene_mean_exp_percentiles))
    }
    dev.off()
    colnames(sample_hb_gene_exp_df) <- c("gene", names(one_gene_mean_exp_percentiles))
    write.csv(sample_hb_gene_exp_df, file=sprintf("%s/%s_hemoglobin_genes.csv", output_dir, sample), quote=FALSE, row.names=FALSE)
    hb_gene_exp_df <- rbind(hb_gene_exp_df, cbind(sample=rep(sample, dim(sample_hb_gene_exp_df)[1]), sample_hb_gene_exp_df))
    
    ## gene percentage of spots with detected expression
    pdf(sprintf("%s/%s_detected_genes.pdf", output_dir, sample))
    percentiles <- plot_density(seurat_obj@assays[[active_assay]]@meta.features$detected, sprintf("Percentage of spots with detected expression\ntotal genes: %d, total spots: %d", nb_genes, nb_spots), "Percentage of spots")
    gene_detected_spot_vector <- c(gene_detected_spot_vector, seurat_obj@assays[[active_assay]]@meta.features$detected)
    gene_stat_vector <- c(gene_stat_vector, rep(sample, length(seurat_obj@assays[[active_assay]]@meta.features$detected)))
    for (one_pct in seq(10, 100, 10)) {
        row_data_df <- seurat_obj@assays[[active_assay]]@meta.features[which(seurat_obj@assays[[active_assay]]@meta.features$detected >= one_pct),]
        row_data_df <- row_data_df[order(row_data_df$detected),]
        one_gene <- rownames(row_data_df[1,])
        one_gene_mean_exp_percentiles <- gene_spatial_mean_exp_plots(seurat_obj, one_gene, active_assay)
    }
    dev.off()
    
    ## gene mean expression
    ### quantiles
    pdf(sprintf("%s/%s_mean_exp_genes_quantiles.pdf", output_dir, sample))
    percentiles <- plot_density(log(seurat_obj@assays[[active_assay]]@meta.features$mean, 2), sprintf("Gene mean expression\ntotal genes: %d, total spots: %d", nb_genes, nb_spots), "log2(Mean expression)")
    gene_mean_expression_vector <- c(gene_mean_expression_vector, log(seurat_obj@assays[[active_assay]]@meta.features$mean, 2))
    mean_exp_quantiles <- quantile(seurat_obj@assays[[active_assay]]@meta.features$mean, probs=c(seq(0, 0.95, 0.05), 0.99, 1))
    quantile_gene_mean_exp_df <- data.frame()
    for(i in 2:length(mean_exp_quantiles)) {
        row_data_df <- seurat_obj@assays[[active_assay]]@meta.features[which(seurat_obj@assays[[active_assay]]@meta.features$mean == mean_exp_quantiles[i]),]
        if (dim(row_data_df)[1] == 0) {
            row_data_df <- seurat_obj@assays[[active_assay]]@meta.features[which(seurat_obj@assays[[active_assay]]@meta.features$mean >= mean_exp_quantiles[i]),]
            row_data_df <- row_data_df[order(row_data_df$mean),]
            one_gene <- rownames(row_data_df[1,])
        } else {
            one_gene <- rownames(row_data_df[1,])
        }
        one_gene_mean_exp_percentiles <- gene_spatial_mean_exp_plots(seurat_obj, one_gene, active_assay)
        quantile_gene_mean_exp_df <- rbind(quantile_gene_mean_exp_df, c(one_gene, one_gene_mean_exp_percentiles))
    }
    dev.off()
    colnames(quantile_gene_mean_exp_df) <- c("gene", names(one_gene_mean_exp_percentiles))
    write.csv(quantile_gene_mean_exp_df, file=sprintf("%s/%s_mean_exp_genes_quantiles.csv", output_dir, sample), quote=FALSE, row.names=FALSE)
    
    ### top 20
    pdf(sprintf("%s/%s_mean_exp_genes_top20.pdf", output_dir, sample))
    percentiles <- plot_density(log(seurat_obj@assays[[active_assay]]@meta.features$mean, 2), sprintf("Gene mean expression\ntotal genes: %d, total spots: %d", nb_genes, nb_spots), "log2(Mean expression)")
    row_data_df <- seurat_obj@assays[[active_assay]]@meta.features
    row_data_df <- row_data_df[order(row_data_df$mean_rank),]
    top20_gene_mean_exp_df <- data.frame()
    for(one_gene in rownames(row_data_df)[1:20]) {
        one_gene_mean_exp_percentiles <- gene_spatial_mean_exp_plots(seurat_obj, one_gene, active_assay)
        top20_gene_mean_exp_df <- rbind(top20_gene_mean_exp_df, c(one_gene, one_gene_mean_exp_percentiles))
    }
    dev.off()
    colnames(top20_gene_mean_exp_df) <- c("gene", names(one_gene_mean_exp_percentiles))
    write.csv(top20_gene_mean_exp_df, file=sprintf("%s/%s_mean_exp_genes_top20.csv", output_dir, sample), quote=FALSE, row.names=FALSE)
    
    
    
    detach("package:Seurat", unload=TRUE)
    detach("package:SeuratObject", unload=TRUE)
}

# dataset plots
## sum of counts, number of expressed features, percentage of mitochondrial gene counts, and percentage of hemoglobin gene counts per spot distributions
df2ggplot <- data.frame(sample=spot_stat_vector, count_sum=count_sum_vector, detected=detected_vector, mito_pct=mito_percent_vector, hb_pct=hb_percent_vector)
write.csv(df2ggplot, file=sprintf("%s/spot_statistic_distributions.csv", dataset_output_dir), quote=FALSE, row.names=FALSE)
pdf(sprintf("%s/spot_statistic_distributions.pdf", dataset_output_dir))
group_density_plots(df2ggplot, "count_sum", "sample", NA, "Sum of counts per spot", "Sum of counts")
group_density_plots(df2ggplot, "detected", "sample", NA, "Detected genes per spot", "Detected genes")
group_density_plots(df2ggplot, "mito_pct", "sample", NA, "Percentage of mitochondrial gene counts per spot", "Percentage of mitochondrial gene counts")
group_density_plots(df2ggplot, "hb_pct", "sample", NA, "Percentage of hemoglobin gene counts per spot", "Percentage of hemoglobin gene counts")
dev.off()

## discarded spots: sum of counts, number of expressed features, percentage of mitochondrial gene counts, and percentage of hemoglobin gene counts per spot
df2ggplot <- data.frame(sample=sample_discarded_vector, metric=metric_vector, threshold=factor(threshold_vector), count=discarded_spot_count_vector)
write.csv(df2ggplot, file=sprintf("%s/spot_statistic_discarded_spots.csv", dataset_output_dir), quote=FALSE, row.names=FALSE)
### sum of counts
qc_statistic <- "sum"
thresholds <- c(1000, 2000, 3000, 4000, 5000)
pdf(sprintf("%s/spot_statistic_discarded_spots_%s.pdf", dataset_output_dir, qc_statistic))
qc_statistic_thresholds_discarded_spots(df2ggplot, "sample", qc_statistic, "threshold", "Sum of counts", "Sample", "Threshold", thresholds, plot_qc_list)
dev.off()
### number of expressed features
qc_statistic <- "detected"
thresholds <- c(500, 1000, 1500, 2000, 2500)
pdf(sprintf("%s/spot_statistic_discarded_spots_%s.pdf", dataset_output_dir, qc_statistic))
qc_statistic_thresholds_discarded_spots(df2ggplot, "sample", qc_statistic, "threshold", "Detected genes", "Sample", "Threshold", thresholds, plot_qc_list)
dev.off()
### percentage of mitochondrial gene counts
qc_statistic <- "subsets_mito_percent"
thresholds <- c(5, 10, 15, 20, 25)
pdf(sprintf("%s/spot_statistic_discarded_spots_%s.pdf", dataset_output_dir, qc_statistic))
qc_statistic_thresholds_discarded_spots(df2ggplot, "sample", qc_statistic, "threshold", "Percentage of mitochondrial gene counts", "Sample", "Threshold", thresholds, plot_qc_list)
dev.off()
### percentage of hemoglobin gene counts
qc_statistic <- "subsets_hb_percent"
thresholds <- c(3, 5, 10, 15)
pdf(sprintf("%s/spot_statistic_discarded_spots_%s.pdf", dataset_output_dir, qc_statistic))
qc_statistic_thresholds_discarded_spots(df2ggplot, "sample", qc_statistic, "threshold", "Percentage of hemoglobin gene counts", "Sample", "Threshold", thresholds, plot_qc_list)
dev.off()
### all qc statistic discarded spots
qc_statistic <- "discard"
pdf(sprintf("%s/spot_statistic_discarded_spots_%s.pdf", dataset_output_dir, qc_statistic))
qc_statistic_thresholds_discarded_spots(df2ggplot, "sample", qc_statistic, NA, "All qc statistic discarded spots", "Sample", NA, NA, plot_qc_list)
dev.off()

### high values of sum of read counts without an increase of number of expressed genes
#### number of discarded spots
##### loess curve slope only
qc_statistic <- "loess_slope"
thresholds <- c(0, 0.01, 0.02)
pdf(sprintf("%s/spot_statistic_discarded_spots_%s.pdf", dataset_output_dir, qc_statistic))
qc_statistic_thresholds_discarded_spots(df2ggplot, "sample", qc_statistic, "threshold", "Loess slope", "Sample", "Threshold", thresholds, plot_qc_list)
dev.off()
##### filter spots based on loess curve slope and if detected < fitted
qc_statistic <- "loess_slope_fitted"
thresholds <- c(0, 0.01, 0.02)
pdf(sprintf("%s/spot_statistic_discarded_spots_%s.pdf", dataset_output_dir, qc_statistic))
qc_statistic_thresholds_discarded_spots(df2ggplot, "sample", qc_statistic, "threshold", "Loess slope", "Sample", "Threshold", thresholds, plot_qc_list)
dev.off()


## library size and number of expressed genes per spot distributions
pdf(sprintf("%s/gene_expression_statistics.pdf", dataset_output_dir))
df2ggplot <- data.frame(sample=gene_stat_vector, mean=gene_mean_expression_vector, detected=gene_detected_spot_vector)
#df2ggplot <- data.frame(sample=gene_stat_vector, mean=gene_mean_expression_vector)
write.csv(df2ggplot, file=sprintf("%s/gene_expression_statistics.csv", dataset_output_dir), quote=FALSE, row.names=FALSE)
density_plot <- ggplot(df2ggplot, aes(x=mean, color=sample)) +
    geom_line(stat="density") +
    labs(title="Gene mean expression", x="Mean expression", y="Density", color="Sample") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50"))
print(density_plot)
density_plot <- ggplot(df2ggplot, aes(x=detected, color=sample)) +
    geom_line(stat="density") +
    labs(title="Gene percentage of spots with detected expression", x="Percentage of spots", y="Density", color="Sample") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50"))
print(density_plot)
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
                                      
















