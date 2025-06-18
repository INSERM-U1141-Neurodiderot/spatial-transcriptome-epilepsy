.libPaths(c("/home/christophe.lepriol/NeuroDev_ADD/R/r_4.1.0", .libPaths()))
library(SpatialExperiment)
library(ggspavis)
library(Seurat)
library(scater)
library(patchwork)
library(scran)


##############
# Parameters #
##############
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"
genome_name <- "NCBIRefSeq108_NCBIRefSeq108GTF"
samples <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")
#samples <- c("A_L1_S1", "A_L2_S5")
sample <- "A_L1_S1"

for (sample in samples) {
    print(sprintf("sample: %s", sample))
    output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s/10-normalization", work_dir, genome_name, sample)
    if (! dir.exists(output_dir)) {
        dir.create(output_dir, recursive=TRUE, mode="0775")
    }

    # load data
    ## expression matrices
    ### Space Ranger
    #h5_file <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs/raw_feature_bc_matrix.h5", work_dir, genome_name, sample)
    #h5_data <- Read10X_h5(h5_file)
    ### ST Pipeline
    st_pipeline_matrix_file <- sprintf("%s/10-ST_analysis/00-ST_Pipeline/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s/10-Pipeline/%s_stdata.tsv", work_dir, genome_name, sample, sample)
    st_pipeline_matrix <- read.table(st_pipeline_matrix_file, sep="\t", header=TRUE, quote="", row.names=1)

    ## Space Ranger output directory for tissue positions list and image
    space_ranger_output_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs", work_dir, genome_name, sample)

    # build SpatialExperiment object
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

    spe <- SpatialExperiment(assay=list(counts=t(st_pipeline_matrix)), rowData=gene_data, colData=spatial_coordinates_df[rownames(st_pipeline_matrix),], imgData=img, spatialCoordsNames=c("pxl_col_in_fullres", "pxl_row_in_fullres"), sample_id=sample)

    # QC
    ## keep spot over tissue
    spe_in_tissue <- spe[, colData(spe)$in_tissue==1]
    ## calculate per-spot QC metrics and store in colData
    ### raw counts
    assay_type <- "counts"
    spe_in_tissue <- addPerCellQC(spe_in_tissue, assay.type=assay_type)
    for (one_stat in c("sum", "detected", "total")) {
        colnames(colData(spe_in_tissue))[colnames(colData(spe_in_tissue)) == one_stat] <- sprintf("%s.%s", assay_type, one_stat)
    }
    head(colData(spe_in_tissue))

    # normalization
    ## methods from scater and scran
    ### quick clustering for pool-based size factors
    set.seed(123)
    qclus <- quickCluster(spe_in_tissue)
    table(qclus)
    ### calculate size factors and store in object
    spe_in_tissue <- computeSumFactors(spe_in_tissue, cluster=qclus)
    summary(sizeFactors(spe_in_tissue))
    hist(sizeFactors(spe_in_tissue), breaks=20)
    ### calculate logcounts (log-transformed normalized counts) and store in object
    spe_in_tissue <- logNormCounts(spe_in_tissue)
    ### check
    assayNames(spe_in_tissue)
    dim(counts(spe_in_tissue))
    dim(logcounts(spe_in_tissue))
    
    ### calculate per-spot QC metrics and store in colData
    #### normalized counts
    assay_type <- "logcounts"
    spe_in_tissue <- addPerCellQC(spe_in_tissue, assay.type=assay_type)
    for (one_stat in c("sum", "detected", "total")) {
        colnames(colData(spe_in_tissue))[colnames(colData(spe_in_tissue)) == one_stat] <- sprintf("%s.%s", assay_type, one_stat)
    }
    head(colData(spe_in_tissue))

    # plots
    pdf(sprintf("%s/raw_normalized_counts.pdf", output_dir))
    ## total counts and number of detected genes
    hist(colData(spe_in_tissue)$counts.sum, breaks=20, main="Sum of counts per spot", xlab="Sum of raw counts")
    hist(colData(spe_in_tissue)$counts.detected, breaks=20, main="Detected genes per spot", xlab="Detected genes - raw counts")
    hist(colData(spe_in_tissue)$logcounts.sum, breaks=20, main="Sum of counts per spot", xlab="Sum of normalized counts")
    hist(colData(spe_in_tissue)$logcounts.detected, breaks=20, main="Detected genes per spot", xlab="Detected genes - normalized counts")
    ## splot plots
    ### raw counts
    print(plotSpots(spe_in_tissue, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="counts.sum", palette=c("white", "black")))
    print(plotVisium(spe_in_tissue, fill="counts.sum", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    print(plotSpots(spe_in_tissue, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="counts.detected", palette=c("white", "black")))
    print(plotVisium(spe_in_tissue, fill="counts.detected", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    ### normalized counts
    print(plotSpots(spe_in_tissue, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="logcounts.sum", palette=c("white", "black")))
    print(plotVisium(spe_in_tissue, fill="logcounts.sum", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    print(plotSpots(spe_in_tissue, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="logcounts.detected", palette=c("white", "black")))
    print(plotVisium(spe_in_tissue, fill="logcounts.detected", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    ### raw and normalized counts
    p1 <- plotVisium(spe_in_tissue, fill="counts.sum", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    p2 <- plotVisium(spe_in_tissue, fill="logcounts.sum", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    print(wrap_plots(p1, p2) + plot_annotation(title="Sum of counts"))
    p1 <- plotVisium(spe_in_tissue, fill="counts.detected", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    p2 <- plotVisium(spe_in_tissue, fill="logcounts.detected", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    print(wrap_plots(p1, p2) + plot_annotation(title="Detected genes"))
    #dev.off()
    
    ## create Seurat object
    seurat_in_tissue <- CreateSeuratObject(counts=counts(spe_in_tissue))
    seurat_in_tissue <- SCTransform(seurat_in_tissue, assay = "RNA")
    ## conversion to SingleCellExperiment object
    #sce_in_tissue <- as.SingleCellExperiment(seurat_in_tissue, assay=c("RNA", "SCT"))
    ## --> same as seurat_in_tissue@assays$SCT@counts
    ## create new SpatialExperiment object with normalized counts to SingleCellExperiment object
    #assays(spe_in_tissue, i="SCT") <- counts(sce_in_tissue) # very long, stop before finished
    #assay(spe_in_tissue, i="SCT", counts(sce_in_tissue)) # not tested
    #assay(spe_in_tissue, i="SCT") <- counts(sce_in_tissue) # very long, stop before finished
    #assay(spe_in_tissue, "SCT") <- assay(sce_in_tissue, i="counts") # very long, stop before finished: max RAM reached ? 9:04
    #norm_counts <- assay(sce_in_tissue, i="counts")
    norm_counts <- seurat_in_tissue@assays$SCT@counts
    ## '_' in row names was replaced by '-', example: 'X--ambiguous.L3mbtl3.LOC120097350.' in norm_counts and 'X__ambiguous.L3mbtl3.LOC120097350.' in assay(spe_in_tissue, "counts")
    rownames(norm_counts)[grep("-", rownames(norm_counts))] <- gsub("-", "_", rownames(norm_counts)[grep("-", rownames(norm_counts))])
    ## remove genes not in normalized counts
    spe_in_tissue_norm <- spe_in_tissue[rownames(norm_counts),]
    assay(spe_in_tissue_norm, i="SCT") <- as.matrix(norm_counts)


    # identify mitochondrial genes
    #is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe_in_tissue_norm)$gene_name)
    #table(is_mito)
    #rowData(spe)$gene_name[is_mito]

    ### calculate per-spot QC metrics and store in colData
    #### normalized counts
    assay_type <- "SCT"
    spe_in_tissue_norm <- addPerCellQC(spe_in_tissue_norm, assay.type=assay_type)
    for (one_stat in c("sum", "detected", "total")) {
        colnames(colData(spe_in_tissue_norm))[colnames(colData(spe_in_tissue_norm)) == one_stat] <- sprintf("%s.%s", assay_type, one_stat)
    }
    head(colData(spe_in_tissue_norm))

    # plots
    #pdf(sprintf("%s/raw_normalized_counts.pdf", output_dir))
    ## total counts and number of detected genes
    hist(colData(spe_in_tissue_norm)$counts.sum, breaks=20, main="Sum of counts per spot", xlab="Sum of raw counts")
    hist(colData(spe_in_tissue_norm)$counts.detected, breaks=20, main="Detected genes per spot", xlab="Detected genes - raw counts")
    hist(colData(spe_in_tissue_norm)$SCT.sum, breaks=20, main="Sum of counts per spot", xlab="Sum of normalized counts")
    hist(colData(spe_in_tissue_norm)$SCT.detected, breaks=20, main="Detected genes per spot", xlab="Detected genes - normalized counts")
    ## splot plots
    ### raw counts
    print(plotSpots(spe_in_tissue_norm, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="counts.sum", palette=c("white", "black")))
    print(plotVisium(spe_in_tissue_norm, fill="counts.sum", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    print(plotSpots(spe_in_tissue_norm, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="counts.detected", palette=c("white", "black")))
    print(plotVisium(spe_in_tissue_norm, fill="counts.detected", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    ### normalized counts
    print(plotSpots(spe_in_tissue_norm, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="SCT.sum", palette=c("white", "black")))
    print(plotVisium(spe_in_tissue_norm, fill="SCT.sum", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    print(plotSpots(spe_in_tissue_norm, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="SCT.detected", palette=c("white", "black")))
    print(plotVisium(spe_in_tissue_norm, fill="SCT.detected", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    ### raw and normalized counts
    p1 <- plotVisium(spe_in_tissue_norm, fill="counts.sum", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    p2 <- plotVisium(spe_in_tissue_norm, fill="SCT.sum", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    print(wrap_plots(p1, p2) + plot_annotation(title="Sum of counts"))
    p1 <- plotVisium(spe_in_tissue_norm, fill="counts.detected", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    p2 <- plotVisium(spe_in_tissue_norm, fill="SCT.detected", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    print(wrap_plots(p1, p2) + plot_annotation(title="Detected genes"))
    ### logcounts and SCT counts
    p1 <- plotVisium(spe_in_tissue_norm, fill="logcounts.sum", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    p2 <- plotVisium(spe_in_tissue_norm, fill="SCT.sum", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    print(wrap_plots(p1, p2) + plot_annotation(title="Sum of counts"))
    p1 <- plotVisium(spe_in_tissue_norm, fill="logcounts.detected", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    p2 <- plotVisium(spe_in_tissue_norm, fill="SCT.detected", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    print(wrap_plots(p1, p2) + plot_annotation(title="Detected genes"))
    dev.off()
}

