.libPaths(c("/home/christophe.lepriol/NeuroDev_ADD/R/r_4.1.0", .libPaths()))
library(SpatialExperiment)
library(ggspavis)
library(Seurat)
library(scater)
library(patchwork)
library(scran) # for buildSNNGraph() function


##############
# Parameters #
##############
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"
genome_name <- "NCBIRefSeq108_NCBIRefSeq108GTF"
#samples <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")
sample <- "A_L1_S1"

print(sprintf("sample: %s", sample))
#output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s/10-normalization", work_dir, genome_name, sample)
#if (! dir.exists(output_dir)) {
#    dir.create(output_dir, recursive=TRUE, mode="0775")
#}

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
## matrix: convert coordinates to barcodes
barcodes <- unlist(lapply(rownames(st_pipeline_matrix), function(x, coord2barcode=spatial_coordinates_df) {
    coordinates <- unlist(strsplit(x, "x"))
    return(coord2barcode[which(coord2barcode$array_col==as.integer(coordinates[1])-1 & coord2barcode$array_row==as.integer(coordinates[2])-1), "barcode"])
}))
rownames(st_pipeline_matrix) <- barcodes
gene_data <- data.frame(gene_name=colnames(st_pipeline_matrix))
## image data
img <- readImgData(imageSources=file.path(space_ranger_output_dir, "spatial", "tissue_lowres_image.png"), scaleFactors=file.path(space_ranger_output_dir, "spatial", "scalefactors_json.json"), sample_id=sample)
spe <- SpatialExperiment(assay=list(counts=t(st_pipeline_matrix)), rowData=gene_data, colData=spatial_coordinates_df[rownames(st_pipeline_matrix),], imgData=img, spatialCoordsNames=c("pxl_col_in_fullres", "pxl_row_in_fullres"), sample_id=sample)

# keep spot over tissue
spe_in_tissue <- spe[, colData(spe)$in_tissue==1]
# normalization
## create Seurat object
seurat_in_tissue <- CreateSeuratObject(counts=counts(spe_in_tissue))
seurat_in_tissue <- SCTransform(seurat_in_tissue, assay = "RNA")
norm_counts <- seurat_in_tissue@assays$SCT@counts
## '_' in row names was replaced by '-', example: 'X--ambiguous.L3mbtl3.LOC120097350.' in norm_counts and 'X__ambiguous.L3mbtl3.LOC120097350.' in assay(spe_in_tissue, "counts")
rownames(norm_counts)[grep("-", rownames(norm_counts))] <- gsub("-", "_", rownames(norm_counts)[grep("-", rownames(norm_counts))])
## remove genes not in normalized counts
spe_in_tissue_norm <- spe_in_tissue[rownames(norm_counts),]
assay(spe_in_tissue_norm, i="SCT") <- as.matrix(norm_counts)


# feature selection
## find variable genes
#variable_features <- FindVariableFeatures(seurat_in_tissue@assays$SCT)
seurat_in_tissue <- FindVariableFeatures(seurat_in_tissue)
VariableFeaturePlot(seurat_in_tissue, assay="SCT")
variable_genes <- seurat_in_tissue@assays$SCT@var.features

# dimensionality reduction
## apply PCA to the set of top HVGs and retain the top 50 PCs
### '_' in row names was replaced by '-', example: 'X--ambiguous.L3mbtl3.LOC120097350.' in norm_counts and 'X__ambiguous.L3mbtl3.LOC120097350.' in assay(spe_in_tissue, "counts")
variable_genes[grep("-", variable_genes)] <- gsub("-", "_", variable_genes[grep("-", variable_genes)])
set.seed(123)
spe_in_tissue_norm <- runPCA(spe_in_tissue_norm, exprs_values="SCT", subset_row=variable_genes)
reducedDimNames(spe_in_tissue_norm)
dim(reducedDim(spe_in_tissue_norm, "PCA"))

## UMAP on the set of top 50 PCs and retain the top 2 UMAP components, which will be used for visualization purposes
set.seed(123)
spe_in_tissue_norm <- runUMAP(spe_in_tissue_norm, dimred="PCA")
reducedDimNames(spe_in_tissue_norm)
dim(reducedDim(spe_in_tissue_norm, "UMAP"))
colnames(reducedDim(spe_in_tissue_norm, "UMAP")) <- paste0("UMAP", 1:2)

## visualizations
### plot top 2 PCA dimensions
plotDimRed(spe_in_tissue_norm, type="PCA")
### plot top 2 UMAP dimensions
plotDimRed(spe_in_tissue_norm, type="UMAP")


# clustering on HVGs
## graph-based clustering
set.seed(123)
k <- 10
g <- buildSNNGraph(spe_in_tissue_norm, k=k, use.dimred="PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
table(clus)
### store cluster labels in column 'label' in colData
colLabels(spe_in_tissue_norm) <- factor(clus)

## visualizations
### in spatial coordinates on the tissue slide
plotSpots(spe_in_tissue_norm, annotate="label", palette="libd_layer_colors")
plotSpots(spe_in_tissue_norm, annotate="label")
plotSpots(spe_in_tissue_norm)
### in reduced dimension spaces




