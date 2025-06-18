.libPaths(c(.libPaths(), "/home/christophe.lepriol/NeuroDev_ADD/R/r_4.1.0"))
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
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
timepoint <- 20
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
spot_hb_gene_count_pct_threshold <- 100
## normalization, feature selection, dimension reduction, clustering
merge_before_norm <- FALSE
normalization_method_vector <- c("SCTransform", "SCTransform_v2", "LogNormalize")
sct_variables_to_regress <- "qc_mito_percent_RNA"
scale_data_variables_to_regress <- "qc_mito_percent_RNA"
variables_to_regress_in_path <- "mito"
scale_data_split_by <- "orig.ident"
scale_data_split_by_in_path <- "orig"
integration_group <- "orig.ident"
integration_group_by_in_path <- "orig"
cluster_markers <- TRUE
DE_analysis <- TRUE
conserved_markers <- FALSE

# number of samples per page when using SpatialFeaturePlot() function
nb_SpatialFeaturePlot_per_page <- 4
# number of samples per page when using SpatialFeaturePlot() function
nb_cluster_markers_to_plot <- 4

nb_DEGs_to_plot <- 30
nb_DEGs_per_page <- 3

#metadata_categories <- c("lame", "zone", "condition", "orig.ident")

merge_before_norm_in_path <- ifelse(merge_before_norm, "Merge_norm", "Norm_merge")
qc_output_dir <- sprintf("%s/%s/Norm_%s-scale_%s-int_%s/QC-sum%d_det%s_mito%d_hb%d_OLD", output_dir, merge_before_norm_in_path, variables_to_regress_in_path, scale_data_split_by_in_path, integration_group_by_in_path, spot_sum_of_counts_threshold, spot_detected_genes_threshold, spot_mito_gene_count_pct_threshold, spot_hb_gene_count_pct_threshold)
qc_reclustering_marker_output_dir <- sprintf("%s/rec_mark", qc_output_dir)

#############
# Functions #
#############
source(sprintf("%s/20-pipelines-Seurat-functions.R", src_dir))
source(sprintf("%s/6X-integration_functions.R", src_dir))


############
# Analysis #
############

# load data
print("load data")
library(SpatialExperiment)
library(ggspavis)
st_pipeline_dir <- sprintf("%s/10-ST_analysis/00-ST_Pipeline/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples", work_dir, genome_name)
space_ranger_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples", work_dir, genome_name)
spe_objects <- lapply(sections, load_ST_pipeline_data, st_dir=st_pipeline_dir, sr_dir=space_ranger_dir)
names(spe_objects) <- sections

# quality control
library(scater) # Loading required package: scuttle
if (! dir.exists(qc_reclustering_marker_output_dir)) {
    dir.create(qc_reclustering_marker_output_dir, recursive=TRUE, mode="0775")
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
mt_genes <- gtf_gene_regex(gtf_GRanges_gene, "gene_id", "mt-", qc_reclustering_marker_output_dir)
#mitochondrial_genes <- gtf_gene_regex(gtf_GRanges_gene, "description", "mitochondrial", qc_reclustering_marker_output_dir)
#mitochondrial_geneset <- c(mt_genes, mitochondrial_genes)
mitochondrial_geneset <- mt_genes
mitochondrial_geneset_spe <- gsub("-", ".", mitochondrial_geneset) # '-' were replaced by '.' during SpatialExperiment object creation
## hemoglobin genes
#hemoglobin_genes <- gtf_gene_regex(gtf_GRanges_gene, "description", "(hemo|neuro|myo|cyto|hapto| )globin", qc_reclustering_marker_output_dir)
hemoglobin_genes <- gtf_gene_regex(gtf_GRanges_gene, "description", "hemoglobin", qc_reclustering_marker_output_dir)
### add LOC120093065: https://www.ncbi.nlm.nih.gov/gene/120093065
hemoglobin_geneset <- c(hemoglobin_genes, "LOC120093065")
hemoglobin_geneset_spe <- gsub("-", ".", hemoglobin_geneset) # '-' were replaced by '.' during SpatialExperiment object creation
## globin genes
#globin_genes <- gtf_gene_regex(gtf_GRanges_gene, "description", "globin", qc_reclustering_marker_output_dir)
rm(gtf_GRanges)
rm(gtf_GRanges_gene)
gc()
detach("package:rtracklayer", unload=TRUE)

## add QC statistics
qc_stats_output_dir <- sprintf("%s/00-QC_stats", qc_reclustering_marker_output_dir)
if (! dir.exists(qc_stats_output_dir)) {
    dir.create(qc_stats_output_dir, recursive=TRUE, mode="0775")
}
raw_column_prefix <- "raw"
spe_objects <- lapply(spe_objects, add_qc_stats, mito_genes=mitochondrial_geneset_spe, hb_genes=hemoglobin_geneset_spe, col_prefix=raw_column_prefix, out_dir=qc_stats_output_dir, out_name="QC_stats")
lapply(spe_objects, qc_stat_plots, col_prefix=raw_column_prefix, out_dir=qc_stats_output_dir, out_name="QC_stats")

## Seurat feature stat plots and gene spatial plots
### SpatialExperiment to Seurat object conversion
#### detach SpatialExperiment package to avoid conflict with Seurat package
detach("package:ggspavis", unload=TRUE) # also unload ggspavis package which imports ‘SpatialExperiment’ namespace
detach("package:SpatialExperiment", unload=TRUE)
library(Seurat) # also attaching SeuratObject package
seurat_objects <- lapply(spe_objects, SpatialExperiment_Seurat_conversion, sr_dir=space_ranger_dir, add_image=load_image) # Loading required package: SpatialExperiment
#### SpatialExperiment package is required for Seurat object creation, if it is not loaded before creating Seurat object, then it is automatically loaded during Seurat object creation
#### detach SpatialExperiment package to avoid conflict with Seurat package
detach("package:SpatialExperiment", unload=TRUE)
### feature stat plots
seurat_objects <- lapply(seurat_objects, addPercentageFeatureSet2Seurat, mito_genes=mitochondrial_geneset_spe, hb_genes=hemoglobin_geneset_spe, assay2use="RNA", col_prefix=raw_column_prefix)
library(RColorBrewer)
lapply(seurat_objects, seurat_all_feature_stat_plots, assay2use="RNA", col_prefix=raw_column_prefix, plot_title="Raw counts", out_dir=qc_stats_output_dir, out_name="QC_stats") 
detach("package:RColorBrewer")
### mitochondrial and hemoglobin gene spatial plots
#### hemoglobin genes
df_list <- lapply(seurat_objects, qc_spatial_mean_exp_plots, gene_vec=hemoglobin_geneset_spe, mean_order=TRUE, out_dir=qc_stats_output_dir, out_name="hemoglobin_genes")
#### mitochondrial genes
df_list <- lapply(seurat_objects, qc_spatial_mean_exp_plots, gene_vec=mitochondrial_geneset_spe, mean_order=TRUE, out_dir=qc_stats_output_dir, out_name="mitochondrial_genes")
rm(seurat_objects)
gc()
### detach Seurat and SeuratObject packages to avoid conflict with SpatialExperiment package
detach("package:Seurat", unload=TRUE)
detach("package:SeuratObject", unload=TRUE)

## QC filtering
library(SpatialExperiment)
library(ggspavis)
qc_filtering_output_dir <- sprintf("%s/10-QC_filtering", qc_reclustering_marker_output_dir)
if (! dir.exists(qc_filtering_output_dir)) {
    dir.create(qc_filtering_output_dir, recursive=TRUE, mode="0775")
}
qc_column_prefix <- "qc"
spe_objects <- lapply(spe_objects, qc_filter, raw_col_prefix=raw_column_prefix, qc_col_prefix=qc_column_prefix, sum_thres=spot_sum_of_counts_threshold, dectected_thres=spot_detected_genes_threshold, mito_pct_thres=spot_mito_gene_count_pct_threshold, hb_pct_thres=spot_hb_gene_count_pct_threshold, out_dir=qc_filtering_output_dir, out_name="filtered_spots")
### QC statistics after filerting
lapply(spe_objects, function(x, out_dir=qc_filtering_output_dir, out_name="filtered_QC_stats") {
    write.csv(colData(x), file=sprintf("%s/%s_%s.csv", out_dir, x@int_metadata$imgData$sample_id, out_name), quote=FALSE, row.names=TRUE)
})
lapply(spe_objects, qc_stat_plots, col_prefix=raw_column_prefix, out_dir=qc_filtering_output_dir, out_name="filtered_QC_stats")
### feature statistics after filtering
spe_objects <- lapply(spe_objects, addFeatureStats2SpatialExperiment)
### Seurat feature stat plots and gene spatial plots after filtering
#### SpatialExperiment to Seurat object conversion
##### detach SpatialExperiment package to avoid conflict with Seurat package
detach("package:ggspavis", unload=TRUE) # also unload ggspavis package which imports ‘SpatialExperiment’ namespace
detach("package:SpatialExperiment", unload=TRUE)
library(Seurat) # also attaching SeuratObject package
seurat_objects <- lapply(spe_objects, SpatialExperiment_Seurat_conversion, sr_dir=space_ranger_dir, add_image=load_image) # does not load SpatialExperiment package !
#detach("package:SpatialExperiment", unload=TRUE)
rm(spe_objects)
gc()
### feature stat plots
seurat_objects <- lapply(seurat_objects, addPercentageFeatureSet2Seurat, mito_genes=mitochondrial_geneset_spe, hb_genes=hemoglobin_geneset_spe, assay2use="RNA", col_prefix=qc_column_prefix)
library(RColorBrewer)
lapply(seurat_objects, seurat_all_feature_stat_plots, assay2use="RNA", col_prefix=qc_column_prefix, plot_title="Raw counts after QC", out_dir=qc_stats_output_dir, out_name="QC_stats_filtering") 
detach("package:RColorBrewer")
#### mitochondrial and hemoglobin gene spatial plots after filerting
##### hemoglobin genes
df_list <- lapply(seurat_objects, qc_spatial_mean_exp_plots, gene_vec=hemoglobin_geneset_spe, mean_order=TRUE, out_dir=qc_filtering_output_dir, out_name="filtered_hemoglobin_genes")
##### mitochondrial genes
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


library(cluster) # silhouette() function

# merge datasets
if (merge_before_norm) {
    seurat_object_merge <- merge(seurat_objects[[1]], y=unlist(seurat_objects)[2:nb_sections], add.cell.ids=sections)
    rm(seurat_objects)
    gc()
    ## orig.ident, lame, zone and condition are not factors in the merged object meta data
    seurat_object_merge@meta.data$orig.ident <- as.factor(seurat_object_merge@meta.data$orig.ident)
    seurat_object_merge@meta.data$lame <- as.factor(seurat_object_merge@meta.data$lame)
    seurat_object_merge@meta.data$zone <- as.factor(seurat_object_merge@meta.data$zone)
    seurat_object_merge@meta.data$condition <- as.factor(seurat_object_merge@meta.data$condition)
    ## add feature statistics to merged Seurat object
    seurat_object_merge <- addFeatureStats2Seurat(seurat_object_merge, "RNA", "counts")
} else {
    seurat_objects <- lapply(seurat_objects, addFeatureStats2Seurat, assay2use="RNA", slot2use="counts")
}



# get global clustering parameter evaluation output
global_clustering_eval_df <- read.csv(sprintf("%s/10-QC_filtering/%s_norm_featnb_int_dims_res_cluster_mean_silhouette.csv", qc_output_dir, output_dirname))
#global_clustering_eval_subdf <- global_clustering_eval_df[which(global_clustering_eval_df$normalization=="SCTransform_v2" & global_clustering_eval_df$integration %in% c("Seurat", "Harmony")),]
#global_clustering_eval_subdf[order(global_clustering_eval_subdf$mean_silhouette_seurat_clusters, decreasing=TRUE),]


# Seurat
integration_method <- "Seurat"
integration_methods <- c(integration_method)
## get top results for global clustering with Seurat integration
global_clustering_eval_subdf <- global_clustering_eval_df[which(global_clustering_eval_df$normalization=="SCTransform_v2" & global_clustering_eval_df$integration == integration_method),]
global_clustering_eval_subdf_order <- global_clustering_eval_subdf[order(global_clustering_eval_subdf$mean_silhouette_seurat_clusters, decreasing=TRUE),]
global_normalization <- global_clustering_eval_subdf_order[1, "normalization"]
global_features <- global_clustering_eval_subdf_order[1, "features"]
#global_integration <- global_clustering_eval_subdf_order[1, "integration"]
global_dimensions <- global_clustering_eval_subdf_order[1, "dimensions"]
global_resolution <- global_clustering_eval_subdf_order[1, "resolution"]
hippocampus_cluster <- "3"

## get top results for reclustering with Seurat integration
reclustering_eval_df <- read.csv(sprintf("%s/rec/10-QC_filtering/%s_norm_featnb_int_dims_res_reclustering_mean_silhouette.csv", qc_output_dir, output_dirname))
reclustering_eval_subdf <- reclustering_eval_df[which(reclustering_eval_df$normalization=="SCTransform_v2" & reclustering_eval_df$integration == integration_method),]
reclustering_eval_subdf_order <- reclustering_eval_subdf[order(reclustering_eval_subdf$mean_silhouette_seurat_clusters, decreasing=TRUE),]
#re_normalization <- reclustering_eval_subdf_order[1, "normalization"]
#re_features <- reclustering_eval_subdf_order[1, "features"]
#re_integration <- reclustering_eval_subdf_order[1, "integration"]
re_dimensions <- reclustering_eval_subdf_order[1, "dimensions"]
re_resolution <- reclustering_eval_subdf_order[1, "resolution"]


rec_marker_qc_filtering_output_dir <- sprintf("%s/%s", qc_filtering_output_dir, integration_method)
print(sprintf("normalization method: %s", global_normalization))
return_only_var_genes <- ifelse(global_normalization %in% c("SCTransform", "SCTransform_v2"), FALSE, NA)
## normalization and feature selection
if (merge_before_norm) {
    returned_list <- seurat_merge_normalization(seurat_object_merge, global_normalization, sct_variables_to_regress, global_features, return_only_var_genes, mitochondrial_geneset_spe, hemoglobin_geneset_spe, rec_marker_qc_filtering_output_dir)
    seurat_object_merge_norm <- returned_list[["seurat_merge"]]
    integration_features_nb <- returned_list[["features_nb"]]
    integration_normalization_method <- returned_list[["norm_method"]]
    default_assay <- returned_list[["assay"]]
    selection_method <- returned_list[["selec"]]
    normalization_output_dir <- returned_list[["out_dir"]]
    normalization_output_name <- returned_list[["out_name"]]
    rm(returned_list)
    gc()

    ### HVG plots
    hvg_df <- hvg_plots(seurat_object_merge_norm, default_assay, selection_method, global_features, normalization_output_dir, "HVGs")
} else {
    returned_list <- seurat_list_normalization(seurat_objects, global_normalization, sct_variables_to_regress, global_features, return_only_var_genes, mitochondrial_geneset_spe, hemoglobin_geneset_spe, rec_marker_qc_filtering_output_dir)
    seurat_objects_norm <- returned_list[["seurat_list"]]
    integration_features_nb <- returned_list[["features_nb"]]
    integration_normalization_method <- returned_list[["norm_method"]]
    default_assay <- returned_list[["assay"]]
    selection_method <- returned_list[["selec"]]
    normalization_output_dir <- returned_list[["out_dir"]]
    normalization_output_name <- returned_list[["out_name"]]
    rm(returned_list)
    gc()

    ### HVG plots
    hvg_df_list <- lapply(seurat_objects_norm, hvg_plots, assay2use=default_assay, selec=selection_method, nb_hvgs=global_features, out_dir=normalization_output_dir, out_name=sprintf("%s_HVGs", normalization_output_name))

    ### merge datasets
    seurat_object_merge_norm <- merge(seurat_objects_norm[[1]], y=unlist(seurat_objects_norm)[2:nb_sections], add.cell.ids=sections)
    #### orig.ident, lame, zone and condition are not factors in the merged object meta data
    seurat_object_merge_norm@meta.data$orig.ident <- as.factor(seurat_object_merge_norm@meta.data$orig.ident)
    seurat_object_merge_norm@meta.data$lame <- as.factor(seurat_object_merge_norm@meta.data$lame)
    seurat_object_merge_norm@meta.data$zone <- as.factor(seurat_object_merge_norm@meta.data$zone)
    seurat_object_merge_norm@meta.data$condition <- as.factor(seurat_object_merge_norm@meta.data$condition)
    #### set variable features: after merging datasets, variable features are not set in the merged Seurat object, which makes RunPCA() fails: Error in PrepDR(object = object, features = features, verbose = verbose): Variable features haven't been set. Run FindVariableFeatures() or provide a vector of feature names.
    ##### get the union of the HVGs of each dataset
    if (global_normalization %in% c("SCTransform", "SCTransform_v2")) {
        seurat_object_merge_norm_hvgs <- get_HVGs_merged_Seurat_SCT(seurat_object_merge_norm, global_features)
    } else {
        seurat_object_merge_norm_hvgs <- get_HVGs_merged_Seurat_LogNormalize(seurat_object_merge_norm, seurat_objects_norm, global_features)
    }
    VariableFeatures(seurat_object_merge_norm) <- seurat_object_merge_norm_hvgs
    rm(seurat_objects_norm)
    gc()

    ### QC plots after normalization
    #### feature stat plots
    default_assay <- DefaultAssay(seurat_object_merge_norm)
    ## add feature statistics to merged Seurat object
    feature_stat_slot <- ifelse(global_normalization == "LogNormalize", "data", "counts")
    seurat_object_merge_norm <- addFeatureStats2Seurat(seurat_object_merge_norm, default_assay, feature_stat_slot)
    norm_column_prefix <- "norm"
    seurat_object_merge_norm <- addPercentageFeatureSet2Seurat(seurat_object_merge_norm, mitochondrial_geneset_spe, hemoglobin_geneset_spe, default_assay, norm_column_prefix)
    library(RColorBrewer)
    seurat_all_feature_stat_plots(seurat_object_merge_norm, default_assay, norm_column_prefix, "Normalized counts", normalization_output_dir, "QC_stats_normalization") 
    detach("package:RColorBrewer")
    #### mitochondrial and hemoglobin gene spatial plots after filerting
    ##### hemoglobin genes
    df_list <- qc_spatial_mean_exp_plots(seurat_object_merge_norm, hemoglobin_geneset_spe, TRUE, normalization_output_dir, "normalized_hemoglobin_genes")
    ##### mitochondrial genes
    df_list <- qc_spatial_mean_exp_plots(seurat_object_merge_norm, mitochondrial_geneset_spe, TRUE, normalization_output_dir, "normalized_mitochondrial_genes")
}

## integration
integration_output_dir <- sprintf("%s/%s", normalization_output_dir, gsub(" ", "_", integration_method))
if (! dir.exists(integration_output_dir)) {
    dir.create(integration_output_dir, recursive=TRUE, mode="0775")
}

if (integration_method == "Seurat") {
    ### Seurat integration
    ## split merged Seurat object into list of Seurat objects for integration methods
    seurat_objects_norm <- SplitObject(seurat_object_merge_norm, split.by="orig.ident")
    ### after splitting, all the objects still have all the images: only keep the image of the corresponding sample
    for (one_sample in names(seurat_objects_norm)) {
        seurat_objects_norm[[one_sample]]@images <- list(seurat_objects_norm[[one_sample]]@images[[one_sample]])
        names(seurat_objects_norm[[one_sample]]@images) <- one_sample
    }
    ### after splitting, all the objects have 4 models in SCTModel.list slot but only one of them has cells: only keep this one
    for (one_sample in names(seurat_objects_norm)) {
        for (one_model in names(seurat_objects_norm[[one_sample]]@assays$SCT@SCTModel.list)) {
            if (dim(seurat_objects_norm[[one_sample]]@assays$SCT@SCTModel.list[[one_model]]@cell.attributes)[1] > 0) {
                seurat_objects_norm[[one_sample]]@assays$SCT@SCTModel.list <- list(seurat_objects_norm[[one_sample]]@assays$SCT@SCTModel.list[[one_model]])
                names(seurat_objects_norm[[one_sample]]@assays$SCT@SCTModel.list) <- one_model
                break
            }
        }
    }

    #### integration
    integration_output_name <- sprintf("%s_%s", normalization_output_name, integration_method)
    features <- SelectIntegrationFeatures(seurat_objects_norm, nfeatures=integration_features_nb)
    if (global_normalization %in% c("SCTransform", "SCTransform_v2")) {
        seurat_objects_norm <- PrepSCTIntegration(object.list=seurat_objects_norm, anchor.features=features)
    }
    int.anchors <- FindIntegrationAnchors(object.list=seurat_objects_norm, anchor.features=features, normalization.method=integration_normalization_method, dims=1:global_dimensions)
    rm(seurat_objects_norm)
    gc()
    seurat_object_to_scale <- IntegrateData(anchorset=int.anchors, normalization.method=integration_normalization_method, dims=1:global_dimensions)
    rm(int.anchors)
    gc()
    ##### orig.ident, lame, zone and condition are not factors in the integrated object meta data
    seurat_object_to_scale@meta.data$orig.ident <- as.factor(seurat_object_to_scale@meta.data$orig.ident)
    seurat_object_to_scale@meta.data$lame <- as.factor(seurat_object_to_scale@meta.data$lame)
    seurat_object_to_scale@meta.data$zone <- as.factor(seurat_object_to_scale@meta.data$zone)
    seurat_object_to_scale@meta.data$condition <- as.factor(seurat_object_to_scale@meta.data$condition)
} else {
    seurat_object_to_scale <- seurat_object_merge_norm
    rm(seurat_objects_norm)
    gc()
}
rm(seurat_object_merge_norm)
gc()

### scaling
if (global_normalization == "LogNormalize") {
    if (! is.na(scale_data_split_by)) {
        if (!is.na(scale_data_variables_to_regress)) {
            seurat_object_integration <- ScaleData(seurat_object_to_scale, split.by=scale_data_split_by, vars.to.regress=scale_data_variables_to_regress)
            integration_output_dir <- sprintf("%s/%s_%s", integration_output_dir, gsub("[.]", "", scale_data_split_by), paste(gsub("[.]", "", scale_data_variables_to_regress), collapse="-"))
        } else {
            seurat_object_integration <- ScaleData(seurat_object_to_scale, split.by=scale_data_split_by)
            integration_output_dir <- sprintf("%s/%s", integration_output_dir, gsub("[.]", "", scale_data_split_by))
        }
    } else {
        seurat_object_integration <- ScaleData(seurat_object_to_scale)
    }
} else {
    seurat_object_integration <- seurat_object_to_scale
}
rm(seurat_object_to_scale)
gc()

### dimensionality reduction
seurat_object_integration <- RunPCA(seurat_object_integration, verbose=FALSE)
print(sprintf("dimension number: %d", global_dimensions))
dimensions_output_dir <- sprintf("%s/%d_dims", integration_output_dir, global_dimensions)
if (! dir.exists(dimensions_output_dir)) {
    dir.create(dimensions_output_dir, recursive=TRUE, mode="0775")
}
dimensions_output_name <- sprintf("%s_%s_%ddims", normalization_output_name, gsub(" ", "_", integration_method), global_dimensions)

#### Harmony integration and dimensionality reduction
if (integration_method == "Harmony") {
    seurat_object_to_cluster <- RunHarmony(seurat_object_integration, group.by.vars=integration_group, dims=1:global_dimensions)
    reduction2use <- "harmony"
} else {
    seurat_object_to_cluster <- seurat_object_integration
    reduction2use <- "pca"
}
rm(seurat_object_integration)
gc()

seurat_object_clustering <- FindNeighbors(seurat_object_to_cluster, reduction=reduction2use, dims=1:global_dimensions)

## global clustering
print(sprintf("global clustering resolution: %.1f", global_resolution))
resolution_output_dir <- sprintf("%s/res%s", dimensions_output_dir, sub("[.]", "_", global_resolution))
if (! dir.exists(resolution_output_dir)) {
    dir.create(resolution_output_dir, recursive=TRUE, mode="0775")
}
resolution_output_name <- sprintf("%s_res%s", dimensions_output_name, sub("[.]", "_", global_resolution))

seurat_object_clustering <- clustering(seurat_object_clustering, global_resolution, reduction2use, global_dimensions)

write.csv(seurat_object_clustering@meta.data, file=sprintf("%s/%s_metadata.csv", resolution_output_dir, resolution_output_name), quote=FALSE, row.names=TRUE)
### plot clusters onto UMAP or onto the tissue section
clustering_plots(seurat_object_clustering, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integrated datasets: %s, all spots", integration_method), resolution_output_dir, sprintf("%s_clustering", resolution_output_name))

### cluster markers, differentially expressed genes across conditions and conserved cell type markers
if (global_normalization %in% c("SCTransform", "SCTransform_v2")) {
    different_SCT_models <- ifelse(length(seurat_object_clustering@assays$SCT@SCTModel.list) > 1, TRUE, FALSE)
} else {
    different_SCT_models <- NA
}
marker_analysis(seurat_object_clustering, default_assay, different_SCT_models, mitochondrial_geneset_spe, hemoglobin_geneset_spe, "negbinom", "orig.ident", cluster_markers, scale_data_split_by, scale_data_variables_to_regress, load_image, nb_cluster_markers_to_plot, nb_SpatialFeaturePlot_per_page, DE_analysis, nb_DEGs_to_plot, nb_DEGs_per_page, conserved_markers, "seurat_clusters", "condition", resolution_output_dir, resolution_output_name)


## hippocampus reclusteing
### get hippocampus spots
barcodes2use <- rownames(seurat_object_clustering@meta.data[which(seurat_object_clustering@meta.data$seurat_clusters == hippocampus_cluster),])
rm(seurat_object_clustering)
gc()

### reclustering
seurat_object_to_cluster_hippocampus <- subset(seurat_object_to_cluster, cells=barcodes2use)
rm(seurat_object_to_cluster)
gc()

print(sprintf("dimension number: %d", re_dimensions))
seurat_object_to_cluster_hippocampus_graph <- FindNeighbors(seurat_object_to_cluster_hippocampus, reduction=reduction2use, dims=1:re_dimensions)


print(sprintf("reclustering resolution: %.1f", re_resolution))
resolution_output_dir <- sprintf("%s/rec/%d_dims/res%s", dimensions_output_dir, re_dimensions, sub("[.]", "_", re_resolution))
if (! dir.exists(resolution_output_dir)) {
    dir.create(resolution_output_dir, recursive=TRUE, mode="0775")
}
resolution_output_name <- sprintf("%s_rec_%d_dims_res%s", dimensions_output_name, re_dimensions, sub("[.]", "_", re_resolution))

seurat_object_clustering_hippocampus <- clustering(seurat_object_to_cluster_hippocampus_graph, re_resolution, reduction2use, re_dimensions)
write.csv(seurat_object_clustering_hippocampus@meta.data, file=sprintf("%s/%s_metadata.csv", resolution_output_dir, resolution_output_name), quote=FALSE, row.names=TRUE)
rm(seurat_object_to_cluster_hippocampus_graph)
gc()

hippocampus_cluster_mean_silhouette <- mean(seurat_object_clustering_hippocampus$silhouette_seurat_clusters)
hippocampus_cluster_median_silhouette <- median(seurat_object_clustering_hippocampus$silhouette_seurat_clusters)
### plot clusters onto UMAP or onto the tissue section
clustering_plots(seurat_object_clustering_hippocampus, load_image, nb_SpatialFeaturePlot_per_page, "Non-integrated datasets, all spots", resolution_output_dir, sprintf("%s_clustering", resolution_output_name))

### cluster markers, differentially expressed genes across conditions and conserved cell type markers
if (global_normalization %in% c("SCTransform", "SCTransform_v2")) {
    different_SCT_models <- ifelse(length(seurat_object_clustering_hippocampus@assays$SCT@SCTModel.list) > 1, TRUE, FALSE)
} else {
    different_SCT_models <- NA
}
marker_analysis(seurat_object_clustering_hippocampus, default_assay, different_SCT_models, mitochondrial_geneset_spe, hemoglobin_geneset_spe, "negbinom", "orig.ident", cluster_markers, scale_data_split_by, scale_data_variables_to_regress, load_image, nb_cluster_markers_to_plot, nb_SpatialFeaturePlot_per_page, DE_analysis, nb_DEGs_to_plot, nb_DEGs_per_page, conserved_markers, "seurat_clusters", "condition", resolution_output_dir, resolution_output_name)

rm(seurat_object_clustering_hippocampus)
gc()






# GO term analysis
library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
de_results_file <- sprintf("%s/10-DEGs/negbinom/origident/SCTransform_v2_featnb500_Seurat_15dims_res0_1_DEGs_cluster_7.csv", resolution_output_dir)
de_results_df <- read.csv(de_results_file, row.names=1)

## GO term enrichment analysis
de_genes <- rownames(de_results_df[which(de_results_df$p_val_adj < 0.05),])
up_genes <- rownames(de_results_df[which(de_results_df$p_val_adj < 0.05 & de_results_df$avg_log2FC > 0),])
down_genes <- rownames(de_results_df[which(de_results_df$p_val_adj < 0.05 & de_results_df$avg_log2FC < 0),])
universe_genes <- rownames(de_results_df)

up_genes_entrez_df <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Rn.eg.db", drop=FALSE)
down_genes_entrez_df <- bitr(down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Rn.eg.db", drop=FALSE)
universe_genes_entrez_df <- bitr(universe_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Rn.eg.db", drop=FALSE)

up_ego <- enrichGO(gene=up_genes_entrez_df$ENTREZID[!is.na(up_genes_entrez_df$ENTREZID)], universe=universe_genes_entrez_df$ENTREZID[!is.na(universe_genes_entrez_df$ENTREZID)], OrgDb=org.Rn.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05, readable=TRUE)
barplot(up_ego, showCategory=20)

down_ego <- enrichGO(gene=down_genes_entrez_df$ENTREZID[!is.na(down_genes_entrez_df$ENTREZID)], universe=universe_genes_entrez_df$ENTREZID[!is.na(universe_genes_entrez_df$ENTREZID)], OrgDb=org.Rn.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05, readable=TRUE)
barplot(down_ego, showCategory=20)

## GSEA
up_results_df <- de_results_df[which(de_results_df$p_val_adj < 0.05 & de_results_df$avg_log2FC > 0),]
up_results_entrez_df <- merge(up_genes_entrez_df, up_results_df, by.x="SYMBOL", by.y="row.names")
up_geneList <- up_results_entrez_df$avg_log2FC
names(up_geneList) <- up_results_entrez_df$ENTREZID
up_geneList <- sort(up_geneList, decreasing=TRUE)
up_gse <- gseGO(geneList=up_geneList, OrgDb=org.Rn.eg.db, ont="BP", minGSSize=100, maxGSSize=500, pvalueCutoff=0.05, verbose=FALSE)

down_results_df <- de_results_df[which(de_results_df$p_val_adj < 0.05 & de_results_df$avg_log2FC < 0),]
down_results_entrez_df <- merge(down_genes_entrez_df, down_results_df, by.x="SYMBOL", by.y="row.names")
down_geneList <- down_results_entrez_df$avg_log2FC
names(down_geneList) <- down_results_entrez_df$ENTREZID
down_geneList <- sort(down_geneList, decreasing=TRUE)
down_gse <- gseGO(geneList=down_geneList, OrgDb=org.Rn.eg.db, ont="BP", minGSSize=100, maxGSSize=500, pvalueCutoff=0.05, verbose=FALSE)





# Harmony
library(harmony)
integration_method <- "Harmony"
integration_methods <- c(integration_methods, integration_method)
## get top results for global clustering with Seurat integration
global_clustering_eval_subdf <- global_clustering_eval_df[which(global_clustering_eval_df$normalization=="SCTransform_v2" & global_clustering_eval_df$integration == integration_method),]
global_clustering_eval_subdf_order <- global_clustering_eval_subdf[order(global_clustering_eval_subdf$mean_silhouette_seurat_clusters, decreasing=TRUE),]
global_normalization <- global_clustering_eval_subdf_order[1, "normalization"]
global_features <- global_clustering_eval_subdf_order[1, "features"]
global_integration <- global_clustering_eval_subdf_order[1, "integration"]
global_dimensions <- global_clustering_eval_subdf_order[1, "dimensions"]
global_resolution <- global_clustering_eval_subdf_order[1, "resolution"]
## set global clustering parameters
normalization_method <- global_normalization
features_param <- global_features
integration_method <- global_integration
dimensions_nb_param <- global_dimensions
clustering_resolution <- global_resolution
hippocampus_cluster <- "3"

rec_marker_qc_filtering_output_dir <- sprintf("%s/%s", qc_filtering_output_dir, integration_method)
print(sprintf("normalization method: %s", normalization_method))
return_only_var_genes <- ifelse(normalization_method %in% c("SCTransform", "SCTransform_v2"), FALSE, NA)
## normalization and feature selection
if (merge_before_norm) {
    returned_list <- seurat_merge_normalization(seurat_object_merge, normalization_method, sct_variables_to_regress, features_param, return_only_var_genes, mitochondrial_geneset_spe, hemoglobin_geneset_spe, rec_marker_qc_filtering_output_dir)
    seurat_object_merge_norm <- returned_list[["seurat_merge"]]
    integration_features_nb <- returned_list[["features_nb"]]
    integration_normalization_method <- returned_list[["norm_method"]]
    default_assay <- returned_list[["assay"]]
    selection_method <- returned_list[["selec"]]
    normalization_output_dir <- returned_list[["out_dir"]]
    normalization_output_name <- returned_list[["out_name"]]
    rm(returned_list)
    gc()

    ### HVG plots
    hvg_df <- hvg_plots(seurat_object_merge_norm, default_assay, selection_method, features_param, normalization_output_dir, "HVGs")
} else {
    returned_list <- seurat_list_normalization(seurat_objects, normalization_method, sct_variables_to_regress, features_param, return_only_var_genes, mitochondrial_geneset_spe, hemoglobin_geneset_spe, rec_marker_qc_filtering_output_dir)
    seurat_objects_norm <- returned_list[["seurat_list"]]
    integration_features_nb <- returned_list[["features_nb"]]
    integration_normalization_method <- returned_list[["norm_method"]]
    default_assay <- returned_list[["assay"]]
    selection_method <- returned_list[["selec"]]
    normalization_output_dir <- returned_list[["out_dir"]]
    normalization_output_name <- returned_list[["out_name"]]
    rm(returned_list)
    gc()

    ### HVG plots
    hvg_df_list <- lapply(seurat_objects_norm, hvg_plots, assay2use=default_assay, selec=selection_method, nb_hvgs=features_param, out_dir=normalization_output_dir, out_name=sprintf("%s_HVGs", normalization_output_name))

    ### merge datasets
    seurat_object_merge_norm <- merge(seurat_objects_norm[[1]], y=unlist(seurat_objects_norm)[2:nb_sections], add.cell.ids=sections)
    #### orig.ident, lame, zone and condition are not factors in the merged object meta data
    seurat_object_merge_norm@meta.data$orig.ident <- as.factor(seurat_object_merge_norm@meta.data$orig.ident)
    seurat_object_merge_norm@meta.data$lame <- as.factor(seurat_object_merge_norm@meta.data$lame)
    seurat_object_merge_norm@meta.data$zone <- as.factor(seurat_object_merge_norm@meta.data$zone)
    seurat_object_merge_norm@meta.data$condition <- as.factor(seurat_object_merge_norm@meta.data$condition)
    #### set variable features: after merging datasets, variable features are not set in the merged Seurat object, which makes RunPCA() fails: Error in PrepDR(object = object, features = features, verbose = verbose): Variable features haven't been set. Run FindVariableFeatures() or provide a vector of feature names.
    ##### get the union of the HVGs of each dataset
    if (normalization_method %in% c("SCTransform", "SCTransform_v2")) {
        seurat_object_merge_norm_hvgs <- get_HVGs_merged_Seurat_SCT(seurat_object_merge_norm, features_param)
    } else {
        seurat_object_merge_norm_hvgs <- get_HVGs_merged_Seurat_LogNormalize(seurat_object_merge_norm, seurat_objects_norm, features_param)
    }
    VariableFeatures(seurat_object_merge_norm) <- seurat_object_merge_norm_hvgs
    rm(seurat_objects_norm)
    gc()

    ### QC plots after normalization
    #### feature stat plots
    default_assay <- DefaultAssay(seurat_object_merge_norm)
    ## add feature statistics to merged Seurat object
    feature_stat_slot <- ifelse(normalization_method == "LogNormalize", "data", "counts")
    seurat_object_merge_norm <- addFeatureStats2Seurat(seurat_object_merge_norm, default_assay, feature_stat_slot)
    norm_column_prefix <- "norm"
    seurat_object_merge_norm <- addPercentageFeatureSet2Seurat(seurat_object_merge_norm, mitochondrial_geneset_spe, hemoglobin_geneset_spe, default_assay, norm_column_prefix)
    library(RColorBrewer)
    seurat_all_feature_stat_plots(seurat_object_merge_norm, default_assay, norm_column_prefix, "Normalized counts", normalization_output_dir, "QC_stats_normalization") 
    detach("package:RColorBrewer")
    #### mitochondrial and hemoglobin gene spatial plots after filerting
    ##### hemoglobin genes
    df_list <- qc_spatial_mean_exp_plots(seurat_object_merge_norm, hemoglobin_geneset_spe, TRUE, normalization_output_dir, "normalized_hemoglobin_genes")
    ##### mitochondrial genes
    df_list <- qc_spatial_mean_exp_plots(seurat_object_merge_norm, mitochondrial_geneset_spe, TRUE, normalization_output_dir, "normalized_mitochondrial_genes")
}

## integration
integration_output_dir <- sprintf("%s/%s", normalization_output_dir, gsub(" ", "_", integration_method))
if (! dir.exists(integration_output_dir)) {
    dir.create(integration_output_dir, recursive=TRUE, mode="0775")
}

if (integration_method == "Seurat") {
    ### Seurat integration
    ## split merged Seurat object into list of Seurat objects for integration methods
    seurat_objects_norm <- SplitObject(seurat_object_merge_norm, split.by="orig.ident")
    ### after splitting, all the objects still have all the images: only keep the image of the corresponding sample
    for (one_sample in names(seurat_objects_norm)) {
        seurat_objects_norm[[one_sample]]@images <- list(seurat_objects_norm[[one_sample]]@images[[one_sample]])
        names(seurat_objects_norm[[one_sample]]@images) <- one_sample
    }

    #### integration
    integration_output_name <- sprintf("%s_%s", normalization_output_name, integration_method)
    features <- SelectIntegrationFeatures(seurat_objects_norm, nfeatures=integration_features_nb)
    if (normalization_method %in% c("SCTransform", "SCTransform_v2")) {
        seurat_objects_norm <- PrepSCTIntegration(object.list=seurat_objects_norm, anchor.features=features)
    }
    int.anchors <- FindIntegrationAnchors(object.list=seurat_objects_norm, anchor.features=features, normalization.method=integration_normalization_method, dims=1:dimensions_nb_param)
    rm(seurat_objects_norm)
    gc()
    seurat_object_to_scale <- IntegrateData(anchorset=int.anchors, normalization.method=integration_normalization_method, dims=1:dimensions_nb_param)
    rm(int.anchors)
    gc()
    ##### orig.ident, lame, zone and condition are not factors in the integrated object meta data
    seurat_object_to_scale@meta.data$orig.ident <- as.factor(seurat_object_to_scale@meta.data$orig.ident)
    seurat_object_to_scale@meta.data$lame <- as.factor(seurat_object_to_scale@meta.data$lame)
    seurat_object_to_scale@meta.data$zone <- as.factor(seurat_object_to_scale@meta.data$zone)
    seurat_object_to_scale@meta.data$condition <- as.factor(seurat_object_to_scale@meta.data$condition)
} else {
    seurat_object_to_scale <- seurat_object_merge_norm
    rm(seurat_objects_norm)
    gc()
}
rm(seurat_object_merge_norm)
gc()

### scaling
if (normalization_method == "LogNormalize") {
    if (! is.na(scale_data_split_by)) {
        if (!is.na(scale_data_variables_to_regress)) {
            seurat_object_integration <- ScaleData(seurat_object_to_scale, split.by=scale_data_split_by, vars.to.regress=scale_data_variables_to_regress)
            integration_output_dir <- sprintf("%s/%s_%s", integration_output_dir, gsub("[.]", "", scale_data_split_by), paste(gsub("[.]", "", scale_data_variables_to_regress), collapse="-"))
        } else {
            seurat_object_integration <- ScaleData(seurat_object_to_scale, split.by=scale_data_split_by)
            integration_output_dir <- sprintf("%s/%s", integration_output_dir, gsub("[.]", "", scale_data_split_by))
        }
    } else {
        seurat_object_integration <- ScaleData(seurat_object_to_scale)
    }
} else {
    seurat_object_integration <- seurat_object_to_scale
}
rm(seurat_object_to_scale)
gc()

### dimensionality reduction
seurat_object_integration <- RunPCA(seurat_object_integration, verbose=FALSE)
print(sprintf("dimension number: %d", dimensions_nb_param))
dimensions_output_dir <- sprintf("%s/%d_dims", integration_output_dir, dimensions_nb_param)
if (! dir.exists(dimensions_output_dir)) {
    dir.create(dimensions_output_dir, recursive=TRUE, mode="0775")
}
dimensions_output_name <- sprintf("%s_%s_%ddims", normalization_output_name, gsub(" ", "_", integration_method), dimensions_nb_param)

#### Harmony integration and dimensionality reduction
if (integration_method == "Harmony") {
    seurat_object_to_cluster <- RunHarmony(seurat_object_integration, group.by.vars=integration_group, dims=1:dimensions_nb_param)
    reduction2use <- "harmony"
} else {
    seurat_object_to_cluster <- seurat_object_integration
    reduction2use <- "pca"
}
rm(seurat_object_integration)
gc()

seurat_object_clustering <- FindNeighbors(seurat_object_to_cluster, reduction=reduction2use, dims=1:dimensions_nb_param)

## global clustering
print(sprintf("global clustering resolution: %.1f", global_resolution))
resolution_output_dir <- sprintf("%s/res%s", dimensions_output_dir, sub("[.]", "_", global_resolution))
if (! dir.exists(resolution_output_dir)) {
    dir.create(resolution_output_dir, recursive=TRUE, mode="0775")
}
resolution_output_name <- sprintf("%s_res%s", dimensions_output_name, sub("[.]", "_", global_resolution))

seurat_object_clustering <- clustering(seurat_object_clustering, global_resolution, reduction2use, dimensions_nb_param)

write.csv(seurat_object_clustering@meta.data, file=sprintf("%s/%s_metadata.csv", resolution_output_dir, resolution_output_name), quote=FALSE, row.names=TRUE)
### plot clusters onto UMAP or onto the tissue section
clustering_plots(seurat_object_clustering, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integrated datasets: %s, all spots", integration_method), resolution_output_dir, sprintf("%s_clustering", resolution_output_name))


## hippocampus reclusteing
### get hippocampus spots
barcodes2use <- rownames(seurat_object_clustering@meta.data[which(seurat_object_clustering@meta.data$seurat_clusters == hippocampus_cluster),])
rm(seurat_object_clustering)
gc()

### reclustering
seurat_object_to_cluster_hippocampus <- subset(seurat_object_to_cluster, cells=barcodes2use)
rm(seurat_object_to_cluster)
gc()

for (dimensions_nb_reclustering in reclustering_dimension_vector) {
    print(sprintf("dimension number: %d", dimensions_nb_reclustering))
    seurat_object_to_cluster_hippocampus_graph <- FindNeighbors(seurat_object_to_cluster_hippocampus, reduction=reduction2use, dims=1:dimensions_nb_reclustering)

    for (clustering_resolution in reclustering_resolution_vector) {
        print(sprintf("reclustering resolution: %.1f", clustering_resolution))
        resolution_output_dir <- sprintf("%s/rec/%d_dims/res%s", dimensions_output_dir, dimensions_nb_reclustering, sub("[.]", "_", clustering_resolution))
        if (! dir.exists(resolution_output_dir)) {
            dir.create(resolution_output_dir, recursive=TRUE, mode="0775")
        }
        resolution_output_name <- sprintf("%s_rec_%d_dims_res%s", dimensions_output_name, dimensions_nb_reclustering, sub("[.]", "_", clustering_resolution))

        seurat_object_clustering_hippocampus <- clustering(seurat_object_to_cluster_hippocampus_graph, clustering_resolution, reduction2use, dimensions_nb_reclustering)
        write.csv(seurat_object_clustering_hippocampus@meta.data, file=sprintf("%s/%s_metadata.csv", resolution_output_dir, resolution_output_name), quote=FALSE, row.names=TRUE)
        hippocampus_cluster_mean_silhouette <- mean(seurat_object_clustering_hippocampus$silhouette_seurat_clusters)
        hippocampus_cluster_median_silhouette <- median(seurat_object_clustering_hippocampus$silhouette_seurat_clusters)
        ### plot clusters onto UMAP or onto the tissue section
        clustering_plots(seurat_object_clustering_hippocampus, load_image, nb_SpatialFeaturePlot_per_page, "Non-integrated datasets, all spots", resolution_output_dir, sprintf("%s_clustering", resolution_output_name))

        ### get mean silhouette scores: all spots and hippocampus spots
        mean_median_silhouette_df <- data.frame(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_reclustering, resolution=clustering_resolution, mean_silhouette_seurat_clusters_hippocampus=hippocampus_cluster_mean_silhouette, median_silhouette_seurat_clusters_hippocampus=hippocampus_cluster_median_silhouette)
        hippocampus_reclustering_mean_silhouette_df <- rbind(hippocampus_reclustering_mean_silhouette_df, mean_median_silhouette_df)

        ### cluster markers, differentially expressed genes across conditions and conserved cell type markers
        if (normalization_method %in% c("SCTransform", "SCTransform_v2")) {
            different_SCT_models <- ifelse(length(seurat_object_clustering_hippocampus@assays$SCT@SCTModel.list) > 1, TRUE, FALSE)
        } else {
            different_SCT_models <- NA
        }
        marker_analysis(seurat_object_clustering_hippocampus, default_assay, different_SCT_models, mitochondrial_geneset_spe, hemoglobin_geneset_spe, "negbinom", "orig.ident", cluster_markers, scale_data_split_by, scale_data_variables_to_regress, load_image, nb_cluster_markers_to_plot, nb_SpatialFeaturePlot_per_page, DE_analysis, nb_DEGs_to_plot, nb_DEGs_per_page, conserved_markers, "seurat_clusters", "condition", resolution_output_dir, resolution_output_name)

        rm(seurat_object_clustering_hippocampus)
        gc()
    }
    rm(seurat_object_to_cluster_hippocampus_graph)
    gc()
}
detach("package:harmony", unload=TRUE)
rm(seurat_object_to_cluster_hippocampus)
gc()
if (! merge_before_norm) {
    rm(seurat_objects)
    gc()
}



hippocampus_reclustering_mean_silhouette_df$integration <- factor(hippocampus_reclustering_mean_silhouette_df$integration, levels=integration_methods)
hippocampus_reclustering_mean_silhouette_df$features <- as.factor(hippocampus_reclustering_mean_silhouette_df$features)
hippocampus_reclustering_mean_silhouette_df$dimensions <- as.factor(hippocampus_reclustering_mean_silhouette_df$dimensions)
hippocampus_reclustering_mean_silhouette_df$resolution <- as.factor(hippocampus_reclustering_mean_silhouette_df$resolution)
plot_output_name <- sprintf("%s_norm_featnb_int_dims_res_reclustering_mean_silhouette", output_dirname)
write.csv(hippocampus_reclustering_mean_silhouette_df, file=sprintf("%s/%s.csv", qc_filtering_output_dir, plot_output_name), quote=FALSE, row.names=FALSE)
pdf(sprintf("%s/%s.pdf", qc_filtering_output_dir, plot_output_name))
one_category <- "seurat_clusters_hippocampus"
p <- ggplot(hippocampus_reclustering_mean_silhouette_df, aes(x=resolution, y=.data[[sprintf("mean_silhouette_%s", one_category)]], fill=integration)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(dimensions~features) +
    labs(title="Hippocampus cluster reclusterint mean silhouette", x="Resolution", y="Mean silhouette", fill="Integration method") +
    theme_bw() +
    theme(legend.position="bottom") +
    theme(axis.text.x=element_text(angle=60, hjust=1)) +
    theme(panel.border=element_rect(color="grey50"))
print(p)
dev.off()





