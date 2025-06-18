.libPaths(c(.libPaths(), "/home/christophe.lepriol/NeuroDev_ADD/R/r_4.1.0"))
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(cluster) # silhouette() function


##############
# Parameters #
##############
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"
src_dir <- sprintf("%s/20-Count_analysis/10-EpiReg/src", work_dir)
genome_name <- "RefSeq108"
load_image <- TRUE

epireg_visium_metadata_df <- data.frame(sample=c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16"),
                                       condition=factor(c(rep("SE", 4), rep("CTRL", 4), rep("SE", 4), rep("CTRL", 4))), 
                                       time=factor(c(rep(c(5, 10), 4), rep(c(20, 40, 40, 20), 2))))

# sections
sections <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")
timepoint <- NA
#timepoint <- 5

# parameters
## count matrix: Space Ranger or ST Pipeline
input_matrix <- "Space Ranger" # allowed values: 'ST Pipeline' or 'Space Ranger'
## QC
spot_sum_of_counts_threshold <- 1000
spot_detected_genes_threshold <- 500
spot_mito_gene_count_pct_threshold <- 100
spot_hb_gene_count_pct_threshold <- 100
## normalization, feature selection, dimension reduction, clustering
#normalization_method_vector <- c("SCTransform", "SCTransform_v2", "LogNormalize")
normalization_method_vector <- c("SCTransform_v2")
sct_variables_to_regress <- "qc_mito_percent_RNA"
scale_data_variables_to_regress <- "qc_mito_percent_RNA"
variables_to_regress_in_path <- "mito"
scale_data_split_by <- "orig.ident"
scale_data_split_by_in_path <- "orig"
integration_method <- "Seurat"
integration_group <- "orig.ident"
integration_group_by_in_path <- "orig"
seed4Harmony <- NA

# number of samples per page when using SpatialFeaturePlot() function
nb_SpatialFeaturePlot_per_page <- 4
# number of samples per page when using SpatialFeaturePlot() function
nb_cluster_markers_to_plot <- 4

nb_DEGs_to_plot <- 30
nb_DEGs_per_page <- 3

metadata_categories <- c("lame", "zone", "condition", "time", "orig.ident")


#############
# Functions #
#############
source(sprintf("%s/6X-integration_functions.R", src_dir))


############
# Analysis #
############

# process parameters
## sections
if (! is.na(timepoint)) {
    sections <- epireg_visium_metadata_df[which(epireg_visium_metadata_df$time==timepoint), "sample"]
}
nb_sections <- length(sections)

## input matrix
if (input_matrix == "ST Pipeline") {
    input_dirname <- "00-ST_P"
    visium_params_in_path <- TRUE
} else {
    if (input_matrix == "Space Ranger") {
        input_dirname <- "10-SR"
        visium_params_in_path <- FALSE
    }
}
output_dir <- sprintf("%s/20-Count_analysis/10-EpiReg/output/%s/TSO_polyA_R1trim1_ov5_n2_min20/%s", work_dir, input_dirname, genome_name)
if (visium_params_in_path) {
    output_dir <- sprintf("%s/Visium_params", output_dir) 
}
output_dirname <- ifelse(is.na(timepoint), "all", sprintf("D%d", timepoint))
output_dir <- sprintf("%s/20-Integration/%s/sub_markers", output_dir, output_dirname)
if (! dir.exists(output_dir)) {
    dir.create(output_dir, recursive=TRUE, mode="0775")
}

# load data
print("load data")
library(SpatialExperiment)
library(ggspavis)
space_ranger_dir <- sprintf("%s/10-Raw2count/10-SR/output/TSO_polyA_R1trim1_ov5_n2_min20/%s/00-Samples", work_dir, genome_name)
if (input_matrix == "ST Pipeline") {
    st_pipeline_dir <- sprintf("%s/10-Raw2count/00-ST_P/output/TSO_polyA_R1trim1_ov5_n2_min20/%s/Visium_params/00-Samples", work_dir, genome_name)
    spe_objects <- lapply(sections, ST_Pipeline2SpatialExperiment, st_dir=st_pipeline_dir, sr_dir=space_ranger_dir, metadata_df=epireg_visium_metadata_df)
} else {
    if (input_matrix == "Space Ranger") {
        spe_objects <- lapply(sections, Space_Ranger2SpatialExperiment, sr_dir=space_ranger_dir, metadata_df=epireg_visium_metadata_df)
    }
}
names(spe_objects) <- sections

# quality control
library(scater) # Loading required package: scuttle
qc_output_dir <- sprintf("%s/Norm_%s-scale_%s-int_%s/QC-sum%d_det%s_mito%d_hb%d", output_dir, variables_to_regress_in_path, scale_data_split_by_in_path, integration_group_by_in_path, spot_sum_of_counts_threshold, spot_detected_genes_threshold, spot_mito_gene_count_pct_threshold, spot_hb_gene_count_pct_threshold)
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
if (input_matrix == "ST Pipeline") {
    mitochondrial_geneset <- gsub("-", ".", mt_genes) # '-' were replaced by '.' by ST Pipeline
} else {
    mitochondrial_geneset <- mt_genes
}
## hemoglobin genes
#hemoglobin_genes <- gtf_gene_regex(gtf_GRanges_gene, "description", "(hemo|neuro|myo|cyto|hapto| )globin", qc_output_dir)
hemoglobin_genes <- gtf_gene_regex(gtf_GRanges_gene, "description", "hemoglobin", qc_output_dir)
### add LOC120093065: https://www.ncbi.nlm.nih.gov/gene/120093065
hemoglobin_genes <- c(hemoglobin_genes, "LOC120093065")
if (input_matrix == "ST Pipeline") {
    hemoglobin_geneset <- gsub("-", ".", hemoglobin_genes) # '-' were replaced by '.' by ST Pipeline
} else {
    hemoglobin_geneset <- hemoglobin_genes
}
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
raw_column_prefix <- "raw"
spe_objects <- lapply(spe_objects, add_qc_stats, mito_genes=mitochondrial_geneset, hb_genes=hemoglobin_geneset, col_prefix=raw_column_prefix, out_dir=qc_stats_output_dir, out_name="QC_stats")
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
seurat_objects <- lapply(seurat_objects, addPercentageFeatureSet2Seurat, mito_genes=mitochondrial_geneset, hb_genes=hemoglobin_geneset, assay2use="RNA", col_prefix=raw_column_prefix)
library(RColorBrewer)
lapply(seurat_objects, seurat_all_feature_stat_plots, assay2use="RNA", col_prefix=raw_column_prefix, plot_title="Raw counts", out_dir=qc_stats_output_dir, out_name="QC_stats") 
detach("package:RColorBrewer")
### mitochondrial and hemoglobin gene spatial plots
#### hemoglobin genes
df_list <- lapply(seurat_objects, qc_spatial_mean_exp_plots, gene_vec=hemoglobin_geneset, mean_order=TRUE, out_dir=qc_stats_output_dir, out_name="hemoglobin_genes")
#### mitochondrial genes
df_list <- lapply(seurat_objects, qc_spatial_mean_exp_plots, gene_vec=mitochondrial_geneset, mean_order=TRUE, out_dir=qc_stats_output_dir, out_name="mitochondrial_genes")
rm(seurat_objects)
gc()
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
seurat_objects <- lapply(seurat_objects, addPercentageFeatureSet2Seurat, mito_genes=mitochondrial_geneset, hb_genes=hemoglobin_geneset, assay2use="RNA", col_prefix=qc_column_prefix)
library(RColorBrewer)
lapply(seurat_objects, seurat_all_feature_stat_plots, assay2use="RNA", col_prefix=qc_column_prefix, plot_title="Raw counts after QC", out_dir=qc_stats_output_dir, out_name="QC_stats_filtering") 
detach("package:RColorBrewer")
#### mitochondrial and hemoglobin gene spatial plots after filerting
##### hemoglobin genes
df_list <- lapply(seurat_objects, qc_spatial_mean_exp_plots, gene_vec=hemoglobin_geneset, mean_order=TRUE, out_dir=qc_filtering_output_dir, out_name="filtered_hemoglobin_genes")
##### mitochondrial genes
df_list <- lapply(seurat_objects, qc_spatial_mean_exp_plots, gene_vec=mitochondrial_geneset, mean_order=TRUE, out_dir=qc_filtering_output_dir, out_name="filtered_mitochondrial_genes")
detach("package:scater", unload=TRUE)
unloadNamespace("DropletUtils") # also unload DropletUtils package which imports ‘scuttle’ namespace
#detach("package:DropletUtils", unload=TRUE) # also unload DropletUtils package which imports ‘scuttle’ namespace
detach("package:scuttle", unload=TRUE)

# add feature statistics to Seurat objects
seurat_objects <- lapply(seurat_objects, addFeatureStats2Seurat, assay2use="RNA", slot2use="counts")


# parameter settings for whole clustering and subclusterings
## whole clustering
### get whole clustering parameter evaluation output
if (visium_params_in_path) {
    whole_clustering_output_dir <- sprintf("%s/20-Count_analysis/10-EpiReg/output/%s/TSO_polyA_R1trim1_ov5_n2_min20/%s/Visium_params/20-Integration/%s/whole/Norm_%s-scale_%s-int_%s/QC-sum%d_det%s_mito%d_hb%d/10-QC_filtering", work_dir, input_dirname, genome_name, output_dirname, variables_to_regress_in_path, scale_data_split_by_in_path, integration_group_by_in_path, spot_sum_of_counts_threshold, spot_detected_genes_threshold, spot_mito_gene_count_pct_threshold, spot_hb_gene_count_pct_threshold)
} else {
    whole_clustering_output_dir <- sprintf("%s/20-Count_analysis/10-EpiReg/output/%s/TSO_polyA_R1trim1_ov5_n2_min20/%s/20-Integration/%s/whole/Norm_%s-scale_%s-int_%s/QC-sum%d_det%s_mito%d_hb%d/10-QC_filtering", work_dir, input_dirname, genome_name, output_dirname, variables_to_regress_in_path, scale_data_split_by_in_path, integration_group_by_in_path, spot_sum_of_counts_threshold, spot_detected_genes_threshold, spot_mito_gene_count_pct_threshold, spot_hb_gene_count_pct_threshold)
}
whole_clustering_eval_df <- read.csv(sprintf("%s/%s_norm_featnb_int_dims_res_cluster_mean_silhouette.csv", whole_clustering_output_dir, output_dirname))

### set best whole clustering parameter setting
whole_params <- clustering_top_params(whole_clustering_eval_df, "mean_silhouette_seurat_clusters", integration_method)
whole_normalization <- whole_params[["normalization"]]
whole_features <- whole_params[["features"]]
whole_integration <- whole_params[["integration"]]
whole_dimensions <- whole_params[["dimensions"]]
whole_resolution <- whole_params[["resolution"]]
print(sprintf("normalization: %s, feature number: %d, integration: %s, dimension number: %d, resolution: %.2f", whole_normalization, whole_features, whole_integration, whole_dimensions, whole_resolution))

## subclusterings
### set subclusterings parameter settings
subsclusterings_parameter_settings <- list()
#### cluster 0
subclustering_name <- "0"
cluster_0_parameters <- list()
cluster_0_parameters[["Seurat"]]  <- subclustering_name
cluster_0_parameters[["normalization"]] <- "SCTransform_v2"
cluster_0_parameters[["features"]]  <- 250
cluster_0_parameters[["integration"]] <- "Seurat"
cluster_0_parameters[["dimensions"]]  <- 10
cluster_0_parameters[["resolution"]] <- 0.1
subsclusterings_parameter_settings[[subclustering_name]] <- cluster_0_parameters
#### cluster 1
subclustering_name <- "1"
cluster_1_parameters <- list()
cluster_1_parameters[["Seurat"]]  <- subclustering_name
cluster_1_parameters[["normalization"]] <- "SCTransform_v2"
cluster_1_parameters[["features"]]  <- 1500
cluster_1_parameters[["integration"]] <- "Seurat"
cluster_1_parameters[["dimensions"]]  <- 10
cluster_1_parameters[["resolution"]] <- 0.03
subsclusterings_parameter_settings[[subclustering_name]] <- cluster_1_parameters
#### cluster 2
subclustering_name <- "2"
cluster_2_parameters <- list()
cluster_2_parameters[["Seurat"]]  <- subclustering_name
cluster_2_parameters[["normalization"]] <- "SCTransform_v2"
cluster_2_parameters[["features"]]  <- 250
cluster_2_parameters[["integration"]] <- "Seurat"
cluster_2_parameters[["dimensions"]]  <- 15
cluster_2_parameters[["resolution"]] <- 0.5
subsclusterings_parameter_settings[[subclustering_name]] <- cluster_2_parameters
#### cluster 3
subclustering_name <- "3"
cluster_3_parameters <- list()
cluster_3_parameters[["Seurat"]]  <- subclustering_name
cluster_3_parameters[["normalization"]] <- "SCTransform_v2"
cluster_3_parameters[["features"]]  <- 2000
cluster_3_parameters[["integration"]] <- "Seurat"
cluster_3_parameters[["dimensions"]]  <- 10
cluster_3_parameters[["resolution"]] <- 0.05
subsclusterings_parameter_settings[[subclustering_name]] <- cluster_3_parameters
#### cluster 4
subclustering_name <- "4"
cluster_4_parameters <- list()
cluster_4_parameters[["Seurat"]]  <- subclustering_name
cluster_4_parameters[["normalization"]] <- "SCTransform_v2"
cluster_4_parameters[["features"]]  <- 500
cluster_4_parameters[["integration"]] <- "Seurat"
cluster_4_parameters[["dimensions"]]  <- 10
cluster_4_parameters[["resolution"]] <- 0.08
subsclusterings_parameter_settings[[subclustering_name]] <- cluster_4_parameters


# whole workflow: whole dataset normalization to subclustering DE analysis (treated vs control comparisons)
integration_output_dir <- sprintf("%s/%s", qc_filtering_output_dir, integration_method)
if (! dir.exists(integration_output_dir)) {
    dir.create(integration_output_dir, recursive=TRUE, mode="0775")
}

## whole clustering
whole_clustering_output_dir <- sprintf("%s/whole", integration_output_dir)
if (! dir.exists(whole_clustering_output_dir)) {
    dir.create(whole_clustering_output_dir, recursive=TRUE, mode="0775")
}

### set parameters
whole_param_list <- list()
whole_param_list$project_name <- output_dirname
whole_param_list$normalization_method <- whole_normalization
whole_param_list$normalization_variables_to_regress <- sct_variables_to_regress
whole_param_list$feature_nb <- whole_features
whole_param_list$mitochondrial_genes <- mitochondrial_geneset
whole_param_list$hemoglobin_genes <- hemoglobin_geneset
whole_param_list$integration_method <- whole_integration
whole_param_list$dimension_nb <- whole_dimensions
whole_param_list$scale_variables_to_regress <- scale_data_variables_to_regress
whole_param_list$sample_variable <- integration_group
whole_param_list$harmony_seed <- seed4Harmony
whole_param_list$metadata_categories <- metadata_categories
whole_param_list$resolution <- whole_resolution
whole_param_list$plots_per_page <- nb_SpatialFeaturePlot_per_page
whole_param_list$markers_test <- "negbinom"
whole_param_list$markers_latent_variables <- "orig.ident"
whole_param_list$markers_cluster <- FALSE
whole_param_list$markers_cluster_condition <- NA
whole_param_list$markers_cluster_plot_nb <- NA
whole_param_list$markers_DE <- TRUE
whole_param_list$markers_DE_timepoint_analysis <- TRUE
whole_param_list$markers_DE_plot_nb <- nb_DEGs_to_plot
whole_param_list$markers_DE_plots_per_page <- nb_DEGs_per_page
whole_param_list$markers_conserved <- FALSE
whole_param_list$markers_condition <- "condition"
print(whole_param_list)

seurat_object_whole_clustering <- norm2markers_workflow(seurat_objects, whole_param_list, whole_clustering_output_dir)
rm(seurat_objects)
gc()

## subclusterings
subclustering_output_dir <- sprintf("%s/sub", integration_output_dir)
if (! dir.exists(subclustering_output_dir)) {
    dir.create(subclustering_output_dir, recursive=TRUE, mode="0775")
}

seurat_object_subclustering_list <- list()
for (one_cluster in names(subsclusterings_parameter_settings)) {
    one_cluster_subclustering_output_dir <- sprintf("%s/c%s", subclustering_output_dir, one_cluster)
    if (! dir.exists(one_cluster_subclustering_output_dir)) {
        dir.create(one_cluster_subclustering_output_dir, recursive=TRUE, mode="0775")
    }
    
    ### set parameters
    sub_param_list <- list()
    sub_param_list$project_name <- output_dirname
    sub_param_list$normalization_method <- subsclusterings_parameter_settings[[one_cluster]][["normalization"]]
    sub_param_list$normalization_variables_to_regress <- sct_variables_to_regress
    sub_param_list$feature_nb <- subsclusterings_parameter_settings[[one_cluster]][["features"]]
    sub_param_list$mitochondrial_genes <- mitochondrial_geneset
    sub_param_list$hemoglobin_genes <- hemoglobin_geneset
    sub_param_list$integration_method <- subsclusterings_parameter_settings[[one_cluster]][["integration"]]
    sub_param_list$dimension_nb <- subsclusterings_parameter_settings[[one_cluster]][["dimensions"]]
    sub_param_list$scale_variables_to_regress <- scale_data_variables_to_regress
    sub_param_list$sample_variable <- integration_group
    sub_param_list$harmony_seed <- seed4Harmony
    sub_param_list$metadata_categories <- metadata_categories
    sub_param_list$resolution <- subsclusterings_parameter_settings[[one_cluster]][["resolution"]]
    sub_param_list$plots_per_page <- nb_SpatialFeaturePlot_per_page
    sub_param_list$markers_test <- "negbinom"
    sub_param_list$markers_latent_variables <- "orig.ident"
    sub_param_list$markers_cluster <- FALSE
    sub_param_list$markers_cluster_condition <- NA
    sub_param_list$markers_cluster_plot_nb <- NA
    sub_param_list$markers_DE <- TRUE
    sub_param_list$markers_DE_timepoint_analysis <- TRUE
    sub_param_list$markers_DE_plot_nb <- nb_DEGs_to_plot
    sub_param_list$markers_DE_plots_per_page <- nb_DEGs_per_page
    sub_param_list$markers_conserved <- FALSE
    sub_param_list$markers_condition <- "condition"
    print(sub_param_list)
    
    ### subset whole clustering Seurat object to the cluster barcodes
    #### get subcluster spots
    barcodes2use <- rownames(seurat_object_whole_clustering@meta.data[which(seurat_object_whole_clustering@meta.data$seurat_clusters == subsclusterings_parameter_settings[[one_cluster]][[sub_param_list$integration_method]]),])
    #### subset
    seurat_object_whole_clustering_sub <- subset(seurat_object_whole_clustering, cells=barcodes2use)
    seurat_object_whole_clustering_sub@project.name <- sprintf("%s_%s", seurat_object_whole_clustering_sub@project.name, one_cluster)
    
    ### split integrated object into a list of objects per sample
    seurat_objects_sub <- SplitObject(seurat_object_whole_clustering_sub, split.by=sub_param_list$sample_variable)
    rm(seurat_object_whole_clustering_sub)
    gc()
    #### after splitting, all the objects still have all the images: only keep the image of the corresponding sample
    for (one_sample in names(seurat_objects_sub)) {
        seurat_objects_sub[[one_sample]]@images <- list(seurat_objects_sub[[one_sample]]@images[[one_sample]])
        names(seurat_objects_sub[[one_sample]]@images) <- one_sample
    }
    if (whole_param_list$normalization_method %in% c("SCTransform", "SCTransform_v2")) {
        #### after splitting, all the objects have 4 models in SCTModel.list slot but only one of them has cells: only keep this one
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
    ### set default assay to RNA and project name to sample
    for (one_sample in names(seurat_objects_sub)) {
        DefaultAssay(seurat_objects_sub[[one_sample]]) <- "RNA"
        seurat_objects_sub[[one_sample]]@project.name <- one_sample
    }

    ### normalization to marker identification
    seurat_object_subclustering_list[[one_cluster]] <- norm2markers_workflow(seurat_objects_sub, sub_param_list, one_cluster_subclustering_output_dir)
    rm(seurat_objects_sub)
    gc()
}


# Seurat object with clusters and subclusters in meta.data
## get subclusters for each cluster
all_clusters_subcluster_df <- data.frame()
for (one_cluster in levels(seurat_object_whole_clustering@meta.data$seurat_clusters)) {
    one_cluster_subcluster_df <- data.frame(barcode=rownames(seurat_object_subclustering_list[[one_cluster]]@meta.data), seurat_subclusters=sprintf("%s.%s", one_cluster, seurat_object_subclustering_list[[one_cluster]]@meta.data$seurat_clusters))
    all_clusters_subcluster_df <- rbind(all_clusters_subcluster_df, one_cluster_subcluster_df)
}
## add subclusters to whole clustering Seurat object
rownames(all_clusters_subcluster_df) <- all_clusters_subcluster_df$barcode
seurat_object_whole_clustering@meta.data$seurat_subclusters <- all_clusters_subcluster_df[rownames(seurat_object_whole_clustering@meta.data), "seurat_subclusters"]
seurat_object_whole_clustering@meta.data$seurat_subclusters <- as.factor(seurat_object_whole_clustering@meta.data$seurat_subclusters)
## add custom clusters to whole clustering Seurat object: subclusters only for clusters 0 and 3
all_clusters_custom_cluster_df <- data.frame()
### clusters 0 and 3: get subclusters
barcodes <- rownames(seurat_object_whole_clustering@meta.data[which(seurat_object_whole_clustering@meta.data$seurat_clusters %in% c("0", "3")),])
all_clusters_custom_cluster_df <- rbind(all_clusters_custom_cluster_df, data.frame(barcode=barcodes, seurat_custom_clusters=seurat_object_whole_clustering@meta.data[barcodes, "seurat_subclusters"]))
### clusters 1, 2 and 4: get clusters
barcodes <- rownames(seurat_object_whole_clustering@meta.data[which(! seurat_object_whole_clustering@meta.data$seurat_clusters %in% c("0", "3")),])
all_clusters_custom_cluster_df <- rbind(all_clusters_custom_cluster_df, data.frame(barcode=barcodes, seurat_custom_clusters=seurat_object_whole_clustering@meta.data[barcodes, "seurat_clusters"]))
all_clusters_custom_cluster_df <- droplevels(all_clusters_custom_cluster_df)
### add custom clusters to whole clustering Seurat object
rownames(all_clusters_custom_cluster_df) <- all_clusters_custom_cluster_df$barcode
seurat_object_whole_clustering@meta.data$seurat_custom_clusters <- all_clusters_custom_cluster_df[rownames(seurat_object_whole_clustering@meta.data), "seurat_custom_clusters"]
seurat_object_whole_clustering@meta.data$seurat_custom_clusters <- as.factor(seurat_object_whole_clustering@meta.data$seurat_custom_clusters)

## save whole clustering Seurat object and list of subclustering Seurat objects
save(seurat_object_whole_clustering, seurat_object_subclustering_list, file=sprintf("%s/whole_clustering_subclustering_Seurat_objects.RData", integration_output_dir))
#load(sprintf("%s/whole_clustering_subclustering_Seurat_objects.RData", integration_output_dir))

## plot clusters, subclusters and custom clusters
pdf(sprintf("%s/whole_clustering_subclusterings.pdf", integration_output_dir))
for (cluster_col in c("seurat_clusters", "seurat_subclusters", "seurat_custom_clusters")) {
    Idents(seurat_object_whole_clustering) <- seurat_object_whole_clustering@meta.data[[cluster_col]]
    # one page per timepoint
    timepoints <- levels(seurat_object_whole_clustering@meta.data$time)
    for (one_timepoint in timepoints) {
        one_timepoint_metadata_df <- seurat_object_whole_clustering@meta.data[which(seurat_object_whole_clustering@meta.data$time==one_timepoint),]
        one_timepoint_metadata_df <- droplevels(one_timepoint_metadata_df)
        sections <- levels(one_timepoint_metadata_df$orig.ident)
        print(SpatialDimPlot(seurat_object_whole_clustering, images=sections, ncol=2, label=TRUE, label.size=2) + plot_annotation(title=sprintf("%s\ntimepoint: %s", cluster_col, one_timepoint), theme=theme(plot.title=element_text(size=18))))
    }
}
dev.off()

## cluster markers
normalization <- whole_param_list$normalization_method
mito_genes <- whole_param_list$mitochondrial_genes
hb_genes <- whole_param_list$hemoglobin_genes
markers_test <- "negbinom"
markers_latent_vars <- "orig.ident"
cluster_markers <- TRUE
cluster_marker_condition <- "CTRL"
sample_variable <- whole_param_list$sample_variable
scale_var2regress <- whole_param_list$scale_variables_to_regress
cluster_markers_plot_nb <- 30
plots_per_page <- whole_param_list$plots_per_page
DE_analysis <- FALSE
timepoint_analysis <- FALSE
nb_DEGs_to_plot <- NA
nb_DEGs_per_page <- NA
conserved_markers <- FALSE
markers_condition <- "condition"

### build output directory path and output basename
normalization_dir_index <- ifelse(normalization == "SCTransform", "00", "01")
vars2regress_in_dirname <- ifelse(is.na(whole_param_list$normalization_variables_to_regress), "no_vars_to_regress", paste(gsub("[.]", "", whole_param_list$normalization_variables_to_regress), collapse="-"))
whole_clustering_resolution_output_dir <- sprintf("%s/%s-%s/%s/featnb%d/%s/%d_dims/res%s", whole_clustering_output_dir, normalization_dir_index, normalization, vars2regress_in_dirname, whole_param_list$feature_nb, gsub(" ", "_", whole_param_list$integration_method), whole_param_list$dimension_nb, sub("[.]", "_", whole_param_list$resolution))
whole_clustering_resolution_output_name <- sprintf("%s_featnb%d_%s_%ddims_res%s", normalization, whole_param_list$feature_nb, gsub(" ", "_", whole_param_list$integration_method), whole_param_list$dimension_nb, sub("[.]", "_", whole_param_list$resolution))

### cluster marker identification analysis
if (normalization %in% c("SCTransform", "SCTransform_v2")) {
    default_assay <- "SCT"
} else {
    default_assay <- "RNA"
}
marker_analysis(seurat_object_whole_clustering, default_assay, mito_genes, hb_genes, markers_test, markers_latent_vars, cluster_markers, cluster_marker_condition, sample_variable, scale_var2regress, TRUE, cluster_markers_plot_nb, plots_per_page, DE_analysis, timepoint_analysis, nb_DEGs_to_plot, nb_DEGs_per_page, conserved_markers, "seurat_subclusters", markers_condition, whole_clustering_resolution_output_dir, whole_clustering_resolution_output_name)

