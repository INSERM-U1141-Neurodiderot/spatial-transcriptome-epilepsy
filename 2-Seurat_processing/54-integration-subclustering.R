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
integration_group <- "orig.ident"
integration_group_by_in_path <- "orig"
seed4Harmony <- 1
#LIGER_integration <- FALSE
cluster_markers <- FALSE
DE_analysis <- FALSE
DE_timepoint_analysis <- FALSE
conserved_markers <- FALSE

# number of samples per page when using SpatialFeaturePlot() function
nb_SpatialFeaturePlot_per_page <- 4
# number of samples per page when using SpatialFeaturePlot() function
nb_cluster_markers_to_plot <- 4

nb_DEGs_to_plot <- 30
nb_DEGs_per_page <- 3

metadata_categories <- c("lame", "zone", "condition", "time", "orig.ident")

# best whole clustering integration method
whole_clustering_integration <- "Seurat"

# subclustering paramets
subclustering_normalization <- "SCTransform_v2"
subclustering_features_param_vector <- c(250, 500, 1000, 1500, 2000)
subclustering_integration_methods <- c("Harmony", "Seurat")
subclustering_dimension_vector <- c(10, 15, 20)
subclustering_resolution_vector <- c(0.01, 0.03, 0.05, 0.08, 0.1, 0.3, 0.5)


##############
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
output_dir <- sprintf("%s/20-Integration/%s/sub", output_dir, output_dirname)
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

# get whole clustering parameter evaluation output
if (visium_params_in_path) {
    whole_clustering_output_dir <- sprintf("%s/20-Count_analysis/10-EpiReg/output/%s/TSO_polyA_R1trim1_ov5_n2_min20/%s/Visium_params/20-Integration/%s/whole/Norm_%s-scale_%s-int_%s/QC-sum%d_det%s_mito%d_hb%d/10-QC_filtering", work_dir, input_dirname, genome_name, output_dirname, variables_to_regress_in_path, scale_data_split_by_in_path, integration_group_by_in_path, spot_sum_of_counts_threshold, spot_detected_genes_threshold, spot_mito_gene_count_pct_threshold, spot_hb_gene_count_pct_threshold)
} else {
    whole_clustering_output_dir <- sprintf("%s/20-Count_analysis/10-EpiReg/output/%s/TSO_polyA_R1trim1_ov5_n2_min20/%s/20-Integration/%s/whole/Norm_%s-scale_%s-int_%s/QC-sum%d_det%s_mito%d_hb%d/10-QC_filtering", work_dir, input_dirname, genome_name, output_dirname, variables_to_regress_in_path, scale_data_split_by_in_path, integration_group_by_in_path, spot_sum_of_counts_threshold, spot_detected_genes_threshold, spot_mito_gene_count_pct_threshold, spot_hb_gene_count_pct_threshold)
}

whole_clustering_eval_df <- read.csv(sprintf("%s/%s_norm_featnb_int_dims_res_cluster_mean_silhouette.csv", whole_clustering_output_dir, output_dirname))

# print best parameter setting per integration method
whole_params <- clustering_top_params(whole_clustering_eval_df, "mean_silhouette_seurat_clusters", whole_clustering_integration)
whole_normalization <- whole_params[["normalization"]]
whole_features <- whole_params[["features"]]
whole_integration <- whole_params[["integration"]]
whole_dimensions <- whole_params[["dimensions"]]
whole_resolution <- whole_params[["resolution"]]
print(sprintf("normalization: %s, feature number: %d, integration: %s, dimension number: %d, resolution: %.2f", whole_normalization, whole_features, whole_integration, whole_dimensions, whole_resolution))

# whole clustering
whole_clustering_output_dir <- sprintf("%s/%s/whole", qc_filtering_output_dir, whole_integration)
if (! dir.exists(whole_clustering_output_dir)) {
    dir.create(whole_clustering_output_dir, recursive=TRUE, mode="0775")
}

## set parameters
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
whole_param_list$markers_test <- NA
whole_param_list$markers_latent_variables <- NA
whole_param_list$markers_cluster <- cluster_markers
whole_param_list$markers_cluster_plot_nb <- nb_cluster_markers_to_plot
whole_param_list$markers_DE <- DE_analysis
whole_param_list$markers_DE_timepoint_analysis <- DE_timepoint_analysis
whole_param_list$markers_DE_plot_nb <- nb_DEGs_to_plot
whole_param_list$markers_DE_plots_per_page <- nb_DEGs_per_page
whole_param_list$markers_conserved <- conserved_markers
whole_param_list$markers_condition <- "condition"

seurat_object_whole_clustering <- norm2markers_workflow(seurat_objects, whole_param_list, whole_clustering_output_dir)

# subclusterings
subclustering_batch_silhouette_df <- subclustering_cluster_silhouette_df <- data.frame()
for (subclustering_integration_method in subclustering_integration_methods) {
    subclustering_output_dir <- sprintf("%s/%s/sub", qc_filtering_output_dir, subclustering_integration_method)
    if (! dir.exists(subclustering_output_dir)) {
        dir.create(subclustering_output_dir, recursive=TRUE, mode="0775")
    }
    for (one_cluster in levels(seurat_object_whole_clustering@meta.data$seurat_clusters)) {
        one_cluster_subclustering_output_dir <- sprintf("%s/c%s", subclustering_output_dir, one_cluster)
        if (! dir.exists(one_cluster_subclustering_output_dir)) {
            dir.create(one_cluster_subclustering_output_dir, recursive=TRUE, mode="0775")
        }

        ## subset whole clustering Seurat object to the cluster barcodes
        ### get subcluster spots
        barcodes2use <- rownames(seurat_object_whole_clustering@meta.data[which(seurat_object_whole_clustering@meta.data$seurat_clusters == one_cluster),])
        ### subset
        seurat_object_whole_clustering_sub <- subset(seurat_object_whole_clustering, cells=barcodes2use)
        seurat_object_whole_clustering_sub@project.name <- sprintf("%s_%s", seurat_object_whole_clustering_sub@project.name, one_cluster)

        ## approach 1: only re-run dimensionality reduction and clustering using subset integrated assay
        ### input: seurat_object_integration
        subclustering_approach <- "dim"
        subclustering_approach_output_dir <- sprintf("%s/%s", one_cluster_subclustering_output_dir, subclustering_approach)
        if (! dir.exists(subclustering_approach_output_dir)) {
            dir.create(subclustering_approach_output_dir, recursive=TRUE, mode="0775")
        }
        subclustering_approach_output_name <- sprintf("%s_featnb%d_%s_sub", whole_normalization, whole_features, subclustering_integration_method)

        if (subclustering_integration_method == "Seurat") {
            seurat_object_integration_dim <- RunPCA(seurat_object_whole_clustering_sub, verbose=TRUE)
            subclustering_dimension_vector_approach <- subclustering_dimension_vector
            reduction2use <- "pca"
        } else {
            seurat_object_integration_dim <- seurat_object_whole_clustering_sub
            subclustering_dimension_vector_approach <- subclustering_dimension_vector[subclustering_dimension_vector <= whole_dimensions]
            if (subclustering_integration_method == "Harmony") {
                reduction2use <- "harmony"
            } else {
                reduction2use <- "pca"
            }
        }

        for (subclustering_dimension_nb in subclustering_dimension_vector_approach) {
            dimensions_output_dir <- sprintf("%s/%d_dims", subclustering_approach_output_dir, subclustering_dimension_nb)
            if (! dir.exists(dimensions_output_dir)) {
                dir.create(dimensions_output_dir, recursive=TRUE, mode="0775")
            }
            dimensions_output_name <- sprintf("%s_%d_dims", subclustering_approach_output_name, subclustering_dimension_nb)
            sub_silhouette_df <- graph_clustering(seurat_object_integration_dim, subclustering_dimension_nb, reduction2use, subclustering_resolution_vector, metadata_categories, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integration: %s, subclustering approach: %s", subclustering_integration_method, subclustering_approach), dimensions_output_dir, dimensions_output_name)
            subclustering_cluster_silhouette_df <- rbind(subclustering_cluster_silhouette_df, cbind(data.frame(whole_cluster=rep(one_cluster, dim(sub_silhouette_df)[1]), approach=rep(subclustering_approach, dim(sub_silhouette_df)[1]), normalization=rep(whole_normalization, dim(sub_silhouette_df)[1]), features=whole_features, integration=rep(subclustering_integration_method, dim(sub_silhouette_df)[1])), sub_silhouette_df))
        }
        rm(seurat_object_integration_dim)
        gc()

        ## approach 2: re-run SCTransform to get HVGs, no integration
        ### input: seurat_objects for SCTransform, seurat_object_integration for dimensionality reduction and subclustering
        subclustering_approach <- "hvg_dim"
        subclustering_approach_output_dir <- sprintf("%s/%s", one_cluster_subclustering_output_dir, subclustering_approach)
        if (! dir.exists(subclustering_approach_output_dir)) {
            dir.create(subclustering_approach_output_dir, recursive=TRUE, mode="0775")
        }

        if (subclustering_integration_method == "Harmony") {
            subclustering_dimension_vector_approach <- subclustering_dimension_vector[subclustering_dimension_vector <= whole_dimensions]    
        } else {
            subclustering_dimension_vector_approach <- subclustering_dimension_vector
        }

        return_only_var_genes <- ifelse(subclustering_normalization %in% c("SCTransform", "SCTransform_v2"), FALSE, NA)
        silhouette_df_list <- norm_dimred_clustering(seurat_object_whole_clustering_sub, subclustering_normalization, subclustering_features_param_vector, sct_variables_to_regress, return_only_var_genes, mitochondrial_geneset, hemoglobin_geneset, FALSE, subclustering_integration_method, NA, integration_group, subclustering_dimension_vector_approach, metadata_categories, subclustering_resolution_vector, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integration: %s, subclustering approach: %s", subclustering_integration_method, subclustering_approach), subclustering_approach_output_dir)
        subclustering_cluster_silhouette_df <- rbind(subclustering_cluster_silhouette_df, cbind(whole_cluster=rep(one_cluster, dim(silhouette_df_list[["cluster"]])[1]), approach=rep(subclustering_approach, dim(silhouette_df_list[["cluster"]])[1]), silhouette_df_list[["cluster"]]))

        ## approach 3: re-run SCTransform to get HVGs, integration
        ### input: seurat_objects for SCTransform, seurat_object_integration for dimensionality reduction and subclustering
        subclustering_approach <- "hvg_dim_int"
        subclustering_approach_output_dir <- sprintf("%s/%s", one_cluster_subclustering_output_dir, subclustering_approach)
        if (! dir.exists(subclustering_approach_output_dir)) {
            dir.create(subclustering_approach_output_dir, recursive=TRUE, mode="0775")
        }

        return_only_var_genes <- ifelse(subclustering_normalization %in% c("SCTransform", "SCTransform_v2"), FALSE, NA)
        silhouette_df_list <- norm_dimred_clustering(seurat_object_whole_clustering_sub, subclustering_normalization, subclustering_features_param_vector, sct_variables_to_regress, return_only_var_genes, mitochondrial_geneset, hemoglobin_geneset, TRUE, subclustering_integration_method, seed4Harmony, integration_group, subclustering_dimension_vector, metadata_categories, subclustering_resolution_vector, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integration: %s, subclustering approach: %s", subclustering_integration_method, subclustering_approach), subclustering_approach_output_dir)
        subclustering_batch_silhouette_df <- rbind(subclustering_batch_silhouette_df, cbind(whole_cluster=rep(one_cluster, dim(silhouette_df_list[["batch"]])[1]), approach=rep(subclustering_approach, dim(silhouette_df_list[["batch"]])[1]), silhouette_df_list[["batch"]]))
        subclustering_cluster_silhouette_df <- rbind(subclustering_cluster_silhouette_df, cbind(whole_cluster=rep(one_cluster, dim(silhouette_df_list[["cluster"]])[1]), approach=rep(subclustering_approach, dim(silhouette_df_list[["cluster"]])[1]), silhouette_df_list[["cluster"]]))
    }
}
rm(seurat_object_whole_clustering)
rm(seurat_objects)
gc()

# mean silhouette plots for all parameter settings
## clustering
subclustering_cluster_silhouette_df$whole_cluster <- as.factor(subclustering_cluster_silhouette_df$whole_cluster)
subclustering_cluster_silhouette_df$approach <- as.factor(subclustering_cluster_silhouette_df$approach)
subclustering_cluster_silhouette_df$integration <- as.factor(subclustering_cluster_silhouette_df$integration)
subclustering_cluster_silhouette_df$features <- factor(subclustering_cluster_silhouette_df$features, levels=subclustering_features_param_vector)
subclustering_cluster_silhouette_df$dimensions <- as.factor(subclustering_cluster_silhouette_df$dimensions)
subclustering_cluster_silhouette_df$resolution <- as.factor(subclustering_cluster_silhouette_df$resolution)
plot_output_name <- sprintf("%s_subclustering_norm_featnb_int_dims_res_cluster_silhouette", output_dirname)
write.csv(subclustering_cluster_silhouette_df, file=sprintf("%s/%s.csv", qc_filtering_output_dir, plot_output_name), quote=FALSE, row.names=FALSE)
for (one_cluster in levels(subclustering_cluster_silhouette_df$whole_cluster)) {
    one_cluster_df <- subclustering_cluster_silhouette_df[which(subclustering_cluster_silhouette_df$whole_cluster == one_cluster),]
    plot_output_name <- sprintf("%s_%s_subclustering_norm_featnb_int_dims_res_cluster_silhouette", output_dirname, one_cluster)
    pdf(sprintf("%s/%s.pdf", qc_filtering_output_dir, plot_output_name))
    for (one_approach in levels(one_cluster_df$approach)) {
        df2ggplot <- one_cluster_df[which(one_cluster_df$approach == one_approach),]
        p <- ggplot(df2ggplot, aes(x=resolution, y=mean_silhouette_seurat_clusters_sub, fill=integration)) +
            geom_bar(stat="identity", position="dodge") +
            geom_text(aes(label=nb_clusters_sub), position=position_dodge(width=1), size=2.5, vjust=2) +
            facet_grid(dimensions~features) +
            labs(title=sprintf("%s subclustering\nmean silhouette: Seurat clusters\nsubclustering approach: %s", one_cluster, one_approach), x="Resolution", y="Mean silhouette", fill="Integration method") +
            theme_bw() +
            theme(legend.position="bottom") +
            theme(axis.text.x=element_text(angle=60, hjust=1)) +
            theme(panel.border=element_rect(color="grey50"))
        print(p)
    }
    dev.off()
}

## batch
subclustering_batch_silhouette_df$whole_cluster <- as.factor(subclustering_batch_silhouette_df$whole_cluster)
subclustering_batch_silhouette_df$approach <- as.factor(subclustering_batch_silhouette_df$approach)
subclustering_batch_silhouette_df$integration <- as.factor(subclustering_batch_silhouette_df$integration)
subclustering_batch_silhouette_df$features <- factor(subclustering_batch_silhouette_df$features, levels=subclustering_features_param_vector)
subclustering_batch_silhouette_df$dimensions <- as.factor(subclustering_batch_silhouette_df$dimensions)
plot_output_name <- sprintf("%s_subclustering_norm_featnb_int_dims_res_batch_silhouette", output_dirname)
write.csv(subclustering_batch_silhouette_df, file=sprintf("%s/%s.csv", qc_filtering_output_dir, plot_output_name), quote=FALSE, row.names=FALSE)
for (one_cluster in levels(subclustering_batch_silhouette_df$whole_cluster)) {
    df2ggplot <- subclustering_batch_silhouette_df[which(subclustering_batch_silhouette_df$whole_cluster==one_cluster),]
    plot_output_name <- sprintf("%s_%s_subclustering_norm_featnb_int_dims_res_batch_silhouette", output_dirname, one_cluster)
    pdf(sprintf("%s/%s.pdf", qc_filtering_output_dir, plot_output_name))
    ### mean silhouette
    for (one_category in metadata_categories) {
        p <- ggplot(df2ggplot, aes(x=dimensions, y=.data[[sprintf("mean_silhouette_%s", one_category)]], fill=integration)) +
            geom_bar(stat="identity", position="dodge") +
            facet_grid(~features) +
            labs(title=sprintf("%s subclustering\nmean silhouette: %s", one_cluster, one_category), x="Number of dimensions", y="Mean silhouette", fill="Integration method") +
            theme_bw() +
            theme(legend.position="bottom") +
            theme(panel.border=element_rect(color="grey50"))
        print(p)
    }
    for (one_category in metadata_categories) {
        p <- ggplot(df2ggplot, aes(x=integration, y=.data[[sprintf("mean_silhouette_%s", one_category)]], fill=integration)) +
            geom_bar(stat="identity", position="dodge") +
            facet_grid(dimensions~features) +
            labs(title=sprintf("%s subclustering\nmean silhouette: %s", one_cluster, one_category), x="Integration method", y="Mean silhouette", fill="Integration method") +
            theme_bw() +
            theme(legend.position="bottom") +
            theme(axis.text.x=element_text(angle=60, hjust=1)) +
            theme(panel.border=element_rect(color="grey50"))
        print(p)
    }
    dev.off()    
}

