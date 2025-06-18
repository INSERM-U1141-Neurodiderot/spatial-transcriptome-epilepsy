.libPaths(c(.libPaths(), "/home/christophe.lepriol/NeuroDev_ADD/R/r_4.1.0"))
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)


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
#features_param_vector <- c(500, 1000, 1500, 2000, 2500, 3000)
#features_param_vector <- c(500, 1000, 1500, 2000, 2500)
features_param_vector <- c(250, 500, 1000, 1500, 2000)
#features_param_vector <- c(2000)
#dimensions_nb_param_vector <- c(15, 20, 30)
dimensions_nb_param_vector <- c(10, 15, 20)
#dimensions_nb_param_vector <- c(15, 20)
#clustering_resolution_vector <- c(0.1, 0.3, 0.5, 0.8, 1, 1.2)
clustering_resolution_vector <- c(0.01, 0.03, 0.05, 0.08, 0.1, 0.3, 0.5)
#clustering_resolution_vector <- c(0.05)
seed4Harmony <- 1
LIGER_integration <- FALSE
cluster_markers <- FALSE
cluster_marker_condition <- NA
DE_analysis <- FALSE
timepoint_analysis <- FALSE
conserved_markers <- FALSE

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
output_dir <- sprintf("%s/20-Integration/%s/whole", output_dir, output_dirname)
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
lapply(seurat_objects, seurat_all_feature_stat_plots, assay2use="RNA", col_prefix=qc_column_prefix, plot_title="Raw counts after QC", out_dir=qc_filtering_output_dir, out_name="QC_stats_filtering") 
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

# no integration, Seurat and Harmony
integration_methods <- c("no integration", "Harmony", "Seurat")
library(harmony)
library(cluster) # silhouette() function
#library(kBET)
integration_batch_mean_silhouette_df <- data.frame()
#integration_kBET_df <- data.frame()
integration_cluster_mean_silhouette_df <- data.frame()
for (normalization_method in normalization_method_vector) {
    print(sprintf("normalization method: %s", normalization_method))
    return_only_var_genes <- ifelse(normalization_method %in% c("SCTransform", "SCTransform_v2"), FALSE, NA)
    for (features_param in features_param_vector) {
        returned_list <- seurat_list_normalization(seurat_objects, normalization_method, sct_variables_to_regress, features_param, return_only_var_genes, mitochondrial_geneset, hemoglobin_geneset, qc_filtering_output_dir)
        seurat_objects_norm <- returned_list[["seurat_list"]]
        integration_features_nb <- returned_list[["features_nb"]]
        integration_normalization_method <- returned_list[["norm_method"]]
        default_assay <- returned_list[["assay"]]
        selection_method <- returned_list[["selec"]]
        normalization_output_dir <- returned_list[["out_dir"]]
        normalization_output_name <- returned_list[["out_name"]]
        rm(returned_list)
        gc()
        
        if (! LIGER_integration & normalization_method == normalization_method_vector[length(normalization_method_vector)] & features_param == features_param_vector[length(features_param_vector)]) {
            # last normalization: Seurat objects before normalization are not needed any more
            rm(seurat_objects)
            gc()
        }

        ## HVG plots
        hvg_df_list <- lapply(seurat_objects_norm, hvg_plots, assay2use=default_assay, selec=selection_method, nb_hvgs=features_param, plots_per_page=NA, out_dir=normalization_output_dir, out_name=sprintf("%s_HVGs", normalization_output_name))

        ## merge datasets, feature statistics and QC plots after normalization
        seurat_object_merge_norm <- merge_stats_plots_after_norm(seurat_objects_norm, output_dirname, normalization_method, features_param, mitochondrial_geneset, hemoglobin_geneset, nb_SpatialFeaturePlot_per_page, normalization_output_dir, normalization_output_name)
        
        ## integration
        for (integration_method in integration_methods) {
            integration_output_dir <- sprintf("%s/%s", normalization_output_dir, gsub(" ", "_", integration_method))
            if (! dir.exists(integration_output_dir)) {
                dir.create(integration_output_dir, recursive=TRUE, mode="0775")
            }
            ### dimensionality reduction
            for (dimensions_nb_param in dimensions_nb_param_vector) {
                print(sprintf("dimension number: %d", dimensions_nb_param))
                dimensions_output_dir <- sprintf("%s/%d_dims", integration_output_dir, dimensions_nb_param)
                if (! dir.exists(dimensions_output_dir)) {
                    dir.create(dimensions_output_dir, recursive=TRUE, mode="0775")
                }
                dimensions_output_name <- sprintf("%s_%s_%ddims", normalization_output_name, gsub(" ", "_", integration_method), dimensions_nb_param)
                
                if (integration_method == "Seurat") {
                    ### Seurat integration
                    features <- SelectIntegrationFeatures(seurat_objects_norm, nfeatures=integration_features_nb)
                    if (normalization_method %in% c("SCTransform", "SCTransform_v2")) {
                        seurat_objects_norm_prep <- PrepSCTIntegration(object.list=seurat_objects_norm, anchor.features=features)
                    } else {
                        seurat_objects_norm_prep <- seurat_objects_norm
                    }
                    if (dimensions_nb_param == dimensions_nb_param_vector[length(dimensions_nb_param_vector)]) {
                        # last Seurat integration: Seurat objects after normalization are not needed any more
                        rm(seurat_objects_norm)
                        gc()
                    }
                    int.anchors <- FindIntegrationAnchors(object.list=seurat_objects_norm_prep, anchor.features=features, normalization.method=integration_normalization_method, dims=1:dimensions_nb_param)
                    rm(seurat_objects_norm_prep)
                    gc()
                    seurat_object_integration <- IntegrateData(anchorset=int.anchors, normalization.method=integration_normalization_method, dims=1:dimensions_nb_param)
                    rm(int.anchors)
                    gc()
                    #### set project name to output dirname
                    seurat_object_integration@project.name <- output_dirname
                    #### orig.ident, lame, zone, condition and time are not factors in the integrated object meta data
                    seurat_object_integration@meta.data$orig.ident <- as.factor(seurat_object_integration@meta.data$orig.ident)
                    seurat_object_integration@meta.data$lame <- as.factor(seurat_object_integration@meta.data$lame)
                    seurat_object_integration@meta.data$zone <- as.factor(seurat_object_integration@meta.data$zone)
                    seurat_object_integration@meta.data$condition <- as.factor(seurat_object_integration@meta.data$condition)
                    seurat_object_integration@meta.data$time <- as.factor(seurat_object_integration@meta.data$time)
                    
                    #### HVGs after integration
                    seurat_object_integration_hvgs <- VariableFeatures(seurat_object_integration)
                    seurat_object_integration_df <- seurat_object_integration@assays$integrated@SCTModel.list$refmodel@feature.attributes[seurat_object_integration_hvgs,]
                    write.csv(seurat_object_integration_df, file=sprintf("%s/%s_%s_HVGs_residual_mean_variance_after_integration.csv", dimensions_output_dir, seurat_object_integration@project.name, dimensions_output_name), quote=FALSE, row.names=TRUE)
                } else {
                    seurat_object_integration <- seurat_object_merge_norm
                    if (integration_method == "Harmony" & dimensions_nb_param == dimensions_nb_param_vector[length(dimensions_nb_param_vector)]) {
                        # last Harmony integration: merged Seurat object after normalization are not needed any more
                        rm(seurat_object_merge_norm)
                        gc()
                    }
                }

                ### scaling
                if (normalization_method == "LogNormalize") {
                    if (! is.na(scale_data_split_by)) {
                        if (!is.na(scale_data_variables_to_regress)) {
                            seurat_object_integration <- ScaleData(seurat_object_integration, split.by=scale_data_split_by, vars.to.regress=scale_data_variables_to_regress)
                            dimensions_output_dir <- sprintf("%s/%s_%s", dimensions_output_dir, gsub("[.]", "", scale_data_split_by), paste(gsub("[.]", "", scale_data_variables_to_regress), collapse="-"))
                        } else {
                            seurat_object_integration <- ScaleData(seurat_object_integration, split.by=scale_data_split_by)
                            dimensions_output_dir <- sprintf("%s/%s", dimensions_output_dir, gsub("[.]", "", scale_data_split_by))
                        }
                    } else {
                        seurat_object_integration <- ScaleData(seurat_object_integration)
                    }
                }
                
                seurat_object_integration <- RunPCA(seurat_object_integration, verbose=FALSE)
                
                #### Harmony integration and dimensionality reduction
                if (integration_method == "Harmony") {
                    ##### dims.use parameter of RunHarmony() function is not used: build a new DimReduc object with the desired number of dimensions
                    ###### https://github.com/immunogenomics/harmony/issues/82
                    ###### https://github.com/immunogenomics/harmony/issues/151
                    seurat_object_integration@reductions$pca_dims <- seurat_object_integration@reductions$pca
                    seurat_object_integration@reductions$pca_dims@cell.embeddings <- seurat_object_integration@reductions$pca_dims@cell.embeddings[, 1:dimensions_nb_param]
                    ##### Harmony results are not reproducible: set seed
                    ###### https://github.com/immunogenomics/harmony/issues/13
                    ###### https://github.com/immunogenomics/harmony/issues/56    
                    set.seed(seed4Harmony)
                    seurat_object_integration <- RunHarmony(seurat_object_integration, group.by.vars=integration_group, reduction="pca_dims")
                    reduction2use <- "harmony"
                } else {
                    reduction2use <- "pca"
                }
                
                ### compute silhouette scores for lames, zones, conditions and samples
                mean_silhouette_vector <- c()
                for (metadata_category_colname in metadata_categories) {
                    seurat_object_integration <- compute_silhouette(seurat_object_integration, reduction2use, dimensions_nb_param, metadata_category_colname)
                    mean_silhouette <- mean(seurat_object_integration@meta.data[[sprintf("silhouette_%s", metadata_category_colname)]])
                    mean_silhouette_vector <- c(mean_silhouette_vector, mean_silhouette)
                }
                mean_silhouette_df <- as.data.frame(t(mean_silhouette_vector))
                colnames(mean_silhouette_df) <- sprintf("mean_silhouette_%s", metadata_categories)
                params_mean_silhouette_df <- cbind(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_param, mean_silhouette_df)
                integration_batch_mean_silhouette_df <- rbind(integration_batch_mean_silhouette_df, params_mean_silhouette_df)
                
                ### compute kBET
                #if (reduction2use == "pca") {
                #    kBET_pca_param <- FALSE
                #    kBET_dims_param <- NULL
                #} else {
                #    kBET_pca_param <- TRUE
                #    kBET_dims_param <- dimensions_nb_param
                #}
                #kBET_results_df <- compute_kBET(Embeddings(seurat_object_integration@reductions[[reduction2use]]), seurat_object_integration@meta.data$orig.ident, kBET_pca_param, kBET_dims_param)
                #mean_median_kBET_df <- kBET_plots(kBET_results_df, dimensions_output_dir, sprintf("%s_kBET", dimensions_output_name))
                #nb_pct_sample_size <- dim(mean_median_kBET_df)[1]
                #integration_kBET_df <- rbind(integration_kBET_df, cbind(data.frame(normalization=rep(normalization_method, nb_pct_sample_size), features=rep(features_param, nb_pct_sample_size), integration=rep(integration_method, nb_pct_sample_size), dimensions=rep(dimensions_nb_param, nb_pct_sample_size), cluster=rep("all", nb_pct_sample_size)), mean_median_kBET_df))
                
                ## clustering
                seurat_object_integration <- FindNeighbors(seurat_object_integration, reduction=reduction2use, dims=1:dimensions_nb_param)
                for (clustering_resolution in clustering_resolution_vector) {
                    print(sprintf("cluster resolution: %.1f", clustering_resolution))
                    resolution_output_dir <- sprintf("%s/res%s", dimensions_output_dir, sub("[.]", "_", clustering_resolution))
                    if (! dir.exists(resolution_output_dir)) {
                        dir.create(resolution_output_dir, recursive=TRUE, mode="0775")
                    }
                    resolution_output_name <- sprintf("%s_res%s", dimensions_output_name, sub("[.]", "_", clustering_resolution))

                    seurat_object_clustering <- clustering(seurat_object_integration, clustering_resolution, reduction2use, dimensions_nb_param)
                    write.csv(seurat_object_clustering@meta.data, file=sprintf("%s/%s_metadata.csv", resolution_output_dir, resolution_output_name), quote=FALSE, row.names=TRUE)
                    ### plot clusters onto UMAP or onto the tissue section
                    clustering_plots(seurat_object_clustering, metadata_categories, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integrated method: %s", integration_method), resolution_output_dir, sprintf("%s_clustering", resolution_output_name))

                    ### get mean silhouette coefficients
                    cluster_mean_silhouette <- mean(seurat_object_clustering$silhouette_seurat_clusters)
                    cluster_nb <- length(levels(seurat_object_clustering@meta.data$seurat_clusters))
                    params_mean_silhouette_df <- data.frame(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_param, resolution=clustering_resolution, nb_clusters=cluster_nb, mean_silhouette_seurat_clusters=cluster_mean_silhouette)
                    integration_cluster_mean_silhouette_df <- rbind(integration_cluster_mean_silhouette_df, params_mean_silhouette_df)
                    
                    ### kBET
                    #if (reduction2use == "pca") {
                    #    kBET_pca_param <- FALSE
                    #    kBET_dims_param <- NULL
                    #} else {
                    #    kBET_pca_param <- TRUE
                    #    kBET_dims_param <- dimensions_nb_param
                    #}
                    #mean_median_cluster_kBET_df <- clustering_kBET(seurat_object_clustering, resolution_output_dir, resolution_output_name)
                    #nb_pct_sample_size <- dim(mean_cluster_kBET_df2ggplot)[1]
                    #integration_kBET_df <- rbind(integration_kBET_df, cbind(data.frame(normalization=rep(normalization_method, nb_pct_sample_size), features=rep(features_param, nb_pct_sample_size), integration=rep(integration_method, nb_pct_sample_size), dimensions=rep(dimensions_nb_param, nb_pct_sample_size)), mean_median_cluster_kBET_df))
                    
                    ### cluster markers, differentially expressed genes across conditions and conserved cell type markers
                    marker_analysis(seurat_object_clustering, default_assay, mitochondrial_geneset, hemoglobin_geneset, "negbinom", "orig.ident", cluster_markers, cluster_marker_condition, scale_data_split_by, scale_data_variables_to_regress, load_image, nb_cluster_markers_to_plot, nb_SpatialFeaturePlot_per_page, DE_analysis, timepoint_analysis, nb_DEGs_to_plot, nb_DEGs_per_page, conserved_markers, "seurat_clusters", "condition", resolution_output_dir, resolution_output_name)
                    
                    rm(seurat_object_clustering)
                    gc()
                }
                rm(seurat_object_integration)
                gc()
            }
        }
    }
}
detach("package:harmony", unload=TRUE)

# LIGER
if (LIGER_integration) {
    library(rliger)
    library(SeuratWrappers)
    integration_method <- "LIGER"
    integration_methods <- c(integration_methods, integration_method)
    for (normalization_method in normalization_method_vector) {
        print(sprintf("normalization method: %s", normalization_method))
        return_only_var_genes <- ifelse(normalization_method %in% c("SCTransform", "SCTransform_v2"), TRUE, NA) # SCTransform: only return variable genes, otherwise error with RunOptimizeALS() function: Error in dimnames(x) <- dn: length of 'dimnames' [2] not equal to array extent
        for (features_param in features_param_vector) {
            ## normalization and feature selection
            returned_list <- seurat_list_normalization(seurat_objects, normalization_method, sct_variables_to_regress, features_param, return_only_var_genes, mitochondrial_geneset, hemoglobin_geneset, qc_filtering_output_dir)
            seurat_objects_norm <- returned_list[["seurat_list"]]
            integration_features_nb <- returned_list[["features_nb"]]
            integration_normalization_method <- returned_list[["norm_method"]]
            default_assay <- returned_list[["assay"]]
            selection_method <- returned_list[["selec"]]
            normalization_output_dir <- returned_list[["out_dir"]]
            if (normalization_method %in% c("SCTransform", "SCTransform_v2")) {
                normalization_output_name <- sprintf("%s_onlyvar", returned_list[["out_name"]])
            } else {
                if (normalization_method == "LogNormalize") {
                    normalization_output_name <- sprintf("%s_forLIGER", returned_list[["out_name"]]) 
                }
            }
            rm(returned_list)
            gc()

            ### HVG plots
            hvg_df_list <- lapply(seurat_objects_norm, hvg_plots, assay2use=default_assay, selec=selection_method, nb_hvgs=features_param, plots_per_page=NA, out_dir=normalization_output_dir, out_name=sprintf("%s_HVGs", normalization_output_name))

            ### merge datasets, feature statistics and QC plots after normalization
            seurat_object_merge_norm <- merge_stats_plots_after_norm(seurat_objects_norm, output_dirname, normalization_method, features_param, mitochondrial_geneset, hemoglobin_geneset, nb_SpatialFeaturePlot_per_page, normalization_output_dir, normalization_output_name)
            rm(seurat_objects_norm)
            gc()

            integration_output_dir <- sprintf("%s/%s", normalization_output_dir, gsub(" ", "_", integration_method))
            if (! dir.exists(integration_output_dir)) {
                dir.create(integration_output_dir, recursive=TRUE, mode="0775")
            }

            ### scaling
            if (normalization_method == "LogNormalize") {
                if (! is.na(scale_data_split_by)) {
                    if (!is.na(scale_data_variables_to_regress)) {
                        seurat_object_merge_norm <- ScaleData(seurat_object_merge_norm, split.by=scale_data_split_by, vars.to.regress=scale_data_variables_to_regress, do.center=FALSE) # LIGER does not center data when scaling, skip that step
                        integration_output_dir <- sprintf("%s/%s_%s", integration_output_dir, gsub("[.]", "", scale_data_split_by), paste(gsub("[.]", "", scale_data_variables_to_regress), collapse="-"))
                    } else {
                        seurat_object_merge_norm <- ScaleData(seurat_object_merge_norm, split.by=scale_data_split_by, do.center=FALSE) # LIGER does not center data when scaling, skip that step
                        integration_output_dir <- sprintf("%s/%s", integration_output_dir, gsub("[.]", "", scale_data_split_by))
                    }
                } else {
                    seurat_object_merge_norm <- ScaleData(seurat_object_merge_norm, do.center=FALSE) # LIGER does not center data when scaling, skip that step
                }
            }

            ### integration and dimensionality reduction
            for (dimensions_nb_param in dimensions_nb_param_vector) {
                print(sprintf("dimension number: %d", dimensions_nb_param))
                dimensions_output_dir <- sprintf("%s/%d_dims", integration_output_dir, dimensions_nb_param)
                if (! dir.exists(dimensions_output_dir)) {
                    dir.create(dimensions_output_dir, recursive=TRUE, mode="0775")
                }
                dimensions_output_name <- sprintf("%s_%s_%ddims", normalization_output_name, gsub(" ", "_", integration_method), dimensions_nb_param)

                ### integration and dimensionality reduction
                seurat_object_integration <- RunOptimizeALS(seurat_object_merge_norm, k=dimensions_nb_param, lambda=5, split.by=integration_group)
                seurat_object_integration <- RunQuantileNorm(seurat_object_integration, split.by=integration_group)
                reduction2use <- "iNMF"

                seurat_object_integration <- FindNeighbors(seurat_object_integration, reduction=reduction2use, dims=1:dimensions_nb_param)

                ### compute silhouette scores for lames, zones, conditions and samples
                mean_silhouette_vector <- c()
                for (metadata_category_colname in metadata_categories) {
                    seurat_object_integration <- compute_silhouette(seurat_object_integration, reduction2use, dimensions_nb_param, metadata_category_colname)
                    mean_silhouette <- mean(seurat_object_integration@meta.data[[sprintf("silhouette_%s", metadata_category_colname)]])
                    mean_silhouette_vector <- c(mean_silhouette_vector, mean_silhouette)
                }
                mean_silhouette_df <- as.data.frame(t(mean_silhouette_vector))
                colnames(mean_silhouette_df) <- sprintf("mean_silhouette_%s", metadata_categories)
                params_mean_silhouette_df <- cbind(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_param, mean_silhouette_df)
                integration_batch_mean_silhouette_df <- rbind(integration_batch_mean_silhouette_df, params_mean_silhouette_df)

                ## clustering
                for (clustering_resolution in clustering_resolution_vector) {
                    print(sprintf("cluster resolution: %.1f", clustering_resolution))
                    resolution_output_dir <- sprintf("%s/res%s", dimensions_output_dir, sub("[.]", "_", clustering_resolution))
                    if (! dir.exists(resolution_output_dir)) {
                        dir.create(resolution_output_dir, recursive=TRUE, mode="0775")
                    }
                    resolution_output_name <- sprintf("%s_res%s", dimensions_output_name, sub("[.]", "_", clustering_resolution))

                    seurat_object_clustering <- clustering(seurat_object_integration, clustering_resolution, reduction2use, dimensions_nb_param)
                    write.csv(seurat_object_clustering@meta.data, file=sprintf("%s/%s_metadata.csv", resolution_output_dir, resolution_output_name), quote=FALSE, row.names=TRUE)
                    ### plot clusters onto UMAP or onto the tissue section
                    clustering_plots(seurat_object_clustering, metadata_categories, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integrated method: %s,  all spots", integration_method), resolution_output_dir, sprintf("%s_clustering", resolution_output_name))

                    ### get mean silhouette scores
                    cluster_mean_silhouette <- mean(seurat_object_clustering$silhouette_seurat_clusters)
                    cluster_nb <- length(levels(seurat_object_clustering@meta.data$seurat_clusters))
                    params_mean_silhouette_df <- data.frame(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_param, resolution=clustering_resolution, nb_clusters=cluster_nb, mean_silhouette_seurat_clusters=cluster_mean_silhouette)
                    integration_cluster_mean_silhouette_df <- rbind(integration_cluster_mean_silhouette_df, params_mean_silhouette_df)
                    
                    ### cluster markers, differentially expressed genes across conditions and conserved cell type markers
                    marker_analysis(seurat_object_clustering, default_assay, mitochondrial_geneset, hemoglobin_geneset, "negbinom", "orig.ident", cluster_markers, cluster_marker_condition, scale_data_split_by, scale_data_variables_to_regress, load_image, nb_cluster_markers_to_plot, nb_SpatialFeaturePlot_per_page, DE_analysis, timepoint_analysis, nb_DEGs_to_plot, nb_DEGs_per_page, conserved_markers, "seurat_clusters", "condition", resolution_output_dir, resolution_output_name)

                    rm(seurat_object_clustering)
                    gc()
                }
                rm(seurat_object_integration)
                gc()
            }
            rm(seurat_object_merge_norm)
            gc()
        }
    }
    rm(seurat_objects)
    gc()
    detach("package:rliger", unload=TRUE)
    detach("package:SeuratWrappers", unload=TRUE)
}

integration_normalization_levels <- c()
for (one_integration_method in integration_methods) {
    for (one_normalization_method in normalization_method_vector) {
        integration_normalization_levels <- c(integration_normalization_levels, sprintf("%s %s", one_integration_method, one_normalization_method))
    }
}

integration_cluster_mean_silhouette_df$integration <- factor(integration_cluster_mean_silhouette_df$integration, levels=integration_methods)
integration_cluster_mean_silhouette_df$integration_normalization <- sprintf("%s %s", integration_cluster_mean_silhouette_df$integration, integration_cluster_mean_silhouette_df$normalization)
integration_cluster_mean_silhouette_df$integration_normalization <- factor(integration_cluster_mean_silhouette_df$integration_normalization, levels=integration_normalization_levels)
integration_cluster_mean_silhouette_df$features <- factor(integration_cluster_mean_silhouette_df$features, levels=features_param_vector)
integration_cluster_mean_silhouette_df$dimensions <- as.factor(integration_cluster_mean_silhouette_df$dimensions)
integration_cluster_mean_silhouette_df$resolution <- as.factor(integration_cluster_mean_silhouette_df$resolution)
plot_output_name <- sprintf("%s_norm_featnb_int_dims_res_cluster_mean_silhouette", output_dirname)
write.csv(integration_cluster_mean_silhouette_df, file=sprintf("%s/%s.csv", qc_filtering_output_dir, plot_output_name), quote=FALSE, row.names=FALSE)
## modify levels after writing csv file to avoid end of lines in csv file
levels(integration_cluster_mean_silhouette_df$integration_normalization) <- unlist(lapply(levels(integration_cluster_mean_silhouette_df$integration_normalization), function(x) { return(gsub(" SCTransform", "\nSCTransform", x)) }))
levels(integration_cluster_mean_silhouette_df$integration_normalization) <- unlist(lapply(levels(integration_cluster_mean_silhouette_df$integration_normalization), function(x) { return(gsub(" LogNormalize", "\nLogNormalize", x)) }))
pdf(sprintf("%s/%s.pdf", qc_filtering_output_dir, plot_output_name))
p <- ggplot(integration_cluster_mean_silhouette_df, aes(x=resolution, y=mean_silhouette_seurat_clusters, fill=integration)) +
    geom_bar(stat="identity", position="dodge") +
    geom_text(aes(label=nb_clusters), position=position_dodge(width=1), size=1.5, vjust=2) +
    facet_grid(dimensions~features) +
    labs(x="Resolution", y="Silhouette score", fill="Integration method") +
    theme_bw() +
    theme(legend.position="bottom") +
    theme(axis.text.x=element_text(angle=60, hjust=1)) +
    theme(panel.border=element_rect(color="grey50"))
print(p)
dev.off()

integration_batch_mean_silhouette_df$integration <- factor(integration_batch_mean_silhouette_df$integration, levels=integration_methods)
integration_batch_mean_silhouette_df$integration_normalization <- sprintf("%s %s", integration_batch_mean_silhouette_df$integration, integration_batch_mean_silhouette_df$normalization)
integration_batch_mean_silhouette_df$integration_normalization <- factor(integration_batch_mean_silhouette_df$integration_normalization, levels=integration_normalization_levels)
integration_batch_mean_silhouette_df$features <- factor(integration_batch_mean_silhouette_df$features, levels=features_param_vector)
integration_batch_mean_silhouette_df$dimensions <- as.factor(integration_batch_mean_silhouette_df$dimensions)
plot_output_name <- sprintf("%s_norm_featnb_int_dims_res_batch_mean_silhouette", output_dirname)
write.csv(integration_batch_mean_silhouette_df, file=sprintf("%s/%s.csv", qc_filtering_output_dir, plot_output_name), quote=FALSE, row.names=FALSE)
## modify levels after writing csv file to avoid end of lines in csv file
levels(integration_batch_mean_silhouette_df$integration_normalization) <- unlist(lapply(levels(integration_batch_mean_silhouette_df$integration_normalization), function(x) { return(gsub(" SCTransform", "\nSCTransform", x)) }))
levels(integration_batch_mean_silhouette_df$integration_normalization) <- unlist(lapply(levels(integration_batch_mean_silhouette_df$integration_normalization), function(x) { return(gsub(" LogNormalize", "\nLogNormalize", x)) }))
pdf(sprintf("%s/%s.pdf", qc_filtering_output_dir, plot_output_name))
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


