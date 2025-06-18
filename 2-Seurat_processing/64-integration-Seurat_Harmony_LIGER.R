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
src_dir <- sprintf("%s/20-Count_analysis/10-EpiReg/src", work_dir)
genome_name <- "RefSeq108"
load_image <- TRUE

epireg_visium_metadata_df <- data.frame(sample=c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16"),
                                       condition=factor(c(rep("SE", 4), rep("CTRL", 4), rep("SE", 4), rep("CTRL", 4))), 
                                       time=factor(c(rep(c(5, 10), 4), rep(c(20, 40, 40, 20), 2))))

# sections
sections <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")
#sections <- c("C_L1_S3", "C_L4_S15", "D_L1_S4", "D_L4_S16")
#sections <- c("A_L3_S9", "B_L3_S10")
timepoint <- NA
#timepoint <- 40
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
output_dir <- sprintf("%s/20-Count_analysis/10-EpiReg/output/00-ST_P/TSO_polyA_R1trim1_ov5_n2_min20/%s/Visium_params/20-Integration/%s", work_dir, genome_name, output_dirname)
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
#normalization_method_vector <- c("SCTransform", "SCTransform_v2", "LogNormalize")
normalization_method_vector <- c("SCTransform_v2")
sct_variables_to_regress <- "qc_mito_percent_RNA"
scale_data_variables_to_regress <- "qc_mito_percent_RNA"
variables_to_regress_in_path <- "mito"
scale_data_split_by <- "orig.ident"
scale_data_split_by_in_path <- "orig"
integration_group <- "orig.ident"
integration_group_by_in_path <- "orig"
#feature_residual_variance_threshold <- 1.1
#feature_residual_variance_threshold <- NA
#feature_residual_variance_threshold_vector <- c(1, 1.1, 1.3, 1.5, 2, 2.5)
#features_param_vector <- c(1.1, 1.3, 1.5, 2)
#features_param_vector <- c(500, 1000, 1500, 2000, 2500)
#features_param_vector <- c(250, 500, 1000, 1500, 2000)
#features_param_vector <- c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000)
#features_param_vector <- c(500, 1000, 1500, 2000, 2500, 3000)
features_param_vector <- c(500, 1000)
#dimensions_nb_param_vector <- c(10, 15, 20)
#dimensions_nb_param_vector <- c(15, 20, 30)
dimensions_nb_param_vector <- c(15)
#clustering_resolution_vector <- c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2)
#clustering_resolution_vector <- c(0.4, 0.6, 0.8, 1)
#clustering_resolution_vector <- c(0.4, 1)
#clustering_resolution_vector <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)
#clustering_resolution_vector <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)
#clustering_resolution_vector <- c(0.1, 0.3, 0.5, 0.8, 1, 1.2)
clustering_resolution_vector <- c(0.1)
#clustering_resolution_vector <- c(0.1, 0.3)
#features_param <- 500
#dimensions_nb_param <- 20
#clustering_resolution <- 0.1
cluster_markers <- FALSE
DE_analysis <- FALSE
conserved_markers <- FALSE

# number of samples per page when using SpatialFeaturePlot() function
nb_SpatialFeaturePlot_per_page <- 4
# number of samples per page when using SpatialFeaturePlot() function
nb_cluster_markers_to_plot <- 4

nb_DEGs_to_plot <- 30
nb_DEGs_per_page <- 3

metadata_categories <- c("lame", "zone", "condition", "orig.ident")

merge_before_norm_in_path <- ifelse(merge_before_norm, "Merge_norm", "Norm_merge")
qc_output_dir <- sprintf("%s/%s/Norm_%s-scale_%s-int_%s/QC-sum%d_det%s_mito%d_hb%d", output_dir, merge_before_norm_in_path, variables_to_regress_in_path, scale_data_split_by_in_path, integration_group_by_in_path, spot_sum_of_counts_threshold, spot_detected_genes_threshold, spot_mito_gene_count_pct_threshold, spot_hb_gene_count_pct_threshold)

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
st_pipeline_dir <- sprintf("%s/10-Raw2count/00-ST_P/output/TSO_polyA_R1trim1_ov5_n2_min20/%s/Visium_params/00-Samples", work_dir, genome_name)
space_ranger_dir <- sprintf("%s/10-Raw2count/10-SR/output/TSO_polyA_R1trim1_ov5_n2_min20/%s/00-Samples", work_dir, genome_name)
spe_objects <- lapply(sections, load_ST_pipeline_data, st_dir=st_pipeline_dir, sr_dir=space_ranger_dir, metadata_df=epireg_visium_metadata_df)
names(spe_objects) <- sections

# quality control
library(scater) # Loading required package: scuttle
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
    space_ranger_output_dir <- sprintf("%s/10-Raw2count/10-SR/output/TSO_polyA_R1trim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs", work_dir, genome_name, sample)
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

# no integration, Seurat and Harmony
library(harmony)
integration_methods <- c("no integration", "Seurat", "Harmony")
for (normalization_method in normalization_method_vector) {
    print(sprintf("normalization method: %s", normalization_method))
    return_only_var_genes <- ifelse(normalization_method %in% c("SCTransform", "SCTransform_v2"), FALSE, NA)
    for (features_param in features_param_vector) {
        ## normalization and feature selection
        if (merge_before_norm) {
            returned_list <- seurat_merge_normalization(seurat_object_merge, normalization_method, sct_variables_to_regress, features_param, return_only_var_genes, mitochondrial_geneset_spe, hemoglobin_geneset_spe, qc_filtering_output_dir)
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
            returned_list <- seurat_list_normalization(seurat_objects, normalization_method, sct_variables_to_regress, features_param, return_only_var_genes, mitochondrial_geneset_spe, hemoglobin_geneset_spe, qc_filtering_output_dir)
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
            #### orig.ident, lame, zone, condition and time are not factors in the merged object meta data
            seurat_object_merge_norm@meta.data$orig.ident <- as.factor(seurat_object_merge_norm@meta.data$orig.ident)
            seurat_object_merge_norm@meta.data$lame <- as.factor(seurat_object_merge_norm@meta.data$lame)
            seurat_object_merge_norm@meta.data$zone <- as.factor(seurat_object_merge_norm@meta.data$zone)
            seurat_object_merge_norm@meta.data$condition <- as.factor(seurat_object_merge_norm@meta.data$condition)
            seurat_object_merge_norm@meta.data$time <- as.factor(seurat_object_merge_norm@meta.data$time)
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
                    if (normalization_method %in% c("SCTransform", "SCTransform_v2")) {
                        seurat_objects_norm <- PrepSCTIntegration(object.list=seurat_objects_norm, anchor.features=features)
                    }
                    int.anchors <- FindIntegrationAnchors(object.list=seurat_objects_norm, anchor.features=features, normalization.method=integration_normalization_method, dims=1:dimensions_nb_param)
                    rm(seurat_objects_norm)
                    gc()
                    seurat_object_to_scale <- IntegrateData(anchorset=int.anchors, normalization.method=integration_normalization_method, dims=1:dimensions_nb_param)
                    rm(int.anchors)
                    gc()
                    ##### orig.ident, lame, zone, condition and time are not factors in the integrated object meta data
                    seurat_object_to_scale@meta.data$orig.ident <- as.factor(seurat_object_to_scale@meta.data$orig.ident)
                    seurat_object_to_scale@meta.data$lame <- as.factor(seurat_object_to_scale@meta.data$lame)
                    seurat_object_to_scale@meta.data$zone <- as.factor(seurat_object_to_scale@meta.data$zone)
                    seurat_object_to_scale@meta.data$condition <- as.factor(seurat_object_to_scale@meta.data$condition)
                    seurat_object_to_scale@meta.data$time <- as.factor(seurat_object_to_scale@meta.data$time)
                } else {
                    seurat_object_to_scale <- seurat_object_merge_norm
                }

                ### scaling
                if (normalization_method == "LogNormalize") {
                    if (! is.na(scale_data_split_by)) {
                        if (!is.na(scale_data_variables_to_regress)) {
                            seurat_object_integration <- ScaleData(seurat_object_to_scale, split.by=scale_data_split_by, vars.to.regress=scale_data_variables_to_regress)
                            dimensions_output_dir <- sprintf("%s/%s_%s", dimensions_output_dir, gsub("[.]", "", scale_data_split_by), paste(gsub("[.]", "", scale_data_variables_to_regress), collapse="-"))
                        } else {
                            seurat_object_integration <- ScaleData(seurat_object_to_scale, split.by=scale_data_split_by)
                            dimensions_output_dir <- sprintf("%s/%s", dimensions_output_dir, gsub("[.]", "", scale_data_split_by))
                        }
                    } else {
                        seurat_object_integration <- ScaleData(seurat_object_to_scale)
                    }
                } else {
                    seurat_object_integration <- seurat_object_to_scale
                }
                rm(seurat_object_to_scale)
                gc()
                
                seurat_object_integration <- RunPCA(seurat_object_integration, verbose=FALSE)
                
                #### Harmony integration and dimensionality reduction
                if (integration_method == "Harmony") {
                    seurat_object_to_cluster <- RunHarmony(seurat_object_integration, group.by.vars=integration_group, dims=1:dimensions_nb_param)
                    reduction2use <- "harmony"
                } else {
                    seurat_object_to_cluster <- seurat_object_integration
                    reduction2use <- "pca"
                }
                
                seurat_object_to_cluster <- FindNeighbors(seurat_object_to_cluster, reduction=reduction2use, dims=1:dimensions_nb_param)

                ### compute silhouette scores for lames, zones, conditions and samples
                mean_silhouette_vector <- c()
                median_silhouette_vector <- c()
                for (metadata_category_colname in metadata_categories) {
                    seurat_object_to_cluster <- compute_silhouette(seurat_object_to_cluster, reduction2use, dimensions_nb_param, metadata_category_colname)
                    mean_silhouette <- mean(seurat_object_to_cluster@meta.data[[sprintf("silhouette_%s", metadata_category_colname)]])
                    mean_silhouette_vector <- c(mean_silhouette_vector, mean_silhouette)
                    median_silhouette <- median(seurat_object_to_cluster@meta.data[[sprintf("silhouette_%s", metadata_category_colname)]])
                    median_silhouette_vector <- c(median_silhouette_vector, median_silhouette)
                }
                mean_silhouette_df <- as.data.frame(t(mean_silhouette_vector))
                colnames(mean_silhouette_df) <- sprintf("mean_silhouette_%s", metadata_categories)
                median_silhouette_df <- as.data.frame(t(median_silhouette_vector))
                colnames(median_silhouette_df) <- sprintf("median_silhouette_%s", metadata_categories)
                mean_median_silhouette_df <- cbind(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_param, mean_silhouette_df, median_silhouette_df)
                integration_batch_mean_silhouette_df <- rbind(integration_batch_mean_silhouette_df, mean_median_silhouette_df)

                ## clustering
                for (clustering_resolution in clustering_resolution_vector) {
                    print(sprintf("cluster resolution: %.1f", clustering_resolution))
                    resolution_output_dir <- sprintf("%s/res%s", dimensions_output_dir, sub("[.]", "_", clustering_resolution))
                    if (! dir.exists(resolution_output_dir)) {
                        dir.create(resolution_output_dir, recursive=TRUE, mode="0775")
                    }
                    resolution_output_name <- sprintf("%s_res%s", dimensions_output_name, sub("[.]", "_", clustering_resolution))

                    seurat_object_clustering <- clustering(seurat_object_to_cluster, clustering_resolution, reduction2use, dimensions_nb_param)
                    write.csv(seurat_object_clustering@meta.data, file=sprintf("%s/%s_metadata.csv", resolution_output_dir, resolution_output_name), quote=FALSE, row.names=TRUE)
                    cluster_mean_silhouette <- mean(seurat_object_clustering$silhouette_seurat_clusters)
                    cluster_median_silhouette <- median(seurat_object_clustering$silhouette_seurat_clusters)
                    ### plot clusters onto UMAP or onto the tissue section
                    clustering_plots(seurat_object_clustering, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integrated method: %s", integration_method), resolution_output_dir, sprintf("%s_clustering", resolution_output_name))

                    ### hippocampus spots
                    barcodes2use <- all_sections_hippocampus_barcodes_Seurat_merge
                    seurat_object_clustering_hippocampus <- subset(seurat_object_clustering, cells=barcodes2use)
                    hippocampus_cluster_mean_silhouette <- mean(seurat_object_clustering_hippocampus$silhouette_seurat_clusters)
                    hippocampus_cluster_median_silhouette <- median(seurat_object_clustering_hippocampus$silhouette_seurat_clusters)
                    #### plot clusters onto UMAP or onto the tissue section
                    clustering_plots(seurat_object_clustering_hippocampus, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integrated method: %s,  hippocampus spots", integration_method), resolution_output_dir, sprintf("%s_clustering_hippocampus", resolution_output_name))

                    ### get mean silhouette scores: all spots and hippocampus spots
                    mean_median_silhouette_df <- data.frame(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_param, resolution=clustering_resolution, mean_silhouette_seurat_clusters=cluster_mean_silhouette, mean_silhouette_seurat_clusters_hippocampus=hippocampus_cluster_mean_silhouette, median_silhouette_seurat_clusters=cluster_median_silhouette, median_silhouette_seurat_clusters_hippocampus=hippocampus_cluster_median_silhouette)
                    integration_cluster_mean_silhouette_df <- rbind(integration_cluster_mean_silhouette_df, mean_median_silhouette_df)
                    
                    ### cluster markers, differentially expressed genes across conditions and conserved cell type markers
                    if (normalization_method %in% c("SCTransform", "SCTransform_v2")) {
                        different_SCT_models <- ifelse(length(seurat_object_clustering@assays$SCT@SCTModel.list) > 1, TRUE, FALSE)
                    } else {
                        different_SCT_models <- NA
                    }
                    marker_analysis(seurat_object_clustering, default_assay, different_SCT_models, mitochondrial_geneset_spe, hemoglobin_geneset_spe, "negbinom", "orig.ident", cluster_markers, scale_data_split_by, scale_data_variables_to_regress, load_image, nb_cluster_markers_to_plot, nb_SpatialFeaturePlot_per_page, DE_analysis, nb_DEGs_to_plot, nb_DEGs_per_page, conserved_markers, "seurat_clusters", "condition", resolution_output_dir, resolution_output_name)
                    
                    #rm(seurat_object_clustering)
                    rm(seurat_object_clustering_hippocampus)
                    gc()
                }
                rm(seurat_object_to_cluster)
                gc()
            }
            rm(seurat_object_integration)
            gc()
        }
        rm(seurat_object_merge_norm)
        gc()
    }
}
detach("package:harmony", unload=TRUE)


# LIGER
library(rliger)
library(SeuratWrappers)
integration_method <- "LIGER"
integration_methods <- c(integration_methods, integration_method)
for (normalization_method in normalization_method_vector) {
    print(sprintf("normalization method: %s", normalization_method))
    return_only_var_genes <- ifelse(normalization_method %in% c("SCTransform", "SCTransform_v2"), TRUE, NA) # SCTransform: only return variable genes, otherwise error with RunOptimizeALS() function: Error in dimnames(x) <- dn: length of 'dimnames' [2] not equal to array extent
    for (features_param in features_param_vector) {
        ## normalization and feature selection
        if (merge_before_norm) {
            returned_list <- seurat_merge_normalization(seurat_object_merge, normalization_method, sct_variables_to_regress, features_param, return_only_var_genes, mitochondrial_geneset_spe, hemoglobin_geneset_spe, qc_filtering_output_dir)
            seurat_object_merge_norm <- returned_list[["seurat_merge"]]
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
            hvg_df <- hvg_plots(seurat_object_merge_norm, default_assay, selection_method, features_param, normalization_output_dir, "HVGs")
        } else {
            returned_list <- seurat_list_normalization(seurat_objects, normalization_method, sct_variables_to_regress, features_param, return_only_var_genes, mitochondrial_geneset_spe, hemoglobin_geneset_spe, qc_filtering_output_dir)
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
        
        integration_output_dir <- sprintf("%s/%s", normalization_output_dir, gsub(" ", "_", integration_method))
        if (! dir.exists(integration_output_dir)) {
            dir.create(integration_output_dir, recursive=TRUE, mode="0775")
        }
        
        ### scaling
        if (normalization_method == "LogNormalize") {
            if (! is.na(scale_data_split_by)) {
                if (!is.na(scale_data_variables_to_regress)) {
                    seurat_object_integration <- ScaleData(seurat_object_merge_norm, split.by=scale_data_split_by, vars.to.regress=scale_data_variables_to_regress, do.center=FALSE) # LIGER does not center data when scaling, skip that step
                    integration_output_dir <- sprintf("%s/%s_%s", integration_output_dir, gsub("[.]", "", scale_data_split_by), paste(gsub("[.]", "", scale_data_variables_to_regress), collapse="-"))
                } else {
                    seurat_object_integration <- ScaleData(seurat_object_merge_norm, split.by=scale_data_split_by, do.center=FALSE) # LIGER does not center data when scaling, skip that step
                    integration_output_dir <- sprintf("%s/%s", integration_output_dir, gsub("[.]", "", scale_data_split_by))
                }
            } else {
                seurat_object_integration <- ScaleData(seurat_object_merge_norm, do.center=FALSE) # LIGER does not center data when scaling, skip that step
            }
        } else {
            seurat_object_integration <- seurat_object_merge_norm
        }
        rm(seurat_object_merge_norm)
        gc()
        
        ### integration and dimensionality reduction
        for (dimensions_nb_param in dimensions_nb_param_vector) {
            print(sprintf("dimension number: %d", dimensions_nb_param))
            dimensions_output_dir <- sprintf("%s/%d_dims", integration_output_dir, dimensions_nb_param)
            if (! dir.exists(dimensions_output_dir)) {
                dir.create(dimensions_output_dir, recursive=TRUE, mode="0775")
            }
            dimensions_output_name <- sprintf("%s_%s_%ddims", normalization_output_name, gsub(" ", "_", integration_method), dimensions_nb_param)
            
            ### integration and dimensionality reduction
            seurat_object_to_cluster <- RunOptimizeALS(seurat_object_integration, k=dimensions_nb_param, lambda=5, split.by=integration_group)
            seurat_object_to_cluster <- RunQuantileNorm(seurat_object_to_cluster, split.by=integration_group)
            reduction2use <- "iNMF"

            seurat_object_to_cluster <- FindNeighbors(seurat_object_to_cluster, reduction=reduction2use, dims=1:dimensions_nb_param)

            ### compute silhouette scores for lames, zones, conditions and samples
            mean_silhouette_vector <- c()
            median_silhouette_vector <- c()
            for (metadata_category_colname in metadata_categories) {
                seurat_object_to_cluster <- compute_silhouette(seurat_object_to_cluster, reduction2use, dimensions_nb_param, metadata_category_colname)
                mean_silhouette <- mean(seurat_object_to_cluster@meta.data[[sprintf("silhouette_%s", metadata_category_colname)]])
                mean_silhouette_vector <- c(mean_silhouette_vector, mean_silhouette)
                median_silhouette <- median(seurat_object_to_cluster@meta.data[[sprintf("silhouette_%s", metadata_category_colname)]])
                median_silhouette_vector <- c(median_silhouette_vector, median_silhouette)
            }
            mean_silhouette_df <- as.data.frame(t(mean_silhouette_vector))
            colnames(mean_silhouette_df) <- sprintf("mean_silhouette_%s", metadata_categories)
            median_silhouette_df <- as.data.frame(t(median_silhouette_vector))
            colnames(median_silhouette_df) <- sprintf("median_silhouette_%s", metadata_categories)
            mean_median_silhouette_df <- cbind(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_param, mean_silhouette_df, median_silhouette_df)
            integration_batch_mean_silhouette_df <- rbind(integration_batch_mean_silhouette_df, mean_median_silhouette_df)
            
            ## clustering
            for (clustering_resolution in clustering_resolution_vector) {
                print(sprintf("cluster resolution: %.1f", clustering_resolution))
                resolution_output_dir <- sprintf("%s/res%s", dimensions_output_dir, sub("[.]", "_", clustering_resolution))
                if (! dir.exists(resolution_output_dir)) {
                    dir.create(resolution_output_dir, recursive=TRUE, mode="0775")
                }
                resolution_output_name <- sprintf("%s_res%s", dimensions_output_name, sub("[.]", "_", clustering_resolution))

                seurat_object_clustering <- clustering(seurat_object_to_cluster, clustering_resolution, reduction2use, dimensions_nb_param)
                write.csv(seurat_object_clustering@meta.data, file=sprintf("%s/%s_metadata.csv", resolution_output_dir, resolution_output_name), quote=FALSE, row.names=TRUE)
                cluster_mean_silhouette <- mean(seurat_object_clustering$silhouette_seurat_clusters)
                cluster_median_silhouette <- median(seurat_object_clustering$silhouette_seurat_clusters)
                ### plot clusters onto UMAP or onto the tissue section
                clustering_plots(seurat_object_clustering, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integrated method: %s,  all spots", integration_method), resolution_output_dir, sprintf("%s_clustering", resolution_output_name))

                ### hippocampus spots
                seurat_object_clustering_hippocampus <- subset(seurat_object_clustering, cells=all_sections_hippocampus_barcodes_Seurat_merge)
                hippocampus_cluster_mean_silhouette <- mean(seurat_object_clustering_hippocampus$silhouette_seurat_clusters)
                hippocampus_cluster_median_silhouette <- median(seurat_object_clustering_hippocampus$silhouette_seurat_clusters)
                #### plot clusters onto UMAP or onto the tissue section
                clustering_plots(seurat_object_clustering_hippocampus, load_image, nb_SpatialFeaturePlot_per_page, sprintf("Integrated method: %s,  hippocampus spots", integration_method), resolution_output_dir, sprintf("%s_clustering_hippocampus", resolution_output_name))

                ### get mean silhouette scores: all spots and hippocampus spots
                mean_silhouette_df <- data.frame(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_param, resolution=clustering_resolution, mean_silhouette_seurat_clusters=cluster_mean_silhouette, mean_silhouette_seurat_clusters_hippocampus=hippocampus_cluster_mean_silhouette)
                mean_median_silhouette_df <- data.frame(normalization=normalization_method, features=features_param, integration=integration_method, dimensions=dimensions_nb_param, resolution=clustering_resolution, mean_silhouette_seurat_clusters=cluster_mean_silhouette, mean_silhouette_seurat_clusters_hippocampus=hippocampus_cluster_mean_silhouette, median_silhouette_seurat_clusters=cluster_median_silhouette, median_silhouette_seurat_clusters_hippocampus=hippocampus_cluster_median_silhouette)
                integration_cluster_mean_silhouette_df <- rbind(integration_cluster_mean_silhouette_df, mean_median_silhouette_df)
                
                ### cluster markers, differentially expressed genes across conditions and conserved cell type markers
                if (normalization_method %in% c("SCTransform", "SCTransform_v2")) {
                    different_SCT_models <- ifelse(length(seurat_object_clustering@assays$SCT@SCTModel.list) > 1, TRUE, FALSE)
                } else {
                    different_SCT_models <- NA
                }
                marker_analysis(seurat_object_clustering, default_assay, different_SCT_models, mitochondrial_geneset_spe, hemoglobin_geneset_spe, "wilcox", NA, cluster_markers, scale_data_split_by, scale_data_variables_to_regress, load_image, nb_cluster_markers_to_plot, nb_SpatialFeaturePlot_per_page, DE_analysis, nb_DEGs_to_plot, nb_DEGs_per_page, conserved_markers, "seurat_clusters", "condition", resolution_output_dir, resolution_output_name)
                
                rm(seurat_object_clustering)
                rm(seurat_object_clustering_hippocampus)
                gc()
            }
            rm(seurat_object_to_cluster)
            gc()
        }
        rm(seurat_object_integration)
        gc()
    }
}
detach("package:rliger", unload=TRUE)
detach("package:SeuratWrappers", unload=TRUE)
if (! merge_before_norm) {
    rm(seurat_objects)
    gc()
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
for (one_category in c("seurat_clusters", "seurat_clusters_hippocampus")) {
    p <- ggplot(integration_cluster_mean_silhouette_df, aes(x=resolution, y=.data[[sprintf("mean_silhouette_%s", one_category)]], fill=integration_normalization)) +
    #p <- ggplot(integration_cluster_mean_silhouette_df, aes(x=resolution, y=!!sym(sprintf("mean_silhouette_%s", one_category)), fill=integration)) +
        geom_bar(stat="identity", position="dodge") +
        facet_grid(dimensions~features) +
        labs(title=sprintf("Mean silhouette: %s", one_category), x="Resolution", y="Mean silhouette", fill="Integration method") +
        theme_bw() +
        theme(legend.position="bottom") +
        theme(axis.text.x=element_text(angle=60, hjust=1)) +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
}
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




