.libPaths(c(.libPaths(), "/home/christophe.lepriol/NeuroDev_ADD/R/r_4.2.0"))
#.libPaths(c(.libPaths(), "/home/christophe.lepriol/NeuroDev_ADD/R/r_4.1.0"))
library(ggplot2)
library(patchwork)
library(cowplot)
library(Seurat)
library(SpatialExperiment)
library(ggspavis)


##############
# Parameters #
##############
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"
src_dir <- sprintf("%s/20-Count_analysis/10-EpiReg/src", work_dir)
genome_name <- "RefSeq108"

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
#cluster_markers <- FALSE
#DE_analysis <- FALSE
#DE_timepoint_analysis <- FALSE
#conserved_markers <- FALSE
markers_test <- "negbinom"
markers_latent_variables <- "orig.ident"

## subclusterings
subclustering_integration <- "Seurat"
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
cluster_0_parameters[["subclusters"]] <- 4
subsclusterings_parameter_settings[[subclustering_name]] <- cluster_0_parameters
#### cluster 3
subclustering_name <- "3"
cluster_3_parameters <- list()
cluster_3_parameters[["Seurat"]]  <- subclustering_name
cluster_3_parameters[["normalization"]] <- "SCTransform_v2"
cluster_3_parameters[["features"]]  <- 2000
cluster_3_parameters[["integration"]] <- "Seurat"
cluster_3_parameters[["dimensions"]]  <- 10
cluster_3_parameters[["resolution"]] <- 0.05
cluster_3_parameters[["subclusters"]] <- 4
subsclusterings_parameter_settings[[subclustering_name]] <- cluster_3_parameters


#############
# Functions #
#############
source(sprintf("%s/6X-integration_functions.R", src_dir))


############
# Analysis #
############

# process parameters
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
## variables to regress
vars2regress_in_dirname <- ifelse(is.na(sct_variables_to_regress), "no_vars_to_regress", paste(gsub("[.]", "", sct_variables_to_regress), collapse="-"))
## markers latent variables
markers_latent_variables_in_path <- ifelse(length(markers_latent_variables > 1), paste(gsub("[.]", "", markers_latent_variables), collapse="-"), ifelse(is.na(markers_latent_variables), "no_var", gsub("[.]", "", markers_latent_variables)))

## output directory
output_dirname <- ifelse(is.na(timepoint), "all", sprintf("D%d", timepoint))

# plot GO term p-value per cluster on tissue section
## load whole clustering Seurat object and list of subclustering Seurat objects
# integration_output_dir <- sprintf("%s/20-Count_analysis/10-EpiReg/output/%s/TSO_polyA_R1trim1_ov5_n2_min20/%s/20-Integration/%s/sub_markers_20230610/Norm_%s-scale_%s-int_%s/QC-sum%d_det%s_mito%d_hb%d/10-QC_filtering/%s", work_dir, input_dirname, genome_name, output_dirname, variables_to_regress_in_path, scale_data_split_by_in_path, integration_group_by_in_path, spot_sum_of_counts_threshold, spot_detected_genes_threshold, spot_mito_gene_count_pct_threshold, spot_hb_gene_count_pct_threshold, integration_method)
integration_output_dir <- sprintf("%s/20-Count_analysis/10-EpiReg/output/%s/TSO_polyA_R1trim1_ov5_n2_min20/%s/20-Integration/%s/sub_markers/Norm_%s-scale_%s-int_%s/QC-sum%d_det%s_mito%d_hb%d/10-QC_filtering/%s", work_dir, input_dirname, genome_name, output_dirname, variables_to_regress_in_path, scale_data_split_by_in_path, integration_group_by_in_path, spot_sum_of_counts_threshold, spot_detected_genes_threshold, spot_mito_gene_count_pct_threshold, spot_hb_gene_count_pct_threshold, integration_method)
load(sprintf("%s/whole_clustering_subclustering_Seurat_objects.RData", integration_output_dir))
## remove unnecessary seurat_object_subclustering_list object
rm(seurat_object_subclustering_list)
gc()

# Seurat to SpatialExperiment object conversion
## code adapted from: https://github.com/drighelli/SpatialExperiment/issues/115
sce <- Seurat::as.SingleCellExperiment(seurat_object_whole_clustering)
sample_ids <- names(seurat_object_whole_clustering@images)
### create SpatialExperiment object
#### detach Seurat and SeuratObject packages to avoid conflict with SpatialExperiment package
#detach("package:Seurat", unload=TRUE)
#detach("package:SeuratObject", unload=TRUE)
#library(SpatialExperiment)
#library(ggspavis)
#### spatial coordinates
#spatialCoords <- do.call("rbind", lapply(sample_ids, function(x, seurat_obj=seurat_object_whole_clustering) {
#    as.matrix(seurat_obj@images[[x]]@coordinates[, c("imagecol", "imagerow")])
#}))
#### image data
imgData <- do.call("rbind", lapply(sample_ids, function(x, seurat_obj=seurat_object_whole_clustering) {
    img <- SpatialImage(x=as.raster(seurat_obj@images[[x]]@image))
    imgDFrame <- DataFrame(sample_id=x, image_id=x, data=I(list(img)), scaleFactor=seurat_obj@images[[x]]@scale.factors$lowres)
    return(imgDFrame)
}))
spe <- SpatialExperiment(assays=assays(sce), rowData=rowData(sce), colData=colData(sce), metadata=metadata(sce), reducedDims=reducedDims(sce), altExps=altExps(sce), sample_id=as.character(colData(sce)$orig.ident), spatialCoordsNames=c("pxl_col_in_fullres", "pxl_row_in_fullres"), imgData=imgData)
#### subset: only a single representative sample per timepoint (one of SE samples)
sample_ids <- c("A_L1_S1", "A_L2_S5", "C_L1_S3", "C_L2_S7")
sample_id_times <- c("5", "10", "20", "40")
names(sample_id_times) <- sample_ids
spe_sub <- spe[, colData(spe)$sample_id %in% sample_ids]
all_barcodes <- rownames(colData(spe_sub))

# output directory
sub_markers_output_dir <- sprintf("%s/20-Count_analysis/10-EpiReg/output/%s/TSO_polyA_R1trim1_ov5_n2_min20/%s/20-Integration/%s/sub_markers", work_dir, input_dirname, genome_name, output_dirname)
qc_filtering_output_subdir <- sprintf("Norm_%s-scale_%s-int_%s/QC-sum%d_det%s_mito%d_hb%d/10-QC_filtering", variables_to_regress_in_path, scale_data_split_by_in_path, integration_group_by_in_path, spot_sum_of_counts_threshold, spot_detected_genes_threshold, spot_mito_gene_count_pct_threshold, spot_hb_gene_count_pct_threshold)

# get all enriched GO terms for all time.(sub)cluster combinations: clusters 1, 2 and 4 and subclusters of clusters 0 and 3
all_times_GO_ORA_results2compare_df <- data.frame()
all_times_GO_ORA_simp_results2compare_df <- data.frame()
## get whole clustering GO analysis results for clusters 1, 2, and 4
### parameter setting: get whole clustering parameter evaluation output
if (visium_params_in_path) {
    qc_filtering_output_dir <- sprintf("%s/20-Count_analysis/10-EpiReg/output/%s/TSO_polyA_R1trim1_ov5_n2_min20/%s/Visium_params/20-Integration/%s/whole/Norm_%s-scale_%s-int_%s/QC-sum%d_det%s_mito%d_hb%d/10-QC_filtering", work_dir, input_dirname, genome_name, output_dirname, variables_to_regress_in_path, scale_data_split_by_in_path, integration_group_by_in_path, spot_sum_of_counts_threshold, spot_detected_genes_threshold, spot_mito_gene_count_pct_threshold, spot_hb_gene_count_pct_threshold)
} else {
    qc_filtering_output_dir <- sprintf("%s/20-Count_analysis/10-EpiReg/output/%s/TSO_polyA_R1trim1_ov5_n2_min20/%s/20-Integration/%s/whole/Norm_%s-scale_%s-int_%s/QC-sum%d_det%s_mito%d_hb%d/10-QC_filtering", work_dir, input_dirname, genome_name, output_dirname, variables_to_regress_in_path, scale_data_split_by_in_path, integration_group_by_in_path, spot_sum_of_counts_threshold, spot_detected_genes_threshold, spot_mito_gene_count_pct_threshold, spot_hb_gene_count_pct_threshold)
}
whole_clustering_eval_df <- read.csv(sprintf("%s/%s_norm_featnb_int_dims_res_cluster_mean_silhouette.csv", qc_filtering_output_dir, output_dirname))
### set best whole clustering parameter setting
whole_params <- clustering_top_params(whole_clustering_eval_df, "mean_silhouette_seurat_clusters", integration_method)
whole_normalization <- whole_params[["normalization"]]
whole_features <- whole_params[["features"]]
whole_integration <- whole_params[["integration"]]
whole_dimensions <- whole_params[["dimensions"]]
whole_resolution <- whole_params[["resolution"]]
whole_nb_clusters <- 5
print(sprintf("normalization: %s, feature number: %d, integration: %s, dimension number: %d, resolution: %.2f", whole_normalization, whole_features, whole_integration, whole_dimensions, whole_resolution))
### build path to per-timepoint DE analysis
whole_clustering_output_dir <- sprintf("%s/%s/%s/whole", sub_markers_output_dir, qc_filtering_output_subdir, whole_integration)
method_dir_index <- ifelse(whole_normalization == "SCTransform", "00", "01")
normalization_output_dir <- sprintf("%s/%s-%s/%s/featnb%s", whole_clustering_output_dir, method_dir_index, whole_normalization, vars2regress_in_dirname, whole_features)
resolution_output_dir <- sprintf("%s/%s/%d_dims/res%s", normalization_output_dir, gsub(" ", "_", whole_integration), whole_dimensions, sub("[.]", "_", whole_resolution))
comp_output_dir <- sprintf("%s/10-DEGs/%s/%s/comp", resolution_output_dir, markers_test, markers_latent_variables_in_path)
### build per-timepoint DE analysis output basename
GO_output_name <- sprintf("%s_featnb%d_%s_%ddims_res%s_DEGs_all_times_all_clusters", whole_normalization, whole_features, gsub(" ", "_", whole_integration), whole_dimensions, sub("[.]", "_", whole_resolution))
### read whole clustering GO analysis results
all_times_all_clusters_GO_ORA_df <- read.table(sprintf("%s/%s_GO_ORA.tsv", comp_output_dir, GO_output_name), header=TRUE, sep="\t", quote="")
all_times_all_clusters_GO_ORA_df$time <- as.factor(all_times_all_clusters_GO_ORA_df$time)
all_times_all_clusters_GO_ORA_df$cluster <- as.factor(all_times_all_clusters_GO_ORA_df$cluster)
all_times_all_clusters_GO_ORA_df$category <- as.factor(all_times_all_clusters_GO_ORA_df$category)
#### only GO analysis results for clusters 1, 2 and 4
all_times_clusters_124_GO_ORA_df <- all_times_all_clusters_GO_ORA_df[which(all_times_all_clusters_GO_ORA_df$cluster %in% c("1", "2", "4")),]
all_times_clusters_124_GO_ORA_df <- droplevels(all_times_clusters_124_GO_ORA_df)
all_times_GO_ORA_results2compare_df <- rbind(all_times_GO_ORA_results2compare_df, all_times_clusters_124_GO_ORA_df)
### simplified GO ORA, cutoff: 0.8
all_times_all_clusters_GO_ORA_simp_df <- read.table(file=sprintf("%s/%s_GO_ORA_simp.tsv", comp_output_dir, GO_output_name), header=TRUE, sep="\t", quote="")
all_times_all_clusters_GO_ORA_simp_df$time <- as.factor(all_times_all_clusters_GO_ORA_simp_df$time)
all_times_all_clusters_GO_ORA_simp_df$cluster <- as.factor(all_times_all_clusters_GO_ORA_simp_df$cluster)
all_times_all_clusters_GO_ORA_simp_df$category <- as.factor(all_times_all_clusters_GO_ORA_simp_df$category)
all_times_all_clusters_GO_ORA_simp_df <- all_times_all_clusters_GO_ORA_simp_df[which(all_times_all_clusters_GO_ORA_simp_df$similarity_cutoff == 0.8),]
all_times_all_clusters_GO_ORA_simp_df <- droplevels(all_times_all_clusters_GO_ORA_simp_df)
#### only GO analysis results for clusters 1, 2 and 4
all_times_clusters_124_GO_ORA_simp_df <- all_times_all_clusters_GO_ORA_simp_df[which(all_times_all_clusters_GO_ORA_simp_df$cluster %in% c("1", "2", "4")),]
all_times_clusters_124_GO_ORA_simp_df <- droplevels(all_times_clusters_124_GO_ORA_simp_df)
all_times_GO_ORA_simp_results2compare_df <- rbind(all_times_GO_ORA_simp_results2compare_df, all_times_clusters_124_GO_ORA_simp_df)

## get subclustering GO analysis results for clusters 0 and 3
for (one_cluster in c("0", "3")) {
    ### set subclustering parameters
    sub_normalization <- subsclusterings_parameter_settings[[one_cluster]][["normalization"]]
    sub_features <- subsclusterings_parameter_settings[[one_cluster]][["features"]]
    sub_integration <- subsclusterings_parameter_settings[[one_cluster]][["integration"]]
    sub_dimensions <- subsclusterings_parameter_settings[[one_cluster]][["dimensions"]]
    sub_resolution <- subsclusterings_parameter_settings[[one_cluster]][["resolution"]]
    sub_nb_clusters <- subsclusterings_parameter_settings[[one_cluster]][["subclusters"]]
    ### build path to per-timepoint DE analysis
    subclustering_output_dir <- sprintf("%s/%s/%s/sub/c%s", sub_markers_output_dir, qc_filtering_output_subdir, subclustering_integration, one_cluster)
    method_dir_index <- ifelse(sub_normalization == "SCTransform", "00", "01")
    normalization_output_dir <- sprintf("%s/%s-%s/%s/featnb%s", subclustering_output_dir, method_dir_index, sub_normalization, vars2regress_in_dirname, sub_features)
    resolution_output_dir <- sprintf("%s/%s/%d_dims/res%s", normalization_output_dir, gsub(" ", "_", sub_integration), sub_dimensions, sub("[.]", "_", sub_resolution))
    comp_output_dir <- sprintf("%s/10-DEGs/%s/%s/comp", resolution_output_dir, markers_test, markers_latent_variables_in_path)
    ### build per-timepoint DE analysis output basename
    GO_output_name <- sprintf("%s_featnb%d_%s_%ddims_res%s_DEGs_all_times_all_clusters", sub_normalization, sub_features, gsub(" ", "_", sub_integration), sub_dimensions, sub("[.]", "_", sub_resolution))
    #GO_output_name <- sprintf("%s_featnb%d_%s_%ddims_res%s_DEGs_all_times_all_subclusters", sub_normalization, sub_features, gsub(" ", "_", sub_integration), sub_dimensions, sub("[.]", "_", sub_resolution))
    
    ### read subclustering GO analysis results
    #all_times_all_clusters_GO_ORA_df <- read.table(sprintf("%s/%s_GO_ORA.tsv", comp_output_dir, GO_output_name), header=TRUE, sep="\t", quote="")
    all_times_all_clusters_GO_ORA_df <- read.table(sprintf("%s/%s_GO_ORA.tsv", comp_output_dir, GO_output_name), header=TRUE, sep="\t", quote="")
    all_times_all_clusters_GO_ORA_df$time <- as.factor(all_times_all_clusters_GO_ORA_df$time)
    all_times_all_clusters_GO_ORA_df$cluster <- as.factor(all_times_all_clusters_GO_ORA_df$cluster)
    all_times_all_clusters_GO_ORA_df$category <- as.factor(all_times_all_clusters_GO_ORA_df$category)
    #### replace cluster column by subcluster
    all_times_all_clusters_GO_ORA_df$cluster <- sprintf("%s_%s", one_cluster, all_times_all_clusters_GO_ORA_df$cluster)
    all_times_all_clusters_GO_ORA_df$cluster <- as.factor(all_times_all_clusters_GO_ORA_df$cluster)
    all_times_GO_ORA_results2compare_df <- rbind(all_times_GO_ORA_results2compare_df, all_times_all_clusters_GO_ORA_df)
    
    ### simplified GO ORA, cutoff: 0.8
    all_times_all_clusters_GO_ORA_simp_df <- read.table(file=sprintf("%s/%s_GO_ORA_simp.tsv", comp_output_dir, GO_output_name), header=TRUE, sep="\t", quote="")
    #all_times_all_clusters_GO_ORA_simp_df <- read.table(file=sprintf("%s/%s_GO_ORA_simp.tsv", comp_output_dir, GO_output_name), header=TRUE, sep="\t", quote="", row.names=1)
    all_times_all_clusters_GO_ORA_simp_df$time <- as.factor(all_times_all_clusters_GO_ORA_simp_df$time)
    all_times_all_clusters_GO_ORA_simp_df$cluster <- as.factor(all_times_all_clusters_GO_ORA_simp_df$cluster)
    all_times_all_clusters_GO_ORA_simp_df$category <- as.factor(all_times_all_clusters_GO_ORA_simp_df$category)
    all_times_all_clusters_GO_ORA_simp_df <- all_times_all_clusters_GO_ORA_simp_df[which(all_times_all_clusters_GO_ORA_simp_df$similarity_cutoff == 0.8),]
    all_times_all_clusters_GO_ORA_simp_df <- droplevels(all_times_all_clusters_GO_ORA_simp_df)
    #### replace cluster column by subcluster
    all_times_all_clusters_GO_ORA_simp_df$cluster <- sprintf("%s_%s", one_cluster, all_times_all_clusters_GO_ORA_simp_df$cluster)
    all_times_all_clusters_GO_ORA_simp_df$cluster <- as.factor(all_times_all_clusters_GO_ORA_simp_df$cluster)
    all_times_GO_ORA_simp_results2compare_df <- rbind(all_times_GO_ORA_simp_results2compare_df, all_times_all_clusters_GO_ORA_simp_df)
}

# replace '_' by '.' in subcluster names to fit with 'seurat_custom_clusters' metadata column
all_times_GO_ORA_results2compare_df$cluster <- as.factor(sub("_", ".", all_times_GO_ORA_results2compare_df$cluster))
#all_times_GO_ORA_simp_results2compare_df$cluster <- as.factor(sub("_", ".", #all_times_GO_ORA_simp_results2compare_df$cluster))

custom_cluster_comp_output_dir <- sprintf("%s/%s/%s/custom/10-DEGs/%s/%s/dynamics", sub_markers_output_dir, qc_filtering_output_subdir, subclustering_integration, markers_test, markers_latent_variables_in_path)
if (! dir.exists(custom_cluster_comp_output_dir)) {
    dir.create(custom_cluster_comp_output_dir, recursive=TRUE, mode="0775")
}

# get mean p-value and number of time and cluster combinations with significant p-value per GO term
mean_pval_df <- aggregate(all_times_GO_ORA_results2compare_df$p.adjust, by=list(all_times_GO_ORA_results2compare_df$ID, all_times_GO_ORA_results2compare_df$Description), mean)
colnames(mean_pval_df) <- c("ID", "Description", "mean_pval")
time_cluster_category_count_df <- aggregate(all_times_GO_ORA_results2compare_df$p.adjust, by=list(all_times_GO_ORA_results2compare_df$ID, all_times_GO_ORA_results2compare_df$Description), length)
colnames(time_cluster_category_count_df) <- c("ID", "Description", "count")
mean_pval_time_cluster_category_count_df <- merge(mean_pval_df, time_cluster_category_count_df, by=c("ID", "Description"))


# dynamics sets of selected GO terms based on custom cluster compareCluster() per time results
all_times_dynamics_go_terms <- list()
all_times_dynamics_go_terms[["1"]] <- c("immune effector process", "regulation of immune response", "leukocyte migration", "phagocytosis", "cytokine production", "antigen processing and presentation", "leukocyte proliferation", "lymphocyte proliferation", "mononuclear cell proliferation", "regulation of lymphocyte proliferation", "antigen processing and presentation of peptide antigen via MHC class II", "humoral immune response", "regulation of cytokine production", "integrin-mediated signaling pathway", "innate immune response", "cytokine-mediated signaling pathway", "leukocyte homeostasis", "neutrophil homeostasis", "myeloid leukocyte migration", "interleukin-10 production", "regulation of interleukin-10 production", "phagocytosis, engulfment", "interleukin-1 production", "regulation of interleukin-1 production", "chemokine production", "regulation of chemokine production", "interferon−gamma production", "regulation of interferon−gamma production", "interleukin-8 production", "regulation of interleukin-8 production", "leukocyte apoptotic process", "antigen processing and presentation of peptide antigen", "positive regulation of immune system process", "B cell mediated immunity", "immunoglobulin mediated immune response", "antigen processing and presentation of exogenous peptide antigen", "regulation of leukocyte apoptotic proces", "antigen processing and presentation of exogenous antigen", "inflammatory response", "production of molecular mediator involved in inflammatory response", "inflammatory response to wounding", "blood vessel morphogenesis", "blood vessel development", "inflammatory response to antigenic stimulus", "amyloid-beta clearance", "response to wounding", "response to mechanical stimulus", "regeneration", "cellular response to external stimulus", "positive regulation of neuron death", "apoptotic signaling pathway", "dendritic cell differentiation", "regulation of neuron death", "gliogenesis", "glial cell activation", "glial cell differentiation", "astrocyte development", "endocytosis", "positive regulation of peptidyl−tyrosine phosphorylation", "regulation of peptidase activity")
all_times_dynamics_go_terms[["2"]] <- c("axonogenesis", "axon development", "neural nucleus development", "modulation of chemical synaptic transmission", "regulation of membrane potential", "synapse organization", "vesicle-mediated transport in synapse", "regulation of trans-synaptic signaling", "synaptic vesicle cycle", "glutamate receptor signaling pathway", "ionotropic glutamate receptor signaling pathway", "excitatory postsynaptic potential", "calcium-mediated signaling", "second-messenger-mediated signaling", "synapse maturation", "regulation of neurotransmitter receptor activity", "regulation of transmembrane transport", "synapse maturation", "potassium ion transmembrane transport")
all_times_dynamics_go_terms[["3"]] <- c("axon ensheathment", "ensheathment of neurons", "myelination", "oligodendrocyte differentiation")
all_times_dynamics_go_terms[["4"]] <- c("dendrite development", "axo-dendritic transport", "neuron migration")
all_times_dynamics_go_terms[["5"]] <- c("synapse pruning")
all_times_dynamics_go_terms[["6"]] <- c("oxidative phosphorylation", "generation of precursor metabolites and energy", "NADH dehydrogenase complex assembly", "mitochondrial respiratory chain complex I assembly", "proton transmembrane transport", "mitochondrion organization", "mitochondrial respiratory chain complex assembly", "electron transport chain")
all_times_dynamics_go_terms[["7"]] <- c("NADH regeneration", "canonical glycolysis", "glucose catabolic process to pyruvate")
all_times_dynamics_go_terms[["8"]] <- c("ribosome assembly", "ribosomal small subunit assembly", "ribosome biogenesis", "ribosomal small subunit biogenesis", "translational elongation", "positive regulation of translation")

all_dynamics_all_times_GO_ORA_results2compare_df <- data.frame()
for (one_dynamics in names(all_times_dynamics_go_terms)) {
    one_dynamics_go_terms <- all_times_dynamics_go_terms[[one_dynamics]]
    # get ID corresponding to description
    one_dynamics_go_ids <- c()
    for (one_go_term in one_dynamics_go_terms) {
        one_dynamics_go_ids <- c(one_dynamics_go_ids, unique(all_times_GO_ORA_results2compare_df[which(all_times_GO_ORA_results2compare_df$Description == one_go_term), "ID"]))
    }
    # spatial plots
    one_dynamics_custom_cluster_comp_spatial_output_name <- sprintf("all_times_dynamics_%s_compareCluster_GO_ORA_form_simp_Rel_0_8_compareClusterResult_spatial_GO_pval", one_dynamics)
    one_dynamics_all_times_GO_ORA_results2compare_df <- go_pval_time_spatial_plots_2(spe_sub, all_times_GO_ORA_results2compare_df, one_dynamics_go_ids, custom_cluster_comp_output_dir, one_dynamics_custom_cluster_comp_spatial_output_name)
    all_dynamics_all_times_GO_ORA_results2compare_df <- rbind(all_dynamics_all_times_GO_ORA_results2compare_df, cbind(dynamic=rep(one_dynamics, dim(one_dynamics_all_times_GO_ORA_results2compare_df)[1]), one_dynamics_all_times_GO_ORA_results2compare_df))
}
## write article supplementary table
write.table(all_dynamics_all_times_GO_ORA_results2compare_df, file=sprintf("%s/all_times_all_dynamics_compareCluster_GO_ORA_form_simp_Rel_0_8_compareClusterResult_spatial_GO_pval.tsv", custom_cluster_comp_output_dir), sep="\t", quote=FALSE, row.names=FALSE)

# dynamics figure
## 2 pages: 3 and 2 dynamics
pdf(sprintf("%s/all_times_dynamics_compareCluster_GO_ORA_form_simp_Rel_0_8_compareClusterResult_spatial_GO_pval_figure_3_2.pdf", custom_cluster_comp_output_dir))
### 3 dynamics: immune respons (regulation of immune response), synapse (modulation of chemical synaptic transmission), synapse pruning
sample_ids <- c()
sample_id_times <- c()
spe_sub_go_pval_1 <- spe2speGOpval(spe_sub, "GO:0050776", all_times_GO_ORA_results2compare_df, "1")
sample_ids <- c(sample_ids, spe_sub_go_pval_1@int_metadata$imgData$sample_id)
sample_id_times <- c(sample_id_times, c("5", "10", "20", "40"))
spe_sub_go_pval_2 <- spe2speGOpval(spe_sub, "GO:0050804", all_times_GO_ORA_results2compare_df, "2")
sample_ids <- c(sample_ids, spe_sub_go_pval_2@int_metadata$imgData$sample_id)
sample_id_times <- c(sample_id_times, rep("", 4)) # only plot times for the first GO ID
spe_sub_go_pval_3 <- spe2speGOpval(spe_sub, "GO:0098883", all_times_GO_ORA_results2compare_df, "3")
sample_ids <- c(sample_ids, spe_sub_go_pval_3@int_metadata$imgData$sample_id)
sample_id_times <- c(sample_id_times, rep("", 4)) # only plot times for the first GO ID

spe_sub_go_pval <- cbind(spe_sub_go_pval_1, spe_sub_go_pval_2, spe_sub_go_pval_3)
names(sample_id_times) <- sample_ids
min_pval <- min(colData(spe_sub_go_pval)$GO_pvalue2plot, na.rm=TRUE)
max_pval <- max(colData(spe_sub_go_pval)$GO_pvalue2plot, na.rm=TRUE)

p <- plotVisium(spe_sub_go_pval, fill="GO_pvalue2plot", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
p <- p + scale_fill_gradientn(colours=c("blue", "grey50", "red"), na.value="white", values=scales::rescale(c(min_pval, -1, 0, 1, max_pval))) +
    theme(legend.position="bottom", legend.key.size=unit(1, "lines"), legend.title=element_text(size=8), legend.text=element_text(size=7)) 
## modify facet number of columns
p$facet$params$ncol <- 4
## modify facet titles
p$facet$params$labeller <- labeller(sample_id=sample_id_times)
## modify point size
p$layers[[length(p$layers)]]$aes_params$size <- 0.5
## modify legend title
p$guides$fill$title <- "-1^{0,1}*log10(p-value)"
print(p)

### 2 dynamics: cellular metabolism (mitochondrial respiratory chain complex I assembly), translation (ribosome assembly)
sample_ids <- c()
sample_id_times <- c()
spe_sub_go_pval_1 <- spe2speGOpval(spe_sub, "GO:0032981", all_times_GO_ORA_results2compare_df, "1")
sample_ids <- c(sample_ids, spe_sub_go_pval_1@int_metadata$imgData$sample_id)
sample_id_times <- c(sample_id_times, c("5", "10", "20", "40"))
spe_sub_go_pval_2 <- spe2speGOpval(spe_sub, "GO:0042255", all_times_GO_ORA_results2compare_df, "2")
sample_ids <- c(sample_ids, spe_sub_go_pval_2@int_metadata$imgData$sample_id)
sample_id_times <- c(sample_id_times, rep("", 4)) # only plot times for the first GO ID

spe_sub_go_pval <- cbind(spe_sub_go_pval_1, spe_sub_go_pval_2)
names(sample_id_times) <- sample_ids
min_pval <- min(colData(spe_sub_go_pval)$GO_pvalue2plot, na.rm=TRUE)
max_pval <- max(colData(spe_sub_go_pval)$GO_pvalue2plot, na.rm=TRUE)

p <- plotVisium(spe_sub_go_pval, fill="GO_pvalue2plot", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
p <- p + scale_fill_gradientn(colours=c("blue", "grey50", "red"), na.value="white", values=scales::rescale(c(min_pval, -1, 0, 1, max_pval))) +
    theme(legend.position="bottom", legend.key.size=unit(1, "lines"), legend.title=element_text(size=8), legend.text=element_text(size=7)) 
## modify facet number of columns
p$facet$params$ncol <- 4
## modify facet titles
p$facet$params$labeller <- labeller(sample_id=sample_id_times)
## modify point size
p$layers[[length(p$layers)]]$aes_params$size <- 0.5
## modify legend title
p$guides$fill$title <- "-1^{0,1}*log10(p-value)"
print(p)
dev.off()

## a single page: 5 dynamics
### 5 dynamics: immune respons (regulation of immune response), synapse (modulation of chemical synaptic transmission), synapse pruning, cellular metabolism (mitochondrial respiratory chain complex I assembly), translation (ribosome assembly)
pdf(sprintf("%s/all_times_dynamics_compareCluster_GO_ORA_form_simp_Rel_0_8_compareClusterResult_spatial_GO_pval_figure_5.pdf", custom_cluster_comp_output_dir), height=14, width=10)
sample_ids <- c()
sample_id_times <- c()
spe_sub_go_pval_1 <- spe2speGOpval(spe_sub, "GO:0050776", all_times_GO_ORA_results2compare_df, "1")
sample_ids <- c(sample_ids, spe_sub_go_pval_1@int_metadata$imgData$sample_id)
sample_id_times <- c(sample_id_times, c("5", "10", "20", "40"))
spe_sub_go_pval_2 <- spe2speGOpval(spe_sub, "GO:0050804", all_times_GO_ORA_results2compare_df, "2")
sample_ids <- c(sample_ids, spe_sub_go_pval_2@int_metadata$imgData$sample_id)
sample_id_times <- c(sample_id_times, rep("", 4)) # only plot times for the first GO ID
spe_sub_go_pval_3 <- spe2speGOpval(spe_sub, "GO:0098883", all_times_GO_ORA_results2compare_df, "3")
sample_ids <- c(sample_ids, spe_sub_go_pval_3@int_metadata$imgData$sample_id)
sample_id_times <- c(sample_id_times, rep("", 4)) # only plot times for the first GO ID
spe_sub_go_pval_4 <- spe2speGOpval(spe_sub, "GO:0032981", all_times_GO_ORA_results2compare_df, "4")
sample_ids <- c(sample_ids, spe_sub_go_pval_4@int_metadata$imgData$sample_id)
sample_id_times <- c(sample_id_times, rep("", 4)) # only plot times for the first GO ID
spe_sub_go_pval_5 <- spe2speGOpval(spe_sub, "GO:0042255", all_times_GO_ORA_results2compare_df, "5")
sample_ids <- c(sample_ids, spe_sub_go_pval_5@int_metadata$imgData$sample_id)
sample_id_times <- c(sample_id_times, rep("", 4)) # only plot times for the first GO ID

spe_sub_go_pval <- cbind(spe_sub_go_pval_1, spe_sub_go_pval_2, spe_sub_go_pval_3, spe_sub_go_pval_4, spe_sub_go_pval_5)
names(sample_id_times) <- sample_ids
min_pval <- min(colData(spe_sub_go_pval)$GO_pvalue2plot, na.rm=TRUE)
max_pval <- max(colData(spe_sub_go_pval)$GO_pvalue2plot, na.rm=TRUE)

p <- plotVisium(spe_sub_go_pval, fill="GO_pvalue2plot", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
p <- p + scale_fill_gradientn(colours=c("blue", "grey50", "red"), na.value="white", values=scales::rescale(c(min_pval, -1, 0, 1, max_pval))) +
    theme(legend.position="bottom", legend.key.size=unit(1, "lines"), legend.title=element_text(size=8), legend.text=element_text(size=7)) 
## modify facet number of columns
p$facet$params$ncol <- 4
## modify facet titles
p$facet$params$labeller <- labeller(sample_id=sample_id_times)
## modify point size
p$layers[[length(p$layers)]]$aes_params$size <- 0.8
## modify legend title
p$guides$fill$title <- "-1^{0,1}*log10(p-value)"
print(p)
dev.off()





# other GO term sets
other_GO_output_dir <- sprintf("%s/%s/%s/custom/10-DEGs/%s/%s/sets", sub_markers_output_dir, qc_filtering_output_subdir, subclustering_integration, markers_test, markers_latent_variables_in_path)
if (! dir.exists(other_GO_output_dir)) {
    dir.create(other_GO_output_dir, recursive=TRUE, mode="0775")
}
all_times_other_go_terms <- list()
all_times_other_go_terms[["synapse_EI"]] <- c("excitatory postsynaptic potential", "modulation of excitatory postsynaptic potential", "positive regulation of excitatory postsynaptic potential", "inhibitory postsynaptic potential", "inhibitory synapse assembly")
all_times_other_go_terms[["neuron"]] <- grep("neuron", unique(all_times_GO_ORA_results2compare_df$Description), value=TRUE)
all_times_other_go_terms[["neurogenesis"]] <- grep("neurogenesis", unique(all_times_GO_ORA_results2compare_df$Description), value=TRUE)
all_times_other_go_terms[["axon"]] <- grep("axon", unique(all_times_GO_ORA_results2compare_df$Description), value=TRUE)
all_times_other_go_terms[["myelin"]]  <- grep("myelin", unique(all_times_GO_ORA_results2compare_df$Description), value=TRUE)
all_times_other_go_terms[["oligodendrocyte"]]  <- grep("oligodendrocyte", unique(all_times_GO_ORA_results2compare_df$Description), value=TRUE)
all_times_other_go_terms[["interleukin"]]  <- grep("interleukin", unique(all_times_GO_ORA_results2compare_df$Description), value=TRUE)
all_times_other_go_terms[["TNF"]]  <- grep("tumor necrosis factor", unique(all_times_GO_ORA_results2compare_df$Description), value=TRUE)
all_times_other_go_terms[["interferon"]]  <- grep("interferon", unique(all_times_GO_ORA_results2compare_df$Description), value=TRUE)



all_set_all_times_GO_ORA_results2compare_df <- data.frame()
for (one_set in names(all_times_other_go_terms)) {
    one_set_go_terms <- all_times_other_go_terms[[one_set]]
    # get ID corresponding to description
    one_set_go_ids <- c()
    for (one_go_term in one_set_go_terms) {
        one_set_go_ids <- c(one_set_go_ids, unique(all_times_GO_ORA_results2compare_df[which(all_times_GO_ORA_results2compare_df$Description == one_go_term), "ID"]))
    }
    # spatial plots
    one_set_custom_cluster_comp_spatial_output_name <- sprintf("all_times_%s_compareCluster_GO_ORA_form_simp_Rel_0_8_compareClusterResult_spatial_GO_pval", one_set)
    one_set_all_times_GO_ORA_results2compare_df <- go_pval_time_spatial_plots_2(spe_sub, all_times_GO_ORA_results2compare_df, one_set_go_ids, other_GO_output_dir, one_set_custom_cluster_comp_spatial_output_name)
    all_set_all_times_GO_ORA_results2compare_df <- rbind(all_set_all_times_GO_ORA_results2compare_df, cbind(set=rep(one_set, dim(one_set_all_times_GO_ORA_results2compare_df)[1]), one_set_all_times_GO_ORA_results2compare_df))
}
## write article supplementary table
write.table(all_set_all_times_GO_ORA_results2compare_df, file=sprintf("%s/all_times_all_sets_compareCluster_GO_ORA_form_compareClusterResult_spatial_GO_pval.tsv", other_GO_output_dir), sep="\t", quote=FALSE, row.names=FALSE)




















#######
# OLD #
#######

# GO term p-value spatial plots
## per time
for (one_time in levels(colData(spe_sub)$time)) {
    print(sprintf("time: %s", one_time))
    one_time_custom_cluster_comp_output_dir <- sprintf("%s/D%s", custom_cluster_comp_output_dir, one_time)
    one_time_custom_cluster_comp_output_name <- sprintf("D%s_top100_compareCluster_GO_ORA_form_simp_Rel_0_8_compareClusterResult_spatial_GO_pval", one_time)
    one_time_compareCluster_GO_ORA_simp_file <- sprintf("%s/DEGs_D%s_compareCluster_GO_ORA_form_simp_Rel_0_8_compareClusterResult.tsv", one_time_custom_cluster_comp_output_dir, one_time)
    one_time_compareCluster_GO_ORA_simp_df <- read.table(one_time_compareCluster_GO_ORA_simp_file, sep="\t", header=TRUE, quote="")
    go_ids2plot <- unique(one_time_compareCluster_GO_ORA_simp_df$ID)[1:100]
    go_pval_time_spatial_plots_2(spe_sub, all_times_GO_ORA_results2compare_df, go_ids2plot, mean_pval_time_cluster_category_count_df, one_time_custom_cluster_comp_output_dir, one_time_custom_cluster_comp_output_name)
}

## all times
all_times_go_ids <- c()
for (one_time in levels(colData(spe_sub)$time)) {
    print(sprintf("time: %s", one_time))
    one_time_custom_cluster_comp_output_dir <- sprintf("%s/D%s", custom_cluster_comp_output_dir, one_time)
    one_time_custom_cluster_comp_output_name <- sprintf("D%s_top100_compareCluster_GO_ORA_form_simp_Rel_0_8_compareClusterResult_spatial_GO_pval", one_time)
    one_time_compareCluster_GO_ORA_simp_file <- sprintf("%s/DEGs_D%s_compareCluster_GO_ORA_form_simp_Rel_0_8_compareClusterResult.tsv", one_time_custom_cluster_comp_output_dir, one_time)
    one_time_compareCluster_GO_ORA_simp_df <- read.table(one_time_compareCluster_GO_ORA_simp_file, sep="\t", header=TRUE, quote="")
    go_ids2plot <- unique(one_time_compareCluster_GO_ORA_simp_df$ID)[1:100]
    all_times_go_ids <- c(all_times_go_ids, go_ids2plot)
}
all_times_go_ids <- unique(all_times_go_ids)
all_times_go_ids_mean_pval_time_cluster_category_count_df <- mean_pval_time_cluster_category_count_df[which(mean_pval_time_cluster_category_count_df$ID %in% all_times_go_ids),]
all_times_go_ids_mean_pval_time_cluster_category_count_df_ordered <- all_times_go_ids_mean_pval_time_cluster_category_count_df[order(all_times_go_ids_mean_pval_time_cluster_category_count_df$mean_pval),]
go_ids2plot <- all_times_go_ids_mean_pval_time_cluster_category_count_df_ordered$ID[1:100]
all_times_custom_cluster_comp_output_name <- "all_times_top100_compareCluster_GO_ORA_form_simp_Rel_0_8_compareClusterResult_spatial_GO_pval"
go_pval_time_spatial_plots_2(spe_sub, all_times_GO_ORA_results2compare_df, go_ids2plot, mean_pval_time_cluster_category_count_df, custom_cluster_comp_output_dir, all_times_custom_cluster_comp_output_name)







### order terms according to mean p-value
mean_pval_df <- aggregate(all_times_GO_ORA_simp_results2compare_df$p.adjust, by=list(all_times_GO_ORA_simp_results2compare_df$ID, all_times_GO_ORA_simp_results2compare_df$Description), mean)
colnames(mean_pval_df) <- c("ID", "Description", "mean_pval")
mean_pval_df <- mean_pval_df[order(mean_pval_df$mean_pval),]
go_order_output_name <- sprintf("%s_mean_pval", custom_cluster_comp_output_name)
write.csv(mean_pval_df, file=sprintf("%s/%s.csv", custom_cluster_comp_output_dir, go_order_output_name), quote=FALSE, row.names=FALSE)
terms2plot_order <- mean_pval_df$ID
go_pval_time_spatial_plots(spe_sub, all_times_GO_ORA_simp_results2compare_df, terms2plot_order[1:100], mean_pval_df, "mean_pval", custom_cluster_comp_output_dir, go_order_output_name)


### order terms according to number of time and cluster combination number
time_cluster_category_count_df <- aggregate(all_times_all_clusters_GO_ORA_simp_df$p.adjust, by=list(all_times_all_clusters_GO_ORA_simp_df$ID, all_times_all_clusters_GO_ORA_simp_df$Description), length)
colnames(time_cluster_category_count_df) <- c("ID", "Description", "count")
time_cluster_category_count_df <- time_cluster_category_count_df[order(time_cluster_category_count_df$count, decreasing=TRUE),]
go_order_output_name <- sprintf("%s_count", subclustering_comp_output_name)
write.csv(time_cluster_category_count_df, file=sprintf("%s/%s.csv", custom_cluster_comp_output_dir, go_order_output_name), quote=FALSE, row.names=FALSE)
terms2plot_order <- time_cluster_category_count_df$ID
go_pval_time_spatial_plots(spe_sub, all_times_all_clusters_GO_ORA_simp_df, terms2plot_order[1:100], time_cluster_category_count_df, "count", custom_cluster_comp_output_dir, go_order_output_name)




# theme sets of selected GO terms based on custom cluster compareCluster() per time results
all_times_selected_go_terms <- list()
## D5
one_time <- "5"
one_time_selected_go_terms <- list()
one_time_selected_go_terms[["immune"]] <- c("immune effector process", "regulation of immune response", "leukocyte migration", "leukocyte cell-cell adhesion", "phagocytosis", "cytokine production", "antigen processing and presentation", "antigen processing and presentation of peptide", "leukocyte proliferation", "lymphocyte proliferation", "mononuclear cell proliferation", "regulation of lymphocyte proliferation", "antigen processing and presentation of peptide antigen via MHC class II", "humoral immune response", "regulation of cytokine production", "integrin-mediated signaling pathway", "innate immune response", "cytokine-mediated signaling pathway", "leukocyte homeostasis", "neutrophil homeostasis", "myeloid leukocyte migration", "interleukin-10 production", "regulation of interleukin-10 production", "phagocytosis, engulfment", "interleukin-1 production", "regulation of interleukin-1 production", "chemokine production", "regulation of chemokine production", "interferon-gamma production", "regulation of interferon-gamma production", "interleukin-8 production", "regulation of interleukin-8 production", "leukocyte apoptotic process")
one_time_selected_go_terms[["inflammation"]] <- c("inflammatory response", "production of molecular mediator involved in inflammatory response")
one_time_selected_go_terms[["synapse"]] <- c("modulation of chemical synaptic transmission", "regulation of trans−synaptic signaling", "regulation of membrane potential", "synapse organization", "endocytosis", "protein localization to synapse", "synapse pruning", "glutamate receptor signaling pathway", "ionotropic glutamate receptor signaling pathway", "positive regulation of peptidyl-tyrosine phosphorylation", "vesicle-mediated transport in synapse")
one_time_selected_go_terms[["response"]] <- c("response to wounding", "response to mechanical stimulus", "regeneration", "cellular response to external stimulus")
one_time_selected_go_terms[["translation"]] <- c("cytoplasmic translation", "ribosome assembly", "regulation of peptidase activity", "ribosomal small subunit assembly")
one_time_selected_go_terms[["neuron"]] <- c("axon ensheathment", "ensheathment of neurons", "axonogenesis", "positive regulation of neuron death", "apoptotic signaling pathway", "axon development", "dendritic cell differentiation")
one_time_selected_go_terms[["glia"]] <- c("gliogenesis", "glial cell activation")
one_time_selected_go_terms[["others"]] <- c("amyloid-beta clearance")
all_times_selected_go_terms[[one_time]] <- one_time_selected_go_terms
## D10
one_time <- "10"
one_time_selected_go_terms <- list()
one_time_selected_go_terms[["mitochondrion"]] <- c("aerobic respiration", "oxidative phosphorylation", "generation of precursor metabolites and energy", "NADH dehydrogenase complex assembly", "mitochondrial respiratory chain complex I assembly", "proton transmembrane transport", "mitochondrion organization", "mitochondrial respiratory chain complex assembly")
one_time_selected_go_terms[["immune"]] <- c("immune effector process", "antigen processing and presentation of peptide antigen", "humoral immune response", "positive regulation of immune system process", "B cell mediated immunity", "immunoglobulin mediated immune response", "antigen processing and presentation of exogenous peptide antigen", "innate immune response", "phagocytosis", "regulation of cytokine production", "neutrophil homeostasis")
one_time_selected_go_terms[["inflammation"]] <- c("inflammatory response", "inflammatory response to wounding")
one_time_selected_go_terms[["translation"]] <- c("cytoplasmic translation", "regulation of peptidase activity", "ribosomal small subunit assembly")
one_time_selected_go_terms[["synapse"]] <- c("regulation of membrane potential", "modulation of chemical synaptic transmission", "regulation of trans-synaptic signaling", "synapse pruning", "synapse organization", "endocytosis", "synaptic vesicle cycle")
one_time_selected_go_terms[["glia"]] <- c("gliogenesis")
one_time_selected_go_terms[["response"]] <- c("response to wounding")
one_time_selected_go_terms[["neuron"]] <- c("neural nucleus development", "positive regulation of neuron death", "myelination")
all_times_selected_go_terms[[one_time]] <- one_time_selected_go_terms
## D20
one_time <- "20"
one_time_selected_go_terms <- list()
one_time_selected_go_terms[["immune"]] <- c("immune effector process", "cytokine production", "regulation of immune response", "antigen processing and presentation", "phagocytosis", "antigen processing and presentation of peptide antigen", "antigen processing and presentation of peptide antigen via MHC class II", "extracellular matrix organization", "extracellular structure organization", "cytokine-mediated signaling pathway", "leukocyte homeostasis", "chemokine production", "regulation of chemokine production", "regulation of cytokine production", "regulation of leukocyte apoptotic process", "leukocyte apoptotic process")
one_time_selected_go_terms[["inflammation"]] <- c("inflammatory response", "blood vessel morphogenesis", "blood vessel development", "production of molecular mediator involved in inflammatory response")
one_time_selected_go_terms[["translation"]] <- c("cytoplasmic translation", "regulation of peptidase activity")
one_time_selected_go_terms[["response"]] <- c("response to mechanical stimulus", "regeneration")
one_time_selected_go_terms[["mitochondrion"]] <- c("oxidative phosphorylation", "NADH dehydrogenase complex assembly", "mitochondrial respiratory chain complex I assembly", "generation of precursor metabolites and energy", "mitochondrial respiratory chain complex assembly", "NADH regeneration", "canonical glycolysis", "glucose catabolic process to pyruvate", "response to oxygen levels")
one_time_selected_go_terms[["synapse"]] <- c("endocytosis", "positive regulation of peptidyl-tyrosine phosphorylation", "modulation of chemical synaptic transmission", "regulation of membrane potential", "regulation of trans-synaptic signaling, synapse organization")
one_time_selected_go_terms[["glia"]] <- c("glial cell differentiation", "gliogenesis")
one_time_selected_go_terms[["neuron"]] <- c("regulation of neuron death")
all_times_selected_go_terms[[one_time]] <- one_time_selected_go_terms
## D40
one_time <- "40"
one_time_selected_go_terms <- list()
one_time_selected_go_terms[["translation"]] <- c("cytoplasmic translation", "ribosome biogenesis", "ribosome assembly", "ribosomal small subunit biogenesis", "ribosomal small subunit assembly", "translational elongation", "positive regulation of translation", "protein depolymerization")
one_time_selected_go_terms[["mitochondrion"]] <- c("oxidative phosphorylation", "NADH dehydrogenase complex assembly", "mitochondrial respiratory chain complex I assembly", "proton transmembrane transport", "mitochondrial respiratory chain complex assembly", "potassium ion transmembrane transport", "electron transport chain", "generation of precursor metabolites and energy")
one_time_selected_go_terms[["synapse"]] <- c("modulation of chemical synaptic transmission", "regulation of trans-synaptic signaling", "synapse organization", "regulation of membrane potential, protein localization to synapse", "excitatory postsynaptic potential", "synaptic vesicle cycle", "calcium-mediated signaling", "vesicle-mediated transport in synapse", "second-messenger-mediated signaling", "synapse maturation", "regulation of neurotransmitter receptor activity", "regulation of transmembrane transport", "maintenance of synapse structure", "regulation of synapse maturation", "calcium ion transport", "receptor localization to synapse", "calcium ion export across plasma membrane", "spontaneous synaptic transmission", "glutamate receptor signaling pathway")
one_time_selected_go_terms[["neuron"]] <- c("dendrite development", "axon ensheathment", "ensheathment of neurons", "axo-dendritic transport", "axon development", "axonogenesis", "neuron migration")
one_time_selected_go_terms[["immune"]] <- c("antigen processing and presentation of exogenous antigen", "humoral immune response", "antigen processing and presentation of exogenous peptide antigen", "lymphocyte proliferation")
one_time_selected_go_terms[["glia"]] <- c("oligodendrocyte differentiation", "astrocyte development")
one_time_selected_go_terms[["others"]] <- c("liver regeneration")
one_time_selected_go_terms[["inflammation"]] <- c("inflammatory response to antigenic stimulus", "inflammatory response")
all_times_selected_go_terms[[one_time]] <- one_time_selected_go_terms

## per time
for (one_time in names(all_times_selected_go_terms)) {
    print(sprintf("time: %s", one_time))
    one_time_selected_go_terms <- all_times_selected_go_terms[[one_time]]
    one_time_spatial_plot_output_dir <- sprintf("%s/D%s/spatial_plots/", custom_cluster_comp_output_dir, one_time)
    if (! dir.exists(one_time_spatial_plot_output_dir)) {
        dir.create(one_time_spatial_plot_output_dir, recursive=TRUE, mode="0775")
    }

    for (one_theme in names(one_time_selected_go_terms)) {
        print(sprintf("theme: %s", one_theme))
        one_time_one_theme_custom_cluster_comp_spatial_output_name <- sprintf("D%s_%s_compareCluster_GO_ORA_form_simp_Rel_0_8_compareClusterResult_spatial_GO_pval", one_time, one_theme)
        # get ID corresponding to description
        selected_go_ids <- c()
        for (one_go_term in one_time_selected_go_terms[[one_theme]]) {
            selected_go_ids <- c(selected_go_ids, unique(all_times_GO_ORA_simp_results2compare_df[which(all_times_GO_ORA_simp_results2compare_df$Description == one_go_term), "ID"]))
        }
        # spatial plots
        go_pval_time_spatial_plots_2(spe_sub, all_times_GO_ORA_results2compare_df, selected_go_ids, mean_pval_time_cluster_category_count_df, one_time_spatial_plot_output_dir, one_time_one_theme_custom_cluster_comp_spatial_output_name)
    }
}

## all times
### get all themes
all_themes <- c()
for (one_time in names(all_times_selected_go_terms)) {
    all_themes <- c(all_themes, names(all_times_selected_go_terms[[one_time]]))
}
all_themes <- unique(all_themes)
for (one_theme in all_themes) {
    print(sprintf("theme: %s", one_theme))
    #### get all terms across all times corresponding to one theme
    one_theme_go_terms <- c()
    for (one_time in names(all_times_selected_go_terms)) {
        one_time_selected_go_terms <- all_times_selected_go_terms[[one_time]]
        if (one_theme %in% names(one_time_selected_go_terms)) {
            one_theme_go_terms <- c(one_theme_go_terms, one_time_selected_go_terms[[one_theme]])
        }
    }
    one_theme_go_terms <- unique(one_theme_go_terms)
    print(one_theme_go_terms)
    
    one_theme_custom_cluster_comp_spatial_output_name <- sprintf("all_times_%s_compareCluster_GO_ORA_form_simp_Rel_0_8_compareClusterResult_spatial_GO_pval", one_theme)
    # get ID corresponding to description
    one_theme_go_ids <- c()
    for (one_go_term in one_theme_go_terms) {
        one_theme_go_ids <- c(one_theme_go_ids, unique(all_times_GO_ORA_simp_results2compare_df[which(all_times_GO_ORA_simp_results2compare_df$Description == one_go_term), "ID"]))
    }
    # spatial plots
    go_pval_time_spatial_plots_2(spe_sub, all_times_GO_ORA_results2compare_df, one_theme_go_ids, mean_pval_time_cluster_category_count_df, custom_cluster_comp_output_dir, one_theme_custom_cluster_comp_spatial_output_name)
}




go_pval_time_spatial_plots <- function(spe_obj, go_df, go_terms, stat_df, stat_col, out_dir, out_name) {
    pdf(sprintf("%s/%s.pdf", out_dir, out_name), width=12)
    for (one_go_id in go_terms) {
        print(sprintf("GO ID: %s", one_go_id))
        one_go_id_stat <- stat_df[which(stat_df$ID == one_go_id), stat_col]
        one_go_id_description <- stat_df[which(stat_df$ID == one_go_id), "Description"]
        
        one_go_id_go_df <- go_df[which(go_df$ID == one_go_id),]
        one_go_id_go_df <- droplevels(one_go_id_go_df)
        go_term_pval_df <- data.frame()
        # get p-value per time and cluster
        for (one_time in levels(one_go_id_go_df$time)) {
            one_time_one_go_id_go_df <- one_go_id_go_df[which(one_go_id_go_df$time == one_time),]
            one_time_one_go_id_go_df <- droplevels(one_time_one_go_id_go_df)
            for (one_cluster in levels(one_time_one_go_id_go_df$cluster)) {
                one_cluster_one_time_one_go_id_go_df <- one_time_one_go_id_go_df[which(one_time_one_go_id_go_df$cluster == one_cluster),]
                if (dim(one_cluster_one_time_one_go_id_go_df)[1] == 2) {
                    # significant p-value for both down and up regualtion: take the most significant p-value
                    min_pvalue <- min(one_cluster_one_time_one_go_id_go_df$p.adjust)
                    one_pvalue <- one_cluster_one_time_one_go_id_go_df[which(one_cluster_one_time_one_go_id_go_df$p.adjust == min_pvalue), "p.adjust"]
                    one_category <- one_cluster_one_time_one_go_id_go_df[which(one_cluster_one_time_one_go_id_go_df$p.adjust == min_pvalue), "category"]
                } else {
                    if (dim(one_cluster_one_time_one_go_id_go_df)[1] == 1) {
                        # significant p-value only for down or up regualtion
                        one_pvalue <- one_cluster_one_time_one_go_id_go_df$p.adjust
                        one_category <- one_cluster_one_time_one_go_id_go_df$category
                    }
                }
                # get barcodes corresponding to time and subcluster
                barcodes <- rownames(colData(spe_obj)[which(colData(spe_obj)$time == one_time & colData(spe_obj)$seurat_custom_clusters == one_cluster),])
                go_term_pval_df <- rbind(go_term_pval_df, data.frame(barcode=barcodes, GO_pvalue=rep(one_pvalue, length(barcodes)), category=rep(one_category, length(barcodes))))
            }
        }
        # log10
        go_term_pval_df$GO_pvalue <- -log(go_term_pval_df$GO_pvalue, 10)
        # log10(p-value) to plot: negative for downregulation, positive for upregulation
        go_term_pval_df$GO_pvalue2plot <- go_term_pval_df$GO_pvalue * ifelse(go_term_pval_df$category == "down", -1, 1)
        na_barcodes <- all_barcodes[! all_barcodes %in% go_term_pval_df$barcode]
        na_barcodes_nb <- length(na_barcodes)
        go_term_pval_df <- rbind(go_term_pval_df, data.frame(barcode=na_barcodes, GO_pvalue=rep(NA, na_barcodes_nb), category=rep(NA, na_barcodes_nb), GO_pvalue2plot=rep(NA, na_barcodes_nb)))
        ## add 'GO_pvalue' column to scolData(spe_obj): order go_term_pval_df
        rownames(go_term_pval_df) <- go_term_pval_df$barcode
        go_term_pval_df <- go_term_pval_df[rownames(colData(spe_obj)),]
        colData(spe_obj)$GO_pvalue2plot <- go_term_pval_df$GO_pvalue2plot
        
        if (stat_col == "mean_pval") {
            plot_title <- sprintf("%s, %s, mean p-value: %.2e", one_go_id, one_go_id_description, one_go_id_stat)
        } else {
            if (stat_col == "count") {
                plot_title <- sprintf("%s, %s, time.subcluser number: %d", one_go_id, one_go_id_description, one_go_id_stat)
            }
        }

        #p <- plotSpots(spe_obj, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="GO_pvalue2plot", size=0.3)
        #p <- p + labs(title=plot_title) +
        #    theme(legend.position="bottom")
        ## modify facet number of columns
        #p$facet$params$ncol <- 4
        ## modify facet titles
        #p$facet$params$labeller <- labeller(sample_id=sample_id_times)
        #print(p)

        p <- plotVisium(spe_obj, fill="GO_pvalue2plot", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
        p <- p + scale_fill_gradient2(name="log10(p-value)", low="blue", high="red", mid="grey50", na.value="white") +
        #p <- p + scale_fill_viridis_c(na.value="white") +
            labs(title=plot_title) +
            theme(legend.position="bottom", legend.key.size=unit(1, "lines"))
            #guides(fill=guide_legend(title="log10(p-value)"))
        ## modify facet number of columns
        p$facet$params$ncol <- 4
        ## modify facet titles
        p$facet$params$labeller <- labeller(sample_id=sample_id_times)
        ## modify point size
        p$layers[[5]]$aes_params$size <- 1
        ## modify legend title
        p$labels$fill <- "log10(p-value)"
        print(p)
    }
    dev.off()
}


