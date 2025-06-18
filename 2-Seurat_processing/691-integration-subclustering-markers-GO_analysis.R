.libPaths(c(.libPaths(), "/home/christophe.lepriol/NeuroDev_ADD/R/r_4.2.0"))
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(GOSemSim)


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
sct_variables_to_regress <- "qc_mito_percent_RNA"
scale_data_variables_to_regress <- "qc_mito_percent_RNA"
variables_to_regress_in_path <- "mito"
scale_data_split_by <- "orig.ident"
scale_data_split_by_in_path <- "orig"
integration_method <- "Seurat"
integration_group <- "orig.ident"
integration_group_by_in_path <- "orig"
markers_test <- "negbinom"
markers_latent_variables <- "orig.ident"

# parameter settings for subclusterings
subclustering_integration <- "Seurat"
## set subclusterings parameter settings
subsclusterings_parameter_settings <- list()
### cluster 0
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
### cluster 1
subclustering_name <- "1"
cluster_1_parameters <- list()
cluster_1_parameters[["Seurat"]]  <- subclustering_name
cluster_1_parameters[["normalization"]] <- "SCTransform_v2"
cluster_1_parameters[["features"]]  <- 1500
cluster_1_parameters[["integration"]] <- "Seurat"
cluster_1_parameters[["dimensions"]]  <- 10
cluster_1_parameters[["resolution"]] <- 0.03
cluster_1_parameters[["subclusters"]] <- 4
subsclusterings_parameter_settings[[subclustering_name]] <- cluster_1_parameters
### cluster 2
subclustering_name <- "2"
cluster_2_parameters <- list()
cluster_2_parameters[["Seurat"]]  <- subclustering_name
cluster_2_parameters[["normalization"]] <- "SCTransform_v2"
cluster_2_parameters[["features"]]  <- 250
cluster_2_parameters[["integration"]] <- "Seurat"
cluster_2_parameters[["dimensions"]]  <- 15
cluster_2_parameters[["resolution"]] <- 0.5
cluster_2_parameters[["subclusters"]] <- 11
subsclusterings_parameter_settings[[subclustering_name]] <- cluster_2_parameters
### cluster 3
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
### cluster 4
subclustering_name <- "4"
cluster_4_parameters <- list()
cluster_4_parameters[["Seurat"]]  <- subclustering_name
cluster_4_parameters[["normalization"]] <- "SCTransform_v2"
cluster_4_parameters[["features"]]  <- 500
cluster_4_parameters[["integration"]] <- "Seurat"
cluster_4_parameters[["dimensions"]]  <- 10
cluster_4_parameters[["resolution"]] <- 0.08
cluster_4_parameters[["subclusters"]] <- 2
subsclusterings_parameter_settings[[subclustering_name]] <- cluster_4_parameters

# analysis to perform
## data
whole_clustering_analysis <- TRUE
subclustering_analysis <- TRUE
## methods
go_analysis <- TRUE
gsea <- FALSE
go_comparison_analysis <- TRUE
similarity_method <- "Rel"
similarity_cutoff <- 0.8


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

# output directory
output_dir <- sprintf("%s/20-Count_analysis/10-EpiReg/output/%s/TSO_polyA_R1trim1_ov5_n2_min20/%s", work_dir, input_dirname, genome_name)
if (visium_params_in_path) {
    output_dir <- sprintf("%s/Visium_params", output_dir) 
}
output_dirname <- ifelse(is.na(timepoint), "all", sprintf("D%d", timepoint))
output_dir <- sprintf("%s/20-Integration/%s", output_dir, output_dirname)
sub_markers_output_dir <- sprintf("%s/sub_markers", output_dir)
## quality control
qc_filtering_output_subdir <- sprintf("Norm_%s-scale_%s-int_%s/QC-sum%d_det%s_mito%d_hb%d/10-QC_filtering", variables_to_regress_in_path, scale_data_split_by_in_path, integration_group_by_in_path, spot_sum_of_counts_threshold, spot_detected_genes_threshold, spot_mito_gene_count_pct_threshold, spot_hb_gene_count_pct_threshold)

# preprare GO data for measuring semantic similarity for both GO analysis and GO term comparison analysis
rnGO <- godata(OrgDb="org.Rn.eg.db", ont="BP")

# whole clustering
if (whole_clustering_analysis) {
    ## parameter setting: get whole clustering parameter evaluation output
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
    whole_nb_clusters <- 5
    print(sprintf("normalization: %s, feature number: %d, integration: %s, dimension number: %d, resolution: %.2f", whole_normalization, whole_features, whole_integration, whole_dimensions, whole_resolution))
    
    ## build path to per-timepoint DE analysis
    whole_clustering_output_dir <- sprintf("%s/%s/%s/whole", sub_markers_output_dir, qc_filtering_output_subdir, whole_integration)
    method_dir_index <- ifelse(whole_normalization == "SCTransform", "00", "01")
    normalization_output_dir <- sprintf("%s/%s-%s/%s/featnb%s", whole_clustering_output_dir, method_dir_index, whole_normalization, vars2regress_in_dirname, whole_features)
    resolution_output_dir <- sprintf("%s/%s/%d_dims/res%s", normalization_output_dir, gsub(" ", "_", whole_integration), whole_dimensions, sub("[.]", "_", whole_resolution))
    DEG_output_dir <- sprintf("%s/10-DEGs/%s/%s", resolution_output_dir, markers_test, markers_latent_variables_in_path)
    ## build per-timepoint DE analysis output basename
    DEG_output_name <- sprintf("%s_featnb%d_%s_%ddims_res%s_DEGs", whole_normalization, whole_features, gsub(" ", "_", whole_integration), whole_dimensions, sub("[.]", "_", whole_resolution))    
    ## GO analysis for all clusters and all times and GO term comparison per time and cluster
    GO_analysis_all_clusters_all_times(epireg_visium_metadata_df, whole_nb_clusters, go_analysis, similarity_method, rnGO, similarity_cutoff, gsea, go_comparison_analysis, DEG_output_dir, DEG_output_name)
}

# subclusterings
if (subclustering_analysis) {
    for (one_cluster in names(subsclusterings_parameter_settings)) {
        ## set subclustering parameters
        sub_normalization <- subsclusterings_parameter_settings[[one_cluster]][["normalization"]]
        sub_features <- subsclusterings_parameter_settings[[one_cluster]][["features"]]
        sub_integration <- subsclusterings_parameter_settings[[one_cluster]][["integration"]]
        sub_dimensions <- subsclusterings_parameter_settings[[one_cluster]][["dimensions"]]
        sub_resolution <- subsclusterings_parameter_settings[[one_cluster]][["resolution"]]
        sub_nb_clusters <- subsclusterings_parameter_settings[[one_cluster]][["subclusters"]]

        ## build path to per-timepoint DE analysis
        subclustering_output_dir <- sprintf("%s/%s/%s/sub/c%s", sub_markers_output_dir, qc_filtering_output_subdir, subclustering_integration, one_cluster)
        method_dir_index <- ifelse(sub_normalization == "SCTransform", "00", "01")
        normalization_output_dir <- sprintf("%s/%s-%s/%s/featnb%s", subclustering_output_dir, method_dir_index, sub_normalization, vars2regress_in_dirname, sub_features)
        resolution_output_dir <- sprintf("%s/%s/%d_dims/res%s", normalization_output_dir, gsub(" ", "_", sub_integration), sub_dimensions, sub("[.]", "_", sub_resolution))
        DEG_output_dir <- sprintf("%s/10-DEGs/%s/%s", resolution_output_dir, markers_test, markers_latent_variables_in_path)

        ## build per-timepoint DE analysis output basename
        DEG_output_name <- sprintf("%s_featnb%d_%s_%ddims_res%s_DEGs", sub_normalization, sub_features, gsub(" ", "_", sub_integration), sub_dimensions, sub("[.]", "_", sub_resolution))

        ## GO analysis for all clusters and all times and GO term comparison per time and cluster
        GO_analysis_all_clusters_all_times(epireg_visium_metadata_df, sub_nb_clusters, go_analysis, similarity_method, rnGO, similarity_cutoff, gsea, go_comparison_analysis, DEG_output_dir, DEG_output_name)
    }
}

# GO term comparison for clusters 1, 2 and 4 and subclusters of clusters 0 and 3: time comparison per (sub)cluster and (sub)cluster comparison per time
all_times_DE_results2compare_df <- data.frame()
all_times_GO_results2compare_df <- data.frame()
## get whole clustering DE analysis results for clusters 1, 2, and 4
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
DEG_output_dir <- sprintf("%s/10-DEGs/%s/%s", resolution_output_dir, markers_test, markers_latent_variables_in_path)
### build per-timepoint DE analysis output basename
DEG_output_name <- sprintf("%s_featnb%d_%s_%ddims_res%s_DEGs", whole_normalization, whole_features, gsub(" ", "_", whole_integration), whole_dimensions, sub("[.]", "_", whole_resolution))
### read whole clustering DE analysis results
all_times_all_clusters_DE_results_df <- read.csv(sprintf("%s/comp/%s_all_times_all_clusters.csv", DEG_output_dir, DEG_output_name))
all_times_all_clusters_DE_results_df$time <- as.factor(all_times_all_clusters_DE_results_df$time)
all_times_all_clusters_DE_results_df$cluster <- as.factor(all_times_all_clusters_DE_results_df$cluster)
#### only DE analysis results for clusters 1, 2 and 4
all_times_clusters_124_DE_results_df <- all_times_all_clusters_DE_results_df[which(all_times_all_clusters_DE_results_df$cluster %in% c("1", "2", "4")),]
all_times_clusters_124_DE_results_df <- droplevels(all_times_clusters_124_DE_results_df)
all_times_DE_results2compare_df <- rbind(all_times_DE_results2compare_df, all_times_clusters_124_DE_results_df)

### read whole clustering GO analysis results
all_times_all_clusters_GO_results_df <- read.table(sprintf("%s/comp/%s_all_times_all_clusters_GO_ORA.tsv", DEG_output_dir, DEG_output_name), header=TRUE, sep="\t", quote="")
all_times_all_clusters_GO_results_df$time <- as.factor(all_times_all_clusters_GO_results_df$time)
all_times_all_clusters_GO_results_df$cluster <- as.factor(all_times_all_clusters_GO_results_df$cluster)
all_times_all_clusters_GO_results_df$category <- as.factor(all_times_all_clusters_GO_results_df$category)
#### only DE analysis results for clusters 1, 2 and 4
all_times_clusters_124_GO_results_df <- all_times_all_clusters_GO_results_df[which(all_times_all_clusters_GO_results_df$cluster %in% c("1", "2", "4")),]
all_times_clusters_124_GO_results_df <- droplevels(all_times_clusters_124_GO_results_df)
all_times_GO_results2compare_df <- rbind(all_times_GO_results2compare_df, all_times_clusters_124_GO_results_df)

## get subclustering DE analysis results for clusters 0 and 3
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
    DEG_output_dir <- sprintf("%s/10-DEGs/%s/%s", resolution_output_dir, markers_test, markers_latent_variables_in_path)
    ### build per-timepoint DE analysis output basename
    DEG_output_name <- sprintf("%s_featnb%d_%s_%ddims_res%s_DEGs", sub_normalization, sub_features, gsub(" ", "_", sub_integration), sub_dimensions, sub("[.]", "_", sub_resolution))
    ### read whole clustering DE analysis results
    all_times_all_clusters_DE_results_df <- read.csv(sprintf("%s/comp/%s_all_times_all_clusters.csv", DEG_output_dir, DEG_output_name))
    all_times_all_clusters_DE_results_df$time <- as.factor(all_times_all_clusters_DE_results_df$time)
    #### update cluster column to the format 'cluster.subcluster'
    all_times_all_clusters_DE_results_df$cluster <- sprintf("%s_%s", one_cluster, all_times_all_clusters_DE_results_df$cluster) # compareCluster does not handle the 'cluster.subcluster' format correctly (the 'cluster' column is only set to cluster value and the 'category' column is set to subcluster value) --> use 'cluster_subscluster' format instead
    all_times_all_clusters_DE_results_df$cluster <- as.factor(all_times_all_clusters_DE_results_df$cluster)
    all_times_DE_results2compare_df <- rbind(all_times_DE_results2compare_df, all_times_all_clusters_DE_results_df)
    
    ### read whole clustering GO analysis results
    all_times_all_clusters_GO_results_df <- read.table(sprintf("%s/comp/%s_all_times_all_clusters_GO_ORA.tsv", DEG_output_dir, DEG_output_name), header=TRUE, sep="\t", quote="")
    all_times_all_clusters_GO_results_df$time <- as.factor(all_times_all_clusters_GO_results_df$time)
    #### update cluster column to the format 'cluster.subcluster'
    all_times_all_clusters_GO_results_df$cluster <- sprintf("%s.%s", one_cluster, all_times_all_clusters_GO_results_df$cluster)
    all_times_all_clusters_GO_results_df$cluster <- as.factor(all_times_all_clusters_GO_results_df$cluster)
    all_times_GO_results2compare_df <- rbind(all_times_GO_results2compare_df, all_times_all_clusters_GO_results_df)
}

markers_latent_variables_in_path <- ifelse(length(markers_latent_variables > 1), paste(gsub("[.]", "", markers_latent_variables), collapse="-"), ifelse(is.na(markers_latent_variables), "no_var", gsub("[.]", "", markers_latent_variables)))
go_term_comparison_output_dir <- sprintf("%s/%s/%s/custom/10-DEGs/%s/%s", sub_markers_output_dir, qc_filtering_output_subdir, subclustering_integration, markers_test, markers_latent_variables_in_path)
if (! dir.exists(go_term_comparison_output_dir)) {
    dir.create(go_term_comparison_output_dir, recursive=TRUE, mode="0775")
}
write.csv(all_times_DE_results2compare_df, file=sprintf("%s/all_times_all_clusters_DE_analysis.csv", go_term_comparison_output_dir), quote=FALSE, row.names=FALSE)
write.table(all_times_GO_results2compare_df, file=sprintf("%s/all_times_all_clusters_GO_ORA_analysis.tsv", go_term_comparison_output_dir), sep="\t", quote=FALSE, row.names=FALSE)

## DEG counts per time and cluster (custom clusters)
all_times_DE_results2compare_only_DEGs_df <- all_times_DE_results2compare_df[which(all_times_DE_results2compare_df$p_val_adj < 0.05),]
up_DEG_count_df <- aggregate(all_times_DE_results2compare_only_DEGs_df$avg_log2FC, by=list(all_times_DE_results2compare_only_DEGs_df$time, all_times_DE_results2compare_only_DEGs_df$cluster), function(x) { length(x[x > 0]) })
colnames(up_DEG_count_df) <- c("time", "cluster", "count")
up_DEG_count_df <- cbind(category=rep("up", dim(up_DEG_count_df)[1]), up_DEG_count_df)
up_DEG_count_df$category <- as.factor(up_DEG_count_df$category)
up_DEG_count_df$cluster <- as.factor(up_DEG_count_df$cluster)
down_DEG_count_df <- aggregate(all_times_DE_results2compare_only_DEGs_df$avg_log2FC, by=list(all_times_DE_results2compare_only_DEGs_df$time, all_times_DE_results2compare_only_DEGs_df$cluster), function(x) { length(x[x < 0]) })
colnames(down_DEG_count_df) <- c("time", "cluster", "count")
down_DEG_count_df <- cbind(category=rep("down", dim(down_DEG_count_df)[1]), down_DEG_count_df)
down_DEG_count_df$category <- as.factor(down_DEG_count_df$category)
down_DEG_count_df$cluster <- as.factor(down_DEG_count_df$cluster)
DE_count_df <- rbind(up_DEG_count_df, down_DEG_count_df)
write.csv(DE_count_df, file=sprintf("%s/all_times_all_clusters_DEG_counts.csv", go_term_comparison_output_dir), quote=FALSE, row.names=FALSE)

## GO term counts per time and cluster (custom clusters)
GO_count_df <- aggregate(all_times_GO_results2compare_df$ID, by=list(all_times_GO_results2compare_df$time, all_times_GO_results2compare_df$cluster, all_times_GO_results2compare_df$category), length)
colnames(GO_count_df) <- c("time", "cluster", "category", "count")
GO_count_df <- GO_count_df[, colnames(DE_count_df)]
write.csv(GO_count_df, file=sprintf("%s/all_times_all_clusters_GO_ORA_counts.csv", go_term_comparison_output_dir), quote=FALSE, row.names=TRUE)

### DEG and GO term counts per time and cluster (custom clusters) plot
levels(DE_count_df$cluster) <- sub("_", ".", levels(DE_count_df$cluster))
DE_GO_count_df <- rbind(cbind(category2=rep("DE", dim(DE_count_df)[1]), DE_count_df), cbind(category2=rep("GO", dim(GO_count_df)[1]), GO_count_df))
pdf(sprintf("%s/all_times_all_clusters_DEG_GO_ORA_counts.pdf", go_term_comparison_output_dir))
p <- ggplot(DE_GO_count_df, aes(x=time, y=count, fill=cluster)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(category2~category, scales="free_y") +
    labs(x="Time", y="Count", fill="Cluster") +
    theme_bw() +
    theme(legend.position="bottom") +
    theme(panel.border=element_rect(color="grey50")) +
    guides(fill=guide_legend(nrow=1, byrow=TRUE))
print(p)
dev.off()

## GO term comparisons
GO_term_comparison_analysis(all_times_DE_results2compare_df, similarity_method, rnGO, similarity_cutoff, gsea, go_term_comparison_output_dir, "DEGs")

## build article supplementary table
all_times_GO_term_comparison_df <- data.frame()
for (one_time in levels(epireg_visium_metadata_df$time)) {
    one_time_output_dir <- sprintf("%s/D%s", go_term_comparison_output_dir, one_time)
    one_time_output_name <- sprintf("DEGs_D%s_compareCluster_GO_ORA_form_simp_%s_%s", one_time, similarity_method, sub("[.]", "_", similarity_cutoff))
    one_time_GO_term_comparison_df <- read.table(sprintf("%s/%s_compareClusterResult.tsv", one_time_output_dir, one_time_output_name), header=TRUE, sep="\t", quote="")
    all_times_GO_term_comparison_df <- rbind(all_times_GO_term_comparison_df, cbind(time=rep(one_time, dim(one_time_GO_term_comparison_df)[1]), one_time_GO_term_comparison_df))
}
write.table(all_times_GO_term_comparison_df, file=sprintf("%s/all_times_compareCluster_GO_ORA_form_simp_%s_%s_compareClusterResult.tsv", go_term_comparison_output_dir, similarity_method, sub("[.]", "_", similarity_cutoff)), sep="\t", quote=FALSE, row.names=FALSE)

