.libPaths(c("/home/christophe.lepriol/NeuroDev_ADD/R/r_4.1.0", .libPaths()))
library(Seurat)
library(scater)
library(patchwork)
library(scran) # for buildSNNGraph() function
library(intrinsicDimension) # for maxLikGlobalDimEst() function
library(bluster) # for bootstrapStability() function
library(tidyr) # for pivot_longer() function
library(tibble) # for rownames_to_column() function
library(cluster) # silhouette() function
library(dplyr) # for mutate() function
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(pheatmap)
#library(ggrepel)


##############
# Parameters #
##############
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"
src_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/src", work_dir)
genome_name <- "NCBIRefSeq108_NCBIRefSeq108GTF"
#sample_vector <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")
#sample <- "A_L1_S1"


#############
# Functions #
#############
source(sprintf("%s/20-pipelines-Seurat-functions.R", src_dir))


############
# Analysis #
############

# multiple parameter evaluation
## SCTransform
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"
genome_name <- "NCBIRefSeq108_NCBIRefSeq108GTF"
#sample_vector <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")
sample_vector <- c("D_L4_S16")
# pipeline parameters
## QC
lib_size_threshold <- 2000
exp_genes_threshold <- 1000
## normalization
#normalization_method <- "SCTransform" # either "SCTransform" or normalization.method parameter of NormalizeData() function, one of the following values "LogNormalize", "CLR", "RC"
normalization_method_vector <- c("SCTransform", "LogNormalize")
## clustering
feature_selection_method <- "HVGs" # either "HVGs" or "SVGs"
#features_nb <- 2000
features_nb_vector <- c(500, 1000, 3000, 5000)
feature_residual_variance_threshold <- NA
#dimensions_nb_param <- "estimation" # either "estimation" or number of PCs
dimensions_nb_param_vector <- c(10, 15, 20)
#nearest_neighbors_nb <- 20 # k.param parameter of FindNeighbors() function
nearest_neighbors_nb_vector <- c(20, 30, 40)
clustering_resolution_vector <- c(0.6, 0.8, 1)
#clustering_method <- "Louvain" # algorithm parameter of FindClusters() function, one of the following values: "Louvain", ...
clustering_method_vector <- c(1) # algorithm parameter of FindClusters() function, one of the following values: "Louvain", ...
cluster_analysis <- TRUE
cluster_analysis_spot_pct <- 0.8
cluster_analysis_iterations <- 20
## marker genes
markers_test <- "wilcox"  # test.use parameter of FindAllMarkers() function, one of the following values: "wilcox", ...
min_logFC <- 0.25 # logfc.threshold parameter of FindAllMarkers() function
min_pct_spots <- 0.1 # min.pct parameter of FindAllMarkers() function

all_samples_df2ggplot <- data.frame()
all_samples_top_df <- data.frame()
#output_basename <- sprintf("%s_%s_%s", normalization_output[["basename"]], feature_selection_method, clustering_output[["basename"]])
#plot_title <- sprintf("Normalization: %s\nFeatures: %s\nClustering: PCs=%s, %dNN, %s", normalization_method, feature_selection_method, dimensions_nb_param, nearest_neighbors_nb, clustering_method)
for (sample in sample_vector) {
    print(sprintf("sample: %s", sample))
    output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s", work_dir, genome_name, sample)
    if (! dir.exists(output_dir)) {
        dir.create(output_dir, recursive=TRUE, mode="0775")
    }
    sample_df2ggplot <- data.frame()
    sample_plot_list <- list()

    # Load data
    library(SpatialExperiment)
    # expression matrices
    ## Space Ranger
    #h5_file <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs/raw_feature_bc_matrix.h5", work_dir, genome_name, sample)
    #h5_data <- Read10X_h5(h5_file)
    ## ST Pipeline
    st_pipeline_matrix_file <- sprintf("%s/10-ST_analysis/00-ST_Pipeline/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s/10-Pipeline/%s_stdata.tsv", work_dir, genome_name, sample, sample)
    st_pipeline_matrix <- read.table(st_pipeline_matrix_file, sep="\t", header=TRUE, quote="", row.names=1)
    # Space Ranger output directory for tissue positions list and image
    space_ranger_output_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs", work_dir, genome_name, sample)
    spe <- load_ST_pipeline_data(st_pipeline_matrix_file, space_ranger_output_dir)
    
    ## QC
    spe_seurat_qc <- Seurat_pipeline(spe, lib_size_threshold, exp_genes_threshold, space_ranger_output_dir, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, output_dir, sample)
    qc_output_basename <- sprintf("QC-libsize%d_exp%d", lib_size_threshold, exp_genes_threshold)
    qc_output_dir <- sprintf("%s/00-%s", output_dir, qc_output_basename)
    
    ## normalization, clustering, marker genes
    for (normalization_method in normalization_method_vector) {
        print(sprintf("normalization method: %s", normalization_method))
        for (features_nb in features_nb_vector) {
            print(sprintf("features number: %d", features_nb))
            ## normalization, feature selection
            spe_seurat_features <- Seurat_pipeline(spe_seurat_qc, NA, NA, NA, normalization_method, feature_selection_method, features_nb, feature_residual_variance_threshold, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, qc_output_dir, sample)
            normalization_output <- get_normalization_subdir_basename(normalization_method)
            feature_selection_output <- get_features_subdir_basename(feature_selection_method, features_nb, feature_residual_variance_threshold)
            features_output_dir <- sprintf("%s/00-%s/%s/%s/%s", output_dir, qc_output_basename, normalization_output[["subdir"]], feature_selection_output[["method_subdir"]], feature_selection_output[["param_subdir"]])

            ## clustering, marker genes
            for (dimensions_nb_param in dimensions_nb_param_vector) {
                print(sprintf("dimension number param: %s", dimensions_nb_param))
                for (nearest_neighbors_nb in nearest_neighbors_nb_vector) {
                    print(sprintf("number of nearest neighbors: %d", nearest_neighbors_nb))
                    spe_seurat_knn <- Seurat_pipeline(spe_seurat_features, NA, NA, NA, NA, NA, NA, NA, dimensions_nb_param, nearest_neighbors_nb, NA, NA, NA, NA, NA, NA, NA, NA, features_output_dir, sample)
                    clustering_output <- get_clustering_subdir_basename(dimensions_nb_param, nearest_neighbors_nb, NA, NA, sample)
                    clustering_output_dir <- sprintf("%s/%s/%s", features_output_dir, clustering_output[["dimensions_subdir"]], clustering_output[["knn_subdir"]])
                    clustering_output_basename <- clustering_output[["basename"]]

                    for (clustering_resolution in clustering_resolution_vector) {
                        print(sprintf("cluster resolution: %.1f", clustering_resolution))
                        mean_silhouette_score_vector <- c()
                        for (clustering_method in clustering_method_vector) {
                            print(sprintf("clustering method number: %d", clustering_method))
                            ### clustering and marker genes can not be executed in a single Seurat_pipeline() call because of out_name paramter set to different values
                            ### clustering
                            bootstrap_nb_spots <- cluster_analysis_spot_pct*ncol(spe_seurat_knn)
                            spe_seurat_clustering <- Seurat_pipeline(spe_seurat_knn, NA, NA, NA, NA, NA, NA, NA, NA, NA, clustering_resolution, clustering_method, cluster_analysis, bootstrap_nb_spots, cluster_analysis_iterations, NA, NA, NA, clustering_output_dir, clustering_output_basename)
                            ### marker genes
                            clustering_output <- get_clustering_subdir_basename(NA, NA, clustering_resolution, clustering_method, clustering_output_basename)
                            clustering_method_output_dir <- sprintf("%s/%s/%s", clustering_output_dir, clustering_output[["resolution_subdir"]], clustering_output[["clustering_subdir"]])
                            spe_seurat <- Seurat_pipeline(spe_seurat_clustering, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, markers_test, min_logFC, min_pct_spots, clustering_method_output_dir, sample)
                            ### get mean silhouette scores
                            silhouette_file <- sprintf("%s/%s_analysis_silhouette.csv", clustering_method_output_dir, clustering_output[["basename"]])
                            silhouette_df <- read.csv(silhouette_file, quote="")
                            mean_silhouette <- mean(silhouette_df$silhouette_score)
                            mean_silhouette_score_vector <- c(mean_silhouette_score_vector, mean_silhouette)
                            ### get cluster labels on tissue section plot
                            spatial_dimplot <- SpatialDimPlot(spe_seurat, label=TRUE, label.size=3) + labs(title=sprintf("Normalization method: %s, number of features: %d\nNumber of dimensions: %s, number of nearest neighbors: %d\nClustering resolution: %.1f, clustering method number: %d\nMean silhouette score: %f", normalization_method, features_nb, dimensions_nb_param, nearest_neighbors_nb, clustering_resolution, clustering_method, mean_silhouette))
                            param_setting <- paste(c(normalization_method, features_nb, dimensions_nb_param, nearest_neighbors_nb, sub("[.]", "_", clustering_resolution), clustering_method), collapse="_")
                            sample_plot_list[[param_setting]] <- spatial_dimplot
                        }
                        df2ggplot <- data.frame(normalization=rep(normalization_method, length(clustering_method_vector)), features_nb=rep(features_nb, length(clustering_method_vector)), dimensions=rep(dimensions_nb_param, length(clustering_method_vector)), knn=rep(nearest_neighbors_nb, length(clustering_method_vector)), resolution=rep(clustering_resolution, length(clustering_method_vector)), clustering_method=clustering_method_vector, silhouette=mean_silhouette_score_vector)
                        sample_df2ggplot <- rbind(sample_df2ggplot, df2ggplot)
                    }
                }
            }
        }
    }
    # silhouette score plots for all parameter settings
    sample_df2ggplot$normalization <- as.factor(sample_df2ggplot$normalization)
    sample_df2ggplot$features_nb <- as.factor(sample_df2ggplot$features_nb)
    sample_df2ggplot$dimensions <- as.factor(sample_df2ggplot$dimensions)
    sample_df2ggplot$knn <- as.factor(sample_df2ggplot$knn)
    sample_df2ggplot$resolution <- as.factor(sample_df2ggplot$resolution)
    sample_df2ggplot$clustering_method <- as.factor(sample_df2ggplot$clustering_method)
    sample_df2ggplot$rank <- rank(-sample_df2ggplot$silhouette)
    sample_df2ggplot$param_setting <- apply(sample_df2ggplot, 1, function(x) { paste(c(x["normalization"], sub("[.]", "_", x["features_nb"]), x["dimensions"], x["knn"], ifelse(x["resolution"]==1, x["resolution"], sub("[.]", "_", x["resolution"])), x["clustering_method"]), collapse="_") })
    sample_df2ggplot$param_setting <- as.factor(sample_df2ggplot$param_setting)
    write.csv(sample_df2ggplot, file=sprintf("%s/%s-%s_mean_silhouette.csv", qc_output_dir, sample, qc_output_basename), quote=FALSE, row.names=FALSE)
    
    pdf(sprintf("%s/%s-%s_mean_silhouette.pdf", qc_output_dir, sample, qc_output_basename))
    
    sample_df2ggplot$param_setting <- as.character(sample_df2ggplot$param_setting)
    p <- ggplot(sample_df2ggplot, aes(x=param_setting, y=silhouette)) +
        geom_line(group=1) +
        geom_point() +
        labs(x="Parameter setting", y="Mean silhouette score") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=90, hjust=1, colour=ifelse(sample_df2ggplot$rank < 10, "red", "black"), size=6)) +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
    
    p <- ggplot(sample_df2ggplot, aes(x=normalization, y=silhouette, fill=features_nb)) +
        geom_boxplot() +
        facet_wrap(~dimensions) +
        labs(x="Normalization method", y="Mean silhouette score", fill="Number of\nfeatures") +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50"))
    print(p) 
    
    for (normalization_method in normalization_method_vector) {
        p <- ggplot(sample_df2ggplot[which(sample_df2ggplot$normalization==normalization_method),], aes(x=features_nb, y=silhouette, fill=dimensions)) +
            geom_boxplot() +
            facet_wrap(~knn) +
            labs(title=sprintf("Normalization: %s", normalization_method), x="Number of features", y="Mean silhouette score", fill="Number of\ndimensions") +
            theme_bw() +
            theme(panel.border=element_rect(color="grey50"))
        print(p)

        for (features_nb in features_nb_vector) {
            p <- ggplot(sample_df2ggplot[which(sample_df2ggplot$normalization=="SCTransform" & sample_df2ggplot$features_nb==features_nb),], aes(x=dimensions, y=silhouette, fill=knn)) +
                geom_boxplot() +
                facet_wrap(~resolution) +
                labs(title=sprintf("Normalization: %s\nnumber of features: %d", normalization_method, features_nb), x="Number of dimensions", y="Mean silhouette score", fill="Number of nearest\nneighbors") +
                theme_bw() +
                theme(panel.border=element_rect(color="grey50"))
            print(p)

            for (dimensions_nb_param in dimensions_nb_param_vector) {
                p <- ggplot(sample_df2ggplot[which(sample_df2ggplot$normalization=="SCTransform" & sample_df2ggplot$features_nb==features_nb & sample_df2ggplot$dimensions==dimensions_nb_param),], aes(x=knn, y=silhouette, fill=resolution)) +
                    geom_bar(stat="identity", position="dodge") +
                    facet_wrap(~clustering_method) +
                    labs(title=sprintf("Normalization: %s\nnumber of features: %d\nnumber of dimensions: %s", normalization_method, features_nb, dimensions_nb_param), x="Number of nearest\nneighbors", y="Mean silhouette score", fill="Resolution") +
                    theme_bw() +
                    theme(panel.border=element_rect(color="grey50"))
                print(p)
            }
        }
    }
    
    # number of times in the top 10 parameter settings
    top_df <- data.frame()
    for (top in c(3, 10, 20)) {
        df2ggplot <- top_counts(sample_df2ggplot, "normalization", top)
        normalization_plot <- top_counts_barplot(df2ggplot, "Normalization method", "Normalization\nmethod", top)
        top_df <- rbind(top_df, cbind(parameter=rep("normalization method", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        df2ggplot <- top_counts(sample_df2ggplot, "features_nb", top)
        features_nb_plot <- top_counts_barplot(df2ggplot, "Number of features", "Number of\nfeatures", top)
        top_df <- rbind(top_df, cbind(parameter=rep("number of features", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        df2ggplot <- top_counts(sample_df2ggplot, "dimensions", top)
        dimensions_plot <- top_counts_barplot(df2ggplot, "Number of dimensions", "Number of\ndimensions", top)
        top_df <- rbind(top_df, cbind(parameter=rep("number of dimensions", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        df2ggplot <- top_counts(sample_df2ggplot, "knn", top)
        knn_plot <- top_counts_barplot(df2ggplot, "Number of nearest neighbors", "Number of nearest\nneighbors", top)
        top_df <- rbind(top_df, cbind(parameter=rep("number of nearest neighbors", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        df2ggplot <- top_counts(sample_df2ggplot, "resolution", top)
        resolution_plot <- top_counts_barplot(df2ggplot, "Clustering resolution", "Clustering\nresolution", top)
        top_df <- rbind(top_df, cbind(parameter=rep("resolution", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        df2ggplot <- top_counts(sample_df2ggplot, "clustering_method", top)
        clustering_method_plot <- top_counts_barplot(df2ggplot, "Clustering method", "Clustering\nmethod", top)
        top_df <- rbind(top_df, cbind(parameter=rep("clustering method", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        multiplot <- ggdraw() +
          draw_plot(normalization_plot, 0, 0.5, 0.33, 0.5) +
          draw_plot(features_nb_plot, 0.33, 0.5, 0.33, 0.5) +
          draw_plot(dimensions_plot, 0.66, 0.5, 0.33, 0.5) +
          draw_plot(knn_plot, 0, 0, 0.33, 0.5) +
          draw_plot(resolution_plot, 0.33, 0, 0.33, 0.5) +
          draw_plot(clustering_method_plot, 0.66, 0, 0.33, 0.5)
        print(multiplot)
    }
    dev.off()
    write.csv(top_df, file=sprintf("%s/%s-%s_mean_silhouette_top_parameter_settings.csv", qc_output_dir, sample, qc_output_basename), quote=FALSE, row.names=FALSE)
    
    all_samples_df2ggplot <- rbind(all_samples_df2ggplot, cbind(sample=rep(sample, dim(sample_df2ggplot)[1]), sample_df2ggplot))
    all_samples_top_df <- rbind(all_samples_top_df, cbind(sample=rep(sample, dim(top_df)[1]), top_df))
    
    # cluster labels on tissue section plot in mean silhouette score decreasing order
    pdf(sprintf("%s/%s-%s_clusters_tissue_section.pdf", qc_output_dir, sample, qc_output_basename))
    for (one_param_setting in as.character(sample_df2ggplot[order(sample_df2ggplot$rank), "param_setting"])) {
        print(one_param_setting)
        print(sample_plot_list[[one_param_setting]])
    }
    dev.off()
}
all_samples_df2ggplot$normalization <- as.factor(all_samples_df2ggplot$normalization)
all_samples_df2ggplot$features_nb <- as.factor(all_samples_df2ggplot$features_nb)
all_samples_df2ggplot$dimensions <- as.factor(all_samples_df2ggplot$dimensions)
all_samples_df2ggplot$knn <- as.factor(all_samples_df2ggplot$knn)
all_samples_df2ggplot$resolution <- as.factor(all_samples_df2ggplot$resolution)
all_samples_df2ggplot$clustering_method <- as.factor(all_samples_df2ggplot$clustering_method)


output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/10-Dataset", work_dir, genome_name)
qc_output_dir <- sprintf("%s/00-%s", output_dir, qc_output_basename)
if (! dir.exists(qc_output_dir)) {
    dir.create(qc_output_dir, recursive=TRUE, mode="0775")
}
write.csv(all_samples_df2ggplot, file=sprintf("%s/%s_mean_silhouette.csv", qc_output_dir, qc_output_basename), quote=FALSE, row.names=FALSE)
pdf(sprintf("%s/%s_mean_silhouette.pdf", qc_output_dir, qc_output_basename))
p <- ggplot(all_samples_df2ggplot, aes(x=normalization, y=silhouette, fill=normalization)) +
    geom_boxplot() +
    labs(x="Normalization method", y="Mean silhouette score", fill="Normalization\nmethod") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50"))
print(p)

p <- ggplot(all_samples_df2ggplot, aes(x=features_nb, y=silhouette, fill=features_nb)) +
    geom_boxplot() +
    facet_wrap(~normalization) +
    labs(x="Number of features", y="Mean silhouette score", fill="Number of\nfeatures") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50"))
print(p)

p <- ggplot(all_samples_df2ggplot, aes(x=dimensions, y=silhouette, fill=dimensions)) +
    geom_boxplot() +
    facet_grid(features_nb~normalization) +
    labs(x="Number of dimensions", y="Mean silhouette score", fill="Number of\ndimensions") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50"))
print(p)
for (normalization_method in normalization_method_vector) {
    p <- ggplot(all_samples_df2ggplot[which(all_samples_df2ggplot$normalization==normalization_method),], aes(x=knn, y=silhouette, fill=knn)) +
        geom_boxplot() +
        facet_grid(dimensions~features_nb) +
        labs(title=sprintf("Normalization: %s", normalization_method), x="Number of nearest neighbors", y="Mean silhouette score", fill="Number of\nnearest neighbors") +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
    
    for (features_nb in features_nb_vector) {
        p <- ggplot(all_samples_df2ggplot[which(all_samples_df2ggplot$normalization=="SCTransform" & all_samples_df2ggplot$features_nb==features_nb),], aes(x=resolution, y=silhouette, fill=resolution)) +
            geom_boxplot() +
            facet_grid(knn~dimensions) +
            labs(title=sprintf("Normalization: %s\nnumber of features: %d", normalization_method, features_nb), x="Resolution", y="Mean silhouette score", fill="Resolution") +
            theme_bw() +
            theme(panel.border=element_rect(color="grey50"))
        print(p)
        
        for (dimensions_nb_param in dimensions_nb_param_vector) {
            p <- ggplot(all_samples_df2ggplot[which(all_samples_df2ggplot$normalization=="SCTransform" & all_samples_df2ggplot$features_nb==features_nb & all_samples_df2ggplot$dimensions==dimensions_nb_param),], aes(x=clustering_method, y=silhouette, fill=clustering_method)) +
                geom_bar(stat="identity") +
                facet_grid(resolution~knn) +
                labs(title=sprintf("Normalization: %s\nnumber of features: %d\nnumber of dimensions: %s", normalization_method, features_nb, dimensions_nb_param), x="Clustering method", y="Mean silhouette score", fill="Clustering\nmethod") +
                theme_bw() +
                theme(panel.border=element_rect(color="grey50"))
            print(p)
        }
    }
}
dev.off()









# multiple parameter evaluation: feature selection based on residual variance
## SCTransform
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"
genome_name <- "NCBIRefSeq108_NCBIRefSeq108GTF"
#sample_vector <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")
sample_vector <- c("D_L2_S8", "D_L3_S12")
# pipeline parameters
## QC
lib_size_threshold <- 2000
exp_genes_threshold <- 1000
## normalization
#normalization_method <- "SCTransform" # either "SCTransform" or normalization.method parameter of NormalizeData() function, one of the following values "LogNormalize", "CLR", "RC"
normalization_method_vector <- c("SCTransform") # TODO: can't use feature residual variance threshold with LogNormalize
## clustering
feature_selection_method <- "HVGs" # either "HVGs" or "SVGs"
features_nb <- NA
feature_residual_variance_threshold_vector <- c(1, 1.1, 1.3, 1.5, 2, 2.5)
#dimensions_nb_param <- "estimation" # either "estimation" or number of PCs
dimensions_nb_param_vector <- c(0, 10, 15, 20)
#nearest_neighbors_nb <- 20 # k.param parameter of FindNeighbors() function
nearest_neighbors_nb_vector <- c(20)
clustering_resolution_vector <- c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2)
#clustering_method <- "Louvain" # algorithm parameter of FindClusters() function, one of the following values: "Louvain", ...
clustering_method_vector <- c(1) # algorithm parameter of FindClusters() function, one of the following values: "Louvain", ...
cluster_analysis <- TRUE
cluster_analysis_spot_pct <- 0.8
cluster_analysis_iterations <- 20
## marker genes
markers_test <- "wilcox"  # test.use parameter of FindAllMarkers() function, one of the following values: "wilcox", ...
min_logFC <- 0.25 # logfc.threshold parameter of FindAllMarkers() function
min_pct_spots <- 0.1 # min.pct parameter of FindAllMarkers() function
## output basename
qc_output_basename <- sprintf("QC-libsize%d_exp%d", lib_size_threshold, exp_genes_threshold)
output_basename <- sprintf("%s-features_residual_variance", qc_output_basename)

all_samples_df2ggplot <- data.frame()
all_samples_feature_residual_variance_df <- data.frame()
all_samples_dimensions_estimation_df <- data.frame()
for (sample in sample_vector) {
    print(sprintf("sample: %s", sample))
    output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s", work_dir, genome_name, sample)
    if (! dir.exists(output_dir)) {
        dir.create(output_dir, recursive=TRUE, mode="0775")
    }
    sample_output_basename <- sprintf("%s-%s", sample, output_basename)
    sample_df2ggplot <- data.frame()
    sample_plot_list <- list()

    # Load data
    library(SpatialExperiment)
    # expression matrices
    ## Space Ranger
    #h5_file <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs/raw_feature_bc_matrix.h5", work_dir, genome_name, sample)
    #h5_data <- Read10X_h5(h5_file)
    ## ST Pipeline
    st_pipeline_matrix_file <- sprintf("%s/10-ST_analysis/00-ST_Pipeline/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s/10-Pipeline/%s_stdata.tsv", work_dir, genome_name, sample, sample)
    st_pipeline_matrix <- read.table(st_pipeline_matrix_file, sep="\t", header=TRUE, quote="", row.names=1)
    # Space Ranger output directory for tissue positions list and image
    space_ranger_output_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs", work_dir, genome_name, sample)
    spe <- load_ST_pipeline_data(st_pipeline_matrix_file, space_ranger_output_dir)
    
    ## hippocampus manual annotations
    manual_annotations_file <- sprintf("%s/manual_annotations.csv", space_ranger_output_dir)
    manual_annotations_df <- read.csv(file=manual_annotations_file)
    
    ## QC
    spe_seurat_qc <- Seurat_pipeline(spe, lib_size_threshold, exp_genes_threshold, space_ranger_output_dir, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, output_dir, sample)
    qc_output_dir <- sprintf("%s/00-%s", output_dir, qc_output_basename)
    
    ## normalization, clustering, marker genes
    for (normalization_method in normalization_method_vector) {
        print(sprintf("normalization method: %s", normalization_method))
        
        nb_features_vector <- c()
        for (feature_residual_variance_threshold in feature_residual_variance_threshold_vector) {
            print(sprintf("feature residual variance threshold: %.1f", feature_residual_variance_threshold))
            ## normalization, feature selection
            spe_seurat_features <- Seurat_pipeline(spe_seurat_qc, NA, NA, NA, normalization_method, feature_selection_method, features_nb, feature_residual_variance_threshold, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, qc_output_dir, sample)
            nb_features <- length(spe_seurat_features@assays$SCT@var.features)
            nb_features_vector <- c(nb_features_vector, nb_features)
            
            normalization_output <- get_normalization_subdir_basename(normalization_method)
            feature_selection_output <- get_features_subdir_basename(feature_selection_method, features_nb, feature_residual_variance_threshold)
            features_output_dir <- sprintf("%s/00-%s/%s/%s/%s", output_dir, qc_output_basename, normalization_output[["subdir"]], feature_selection_output[["method_subdir"]], feature_selection_output[["param_subdir"]])

            ## clustering, marker genes
            nb_dimensions_vector <- feature_resvar_vector <- c()
            for (dimensions_nb_param in dimensions_nb_param_vector) {
                if (dimensions_nb_param == 0) {
                    dimensions_nb_param <- "estimation"
                }
                print(sprintf("dimension number param: %s", dimensions_nb_param))
                for (nearest_neighbors_nb in nearest_neighbors_nb_vector) {
                    print(sprintf("number of nearest neighbors: %d", nearest_neighbors_nb))
                    spe_seurat_knn <- Seurat_pipeline(spe_seurat_features, NA, NA, NA, NA, NA, NA, NA, dimensions_nb_param, nearest_neighbors_nb, NA, NA, NA, NA, NA, NA, NA, NA, features_output_dir, sample)
                    clustering_output <- get_clustering_subdir_basename(dimensions_nb_param, nearest_neighbors_nb, NA, NA, sample)
                    clustering_output_dir <- sprintf("%s/%s/%s", features_output_dir, clustering_output[["dimensions_subdir"]], clustering_output[["knn_subdir"]])
                    clustering_output_basename <- clustering_output[["basename"]]

                    for (clustering_resolution in clustering_resolution_vector) {
                        print(sprintf("cluster resolution: %.1f", clustering_resolution))
                        mean_silhouette_score_vector <- hippocampus_mean_silhouette_score_vector <- c()
                        for (clustering_method in clustering_method_vector) {
                            print(sprintf("clustering method number: %d", clustering_method))
                            ### clustering and marker genes can not be executed in a single Seurat_pipeline() call because of out_name paramter set to different values
                            ### clustering
                            bootstrap_nb_spots <- cluster_analysis_spot_pct*ncol(spe_seurat_knn)
                            spe_seurat_clustering <- Seurat_pipeline(spe_seurat_knn, NA, NA, NA, NA, NA, NA, NA, NA, NA, clustering_resolution, clustering_method, cluster_analysis, bootstrap_nb_spots, cluster_analysis_iterations, NA, NA, NA, clustering_output_dir, clustering_output_basename)
                            ### marker genes
                            clustering_output <- get_clustering_subdir_basename(NA, NA, clustering_resolution, clustering_method, clustering_output_basename)
                            clustering_method_output_dir <- sprintf("%s/%s/%s", clustering_output_dir, clustering_output[["resolution_subdir"]], clustering_output[["clustering_subdir"]])
                            spe_seurat <- Seurat_pipeline(spe_seurat_clustering, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, markers_test, min_logFC, min_pct_spots, clustering_method_output_dir, sample)
                            ### get mean silhouette scores
                            silhouette_file <- sprintf("%s/%s_analysis_silhouette.csv", clustering_method_output_dir, clustering_output[["basename"]])
                            silhouette_df <- read.csv(silhouette_file, quote="")
                            mean_silhouette <- mean(silhouette_df$silhouette_score)
                            mean_silhouette_score_vector <- c(mean_silhouette_score_vector, mean_silhouette)
                            silhouette_annotations_df <- merge(silhouette_df, manual_annotations_df, by.x="barcode", by.y="Barcode")
                            hippocampus_mean_silhouette <- mean(silhouette_annotations_df[which(silhouette_annotations_df$manual.annotations=="hippocampus"), "silhouette_score"])
                            hippocampus_mean_silhouette_score_vector <- c(hippocampus_mean_silhouette_score_vector, hippocampus_mean_silhouette)
                            
                            ### get cluster labels on tissue section plot
                            spatial_dimplot <- SpatialDimPlot(spe_seurat, label=TRUE, label.size=3) + labs(title=sprintf("Normalization method: %s, feature residual variance threshold: %.1f\nNumber of dimensions: %s, number of nearest neighbors: %d\nClustering resolution: %.1f, clustering method number: %d\nMean silhouette score - all spots: %f, hippocampus spots: %f", normalization_method, feature_residual_variance_threshold, dimensions_nb_param, nearest_neighbors_nb, clustering_resolution, clustering_method, mean_silhouette, hippocampus_mean_silhouette))
                            param_setting <- paste(c(normalization_method, sub("[.]", "_", feature_residual_variance_threshold), dimensions_nb_param, nearest_neighbors_nb, sub("[.]", "_", clustering_resolution), clustering_method), collapse="_")
                            sample_plot_list[[param_setting]] <- spatial_dimplot
                            
                        }
                        df2ggplot <- data.frame(normalization=rep(normalization_method, length(clustering_method_vector)), features_resvar=rep(feature_residual_variance_threshold, length(clustering_method_vector)), dimensions=rep(dimensions_nb_param, length(clustering_method_vector)), knn=rep(nearest_neighbors_nb, length(clustering_method_vector)), resolution=rep(clustering_resolution, length(clustering_method_vector)), clustering_method=clustering_method_vector, silhouette=mean_silhouette_score_vector, silhouette_hippocampus=hippocampus_mean_silhouette_score_vector)
                        sample_df2ggplot <- rbind(sample_df2ggplot, df2ggplot)
                    }
                }
                if (dimensions_nb_param == "estimation") {
                    nb_dimensions <- spe_seurat_knn@commands$FindNeighbors.SCT.pca@params$dims[length(spe_seurat_knn@commands$FindNeighbors.SCT.pca@params$dims)]
                    nb_dimensions_vector <- c(nb_dimensions_vector, nb_dimensions)
                    feature_resvar_vector <- c(feature_resvar_vector, feature_residual_variance_threshold)
                }
            }
        }
        sample_feature_residual_variance_df <- data.frame(sample=rep(sample, length(feature_residual_variance_threshold_vector)), threshold=feature_residual_variance_threshold_vector, nb_features=nb_features_vector)
        sample_dimensions_estimation_df <- data.frame(sample=rep(sample, length(nb_dimensions_vector)), feature_resvar=feature_resvar_vector, nb_dimensions=nb_dimensions_vector)
    }
    # silhouette score plots for all parameter settings
    sample_df2ggplot$normalization <- as.factor(sample_df2ggplot$normalization)
    sample_df2ggplot$features_resvar <- as.factor(sample_df2ggplot$features_resvar)
    sample_df2ggplot$dimensions <- as.factor(sample_df2ggplot$dimensions)
    sample_df2ggplot$knn <- as.factor(sample_df2ggplot$knn)
    sample_df2ggplot$resolution <- as.factor(sample_df2ggplot$resolution)
    sample_df2ggplot$clustering_method <- as.factor(sample_df2ggplot$clustering_method)
    sample_df2ggplot$rank <- rank(-sample_df2ggplot$silhouette)
    sample_df2ggplot$rank_hippocampus <- rank(-sample_df2ggplot$silhouette_hippocampus)
    sample_df2ggplot$param_setting <- apply(sample_df2ggplot, 1, function(x) { paste(c(x["normalization"], sub("[.]", "_", x["features_resvar"]), x["dimensions"], x["knn"], sub("[.]", "_", x["resolution"]), x["clustering_method"]), collapse="_") })
    sample_df2ggplot$param_setting <- as.factor(sample_df2ggplot$param_setting)
    write.csv(sample_df2ggplot, file=sprintf("%s/%s_mean_silhouette.csv", qc_output_dir, sample_output_basename), quote=FALSE, row.names=FALSE)
    
    #sample_df2ggplot$param_setting <- as.character(sample_df2ggplot$param_setting)
    mean_silhouette_plots(sample_df2ggplot, "silhouette", "rank", qc_output_dir, sprintf("%s_mean_silhouette", sample_output_basename))
    mean_silhouette_plots_top(sample_df2ggplot, "rank", qc_output_dir, sprintf("%s_mean_silhouette_top_parameter_settings", sample_output_basename))  
    mean_silhouette_plots(sample_df2ggplot, "silhouette_hippocampus", "rank_hippocampus", qc_output_dir, sprintf("%s_mean_silhouette-hippocampus", sample_output_basename))
    mean_silhouette_plots_top(sample_df2ggplot, "rank_hippocampus", qc_output_dir, sprintf("%s_mean_silhouette-hippocampus_top_parameter_settings", sample_output_basename))
    
    # cluster labels on tissue section plot in mean silhouette score decreasing order
    pdf(sprintf("%s/%s_clusters_tissue_section.pdf", qc_output_dir, sample_output_basename))
    for (one_param_setting in as.character(sample_df2ggplot[order(sample_df2ggplot$rank_hippocampus), "param_setting"])) {
        print(one_param_setting)
        all_spots_rank <- sample_df2ggplot[which(sample_df2ggplot$param_setting==one_param_setting), "rank"]
        hippocampus_spots_rank <- sample_df2ggplot[which(sample_df2ggplot$param_setting==one_param_setting), "rank_hippocampus"]
        one_tissue_section_plot <- sample_plot_list[[one_param_setting]]
        one_tissue_section_plot <- one_tissue_section_plot + labs(title=sprintf("%s\nMean silhouette score rank- all spots: %.1f, hippocampus spots: %.1f", one_tissue_section_plot$labels$title, all_spots_rank, hippocampus_spots_rank))
        print(one_tissue_section_plot)
    }
    dev.off()
    
    all_samples_df2ggplot <- rbind(all_samples_df2ggplot, cbind(sample=rep(sample, dim(sample_df2ggplot)[1]), sample_df2ggplot))
    all_samples_feature_residual_variance_df <- rbind(all_samples_feature_residual_variance_df, sample_feature_residual_variance_df)
    all_samples_dimensions_estimation_df <- rbind(all_samples_dimensions_estimation_df, sample_dimensions_estimation_df)
}
all_samples_df2ggplot$normalization <- as.factor(all_samples_df2ggplot$normalization)
all_samples_df2ggplot$features_nb <- as.factor(all_samples_df2ggplot$features_nb)
all_samples_df2ggplot$dimensions <- as.factor(all_samples_df2ggplot$dimensions)
all_samples_df2ggplot$knn <- as.factor(all_samples_df2ggplot$knn)
all_samples_df2ggplot$resolution <- as.factor(all_samples_df2ggplot$resolution)
all_samples_df2ggplot$clustering_method <- as.factor(all_samples_df2ggplot$clustering_method)


output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/10-Dataset", work_dir, genome_name)
qc_output_dir <- sprintf("%s/00-%s", output_dir, qc_output_basename)
features_output_dir <- sprintf("%s/%s/%s", qc_output_dir, normalization_output[["subdir"]], feature_selection_output[["method_subdir"]])
if (! dir.exists(features_output_dir)) {
    dir.create(features_output_dir, recursive=TRUE, mode="0775")
}
features_output_basename <- sprintf("%s-%s", normalization_output[["basename"]], feature_selection_output[["method_basename"]])
# number of features according to residual variance threshold
all_samples_feature_residual_variance_df$sample <- as.factor(all_samples_feature_residual_variance_df$sample)
all_samples_feature_residual_variance_df$threshold <- as.factor(all_samples_feature_residual_variance_df$threshold)
write.csv(all_samples_feature_residual_variance_df, file=sprintf("%s/%s_feature_residual_variance_features_number.csv", features_output_dir, features_output_basename), quote=FALSE, row.names=FALSE)
pdf(sprintf("%s/%s_feature_residual_variance_features_number.pdf", features_output_dir, features_output_basename))
p <- ggplot(all_samples_feature_residual_variance_df, aes(x=sample, y=nb_features, fill=sample)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(~threshold) +
    labs(title="Number of features according to residual variance threshold", x="Sample", y="Number of features", fill="Sample") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=60, hjust=1, size=5)) +
    theme(panel.border=element_rect(color="grey50"))
print(p)

p <- ggplot(all_samples_feature_residual_variance_df, aes(x=sample, y=nb_features, fill=threshold)) +
    geom_bar(stat="identity", position="dodge") +
    labs(title="Number of features according to residual variance threshold", x="Sample", y="Number of features", fill="Threshold") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    theme(panel.border=element_rect(color="grey50"))
print(p)

p <- ggplot(all_samples_feature_residual_variance_df, aes(x=threshold, y=nb_features, fill=threshold)) +
    geom_boxplot() +
    labs(title="Number of features according to residual variance threshold", x="Threshold", y="Number of features", fill="Threshold") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50"))
print(p)

p <- ggplot(all_samples_feature_residual_variance_df, aes(x=threshold, y=nb_features, fill=sample)) +
    geom_bar(stat="identity", position="dodge") +
    labs(title="Number of features according to residual variance threshold", x="Threshold", y="Number of features", fill="Sample") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50"))
print(p)
dev.off()

# number of dimensions according to residual variance threshold
all_samples_dimensions_estimation_df$sample <- as.factor(all_samples_dimensions_estimation_df$sample)
write.csv(all_samples_feature_residual_variance_df, file=sprintf("%s/%s_feature_residual_variance_dimensions_number.csv", features_output_dir, features_output_basename), quote=FALSE, row.names=FALSE)
pdf(sprintf("%s/%s_feature_residual_variance_dimensions_number.pdf", features_output_dir, features_output_basename))
p <- ggplot(all_samples_dimensions_estimation_df, aes(x=feature_resvar, y=nb_dimensions, fill=feature_resvar)) +
    geom_boxplot() +
    facet_wrap(~feature_resvar) +
    labs(title="Estimated number of dimensions according to residual variance threshold", x="Feature residual variance threshold", y="Number of dimensions", fill="Threshold") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50"))
print(p)
p <- ggplot(all_samples_dimensions_estimation_df, aes(x=sample, y=nb_dimensions, fill=sample)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(~feature_resvar) +
    labs(title="Estimated number of dimensions according to residual variance threshold", x="Feature residual variance threshold", y="Number of dimensions", fill="Threshold") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50"))
print(p)
dev.off()









# Machado et al rat hippocampus marker genes
machado_marker_output_dir <- sprintf("%s/10-Machado_et_al", clustering_method_output_dir)
if (! dir.exists(machado_marker_output_dir)) {
    dir.create(machado_marker_output_dir, recursive=TRUE, mode="0775")
}
markers_per_page <- 6
pdf(sprintf("%s/%s_CA1_Machado_et_al_marker_genes.pdf", machado_marker_output_dir, clustering_output[["basename"]]))
CA1_marker_genes <- c("Cbln1", "Spata18", "Htr5b", "Epgn", "Wnt3a", "Ldb2", "Ush2a", "Ascl1", "C1qtnf7", "Wfs1", "Drd5", "Gabrr2", "Diaph3", "Ndst4", "Satb2", "Pex5l", "Htr1b", "Cdh6", "RGD1563354", "Cd3e", "Galr1", "Arhgap12", "Olfm3", "Lyzl4", "Ackr2", "Mybpc2", "Cdk15", "Crabp1", "AY172581.19", "AABR07061707.1", "Met", "AABR07044001.4", "AABR07070872.1", "Aox3")
markers2print <- CA1_marker_genes[CA1_marker_genes %in% rownames(spe_seurat_clustering@assays$SCT@counts)]
nb_markers2print <- length(markers2print)
for (i in seq(1, nb_markers2print, markers_per_page)) {
    if ((i+markers_per_page-1) <=  nb_markers2print) {
        print(SpatialFeaturePlot(spe_seurat, features=markers2print[i:(i+markers_per_page-1)], alpha=c(0.1,1)) + plot_annotation(title="CA1 marker genes", theme=theme(plot.title=element_text(size=14))))
    } else {
        print(SpatialFeaturePlot(spe_seurat, features=markers2print[i:nb_markers2print], alpha=c(0.1,1)) + plot_annotation(title="CA1 marker genes", theme=theme(plot.title=element_text(size=14))))
    }
}
dev.off()

pdf(sprintf("%s/%s_CA2_Machado_et_al_marker_genes.pdf", machado_marker_output_dir, clustering_output[["basename"]]))
CA2_marker_genes <- c("Casr", "Rtl3", "Cacng5", "Gabra6", "Amigo2", "Papln", "Sfrp2", "Glra1", "Glipr2", "Ccdc3", "Cyp1b1", "Fam72a", "Spp1", "AABR07001734.1")
markers2print <- CA2_marker_genes[CA2_marker_genes %in% rownames(spe_seurat_clustering@assays$SCT@counts)]
nb_markers2print <- length(markers2print)
for (i in seq(1, nb_markers2print, markers_per_page)) {
    if ((i+markers_per_page-1) <=  nb_markers2print) {
        print(SpatialFeaturePlot(spe_seurat, features=markers2print[i:(i+markers_per_page-1)], alpha=c(0.1,1)) + plot_annotation(title="CA2 marker genes", theme=theme(plot.title=element_text(size=14))))
    } else {
        print(SpatialFeaturePlot(spe_seurat, features=markers2print[i:nb_markers2print], alpha=c(0.1,1)) + plot_annotation(title="CA2 marker genes", theme=theme(plot.title=element_text(size=14))))
    }
}
dev.off()

pdf(sprintf("%s/%s_CA3_Machado_et_al_marker_genes.pdf", machado_marker_output_dir, clustering_output[["basename"]]))
CA3_marker_genes <- c("Mx1", "Arhgap36", "Gjd2", "Bhlhe23", "Col8a2", "Hcrtr2", "Cyp26b1", "Fxyd7", "Fermt1", "RT1-M6-1", "Gpr88", "Vgll3", "Frem3", "Olr1353", "Tacr2")
markers2print <- CA3_marker_genes[CA3_marker_genes %in% rownames(spe_seurat_clustering@assays$SCT@counts)]
nb_markers2print <- length(markers2print)
for (i in seq(1, nb_markers2print, markers_per_page)) {
    if ((i+markers_per_page-1) <=  nb_markers2print) {
        print(SpatialFeaturePlot(spe_seurat, features=markers2print[i:(i+markers_per_page-1)], alpha=c(0.1,1)) + plot_annotation(title="CA3 marker genes", theme=theme(plot.title=element_text(size=14))))
    } else {
        print(SpatialFeaturePlot(spe_seurat, features=markers2print[i:nb_markers2print], alpha=c(0.1,1)) + plot_annotation(title="CA3 marker genes", theme=theme(plot.title=element_text(size=14))))
    }
}
dev.off()

pdf(sprintf("%s/%s_DG_Machado_et_al_marker_genes.pdf", machado_marker_output_dir, clustering_output[["basename"]]))
DG_marker_genes <- c("RGD1304622", "Ros1", "Pim1", "Slc16a10", "Brs3", "Rsph10b", "Gabrr3", "Cldn1", "Il1rap", "Tgfb2", "Galnt10", "Brinp3", "Ppl", "Btg2", "Ccdc85a", "Ptchd1", "Prox1", "Rgs9", "Tmigd1", "Rgs13", "Stxbp6", "Kcnj2", "Lratd2", "Vangl2", "Trhr", "Tmco5a", "Nhlh1", "Kitlg", "Itgb4", "Wnt16", "Lrrc3b", "Fsip1", "Trpc6", "Fam163b", "Itga11", "Spo11", "Stk32c", "Cntn5", "Ccn4", "C1qtnf6", "Zfp467", "Mgst1", "Frmd4b", "Stk26", "Klhl38", "Ccdc33", "Zbtb8a", "Pdzd3", "Tp53i11", "Dach1", "Penk", "Amotl1", "Grik3", "Fut4", "Ddit4l2", "Dgkh", "Hdc", "Styk1", "Hspb2", "Gnal", "Plau", "Otog", "Zdhhc22", "Fam167a", "Filip1", "Adamts18", "Fst", "Unc5d", "Kiss1r", "Sema5a", "Perp", "RGD1562914", "Sox15", "Trdn", "Msr1", "Tmem164", "Zcchc12", "Grm2", "Rps6ka2", "Akap7", "Arg1", "Dock11", "Pde7b", "Polm", "Aebp1", "Adra2b", "Dsp", "Smoc2", "Syt10", "Cdh13", "Rfx3", "Zfp423", "Rbm20", "Glis3", "Synpo2", "Slc35f4", "Pgm5", "Ankdd1a", "Ptgfrn", "Rreb1", "Tshz1", "Gabrd", "Plekhg4", "Vwa5b1", "Eya4", "Tiam2", "Tcerg1l", "Fgfr4", "Cartpt", "Mpp7", "Armc4", "Mkx", "AABR07008724.1", "Nrg2", "Npr3", "Ablim3", "Adamts19", "Colq", "Ntf3", "Kcna5", "Actn3", "Sipa1l2", "Ryr1", "Efna1", "Slc4a11", "Tiam1", "Lrrtm4", "Fam216b", "Itga2b", "LOC103693984", "Olr1462", "Ctxn3", "Col11a1", "Adamts9", "Sned1", "Drd1", "Igsf3", "Fndc3b", "Ghsr", "Plet1", "Fhdc1", "Ccdc88b", "Atp8b1", "Rgs22", "Pdyn", "Tnfaip8", "Mgam", "Disp3", "Cdc42bpg", "Plppr3", "Fat4", "Slc9a5", "Trabd2b", "Pitpnm2", "LOC500877", "Gypc", "Usp13", "Stk32b", "Cdh23", "Svep1", "Plk5", "Btbd16", "Faap100", "Pde1b", "Tyms", "Calml4", "Plekha2", "Piezo2", "Ahcyl2", "Niban3", "Lnp1", "Slc44a5", "4933403O08Rik", "Rgs10", "Zmat4", "Rfx2", "Gpc6", "Msx3", "C1ql2", "Kcnk13", "E2f2", "Klhdc8b", "Fam160a1", "AABR07017624.2", "Traf1", "Krt71", "Cebpd", "AABR07013918.1", "AABR07070161.1", "Marveld2", "Pkd1l3", "Npnt", "AABR07068161.1", "Ror2", "Cracr2a", "LOC102548534", "Prelid2", "AABR07018116.1", "Srpk3", "AABR07041096.1", "AABR07058124.4", "AABR07012826.1", "Gcnt4", "RGD1564664", "SNORA63", "Unc13c", "Susd1", "AABR07013140.1", "7SK", "U6", "Tnc", "Oas3", "Adcy1", "AABR07060341.1", "AABR07070161.4", "Gpc3", "AABR07060341.2", "AABR07007717.3", "Pxdn", "Kcnc4", "Htra4", "AABR07060341.3", "AABR07061230.1", "U4", "Kcna3", "AABR07061237.1")
markers2print <- DG_marker_genes[DG_marker_genes %in% rownames(spe_seurat_clustering@assays$SCT@counts)]
nb_markers2print <- length(markers2print)
for (i in seq(1, nb_markers2print, markers_per_page)) {
    if ((i+markers_per_page-1) <=  nb_markers2print) {
        print(SpatialFeaturePlot(spe_seurat, features=markers2print[i:(i+markers_per_page-1)], alpha=c(0.1,1)) + plot_annotation(title="DG marker genes", theme=theme(plot.title=element_text(size=14))))
    } else {
        print(SpatialFeaturePlot(spe_seurat, features=markers2print[i:nb_markers2print], alpha=c(0.1,1)) + plot_annotation(title="DG marker genes", theme=theme(plot.title=element_text(size=14))))
    }
}
dev.off()




# manual annotations
## test STutility
library(STutility)
library(hdf5r)

work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"
genome_name <- "NCBIRefSeq108_NCBIRefSeq108GTF"
sample <- "A_L1_S1"
space_ranger_output_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs", work_dir, genome_name, sample)

h5_file <- sprintf("%s/raw_feature_bc_matrix.h5", space_ranger_output_dir)
tissue_positions_list_file <- sprintf("%s/spatial/tissue_positions_list.csv", space_ranger_output_dir)
he_image_file <- sprintf("%s/spatial/tissue_hires_image.png", space_ranger_output_dir)
scaling_factors_file <- sprintf("%s/spatial/scalefactors_json.json", space_ranger_output_dir)
df2input <- data.frame(samples=h5_file, spotfiles=tissue_positions_list_file, imgs=he_image_file, json=scaling_factors_file)
se <- InputFromTable(infotable=df2input)
se <- LoadImages(se, time.resolve = FALSE, verbose = TRUE)
options(browser="firefox")
options(shiny.port=9870)
se <- ManualAnnotation(se, type=NULL, res=1000, verbose=FALSE)


## Loupe Browser
sample <- "A_L1_S1"
space_ranger_output_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs", work_dir, genome_name, sample)
manual_annotations_file <- sprintf("%s/manual_annotations.csv", space_ranger_output_dir)
manual_annotations_df <- read.csv(file=manual_annotations_file)
silhouette_file <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/NCBIRefSeq108_NCBIRefSeq108GTF/00-Visium_recommended/00-Samples/A_L1_S1/00-QC-libsize2000_exp1000/00-Normalization-SCTransform/00-Features-HVGs/residualvariance1/10-10_PCs/k20/resolution0_4/00-Louvain/A_L1_S1-Clustering-10PCs_k20_resolution0_4_Louvain_analysis_silhouette.csv"
silhouette_df <- read.csv(silhouette_file, quote="")
silhouette_annotations_df <- merge(silhouette_df, manual_annotations_df, by.x="barcode", by.y="Barcode")
mean(silhouette_df$silhouette_score)
mean(silhouette_annotations_df[which(silhouette_annotations_df$manual.annotations=="hippocampus"), "silhouette_score"])

qc_output_basename <- sprintf("QC-libsize%d_exp%d", lib_size_threshold, exp_genes_threshold)
sample_vector <- c("A_L1_S1")
all_samples_df2ggplot <- data.frame()
all_samples_top_df <- data.frame()
for (sample in sample_vector) {
    sample_df2ggplot <- data.frame()
    sample_plot_list <- list()
    
    print(sprintf("sample: %s", sample))
    output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s", work_dir, genome_name, sample)
    qc_output_dir <- sprintf("%s/00-%s", output_dir, qc_output_basename)
    space_ranger_output_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs", work_dir, genome_name, sample)
    manual_annotations_file <- sprintf("%s/manual_annotations.csv", space_ranger_output_dir)
    manual_annotations_df <- read.csv(file=manual_annotations_file)
    
    ## normalization, clustering, marker genes
    for (normalization_method in normalization_method_vector) {
        print(sprintf("normalization method: %s", normalization_method))
        for (feature_residual_variance_threshold in feature_residual_variance_threshold_vector) {
            print(sprintf("feature residual variance threshold: %.1f", feature_residual_variance_threshold))
            normalization_output <- get_normalization_subdir_basename(normalization_method)
            feature_selection_output <- get_features_subdir_basename(feature_selection_method, features_nb, feature_residual_variance_threshold)
            features_output_dir <- sprintf("%s/00-%s/%s/%s/%s", output_dir, qc_output_basename, normalization_output[["subdir"]], feature_selection_output[["method_subdir"]], feature_selection_output[["param_subdir"]])

            ## clustering, marker genes
            for (dimensions_nb_param in dimensions_nb_param_vector) {
                if (dimensions_nb_param == 0) {
                    dimensions_nb_param <- "estimation"
                }
                print(sprintf("dimension number param: %s", dimensions_nb_param))
                for (nearest_neighbors_nb in nearest_neighbors_nb_vector) {
                    print(sprintf("number of nearest neighbors: %d", nearest_neighbors_nb))
                    clustering_output <- get_clustering_subdir_basename(dimensions_nb_param, nearest_neighbors_nb, NA, NA, sample)
                    clustering_output_dir <- sprintf("%s/%s/%s", features_output_dir, clustering_output[["dimensions_subdir"]], clustering_output[["knn_subdir"]])
                    clustering_output_basename <- clustering_output[["basename"]]

                    for (clustering_resolution in clustering_resolution_vector) {
                        print(sprintf("cluster resolution: %.1f", clustering_resolution))

                        hippocampus_mean_silhouette_score_vector <- c()
                        for (clustering_method in clustering_method_vector) {
                            print(sprintf("clustering method number: %d", clustering_method))
                            ### marker genes
                            clustering_output <- get_clustering_subdir_basename(NA, NA, clustering_resolution, clustering_method, clustering_output_basename)
                            clustering_method_output_dir <- sprintf("%s/%s/%s", clustering_output_dir, clustering_output[["resolution_subdir"]], clustering_output[["clustering_subdir"]])
                            ### get mean silhouette scores only for hippocampus spots
                            silhouette_file <- sprintf("%s/%s_analysis_silhouette.csv", clustering_method_output_dir, clustering_output[["basename"]])
                            silhouette_df <- read.csv(silhouette_file, quote="")
                            silhouette_annotations_df <- merge(silhouette_df, manual_annotations_df, by.x="barcode", by.y="Barcode")
                            hippocampus_mean_silhouette <- mean(silhouette_annotations_df[which(silhouette_annotations_df$manual.annotations=="hippocampus"), "silhouette_score"])
                            hippocampus_mean_silhouette_score_vector <- c(hippocampus_mean_silhouette_score_vector, hippocampus_mean_silhouette)
                        }
                        df2ggplot <- data.frame(normalization=rep(normalization_method, length(clustering_method_vector)), features_resvar=rep(feature_residual_variance_threshold, length(clustering_method_vector)), dimensions=rep(dimensions_nb_param, length(clustering_method_vector)), knn=rep(nearest_neighbors_nb, length(clustering_method_vector)), resolution=rep(clustering_resolution, length(clustering_method_vector)), clustering_method=clustering_method_vector, hippocampus_silhouette=hippocampus_mean_silhouette_score_vector)
                        sample_df2ggplot <- rbind(sample_df2ggplot, df2ggplot)
                    }
                }

            }
        }
    }
    sample_output_basename <- sprintf("%s-%s-features_residual_variance-hippocampus", sample, qc_output_basename)
    
    # silhouette score plots for all parameter settings
    sample_df2ggplot$normalization <- as.factor(sample_df2ggplot$normalization)
    sample_df2ggplot$features_resvar <- as.factor(sample_df2ggplot$features_resvar)
    sample_df2ggplot$dimensions <- as.factor(sample_df2ggplot$dimensions)
    sample_df2ggplot$knn <- as.factor(sample_df2ggplot$knn)
    sample_df2ggplot$resolution <- as.factor(sample_df2ggplot$resolution)
    sample_df2ggplot$clustering_method <- as.factor(sample_df2ggplot$clustering_method)
    sample_df2ggplot$rank <- rank(-sample_df2ggplot$hippocampus_silhouette)
    sample_df2ggplot$param_setting <- apply(sample_df2ggplot, 1, function(x) { paste(c(x["normalization"], sub("[.]", "_", x["features_resvar"]), x["dimensions"], x["knn"], sub("[.]", "_", x["resolution"]), x["clustering_method"]), collapse="_") })
    sample_df2ggplot$param_setting <- as.factor(sample_df2ggplot$param_setting)
    write.csv(sample_df2ggplot, file=sprintf("%s/%s_mean_silhouette.csv", qc_output_dir, sample_output_basename), quote=FALSE, row.names=FALSE)
    
    pdf(sprintf("%s/%s_mean_silhouette.pdf", qc_output_dir, sample_output_basename))
    
    sample_df2ggplot$param_setting <- as.character(sample_df2ggplot$param_setting)
    p <- ggplot(sample_df2ggplot, aes(x=param_setting, y=hippocampus_silhouette)) +
        geom_line(group=1) +
        geom_point() +
        labs(x="Parameter setting", y="Mean silhouette score") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=90, hjust=1, colour=ifelse(sample_df2ggplot$rank < 10, "red", "black"), size=6)) +
        theme(panel.border=element_rect(color="grey50"))
    print(p)
    
    p <- ggplot(sample_df2ggplot, aes(x=normalization, y=hippocampus_silhouette, fill=features_resvar)) +
        geom_boxplot() +
        facet_wrap(~dimensions) +
        labs(x="Normalization method", y="Mean silhouette score", fill="Residual variance threshold\nfor feature selection") +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50"))
    print(p) 
    
    for (normalization_method in normalization_method_vector) {
        p <- ggplot(sample_df2ggplot[which(sample_df2ggplot$normalization==normalization_method),], aes(x=features_resvar, y=hippocampus_silhouette, fill=dimensions)) +
            geom_boxplot() +
            facet_wrap(~knn) +
            labs(title=sprintf("Normalization: %s", normalization_method), x="Residual variance threshold for feature selection", y="Mean silhouette score", fill="Number of\ndimensions") +
            theme_bw() +
            theme(panel.border=element_rect(color="grey50"))
        print(p)

        for (feature_residual_variance_threshold in feature_residual_variance_threshold_vector) {
            p <- ggplot(sample_df2ggplot[which(sample_df2ggplot$normalization=="SCTransform" & sample_df2ggplot$features_resvar==feature_residual_variance_threshold),], aes(x=dimensions, y=hippocampus_silhouette, fill=knn)) +
                geom_boxplot() +
                facet_wrap(~resolution) +
                labs(title=sprintf("Normalization: %s\nresidual variance threshold for feature selection: %.1f", normalization_method, feature_residual_variance_threshold), x="Number of dimensions", y="Mean silhouette score", fill="Number of nearest\nneighbors") +
                theme_bw() +
                theme(panel.border=element_rect(color="grey50"))
            print(p)

            for (dimensions_nb_param in dimensions_nb_param_vector) {
                if (dimensions_nb_param == 0) {
                    dimensions_nb_param <- "estimation"
                }
                p <- ggplot(sample_df2ggplot[which(sample_df2ggplot$normalization=="SCTransform" & sample_df2ggplot$features_resvar==feature_residual_variance_threshold & sample_df2ggplot$dimensions==dimensions_nb_param),], aes(x=knn, y=hippocampus_silhouette, fill=resolution)) +
                    geom_bar(stat="identity", position="dodge") +
                    facet_wrap(~clustering_method) +
                    labs(title=sprintf("Normalization: %s\nresidual variance threshold for feature selection: %.1f\nnumber of dimensions: %s", normalization_method, feature_residual_variance_threshold, dimensions_nb_param), x="Number of nearest\nneighbors", y="Mean silhouette score", fill="Resolution") +
                    theme_bw() +
                    theme(panel.border=element_rect(color="grey50"))
                print(p)
            }
        }
    }
    
    # number of times in the top 10 parameter settings
    top_df <- data.frame()
    for (top in c(3, 10, 20)) {
        df2ggplot <- top_counts(sample_df2ggplot, "normalization", top)
        normalization_plot <- top_counts_barplot(df2ggplot, "Normalization method", "Normalization\nmethod", top)
        top_df <- rbind(top_df, cbind(parameter=rep("normalization method", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        df2ggplot <- top_counts(sample_df2ggplot, "features_resvar", top)
        features_resvar_plot <- top_counts_barplot(df2ggplot, "Feature residual variance threshold", "Feature residual\nvariance threshold", top)
        top_df <- rbind(top_df, cbind(parameter=rep("feature residual variance threshold", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        df2ggplot <- top_counts(sample_df2ggplot, "dimensions", top)
        dimensions_plot <- top_counts_barplot(df2ggplot, "Number of dimensions", "Number of\ndimensions", top)
        top_df <- rbind(top_df, cbind(parameter=rep("number of dimensions", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        df2ggplot <- top_counts(sample_df2ggplot, "knn", top)
        knn_plot <- top_counts_barplot(df2ggplot, "Number of nearest neighbors", "Number of nearest\nneighbors", top)
        top_df <- rbind(top_df, cbind(parameter=rep("number of nearest neighbors", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        df2ggplot <- top_counts(sample_df2ggplot, "resolution", top)
        resolution_plot <- top_counts_barplot(df2ggplot, "Clustering resolution", "Clustering\nresolution", top)
        top_df <- rbind(top_df, cbind(parameter=rep("resolution", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        df2ggplot <- top_counts(sample_df2ggplot, "clustering_method", top)
        clustering_method_plot <- top_counts_barplot(df2ggplot, "Clustering method", "Clustering\nmethod", top)
        top_df <- rbind(top_df, cbind(parameter=rep("clustering method", dim(df2ggplot)[1]), df2ggplot, top=rep(top, dim(df2ggplot)[1])))
        multiplot <- ggdraw() +
          draw_plot(normalization_plot, 0, 0.5, 0.33, 0.5) +
          draw_plot(features_resvar_plot, 0.33, 0.5, 0.33, 0.5) +
          draw_plot(dimensions_plot, 0.66, 0.5, 0.33, 0.5) +
          draw_plot(knn_plot, 0, 0, 0.33, 0.5) +
          draw_plot(resolution_plot, 0.33, 0, 0.33, 0.5) +
          draw_plot(clustering_method_plot, 0.66, 0, 0.33, 0.5)
        print(multiplot)
    }
    dev.off()
    write.csv(top_df, file=sprintf("%s/%s_mean_silhouette_top_parameter_settings.csv", qc_output_dir, sample_output_basename), quote=FALSE, row.names=FALSE)
    
    
    # mean silhouette score ranks: all spots vs hippocampus spots
    sample_output_basename <- sprintf("%s-%s-features_residual_variance", sample, qc_output_basename)
    sample_df2ggplot_all_spots <- read.csv(file=sprintf("%s/%s_mean_silhouette.csv", qc_output_dir, sample_output_basename), quote="")
    colnames(sample_df2ggplot_all_spots)[which(colnames(sample_df2ggplot_all_spots)=="rank")] <- "all_spots_rank"
    colnames(sample_df2ggplot)[which(colnames(sample_df2ggplot)=="rank")] <- "hippocampus_rank"
    
    ranks_df2ggplot <- merge(sample_df2ggplot_all_spots[, c("param_setting", "all_spots_rank")], sample_df2ggplot[, c("param_setting", "hippocampus_rank")], by="param_setting")
    sample_output_basename <- sprintf("%s-%s-features_residual_variance-hippocampus", sample, qc_output_basename)
    pdf(sprintf("%s/%s_mean_silhouette_ranks.pdf", qc_output_dir, sample_output_basename))
    p <- ggplot(ranks_df2ggplot, aes(x=all_spots_rank, y=hippocampus_rank)) +
        geom_point() +
        geom_label_repel(data=subset(ranks_df2ggplot, hippocampus_rank <= 5), aes(x=all_spots_rank, y=hippocampus_rank, label=param_setting), max.overlaps=Inf) +
        geom_label_repel(data=subset(ranks_df2ggplot, all_spots_rank <= 5), aes(x=all_spots_rank, y=hippocampus_rank, label=param_setting), max.overlaps=Inf) +
        labs(title="Mean silhouette score rank", x="All spots rank", y="Hippocampus spots ranks") +
        theme_bw() +
        theme(panel.border=element_rect(color="grey50"))
    print(p) 
    dev.off()
    
    all_samples_df2ggplot <- rbind(all_samples_df2ggplot, cbind(sample=rep(sample, dim(sample_df2ggplot)[1]), sample_df2ggplot))
    all_samples_top_df <- rbind(all_samples_top_df, cbind(sample=rep(sample, dim(top_df)[1]), top_df))
}















#############################################
# dev version of Seurat_pipeline() function #
#############################################

library(SpatialExperiment)

# parameters

spe_seurat <- spe_seurat_qc
libsize <- NA
exp_genes <- NA
norm_method <- normalization_method
features_method <- feature_selection_method
features_nb <- features_nb
features_resvar <- feature_residual_variance_threshold
dim_param <- dimensions_nb_param
knn <- nearest_neighbors_nb
clust_method <- clustering_method
clust_analysis <- cluster_analysis
test <- markers_test
logFC <- min_logFC
pct <- min_pct_spots
out_dir <- qc_output_dir

# function

# parameters
## spe_seurat: either a SpatialExperiment object for QC() function or a Seurat for all the other steps
## if libsize and exp_genes are set to NA, then do not perform QC
## if norm_method, features_method and (features_nb and features_resvar) are set to NA, then do not perform normalization and feature selection
## if dim_param, knn, clust_method and clust_analysis parameters are set to NA, then do not perform clustering
## if test, logFC and pct are set to NA, then do not perform marker genes

# QC
if (!is.na(libsize) & !is.na(exp_genes)) {
    qc_out_basename <- sprintf("QC-libsize%d_exp%d", libsize, exp_genes)
    qc_out_dir <- sprintf("%s/00-%s", out_dir, qc_out_basename)
    if (! dir.exists(qc_output_dir)) {
        dir.create(qc_output_dir, recursive=TRUE, mode="0775")
    }
    spe_seurat <- QC(spe_seurat, libsize, exp_genes, space_ranger_out_dir, qc_out_dir, qc_out_basename)
} else {
    qc_out_dir <- out_dir
}

# normalization and feature selection
if (!is.na(norm_method) & !is.na(features_method) & !(is.na(features_nb) & is.na(features_resvar))) {
    norm_out_name <- sprintf("Normalization-%s", norm_method)
    if (norm_method == "SCTransform") {
        norm_method_dir_index <- "00"
    } else {
        if (norm_method == "LogNormalize") {
            norm_method_dir_index <- "10"
        } else {
            if (norm_method == "CLR") {
                norm_method_dir_index <- "20"
            } else {
                if (norm_method == "RC") {
                    norm_method_dir_index <- "30"
                }
            }
        }
    }
    norm_out_dir <- sprintf("%s/%s-%s", qc_out_dir, norm_method_dir_index, norm_out_name)
    if (! dir.exists(norm_out_dir)) {
        dir.create(norm_out_dir, recursive=TRUE, mode="0775")
    }

    if (features_method == "HVGs") {
        features_method_dir_index <- "00"
        features_method_out_basename <- sprintf("Features-%s", features_method)
        if (!is.na(features_nb) & is.na(features_resvar)) {
            features_param_out_basename <- sprintf("nb%d", features_nb)
        } else {
            if (!is.na(features_resvar) & is.na(features_nb)) {
                features_resvar_in_dir <- sub("[.]", "_", features_resvar)
                features_param_out_basename <- sprintf("residualvariance%s", features_resvar_in_dir)
            }
        }
        features_out_basename <- sprintf("%s_%s", features_method_out_basename, features_param_out_basename)
    } else {
        if (features_method == "SVGs") {
            features_method_dir_index <- "10"
            features_method_out_basename <- sprintf("00-Features-%s", features_method)
            # TODO
        }
    }

    features_out_dir <- sprintf("%s/%s-%s/%s", norm_out_dir, features_method_dir_index, features_method_out_basename, features_param_out_basename)
    if (! dir.exists(features_out_dir)) {
        dir.create(features_out_dir, recursive=TRUE, mode="0775")
    }

    spe_seurat <- normalization_feature_selection(spe_seurat, norm_method, features_method, features_nb, features_resvar, norm_out_dir, norm_out_name, features_out_dir, features_out_basename)
    
        ##########################################
        # dev version of marker_genes() function #
        ##########################################
        
        # parameters

        seurat_obj <- spe_seurat
        norm_method <- norm_method
        features_method <- features_method
        features_nb <- features_nb
        features_resvar <- features_resvar
        norm_out_dir <- norm_out_dir
        norm_out_name <- norm_out_name
        features_out_dir <- features_out_dir
        features_out_name <- features_out_basename

        # function

        if (norm_method == "SCTransform") {
            if (features_method == "HVGs") {
                if (! is.na(features_nb) & is.na(features_resvar)) {
                    seurat_obj <- SCTransform(seurat_obj, assay="RNA", variable.features.n=features_nb)
                } else {
                    if (! is.na(features_resvar) & is.na(features_nb)) {
                        seurat_obj <- SCTransform(seurat_obj, assay="RNA", variable.features.rv.th=features_resvar)
                    }
                }
            }
            features_suffix2plot <- "SCT"
            selection_method <- "sct"
            assay_name <- "SCT"
        } else {
            # NormalizeData() function
            seurat_obj <- NormalizeData(seurat_obj, normalization.method=norm_method)
            # compute sum of counts and number of expressed genes after normalization
            counts <- as.matrix(seurat_obj@assays$RNA@data)
            nCount <- colSums(counts)
            seurat_obj$nCount_norm <- nCount
            nFeature <- apply(counts, 2, function(x) { length(x[x>0]) })
            seurat_obj$nFeature_norm <- nFeature
            features_suffix2plot <- "norm"
            selection_method <- "vst"
            assay_name <- "RNA"

            if (features_method == "HVGs") {
                if (! is.na(features_nb) & is.na(features_resvar)) {
                    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method="vst", nfeatures=features_nb)
                }
            }
        }

        # plots
        ## normalization
        pdf(sprintf("%s/%s.pdf", norm_out_dir, norm_out_name))
        ## counts per spot
        ### raw counts
        raw_plots <- normalization_plots(seurat_obj, "nCount", "RNA")
        print((wrap_plots(raw_plots[["violin"]], raw_plots[["spatial"]])) + plot_annotation(title="Raw counts"))
        ### normalized counts
        norm_plots <- normalization_plots(seurat_obj, "nCount", features_suffix2plot)
        print((wrap_plots(norm_plots[["violin"]], norm_plots[["spatial"]])) + plot_annotation(title="Normalized counts"))
        ## number of expressed genes per spot
        ### raw counts
        raw_plots <- normalization_plots(seurat_obj, "nFeature", "RNA")
        print((wrap_plots(raw_plots[["violin"]], raw_plots[["spatial"]])) + plot_annotation(title="Raw counts"))
        ### normalized counts
        norm_plots <- normalization_plots(seurat_obj, "nFeature", features_suffix2plot)
        print((wrap_plots(norm_plots[["violin"]], norm_plots[["spatial"]])) + plot_annotation(title="Normalized counts"))

        if (norm_method == "SCTransform") {
            # add comparison with log-normalization: https://satijalab.org/seurat/articles/spatial_vignette.html
            ## rerun normalization to store sctranform residuals for all genes
            seurat_obj_norm_comp <- SCTransform(seurat_obj, assay = "RNA", return.only.var.genes=FALSE, verbose=FALSE)
            ## also run standard log normalization for comparison
            seurat_objnorm_comp <- NormalizeData(seurat_obj_norm_comp, assay="RNA", verbose=FALSE)
            ## compute the correlation between the log normalized data and sctransform residuals with the number of UMIs
            seurat_obj_norm_comp <- GroupCorrelation(seurat_obj_norm_comp, group.assay="RNA", assay="RNA", slot="data", do.plot=FALSE)
            seurat_obj_norm_comp <- GroupCorrelation(seurat_obj_norm_comp, group.assay="RNA", assay="SCT", slot="scale.data", do.plot=FALSE)
            seurat_obj <- GroupCorrelation(seurat_obj, group.assay="RNA", assay="SCT", slot="scale.data", do.plot=FALSE)
            p1 <- GroupCorrelationPlot(seurat_obj_norm_comp, assay="RNA", cor="nCount_RNA_cor") + ggtitle("Log normalization") + theme(plot.title=element_text(hjust=0.5))
            p2 <- GroupCorrelationPlot(seurat_obj_norm_comp, assay="SCT", cor="nCount_RNA_cor") + ggtitle("SCTransform Normalization\nscale data for all genes") + theme(plot.title=element_text(hjust=0.5))
            p3 <- GroupCorrelationPlot(seurat_obj, assay="SCT", cor="nCount_RNA_cor") + ggtitle("SCTransform Normalization\n scale data only for the\nvariable genes") + theme(plot.title=element_text(hjust=0.5))
            print(p1 + p2)
            print(p1 + p3)
            print(p2 + p3)
        }
        dev.off()
        
        ## feature selection
        # Identify the 10 most highly variable genes
        top10 <- head(VariableFeatures(seurat_obj), 10)
        # plot variable features with and without labels
        pdf(sprintf("%s/%s.pdf", features_out_dir, features_out_name))
        plot1 <- VariableFeaturePlot(seurat_obj, selection.method=selection_method, assay=assay_name)
        plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
        print(plot2)
        dev.off()
        
        # scale data
        if (norm_method != "SCTransform" & features_method == "HVGs") {
            # linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA
            all.genes <- rownames(seurat_obj)
            seurat_obj <- ScaleData(seurat_obj, features=all.genes)
        }
}

# clustering
if (!is.na(dim_param) & !is.na(knn) & !is.na(clust_method) & !is.na(clust_analysis)) {
    if ("SCT" %in% names(spe_seurat@assays)) {
        spe_seurat_assay="SCT"
    } else {
        # NormalizeData() function
        spe_seurat_assay="RNA"
    }
    spe_seurat <- RunPCA(spe_seurat, assay=spe_seurat_assay, verbose=FALSE)

    if (dim_param == "estimation") {
        # estimate intrinsic dimension
        dim_estimate <- maxLikGlobalDimEst(spe_seurat@reductions$pca@cell.embeddings, k=10)
        dim_nb <- round(dim_estimate$dim.est)+5
        df2ggplot <- data.frame(PC=1:length(spe_seurat@reductions$pca), stdev=spe_seurat@reductions$pca@stdev)
        write.csv(df2ggplot, file=sprintf("%s/Clustering-PCs_stdev.csv", features_out_dir), quote=FALSE, row.names=TRUE)
        pdf(sprintf("%s/Clustering-PCs_stdev.pdf", features_out_dir))
        p <- ggplot(df2ggplot, aes(PC, stdev)) +
            geom_point() +
            geom_point() +
            geom_vline(xintercept = round(dim_estimate$dim.est), color = "blue") +
            geom_vline(xintercept = dim_nb, color = "red") +
            theme_bw() +
            labs(x = "Principal components", y = "Standard deviation")
        print(p)
        dev.off()
        dim_subdir <- "00-maxLikGlobalDimEst_p5"
    } else {
        dim_nb <- dim_param
        dim_subdir <- sprintf("10-%d_PCs", dim_nb)
    }

    # clustering
    if (clust_method == "Louvain") {
        clust_method_subdir <- sprintf("00-%s", clust_method)
    }

    clust_out_dir <- sprintf("%s/%s/k%d/%s", features_out_dir, dim_subdir, knn, clust_method_subdir)
    clust_out_name <- sprintf("Clustering-nbPCs%d_k%d_%s", dim_nb, knn, clust_method)
    if (! dir.exists(clust_out_dir)) {
        dir.create(clust_out_dir, recursive=TRUE, mode="0775")
    }

    spe_seurat <- clustering(spe_seurat, dim_nb, knn, clust_method, clust_analysis, clust_out_dir, clust_out_name)
}

# marker genes
if (!is.na(test) & !is.na(logFC) & !is.na(pct)) {
    if (test == "wilcox") {
        mark_test_dir_index <- "00"

    } else {
        if (mark_test == "bimod") {
            mark_test_dir_index <- "10"
        }
    }
    mark_test_out_basename <- sprintf("Markers-%s", test)
    mark_param_out_basename <- sprintf("logFC%s_minpct%s", sub("[.]", "_", logFC), sub("[.]", "_", pct))
    mark_out_basename <- sprintf("%s_%s", mark_test_out_basename, mark_param_out_basename)

    mark_out_dir <- sprintf("%s/%s-%s/%s", clust_out_dir, mark_test_dir_index, mark_test_out_basename, mark_param_out_basename)
    if (! dir.exists(mark_out_dir)) {
        dir.create(mark_out_dir, recursive=TRUE, mode="0775")
    }

    spe_seurat <- marker_genes(spe_seurat, test, logFC, pct, mark_out_dir, mark_out_basename)
}

    ##########################################
    # dev version of marker_genes() function #
    ##########################################

    # parameters

seurat_obj <- spe_seurat
test <- test
logFC <- logFC
pct <- pct
out_dir <- mark_out_dir
out_name <- mark_out_basename

    # function

library(dplyr)
# Finding differentially expressed features (cluster biomarkers): https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
## find markers for every cluster compared to all remaining spots, report only the positive ones
markers <- FindAllMarkers(seurat_obj, test.use=test, logfc.threshold=logFC, min.pct=pct, only.pos=TRUE)
write.csv(markers, file=sprintf("%s/%s.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
### visualizations
pdf(sprintf("%s/%s.pdf", out_dir, out_name))
#### expression heatmap
markers <- droplevels(markers) # removes clusters without markers to avoid errors
top10_markers <- markers %>%
    group_by(cluster) %>%
    top_n(n=10, wt=avg_log2FC)
DoHeatmap(seurat_obj, features=top10_markers$gene) + NoLegend()
#### get top6 markers for each cluster
top6_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(n=6, order_by=avg_log2FC)
for (one_cluster in levels(markers$cluster)) {
    one_cluster_markers <- top6_markers[which(top6_markers$cluster==one_cluster), "gene", drop=TRUE]
    #### violin plot
    print(VlnPlot(seurat_obj, features=one_cluster_markers) + plot_annotation(title=sprintf("Cluster: %s", one_cluster), theme=theme(plot.title=element_text(size=16))))
    #### UMAP reduced dimensions
    print(FeaturePlot(seurat_obj, features=one_cluster_markers) + plot_annotation(title=sprintf("Cluster: %s", one_cluster), theme=theme(plot.title=element_text(size=16))))
    #### spatial coordinates
    print(SpatialFeaturePlot(seurat_obj, features=one_cluster_markers, alpha=c(0.1,1)) + plot_annotation(title=sprintf("Cluster: %s", one_cluster), theme=theme(plot.title=element_text(size=14))))
}
dev.off()
    











