.libPaths(c(.libPaths(), "/home/christophe.lepriol/NeuroDev_ADD/R/r_4.1.0"))
library(ggplot2)
library(cowplot)


##############
# Parameters #
##############
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"
genome_name <- "NCBIRefSeq108_NCBIRefSeq108GTF"
load_image <- TRUE

epireg_visium_metadata_df <- data.frame(sample=c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16"),
                                       condition=factor(c(rep("SE", 4), rep("CTRL", 4), rep("SE", 4), rep("CTRL", 4))), 
                                       time=factor(c(rep(c(5, 10), 4), rep(c(20, 40, 40, 20), 2))))

# sections
sections <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")
#sections <- c("C_L1_S3", "C_L4_S15", "D_L1_S4", "D_L4_S16")
#sections <- c("A_L3_S9", "B_L3_S10")
timepoint <- NA
#timepoint <- 20
if (! is.na(timepoint)) {
    sections <- epireg_visium_metadata_df[which(epireg_visium_metadata_df$time==timepoint), "sample"]
}
nb_sections <- length(sections)

# output directory name
if (! is.na(timepoint)) {
    output_dirname <- sprintf("%dd_samples", timepoint)
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
#nb_integration_features <- 2000
#feature_residual_variance_threshold <- 1.1
#feature_residual_variance_threshold <- NA
#feature_residual_variance_threshold_vector <- c(1, 1.1, 1.3, 1.5, 2, 2.5)
#features_param_vector <- c(1.1, 1.3, 1.5, 2)
features_param_vector <- c(500, 1000, 1500, 2000, 2500)
#features_param_vector <- c(500, 1000)
#dimensions_nb_param_vector <- c(10, 15, 20)
dimensions_nb_param_vector <- c(20)
#clustering_resolution_vector <- c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2)
#clustering_resolution_vector <- c(0.4, 0.6, 0.8, 1)
#clustering_resolution_vector <- c(0.4, 1)
clustering_resolution_vector <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)
cluster_markers <- TRUE
DE_analysis <- TRUE
conserved_markers <- FALSE
integration_method <- "Seurat"



############
# Analysis #
############

integration_cluster_mean_silhouette_df <- data.frame()

# non-integrated datasets
for (features_param in features_param_vector) {
    ## normalization
    if (features_param > 10) {
        ### features_param is a number of features
        print(sprintf("features param: %d", features_param))
        normalization_output_dir <- sprintf("%s/00-SCTransform/featurenumber_%d", output_dir, features_param)
    } else {
        ### features_param is a residual variance threshold
        print(sprintf("features param: %.1f", features_param))
        normalization_output_dir <- sprintf("%s/00-SCTransform/residualvariance%s", output_dir, sub("[.]", "_", features_param))
        if (! dir.exists(normalization_output_dir)) {
            dir.create(normalization_output_dir, recursive=TRUE, mode="0775")
        }
    }
    
    for (dimensions_nb_param in dimensions_nb_param_vector) {
        ## dimensionality reduction and clustering
        print(sprintf("dimension number: %d", dimensions_nb_param))
        dimensions_output_dir <- sprintf("%s/%d_PCs", normalization_output_dir, dimensions_nb_param)
        if (! dir.exists(dimensions_output_dir)) {
            dir.create(dimensions_output_dir, recursive=TRUE, mode="0775")
        }
        
        cluster_mean_silhouette_vector <- c()
        for (clustering_resolution in clustering_resolution_vector) {
            ## dimensionality reduction and clustering
            print(sprintf("cluster resolution: %.1f", clustering_resolution))
            resolution_output_dir <- sprintf("%s/resolution%s/00-No_integration", dimensions_output_dir, sub("[.]", "_", clustering_resolution))
            if (! dir.exists(resolution_output_dir)) {
                dir.create(resolution_output_dir, recursive=TRUE, mode="0775")
            }
            
            ## silhouette
            ### cluster purity
            metadata_file <- sprintf("%s/metadata.csv", resolution_output_dir)
            metadata_df <- read.csv(metadata_file, quote="", row.names=1)
            category_colname <- "seurat_clusters"
            mean_silhouette <- mean(metadata_df[[sprintf("silhouette_%s", category_colname)]])
            cluster_mean_silhouette_vector <- c(cluster_mean_silhouette_vector, mean_silhouette)
        }
        df2ggplot <- data.frame(integration=rep("no integration", length(clustering_resolution_vector)), normalization=rep("SCTransform", length(clustering_resolution_vector)), features_param=rep(features_param, length(clustering_resolution_vector)), dimensions=rep(dimensions_nb_param, length(clustering_resolution_vector)), resolution=clustering_resolution_vector, silhouette_clusters=cluster_mean_silhouette_vector)
        integration_cluster_mean_silhouette_df <- rbind(integration_cluster_mean_silhouette_df, df2ggplot)
    }
}


# integrated datasets
integration_features_param_vector <- c(500, 1000, 1500, 2000)
## need to set maxSize to avoid error when identifying conserved markers: Error in getGlobalsAndPackages(expr, envir = envir, globals = globals): The total size of the 3 globals exported for future expression (‘FUN()’) is 1.38 GiB.. This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize').
#options(future.globals.maxSize=2000 * 1024^2) # set allowed size to 2K MiB
for (features_param in integration_features_param_vector) {    
    ## normalization
    if (features_param > 10) {
        ### features_param is a number of features
        print(sprintf("features param: %d", features_param))
        normalization_output_dir <- sprintf("%s/00-SCTransform/featurenumber_%d", output_dir, features_param)
        if (! dir.exists(normalization_output_dir)) {
            dir.create(normalization_output_dir, recursive=TRUE, mode="0775")
        }
    } else {
        ### features_param is a residual variance threshold
        print(sprintf("features param: %.1f", features_param))
        normalization_output_dir <- sprintf("%s/00-SCTransform/residualvariance%s", output_dir, sub("[.]", "_", features_param))
        if (! dir.exists(normalization_output_dir)) {
            dir.create(normalization_output_dir, recursive=TRUE, mode="0775")
        }
    }
    
    ## integration
    for (dimensions_nb_param in dimensions_nb_param_vector) {
        ## dimensionality reduction and clustering
        print(sprintf("dimension number: %d", dimensions_nb_param))
        dimensions_output_dir <- sprintf("%s/%d_PCs", normalization_output_dir, dimensions_nb_param)
        if (! dir.exists(dimensions_output_dir)) {
            dir.create(dimensions_output_dir, recursive=TRUE, mode="0775")
        }

        cluster_mean_silhouette_vector <- c()
        for (clustering_resolution in clustering_resolution_vector) {
            ## dimensionality reduction and clustering
            print(sprintf("cluster resolution: %.1f", clustering_resolution))
            resolution_output_dir <- sprintf("%s/resolution%s/10-Integration", dimensions_output_dir, sub("[.]", "_", clustering_resolution))
            if (! dir.exists(resolution_output_dir)) {
                dir.create(resolution_output_dir, recursive=TRUE, mode="0775")
            }

            ## silhouette
            ### cluster purity
            metadata_file <- sprintf("%s/metadata.csv", resolution_output_dir)
            metadata_df <- read.csv(metadata_file, quote="", row.names=1)
            category_colname <- "seurat_clusters"
            mean_silhouette <- mean(metadata_df[[sprintf("silhouette_%s", category_colname)]])
            cluster_mean_silhouette_vector <- c(cluster_mean_silhouette_vector, mean_silhouette)
        }
        df2ggplot <- data.frame(integration=rep(integration_method, length(clustering_resolution_vector)), normalization=rep("SCTransform", length(clustering_resolution_vector)), features_param=rep(features_param, length(clustering_resolution_vector)), dimensions=rep(dimensions_nb_param, length(clustering_resolution_vector)), resolution=clustering_resolution_vector, silhouette_clusters=cluster_mean_silhouette_vector)
        integration_cluster_mean_silhouette_df <- rbind(integration_cluster_mean_silhouette_df, df2ggplot)
    }
}
write.csv(integration_cluster_mean_silhouette_df, file=sprintf("%s/00-SCTransform/featurenumber_cluster_mean_silhouette.csv", output_dir), quote=FALSE, row.names=FALSE)
pdf(sprintf("%s/00-SCTransform/featurenumber_cluster_mean_silhouette.pdf", output_dir))
integration_cluster_mean_silhouette_df$integration <- factor(integration_cluster_mean_silhouette_df$integration, levels=c("no integration", "Seurat"))
integration_cluster_mean_silhouette_df$features_param <- as.factor(integration_cluster_mean_silhouette_df$features_param)
integration_cluster_mean_silhouette_df$resolution <- as.factor(integration_cluster_mean_silhouette_df$resolution)
p <- ggplot(integration_cluster_mean_silhouette_df, aes(x=resolution, y=silhouette_clusters, fill=resolution)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(integration~features_param) +
    labs(title="Cluster mean silhouette", x="Number of features", y="Integration method", fill="Resolution") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=60, hjust=1)) +
    theme(panel.border=element_rect(color="grey50"))
print(p)

p <- ggplot(integration_cluster_mean_silhouette_df, aes(x=features_param, y=silhouette_clusters, fill=features_param)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(integration~resolution) +
    labs(title="Cluster mean silhouette", x="Resolution", y="Integration method", fill="Number of\nfeatures") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=60, hjust=1)) +
    theme(panel.border=element_rect(color="grey50"))
print(p)
dev.off()


