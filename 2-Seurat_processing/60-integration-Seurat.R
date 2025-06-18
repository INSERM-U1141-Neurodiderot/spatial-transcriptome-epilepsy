.libPaths(c(.libPaths(), "/home/christophe.lepriol/NeuroDev_ADD/R/r_4.1.0"))
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(cluster) # silhouette() function
library(cowplot)
#library(kBET)
library(harmony)

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
#timepoint <- NA
timepoint <- 20
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
#features_param_vector <- c(500, 1000, 1500, 2000, 2500)
features_param_vector <- c(500, 1000)
#dimensions_nb_param_vector <- c(10, 15, 20)
dimensions_nb_param_vector <- c(15, 20)
#dimensions_nb_param_vector <- c(20)
#clustering_resolution_vector <- c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2)
#clustering_resolution_vector <- c(0.4, 0.6, 0.8, 1)
#clustering_resolution_vector <- c(0.4, 1)
#clustering_resolution_vector <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)
clustering_resolution_vector <- c(0.1, 0.2)
cluster_markers <- FALSE
DE_analysis <- FALSE
conserved_markers <- FALSE
#integration_method <- "Seurat"

# number of samples per page when using SpatialFeaturePlot() function
nb_SpatialFeaturePlot_per_page <- 2
# number of samples per page when using SpatialFeaturePlot() function
nb_cluster_markers_to_plot <- 4

nb_DEGs_to_plot <- 30
nb_DEGs_per_page <- 3



#############
# Functions #
#############
load_data <- function(sample, st_dir, sr_dir, add_image) {
    # expression matrices
    ## ST Pipeline
    st_pipeline_matrix_file <- sprintf("%s/%s/10-Pipeline/%s_stdata.tsv", st_dir, sample, sample)
    st_pipeline_matrix <- read.table(st_pipeline_matrix_file, sep="\t", header=TRUE, quote="", row.names=1)
    # Space Ranger output directory for tissue positions list and image
    space_ranger_output_dir <- sprintf("%s/%s/10-Pipeline/outs", sr_dir, sample)

    # build SpatialExperiment object
    spatial_coordinates_file <- file.path(space_ranger_output_dir, "spatial", "tissue_positions_list.csv")
    spatial_coordinates_df <- read.csv(spatial_coordinates_file, header=FALSE, quote="")
    colnames(spatial_coordinates_df) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    rownames(spatial_coordinates_df) <- spatial_coordinates_df$barcode
    nb_spots <- dim(spatial_coordinates_df)[1]
    ## add lame, zone and condition
    sample_strsplit <- unlist(strsplit(sample, "_"))
    zone <- sample_strsplit[1]
    lame <- sample_strsplit[2]
    condition <- ifelse(zone=="A" | zone=="C", "SE", "CTRL")
    spatial_coordinates_df <- cbind(spatial_coordinates_df, lame=rep(lame, nb_spots), zone=rep(zone, nb_spots), condition=rep(condition, nb_spots))
    spatial_coordinates_df$zone <- as.factor(spatial_coordinates_df$zone)
    spatial_coordinates_df$lame <- as.factor(spatial_coordinates_df$lame)
    spatial_coordinates_df$condition <- as.factor(spatial_coordinates_df$condition)
    ## matrix: convert coordinates to barcodes
    barcodes <- unlist(lapply(rownames(st_pipeline_matrix), function(x, coord2barcode=spatial_coordinates_df) {
        coordinates <- unlist(strsplit(x, "x"))
        return(coord2barcode[which(coord2barcode$array_col==as.integer(coordinates[1])-1 & coord2barcode$array_row==as.integer(coordinates[2])-1), "barcode"])
    }))
    rownames(st_pipeline_matrix) <- barcodes
    
    # create Seurat object
    seurat_obj <- CreateSeuratObject(counts=t(st_pipeline_matrix), project=sample, meta.data=spatial_coordinates_df)
    ## add image
    if (add_image) {
        img <- Read10X_Image(image.dir=file.path(space_ranger_output_dir, "spatial"))
        img@assay <- Assays(seurat_obj) # set assay identical to Seurat object
        Key(img) <- sprintf("%s_", gsub("_", "", sample)) # add key for SpatialFeaturePlot() and SpatialPlot() functions: Keys should be one or more alphanumeric characters followed by an underscore
        img_list <- list(img) # the image must be stored in a list
        names(img_list) <- sample # change image name in the list: used as title by SpatialFeaturePlot() function for example
        seurat_obj@images <- img_list
    }
    ## keep spot over tissue
    seurat_obj <- seurat_obj[, seurat_obj$in_tissue==1]
    
    return(seurat_obj)
}

compute_silhouette <- function(seurat_obj, metadata_colname, dist_matrix) {
    categories <- seurat_obj@meta.data[[metadata_colname]]
    silhouette_obj <- silhouette(as.numeric(categories), dist=dist_matrix)
    #seurat_obj@meta.data$silhouette <- silhouette_obj[,3]
    #colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data)=="silhouette")] <- sprintf("silhouette_%s", metadata_colname)
    seurat_obj@meta.data[[sprintf("silhouette_%s", metadata_colname)]] <- silhouette_obj[,3]
    return(seurat_obj)
}


############
# Analysis #
############

integration_cluster_mean_silhouette_df <- data.frame()

# load data
print("load data")
st_pipeline_dir <- sprintf("%s/10-ST_analysis/00-ST_Pipeline/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples", work_dir, genome_name)
space_ranger_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples", work_dir, genome_name)
seurat_objects <- lapply(sections, load_data, st_dir=st_pipeline_dir, sr_dir=space_ranger_dir, add_image=load_image)
names(seurat_objects) <- sections

# quality control
## merge datasets
print("merge datasets")
seurat_objects_merge <- merge(seurat_objects[[1]], y=unlist(seurat_objects)[2:nb_sections], add.cell.ids=sections)
if (load_image) {
    print("quality control")
    ## statistics: number of counts, number of features
    pdf(sprintf("%s/QC.pdf", output_dir))
    VlnPlot(seurat_objects_merge, features=c("nCount_RNA", "nFeature_RNA"), pt.size=0.1, ncol=2) + NoLegend()
    ### SpatialFeaturePlot: impossible to change legend parameters due to combine=TRUE
    #### SpatialFeaturePlot(seurat_objects_merge, features=c("nCount_RNA", "nFeature_RNA")) + theme(legend.title=element_text(size=8), legend.text=element_text(size=6)): the parameters are only applied to the last of plot
    #### only solution to make the legend readable: plot only 2 samples per page
    for (i in 1:(nb_sections/nb_SpatialFeaturePlot_per_page)) {
        print(SpatialFeaturePlot(seurat_objects_merge, features=c("nCount_RNA", "nFeature_RNA"), images=names(seurat_objects_merge@images)[((i-1)*nb_SpatialFeaturePlot_per_page+1):(i*nb_SpatialFeaturePlot_per_page)]))
    }
    ## filter
    ### select all spots with less than 25% mitochondrial reads, less than 20% hb-reads and more than 500 detected genes
    seurat_objects_merge <- seurat_objects_merge[, seurat_objects_merge$nFeature_RNA > 500]
    for (i in 1:(nb_sections/nb_SpatialFeaturePlot_per_page)) {
        print(SpatialFeaturePlot(seurat_objects_merge, features=c("nCount_RNA", "nFeature_RNA"), images=names(seurat_objects_merge@images)[((i-1)*nb_SpatialFeaturePlot_per_page+1):(i*nb_SpatialFeaturePlot_per_page)]))
    }
    dev.off()
}

# non-integrated datasets
for (features_param in features_param_vector) {
    ## normalization
    if (features_param > 10) {
        ### features_param is a number of features
        print(sprintf("features param: %d", features_param))
        normalization_output_dir <- sprintf("%s/00-SCTransform/featurenumber_%d", output_dir, features_param)
        if (! dir.exists(normalization_output_dir)) {
            dir.create(normalization_output_dir, recursive=TRUE, mode="0775")
        }
        seurat_objects_merge <- SCTransform(seurat_objects_merge, assay="RNA", verbose=TRUE, method="poisson", variable.features.n=features_param, return.only.var.genes=FALSE)
    } else {
        ### features_param is a residual variance threshold
        print(sprintf("features param: %.1f", features_param))
        normalization_output_dir <- sprintf("%s/00-SCTransform/residualvariance%s", output_dir, sub("[.]", "_", features_param))
        if (! dir.exists(normalization_output_dir)) {
            dir.create(normalization_output_dir, recursive=TRUE, mode="0775")
        }
        seurat_objects_merge <- SCTransform(seurat_objects_merge, assay="RNA", verbose=TRUE, method="poisson", variable.features.n=NULL, variable.features.rv.th=features_param, return.only.var.genes=FALSE)
    }
    
    for (dimensions_nb_param in dimensions_nb_param_vector) {
        ## dimensionality reduction and clustering
        print(sprintf("dimension number: %d", dimensions_nb_param))
        dimensions_output_dir <- sprintf("%s/00-No_integration/%d_PCs", normalization_output_dir, dimensions_nb_param)
        if (! dir.exists(dimensions_output_dir)) {
            dir.create(dimensions_output_dir, recursive=TRUE, mode="0775")
        }
        seurat_objects_merge <- RunPCA(seurat_objects_merge, assay="SCT", verbose=FALSE)
        seurat_objects_merge <- FindNeighbors(seurat_objects_merge, reduction="pca", dims=1:dimensions_nb_param)
        
        cluster_mean_silhouette_vector <- c()
        for (clustering_resolution in clustering_resolution_vector) {
            ## dimensionality reduction and clustering
            print(sprintf("cluster resolution: %.1f", clustering_resolution))
            resolution_output_dir <- sprintf("%s/resolution%s", dimensions_output_dir, sub("[.]", "_", clustering_resolution))
            if (! dir.exists(resolution_output_dir)) {
                dir.create(resolution_output_dir, recursive=TRUE, mode="0775")
            }
            seurat_objects_merge <- FindClusters(seurat_objects_merge, resolution=clustering_resolution, verbose=FALSE)
            seurat_objects_merge <- RunUMAP(seurat_objects_merge, reduction="pca", dims=1:dimensions_nb_param)
            
            ## silhouette
            distance_matrix <- dist(Embeddings(seurat_objects_merge$pca))
            ### cluster purity, lame mixing, zone mixing, condition, sample
            for (category_colname in c("seurat_clusters", "lame", "zone", "condition", "orig.ident")) {
                seurat_objects_merge@meta.data[[category_colname]] <- as.factor(seurat_objects_merge@meta.data[[category_colname]])
                if (length(levels(seurat_objects_merge@meta.data[[category_colname]])) > 1) {
                    seurat_objects_merge <- compute_silhouette(seurat_objects_merge, category_colname, distance_matrix)
                    mean_silhouette <- mean(seurat_objects_merge@meta.data[[sprintf("silhouette_%s", category_colname)]])
                    print(sprintf("mean silhouette for %s: %f", category_colname, mean_silhouette))
                    if (category_colname == "seurat_clusters") {
                        cluster_mean_silhouette_vector <- c(cluster_mean_silhouette_vector, mean_silhouette)
                    }
                }
            }
            write.csv(seurat_objects_merge@meta.data, file=sprintf("%s/metadata.csv", resolution_output_dir), quote=FALSE, row.names=TRUE)
            
            ## kBET
            ### clusters
            #kBET_values <- kBET(Embeddings(seurat_objects_merge$pca), seurat_objects_merge@meta.data$seurat_clusters)
            ### lame
            #### 50 PCs
            #kBET_values_lame <- kBET(Embeddings(seurat_objects_merge$pca), seurat_objects_merge@meta.data$lame)
            #### 20 PCs
            #kBET_values_lame_20PCs <- kBET(Embeddings(seurat_objects_merge$pca)[,1:20], seurat_objects_merge@meta.data$lame)
            
            
            

            ## plot clusters onto UMAP or onto the tissue section
            pdf(sprintf("%s/clustering.pdf", resolution_output_dir))
            dimplot_list <- list()
            for (category_colname in c("seurat_clusters", "lame", "zone", "condition", "orig.ident")) {
                p <- DimPlot(seurat_objects_merge, reduction="umap", group.by=category_colname, label=TRUE)
                if (length(levels(seurat_objects_merge@meta.data[[category_colname]])) > 1) {
                    p <- p + ggtitle(sprintf("Mean silhouette for %s: %f", category_colname, mean(seurat_objects_merge@meta.data[[sprintf("silhouette_%s", category_colname)]]))) + theme(plot.title=element_text(size=14))
                }
                dimplot_list[[category_colname]] <- p
            }
            print(dimplot_list[["seurat_clusters"]] / dimplot_list[["lame"]] + plot_annotation(title="Non-integrated datasets", theme=theme(plot.title=element_text(size=16)))) 
            print(dimplot_list[["seurat_clusters"]] / dimplot_list[["zone"]] + plot_annotation(title="Non-integrated datasets", theme=theme(plot.title=element_text(size=16))))
            print(dimplot_list[["seurat_clusters"]] / dimplot_list[["condition"]] + plot_annotation(title="Non-integrated datasets", theme=theme(plot.title=element_text(size=16))))
            print(dimplot_list[["seurat_clusters"]] / dimplot_list[["orig.ident"]] + plot_annotation(title="Non-integrated datasets", theme=theme(plot.title=element_text(size=16))))
            if (load_image) {
                for (i in 1:(nb_sections/nb_SpatialFeaturePlot_per_page)) {
                    print(SpatialDimPlot(seurat_objects_merge, images=names(seurat_objects_merge@images)[((i-1)*nb_SpatialFeaturePlot_per_page+1):(i*nb_SpatialFeaturePlot_per_page)], label=TRUE, label.size=2) + plot_annotation(title="Non-integrated datasets", theme=theme(plot.title=element_text(size=18))))
                }
                for (i in 1:nb_sections) {
                    print(SpatialDimPlot(seurat_objects_merge, images=names(seurat_objects_merge@images)[i], label=TRUE, label.size=3) + plot_annotation(title=sprintf("Non-integrated datasets: %s", names(seurat_objects_merge@images)[i]), theme=theme(plot.title=element_text(size=16))))
                }
            }
            dev.off()
        }
        df2ggplot <- data.frame(integration=rep("no integration", length(clustering_resolution_vector)), normalization=rep("SCTransform", length(clustering_resolution_vector)), features_param=rep(features_param, length(clustering_resolution_vector)), dimensions=rep(dimensions_nb_param, length(clustering_resolution_vector)), resolution=clustering_resolution_vector, silhouette_clusters=cluster_mean_silhouette_vector)
        integration_cluster_mean_silhouette_df <- rbind(integration_cluster_mean_silhouette_df, df2ggplot)
    }
}
rm(seurat_objects_merge)
rm(distance_matrix)
gc()


# integrated datasets
## need to set maxSize to avoid error when identifying conserved markers: Error in getGlobalsAndPackages(expr, envir = envir, globals = globals): The total size of the 3 globals exported for future expression (‘FUN()’) is 1.38 GiB.. This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize').
#options(future.globals.maxSize=2000 * 1024^2) # set allowed size to 2K MiB
for (features_param in features_param_vector) {
    for (integration_method in c("Seurat", "Harmony")) {
        # load data
        print("load data")
        st_pipeline_dir <- sprintf("%s/10-ST_analysis/00-ST_Pipeline/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples", work_dir, genome_name)
        space_ranger_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples", work_dir, genome_name)
        seurat_objects <- lapply(sections, load_data, st_dir=st_pipeline_dir, sr_dir=space_ranger_dir, add_image=load_image)
        names(seurat_objects) <- sections
        
        ## normalization
        if (features_param > 10) {
            ### features_param is a number of features
            print(sprintf("features param: %d", features_param))
            normalization_output_dir <- sprintf("%s/00-SCTransform/featurenumber_%d", output_dir, features_param)
            if (! dir.exists(normalization_output_dir)) {
                dir.create(normalization_output_dir, recursive=TRUE, mode="0775")
            }
            
            if (integration_method == "Seurat") {
                # Seurat integration methods require a list of Seurat object
                seurat_objects <- lapply(seurat_objects, FUN=SCTransform, variable.features.n=features_param, return.only.var.genes=FALSE)
            } else {
                if (integration_method == "Harmony") {
                    # Harmony integration function requires a merged Seurat object
                    seurat_objects_merge <- merge(seurat_objects[[1]], y=unlist(seurat_objects)[2:nb_sections], add.cell.ids=sections)
                    seurat_objects_merge <- SCTransform(seurat_objects_merge, variable.features.n=features_param, return.only.var.genes=FALSE)
                    rm(seurat_objects)
                    gc()
                }
            }
            integration_features_nb <- features_param
        } else {
            ### features_param is a residual variance threshold
            print(sprintf("features param: %.1f", features_param))
            normalization_output_dir <- sprintf("%s/00-SCTransform/residualvariance%s", output_dir, sub("[.]", "_", features_param))
            if (! dir.exists(normalization_output_dir)) {
                dir.create(normalization_output_dir, recursive=TRUE, mode="0775")
            }
            
            if (integration_method == "Seurat") {
                # Seurat integration methods require a list of Seurat object
                seurat_objects <- lapply(seurat_objects, FUN=SCTransform, variable.features.n=NULL, variable.features.rv.th=features_param, return.only.var.genes=FALSE)
                ## get minimum number of variable features
                min_variable_features <- length(seurat_objects[[sections[1]]]@assays$SCT@var.features)
                for (sample in sections[2:length(sections)]) {
                    if (length(seurat_objects[[sample]]@assays$SCT@var.features) < min_variable_features) {
                        min_variable_features <- length(seurat_objects[[sample]]@assays$SCT@var.features)
                    }
                }
            } else {
                if (integration_method == "Harmony") {
                    # Harmony integration function requires a merged Seurat object
                    seurat_objects_merge <- merge(seurat_objects[[1]], y=unlist(seurat_objects)[2:nb_sections], add.cell.ids=sections)
                    seurat_objects_merge <- SCTransform(seurat_objects_merge, variable.features.n=NULL, variable.features.rv.th=features_param, return.only.var.genes=FALSE)
                    ## TODO: ## get minimum number of variable features
                    rm(seurat_objects)
                    gc()
                }
            }
            integration_features_nb <- min_variable_features
        }
        
        ## integration: Seurat
        if (integration_method == "Seurat") {
            integration_output_dir <- sprintf("%s/10-Seurat", normalization_output_dir)
            if (! dir.exists(integration_output_dir)) {
                dir.create(integration_output_dir, recursive=TRUE, mode="0775")
            }
            features <- SelectIntegrationFeatures(seurat_objects, nfeatures=integration_features_nb)
            seurat_objects <- PrepSCTIntegration(object.list=seurat_objects, anchor.features=features)
            int.anchors <- FindIntegrationAnchors(object.list=seurat_objects, normalization.method="SCT", anchor.features=features)
            rm(seurat_objects)
            gc()
            seurat_objects_integrated <- IntegrateData(anchorset=int.anchors, normalization.method="SCT")
            rm(int.anchors)
            gc()
            reduction2use <- "pca"
        }
        
        if (integration_method == "Harmony") {
            seurat_objects_merge <- RunPCA(seurat_objects_merge, verbose=FALSE)
        }
        
        for (dimensions_nb_param in dimensions_nb_param_vector) {
            ## integration: Harmony
            if (integration_method == "Harmony") {
                integration_output_dir <- sprintf("%s/11-Harmony", normalization_output_dir)
                if (! dir.exists(integration_output_dir)) {
                    dir.create(integration_output_dir, recursive=TRUE, mode="0775")
                }
                seurat_objects_integrated <- RunHarmony(seurat_objects_merge, group.by.vars="lame", dims=1:dimensions_nb_param)
                reduction2use <- "harmony"
            }
            
            ## dimensionality reduction and clustering
            print(sprintf("dimension number: %d", dimensions_nb_param))
            dimensions_output_dir <- sprintf("%s/%d_PCs", integration_output_dir, dimensions_nb_param)
            if (! dir.exists(dimensions_output_dir)) {
                dir.create(dimensions_output_dir, recursive=TRUE, mode="0775")
            }

            if (integration_method == "Seurat") {
                seurat_objects_integrated <- RunPCA(seurat_objects_integrated, assay="integrated", verbose=FALSE)
            }
            seurat_objects_integrated <- FindNeighbors(seurat_objects_integrated, reduction=reduction2use, dims=1:dimensions_nb_param)

            cluster_mean_silhouette_vector <- c()
            for (clustering_resolution in clustering_resolution_vector) {
                ## dimensionality reduction and clustering
                print(sprintf("cluster resolution: %.1f", clustering_resolution))
                resolution_output_dir <- sprintf("%s/resolution%s", dimensions_output_dir, sub("[.]", "_", clustering_resolution))
                if (! dir.exists(resolution_output_dir)) {
                    dir.create(resolution_output_dir, recursive=TRUE, mode="0775")
                }

                # set default assay to integrated before clustering
                if (integration_method == "Seurat") {
                    if (DefaultAssay(seurat_objects_integrated) != "integrated") {
                        DefaultAssay(seurat_objects_integrated) <- "integrated"
                    }
                }
                seurat_objects_integrated <- FindClusters(seurat_objects_integrated, resolution=clustering_resolution, verbose=FALSE)
                seurat_objects_integrated <- RunUMAP(seurat_objects_integrated, reduction=reduction2use, dims=1:dimensions_nb_param)

                ## silhouette
                distance_matrix <- dist(Embeddings(seurat_objects_integrated[[reduction2use]]))
                ### cluster purity, lame mixing, zone mixing, condition, sample
                for (category_colname in c("seurat_clusters", "lame", "zone", "condition", "orig.ident")) {
                    seurat_objects_integrated@meta.data[[category_colname]] <- as.factor(seurat_objects_integrated@meta.data[[category_colname]])
                    if (length(levels(seurat_objects_integrated@meta.data[[category_colname]])) > 1) {
                        seurat_objects_integrated <- compute_silhouette(seurat_objects_integrated, category_colname, distance_matrix)
                        mean_silhouette <- mean(seurat_objects_integrated@meta.data[[sprintf("silhouette_%s", category_colname)]])
                        print(sprintf("mean silhouette for %s: %f", category_colname, mean_silhouette))
                        if (category_colname == "seurat_clusters") {
                            cluster_mean_silhouette_vector <- c(cluster_mean_silhouette_vector, mean_silhouette)
                        }
                    }
                }
                write.csv(seurat_objects_integrated@meta.data, file=sprintf("%s/metadata.csv", resolution_output_dir), quote=FALSE, row.names=TRUE)

                ## kBET
                ### clusters
                #kBET_values <- kBET(Embeddings(seurat_objects_merge$pca), seurat_objects_merge@meta.data$seurat_clusters)
                ### lame
                #### 50 PCs
                #kBET_values_lame <- kBET(Embeddings(seurat_objects_merge$pca), seurat_objects_merge@meta.data$lame)
                #### 20 PCs
                #kBET_values_inetgrated_lame_20PCs <- kBET(Embeddings(seurat_objects_integrated$pca)[,1:20], seurat_objects_integrated@meta.data$lame)


                ## plot clusters onto UMAP or onto the tissue section
                pdf(sprintf("%s/clustering.pdf", resolution_output_dir))
                integration_dimplot_list <- list()
                for (category_colname in c("seurat_clusters", "lame", "zone", "condition", "orig.ident")) {
                    p <- DimPlot(seurat_objects_integrated, reduction="umap", group.by=category_colname, label=TRUE)
                    if (length(levels(seurat_objects_integrated@meta.data[[category_colname]])) > 1) {
                        p <- p + ggtitle(sprintf("Mean silhouette for %s: %f", category_colname, mean(seurat_objects_integrated@meta.data[[sprintf("silhouette_%s", category_colname)]]))) + theme(plot.title=element_text(size=14))
                    }
                    integration_dimplot_list[[category_colname]] <- p
                }
                print(integration_dimplot_list[["seurat_clusters"]] / integration_dimplot_list[["lame"]] + plot_annotation(title="Integrated datasets", theme=theme(plot.title=element_text(size=16)))) 
                print(integration_dimplot_list[["seurat_clusters"]] / integration_dimplot_list[["zone"]] + plot_annotation(title="Integrated datasets", theme=theme(plot.title=element_text(size=16))))
                print(integration_dimplot_list[["seurat_clusters"]] / integration_dimplot_list[["condition"]] + plot_annotation(title="Integrated datasets", theme=theme(plot.title=element_text(size=16))))
                print(integration_dimplot_list[["seurat_clusters"]] / integration_dimplot_list[["orig.ident"]] + plot_annotation(title="Integrated datasets", theme=theme(plot.title=element_text(size=16))))
                if (load_image) {
                    for (i in 1:(nb_sections/nb_SpatialFeaturePlot_per_page)) {
                        print(SpatialDimPlot(seurat_objects_integrated, images=names(seurat_objects_integrated@images)[((i-1)*nb_SpatialFeaturePlot_per_page+1):(i*nb_SpatialFeaturePlot_per_page)], label=TRUE, label.size=2) + plot_annotation(title="Integrated datasets", theme=theme(plot.title=element_text(size=18))))
                    }
                    for (i in 1:nb_sections) {
                        print(SpatialDimPlot(seurat_objects_integrated, images=names(seurat_objects_integrated@images)[i], label=TRUE, label.size=3) + plot_annotation(title=sprintf("Integrated datasets: %s", names(seurat_objects_integrated@images)[i]), theme=theme(plot.title=element_text(size=16))))
                    }
                }
                dev.off()

                if (cluster_markers) {
                    # find markers for every cluster compared to all remaining spots, report only the positive ones
                    markers_output_dir <- sprintf("%s/00-Cluster_markers", resolution_output_dir)
                    if (! dir.exists(markers_output_dir)) {
                        dir.create(markers_output_dir, recursive=TRUE, mode="0775")
                    }
                    default_assay <- "SCT"
                    DefaultAssay(seurat_objects_integrated) <- default_assay
                    markers <- FindAllMarkers(seurat_objects_integrated, logfc.threshold=0.25, min.pct=0.25, only.pos=TRUE)
                    write.csv(markers, file=sprintf("%s/cluster_marker_genes.csv", markers_output_dir), quote=FALSE, row.names=TRUE)
                    pdf(sprintf("%s/cluster_marker_genes.pdf", markers_output_dir))
                    ## expression heatmap
                    top10_markers <- markers %>%
                        group_by(cluster) %>%
                        top_n(n=10, wt=avg_log2FC)
                    print(DoHeatmap(seurat_objects_integrated, features=top10_markers$gene) + NoLegend())
                    ## get top markers for each cluster
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
                }

                if (DE_analysis) {
                    # differentially expressed genes across conditions
                    DE_output_dir <- sprintf("%s/10-DEGs", resolution_output_dir)
                    if (! dir.exists(DE_output_dir)) {
                        dir.create(DE_output_dir, recursive=TRUE, mode="0775")
                    }
                    default_assay <- "SCT"
                    DefaultAssay(seurat_objects_integrated) <- default_assay
                    ## average expression plots
                    theme_set(theme_cowplot())
                    Idents(seurat_objects_integrated) <- seurat_objects_integrated$seurat_clusters
                    seurat_objects_integrated.all_clusters.avg <- data.frame()
                    pdf(sprintf("%s/DEGs_all_clusters_average_expression.pdf", DE_output_dir))
                    for (one_cluster in levels(seurat_objects_integrated$seurat_clusters)) {
                        print(sprintf("cluster: %s", one_cluster))
                        seurat_objects_integrated.one_cluster <- subset(seurat_objects_integrated, idents=one_cluster)
                        ## set class to conditions
                        Idents(seurat_objects_integrated.one_cluster) <- "condition"
                        seurat_objects_integrated.one_cluster.avg <- as.data.frame(log1p(AverageExpression(seurat_objects_integrated.one_cluster, verbose=FALSE)[[default_assay]]))
                        seurat_objects_integrated.one_cluster.avg$gene <- rownames(seurat_objects_integrated.one_cluster.avg)
                        seurat_objects_integrated.one_cluster.avg$diff <- seurat_objects_integrated.one_cluster.avg$SE - seurat_objects_integrated.one_cluster.avg$CTRL
                        seurat_objects_integrated.all_clusters.avg <- rbind(seurat_objects_integrated.all_clusters.avg, cbind(seurat_objects_integrated.one_cluster.avg, cluster=rep(one_cluster, dim(seurat_objects_integrated.one_cluster.avg)[1])))
                        genes.to.label <- c(rownames(seurat_objects_integrated.one_cluster.avg[order(seurat_objects_integrated.one_cluster.avg$diff, decreasing=TRUE),])[1:10], tail(rownames(seurat_objects_integrated.one_cluster.avg[order(seurat_objects_integrated.one_cluster.avg$diff, decreasing=TRUE),]), n=10))
                        p <- ggplot(seurat_objects_integrated.one_cluster.avg, aes(CTRL, SE)) + geom_point() + ggtitle(sprintf("Cluster %s", one_cluster))
                        p <- LabelPoints(p, points=genes.to.label, repel=TRUE)
                        print(p)
                    }
                    dev.off()
                    write.csv(seurat_objects_integrated.all_clusters.avg, file=sprintf("%s/DEGs_all_clusters_average_expression.csv", DE_output_dir), quote=FALSE, row.names=TRUE)

                    ## differentially expressed genes in different conditions for cells of the same type
                    ### column in metadata table to hold both the cell type and stimulation information and switch the current ident to that column
                    seurat_objects_integrated$seurat_clusters.condition <- paste(Idents(seurat_objects_integrated), seurat_objects_integrated$condition, sep="_")
                    Idents(seurat_objects_integrated) <- seurat_objects_integrated$seurat_clusters.condition
                    ### differentially expressed genes between SE and CTRL spots
                    seurat_objects_integrated.DE.all_clusters <- data.frame()
                    for (one_cluster in levels(seurat_objects_integrated$seurat_clusters)) {
                        print(sprintf("cluster: %s", one_cluster))
                        seurat_objects_integrated.DE.one_cluster <- FindMarkers(seurat_objects_integrated, ident.1=sprintf("%s_SE", one_cluster), ident.2=sprintf("%s_CTRL", one_cluster), verbose=FALSE)
                        write.csv(seurat_objects_integrated.DE.one_cluster, file=sprintf("%s/DEGs_cluster_%s.csv", DE_output_dir, one_cluster), quote=FALSE, row.names=TRUE)
                        seurat_objects_integrated.DE.all_clusters <- rbind(seurat_objects_integrated.DE.all_clusters, cbind(seurat_objects_integrated.DE.one_cluster, cluster=rep(one_cluster, dim(seurat_objects_integrated.DE.one_cluster)[1])))
                        #### visualizations
                        pdf(sprintf("%s/DEGs_cluster_%s.pdf", DE_output_dir, one_cluster))
                        for(i in seq(1, nb_DEGs_to_plot, nb_DEGs_per_page)) {
                            print(FeaturePlot(seurat_objects_integrated, features=rownames(seurat_objects_integrated.DE.one_cluster)[i:(i+nb_DEGs_per_page-1)], split.by="condition", max.cutoff=3, cols=c("grey", "red")))
                            plots <- VlnPlot(seurat_objects_integrated, features=rownames(seurat_objects_integrated.DE.one_cluster)[i:(i+nb_DEGs_per_page-1)], split.by="condition", group.by="seurat_clusters", pt.size=0, combine=FALSE)
                            print(wrap_plots(plots=plots, ncol=1))
                        }
                        dev.off()
                    }
                    write.csv(seurat_objects_integrated.DE.all_clusters, file=sprintf("%s/DEGs_all_clusters.csv", DE_output_dir), quote=FALSE, row.names=TRUE)
                    seurat_objects_integrated.DE.all_clusters.ordered <- seurat_objects_integrated.DE.all_clusters[order(seurat_objects_integrated.DE.all_clusters$p_val_adj),]
                }

                if (conserved_markers) {
                    # conserved cell type markers
                    DefaultAssay(seurat_objects_integrated) <- "SCT"
                    conserved_markers_output_dir <- sprintf("%s/20-Conserved_markers", resolution_output_dir)
                    if (! dir.exists(conserved_markers_output_dir)) {
                        dir.create(conserved_markers_output_dir, recursive=TRUE, mode="0775")
                    }
                    pdf(sprintf("%s/conserved_markers_umap_dotplots.pdf", conserved_markers_output_dir))
                    for (one_cluster in levels(Idents(seurat_objects_integrated))) {
                        print(sprintf("cluster: %s", one_cluster))
                        cluster_conserved_markers <- FindConservedMarkers(seurat_objects_integrated, ident.1=one_cluster, grouping.var="condition", verbose=FALSE)
                        write.csv(cluster_conserved_markers, file=sprintf("%s/conserved_markers_cluster_%s.csv", conserved_markers_output_dir, one_cluster), quote=FALSE, row.names=TRUE)
                        print(FeaturePlot(seurat_objects_integrated, features=rownames(cluster_conserved_markers)[1:9], min.cutoff="q9") + plot_annotation(title=sprintf("Cluster %s", one_cluster), theme=theme(plot.title=element_text(size=16))))
                        print(DotPlot(seurat_objects_integrated, features=rownames(cluster_conserved_markers)[1:30], cols=c("blue", "red"), dot.scale=8, split.by="condition") + RotatedAxis() + plot_annotation(title=sprintf("Cluster %s", one_cluster), theme=theme(plot.title=element_text(size=16))))
                    }
                    dev.off()
                }
            }
            df2ggplot <- data.frame(integration=rep(integration_method, length(clustering_resolution_vector)), normalization=rep("SCTransform", length(clustering_resolution_vector)), features_param=rep(features_param, length(clustering_resolution_vector)), dimensions=rep(dimensions_nb_param, length(clustering_resolution_vector)), resolution=clustering_resolution_vector, silhouette_clusters=cluster_mean_silhouette_vector)
            integration_cluster_mean_silhouette_df <- rbind(integration_cluster_mean_silhouette_df, df2ggplot)
        }
        
        if (integration_method == "Harmony") {
            rm(seurat_objects_merge)
            gc()
        }
        
        rm(seurat_objects_integrated)
        gc()
    }
}
write.csv(integration_cluster_mean_silhouette_df, file=sprintf("%s/00-SCTransform/cluster_mean_silhouette.csv", output_dir), quote=FALSE, row.names=FALSE)












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




