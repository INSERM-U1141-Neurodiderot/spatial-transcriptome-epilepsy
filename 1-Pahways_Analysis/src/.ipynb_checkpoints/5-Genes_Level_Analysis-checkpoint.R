.libPaths(c("usr/local/lib/R/site-library", .libPaths()))
library(Seurat)
library(dplyr)
library(tidyverse)
library(msigdbr)
library(ggplot2)
library(rlang)
library(future)
library(pbapply)
library(Matrix)
library(SeuratData)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(AnnotationDbi)
library(tibble)
library(clusterProfiler)

##############
# Parameters #
##############
work_dir <- "/home/ronan.jouanard/NeuroDev_ADD/Share_ronan.jouanard/EpiReg_suite"
old_work_dir <- "/home/ronan.jouanard/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"

input_model_organism <- "rat"
output_model_organism <- "human"

# For human
output_dir_human <- sprintf("/home/ronan.jouanard/NeuroDev_ADD/Share_ronan.jouanard/EpiReg_suite/1-Pathways_Analysis/output/%s", output_model_organism)

# For rat
output_dir_rat <- sprintf("/home/ronan.jouanard/NeuroDev_ADD/Share_ronan.jouanard/EpiReg_suite/1-Pathways_Analysis/output/%s", input_model_organism)

# Parameters
seurat_custom_cluster.v <- c("0.0","0.1","0.2","0.3","1","2","3.0","3.1","3.2","3.3","4")
time_points.v <- c("5","10","20","40")
pathway_collections <- c("HALLMARK", "REACTOME", "GOBP", "KEGG")
col_clusters <- c("1" = "#cccccc", "2" = "#f768a1", "4" = "#bae1ff", "0.0" = "#238b45", "0.1" = "#74c476", "0.2" = "#bae4b3", "0.3" = "#edf9e9", "3.0" = "#fdf498", "3.1" = "#807dba", "3.2" = "#4a1486", "3.3" = "#ac8fff")
time_points_col <- c("5" = "#feedde", "10" = "#fdbe85", "20" = "#fd8d3c", "40" = "#d94701")

# Gene expression Levels - Z-score design parameters
breaks_gradient <- c(-4, 0, 4)
colors_gradient <- c("blue", "white", "red")
color_gradient_mapping <- colorRamp2(breaks_gradient, colors_gradient)

##############################
# Functions for the analysis #
##############################

PrepSCTFindMarkers.test <- function(object, assay = "SCT", verbose = TRUE) {
  if (verbose && nbrOfWorkers() == 1) {
    my.lapply <- pblapply
  } else {
    my.lapply <- future_lapply
  }
  if (length(x = levels(x = object[[assay]])) == 1) {
    if (verbose) {
      message("Only one SCT model is stored - skipping recalculating corrected counts")
    }
    return(object)
  }
  observed_median_umis <- lapply(
    X = SCTResults(object = object[[assay]], slot = "cell.attributes"),
    FUN = function(x) median(x[, "umi"])
  )
  model.list <- slot(object = object[[assay]], name = "SCTModel.list")
  median_umi.status <- lapply(X = model.list,
                              FUN = function(x) { return(tryCatch(
                                expr = slot(object = x, name = 'median_umi'),
                                error = function(...) {return(NULL)})
                              )})
  if (any(is.null(x = unlist(x = median_umi.status)))){
    # For old SCT objects  median_umi is set to median umi as calculated from observed UMIs
    slot(object = object[[assay]], name = "SCTModel.list") <- lapply(X = model.list,
                                                                     FUN = UpdateSlots)
    SCTResults(object = object[[assay]], slot = "median_umi") <- observed_median_umis

  }
  model_median_umis <- SCTResults(object = object[[assay]], slot = "median_umi")
  min_median_umi <- min(unlist(x = observed_median_umis), na.rm = TRUE)
  if (all(unlist(x = model_median_umis) > min_median_umi)){
    if (verbose){
      message("Minimum UMI unchanged. Skipping re-correction.")
    }
    return(object)
  }
  if (verbose) {
    message(paste0("Found ",
                   length(x = levels(x = object[[assay]])),
                   " SCT models.",
                   " Recorrecting SCT counts using minimum median counts: ",
                   min_median_umi))
  }
  umi.assay <- unique(
    x = unlist(
      x = SCTResults(object = object[[assay]], slot = "umi.assay")
    )
  )
  if (length(x = umi.assay) > 1) {
    stop("Multiple UMI assays are used for SCTransform: ",
         paste(umi.assay, collapse = ", ")
    )
  }
  umi.layers <- Layers(object = object, assay = umi.assay, search = 'counts')
  if (length(x = umi.layers) > 1) {
    object[[umi.assay]] <- JoinLayers(
      object = object[[umi.assay]],
      layers = "counts", new = "counts")
  }
  raw_umi <- GetAssayData(object = object, assay = umi.assay, slot = "counts")
  corrected_counts <- Matrix(
    nrow = nrow(x = raw_umi),
    ncol = ncol(x = raw_umi),
    data = 0,
    dimnames = dimnames(x = raw_umi),
    sparse = TRUE
  )
  cell_attr <- SCTResults(object = object[[assay]], slot = "cell.attributes")
  model_pars_fit <- lapply(
    X = SCTResults(object = object[[assay]], slot = "feature.attributes"),
    FUN = function(x) x[, c("theta", "(Intercept)", "log_umi")]
  )
  arguments <- SCTResults(object = object[[assay]], slot = "arguments")
  model_str <- SCTResults(object = object[[assay]], slot = "model")
  set_median_umi <- rep(min_median_umi, length(levels(x = object[[assay]])))
  names(set_median_umi) <- levels(x = object[[assay]])
  set_median_umi <- as.list(set_median_umi)
  all_genes <- rownames(x = object[[assay]])
  # correct counts
  my.correct_counts <- function(model_name){
    model_genes <- rownames(x = model_pars_fit[[model_name]])
    x <- list(
      model_str = model_str[[model_name]],
      arguments = arguments[[model_name]],
      model_pars_fit = as.matrix(x = model_pars_fit[[model_name]]),
      cell_attr = cell_attr[[model_name]]
    )
    cells <- rownames(x = cell_attr[[model_name]])
    
    # Subset raw_umi only on common genes and cells
    common_genes <- intersect(all_genes, rownames(raw_umi))
    umi <- raw_umi[common_genes, cells]

    umi_corrected <- correct_counts(
      x = x,
      umi = umi,
      verbosity = 0,
      scale_factor = min_median_umi
    )
    missing_features <- setdiff(x = all_genes, y = rownames(x = umi_corrected))
    corrected_counts.list <- NULL
    gc(verbose = FALSE)
    empty <- SparseEmptyMatrix(nrow = length(x = missing_features), ncol = ncol(x = umi_corrected))
    rownames(x = empty) <- missing_features
    colnames(x = umi_corrected) <- colnames(x = umi_corrected)

    umi_corrected <- rbind(umi_corrected, empty)[all_genes,]

    return(umi_corrected)
  }
  corrected_counts.list <- my.lapply(X = levels(x = object[[assay]]),
                                     FUN = my.correct_counts)
  names(x = corrected_counts.list) <- levels(x = object[[assay]])

  corrected_counts <- do.call(what = MergeSparseMatrices, args = corrected_counts.list)
  corrected_counts <- as.sparse(x = corrected_counts)
  corrected_data <- log1p(x = corrected_counts)
  suppressWarnings({object <- SetAssayData(object = object,
                                           assay = assay,
                                           slot = "counts",
                                           new.data = corrected_counts)})
  suppressWarnings({object <- SetAssayData(object = object,
                                           assay = assay,
                                           slot = "data",
                                           new.data = corrected_data)})
  SCTResults(object = object[[assay]], slot = "median_umi") <- set_median_umi
  return(object)
}

MergeSparseMatrices <- function(...) {

  colname.new <- character()
  rowname.new <- character()
  x <- vector()
  i <- numeric()
  j <- numeric()

  for (mat in list(...)) {
    colname.old <- colnames(x = mat)
    rowname.old <- rownames(x = mat)

    # does not check if there are overlapping cells
    colname.new <- union(x = colname.new, y = colname.old)
    rowname.new <- union(x = rowname.new, y = rowname.old)

    colindex.new <- match(x = colname.old, table = colname.new)
    rowindex.new <- match(x = rowname.old, table = rowname.new)

    ind <- summary(object = mat)
    # Expand the list of indices and x
    i <- c(i, rowindex.new[ind[,1]])
    j <- c(j, colindex.new[ind[,2]])
    x <- c(x, ind[,3])
  }

  merged.mat <- sparseMatrix(i=i,
                             j=j,
                             x=x,
                             dims=c(length(rowname.new), length(colname.new)),
                             dimnames=list(rowname.new, colname.new))
  return (merged.mat)
}


correct_counts <- function(x, umi, cell_attr = x$cell_attr, scale_factor = NA, verbosity = 2,
                           verbose = NULL, show_progress = NULL) {
  # Take care of deprecated arguments
  if (!is.null(verbose)) {
    warning("The 'verbose' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)", immediate. = TRUE, call. = FALSE)
    verbosity <- as.numeric(verbose)
  }
  if (!is.null(show_progress)) {
    warning("The 'show_progress' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)", immediate. = TRUE, call. = FALSE)
    if (show_progress) {
      verbosity <- 2
    } else {
      verbosity <- min(verbosity, 1)
    }
  }

  regressor_data_orig <- model.matrix(as.formula(gsub('^y', '', x$model_str)), cell_attr)
  # when correcting, set all latent variables to median values
  cell_attr[, x$arguments$latent_var] <- apply(cell_attr[, x$arguments$latent_var, drop=FALSE], 2, function(x) rep(median(x), length(x)))

  if (!is.na(scale_factor) && !is.numeric(scale_factor)){
    stop("`scale_factor` should be numeric")
  }
  if (!is.na(scale_factor)){
    if (verbosity>0){
      message(paste("Setting log_umi for correcting counts to", scale_factor))
    }
    cell_attr[, "log_umi"] <- log10(scale_factor)
  }
  regressor_data <- model.matrix(as.formula(gsub('^y', '', x$model_str)), cell_attr)

  genes <- rownames(umi)[rownames(umi) %in% rownames(x$model_pars_fit)]
  bin_size <- x$arguments$bin_size
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (verbosity > 0) {
    message('Computing corrected UMI count matrix')
  }
  if (verbosity > 1) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  #corrected_data <- matrix(NA_real_, length(genes), nrow(regressor_data), dimnames = list(genes, rownames(regressor_data)))
  corrected_data <- list()
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    coefs <- x$model_pars_fit[genes_bin, -1, drop=FALSE]
    theta <- x$model_pars_fit[genes_bin, 1]
    # get pearson residuals
    mu <- exp(tcrossprod(coefs, regressor_data_orig))
    variance <- mu + mu^2 / theta
    y <- as.matrix(umi[genes_bin, , drop=FALSE])
    pearson_residual <- (y - mu) / sqrt(variance)
    # generate output
    mu <- exp(tcrossprod(coefs, regressor_data))
    variance <- mu + mu^2 / theta
    y.res <- mu + pearson_residual * sqrt(variance)
    y.res <- round(y.res, 0)
    y.res[y.res < 0] <- 0
    corrected_data[[length(corrected_data) + 1]] <- make.sparse(mat = y.res)
    if (verbosity > 1) {
      setTxtProgressBar(pb, i)
    }
  }
  if (verbosity > 1) {
    close(pb)
  }
  corrected_data <- do.call(what = rbind, args = corrected_data)

  return(corrected_data)
}

reverse_regression <- function(pearson_residual, theta, coefs, data) {
  mu <- exp(data %*% coefs)[, 1]
  variance <- mu + mu^2 / theta
  return(mu + pearson_residual * sqrt(variance))
}

#' @return A dgCMatrix
make.sparse <- function(mat){
  mat <- as(object = mat, Class = "Matrix")
  return (as(object = as(object = as(object = mat, Class = "dMatrix"), Class = "generalMatrix"), Class = "CsparseMatrix"))
}

Heatmap_structure <- function(data_scaled, separate_conditions = FALSE, cluster_hierarchisation = FALSE, time_ordered = FALSE) {
    # Removing X prefix in colnames
    colnames(data_scaled) <- gsub("X", "", colnames(data_scaled))
    
    # Clusters, condition, and time levels extraction
    condition <- colnames(data_scaled) %>% strsplit("_") %>% sapply(function(x) x[1])
    time <- colnames(data_scaled) %>% strsplit("_") %>% sapply(function(x) x[2]) %>% as.numeric() %>% as.factor()
    cluster <- colnames(data_scaled) %>% strsplit("_") %>% sapply(function(x) x[3])
    
    names(condition) <- colnames(data_scaled)
    names(time) <- colnames(data_scaled)
    names(cluster) <- colnames(data_scaled)
    
    # Create a data frame for reordering
    df <- data.frame(condition, time, cluster, combined_names = colnames(data_scaled), stringsAsFactors = FALSE)
    
    # Separate SE and CTRL samples if required
    if (separate_conditions) {
        se_df <- df[df$condition == "SE", ]
        ctrl_df <- df[df$condition == "CTRL", ]
        
        # Reorder based on specified conditions
        if ((cluster_hierarchisation == FALSE) & (time_ordered == FALSE)) {
            se_ordered <- se_df
            ctrl_ordered <- ctrl_df
        } else if ((cluster_hierarchisation == TRUE) & (time_ordered == FALSE)) {
            se_ordered <- se_df[order(se_df$cluster), ]
            ctrl_ordered <- ctrl_df[order(ctrl_df$cluster), ]
        } else if ((cluster_hierarchisation == TRUE) & (time_ordered == TRUE)) {
            se_ordered <- se_df[order(se_df$cluster, se_df$time), ]
            ctrl_ordered <- ctrl_df[order(ctrl_df$cluster, ctrl_df$time), ]
        } else {
            se_ordered <- se_df[order(se_df$time), ]
            ctrl_ordered <- ctrl_df[order(ctrl_df$time), ]
        }
        
        # Combine the ordered data frames with a white space in between
        df_ordered <- rbind(ctrl_ordered, df[0,], se_ordered)
    } else {
        # Reorder based on specified conditions
        if ((cluster_hierarchisation == FALSE) & (time_ordered == FALSE)) {
            df_ordered <- df[order(df$condition), ]
        } else if ((cluster_hierarchisation == TRUE) & (time_ordered == FALSE)) {
            df_ordered <- df[order(df$condition, df$cluster), ]
        } else if ((cluster_hierarchisation == TRUE) & (time_ordered == TRUE)) {
            df_ordered <- df[order(df$condition, df$cluster, df$time), ]
        } else if ((cluster_hierarchisation == FALSE) & (time_ordered == TRUE)) {
            # Custom ordering: CTRL and SE side by side block by block by time
            df_ordered <- df[order(df$time, df$condition), ]
        } else {
            df_ordered <- df[order(df$condition, df$time), ]
        }
    }
    
    # Extract ordered columns
    ordered_columns <- df_ordered$combined_names
    
    # Reorder data_scaled based on the ordered columns
    data_scaled_order <- data_scaled[, ordered_columns, drop = FALSE]
    
    # Get clusters and time orders
    clusters_order <- df_ordered$cluster
    times_order <- df_ordered$time
    
    # Get back to character vectors
    clusters_order <- as.character(clusters_order)
    times_order <- as.character(times_order)
    
    return(list(clusters_order, times_order, data_scaled_order))
}

SpatialPathways2GenesHeatmap <- function(seurat_object, clusters, time_points, genes, cluster_colors, time_colors, output_file, separate_conditions = FALSE, cluster_hierarchisation = FALSE, time_ordered = FALSE, row_km = NULL, column_km = NULL, cluster_rows = TRUE) {
  
  # Extract the data matrix from the Seurat object
  data <- GetAssayData(seurat_object, assay = "SCT", slot = "data")
  
  # Subset data for the selected genes
  data <- data[rownames(data) %in% genes, ]
  
  # Remove any rows with NA/NaN/Inf values
  data <- data[apply(data, 1, function(x) all(is.finite(x))), ]
  
  # Create a data frame with metadata
  metadata <- seurat_object@meta.data
  
  # Subset the metadata for the selected clusters and time points
  if (!is.null(clusters)) {
    metadata <- metadata[metadata$seurat_custom_clusters %in% clusters, ]
  }
  
  if (!is.null(time_points)) {
    metadata <- metadata[metadata$time %in% time_points, ]
  }
  
  # Subset the data for the selected cells
  selected_cells <- rownames(metadata)
  data <- data[, selected_cells]
  
  # Combine data and metadata
  combined_data <- data.frame(t(data))
  combined_data$condition <- metadata$condition
  combined_data$time_point <- metadata$time
  combined_data$cluster <- metadata$seurat_custom_clusters
  
  # Calculate the z-score for each gene in each cluster
  data_long <- combined_data %>%
    pivot_longer(cols = -c(condition, time_point, cluster), names_to = "gene", values_to = "expression") %>%
    group_by(gene, condition, time_point, cluster) %>%
    summarize(expression = mean(expression), .groups = 'drop') %>%
    ungroup()
  
  data_wide <- data_long %>%
    pivot_wider(names_from = c(condition, time_point, cluster), values_from = expression) %>%
    column_to_rownames(var = "gene") %>%
    as.matrix()
  
  # Scale the data to z-scores
  data_scaled <- t(scale(t(data_wide)))
  
  # Remove any rows with NA/NaN/Inf values after scaling
  data_scaled <- data_scaled[apply(data_scaled, 1, function(x) all(is.finite(x))), ]
  
  # Adjust the values and colors as needed
  expression_colors <- colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
  
  # Organize heatmap structure
  heatmap_structure <- Heatmap_structure(data_scaled, separate_conditions, cluster_hierarchisation, time_ordered)
  clusters_order <- heatmap_structure[[1]]
  times_order <- heatmap_structure[[2]]
  data_scaled_order <- heatmap_structure[[3]]
  
  # Ensure all clusters and times have defined colors
  cluster_levels <- unique(clusters_order)
  time_levels <- unique(times_order)
  
  cluster_colors <- cluster_colors[names(cluster_colors) %in% cluster_levels]
  time_colors <- time_colors[names(time_colors) %in% time_levels]
  
  # Create column annotations for clusters and time points
  column_annotation <- HeatmapAnnotation(
    Cluster = clusters_order,
    Time = times_order,
    col = list(Cluster = cluster_colors, Time = time_colors),
    gp = gpar(col = "white", lwd = 0.05),
    annotation_name_gp = gpar(fontsize = 9),
    simple_anno_size = unit(3, "mm")
  )
  
  # Generate the heatmap
  pdf(output_file, width = 12, height = 8)
  hm <- Heatmap(data_scaled_order,
                name = "Expression",
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 8),
                top_annotation = column_annotation,
                border = TRUE,
                show_row_dend = FALSE,
                show_column_dend = TRUE,
                show_column_names = TRUE,
                cluster_columns = FALSE,
                column_labels = colnames(data_scaled_order),
                column_dend_height = unit(3, "mm"),
                row_km = row_km,
                column_km = column_km,
                col = expression_colors) # Apply the custom color gradient
  draw(hm)
  
  dev.off()
}

## Genes annotation - Human
# Load Seurat object
seurat_orth_human_file <- file.path(output_dir_human, "seurat_obj/seurat_orth_human.RData")
load(seurat_orth_human_file)

### Initiate new cluster_time idents
seurat_orth_human$cluster.time.condition <- paste(seurat_orth_human$seurat_custom_clusters, seurat_orth_human$time, seurat_orth_human$condition,
    sep = "_")
Idents(seurat_orth_human) <- "cluster.time.condition"

### Find marker genes between SE and CTRL condition in identical cluster, identical time points 
seurat_orth_human <- PrepSCTFindMarkers.test(seurat_orth_human)

list_markers_time_point.cluster_SEvsCTRL <- list() 
for (time in time_points.v) {
    for (cluster in seurat_custom_cluster.v) {
        name_markers <- paste(cluster, time, sep = "_")
        ident_SE <- paste(cluster, time, "SE", sep = "_")
        ident_CTRL <- paste(cluster, time, "CTRL", sep = "_")
        list_markers_time_point.cluster_SEvsCTRL[[name_markers]] <- FindMarkers(seurat_orth_human, assay = "SCT", ident.1 = ident_SE, ident.2 = ident_CTRL, verbose = FALSE)
        }
    }
# Rank marker genes by abs(pct.1 - pct.2) * avg_log2FC * (-log10(p_val_adj))

filtered_markers_list <- list()
# For each time point and each cluster, rank markers by the score and filter
for (name_markers in names(list_markers_time_point.cluster_SEvsCTRL)) {
  # Extract the current markers data frame
  markers_df <- list_markers_time_point.cluster_SEvsCTRL[[name_markers]]
  
  # Apply the transformation to calculate the new metric and add it as a new column
  markers_df$rank_score <- abs(markers_df$pct.1 - markers_df$pct.2) * markers_df$avg_log2FC * (-log10(markers_df$p_val_adj))
  
  # Rank the markers from highest to lowest based on the new metric
  markers_df <- markers_df[order(-markers_df$rank_score), ]
  
  # Update the original list with the new ranked data frame
  list_markers_time_point.cluster_SEvsCTRL[[name_markers]] <- markers_df
  
  # Keep only genes with a rank score higher than 10 and store them in the new list
  filtered_markers_list[[name_markers]] <- markers_df[markers_df$rank_score > 10, ]
}

### Number of marker genes ound per cluster and time point                
# cluster size
cluster_sizes <- table(seurat_orth_human@meta.data$cluster.time.condition) %>% as.list()

# Compute mean cluster sizes between SE and CTRL for each cluster and time point
mean_cluster_sizes <- sapply(unique(gsub("_CTRL|_SE", "", names(cluster_sizes))), function(cluster) {
  unlist(cluster_sizes[paste0(cluster, c("_CTRL", "_SE"))]) %>% mean()
})

# Summary marker genes count 
df_summary_marker_genes <- tibble::tibble(
  cluster = names(filtered_markers_list),
  number_of_genes = sapply(filtered_markers_list, nrow),
  normalized_gene_count = sapply(names(filtered_markers_list), function(name) {
    nrow(filtered_markers_list[[name]]) / mean_cluster_sizes[name]
  })
)

df_summary_marker_genes <- df_summary_marker_genes %>%
  mutate(
    time_point = as.numeric(sub(".*_(\\d+)$", "\\1", cluster)),       # Extract time point
    cluster_number = sub("_.*$", "", cluster)                        # Extract everything before the first underscore
  ) %>%
  arrange(time_point, cluster_number) %>%  # Order by time point and cluster number
  mutate(cluster_ordered = factor(cluster, levels = cluster))

### Genes count

# Plot the bar plot of the number of genes per cluster with custom colors
plot <- ggplot(df_summary_marker_genes, aes(x = cluster_number, y = number_of_genes, fill = cluster_number)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.7) +
  scale_fill_manual(values = col_clusters) +  # Apply custom colors
  labs(
    title = "Number of Highly Differentially Expressed Genes (SE vs. CTRL) per Cluster by Time Point",
    x = "Cluster",
    y = "Genes Count"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ time_point, scales = "free_x")

pdf_file <- file.path(output_dir_human, "Gene_level_analysis/gene_count/gene_count_custom_cluster_by_time.pdf")
pdf(pdf_file, width = 10, height = 10)
print(plot)
dev.off()

### Normalised genes count
plot <- ggplot(df_summary_marker_genes, aes(x = cluster_number, y = normalized_gene_count, fill = cluster_number)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.7) +
  scale_fill_manual(values = col_clusters) +  # Apply custom colors
  labs(
    title = "Normalised (cluster size) number of Highly Differentially Expressed Genes (SE vs. CTRL) per Cluster by Time Point",
    x = "Cluster",
    y = "Genes Count"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ time_point, scales = "free_x")

pdf_file <- file.path(output_dir_human, "Gene_level_analysis/gene_count/normalised_gene_count_custom_cluster_by_time.pdf")
pdf(pdf_file, width = 10, height = 10)
print(plot)
dev.off()



### Markers frequency across clusters and time points
marker_counts <- bind_rows(lapply(names(filtered_markers_list), function(name) {
  filtered_markers_list[[name]] %>% 
    rownames_to_column("gene")
})) %>%
  count(gene, name = "frequency") %>%
  mutate(gene = factor(gene, levels = gene[order(-frequency)]))

# Plot the bar plot of marker frequencies with ordered markers
plot <- ggplot(marker_counts, aes(x = gene, y = frequency)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", alpha = 0.7) +
  labs(
    title = "Frequency of Marker Genes Across Clusters and Time Points",
    x = "Gene",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4))

pdf_file <- file.path(output_dir_human, "Gene_level_analysis/markers_frequency/markers_frequency_cluster_time_points.pdf")
pdf(pdf_file, width = 10, height = 10)
print(plot)
dev.off()

# Most frequent markers
frequent_marker_counts <- marker_counts %>% filter(frequency >= 5)

plot <- ggplot(frequent_marker_counts, aes(x = gene, y = frequency)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", alpha = 0.7) +
  labs(
    title = "Frequency of most frequent marker genes across clusters and time Points",
    x = "Gene",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf_file <- file.path(output_dir_human, "Gene_level_analysis/markers_frequency/top_markers_frequency_cluster_time_points.pdf")
pdf(pdf_file, width = 10, height = 10)
print(plot)
dev.off()

# Select top markers

top_markers.v <- frequent_marker_counts[, "gene"] 

### Expression level analysis - Top markers heatmap

pdf_file <- file.path(output_dir_human, "Gene_level_analysis/top_markers_heatmap/test.pdf")

SpatialPathways2GenesHeatmap(seurat_orth_human, clusters =  seurat_custom_cluster.v, time_points = time_points.v, genes = top_markers.v, cluster_colors = col_clusters, time_colors = time_points_col, output_file = pdf_file, separate_conditions = TRUE, cluster_hierarchisation = TRUE, time_ordered = TRUE)
