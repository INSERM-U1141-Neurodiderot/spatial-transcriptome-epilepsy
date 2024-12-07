.libPaths(c("usr/local/lib/R/site-library", .libPaths()))
library(magrittr)
library(SCPA)
library(Seurat)
library(SingleCellExperiment)
library(SpatialExperiment)
library(ggspavis)
library(msigdbr)
library(withr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(backports)
library(ggpubr)
library(memoise, attach.required=T)
library(ggtree, attach.required=T)
library(gprofiler2)
library(grr)
library(orthogene)
library(RCurl)
library(ComplexHeatmap)
library(tzdb)
library(readr)
library(forcats)
library(lubridate)
library(tidyverse)
library(multicross)
library(circlize)

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

# Spatial visualisation - Q-value design parameters
breaks_gradient <- c(0, 5, 10)
colors_gradient <- c("white", "#85c6ad", "#304d36")
color_gradient_mapping <- colorRamp2(breaks_gradient, colors_gradient)

# SCPA outputs - directory
scpa_human <- sprintf("%s/SCPA_analysis", output_dir_human)

# Spatial visualisation of Q-value - directory
spatial_SCPA_vis_human <- sprintf("%s/SCPA_analysis/visualisation/spatial", output_dir_human)

# Seurat
## Human Gene Annotations
human_seurat_file <- sprintf("%s/seurat_obj/seurat_orth_human.RData", output_dir_human)
load(human_seurat_file)

##############################
# Functions for the analysis #
##############################

get_qval <- function(clusterid, time, pathway_qval_per_condition) {
pattern_cluster <- paste0("X", clusterid)
pattern <- paste0(pattern_cluster, ".*", time)
selected_columns <- grepl(pattern, colnames(pathway_qval_per_condition))
is_single_numeric_value <- is.numeric(pathway_qval_per_condition[, selected_columns]) && length(pathway_qval_per_condition[, selected_columns]) == 1
if (is_single_numeric_value) {
    return(pathway_qval_per_condition[, selected_columns])
    }
    return(NA)}

spe2speGOqval <- function(spe_obj, GOqval_cluster_time)  {
    
    all_barcodes <- rownames(colData(spe_obj))
    
    go_term_qval_df <- data.frame()
    unique_times <- unique(GOqval_cluster_time$time)
    
    complete_pathway_name <- unique(GOqval_cluster_time$complete_pathway_name)
    
    for (one_time in unique_times) {
        one_time_qval_df <- GOqval_cluster_time[GOqval_cluster_time$time == one_time, ]
        unique_clusters <- unique(one_time_qval_df$cluster)
        
        for (one_cluster in unique_clusters) {
            one_cluster_one_time_qval_df <- one_time_qval_df[one_time_qval_df$cluster == one_cluster, ]
            
            if (nrow(one_cluster_one_time_qval_df) == 2) {
                min_qvalue <- min(one_cluster_one_time_qval_df$qval)
                one_qvalue <- one_cluster_one_time_qval_df$qval == min_qvalue
                one_category <- one_cluster_one_time_qval_df$category[one_qvalue]
                one_qvalue <- one_cluster_one_time_qval_df$qval[one_qvalue]
            } else if (nrow(one_cluster_one_time_qval_df) == 1) {
                one_qvalue <- one_cluster_one_time_qval_df$qval
                one_category <- one_cluster_one_time_qval_df$category
            }
            
            barcodes <- rownames(colData(spe_obj)[which(colData(spe_obj)$time == one_time & colData(spe_obj)$seurat_custom_clusters == one_cluster),])
            go_term_qval_df <- rbind(go_term_qval_df, data.frame(barcode=barcodes, GO_qvalue=rep(one_qvalue, length(barcodes)), category=rep(one_category, length(barcodes))))
        }
    }
        
    na_barcodes <- all_barcodes[! all_barcodes %in% go_term_qval_df$barcode]
    na_barcodes_nb <- length(na_barcodes)
    go_term_qval_df <- rbind(go_term_qval_df, data.frame(barcode=na_barcodes, GO_qvalue=rep(NA, na_barcodes_nb), category=rep(NA, na_barcodes_nb)))
    
    rownames(go_term_qval_df) <- go_term_qval_df$barcode
    go_term_qval_df <- go_term_qval_df[rownames(colData(spe_obj)),]
    colData(spe_obj)$GO_qvalue <- go_term_qval_df$GO_qvalue
    
    # Modification of barcodes prefix
    spe_obj@colData@rownames <- sprintf("%s_%s", complete_pathway_name, spe_obj@colData@rownames)
    spe_obj@colData$sample_id <- sprintf("%s_%s", complete_pathway_name, spe_obj@colData$sample_id)
    spe_obj@int_metadata$imgData$sample_id <- sprintf("%s_%s", complete_pathway_name, spe_obj@int_metadata$imgData$sample_id)
    
    return(spe_obj)
}

get_GOqval_cluster_time <- function(save_dir, complete_pathway_name) {
    
    parts_complete_pathway_name <- strsplit(x = complete_pathway_name, split = "_")[[1]]
    
    pathway_collection <- parts_complete_pathway_name[1]
    
    collection_pathway_name <- parts_complete_pathway_name[-1] %>% paste(collapse = "_")
    
    scpa_output_file <- sprintf("%s/%s.csv", save_dir, pathway_collection)
    scpa_output <- read.csv(file = scpa_output_file)
    scpa_output <- scpa_output %>%  column_to_rownames(var = "X")
    
    clusterids <- strsplit(colnames(scpa_output), split = "_")  %>% sapply("[[", 1) %>% sub("X", "", .)
    timepoints <- strsplit(colnames(scpa_output), split = "_")  %>% sapply("[[", 3) %>% sub("X", "", .)

    GOqval_cluster_time <- tibble(cluster = clusterids, time = timepoints)
    GOqval_cluster_time <- GOqval_cluster_time %>% mutate(row = row_number(), .before = cluster)
    
    pathway_qval_given_condition <- scpa_output[which(rownames(scpa_output) == complete_pathway_name),]
    
    # Add the q-values as a new column using map2
    GOqval_cluster_time <- GOqval_cluster_time %>%
    mutate(qval = NA)
    
    for (row_number in GOqval_cluster_time$row) {
    cluster <- GOqval_cluster_time$cluster[row_number]
    time <- GOqval_cluster_time$time[row_number]
    GOqval_cluster_time$qval[row_number] <-  get_qval(clusterid = cluster, time = time, pathway_qval_per_condition = pathway_qval_given_condition)
    }
    
    # Add the q-values as a new column using map2
    GOqval_cluster_time <- GOqval_cluster_time %>%
    mutate(category = NA, .before = "qval")
    
    for (row_number in GOqval_cluster_time$row) {
    if (sign(GOqval_cluster_time$qval[row_number]) == 1) {
        GOqval_cluster_time$category[row_number] <- "up"
        }
    else if (sign(GOqval_cluster_time$qval[row_number]) == -1) {
        GOqval_cluster_time$category[row_number] <- "down"
        }
    else {
        GOqval_cluster_time$category[row_number] <- "0"
    }}
    
    GOqval_cluster_time <- GOqval_cluster_time %>% mutate(complete_pathway_name = complete_pathway_name, .before = "cluster")

return(GOqval_cluster_time)}

plotVisiumQval <- function(spe_obj, scpa_save_dir, sample_ids, complete_pathway_name, name_replicate, spat_vis_dir) {
    
    spe_sub <- spe_obj[, colData(spe_obj)$sample_id %in% sample_ids]
    
    GOqval_cluster_time <- get_GOqval_cluster_time(save_dir = scpa_save_dir, complete_pathway_name = complete_pathway_name)
    
    spe_sub_GOqval <- spe2speGOqval(spe_obj = spe_sub, GOqval_cluster_time = GOqval_cluster_time)
    
    # Rename "GO_qvalue" to "Qval"
    names(colData(spe_sub_GOqval))[names(colData(spe_sub_GOqval)) == "GO_qvalue"] <- "Qval"
    
    min_qval <- min(colData(spe_sub_GOqval)$Qval, na.rm = TRUE)
    max_qval <- max(colData(spe_sub_GOqval)$Qval, na.rm = TRUE)
    
    ### PDF
    if (!dir.exists(spat_vis_dir)) {
        dir.create(spat_vis_dir, recursive = TRUE)
    }
    
    pdf_file_name <- sprintf("%s_%s.pdf", complete_pathway_name, name_replicate)
    pdf_full_path <- file.path(spat_vis_dir, pdf_file_name)
    
    for (sample_id in sample_ids) { 
    
    indices <- grepl(sample_id, colData(spe_sub_GOqval)$sample_id) %>% which()
    
    colData(spe_sub_GOqval)$sample_id[indices] <- sample_id
    }
    
    pdf(pdf_full_path)

    # Plot using the new column name "Qval"
    p <- plotVisium(spe_sub_GOqval, annotate="Qval", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    
    p <- p + scale_fill_gradient2(midpoint = breaks_gradient[2], limit = c(breaks_gradient[1],breaks_gradient[3]), low = colors_gradient[1], mid = colors_gradient[2], high = colors_gradient[3]) 
    
    p <- p + ggtitle(complete_pathway_name)
    
    p <- p + theme(plot.title = element_text(hjust = 0.5, vjust = 5, size = 8), 
                   legend.position="bottom", 
                   legend.key.size=unit(1, "lines"), 
                   legend.title=element_text(size=8), 
                   legend.text=element_text(size=7))

    p$facet$params$ncol <- 4
    p$layers[[length(p$layers)]]$aes_params$size <- 0.5
    
    print(p)
    dev.off()
}

# Spatial Experiment object generation from Seurat object
## Human
sce_human <- as.SingleCellExperiment(seurat_orth_human)
sample_ids_human <- names(seurat_orth_human@images)

imgData_human <- do.call("rbind", lapply(sample_ids_human, function(x, seurat_obj=seurat_orth_human) {
    img <- SpatialImage(x=as.raster(seurat_obj@images[[x]]@image))
    imgDFrame <- DataFrame(sample_id=x, image_id=x, data=I(list(img)), scaleFactor=seurat_obj@images[[x]]@scale.factors$lowres)
    return(imgDFrame)
}))

spe_human <- SpatialExperiment(assays=assays(sce_human), rowData=rowData(sce_human), colData=colData(sce_human), metadata=metadata(sce_human), reducedDims=reducedDims(sce_human), altExps=altExps(sce_human), sample_id=as.character(colData(sce_human)$orig.ident), spatialCoordsNames=c("pxl_col_in_fullres", "pxl_row_in_fullres"), imgData=imgData_human)

## Replicates samples - CTRL
sample_ids_1st_repli_CTRL <- c("B_L1_S2", "B_L2_S6", "D_L1_S4", "D_L2_S8")
sample_ids_2nd_repli_CTRL <- c("B_L3_S10", "B_L4_S14", "D_L4_S16", "D_L3_S12")

## Replicates samples - SE
sample_ids_1st_repli_SE <- c("A_L1_S1", "A_L2_S5", "C_L1_S3", "C_L2_S7")
sample_ids_2nd_repli_SE <- c("A_L3_S9", "A_L4_S13", "C_L4_S15", "C_L3_S11")

## Custom grouping 
sample_ids_custom <- c("B_L1_S2", "B_L4_S14", "D_L1_S4", "D_L3_S12")

# Visualisation Q-value
## Human
### GOBP : REGULATION_OF_NEUROGENESIS
output_dir_spat_GOBP_REGULATION_OF_NEUROGENESIS <- sprintf("%s/GOBP_REGULATION_OF_NEUROGENESIS", spatial_SCPA_vis_human)
if (! dir.exists(output_dir_spat_GOBP_REGULATION_OF_NEUROGENESIS)) {
    dir.create(output_dir_spat_GOBP_REGULATION_OF_NEUROGENESIS, recursive=TRUE, mode="0775")
}
complete_pathway_name <- "GOBP_REGULATION_OF_NEUROGENESIS"

## Replicates samples - CTRL
### 1st replicates
plotVisiumQval(spe_obj = spe_human, scpa_save_dir = scpa_human, sample_ids = sample_ids_1st_repli_CTRL, complete_pathway_name = complete_pathway_name, name_replicate = "1st_replicates_CTRL", spat_vis_dir = output_dir_spat_GOBP_REGULATION_OF_NEUROGENESIS)
### 2nd replicates
plotVisiumQval(spe_obj = spe_human, scpa_save_dir = scpa_human, sample_ids = sample_ids_2nd_repli_CTRL, complete_pathway_name = complete_pathway_name, name_replicate = "2nd_replicates_CTRL", spat_vis_dir = output_dir_spat_GOBP_REGULATION_OF_NEUROGENESIS)
## Replicates samples - SE
### 1st replicates
plotVisiumQval(spe_obj = spe_human, scpa_save_dir = scpa_human, sample_ids = sample_ids_1st_repli_SE, complete_pathway_name = complete_pathway_name, name_replicate = "1st_replicates_SE", spat_vis_dir = output_dir_spat_GOBP_REGULATION_OF_NEUROGENESIS)
### 2nd replicates
plotVisiumQval(spe_obj = spe_human, scpa_save_dir = scpa_human, sample_ids = sample_ids_2nd_repli_SE, complete_pathway_name = complete_pathway_name, name_replicate = "2nd_replicates_SE", spat_vis_dir = output_dir_spat_GOBP_REGULATION_OF_NEUROGENESIS)
## Custom grouping
plotVisiumQval(spe_obj = spe_human, scpa_save_dir = scpa_human, sample_ids = sample_ids_custom, complete_pathway_name = complete_pathway_name, name_replicate = "custom", spat_vis_dir = output_dir_spat_GOBP_REGULATION_OF_NEUROGENESIS)

### GOBP : REGULATION_OF_TRANS_SYNAPTIC_SIGNALING
output_dir_spat_GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING <- sprintf("%s/GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING", spatial_SCPA_vis_human)
if (! dir.exists(output_dir_spat_GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING)) {
    dir.create(output_dir_spat_GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING, recursive=TRUE, mode="0775")
}
complete_pathway_name <- "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING"

## Replicates samples - CTRL
### 1st replicates
plotVisiumQval(spe_obj = spe_human, scpa_save_dir = scpa_human, sample_ids = sample_ids_1st_repli_CTRL, complete_pathway_name = complete_pathway_name, name_replicate = "1st_replicates_CTRL", spat_vis_dir = output_dir_spat_GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING)
### 2nd replicates
plotVisiumQval(spe_obj = spe_human, scpa_save_dir = scpa_human, sample_ids = sample_ids_2nd_repli_CTRL, complete_pathway_name = complete_pathway_name, name_replicate = "2nd_replicates_CTRL", spat_vis_dir = output_dir_spat_GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING)
## Replicates samples - SE
### 1st replicates
plotVisiumQval(spe_obj = spe_human, scpa_save_dir = scpa_human, sample_ids = sample_ids_1st_repli_SE, complete_pathway_name = complete_pathway_name, name_replicate = "1st_replicates_SE", spat_vis_dir = output_dir_spat_GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING)
### 2nd replicates
plotVisiumQval(spe_obj = spe_human, scpa_save_dir = scpa_human, sample_ids = sample_ids_2nd_repli_SE, complete_pathway_name = complete_pathway_name, name_replicate = "2nd_replicates_SE", spat_vis_dir = output_dir_spat_GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING)
## Custom grouping
plotVisiumQval(spe_obj = spe_human, scpa_save_dir = scpa_human, sample_ids = sample_ids_custom, complete_pathway_name = complete_pathway_name, name_replicate = "custom", spat_vis_dir = output_dir_spat_GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING)