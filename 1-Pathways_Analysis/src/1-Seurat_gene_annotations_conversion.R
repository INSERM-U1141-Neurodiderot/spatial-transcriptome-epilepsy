.libPaths(c("usr/local/lib/R/site-library", .libPaths()))
library(Seurat)
library(msigdbr)
library(ggplot2)
library(grr)
library(ggtree)
library(RCurl)
library(gprofiler2)
library(orthogene)
library(SCPA)

##############
# Parameters #
##############
work_dir <- "/home/ronan.jouanard/NeuroDev_ADD/Share_ronan.jouanard/EpiReg_suite"
old_work_dir <- "/home/ronan.jouanard/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"

input_model_organism <- "rat"
output_model_organisms <- c("human", "mouse")

# For human gene annotations conversion
output_dir_human <- sprintf("/home/ronan.jouanard/NeuroDev_ADD/Share_ronan.jouanard/EpiReg_suite/1-Pathways_Analysis/output/%s", output_model_organisms[1])
if (! dir.exists(output_dir_human)) {
    dir.create(output_dir_human, recursive=TRUE, mode="0775")
}

# For mouse gene annotations conversion
output_dir_mouse <- sprintf("/home/ronan.jouanard/NeuroDev_ADD/Share_ronan.jouanard/EpiReg_suite/1-Pathways_Analysis/output/%s", output_model_organisms[2])
if (! dir.exists(output_dir_mouse)) {
    dir.create(output_dir_mouse, recursive=TRUE, mode="0775")
}

################################################
# Function - Seurat gene annotation conversion #
################################################

UpdateSeuratOrth <- function(seurat_obj, input_species, output_species) {
  
  # Function to map genes using orthologs conversion
  map_genes <- function(assay_data, input_species, output_species) {
    convert_orthologs(assay_data, input_species = input_species, output_species = output_species,
                      non121_strategy = "drop_both_species", drop_nonorths = TRUE, 
                      agg_fun = "max", method = "gprofiler")
  }
  
  # Loop through each assay in the Seurat object
  for (assay_name in names(seurat_obj@assays)) {
    cat("Processing assay:", assay_name, "\n")
    
    # Extract the assay data
    assay_data <- seurat_extract(seu_obj = seurat_obj, assay = assay_name)
    
    # Convert orthologs for this assay
    orthologs_data <- map_genes(assay_data, input_species, output_species)
    
    # Create a new Seurat object from the orthologs data
    options(Seurat.object.assay.version = "v3")
    seurat_new <- CreateSeuratObject(counts = orthologs_data, assay = assay_name)
    
    # Update the original Seurat object with the converted gene annotations
    seurat_obj@assays[[assay_name]]@counts <- seurat_new@assays[[assay_name]]@counts
    seurat_obj@assays[[assay_name]]@data <- seurat_new@assays[[assay_name]]@data
    seurat_obj@assays[[assay_name]]@scale.data <- seurat_new@assays[[assay_name]]@scale.data
    seurat_obj@assays[[assay_name]]@meta.features <- seurat_new@assays[[assay_name]]@meta.features
    
    # Modify SCT models if present
    if (assay_name == "SCT") {
      SCTModels <- seurat_obj[["SCT"]]@SCTModel.list
      for (i in seq_along(SCTModels)) {
        model <- SCTModels[[i]]
        features <- rownames(model@feature.attributes)

        cat(sprintf("Converting SCT model number %s\n", i))
          
        # Convert orthologs for SCT model features
        convert_orth_dict <- convert_orthologs(
          gene_df = features,
          input_species = input_species,
          output_species = output_species,
          non121_strategy = "drop_both_species",
          drop_nonorths = TRUE,
          agg_fun = "max",
          method = "gprofiler"
        )
        
        # Extract keys and values
        input_orth <- convert_orth_dict$input_gene
        output_orth <- rownames(convert_orth_dict)  # Convert each element to a string for display
        
        # Create the data frame
        df_orthologs <- data.frame(
          Input_Orthologs = input_orth,
          Output_Orthologs = output_orth,
          stringsAsFactors = FALSE
        )

        # Filter the model features to keep only genes present in the orthologs list
        model_filtered <- model@feature.attributes[rownames(model@feature.attributes) %in% df_orthologs$Input_Orthologs, , drop = FALSE]

        # Rename the row names to the output orthologs
        matched_indices <- match(rownames(model_filtered), df_orthologs$Input_Orthologs)
        new_rownames <- df_orthologs$Output_Orthologs[matched_indices]
        rownames(model_filtered) <- new_rownames
        
        # Ensure there are no underscores in the feature names
        rownames(model_filtered) <- gsub("_", "-", rownames(model_filtered))
        
        # Update the feature attributes in the model
        model@feature.attributes <- model_filtered
        
        # Store the modified model back in the Seurat object
        seurat_obj[["SCT"]]@SCTModel.list[[i]] <- model
      }
    }
  }
  
  return(seurat_obj)
}

# load data
## Seurat object
### Whole clustering

seurat_file <- sprintf("%s/20-Count_analysis/10-EpiReg/output/10-SR/TSO_polyA_R1trim1_ov5_n2_min20/RefSeq108/20-Integration/all/sub_markers/Norm_mito-scale_orig-int_orig/QC-sum1000_det500_mito100_hb100/10-QC_filtering/Seurat/whole_clustering_subclustering_Seurat_objects.RData", old_work_dir)

load(seurat_file)

# Annotations conversion to human 
seurat_orth_human <- UpdateSeuratOrth(seurat_obj = seurat_object_whole_clustering, input_species = "rat", output_species = "human")

# Converted features summary
rna_expm <- seurat_extract(seu_obj = seurat_object_whole_clustering, assay = "RNA")
conversion_genes <- orthogene::convert_orthologs(gene_df = rna_expm, gene_output = "columns", input_species = "rat", output_species = "human", non121_strategy = "drop_both_species")
colnames(conversion_genes)[2] <- "RNA_assay"
sct_counts_names <- seurat_orth_human[["SCT"]]@counts@Dimnames[[1]]
names(sct_counts_names) <- sct_counts_names

# Initialize the new column with NA values
conversion_genes$SCT_assay <- NA

# Match and fill the new column based on conversion_genes$RNA_assay
matches <- match(conversion_genes$RNA_assay, sct_counts_names)
conversion_genes$SCT_assay <- sct_counts_names[matches]

# Save converted features summary
output_dir_summary_human <- sprintf("%s/summary", output_dir_human)
if (! dir.exists(output_dir_summary_human)) {
    dir.create(output_dir_summary_human, recursive=TRUE, mode="0775")
}
conversion_summary_file <- file.path(output_dir_summary_human, sprintf("%s_features_conversion_summary.csv", output_model_organisms[1]))
write.csv(conversion_genes, conversion_summary_file, row.names = FALSE)

# Annotations conversion to mouse 
seurat_orth_mouse <- UpdateSeuratOrth(seurat_obj = seurat_object_whole_clustering, input_species = "rat", output_species = "mouse")

# Converted features summary
rna_expm <- seurat_extract(seu_obj = seurat_object_whole_clustering, assay = "RNA")
conversion_genes <- orthogene::convert_orthologs(gene_df = rna_expm, gene_output = "columns", input_species = "rat", output_species = "mouse", non121_strategy = "drop_both_species")
colnames(conversion_genes)[2] <- "RNA_assay"
sct_counts_names <- seurat_orth_mouse[["SCT"]]@counts@Dimnames[[1]]
names(sct_counts_names) <- sct_counts_names

# Initialize the new column with NA values
conversion_genes$SCT_assay <- NA

# Match and fill the new column based on conversion_genes$RNA_assay
matches <- match(conversion_genes$RNA_assay, sct_counts_names)
conversion_genes$SCT_assay <- sct_counts_names[matches]

# Save converted features summary
output_dir_summary_mouse <- sprintf("%s/summary", output_dir_mouse)
if (! dir.exists(output_dir_summary_mouse)) {
    dir.create(output_dir_summary_mouse, recursive=TRUE, mode="0775")
}
conversion_summary_file <- file.path(output_dir_summary_mouse, sprintf("%s_features_conversion_summary.csv", output_model_organisms[2]))
write.csv(conversion_genes, conversion_summary_file, row.names = FALSE)

# Save seurat objects
## Human
output_dir_seurat_human <- sprintf("%s/seurat_obj", output_dir_human)
if (! dir.exists(output_dir_seurat_human)) {
    dir.create(output_dir_seurat_human, recursive=TRUE, mode="0775")
}

save(seurat_orth_human, file = file.path(output_dir_seurat_human, "seurat_orth_human.RData"))

## Mouse
output_dir_seurat_mouse <- sprintf("%s/seurat_obj", output_dir_mouse)
if (! dir.exists(output_dir_seurat_mouse)) {
    dir.create(output_dir_seurat_mouse, recursive=TRUE, mode="0775")
}

save(seurat_orth_mouse, file = file.path(output_dir_seurat_mouse, "seurat_orth_mouse.RData"))
