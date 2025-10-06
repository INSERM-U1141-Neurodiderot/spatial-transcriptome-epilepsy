library(AUCell)
library(Seurat)
library(msigdbr)
library(ExperimentHub)
library(GSEABase)
library(dplyr)
library(SCENIC)

load("/work/project/fragencode/workspace/umr_1141/results/epilepsie/spatial/seurat_obj/whole_clustering_subclustering_Seurat_objects.RData")

epi_assay <- LayerData(seurat_object_whole_clustering, assay = "RNA", layer = "counts")

go_df <- msigdbr(species = "Rattus norvegicus", collection = "C5")
kegg_df <- msigdbr(species = "Rattus norvegicus", collection = "C2", subcollection = "CP:KEGG_LEGACY")

convert_to_GeneSetCollection <- function(msig_df) {
  # Split by gene set and ensure unique genes in each
  gene_sets <- msig_df %>%
    select(gs_name, gene_symbol) %>%
    distinct() %>%
    group_by(gs_name) %>%
    summarise(genes = list(unique(gene_symbol)), .groups = "drop")
  
  # Create GeneSet list
  gene_set_list <- lapply(1:nrow(gene_sets), function(i) {
    GeneSet(
      geneIds = gene_sets$genes[[i]],
      setName = gene_sets$gs_name[[i]],
      collectionType = ComputedCollection()
    )
  })

  GeneSetCollection(gene_set_list)
}

go_gsc <- convert_to_GeneSetCollection(go_df)
kegg_gsc <- convert_to_GeneSetCollection(kegg_df)

# Combine them into one list of GeneSet objects
all_gene_sets <- c(as.list(go_gsc), as.list(kegg_gsc))

# Create a single GeneSetCollection
merged_gsc <- GeneSetCollection(all_gene_sets)

cells_rankings <- AUCell_buildRankings(epi_assay, plotStats=FALSE)

cells_AUC <- AUCell_calcAUC(merged_gsc, cells_rankings, nCores=15)

save(cells_AUC, file="/home/adufour/work/rds_storage/epilepsie/cells_AUC_epilepsie.RData")