library(SCPA)
library(dplyr)
library(tibble)
library(msigdbr)
library(Seurat)
library(magrittr)
library(purrr)

load("/work/project/fragencode/workspace/umr_1141/results/epilepsie/spatial/seurat_obj/whole_clustering_subclustering_Seurat_objects.RData")

tb1 <- msigdbr("Rattus norvegicus", category = "C2", subcategory = "CP:KEGG")
tb2 <- msigdbr("Rattus norvegicus", category = "C5", subcategory = "GO:BP")
tb3 <- msigdbr("Rattus norvegicus", category = "H")

pathways <- rbind(tb1, tb2, tb3) %>%
  format_pathways()


seurat_object_whole_clustering@meta.data$condition_time <- paste0(seurat_object_whole_clustering@meta.data$condition, "_", seurat_object_whole_clustering@meta.data$time)

epilepsie_list <- SplitObject(seurat_object_whole_clustering, split.by = "condition_time")

##
# Epilepsie vs Control D40
##

scpa_out_d40_CTRL_EPI <- list()
for (i in unique(seurat_object_whole_clustering@meta.data$seurat_custom_clusters)) {
  
  epilepticus <- seurat_extract(epilepsie_list$SE_40, 
                            meta1 = "seurat_custom_clusters", value_meta1 = i)
  
  control <- seurat_extract(epilepsie_list$CTRL_40, 
                          meta1 = "seurat_custom_clusters", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out_d40_CTRL_EPI[[i]] <- compare_pathways(list(epilepticus, control), pathways, parallel = TRUE, cores = 10) %>%
    select(Pathway, qval, adjPval, FC) %>%
    set_colnames(c("Pathway", paste(i, "qval", sep = "_"), paste(i, "adjPval", sep = "_"), paste(i, "FC", sep = "_")))
  
}

saveRDS(scpa_out_d40_CTRL_EPI, file = "/home/adufour/work/rds_storage/epilepsie/SCPA_epilepsie_d40_ctrl_epi.rds")

##
# Epilepsie vs Control D20
##

scpa_out_d20_CTRL_EPI <- list()
for (i in unique(seurat_object_whole_clustering@meta.data$seurat_custom_clusters)) {
  
  epilepticus <- seurat_extract(epilepsie_list$SE_20, 
                            meta1 = "seurat_custom_clusters", value_meta1 = i)

  control <- seurat_extract(epilepsie_list$CTRL_20,
                          meta1 = "seurat_custom_clusters", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out_d20_CTRL_EPI[[i]] <- compare_pathways(list(epilepticus, control), pathways, parallel = TRUE, cores = 10) %>%
    select(Pathway, qval, adjPval, FC) %>%
    set_colnames(c("Pathway", paste(i, "qval", sep = "_"), paste(i, "adjPval", sep = "_"), paste(i, "FC", sep = "_")))
  
}

saveRDS(scpa_out_d20_CTRL_EPI, file = "/home/adufour/work/rds_storage/epilepsie/SCPA_epilepsie_d20_ctrl_epi.rds")

##
# Epilepsie vs Control D10
##

scpa_out_d10_CTRL_EPI <- list()
for (i in unique(seurat_object_whole_clustering@meta.data$seurat_custom_clusters)) {
  
  epilepticus <- seurat_extract(epilepsie_list$SE_10, 
                            meta1 = "seurat_custom_clusters", value_meta1 = i)

  control <- seurat_extract(epilepsie_list$CTRL_10,
                          meta1 = "seurat_custom_clusters", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out_d10_CTRL_EPI[[i]] <- compare_pathways(list(epilepticus, control), pathways, parallel = TRUE, cores = 10) %>%
    select(Pathway, qval, adjPval, FC) %>%
    set_colnames(c("Pathway", paste(i, "qval", sep = "_"), paste(i, "adjPval", sep = "_"), paste(i, "FC", sep = "_")))
  
}

saveRDS(scpa_out_d10_CTRL_EPI, file = "/home/adufour/work/rds_storage/epilepsie/SCPA_epilepsie_d10_ctrl_epi.rds")

##
# Epilepsie vs Control D5
##

scpa_out_d5_CTRL_EPI <- list()
for (i in unique(seurat_object_whole_clustering@meta.data$seurat_custom_clusters)) {
  
  epilepticus <- seurat_extract(epilepsie_list$SE_5, 
                            meta1 = "seurat_custom_clusters", value_meta1 = i)

  control <- seurat_extract(epilepsie_list$CTRL_5,
                          meta1 = "seurat_custom_clusters", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out_d5_CTRL_EPI[[i]] <- compare_pathways(list(epilepticus, control), pathways, parallel = TRUE, cores = 10) %>%
    select(Pathway, qval, adjPval, FC) %>%
    set_colnames(c("Pathway", paste(i, "qval", sep = "_"), paste(i, "adjPval", sep = "_"), paste(i, "FC", sep = "_")))
  
}

saveRDS(scpa_out_d5_CTRL_EPI, file = "/home/adufour/work/rds_storage/epilepsie/SCPA_epilepsie_d5_ctrl_epi.rds")