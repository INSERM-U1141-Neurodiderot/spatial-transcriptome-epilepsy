library(Seurat)
library(CARD)

reference <- readRDS("/work/project/fragencode/workspace/umr_1141/results/astro/big_object/l5_all_v2.rds")
load("/work/project/fragencode/workspace/umr_1141/results/epilepsie/spatial/seurat_obj/whole_clustering_subclustering_Seurat_objects.RData")

counts <- seurat_object_whole_clustering[["RNA"]]$counts

counts_ref <- reference[["RNA"]]$counts

metadata <- reference@meta.data

cell_names <- c()
x_coord <- c()
y_coord <- c()

for (i in names(seurat_object_whole_clustering@images)) {
	coords <- GetTissueCoordinates(seurat_object_whole_clustering, image = i)
	cell_names <- c(cell_names, rownames(coords))
	x_coord <- c(x_coord, coords$imagerow)
	y_coord <- c(y_coord, coords$imagecol)
}

coords <- data.frame(x = x_coord, y = y_coord)
rownames(coords) <- cell_names
coords <- coords[colnames(counts),]

CARD_obj = createCARDObject(
	sc_count = counts_ref,
	sc_meta = metadata,
	spatial_count = counts,
	spatial_location = coords,
	ct.varname = "Description",
	ct.select = unique(metadata$Description),
	sample.varname = "SampleID",
	minCountGene = 100,
	minCountSpot = 5)
	
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
	
saveRDS(CARD_obj, "/home/adufour/work/rds_storage/epilepsie/card_deconv.rds")