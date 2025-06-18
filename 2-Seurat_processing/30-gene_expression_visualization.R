.libPaths(c("/home/christophe.lepriol/NeuroDev_ADD/R/r_4.1.0", .libPaths()))
library(SpatialExperiment)
library(ggspavis)
library(Seurat)
library(scater)
library(patchwork)



##############
# Parameters #
##############
# dataset <- "GSE136913"
dataset <- "GSE140420"
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics"
input_dir <- sprintf("%s/projects/20-R_norvegicus_ref_genome_eval/20-Expression_analysis/input/%s", work_dir, dataset)
output_dir <- sprintf("%s/projects/20-R_norvegicus_ref_genome_eval/20-Expression_analysis/output/20-DE_analysis/%s", work_dir, dataset)
src_dir <- sprintf("%s/projects/20-R_norvegicus_ref_genome_eval/20-Expression_analysis/src", work_dir)

# reference genomes
raw_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/data/ref_genome/Rattus_norvegicus"
ensembl_104 <- list(ref_genome="NCBIRefSeq106_Ensembl104GTF", gene_ID_column_name="ensembl_gene_id", AnnotationDb_keytype="ENSEMBL", gtf_file=sprintf("%s/Rnor_6.0/annotations/Ensembl/release_104/gtf/Rattus_norvegicus.Rnor_6.0.104.NCBIseqname.gtf.gz", raw_dir))
ncbi_refseq_106 <- list(ref_genome="NCBIRefSeq106_NCBIRefSeq106GTF", gene_ID_column_name="entrezgene_accession", AnnotationDb_keytype="SYMBOL", gtf_file=sprintf("%s/Rnor_6.0/annotations/NCBI_RefSeq/release_106/GCF_000001895.5_Rnor_6.0_genomic.gtf.gz", raw_dir))
ncbi_refseq_108 <- list(ref_genome="NCBIRefSeq108_NCBIRefSeq108GTF", gene_ID_column_name="entrezgene_accession", AnnotationDb_keytype="SYMBOL", gtf_file=sprintf("%s/mRatBN7.2/annotations/NCBI_RefSeq/release_108/GCF_015227675.2_mRatBN7.2_genomic.gtf.gz", raw_dir))
ensembl_105 <- list(ref_genome="NCBIRefSeq108_Ensembl105GTF", gene_ID_column_name="ensembl_gene_id", AnnotationDb_keytype="ENSEMBL", gtf_file=sprintf("%s/mRatBN7.2/annotations/Ensembl/release_105/gtf/Rattus_norvegicus.mRatBN7.2.105.NCBIseqname.gtf.gz", raw_dir))
# genome_list <- list(ensembl_104, ncbi_refseq_106, ncbi_refseq_108, ensembl_105)
genome_list <- list(ncbi_refseq_108, ensembl_105)
# genome_list <- list(ensembl_104)
# names(genome_list) <- c("Ensembl 104", "NCBI RefSeq 106", "NCBI RefSeq 108", "Ensembl 105")
names(genome_list) <- c("NCBI RefSeq 108", "Ensembl 105")
# names(genome_list) <- c("Ensembl 104")

# n value in consistency/stringency formula
consistency_stringency_n_value_to_add <- 1

# design and contrast names
## GSE140420
### Region, blocking: Age
design_formula <- "~Age+Region"
condition_name <- "Region"
condition_ref <- NA
second_condition_name <- NA
design_output_name_append <- NULL
### Region and Age: all combinations (edgeR user's guide 3.3.1 Experiments with all combinations of multiple factors/Defining each treatment combination as a group)
#### contrasts: Region
# condition_name <- "Region"
# condition_ref <- NA
# second_condition_name <- "Age"
# design_output_name_append <- sprintf("_%s", condition_name)

# design output name
design_output_name <- sub("~|~0\\+", "", design_formula)
design_output_name <- gsub("\\+", "_", design_output_name)
design_output_name <- gsub("\\.", "", design_output_name)

# contrast names
## GSE140420
### Region, blocking: Age
contrast_output_names <- c("CA1vsCA3DG", "CA3vsCA1DG", "DGvsCA1CA3")
### Region and Age: all combinations (edgeR user's guide 3.3.1 Experiments with all combinations of multiple factors/Defining each treatment combination as a group)
#### contrasts: Region
# contrast_output_names <- c("12months_CA3vsCA1", "12months_DGvsCA1", "12months_CA1vsCA3DG", "12months_DGvsCA3", "12months_CA3vsCA1DG", "12months_DGvsCA1CA3", "5months_CA3vsCA1", "5months_DGvsCA1", "5months_CA1vsCA3DG", "5months_DGvsCA3", "5months_CA3vsCA1DG", "5months_DGvsCA1CA3")

# fold-change threshold
# fc_thresholds <- c(1.25, 1.5)
fc_thresholds <- c(1.25)

if (is.null(design_output_name_append)) {
  design_output_name_append <- ""
}

samples <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")
# sample <- "A_L1_S1"


#############
# Functions #
#############

load_ST_pipeline_matrix_and_normalize <- function(matrix_file, space_ranger_dir, sample) {
    st_pipeline_matrix <- read.table(matrix_file, sep="\t", header=TRUE, quote="", row.names=1)
    
    # build SpatialExperiment object
    spatial_coordinates_file <- file.path(space_ranger_dir, "spatial", "tissue_positions_list.csv")
    spatial_coordinates_df <- read.csv(spatial_coordinates_file, header=FALSE, quote="")
    colnames(spatial_coordinates_df) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    rownames(spatial_coordinates_df) <- spatial_coordinates_df$barcode
    
    ### matrix: convert coordinates to barcodes
    #barcodes <- unlist(lapply(rownames(st_pipeline_matrix), function(x, coord2barcode=visium_v1_coordinates_df) {
    #    coordinates <- unlist(strsplit(x, "x"))
    #    return(coord2barcode[which(coord2barcode$x==as.integer(coordinates[1]) & coord2barcode$y==as.integer(coordinates[2])), "barcode"])
    #}))
    barcodes <- unlist(lapply(rownames(st_pipeline_matrix), function(x, coord2barcode=spatial_coordinates_df) {
            coordinates <- unlist(strsplit(x, "x"))
            return(coord2barcode[which(coord2barcode$array_col==as.integer(coordinates[1])-1 & coord2barcode$array_row==as.integer(coordinates[2])-1), "barcode"])
        }))
    rownames(st_pipeline_matrix) <- barcodes
    
    gene_data <- data.frame(gene_name=colnames(st_pipeline_matrix))
    
    ### image data
    img <- readImgData(imageSources=file.path(space_ranger_dir, "spatial", "tissue_lowres_image.png"), scaleFactors=file.path(space_ranger_output_dir, "spatial", "scalefactors_json.json"), sample_id=sample)
    
    spe <- SpatialExperiment(assay=list(counts=t(st_pipeline_matrix)), rowData=gene_data, colData=spatial_coordinates_df[rownames(st_pipeline_matrix),], imgData=img, spatialCoordsNames=c("pxl_col_in_fullres", "pxl_row_in_fullres"), sample_id=sample)
    
    # QC
    # keep spot over tissue
    spe_in_tissue <- spe[, colData(spe)$in_tissue==1]
    
    # normalization
    ## create Seurat object
    seurat_in_tissue <- CreateSeuratObject(counts=counts(spe_in_tissue))
    seurat_in_tissue <- SCTransform(seurat_in_tissue, assay = "RNA")
    ## conversion to SingleCellExperiment object
    norm_counts <- seurat_in_tissue@assays$SCT@counts
    ## '_' in row names was replaced by '-', example: 'X--ambiguous.L3mbtl3.LOC120097350.' in norm_counts and 'X__ambiguous.L3mbtl3.LOC120097350.' in assay(spe_in_tissue, "counts")
    rownames(norm_counts)[grep("-", rownames(norm_counts))] <- gsub("-", "_", rownames(norm_counts)[grep("-", rownames(norm_counts))])
    ## remove genes not in normalized counts
    spe_in_tissue_norm <- spe_in_tissue[rownames(norm_counts),]
    assay(spe_in_tissue_norm, i="SCT") <- as.matrix(norm_counts)
    return(spe_in_tissue_norm)
}




############
# Analysis #
############

work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"

for (one_sample in samples) {
    print(sprintf("sample: %s", one_sample))

    for(i in 1:(length(names(genome_list))-1)) {
        exp_matrix_1_name <- names(genome_list)[i]
        exp_matrix_1_ref_genome <- genome_list[[exp_matrix_1_name]]$ref_genome
        exp_matrix_1_ID_col <- genome_list[[exp_matrix_1_name]]$gene_ID_column_name

        # load ST Pipeline expression matrix
        exp_matrix_1_st_pipeline_matrix_file <- sprintf("%s/10-ST_analysis/00-ST_Pipeline/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s/10-Pipeline/%s_stdata.tsv", work_dir, exp_matrix_1_ref_genome, one_sample, one_sample)
        exp_matrix_1_st_pipeline_matrix <- read.table(exp_matrix_1_st_pipeline_matrix_file, sep="\t", header=TRUE, quote="", row.names=1)
        ## Space Ranger output directory for tissue positions list and image
        space_ranger_output_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs", work_dir, exp_matrix_1_ref_genome, one_sample)
        exp_matrix_1_spe_norm <- load_ST_pipeline_matrix_and_normalize(exp_matrix_1_st_pipeline_matrix_file, space_ranger_output_dir, one_sample)

        for(j in (i+1):length(names(genome_list))) {
            exp_matrix_2_name <- names(genome_list)[j]
            exp_matrix_2_ref_genome <- genome_list[[exp_matrix_2_name]]$ref_genome
            exp_matrix_2_ID_col <- genome_list[[exp_matrix_2_name]]$gene_ID_column_name

            # load ST Pipeline expression matrix
            exp_matrix_2_st_pipeline_matrix_file <- sprintf("%s/10-ST_analysis/00-ST_Pipeline/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s/10-Pipeline/%s_stdata.tsv", work_dir, exp_matrix_2_ref_genome, one_sample, one_sample)
            exp_matrix_2_st_pipeline_matrix <- read.table(exp_matrix_2_st_pipeline_matrix_file, sep="\t", header=TRUE, quote="", row.names=1)
            ## Space Ranger output directory for tissue positions list and image
            space_ranger_output_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs", work_dir, exp_matrix_2_ref_genome, one_sample)
            exp_matrix_2_spe_norm <- load_ST_pipeline_matrix_and_normalize(exp_matrix_2_st_pipeline_matrix_file, space_ranger_output_dir, one_sample)


            comparison <- sprintf("%s_%s", gsub(" ", "", exp_matrix_1_name), gsub(" ", "", exp_matrix_2_name))
            print(sprintf("comparison: %s", comparison))

            for (one_contrast_output_name in contrast_output_names) {
                print(sprintf("contrast: %s", one_contrast_output_name))

                for (one_fc_threshold in fc_thresholds) {
                    print(sprintf("FC: %s", one_fc_threshold))

                    # output directory
                    ST_output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/comparisons/%s/%s/%s/FC_%s/00-Visium_recommended/00-Samples/%s", work_dir, comparison, design_output_name, one_contrast_output_name, sub("\\.", "_", one_fc_threshold), one_sample)
                    if (! dir.exists(ST_output_dir)) {
                        dir.create(ST_output_dir, recursive=TRUE, mode="0775")
                    }

                    # get genome 1 specific DE genes
                    comparison_FC_output_dir <- sprintf("%s/comparisons/%s/%s/%s/FC_%s/DE_analysis", output_dir, comparison, design_output_name, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
                    comparison_FC_output_basename <- sprintf("%s_%s_%s%s_%s_FC_%s", dataset, comparison, design_output_name, design_output_name_append, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
                    exp_matrix_1_specific_DE_df <- read.csv(sprintf("%s/%s_gene_ID_matching_with_non_DE_%s_specific.csv", comparison_FC_output_dir, comparison_FC_output_basename, gsub(" ", "", exp_matrix_1_name)), quote="")

                    # get mean expression over all spots for specific DE genes
                    ## expression matrix: genome 1
                    exp_matrix_1_specific_DE_genes_exp_matrix_1_id <- unique(exp_matrix_1_specific_DE_df[[sprintf("gene_id_%s", gsub(" ", "", exp_matrix_1_name))]])
                    exp_matrix_1_specific_DE_ST_sample_mean_exp_vector <- c()
                    for (one_gene in exp_matrix_1_specific_DE_genes_exp_matrix_1_id) {
                        if (one_gene %in% rownames(assay(exp_matrix_1_spe_norm, i = "SCT"))) {
                            one_gene_mean_exp <- mean(assay(exp_matrix_1_spe_norm, i="SCT")[one_gene, rownames(colData(exp_matrix_1_spe_norm))])
                        } else {
                            one_gene_mean_exp <- NA
                        }
                        exp_matrix_1_specific_DE_ST_sample_mean_exp_vector <- c(exp_matrix_1_specific_DE_ST_sample_mean_exp_vector, one_gene_mean_exp)
                    }
                    exp_matrix_1_specific_DE_ST_sample_exp_matrix_1_mean_exp_df <- data.frame(gene=exp_matrix_1_specific_DE_genes_exp_matrix_1_id, exp=exp_matrix_1_specific_DE_ST_sample_mean_exp_vector)
                    colnames(exp_matrix_1_specific_DE_ST_sample_exp_matrix_1_mean_exp_df)[colnames(exp_matrix_1_specific_DE_ST_sample_exp_matrix_1_mean_exp_df)=="exp"] <- sprintf("ST_sample_mean_exp_%s", gsub(" ", "", exp_matrix_1_name))
                    exp_matrix_1_specific_DE_ST_sample_mean_exp_df <- merge(exp_matrix_1_specific_DE_df, exp_matrix_1_specific_DE_ST_sample_exp_matrix_1_mean_exp_df, by.x=sprintf("gene_id_%s", gsub(" ", "", exp_matrix_1_name)), by.y="gene")
                    ## expression matrix: genome 2
                    exp_matrix_1_specific_DE_genes_exp_matrix_2_id <- unique(exp_matrix_1_specific_DE_df[[sprintf("gene_id_%s", gsub(" ", "", exp_matrix_2_name))]])
                    exp_matrix_1_specific_DE_ST_sample_mean_exp_vector <- c()
                    for (one_gene in exp_matrix_1_specific_DE_genes_exp_matrix_2_id) {
                        if (one_gene %in% rownames(assay(exp_matrix_2_spe_norm, i = "SCT"))) {
                            one_gene_mean_exp <- mean(assay(exp_matrix_2_spe_norm, i="SCT")[one_gene, rownames(colData(exp_matrix_2_spe_norm))])
                        } else {
                            one_gene_mean_exp <- NA
                        }
                        exp_matrix_1_specific_DE_ST_sample_mean_exp_vector <- c(exp_matrix_1_specific_DE_ST_sample_mean_exp_vector, one_gene_mean_exp)
                    }
                    exp_matrix_1_specific_DE_ST_sample_exp_matrix_2_mean_exp_df <- data.frame(gene=exp_matrix_1_specific_DE_genes_exp_matrix_2_id, exp=exp_matrix_1_specific_DE_ST_sample_mean_exp_vector)
                    colnames(exp_matrix_1_specific_DE_ST_sample_exp_matrix_2_mean_exp_df)[colnames(exp_matrix_1_specific_DE_ST_sample_exp_matrix_2_mean_exp_df)=="exp"] <- sprintf("ST_sample_mean_exp_%s", gsub(" ", "", exp_matrix_2_name))
                    exp_matrix_1_specific_DE_ST_sample_mean_exp_df <- merge(exp_matrix_1_specific_DE_ST_sample_mean_exp_df, exp_matrix_1_specific_DE_ST_sample_exp_matrix_2_mean_exp_df, by.x=sprintf("gene_id_%s", gsub(" ", "", exp_matrix_2_name)), by.y="gene")


                    ## only keep upregulated DE genes
                    ### get DE results
                    exp_matrix_1_ref_genome_output_dir <- sprintf("%s/%s/%s/%s/FC_%s", output_dir, exp_matrix_1_ref_genome, design_output_name, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
                    exp_matrix_1_output_basename <- sprintf("%s_%s_%s%s_%s_FC_%s", dataset, exp_matrix_1_ref_genome, design_output_name, design_output_name_append, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
                    exp_matrix_1_de_results_file <- sprintf("%s/%s_pval.csv", exp_matrix_1_ref_genome_output_dir, exp_matrix_1_output_basename)
                    exp_matrix_1_de_results_df <- read.csv(exp_matrix_1_de_results_file, quote="", row.names=1)
                    print(sprintf("upregulated %s specific DE genes", exp_matrix_1_name))
                    exp_matrix_1_up_DE_genes <- rownames(exp_matrix_1_de_results_df[which(exp_matrix_1_de_results_df$FDR < 0.05 & exp_matrix_1_de_results_df$logFC > 0),])
                    exp_matrix_1_specific_DE_ST_sample_mean_exp_up_df <- exp_matrix_1_specific_DE_ST_sample_mean_exp_df[which(exp_matrix_1_specific_DE_ST_sample_mean_exp_df[[sprintf("gene_id_%s", gsub(" ", "", exp_matrix_1_name))]] %in% exp_matrix_1_up_DE_genes),]

                    # order exp matrix 1 specific DE genes by mean expression
                    exp_matrix_1_specific_DE_ST_sample_mean_exp_up_df <- exp_matrix_1_specific_DE_ST_sample_mean_exp_up_df[order(exp_matrix_1_specific_DE_ST_sample_mean_exp_up_df[[sprintf("ST_sample_mean_exp_%s", gsub(" ", "", exp_matrix_1_name))]], decreasing=TRUE),]
                    pdf(sprintf("%s/%s_%s_%s_specific_DE_gene_mean_exp_up.pdf", ST_output_dir, comparison_FC_output_basename, one_sample, gsub(" ", "", exp_matrix_1_name)))
                    for (i in 1:dim(exp_matrix_1_specific_DE_ST_sample_mean_exp_up_df)[1]) {
                        exp_matrix_1_mean_exp <- exp_matrix_1_specific_DE_ST_sample_mean_exp_up_df[i, sprintf("ST_sample_mean_exp_%s",  gsub(" ", "", exp_matrix_1_name))]
                        exp_matrix_2_mean_exp <- exp_matrix_1_specific_DE_ST_sample_mean_exp_up_df[i, sprintf("ST_sample_mean_exp_%s",  gsub(" ", "", exp_matrix_2_name))]

                        if (! is.na(exp_matrix_1_mean_exp) & ! is.na(exp_matrix_2_mean_exp)) {
                            exp_matrix_1_gene_id <- exp_matrix_1_specific_DE_ST_sample_mean_exp_up_df[i, sprintf("gene_id_%s",  gsub(" ", "", exp_matrix_1_name))]
                            p1 <- plotVisium(exp_matrix_1_spe_norm, fill=exp_matrix_1_gene_id, assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette="navy") 
                            exp_matrix_2_gene_id <- exp_matrix_1_specific_DE_ST_sample_mean_exp_up_df[i, sprintf("gene_id_%s",  gsub(" ", "", exp_matrix_2_name))]
                            p2 <- plotVisium(exp_matrix_2_spe_norm, fill=exp_matrix_2_gene_id, assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette="navy") 
                            print(wrap_plots(p1, p2) + plot_annotation(title=sprintf("genome 1 ID: %s, genome 2 ID: %s", exp_matrix_1_gene_id, exp_matrix_2_gene_id)))
                        }


                    }
                    dev.off()

                    ## only keep downregulated DE genes
                    print(sprintf("downregulated %s specific DE genes", exp_matrix_1_name))
                    exp_matrix_1_down_DE_genes <- rownames(exp_matrix_1_de_results_df[which(exp_matrix_1_de_results_df$FDR < 0.05 & exp_matrix_1_de_results_df$logFC < 0),])
                    exp_matrix_1_specific_DE_ST_sample_mean_exp_down_df <- exp_matrix_1_specific_DE_ST_sample_mean_exp_df[which(exp_matrix_1_specific_DE_ST_sample_mean_exp_df[[sprintf("gene_id_%s", gsub(" ", "", exp_matrix_1_name))]] %in% exp_matrix_1_down_DE_genes),]

                    # order exp matrix 1 specific DE genes by mean expression
                    exp_matrix_1_specific_DE_ST_sample_mean_exp_down_df <- exp_matrix_1_specific_DE_ST_sample_mean_exp_down_df[order(exp_matrix_1_specific_DE_ST_sample_mean_exp_down_df[[sprintf("ST_sample_mean_exp_%s", gsub(" ", "", exp_matrix_1_name))]], decreasing=TRUE),]
                    pdf(sprintf("%s/%s_%s_%s_specific_DE_gene_mean_exp_down.pdf", ST_output_dir, comparison_FC_output_basename, one_sample, gsub(" ", "", exp_matrix_1_name)))
                    for (i in 1:dim(exp_matrix_1_specific_DE_ST_sample_mean_exp_down_df)[1]) {
                        exp_matrix_1_mean_exp <- exp_matrix_1_specific_DE_ST_sample_mean_exp_down_df[i, sprintf("ST_sample_mean_exp_%s",  gsub(" ", "", exp_matrix_1_name))]
                        exp_matrix_2_mean_exp <- exp_matrix_1_specific_DE_ST_sample_mean_exp_down_df[i, sprintf("ST_sample_mean_exp_%s",  gsub(" ", "", exp_matrix_2_name))]

                        if (! is.na(exp_matrix_1_mean_exp) & ! is.na(exp_matrix_2_mean_exp)) {
                            exp_matrix_1_gene_id <- exp_matrix_1_specific_DE_ST_sample_mean_exp_down_df[i, sprintf("gene_id_%s",  gsub(" ", "", exp_matrix_1_name))]
                            p1 <- plotVisium(exp_matrix_1_spe_norm, fill=exp_matrix_1_gene_id, assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette="navy") 
                            exp_matrix_2_gene_id <- exp_matrix_1_specific_DE_ST_sample_mean_exp_down_df[i, sprintf("gene_id_%s",  gsub(" ", "", exp_matrix_2_name))]
                            p2 <- plotVisium(exp_matrix_2_spe_norm, fill=exp_matrix_2_gene_id, assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette="navy") 
                            print(wrap_plots(p1, p2) + plot_annotation(title=sprintf("genome 1 ID: %s, genome 2 ID: %s", exp_matrix_1_gene_id, exp_matrix_2_gene_id)))
                        }


                    }
                    dev.off()


                    # get genome 2 specific DE genes
                    comparison_FC_output_dir <- sprintf("%s/comparisons/%s/%s/%s/FC_%s/DE_analysis", output_dir, comparison, design_output_name, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
                    comparison_FC_output_basename <- sprintf("%s_%s_%s%s_%s_FC_%s", dataset, comparison, design_output_name, design_output_name_append, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
                    exp_matrix_2_specific_DE_df <- read.csv(sprintf("%s/%s_gene_ID_matching_with_non_DE_%s_specific.csv", comparison_FC_output_dir, comparison_FC_output_basename, gsub(" ", "", exp_matrix_2_name)), quote="")

                    # get mean expression over all spots for specific DE genes
                    ## expression matrix: genome 1
                    exp_matrix_2_specific_DE_genes_exp_matrix_1_id <- unique(exp_matrix_2_specific_DE_df[[sprintf("gene_id_%s", gsub(" ", "", exp_matrix_1_name))]])
                    exp_matrix_2_specific_DE_ST_sample_mean_exp_vector <- c()
                    for (one_gene in exp_matrix_2_specific_DE_genes_exp_matrix_1_id) {
                        if (one_gene %in% rownames(assay(exp_matrix_1_spe_norm, i = "SCT"))) {
                            one_gene_mean_exp <- mean(assay(exp_matrix_1_spe_norm, i="SCT")[one_gene, rownames(colData(exp_matrix_1_spe_norm))])
                        } else {
                            one_gene_mean_exp <- NA
                        }
                        exp_matrix_2_specific_DE_ST_sample_mean_exp_vector <- c(exp_matrix_2_specific_DE_ST_sample_mean_exp_vector, one_gene_mean_exp)
                    }
                    exp_matrix_2_specific_DE_ST_sample_exp_matrix_1_mean_exp_df <- data.frame(gene=exp_matrix_2_specific_DE_genes_exp_matrix_1_id, exp=exp_matrix_2_specific_DE_ST_sample_mean_exp_vector)
                    colnames(exp_matrix_2_specific_DE_ST_sample_exp_matrix_1_mean_exp_df)[colnames(exp_matrix_2_specific_DE_ST_sample_exp_matrix_1_mean_exp_df)=="exp"] <- sprintf("ST_sample_mean_exp_%s", gsub(" ", "", exp_matrix_1_name))
                    exp_matrix_2_specific_DE_ST_sample_mean_exp_df <- merge(exp_matrix_2_specific_DE_df, exp_matrix_2_specific_DE_ST_sample_exp_matrix_1_mean_exp_df, by.x=sprintf("gene_id_%s", gsub(" ", "", exp_matrix_1_name)), by.y="gene")
                    ## expression matrix: genome 2
                    exp_matrix_2_specific_DE_genes_exp_matrix_2_id <- unique(exp_matrix_2_specific_DE_df[[sprintf("gene_id_%s", gsub(" ", "", exp_matrix_2_name))]])
                    exp_matrix_2_specific_DE_ST_sample_mean_exp_vector <- c()
                    for (one_gene in exp_matrix_2_specific_DE_genes_exp_matrix_2_id) {
                        if (one_gene %in% rownames(assay(exp_matrix_2_spe_norm, i = "SCT"))) {
                            one_gene_mean_exp <- mean(assay(exp_matrix_2_spe_norm, i="SCT")[one_gene, rownames(colData(exp_matrix_2_spe_norm))])
                        } else {
                            one_gene_mean_exp <- NA
                        }
                        exp_matrix_2_specific_DE_ST_sample_mean_exp_vector <- c(exp_matrix_2_specific_DE_ST_sample_mean_exp_vector, one_gene_mean_exp)
                    }
                    exp_matrix_2_specific_DE_ST_sample_exp_matrix_2_mean_exp_df <- data.frame(gene=exp_matrix_2_specific_DE_genes_exp_matrix_2_id, exp=exp_matrix_2_specific_DE_ST_sample_mean_exp_vector)
                    colnames(exp_matrix_2_specific_DE_ST_sample_exp_matrix_2_mean_exp_df)[colnames(exp_matrix_2_specific_DE_ST_sample_exp_matrix_2_mean_exp_df)=="exp"] <- sprintf("ST_sample_mean_exp_%s", gsub(" ", "", exp_matrix_2_name))
                    exp_matrix_2_specific_DE_ST_sample_mean_exp_df <- merge(exp_matrix_2_specific_DE_ST_sample_mean_exp_df, exp_matrix_2_specific_DE_ST_sample_exp_matrix_2_mean_exp_df, by.x=sprintf("gene_id_%s", gsub(" ", "", exp_matrix_2_name)), by.y="gene")

                    ## only keep upregulated DE genes
                    ### get DE results
                    exp_matrix_2_ref_genome_output_dir <- sprintf("%s/%s/%s/%s/FC_%s", output_dir, exp_matrix_2_ref_genome, design_output_name, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
                    exp_matrix_2_output_basename <- sprintf("%s_%s_%s%s_%s_FC_%s", dataset, exp_matrix_2_ref_genome, design_output_name, design_output_name_append, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
                    exp_matrix_2_de_results_file <- sprintf("%s/%s_pval.csv", exp_matrix_2_ref_genome_output_dir, exp_matrix_2_output_basename)
                    exp_matrix_2_de_results_df <- read.csv(exp_matrix_2_de_results_file, quote="", row.names=1)
                    print(sprintf("upregulated %s specific DE genes", exp_matrix_2_name))
                    exp_matrix_2_up_DE_genes <- rownames(exp_matrix_2_de_results_df[which(exp_matrix_2_de_results_df$FDR < 0.05 & exp_matrix_2_de_results_df$logFC > 0),])
                    exp_matrix_2_specific_DE_ST_sample_mean_exp_up_df <- exp_matrix_2_specific_DE_ST_sample_mean_exp_df[which(exp_matrix_2_specific_DE_ST_sample_mean_exp_df[[sprintf("gene_id_%s", gsub(" ", "", exp_matrix_2_name))]] %in% exp_matrix_2_up_DE_genes),]

                    # order exp matrix 2 specific DE genes by mean expression
                    exp_matrix_2_specific_DE_ST_sample_mean_exp_up_df <- exp_matrix_2_specific_DE_ST_sample_mean_exp_up_df[order(exp_matrix_2_specific_DE_ST_sample_mean_exp_up_df[[sprintf("ST_sample_mean_exp_%s", gsub(" ", "", exp_matrix_1_name))]], decreasing=TRUE),]
                    pdf(sprintf("%s/%s_%s_%s_specific_DE_gene_mean_exp_up.pdf", ST_output_dir, comparison_FC_output_basename, one_sample, gsub(" ", "", exp_matrix_2_name)))
                    for (i in 1:dim(exp_matrix_2_specific_DE_ST_sample_mean_exp_up_df)[1]) {
                        exp_matrix_1_mean_exp <- exp_matrix_2_specific_DE_ST_sample_mean_exp_up_df[i, sprintf("ST_sample_mean_exp_%s",  gsub(" ", "", exp_matrix_1_name))]
                        exp_matrix_2_mean_exp <- exp_matrix_2_specific_DE_ST_sample_mean_exp_up_df[i, sprintf("ST_sample_mean_exp_%s",  gsub(" ", "", exp_matrix_2_name))]

                        if (! is.na(exp_matrix_1_mean_exp) & ! is.na(exp_matrix_2_mean_exp)) {
                            exp_matrix_1_gene_id <- exp_matrix_2_specific_DE_ST_sample_mean_exp_up_df[i, sprintf("gene_id_%s",  gsub(" ", "", exp_matrix_1_name))]
                            p1 <- plotVisium(exp_matrix_1_spe_norm, fill=exp_matrix_1_gene_id, assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette="navy") 
                            exp_matrix_2_gene_id <- exp_matrix_2_specific_DE_ST_sample_mean_exp_up_df[i, sprintf("gene_id_%s",  gsub(" ", "", exp_matrix_2_name))]
                            p2 <- plotVisium(exp_matrix_2_spe_norm, fill=exp_matrix_2_gene_id, assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette="navy") 
                            print(wrap_plots(p1, p2) + plot_annotation(title=sprintf("genome 1 ID: %s, genome 2 ID: %s", exp_matrix_1_gene_id, exp_matrix_2_gene_id)))
                        }
                    }
                    dev.off()

                    ## only keep downregulated DE genes
                    print(sprintf("downregulated %s specific DE genes", exp_matrix_2_name))
                    exp_matrix_2_down_DE_genes <- rownames(exp_matrix_2_de_results_df[which(exp_matrix_2_de_results_df$FDR < 0.05 & exp_matrix_2_de_results_df$logFC < 0),])
                    exp_matrix_2_specific_DE_ST_sample_mean_exp_down_df <- exp_matrix_2_specific_DE_ST_sample_mean_exp_df[which(exp_matrix_2_specific_DE_ST_sample_mean_exp_df[[sprintf("gene_id_%s", gsub(" ", "", exp_matrix_2_name))]] %in% exp_matrix_2_down_DE_genes),]

                    # order exp matrix 2 specific DE genes by mean expression
                    exp_matrix_2_specific_DE_ST_sample_mean_exp_down_df <- exp_matrix_2_specific_DE_ST_sample_mean_exp_down_df[order(exp_matrix_2_specific_DE_ST_sample_mean_exp_down_df[[sprintf("ST_sample_mean_exp_%s", gsub(" ", "", exp_matrix_1_name))]], decreasing=TRUE),]
                    pdf(sprintf("%s/%s_%s_%s_specific_DE_gene_mean_exp_down.pdf", ST_output_dir, comparison_FC_output_basename, one_sample, gsub(" ", "", exp_matrix_2_name)))
                    for (i in 1:dim(exp_matrix_2_specific_DE_ST_sample_mean_exp_down_df)[1]) {
                        exp_matrix_1_mean_exp <- exp_matrix_2_specific_DE_ST_sample_mean_exp_down_df[i, sprintf("ST_sample_mean_exp_%s",  gsub(" ", "", exp_matrix_1_name))]
                        exp_matrix_2_mean_exp <- exp_matrix_2_specific_DE_ST_sample_mean_exp_down_df[i, sprintf("ST_sample_mean_exp_%s",  gsub(" ", "", exp_matrix_2_name))]

                        if (! is.na(exp_matrix_1_mean_exp) & ! is.na(exp_matrix_2_mean_exp)) {
                            exp_matrix_1_gene_id <- exp_matrix_2_specific_DE_ST_sample_mean_exp_down_df[i, sprintf("gene_id_%s",  gsub(" ", "", exp_matrix_1_name))]
                            p1 <- plotVisium(exp_matrix_1_spe_norm, fill=exp_matrix_1_gene_id, assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette="navy") 
                            exp_matrix_2_gene_id <- exp_matrix_2_specific_DE_ST_sample_mean_exp_down_df[i, sprintf("gene_id_%s",  gsub(" ", "", exp_matrix_2_name))]
                            p2 <- plotVisium(exp_matrix_2_spe_norm, fill=exp_matrix_2_gene_id, assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette="navy") 
                            print(wrap_plots(p1, p2) + plot_annotation(title=sprintf("genome 1 ID: %s, genome 2 ID: %s", exp_matrix_1_gene_id, exp_matrix_2_gene_id)))
                        }
                    }
                    dev.off()
                }
            }
        }
    }
}





# genome comparison specific DE genes expression visualization on ST datasets: plots for JOBIM poster
one_sample <- "B_L4_S14"
refseq_genes_to_plot <- c("Rsu1", "Lsm11", "Kcng2", "Lmo3", "Usp22", "Adcy1", "Nlgn3", "Ikzf4")
ensembl_genes_to_plot <- c("ENSRNOG00000017595", "ENSRNOG00000005450", "ENSRNOG00000053640", "ENSRNOG00000047450", "ENSRNOG00000032492", "ENSRNOG00000059479", "ENSRNOG00000003812", "ENSRNOG00000005535")

# output directory
ST_output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/comparisons/%s/%s", work_dir, comparison, design_output_name)
if (! dir.exists(ST_output_dir)) {
    dir.create(ST_output_dir, recursive=TRUE, mode="0775")
}
pdf(sprintf("%s/%s_JOBIM_poster_selected_genes.pdf", ST_output_dir, one_sample))
for (i in 1:length(refseq_genes_to_plot)) {
    p1 <- plotVisium(exp_matrix_1_spe_norm, fill=refseq_genes_to_plot[i], assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette="black")
    p1 <- p1 + theme(legend.title=element_blank())
    p1 <- p1 + theme(strip.text=element_blank())
    p2 <- plotVisium(exp_matrix_2_spe_norm, fill=ensembl_genes_to_plot[i], assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette="black")
    p2 <- p2 + theme(legend.title=element_blank())
    p2 <- p2 + theme(strip.text=element_blank())
    #print(wrap_plots(p1, p2, nrow=2))
    print(p1)
    print(p2)
    
    # p1 <- plotVisium(exp_matrix_1_spe_norm, fill=refseq_genes_to_plot[i], assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette="navy")
    # p1 <- p1 + theme(legend.title=element_blank())
    # p1 <- p1 + theme(strip.text=element_blank())
    # p2 <- plotVisium(exp_matrix_2_spe_norm, fill=ensembl_genes_to_plot[i], assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette="navy")
    # p2 <- p2 + theme(legend.title=element_blank())
    # p2 <- p2 + theme(strip.text=element_blank())
    # print(wrap_plots(p1, p2, nrow=2))
    
    # print(p1)
    # print(p2)
}
dev.off()


for (i in 1:length(refseq_genes_to_plot)) {
    # black
    p1 <- plotVisium(exp_matrix_1_spe_norm, fill=refseq_genes_to_plot[i], assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette="black")
    p1 <- p1 + theme(legend.title=element_blank())
    p1 <- p1 + theme(strip.text=element_blank())
    p2 <- plotVisium(exp_matrix_2_spe_norm, fill=ensembl_genes_to_plot[i], assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette="black")
    p2 <- p2 + theme(legend.title=element_blank())
    p2 <- p2 + theme(strip.text=element_blank())
    #print(wrap_plots(p1, p2, nrow=2))
    png(sprintf("%s/%s_%s_%s_RefSeq_108.png", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i]), width=2048, height=2048, res=450)
    # jpeg(sprintf("%s/%s_%s_%s_RefSeq_108.jpeg", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i]), quality=100)
    # tiff(sprintf("%s/%s_%s_%s_RefSeq_108.tiff", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i]), width=1024, height=1024, res=150)
    print(p1)
    dev.off()
    
    png(sprintf("%s/%s_%s_%s_Ensembl_105.png", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i]), width=2048, height=2048, res=450)
    # jpeg(sprintf("%s/%s_%s_%s_Ensembl_105.jpeg", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i]), quality=100)
    # tiff(sprintf("%s/%s_%s_%s_Ensembl_105.tiff", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i]), width=1024, height=1024, res=150)
    print(p2)
    dev.off()
    
    # yellow
    color <- "yellow"
    p1 <- plotVisium(exp_matrix_1_spe_norm, fill=refseq_genes_to_plot[i], assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette=color)
    p1 <- p1 + theme(legend.title=element_blank())
    p1 <- p1 + theme(strip.text=element_blank())
    p2 <- plotVisium(exp_matrix_2_spe_norm, fill=ensembl_genes_to_plot[i], assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette=color)
    p2 <- p2 + theme(legend.title=element_blank())
    p2 <- p2 + theme(strip.text=element_blank())
    #print(wrap_plots(p1, p2, nrow=2))
    png(sprintf("%s/%s_%s_%s_RefSeq_108_%s.png", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i], color), width=2048, height=2048, res=450)
    # jpeg(sprintf("%s/%s_%s_%s_RefSeq_108.jpeg", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i]), quality=100)
    # tiff(sprintf("%s/%s_%s_%s_RefSeq_108.tiff", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i]), width=1024, height=1024, res=150)
    print(p1)
    dev.off()
    
    png(sprintf("%s/%s_%s_%s_Ensembl_105_%s.png", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i], color), width=2048, height=2048, res=450)
    # jpeg(sprintf("%s/%s_%s_%s_Ensembl_105.jpeg", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i]), quality=100)
    # tiff(sprintf("%s/%s_%s_%s_Ensembl_105.tiff", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i]), width=1024, height=1024, res=150)
    print(p2)
    dev.off()
    
    # green
    color <- "green"
    p1 <- plotVisium(exp_matrix_1_spe_norm, fill=refseq_genes_to_plot[i], assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette=color)
    p1 <- p1 + theme(legend.title=element_blank())
    p1 <- p1 + theme(strip.text=element_blank())
    p2 <- plotVisium(exp_matrix_2_spe_norm, fill=ensembl_genes_to_plot[i], assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette=color)
    p2 <- p2 + theme(legend.title=element_blank())
    p2 <- p2 + theme(strip.text=element_blank())
    #print(wrap_plots(p1, p2, nrow=2))
    png(sprintf("%s/%s_%s_%s_RefSeq_108_%s.png", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i], color), width=2048, height=2048, res=450)
    # jpeg(sprintf("%s/%s_%s_%s_RefSeq_108.jpeg", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i]), quality=100)
    # tiff(sprintf("%s/%s_%s_%s_RefSeq_108.tiff", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i]), width=1024, height=1024, res=150)
    print(p1)
    dev.off()
    
    png(sprintf("%s/%s_%s_%s_Ensembl_105_%s.png", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i], color), width=2048, height=2048, res=450)
    # jpeg(sprintf("%s/%s_%s_%s_Ensembl_105.jpeg", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i]), quality=100)
    # tiff(sprintf("%s/%s_%s_%s_Ensembl_105.tiff", ST_output_dir, one_sample, refseq_genes_to_plot[i], ensembl_genes_to_plot[i]), width=1024, height=1024, res=150)
    print(p2)
    dev.off()
    
    
    
    
    # p1 <- plotVisium(exp_matrix_1_spe_norm, fill=refseq_genes_to_plot[i], assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette="navy")
    # p1 <- p1 + theme(legend.title=element_blank())
    # p1 <- p1 + theme(strip.text=element_blank())
    # p2 <- plotVisium(exp_matrix_2_spe_norm, fill=ensembl_genes_to_plot[i], assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", palette="navy")
    # p2 <- p2 + theme(legend.title=element_blank())
    # p2 <- p2 + theme(strip.text=element_blank())
    # print(wrap_plots(p1, p2, nrow=2))
    
    # print(p1)
    # print(p2)
}








# load data
## expression matrices
### Space Ranger
#h5_file <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs/raw_feature_bc_matrix.h5", work_dir, genome_name, sample)
#h5_data <- Read10X_h5(h5_file)
### ST Pipeline


# identify mitochondrial genes
#is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe_in_tissue_norm)$gene_name)
#table(is_mito)
#rowData(spe)$gene_name[is_mito]

# plots
#pdf(sprintf("%s/raw_normalized_counts.pdf", output_dir))
    ## total counts and number of detected genes
    hist(colData(spe_in_tissue_norm)$counts.sum, breaks=20, main="Sum of counts per spot", xlab="Sum of raw counts")
    hist(colData(spe_in_tissue_norm)$counts.detected, breaks=20, main="Detected genes per spot", xlab="Detected genes - raw counts")
    hist(colData(spe_in_tissue_norm)$SCT.sum, breaks=20, main="Sum of counts per spot", xlab="Sum of normalized counts")
    hist(colData(spe_in_tissue_norm)$SCT.detected, breaks=20, main="Detected genes per spot", xlab="Detected genes - normalized counts")
    ## splot plots
    ### raw counts
    print(plotSpots(spe_in_tissue_norm, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="counts.sum", palette=c("white", "black")))
    print(plotVisium(spe_in_tissue_norm, fill="counts.sum", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    print(plotSpots(spe_in_tissue_norm, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="counts.detected", palette=c("white", "black")))
    print(plotVisium(spe_in_tissue_norm, fill="counts.detected", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    ### normalized counts
    print(plotSpots(spe_in_tissue_norm, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="SCT.sum", palette=c("white", "black")))
    print(plotVisium(spe_in_tissue_norm, fill="SCT.sum", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    print(plotSpots(spe_in_tissue_norm, x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres", annotate="SCT.detected", palette=c("white", "black")))
    print(plotVisium(spe_in_tissue_norm, fill="SCT.detected", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres"))
    ### raw and normalized counts
    p1 <- plotVisium(spe_in_tissue_norm, fill="counts.sum", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    p2 <- plotVisium(spe_in_tissue_norm, fill="SCT.sum", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    print(wrap_plots(p1, p2) + plot_annotation(title="Sum of counts"))
    p1 <- plotVisium(spe_in_tissue_norm, fill="counts.detected", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    p2 <- plotVisium(spe_in_tissue_norm, fill="SCT.detected", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")
    print(wrap_plots(p1, p2) + plot_annotation(title="Sum of counts"))
    dev.off()

plotVisium(spe_in_tissue_norm, fill="Cntnap5c", assay="SCT", x_coord="pxl_col_in_fullres", y_coord="pxl_row_in_fullres")








genome_name <- "NCBIRefSeq106_NCBIRefSeq106GTF"
#samples <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")
sample <- "A_L1_S1"


print(sprintf("sample: %s", sample))
output_dir <- sprintf("%s/20-Data_analysis/10-EpiReg_data/output/00-ST_Pipeline/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s/10-normalization", work_dir, genome_name, sample)
if (! dir.exists(output_dir)) {
    dir.create(output_dir, recursive=TRUE, mode="0775")
}

# load data
## expression matrices
### Space Ranger
#h5_file <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Samples/%s/10-Pipeline/outs/raw_feature_bc_matrix.h5", work_dir, genome_name, sample)
#h5_data <- Read10X_h5(h5_file)
### ST Pipeline
st_pipeline_matrix_file <- sprintf("%s/10-ST_analysis/00-ST_Pipeline/output/10-TSO_polyA_R1hardtrim1_ov5_n2_min20/%s/00-Visium_recommended/00-Samples/%s/10-Pipeline/%s_stdata.tsv", work_dir, genome_name, sample, sample)
st_pipeline_matrix <- read.table(st_pipeline_matrix_file, sep="\t", header=TRUE, quote="", row.names=1)