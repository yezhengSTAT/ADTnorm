#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(Seurat)
library(dplyr)
library(ggplot2)
library(httpgd)

method <- args[1] #"Arcsinh_b5"
run_name = args[2] #"publicData_CITEseq"
# run_name = "public13Dataset_CITEseq"

master_path = "./"
in_path = "/publicData_CITEseq/data/"
out_path = paste0(master_path, "manuscript/results/", run_name)
fig_path = paste0(out_path, "/Figures/")

## load data directly
marker_list = c("CD3", "CD4", "CD8",  "CD14", "CD19", "CD25", "CD45RA", "CD56", "CD127")
adt_data = readRDS(file = paste0(out_path, "/RDS/adt_data_RawCount_", run_name, ".rds"))
adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", run_name, ".rds"))

## RPCA integration methods
# creates a Seurat object based on the CITE-seq data
rna_data = readRDS(file = './results/publicData_CITEseq/RDS/rna_data_common_RawCount_public13Dataset_CITEseq.rds')
citeseq_obj = CreateSeuratObject(counts = t(rna_data))
citeseq_obj = AddMetaData(citeseq_obj, metadata = adt_feature)

cite_list = SplitObject(citeseq_obj, split.by = "study_name")
cite_list <- lapply(X = cite_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000) ##
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = cite_list)
cite_list <- lapply(X = cite_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
immune.anchors <- FindIntegrationAnchors(object.list = cite_list, anchor.features = features, reduction = "rpca") ##
# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 100, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:50)
saveRDS(immune.combined, file = paste0(out_path, "/RDS/citeseq_obj_RPCA_", run_name, ".rds"))

citeNorm = readRDS(file = paste0(out_path, "/RDS/citeseq_obj_RPCA_", run_name, ".rds"))

adt_data_adtnorm = readRDS(file = paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds"))
adt_data_adtnorm[is.na(adt_data_adtnorm)] = 0

rownames(adt_data_adtnorm) = colnames(citeNorm)
adt_assay <- CreateAssay5Object(counts = t(adt_data_adtnorm), data = t(adt_data_adtnorm))
citeNorm[["ADT"]] <- adt_assay
DefaultAssay(citeNorm) <- 'ADT'
VariableFeatures(citeNorm) <- rownames(citeNorm[["ADT"]])
## set scale.data
if(method %in% c("fastMNN_study", "fastMNN_sample", "logCPM", "decontPro")){
    citeNorm = ScaleData(citeNorm)
}else{
    citeNorm <- SetAssayData(citeNorm, assay = "ADT", slot = "scale.data", new.data = t(adt_data_adtnorm))
}
citeNorm <- RunPCA(citeNorm, reduction.name = 'apca')

DefaultAssay(citeNorm) <- 'integrated'
citeNorm <- FindMultiModalNeighbors(
  citeNorm, reduction.list = list("pca", "apca"), 
  dims.list = list(1:15, 1:8), modality.weight.name = "RNA.weight", k.nn = 30
)

citeNorm <- RunUMAP(citeNorm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", min.dist = 0.01)
citeNorm <- FindClusters(citeNorm, graph.name = "wsnn", algorithm = 1, resolution = 0.15, verbose = FALSE)
saveRDS(citeNorm, file = paste0(out_path, "/RDS/citeseq_obj_RPCA_citeNorm_WNN_15_8_withCluster_noScalePCA_nneighbor20_mindist0.01_", run_name, "_", method, ".rds"))

ans_umap = Embeddings(citeNorm, reduction = "wnn.umap") %>% data.frame
colnames(ans_umap) = c("X1", "X2")
adt_feature$seurat_clusters = citeNorm@meta.data$seurat_clusters


pdf(paste0("./manuscript/results/public13Dataset_CITEseq/Figures/WNN/citeseq_obj_RPCA_citeNorm_WNN_15_8_withCluster_noScalePCA_nneighbor20_mindist0.01_", run_name, "_", method, ".pdf"), width = 15, height = 11)

reindex1 = which(adt_feature$cell_type_l1 == "undefined")
reindex2 = which(adt_feature$cell_type_l1 != "undefined")
reindex = c(reindex1, reindex2)

# print(plot_umap_raster(ans_umap[reindex, ], adt_feature[reindex, ], point_size = 0.3, parameter_list = list(
# target_feature = "seurat_clusters", 
# color_design = colorRampPalette(brewer.pal(8, "Spectral"))(length(adt_feature$seurat_clusters %>% unique)),
# method_label = method
# )) #+ scale_color_brewer(palette = "Dark2")
# )


print(plot_umap_raster(ans_umap[reindex, ], adt_feature[reindex, ], point_size = 0.3, parameter_list = list(
target_feature = "sample", 
color_design = colorRampPalette(brewer.pal(8, "Spectral"))(length(adt_feature$sample %>% unique)),
method_label = method
)) #+ scale_color_brewer(palette = "Dark2")
)

print(plot_umap_raster(ans_umap[reindex, ], adt_feature[reindex, ], point_size = 0.3, parameter_list = list(
target_feature = "study_name", 
color_design = colorRampPalette(brewer.pal(8, "Dark2"))(length(adt_feature$study_name %>% unique)),
method_label = method
)))

print(plot_umap_raster(ans_umap[reindex, ], adt_feature[reindex, ], point_size = 0.3, parameter_list = list(
target_feature = "sample_status", 
color_design = colorRampPalette(brewer.pal(8, "Set1"))(length(adt_feature$sample_status %>% unique)),
method_label = method
)))

adt_feature$cell_type_l1 = factor(adt_feature$cell_type_l1, levels = c("B", "CD4 T", "CD8 T", "monocytes", "NK", "DCs", "undefined"))
print(plot_umap_raster(ans_umap[reindex, ], adt_feature[reindex, ], point_size = 0.3, parameter_list = list(
target_feature = "cell_type_l1", 
color_design = c("#B00000", "#3F918B", "#896191", "#FF980A", "#226fa7", "#F781BF", "grey"), 
#color_design = c(colorRampPalette(brewer.pal(8, "Set1"))(6), "grey"),
color_break = c("B", "CD4 T", "CD8 T", "monocytes", "NK", "DCs", "undefined"),
method_label = method
)))

adt_feature$cell_type_l2 = factor(adt_feature$cell_type_l2, levels = c("naive B", "memory B", 
    "naive CD4", "memory CD4", "Treg", 
    "naive CD8", "memory CD8", 
    "classical monocyte", "intermediate monocyte", "non-classical CD16+ monocyte",
    "CD16- NK", "CD16+ NK", 
    "myeloid DC", "plasmacytoid DC", "undefined"))
print(plot_umap_raster(ans_umap[reindex, ], adt_feature[reindex, ], point_size = 0.3, parameter_list = list(
target_feature = "cell_type_l2",
color_design = c(
    "#B00000", "#FF3380", 
    "#9ef954", "#3F918B", "#03500e",
    "#896191", "#350154",  
    "#FF980A", "#A78300", "#DB6400",
    "#1000a0", "#226fa7",
    "#F781BF", "#7B0054", "grey"),
#color_design = c(colorRampPalette(brewer.pal(8, "Dark2"))(14), "grey"),
color_break = c(
    "naive B", "memory B", 
    "naive CD4", "memory CD4", "Treg", 
    "naive CD8", "memory CD8", 
    "classical monocyte", "intermediate monocyte", "non-classical CD16+ monocyte",
    "CD16- NK", "CD16+ NK", 
    "myeloid DC", "plasmacytoid DC", "undefined"),
method_label = method
)))

dev.off()

pdf(paste0("./manuscript/results/public13Dataset_CITEseq/Figures/WNN/citeseq_obj_RPCA_citeNorm_WNN_15_8_withCluster_noScalePCA_nneighbor20_mindist0.01_", run_name, "_", method, "_weights.pdf"), width = 18, height = 9)
VlnPlot(citeNorm, features = "ADT.weight", group.by = 'cell_type_l1', sort = TRUE, pt.size = 0) +
scale_fill_manual(values = c("#B00000", "#3F918B", "#896191", "#FF980A", "#226fa7", "#F781BF", "grey"), 
breaks = c("B", "CD4 T", "CD8 T", "monocytes", "NK", "DCs", "undefined"))
VlnPlot(citeNorm, features = "ADT.weight", group.by = 'cell_type_l2', sort = TRUE, pt.size = 0) +
scale_fill_manual(values = c(
    "#B00000", "#FF3380", 
    "#9ef954", "#3F918B", "#03500e",
    "#896191", "#350154",  
    "#FF980A", "#A78300", "#DB6400",
    "#1000a0", "#226fa7",
    "#F781BF", "#7B0054", "grey"),
breaks = c(
    "naive B", "memory B", 
    "naive CD4", "memory CD4", "Treg", 
    "naive CD8", "memory CD8", 
    "classical monocyte", "intermediate monocyte", "non-classical CD16+ monocyte",
    "CD16- NK", "CD16+ NK", 
    "myeloid DC", "plasmacytoid DC", "undefined"))
VlnPlot(citeNorm, features = "ADT.weight", group.by = 'seurat_clusters', sort = TRUE, pt.size = 0) +
scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Spectral"))(length(adt_feature$seurat_clusters %>% unique)))
dev.off()

