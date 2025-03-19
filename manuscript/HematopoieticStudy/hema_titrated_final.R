library(pacman)
library(ggplot2)
library(viridis)
library(data.table)
library(tidyr)
library(gridExtra)
library(magrittr)
library(umap)
library(ggpubr)
library(SingleCellExperiment)
library(RColorBrewer)
library(ggridges)
library(dplyr)
library(CytofRUV)
library(shapr)
library(fda)
library(pryr)
library(Seurat)
library(systemfonts)
library(shiny)
library(httpgd)
library(ggrepel)
require(magrittr)
require(ggplot2)
require(RColorBrewer)
require(tidyr)
require(ggridges)
require(cytotidyr)
library(harmony)
library(zellkonverter)
library(reticulate)
library(hdf5r)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(ADTnorm)
# p_load(flowStats, flowCore, FlowSOM, ncdfFlow, flowViz, pdfCluster, cluster)

set.seed(20250121)
## =====================
## study name and paths
## =====================
run_name = "hema_titrated_final"
master_path = "./"
out_path = paste0(master_path, "manuscript/results/", run_name)
fig_path = paste0(out_path, "/Figures_", run_name, "/")

## =================
## load titration data
## =================
in_path = "./data/Zhang2024/"
annotation = fread(paste0(in_path, "/GEO/GSE245108_Level3R-Titrated-Titration-Annotation.txt"), header = TRUE) %>% mutate(filename = gsub(".*\\.", "", UID), barcode = gsub("\\..*", "", UID), filename = paste0(sub("_.*", "", filename), "-", sub(".*_", "", filename)))
cd271_obj = Read10X_h5(paste0(in_path, "/GEO/GSE245108_BF21-CD271_filtered_feature_bc_matrix.h5"))
protein_list = rownames(cd271_obj$`Antibody Capture`)

# adt_data = c()
# adt_feature = c()
# for(item_name in paste(rep(c("BF21", "BM27", "WF26", "WM34"), each = 3), rep(c("CD271", "CD34", "TNC"), 4), sep = "-")){
#     obj = Read10X_h5(paste0(in_path, "/GEO/GSE245108_", item_name, "_filtered_feature_bc_matrix.h5"))
#     cell_list = colnames(obj$`Antibody Capture`)[which(colnames(obj$`Antibody Capture`) %in% (annotation %>% dplyr::filter(filename == item_name) %$% barcode))]

#     adt_feature_tmp = annotation %>% dplyr::filter(filename == item_name) %>% data.frame
#     rownames(adt_feature_tmp) = adt_feature_tmp$barcode
#     adt_feature_each = adt_feature_tmp[cell_list, ]
#     adt_data_each = obj$`Antibody Capture`[protein_list, cell_list]
#     print(which(rownames(adt_feature_each) != colnames(adt_data_each)))

#     adt_data = cbind(adt_data, adt_data_each)
#     adt_feature = rbind(adt_feature, adt_feature_each)
# }
# adt_data %>% dim
# adt_feature %>% dim
# colnames(adt_data) %>% head
# rownames(adt_feature) %>% head
# head(adt_feature)
# adt_feature_bp = adt_feature
# adt_feature$donor = factor(adt_feature$filename %>% gsub("-.*", "", .), levels = c("BF21", "BM27", "WF26", "WM34"), labels = c(paste0("Donor", c(1:4))))
# adt_feature$subset = factor(adt_feature$filename %>% gsub(".*-", "", .), levels = c("CD271", "CD34", "TNC"), labels = c("CD271", "CD34", "BMNC"))
# adt_feature$sample = paste0(adt_feature$donor, "_", adt_feature$subset)
# adt_feature$batch = adt_feature$sample

# cell_x_adt = data.frame(t(adt_data))
# cell_x_feature = data.frame(adt_feature)
# colnames(cell_x_adt) = colnames(cell_x_adt) %>% gsub("HuMs\\.", "", .) %>% gsub("Hu\\.", "", .) %>% gsub("Isotype_", "Isotype", .) %>% gsub("_.*", "", .)

# saveRDS(cell_x_adt, paste0(out_path, "/RDS/cell_x_adt_", run_name, ".rds"))
# saveRDS(cell_x_feature, paste0(out_path, "/RDS/cell_x_feature_", run_name, ".rds"))

adt_data = readRDS(paste0(out_path, "/RDS/cell_x_adt_", run_name, ".rds"))
adt_feature = readRDS(paste0(out_path, "/RDS/cell_x_feature_", run_name, ".rds"))

## ================
## QC and ADTnorm
## ================

## check the unique value number for each protein marker for each sample in adt_feature
uniq_num = cbind(adt_data, adt_feature) %>% dplyr::group_by(sample) %>% summarise_all(n_distinct) %>% as.data.frame

# adt_feature$batch = factor(adt_feature$Concentration, levels = c("0.25", "0.5", "1", "2", "4"))
adt_feature$sample = factor(
    adt_feature$sample, 
    levels = c("Donor1_CD271", "Donor1_CD34", "Donor1_BMNC", "Donor2_CD271", "Donor2_CD34", "Donor2_BMNC", "Donor3_CD271", "Donor3_CD34", "Donor3_BMNC", "Donor4_CD271", "Donor4_CD34", "Donor4_BMNC")
)
adt_feature$batch = adt_feature$sample

# res_norm = ADTnorm(customize_landmark = TRUE,
# save_fig = TRUE,
# save_landmark = TRUE,
# cell_x_adt = adt_data,
# cell_x_feature = adt_feature,
# save_outpath = out_path,
# study_name = run_name,
# trimodal_marker = c("CD11a", "CD162", "CD1d", "CD45RB", "CD4", "IgG.Fc", "CD82", "CD47", "CD29", "CD34", "CD61", "CD62P", "CD32"), #c("CD4"),
# bw_smallest_bi = 0.3,
# bw_smallest_tri = 0.3,
# bw_smallest_adjustments = list(CD3 = 0.8, CD4 = 0.8, CD8 = 0.8),
# shoulder_valley_slope = -1,
# marker_to_process = c("IgG.Fc"),
# exclude_zeroes = FALSE,
# bimodal_marker = NULL,
# positive_peak = NULL,
# quantile_clip = 1,
# peak_type = "mode",
# multi_sample_per_batch = FALSE,
# shoulder_valley = TRUE,
# valley_density_adjust = 3,
# landmark_align_type = "negPeak_valley_posPeak",
# midpoint_type = "valley",
# neg_candidate_thres = asinh(2/5 + 1),
# lower_peak_thres = 0.00005,
# brewer_palettes = "Set1",
# detect_outlier_valley = FALSE,
# target_landmark_location = NULL,
# clean_adt_name = FALSE,
# override_landmark = NULL,
# verbose = TRUE)

# saveRDS(res_norm, file = paste0(out_path, "/RDS/adt_data_norm_ADTnorm_sample_manual_", run_name, ".rds"))

## check the unique value number for each protein marker for each sample in adt_feature
uniq_num = cbind(adt_data, adt_feature) %>% dplyr::group_by(sample) %>% summarise_all(n_distinct) %>% as.data.frame
uniq_num[, 2:143] %>% colMeans %>% sort
uniq_num[, 2:143] %>% apply(., 2, median) %>% sort
uniq_num[, 2:143] %>% apply(., 2, max) %>% sort

## only process adt markers that at least have 10 counts per sample per concentration
adt_marker_filter = uniq_num[, 2:143] %>% apply(., 2, max) %>% data.frame %>% dplyr::filter(. >=10) %>% rownames

## =============
## UMAP
## =============
## read in normalized count
# norm_rds_path = "./results/hema_titrated_final/RDS/each_marker/"
# res_norm = matrix(0, nrow = nrow(adt_data), ncol = length(adt_marker_filter))
# colnames(res_norm) = adt_marker_filter
# rownames(res_norm) = rownames(adt_data)

# for(adt_marker_each in adt_marker_filter){
#     tmp = readRDS(paste0(norm_rds_path, "adt_data_norm_ADTnorm_sample_default_", run_name, "_", adt_marker_each, ".rds"))
#     res_norm[, adt_marker_each] = tmp[,1]
# }

# saveRDS(res_norm, file = paste0(out_path, "/RDS/adt_data_norm_ADTnorm_sample_manual_", run_name, ".rds"))

res_norm = readRDS(file = paste0(out_path, "/RDS/adt_data_norm_ADTnorm_sample_manual_", run_name, ".rds"))

# adt_data_select = adt_data[cell_select, ]
# adt_feature_select = adt_feature[cell_select, ]

## add three levels of annotation
label3levels = fread("./data/Zhang2024/Titrated_Datasets/MetaData/scTriangulate-titrated-clusters.txt", header = TRUE)
which(!(label3levels$UID %in% adt_feature$UID))
adt_feature = left_join(adt_feature, label3levels, by = "UID")
head(adt_feature)

## Arcsihn transformation
adt_obj = CreateSeuratObject(counts = t(adt_data[, adt_marker_filter]), data = t(adt_data[, adt_marker_filter]))
adt_obj = SetAssayData(adt_obj, slot = "scale.data", new.data = t(arcsinh_transform(cell_x_adt = adt_data[, adt_marker_filter])))
adt_obj = FindVariableFeatures(adt_obj, nfeature = 2000)
adt_obj = AddMetaData(adt_obj, metadata = adt_feature)
adt_obj = RunPCA(adt_obj)
adt_obj = FindNeighbors(adt_obj, dims = 1:8, reduction = "pca")
adt_obj = FindClusters(adt_obj, resolution = 0.3, algorithm = 2)
adt_obj = RunUMAP(adt_obj, dims = 1:8, reduction = "pca", min.dist = 0.0001, n.neighbors = 50)

DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "seurat_clusters", raster = TRUE) 
p1 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "donor", raster = TRUE) + scale_color_brewer(palette = "Set1")
p2 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "subset", raster = TRUE) + scale_color_brewer(palette = "Set2")
p3 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "Level 1", raster = TRUE) + scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(12))
p4 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "Level 2", raster = TRUE) + scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(27))

ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, align = "hv")
ggsave(paste0(fig_path, "Arcsihn_titrated_umap_pc8_dist0.0001_neighbors50.pdf"), width = 18, height = 10.5)

## CLR transformation
adt_obj = CreateSeuratObject(counts = t(adt_data[, adt_marker_filter]))
adt_obj = NormalizeData(adt_obj, method = "CLR")
adt_obj = ScaleData(adt_obj)
adt_obj = FindVariableFeatures(adt_obj, nfeature = 2000)
adt_obj = AddMetaData(adt_obj, metadata = adt_feature)
adt_obj = RunPCA(adt_obj)
adt_obj = FindNeighbors(adt_obj, dims = 1:8, reduction = "pca")
adt_obj = FindClusters(adt_obj, resolution = 0.3, algorithm = 2)
adt_obj = RunUMAP(adt_obj, dims = 1:8, reduction = "pca", min.dist = 0.0001, n.neighbors = 50)

p1 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "donor", raster = TRUE) + scale_color_brewer(palette = "Set1")
p2 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "subset", raster = TRUE) + scale_color_brewer(palette = "Set2")
p3 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "Level 1", raster = TRUE) + scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(12))
p4 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "Level 2", raster = TRUE) + scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(27))

ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, align = "hv")
ggsave(paste0(fig_path, "CLR_titrated_umap_pc8_dist0.0001_neighbors50.pdf"), width = 18, height = 10.5)

## Harmony
# out <- HarmonyMatrix(
#         data_mat = as.matrix(arcsinh_transform(cell_x_adt = adt_data[, adt_marker_filter])) + 0.0001,
#         meta_data = adt_feature$batch, #adt_feature$sample
#         do_pca = FALSE
# )
# saveRDS(out, file = paste0(out_path, "/RDS/adt_data_norm_Harmony_batch_", run_name, ".rds"))
# out = readRDS(file = paste0(out_path, "/RDS/adt_data_norm_Harmony_batch_", run_name, ".rds"))
out = readRDS(file = paste0(out_path, "/RDS/adt_data_norm_Harmony_donor_", run_name, ".rds"))

adt_obj = CreateSeuratObject(counts = t(adt_data[, adt_marker_filter]), data = t(adt_data[, adt_marker_filter]))
adt_obj = SetAssayData(adt_obj, slot = "scale.data", new.data = t(out))
adt_obj = FindVariableFeatures(adt_obj, nfeature = 2000)
adt_obj = AddMetaData(adt_obj, metadata = adt_feature)
adt_obj = RunPCA(adt_obj)
adt_obj = FindNeighbors(adt_obj, dims = 1:10, reduction = "pca")
adt_obj = FindClusters(adt_obj, resolution = 0.3, algorithm = 2)
adt_obj = RunUMAP(adt_obj, dims = 1:10, reduction = "pca", min.dist = 0.0001, n.neighbors = 50)

p1 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "donor", raster = TRUE) + scale_color_brewer(palette = "Set1")
p2 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "subset", raster = TRUE) + scale_color_brewer(palette = "Set2")
p3 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "Level 1", raster = TRUE) + scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(12))
p4 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "Level 2", raster = TRUE) + scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(27))

ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, align = "hv")
ggsave(paste0(fig_path, "Harmony_donor_titrated_umap_pc8_dist0.0001_neighbors50.pdf"), width = 18, height = 10.5)

## totalVI from the original paper
seurat_obj = readRDS(paste0("./data/Zhang2024/Titrated_Datasets/R-object/titrated-seurat.rds"))

totalVI_input = fread("./data/Zhang2024/Titrated_Datasets/TotalVI-ADTs/ADTs-Titrated-TotalVI-log2-transposed.txt", header = TRUE) %>% data.frame

colnames(totalVI_input)[-1] == gsub("-", "_", rownames(seurat_obj))
rownames(totalVI_input) = totalVI_input$uid
totalVI_input = totalVI_input[, -1]

## the following code prove that original paper use CLR + Scaling to get the scale.data.
# tmp_count = seurat_obj@assays$ADT@counts
# tmp_cpm = NormalizeData(tmp_count, normalization.method = "CLR", margin = 2)
# tmp_scale = ScaleData(tmp_cpm, margin = 2)
library(zellkonverter)
sce <- readH5AD("./data/Zhang2024/Titrated_Datasets/adata_combined_rna_adt_annotated-titrated.h5ad")
sce


adt_obj = CreateSeuratObject(counts = as.matrix(seurat_obj@assays$ADT@counts), data = as.matrix(seurat_obj@assays$ADT@data))
adt_obj = SetAssayData(adt_obj, slot = "scale.data", new.data = as.matrix(t(totalVI_input[colnames(seurat_obj), ])))

adt_obj = FindVariableFeatures(adt_obj, nfeature = 2000)
adt_obj = AddMetaData(adt_obj, metadata = seurat_obj@meta.data)
adt_obj = RunPCA(adt_obj)
adt_obj = FindNeighbors(adt_obj, dims = 1:8, reduction = "pca")
adt_obj = FindClusters(adt_obj, resolution = 0.3, algorithm = 2)
adt_obj = RunUMAP(adt_obj, dims = 1:8, reduction = "pca", min.dist = 0.0001, n.neighbors = 50)
adt_obj@meta.data$donor = rownames(adt_obj@meta.data) %>% gsub(".*\\.", "", .) %>% gsub("_.*", "", .) %>% factor(levels = c("BF21", "BM27", "WF26", "WM34"), labels = c(paste0("Donor", c(1:4))))
adt_obj@meta.data$subset = rownames(adt_obj@meta.data) %>% gsub(".*_", "", .) %>% factor(levels = c("CD271", "CD34", "TNC"), labels = c("CD271", "CD34", "BMNC"))

c(
    adjustedRandIndex(adt_obj@meta.data$seurat_clusters, paste0(adt_obj@meta.data$donor, adt_obj@meta.data$subset)),
    adjustedRandIndex(adt_obj@meta.data$seurat_clusters, adt_obj@meta.data$donor),
    adjustedRandIndex(adt_obj@meta.data$seurat_clusters, adt_obj@meta.data$subset),
    adjustedRandIndex(adt_obj@meta.data$seurat_clusters, adt_obj@meta.data$`Level1`),
    adjustedRandIndex(adt_obj@meta.data$seurat_clusters, adt_obj@meta.data$`Level2`)
)
# 0.26187469 0.13423286 0.07592449 0.23326998 0.23532791
# 0.3347427 0.2301026 0.1064639 0.2361489 0.1855546
#  0.3562089 0.2428871 0.1225801 0.2502295 0.1803498
# 0.062368711 0.005668332 0.100206547 0.533556539 0.474844035
p1 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "donor", raster = TRUE) + scale_color_brewer(palette = "Set1")
p2 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "subset", raster = TRUE) + scale_color_brewer(palette = "Set2")
p3 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "Level1", raster = TRUE) + scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(12))
p4 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "Level2", raster = TRUE) + scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(27))

ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, align = "hv")
ggsave(paste0(fig_path, "totalVI_titrated_umap_pc8_dist0.0001_neighbors50_new.pdf"), width = 18, height = 10.5)

## ADTnorm
adt_obj = CreateSeuratObject(counts = t(adt_data[, adt_marker_filter]), data = t(adt_data[, adt_marker_filter]))
# adt_obj = CreateSeuratObject(counts = t(res_norm), data = t(res_norm))
adt_obj = SetAssayData(adt_obj, slot = "scale.data", new.data = t(res_norm))
adt_obj = FindVariableFeatures(adt_obj, nfeature = 2000)
adt_obj = AddMetaData(adt_obj, metadata = adt_feature)
adt_obj = RunPCA(adt_obj)
adt_obj = FindNeighbors(adt_obj, dims = 1:8, reduction = "pca")
adt_obj = FindClusters(adt_obj, resolution = 0.3, algorithm = 2)
adt_obj = RunUMAP(adt_obj, dims = 1:8, reduction = "pca", min.dist = 0.0001, n.neighbors = 50)

c(
    adjustedRandIndex(adt_obj@meta.data$seurat_clusters, adt_feature$sample),
    adjustedRandIndex(adt_obj@meta.data$seurat_clusters, adt_feature$donor),
    adjustedRandIndex(adt_obj@meta.data$seurat_clusters, adt_feature$subset),
    adjustedRandIndex(adt_obj@meta.data$seurat_clusters, adt_feature$`Level 1`),
    adjustedRandIndex(adt_obj@meta.data$seurat_clusters, adt_feature$`Level 2`)
)

# 0.21222546 0.07141922 0.08759827 0.25287953 0.30561155
# 0.21592994 0.07839284 0.12521183 0.38659873 0.30584253
#  0.10637068 0.03122535 0.10277495 0.44759984 0.38782916
p1 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "donor", raster = TRUE) + scale_color_brewer(palette = "Set1")
p2 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "subset", raster = TRUE) + scale_color_brewer(palette = "Set2")
p3 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "Level 1", raster = TRUE) + scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(12))
p4 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "Level 2", raster = TRUE) + scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(27))

ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, align = "hv")
ggsave(paste0(fig_path, "ADTnorm_titrated_umap_pc8_dist0.0001_neighbors50.pdf"), width = 18, height = 10.5)

## ADTnorm only use the 70 markers recommended.
adtnorm_select_marker = c("CD107a", "CD109", "CD126", "CD138", "CD141", "CD144", "CD207", "CD324", "CD340", "CD49b", "CD57", "CD62E", "CD66b", "EGFR", "IgG.Fc", "LT.bR", "TCR.Vb13.1", "VEGFR.3", "XCR1", "CD10", "CD82", "IgM", "CD106", "CD305", "CD47", "HLA.ABC", "CD19", "CD205", "CD230", "CD5", "CD63", "CD73", "CD13", "CD235ab", "CD29", "CD62L", "CD36", "CD8", "CD117", "CD55", "CD94", "CD271", "CD34", "CD24", "CD61", "CD81", "CD49a", "CD72", "CD11a", "CD46", "CD86", "CD9", "CD123", "CD39", "CD49d", "CD33", "CD43", "CD62P", "HLA.DR", "CD11b", "CD11c", "CD32", "CD244", "CD26", "CD35", "CD41", "CD99", "CLEC12A", "HLA.DR.DP.DQ", "CD4")
adtnorm_select_marker = adtnorm_select_marker[which(adtnorm_select_marker %in% adt_marker_filter)]

adt_obj = CreateSeuratObject(counts = t(adt_data[, adtnorm_select_marker]), data = t(adt_data[, adtnorm_select_marker]))
# adt_obj = CreateSeuratObject(counts = t(res_norm), data = t(res_norm))
adt_obj = SetAssayData(adt_obj, slot = "scale.data", new.data = t(res_norm))
adt_obj = FindVariableFeatures(adt_obj, nfeature = 2000)
adt_obj = AddMetaData(adt_obj, metadata = adt_feature)
adt_obj = RunPCA(adt_obj)
adt_obj = FindNeighbors(adt_obj, dims = 1:30, reduction = "pca")
adt_obj = FindClusters(adt_obj, resolution = 0.3, algorithm = 2)
adt_obj = RunUMAP(adt_obj, dims = 1:30, reduction = "pca", min.dist = 0.0001, n.neighbors = 50)

c(
    adjustedRandIndex(adt_obj@meta.data$seurat_clusters, adt_feature$sample),
    adjustedRandIndex(adt_obj@meta.data$seurat_clusters, adt_feature$donor),
    adjustedRandIndex(adt_obj@meta.data$seurat_clusters, adt_feature$subset),
    adjustedRandIndex(adt_obj@meta.data$seurat_clusters, adt_feature$`Level 1`),
    adjustedRandIndex(adt_obj@meta.data$seurat_clusters, adt_feature$`Level 2`)
)
# [1] 0.09081758 0.01850635 0.12374618 0.46360061 0.35369880

p1 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "donor", raster = TRUE) + scale_color_brewer(palette = "Set1")
p2 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "subset", raster = TRUE) + scale_color_brewer(palette = "Set2")
p3 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "Level 1", raster = TRUE) + scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(12))
p4 = DimPlot(adt_obj, reduction = "umap", label = FALSE, group.by = "Level 2", raster = TRUE) + scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(27))

ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, align = "hv")
ggsave(paste0(fig_path, "ADTnorm_titrated_70recommended_umap_pc30_dist0.0001_neighbors50.pdf"), width = 18, height = 10.5)


## ================
## staining quality
## ================
cell_x_adt = arcsinh_transform(cell_x_adt = adt_data[, adt_marker_filter]) 
cell_x_feature = adt_feature

area_under_curve = function(densityObj, threshold, peak_num){

    if(peak_num > 1){
        index <- which(densityObj$x >= threshold)
        auc <- sum(densityObj$y[index] * diff(c(densityObj$x[index], max(densityObj$x))))
        return(auc)

    }else{
        index <- which(densityObj$x >= threshold)
        auc <- sum(densityObj$y[index] * diff(c(densityObj$x[index], max(densityObj$x))))
        return(min(auc, 1-auc))
    }

}

valley_deep_scaler = function(cell_x_feature, cell_x_adt, adt_marker_select, each_sample, peak_info, valley_info, adjust = 3, peak_num = NA){
    ind = which(cell_x_feature$sample == each_sample)
    # mean_diff = (peak_info[each_sample, ncol(peak_info)] - peak_info[each_sample, 1])
    adt_tmp = cell_x_adt[ind, adt_marker_select]

    density_res = stats::density(
        adt_tmp,
        adjust = adjust , na.rm = TRUE
    )
    x = density_res$x
    y = density_res$y
    if(peak_num >1){
        x_peakR = which(abs(x - peak_info[each_sample, ncol(peak_info)]) == min(abs(x - peak_info[each_sample, ncol(peak_info)])))
        y_peakR = y[x_peakR]
    }else{
        y_peakR = 0
    }
    
    # x_peakL = which(abs(x - peak_info[each_sample, 1]) == min(abs(x - peak_info[each_sample, 1])))
    if(peak_num <=2){
        x_valley = which(abs(x - valley_info[each_sample, 1]) == min(abs(x - valley_info[each_sample, 1])))
    }else{
        x_valley = which(abs(x - valley_info[each_sample, ncol(valley_info)]) == min(abs(x - valley_info[each_sample, ncol(valley_info)])))
    }
    # y_peakL = y[x_peakL]
    if(peak_num > 1){
        y_valley = y[x_valley]
    }else{
        y_valley = 0
    }
    auc_scaler = area_under_curve(density_res, threshold = valley_info[each_sample, 1], peak_num = peak_num)
    deep_scaler = (1 + y_peakR - y_valley) * (1 + auc_scaler)
    return(deep_scaler)
}

two_peak_stain_quality = function(cell_x_feature, cell_x_adt, adt_marker_select, each_sample, peak_info, valley_info, peak_num = 2){

    ind = which(cell_x_feature$sample == each_sample)
    mean_diff = abs(peak_info[each_sample, ncol(peak_info)] - peak_info[each_sample, 1])
    adt_tmp = cell_x_adt[ind, adt_marker_select]
    within_peak_sd = sqrt((sum((adt_tmp[which(adt_tmp < valley_info[each_sample, 1])] - peak_info[each_sample, 1])^2) + sum((adt_tmp[which(adt_tmp > valley_info[each_sample, 1])] - peak_info[each_sample, ncol(peak_info)])^2))/length(adt_tmp))
    
    deep_scaler = valley_deep_scaler(cell_x_feature, cell_x_adt, adt_marker_select, each_sample, peak_info, valley_info, peak_num = peak_num)
    # return(mean_diff/within_peak_sd)
    return((mean_diff*deep_scaler)/within_peak_sd)
}

multi_peak_stain_quality = function(cell_x_feature, cell_x_adt, adt_marker_select, each_sample, peak_info, valley_info, peak_num = 3){

    ind = which(cell_x_feature$sample == each_sample)
    mean_diff = abs(peak_info[each_sample, ncol(peak_info)] - peak_info[each_sample, 1])
    adt_tmp = cell_x_adt[ind, adt_marker_select]
    within_peak_var = 0
    for(col_index in 1:ncol(peak_info)){
        if(col_index == 1){
            ## first peak
            adt_tmp_select = which(adt_tmp < valley_info[each_sample, col_index])
            within_peak_var = within_peak_var + sum((adt_tmp[adt_tmp_select] - peak_info[each_sample, col_index])^2)
        }else if(col_index == ncol(peak_info)){
            ## last peak
            adt_tmp_select = which(adt_tmp > valley_info[each_sample, col_index-1])
            within_peak_var = within_peak_var + sum((adt_tmp[adt_tmp_select] - peak_info[each_sample, col_index])^2)
        }else{
            ## middle peak
            adt_tmp_select = which(adt_tmp < valley_info[each_sample, col_index] & adt_tmp > valley_info[each_sample, col_index-1])
            within_peak_var = within_peak_var + sum((adt_tmp[adt_tmp_select] - peak_info[each_sample, col_index])^2)
        }        
    }
    within_peak_sd = sqrt(within_peak_var/length(adt_tmp))
    
    deep_scaler = valley_deep_scaler(cell_x_feature, cell_x_adt, adt_marker_select, each_sample, peak_info, valley_info, peak_num = peak_num)
    # return(mean_diff/within_peak_sd)
    return((mean_diff*deep_scaler)/within_peak_sd)
}

one_peak_stain_quality = function(cell_x_feature, cell_x_adt, adt_marker_select, each_sample, peak_info, valley_info, peak_num = 1){

    ind = which(cell_x_feature$sample == each_sample)
    adt_tmp = cell_x_adt[ind, adt_marker_select]

    # median_pos = median(adt_tmp[which(adt_tmp > valley_info[each_sample, 1])])
    if(length(peak_info[each_sample,]) == 1){
        mean_diff = abs(valley_info[each_sample, 1] - peak_info[each_sample, 1])
    }else if(is.na(peak_info[each_sample, 1])){
        mean_diff = abs(valley_info[each_sample, 1] - peak_info[each_sample, ncol(peak_info)])
    }else{
        mean_diff = abs(valley_info[each_sample, 1] - peak_info[each_sample, 1])
    }
    within_peak_sd = sqrt(sd(adt_tmp))
    deep_scaler = valley_deep_scaler(cell_x_feature, cell_x_adt, adt_marker_select, each_sample, peak_info, valley_info, peak_num = peak_num)
    # return(mean_diff/within_peak_sd)
    return((mean_diff*deep_scaler)/within_peak_sd)
 

}

## use valley and peak
# peak_num_summary = c()
# peak_sep_summary = c()

# for(adt_marker_select in colnames(cell_x_adt)){ ##
    
#     peak_valley_density = readRDS(paste0(out_path, "/RDS/peak_valley_locations_", adt_marker_select, "_", run_name, ".rds"))
    
#     peak_info = peak_valley_density$peak_landmark_list
#     valley_info = peak_valley_density$valley_landmark_list

#     for(each_sample in unique(adt_feature$sample)){
#         batch_info = each_sample
#         each_peak_info = peak_info[each_sample, ]
#         peak_num = sum(is.na(each_peak_info) == FALSE)
#         peak_num_summary = data.frame(
#             peak_num = peak_num, 
#             batch = batch_info, 
#             sample = each_sample, adt_marker = adt_marker_select) %>% 
#             rbind(peak_num_summary, .)
#         if(peak_num == 1){
#             peak_sep_summary = data.frame(
#                 sep_power = one_peak_stain_quality(adt_feature, cell_x_adt, adt_marker_select, each_sample, peak_info, valley_info, peak_num = peak_num),
#                 adt_marker = adt_marker_select, sample = each_sample, batch = batch_info, peak_num = paste0("# of peak: ", peak_num)) %>% rbind(peak_sep_summary, .)
#         }
#         if(peak_num == 2){
#             peak_sep_summary = data.frame(
#                 sep_power = two_peak_stain_quality(adt_feature, cell_x_adt, adt_marker_select, each_sample, peak_info, valley_info, peak_num = peak_num),
#                 adt_marker = adt_marker_select, sample = each_sample, batch = batch_info, peak_num = paste0("# of peak: ", peak_num)) %>% rbind(peak_sep_summary, .)
#         }
#         if(peak_num > 2){
#             peak_sep_summary = data.frame(
#                 sep_power = multi_peak_stain_quality(adt_feature, cell_x_adt, adt_marker_select, each_sample, peak_info, valley_info, peak_num = peak_num),
#                 adt_marker = adt_marker_select, sample = each_sample, batch = batch_info, peak_num = paste0("# of peak: ", peak_num)) %>% rbind(peak_sep_summary, .)
#         }
#     }
# }
# peak_sep_summary$Donor = factor(gsub("_.*", "", peak_sep_summary$sample), levels = paste0("Donor", c(1:4)))
# peak_sep_summary$Capture = factor(gsub(".*_", "", peak_sep_summary$sample), levels = c("CD271", "CD34", "BMNC"))

# saveRDS(peak_sep_summary, file = paste0(out_path, "/RDS/peak_sep_summary_", run_name, ".rds"))
peak_sep_summary = readRDS(file = paste0(out_path, "/RDS/peak_sep_summary_", run_name, ".rds"))

library(gghalves)

peak_sep_summary %>% 
ggplot(aes(x = Capture, y = sep_power, fill = Donor)) +
geom_half_violin() +
geom_half_point(size = 0.5) +
theme_bw(base_size = 25) +
scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(4, "Set1"))(4)) +
scale_color_brewer(palette = "Set1") +
ylab("Stain Quality Score") +
xlab("") +
theme(legend.position = "top") +
# rotate_x_text(angle = 45) +
scale_y_continuous(trans = "log10") +
geom_hline(yintercept = 2, linewidth =2, linetype = "dashed", color = "grey30")
ggsave(paste0(fig_path, "titrated_stain_quality_score_acrossDonor.pdf"), width = 8, height = 9)

## generate for publication
peak_sep_summary %>% 
ggplot(aes(x = Capture, y = sep_power, fill = Donor)) +
geom_violin() +
geom_jitter(size = 0.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9)) +
theme_bw(base_size = 25) +
scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(4, "Set1"))(4)) +
scale_color_brewer(palette = "Set1") +
ylab("Stain Quality Score") +
xlab("") +
theme(legend.position = "top") +
# rotate_x_text(angle = 45) +
scale_y_continuous(trans = "log10") +
geom_hline(yintercept = 2, linewidth =2, linetype = "dashed", color = "grey30")
ggsave(paste0(fig_path, "titrated_stain_quality_score_acrossDonor_polished.pdf"), width = 8, height = 9)

adt_marker_select = peak_sep_summary %>% dplyr::filter(sep_power > 2) %$% adt_marker %>% table %>% sort %>% data.frame %>% dplyr::filter(Freq > 1) %$% .
adtnorm_select_marker = c("CD107a", "CD109", "CD126", "CD138", "CD141", "CD144", "CD207", "CD324", "CD340", "CD49b", "CD57", "CD62E", "CD66b", "EGFR", "IgG.Fc", "LT.bR", "TCR.Vb13.1", "VEGFR.3", "XCR1", "CD10", "CD82", "IgM", "CD106", "CD305", "CD47", "HLA.ABC", "CD19", "CD205", "CD230", "CD5", "CD63", "CD73", "CD13", "CD235ab", "CD29", "CD62L", "CD36", "CD8", "CD117", "CD55", "CD94", "CD271", "CD34", "CD24", "CD61", "CD81", "CD49a", "CD72", "CD11a", "CD46", "CD86", "CD9", "CD123", "CD39", "CD49d", "CD33", "CD43", "CD62P", "HLA.DR", "CD11b", "CD11c", "CD32", "CD244", "CD26", "CD35", "CD41", "CD99", "CLEC12A", "HLA.DR.DP.DQ", "CD4")

peak_sep_summary %>% dplyr::filter(sep_power > 2) %>%
  ggplot(aes(x = Capture, y = sep_power, fill = Donor)) +
  geom_half_boxplot(outlier.shape = NA) +
  geom_half_point(size = 1) +
  theme_bw(base_size = 25) +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(4, "Set1"))(4)) +
  ylab("Stain Quality Score") +
  xlab("") +
  theme(legend.position = "top") +
  # rotate_x_text(angle = 45) +
  scale_y_continuous(trans = "log10")
ggsave(paste0(fig_path, "titrated_stain_quality_score_acrossCapture_filter.pdf"), width = 8, height = 9)

adt_marker_select = peak_sep_summary %>% dplyr::filter(sep_power > 2) %$% adt_marker %>% table %>% sort %>% data.frame %>% dplyr::filter(Freq > 0) %$% .
adtnorm_select_marker = c("CD107a", "CD109", "CD126", "CD138", "CD141", "CD144", "CD207", "CD324", "CD340", "CD49b", "CD57", "CD62E", "CD66b", "EGFR", "IgG.Fc", "LT.bR", "TCR.Vb13.1", "VEGFR.3", "XCR1", "CD10", "CD82", "IgM", "CD106", "CD305", "CD47", "HLA.ABC", "CD19", "CD205", "CD230", "CD5", "CD63", "CD73", "CD13", "CD235ab", "CD29", "CD62L", "CD36", "CD8", "CD117", "CD55", "CD94", "CD271", "CD34", "CD24", "CD61", "CD81", "CD49a", "CD72", "CD11a", "CD46", "CD86", "CD9", "CD123", "CD39", "CD49d", "CD33", "CD43", "CD62P", "HLA.DR", "CD11b", "CD11c", "CD32", "CD244", "CD26", "CD35", "CD41", "CD99", "CLEC12A", "HLA.DR.DP.DQ", "CD4")
adtnorm_select_marker = adtnorm_select_marker[which(adtnorm_select_marker %in% adt_marker_filter)]
sum((adt_marker_select %in% adtnorm_select_marker))/length(adt_marker_select) * 100
sum((adt_marker_select %in% adtnorm_select_marker))/length(adtnorm_select_marker) * 100

adt_marker_select[which(!(adt_marker_select %in% adtnorm_select_marker))]

adtnorm_select_marker[which(!(adtnorm_select_marker %in% adt_marker_select))]

peak_sep_summary$label = "Others"
peak_sep_summary$label[which(peak_sep_summary$adt_marker %in% adtnorm_select_marker)] = "ADTnorm Recommended"


peak_sep_summary %>% 
    data.frame %>% 
    dplyr::select(adt_marker, Capture, label, sep_power) %>% 
    dplyr::group_by(adt_marker, Capture, label) %>% 
    summarize(sep_power_mean = median(sep_power)) %>% 
    data.frame %>% ggplot(aes(x = label, y = sep_power_mean, fill = Capture)) +
    # geom_half_boxplot(outlier.shape = NA) +
    # geom_half_point(size = 0.5) +
    geom_boxplot() +
theme_bw(base_size = 25) +
scale_fill_brewer(palette = "Set2") +
ylab("Protein Marker Stain Quality Score") +
xlab("") +
theme(legend.position = "top") +
scale_y_continuous(trans = "log10")

ggsave(paste0(fig_path, "titrated_stain_quality_score_ADTnorm_recommended_boxplot.pdf"), width = 6, height = 9)

## ==============================================================================
## Differential detection using proportion of positive cells
## ==============================================================================
# pos_x_adt = matrix(-1, nrow = nrow(cell_x_adt), ncol = ncol(cell_x_adt))
# rownames(pos_x_adt) = rownames(cell_x_adt)
# colnames(pos_x_adt) = colnames(cell_x_adt)
# for(adt_marker_select in colnames(cell_x_adt)){ ##
    
#     peak_valley_density = readRDS(paste0(out_path, "/RDS/peak_valley_locations_", adt_marker_select, "_", run_name, ".rds"))
    
#     # peak_info = peak_valley_density$peak_landmark_list
#     valley_info = peak_valley_density$valley_landmark_list

#     for(each_sample in unique(adt_feature$sample)){

#         each_valley_info = valley_info[each_sample, 1]
#         ind_sample = which(adt_feature$sample == each_sample)
#         pos_x_adt[ind_sample, adt_marker_select][which(cell_x_adt[ind_sample, adt_marker_select] > each_valley_info)] = 1
#     }
# }
# saveRDS(pos_x_adt, file = paste0(out_path, "/RDS/pos_x_adt_", run_name, ".rds"))

pos_x_adt = readRDS(paste0(out_path, "/RDS/pos_x_adt_", run_name, ".rds"))



## Level 3 cell type
combined_data <- cbind(pos_x_adt[, !grepl("Isotype", colnames(pos_x_adt))], adt_feature)
pos_x_adt_df <- as.data.frame(pos_x_adt)
pos_x_adt_df$donor <- adt_feature$donor
pos_x_adt_df$Level3 <- adt_feature$`Level 3`
proportion_data <- pos_x_adt_df %>%
     group_by(donor, Level3) %>%
     summarise(across(everything(), ~ sum(. == 1, na.rm = TRUE) / n() * 100))
proportion_data
proportion_data = proportion_data %>% dplyr::filter(!is.na(Level3))

pos_data <- pos_x_adt_df %>%
     group_by(donor, Level3) %>%
     summarise(across(everything(), ~ sum(. == 1, na.rm = TRUE)))
pos_data

# Define columns
markers <- setdiff(colnames(proportion_data), c("Level3","donor"))

# Prepare empty list to store results
all_results <- list()
# Loop through each cell type in Level3
for (ct in unique(proportion_data$Level3)) {
  # Subset cell type (group1) and the rest (group2)
  group1 <- proportion_data[which(proportion_data$Level3 == ct), ]
  group2 <- proportion_data[which(proportion_data$Level3 != ct), ]
  
  # Store p-values for this cell type
  p_values <- sapply(markers, function(m) {
    test_res <- wilcox.test(group1[[m]], group2[[m]])
    test_res$p.value
  })
  
  # Adjust p-values
  p_adj <- p.adjust(p_values, method = "BH")
  
  # Filter significant markers
  sig_markers <- names(p_adj)[which(p_adj < 0.05)]
  
  if(length(sig_markers) > 0){
    # Save to result list
    all_results[[ct]] <- data.frame(
        cell_type = ct,
        marker = sig_markers,
        p_value = p_values[sig_markers],
        p_adj = p_adj[sig_markers]
    )

  }
}

# Combine all results for all cell types
final_results <- do.call(rbind, all_results)

# Output final results
final_results


# Define the custom order for cell types
custom_celltype_order <- c("HSC-1", "HSC-2", "MPP-1", "MPP-2", "MPP-MEP", "LMPP-1", "LMPP-1-cycling", "Multilin-1", "Multilin-2", "MEP-1", "MEP-2", "MEP-Eryth-1", "MEP-Eryth-2", "ERP-1", "ERP-2", "ERP-3", "ERP-4", "ERP-5", "ERP-6", "ERP-7", "ERP-8", "Erythroblast-1", "Erythroblast-2", "Erythroblast-3", "MK-Platelet", "MKP-early", "MKP-late", "BMCP-1", "BMCP-2",  "LMPP-2", "CLP", "Pro-B-Early", "Pro-B-Early-cycling", "Pro-B-cycling-1", "Pro-B-cycling-2", "Transitional-B-1", "Transitional-B-2", "Pro-B-1", "Pro-B-2", "Pro-B-3", "pre-B", "B Memory-1", "B Memory-2", "Plasma Cell", "T CD4 Naive-1", "T CD4 Naive-2", "T CD4 Activated", "CD4 TCM", "T CD8 Naive", "CD8 TEM", "MAIT", "NK-Mature", "MultiLin-GMP-1", "MultiLin-GMP-2", "MultiLin-GMP-3", "preNeu", "immNeu-1", "immNeu-2", "cMOP", "Mono-1", "Mono-2", "Intermediate Mono-1", "Intermediate Mono-2", "Intermediate Mono-3", "Classical-Mono", "Non-Classical Mono-1", "Non-Classical Mono-2", "Myeloid intermediate 1", "Myeloid intermediate 2", "Myeloid intermediate 3", "MDP-1", "MDP-2", "MDP-3", "MDP-4", "MDP-5", "MDP-6", "pre-DC-1", "pre-DC-2", "pre-DC-3", "cDC1", "cDC2-1", "cDC2-2", "pDC", "ASDC", "Mac", "MSC Fibroblast-1", "MSC Fibroblast-2", "Osteoblast", "Stromal Vascular")  # Example custom order


selected_marker = final_results %>% dplyr::filter(cell_type %in% custom_celltype_order) %$% marker %>% unique
selected_marker


# Load necessary libraries
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# Reshape the data into matrix format
heatmap_data <- proportion_data  %>% dplyr::filter(Level3 %in% custom_celltype_order) %>% #%>% dplyr::select(donor, Level3, all_of(selected_marker))
  pivot_longer(
    cols = -c(donor, Level3),    # Keep donor and Level3 as identifiers
    names_to = "ProteinMarker",  # Name for protein marker column
    values_to = "PositiveProp"          # Name for PositiveProp column
  ) %>%
  pivot_wider(
    names_from = c(donor, Level3), # Use donor and Level3 as columns
    values_from = PositiveProp            # Values are PositiveProps
  )

# Extract matrix data
protein_markers <- heatmap_data$ProteinMarker     # Save the protein markers
heatmap_matrix <- as.matrix(heatmap_data[, -1])   # Convert all but the first column to a matrix
rownames(heatmap_matrix) <- protein_markers       # Assign rownames from protein_markers

# Extract donor and cell type information for the columns
column_info <- strsplit(colnames(heatmap_matrix), "_")
column_info <- do.call(rbind, column_info)        # Convert list to matrix
colnames(column_info) <- c("Donor", "CellType")   # Assign column names
column_info <- as.data.frame(column_info)

# Define the custom order for cell types
column_info$CellType <- factor(column_info$CellType, levels = custom_celltype_order)

# Define colors for cell types and donors
celltype_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(custom_celltype_order)),
  custom_celltype_order
)

donor_colors <- setNames(
  RColorBrewer::brewer.pal(n = length(unique(column_info$Donor)), "Set1"),
  unique(column_info$Donor)
)

# Create annotation for the columns
col_annotation <- HeatmapAnnotation(
  CellType = column_info$CellType,
  Donor = column_info$Donor,
  col = list(
    CellType = celltype_colors,
    Donor = donor_colors
  ),
  annotation_name_side = "left"
)

pdf(paste0(fig_path, "titrated_heatmap_allcelltypeLevel3_allmarker.pdf"), width = 33, height = 40)
# Plot the heatmap
Heatmap(
  heatmap_matrix,
  name = "PositiveProp",                # Legend title
  cluster_rows = TRUE,           # Cluster rows (protein markers)
  cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
  top_annotation = col_annotation, # Add column annotations
  column_order = order(column_info$CellType), # Order columns by custom cell type
  column_title_rot = 90,
  column_split = column_info$CellType, # Split columns by cell type
  row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
  column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
#   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
  col = colorRamp2(c(min(heatmap_matrix, na.rm = TRUE), 
                     median(heatmap_matrix, na.rm = TRUE), 
                     max(heatmap_matrix, na.rm = TRUE)),
                   c("purple", "black", "yellow"))  # Color gradient
)
dev.off()

## aggregate donors

# Reshape the data into matrix format
heatmap_data <- proportion_data  %>% dplyr::filter(Level3 %in% custom_celltype_order) %>% group_by(Level3) %>% summarise(across(everything(), ~ mean(.))) %>%
#%>% dplyr::select(donor, Level3, all_of(selected_marker))
  pivot_longer(
    cols = -c(Level3),    # Keep donor and Level3 as identifiers
    names_to = "ProteinMarker",  # Name for protein marker column
    values_to = "PositiveProp"          # Name for PositiveProp column
  ) %>%
  pivot_wider(
    names_from = c(Level3), # Use donor and Level3 as columns
    values_from = PositiveProp            # Values are PositiveProps
  )

# Extract matrix data
protein_markers <- heatmap_data$ProteinMarker[-1]     # Save the protein markers
heatmap_matrix <- as.matrix(heatmap_data[-1, -1])   # Convert all but the first column to a matrix
rownames(heatmap_matrix) <- protein_markers       # Assign rownames from protein_markers

# Extract donor and cell type information for the columns
column_info <- strsplit(colnames(heatmap_matrix), "_")
column_info <- do.call(rbind, column_info)        # Convert list to matrix
colnames(column_info) <- c("CellType")   # Assign column names
column_info <- as.data.frame(column_info)

# Define the custom order for cell types
column_info$CellType <- factor(column_info$CellType, levels = custom_celltype_order)

# Define colors for cell types and donors
celltype_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(custom_celltype_order)),
  custom_celltype_order
)

# Create annotation for the columns
col_annotation <- HeatmapAnnotation(
  CellType = column_info$CellType,
  col = list(
    CellType = celltype_colors),
  annotation_name_side = "left"
)

pdf(paste0(fig_path, "titrated_heatmap_allcelltypeLevel3_allmarker_aggreDonor.pdf"), width = 33, height = 40)
# Plot the heatmap
Heatmap(
  heatmap_matrix,
  name = "PositiveProp",                # Legend title
  cluster_rows = TRUE,           # Cluster rows (protein markers)
  cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
  top_annotation = col_annotation, # Add column annotations
  column_order = order(column_info$CellType), # Order columns by custom cell type
  column_title_rot = 90,
  gap = unit(0, "cm"),
  column_gap = unit(0, "mm"),       # Remove the white gap between columns
  column_split = column_info$CellType, # Split columns by cell type
  row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
  column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
#   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
  col = colorRamp2(c(min(heatmap_matrix, na.rm = TRUE), 
                     median(heatmap_matrix, na.rm = TRUE), 
                     max(heatmap_matrix, na.rm = TRUE)),
                   c("purple", "black", "yellow"))  # Color gradient
)
dev.off()

heatmap_matrix %>% round(1) %>% write.csv(paste0(out_path, "/titrated_heatmap_allcelltypeLevel3_allmarker_aggreDonor_ADTnorm.csv"), row.names = TRUE)

## ====
# Create the heatmap object for the totalVI data
heatmap_obj <- Heatmap(
  heatmap_matrix,
  name = "PositiveProp",                # Legend title
  cluster_rows = TRUE,           # Cluster rows (protein markers)
  cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
  top_annotation = col_annotation, # Add column annotations
  column_order = order(column_info$CellType), # Order columns by custom cell type
  column_title_rot = 90,
  gap = unit(0, "cm"),
  column_gap = unit(0, "mm"),       # Remove the white gap between columns
  column_split = column_info$CellType, # Split columns by cell type
  row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
  column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
#   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
  col = colorRamp2(c(min(heatmap_matrix, na.rm = TRUE), 
                     median(heatmap_matrix, na.rm = TRUE), 
                     max(heatmap_matrix, na.rm = TRUE)),
                   c("purple", "black", "yellow"))  # Color gradient
)

# Draw the heatmap to get the clustering information
heatmap_drawn <- draw(heatmap_obj)
# Extract the clustered row order
clustered_row_order <- row_order(heatmap_drawn)
# Get the row names in the clustered order
clustered_row_names <- rownames(heatmap_matrix)[clustered_row_order]
# Output the clustered row names
clustered_row_names

## function!!
level3_heatmap = function(input_data, selected_marker, clustered_row_names, fig_path, file_suffix, low_t = -2, mid_t = 0, high_t = 2){

  # Reshape the data into matrix format
  heatmap_data <- input_data %>% dplyr::select(-donor) %>%
  dplyr::group_by(Level3M) %>% 
  dplyr::summarise(across(everything(), ~ mean(.))) %>% 
  dplyr::select(Level3M, all_of(selected_marker)) %>%
    pivot_longer(
      cols = -c(Level3M),    # Keep donor and Level3M as identifiers
      names_to = "ProteinMarker",  # Name for protein marker column
      values_to = "PositiveProp"          # Name for PositiveProp column
    ) %>%
    pivot_wider(
      names_from = c(Level3M), # Use donor and Level3M as columns
      values_from = PositiveProp            # Values are PositiveProps
    )


  # Extract matrix data
  protein_markers <- heatmap_data$ProteinMarker     # Save the protein markers
  heatmap_matrix <- as.matrix(heatmap_data[, -1])   # Convert all but the first column to a matrix
  rownames(heatmap_matrix) <- protein_markers       # Assign rownames from protein_markers

  # Extract donor and cell type information for the columns
  column_info <- strsplit(colnames(heatmap_matrix), "_")
  column_info <- do.call(rbind, column_info)        # Convert list to matrix
  colnames(column_info) <- c("CellType")   # Assign column names
  column_info <- as.data.frame(column_info)

  # Define the custom order for cell types
  column_info$CellType <- factor(column_info$CellType, levels = custom_celltype_order)

  # Define colors for cell types and donors
  celltype_colors <- setNames(
    colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(custom_celltype_order)),
    custom_celltype_order
  )

  # Create annotation for the columns
  col_annotation <- HeatmapAnnotation(
    CellType = column_info$CellType,
    col = list(
      CellType = celltype_colors),
    annotation_name_side = "left"
  )
  heatmap_matrix %>% round(1) %>% write.csv(paste0(out_path, "/titrated_heatmap_allcelltypeLevel3_allmarker_nonscaled_", file_suffix, ".csv"),  row.names = TRUE)

  pdf(paste0(fig_path, "titrated_heatmap_allcelltypeLevel3_allmarker_aggreDonor_", file_suffix, "_nonscaled.pdf"), width = 33, height = 40)
  # Plot the heatmap
  print(Heatmap(
    heatmap_matrix[clustered_row_names, ],
    name = "PositiveProp",                # Legend title
    cluster_rows = TRUE,           # Cluster rows (protein markers)
    cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
    top_annotation = col_annotation, # Add column annotations
    column_order = order(column_info$CellType), # Order columns by custom cell type
    column_title_rot = 90,
    gap = unit(0, "cm"),
    column_gap = unit(0, "mm"),       # Remove the white gap between columns
    column_split = column_info$CellType, # Split columns by cell type
    row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
    column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
  #   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
    col = colorRamp2(c(min(heatmap_matrix, na.rm = TRUE), 
                      median(heatmap_matrix, na.rm = TRUE), 
                      max(heatmap_matrix, na.rm = TRUE)),
                    c("purple", "black", "yellow"))  # Color gradient
  ))
  print(Heatmap(
    heatmap_matrix[clustered_row_names, ],
    name = "PositiveProp",                # Legend title
    cluster_rows = FALSE,           # Cluster rows (protein markers)
    cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
    top_annotation = col_annotation, # Add column annotations
    column_order = order(column_info$CellType), # Order columns by custom cell type
    column_title_rot = 90,
    gap = unit(0, "cm"),
    column_gap = unit(0, "mm"),       # Remove the white gap between columns
    column_split = column_info$CellType, # Split columns by cell type
    row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
    column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
  #   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
    col = colorRamp2(c(0, 3, 6),
                    c("purple", "black", "yellow"))  # Color gradient
  ))
  print(Heatmap(
    heatmap_matrix[clustered_row_names, ],
    name = "PositiveProp",                # Legend title
    cluster_rows = TRUE,           # Cluster rows (protein markers)
    cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
    top_annotation = col_annotation, # Add column annotations
    column_order = order(column_info$CellType), # Order columns by custom cell type
    column_title_rot = 90,
    gap = unit(0, "cm"),
    column_gap = unit(0, "mm"),       # Remove the white gap between columns
    column_split = column_info$CellType, # Split columns by cell type
    row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
    column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
  #   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
    col = colorRamp2(c(0, 3, 6),
                    c("purple", "black", "yellow"))  # Color gradient
  ))
  dev.off()
  # Function to scale rows to [-1, 1]
  # scale_to_range <- function(x) {
  #   (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * 6 - 3
  # }
  scale_to_range <- function(y) {
    x = log2(y/mean(y, na.rm = TRUE))
    return(x)
    # (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * 6 - 3
  }
  # Apply scaling function to each row
  scaled_heatmap_matrix <- t(apply(heatmap_matrix[clustered_row_names, ], 1, scale_to_range))
  scaled_heatmap_matrix %>% round(1) %>% write.csv(paste0(out_path, "/titrated_heatmap_allcelltypeLevel3_allmarker_log2FC_byMean_", file_suffix, ".csv"),  row.names = TRUE)

  pdf(paste0(fig_path, "titrated_heatmap_allcelltypeLevel3_allmarker_aggreDonor_", file_suffix, "_log2FC_byMean.pdf"), width = 33, height = 40)
  # Plot the heatmap
  print(Heatmap(
    scaled_heatmap_matrix[clustered_row_names, ],
    name = "PositiveProp",                # Legend title
    cluster_rows = FALSE,           # Cluster rows (protein markers)
    cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
    top_annotation = col_annotation, # Add column annotations
    column_order = order(column_info$CellType), # Order columns by custom cell type
    column_title_rot = 90,
    gap = unit(0, "cm"),
    column_gap = unit(0, "mm"),       # Remove the white gap between columns
    column_split = column_info$CellType, # Split columns by cell type
    row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
    column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
  #   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
    col = colorRamp2(c(-1, 0, 1),
                    c("purple", "black", "yellow"))  # Color gradient
  ))
  print(Heatmap(
    scaled_heatmap_matrix[clustered_row_names, ],
    name = "PositiveProp",                # Legend title
    cluster_rows = TRUE,           # Cluster rows (protein markers)
    cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
    top_annotation = col_annotation, # Add column annotations
    column_order = order(column_info$CellType), # Order columns by custom cell type
    column_title_rot = 90,
    gap = unit(0, "cm"),
    column_gap = unit(0, "mm"),       # Remove the white gap between columns
    column_split = column_info$CellType, # Split columns by cell type
    row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
    column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
  #   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
    col = colorRamp2(c(-1, 0, 1),
                    c("purple", "black", "yellow"))  # Color gradient
  ))
  dev.off()

  scale_to_range <- function(y) {
    x = log2(y/median(y, na.rm = TRUE))
    return(x)
    # (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * 6 - 3
  }
  # Apply scaling function to each row
  scaled_heatmap_matrix <- t(apply(heatmap_matrix[clustered_row_names, ], 1, scale_to_range))

  pdf(paste0(fig_path, "titrated_heatmap_allcelltypeLevel3_allmarker_aggreDonor_", file_suffix, "_log2FC_byMedian.pdf"), width = 33, height = 40)
  # Plot the heatmap
  print(Heatmap(
    scaled_heatmap_matrix[clustered_row_names, ],
    name = "PositiveProp",                # Legend title
    cluster_rows = FALSE,           # Cluster rows (protein markers)
    cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
    top_annotation = col_annotation, # Add column annotations
    column_order = order(column_info$CellType), # Order columns by custom cell type
    column_title_rot = 90,
    gap = unit(0, "cm"),
    column_gap = unit(0, "mm"),       # Remove the white gap between columns
    column_split = column_info$CellType, # Split columns by cell type
    row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
    column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
  #   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
    col = colorRamp2(c(-1, 0, 1),
                    c("purple", "black", "yellow"))  # Color gradient
  ))
  print(Heatmap(
    scaled_heatmap_matrix[clustered_row_names, ],
    name = "PositiveProp",                # Legend title
    cluster_rows = TRUE,           # Cluster rows (protein markers)
    cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
    top_annotation = col_annotation, # Add column annotations
    column_order = order(column_info$CellType), # Order columns by custom cell type
    column_title_rot = 90,
    gap = unit(0, "cm"),
    column_gap = unit(0, "mm"),       # Remove the white gap between columns
    column_split = column_info$CellType, # Split columns by cell type
    row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
    column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
  #   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
    col = colorRamp2(c(-1, 0, 1),
                    c("purple", "black", "yellow"))  # Color gradient
  ))
  dev.off()

}

totalvi_data = totalVI_input[colnames(seurat_obj), ] %>% data.frame
colnames(totalvi_data) = rownames(adt_obj)
totalvi_data$donor = adt_obj@meta.data$donor
totalvi_data$Level3M = adt_obj@meta.data$Level3M
totalvi_data$Level3M[which(totalvi_data$Level3M == "BMCP")] = "BMCP-1"
totalvi_data$Level3M[which(totalvi_data$Level3M == "Eosinophil-early")] = "BMCP-2"

# unique(totalvi_data$Level3M)[which(!(totalvi_data$Level3M %>% unique %in% custom_celltype_order))]

# custom_celltype_order[which(!(custom_celltype_order %in% unique(totalvi_data$Level3M)))]

selected_marker = rownames(adt_obj)[!grepl("Isotype", rownames(adt_obj))]

level3_heatmap(totalvi_data, selected_marker, clustered_row_names, fig_path, "totalVI_adjust")


clr_data = seurat_obj@assays$ADT@scale.data %>% t %>% data.frame
colnames(clr_data) = rownames(adt_obj)
clr_data$donor = adt_obj@meta.data$donor
clr_data$Level3M = adt_obj@meta.data$Level3M
clr_data$Level3M[which(clr_data$Level3M == "BMCP")] = "BMCP-1"
clr_data$Level3M[which(clr_data$Level3M == "Eosinophil-early")] = "BMCP-2"
level3_heatmap(clr_data, clustered_row_names, clustered_row_names, fig_path, "CLR_adjust")


harmony_data = out %>% data.frame
harmony_data$donor = adt_feature$donor
harmony_data$Level3M = adt_feature$`Level 3`
level3_heatmap(harmony_data, clustered_row_names, clustered_row_names, fig_path, "Harmony_adjust")






## ================
## Level 3 selected a few marker to make figure smaller
## ================
combined_data <- cbind(pos_x_adt[, !grepl("Isotype", colnames(pos_x_adt))], adt_feature)
pos_x_adt_df <- as.data.frame(pos_x_adt)
pos_x_adt_df$donor <- adt_feature$donor
pos_x_adt_df$Level3 <- adt_feature$`Level 3`
proportion_data <- pos_x_adt_df %>%
     group_by(donor, Level3) %>%
     summarise(across(everything(), ~ sum(. == 1, na.rm = TRUE) / n() * 100))
proportion_data
proportion_data = proportion_data %>% dplyr::filter(!is.na(Level3))

pos_data <- pos_x_adt_df %>%
     group_by(donor, Level3) %>%
     summarise(across(everything(), ~ sum(. == 1, na.rm = TRUE)))
pos_data

# Define columns
markers <- setdiff(colnames(proportion_data), c("Level3","donor"))

# Prepare empty list to store results
all_results <- list()
# Loop through each cell type in Level3
for (ct in unique(proportion_data$Level3)) {
  # Subset cell type (group1) and the rest (group2)
  group1 <- proportion_data[which(proportion_data$Level3 == ct), ]
  group2 <- proportion_data[which(proportion_data$Level3 != ct), ]
  
  # Store p-values for this cell type
  p_values <- sapply(markers, function(m) {
    test_res <- wilcox.test(group1[[m]], group2[[m]])
    test_res$p.value
  })
  
  # Adjust p-values
  p_adj <- p.adjust(p_values, method = "BH")
  
  # Filter significant markers
  sig_markers <- names(p_adj)[which(p_adj < 0.05)]
  
  if(length(sig_markers) > 0){
    # Save to result list
    all_results[[ct]] <- data.frame(
        cell_type = ct,
        marker = sig_markers,
        p_value = p_values[sig_markers],
        p_adj = p_adj[sig_markers]
    )

  }
}

# Combine all results for all cell types
final_results <- do.call(rbind, all_results)

# Output final results
final_results


# Define the custom order for cell types
custom_celltype_order <- c("HSC-1", "HSC-2", "MPP-1", "MPP-2", "MPP-MEP", "LMPP-1", "LMPP-1-cycling", "Multilin-1", "Multilin-2", "MEP-1", "MEP-2", "MEP-Eryth-1", "MEP-Eryth-2", "ERP-1", "ERP-2", "ERP-3", "ERP-4", "ERP-5", "ERP-6", "ERP-7", "ERP-8", "Erythroblast-1", "Erythroblast-2", "Erythroblast-3", "MK-Platelet", "MKP-early", "MKP-late", "BMCP-1", "BMCP-2",  "LMPP-2", "CLP", "Pro-B-Early", "Pro-B-Early-cycling", "Pro-B-cycling-1", "Pro-B-cycling-2", "Transitional-B-1", "Transitional-B-2", "Pro-B-1", "Pro-B-2", "Pro-B-3", "pre-B", "B Memory-1", "B Memory-2", "Plasma Cell", "T CD4 Naive-1", "T CD4 Naive-2", "T CD4 Activated", "CD4 TCM", "T CD8 Naive", "CD8 TEM", "MAIT", "NK-Mature", "MultiLin-GMP-1", "MultiLin-GMP-2", "MultiLin-GMP-3", "preNeu", "immNeu-1", "immNeu-2", "cMOP", "Mono-1", "Mono-2", "Intermediate Mono-1", "Intermediate Mono-2", "Intermediate Mono-3", "Classical-Mono", "Non-Classical Mono-1", "Non-Classical Mono-2", "Myeloid intermediate 1", "Myeloid intermediate 2", "Myeloid intermediate 3", "MDP-1", "MDP-2", "MDP-3", "MDP-4", "MDP-5", "MDP-6", "pre-DC-1", "pre-DC-2", "pre-DC-3", "cDC1", "cDC2-1", "cDC2-2", "pDC", "ASDC", "Mac", "MSC Fibroblast-1", "MSC Fibroblast-2", "Osteoblast", "Stromal Vascular")  # Example custom order


selected_marker = final_results %>% dplyr::filter(cell_type %in% custom_celltype_order) %$% marker %>% unique
selected_marker


# Reshape the data into matrix format
heatmap_data <- proportion_data  %>% dplyr::filter(Level3 %in% custom_celltype_order) %>% #%>% dplyr::select(donor, Level3, all_of(selected_marker))
  pivot_longer(
    cols = -c(donor, Level3),    # Keep donor and Level3 as identifiers
    names_to = "ProteinMarker",  # Name for protein marker column
    values_to = "PositiveProp"          # Name for PositiveProp column
  ) %>%
  pivot_wider(
    names_from = c(donor, Level3), # Use donor and Level3 as columns
    values_from = PositiveProp            # Values are PositiveProps
  )

# Extract matrix data
protein_markers <- heatmap_data$ProteinMarker     # Save the protein markers
heatmap_matrix <- as.matrix(heatmap_data[, -1])   # Convert all but the first column to a matrix
rownames(heatmap_matrix) <- protein_markers       # Assign rownames from protein_markers

# Extract donor and cell type information for the columns
column_info <- strsplit(colnames(heatmap_matrix), "_")
column_info <- do.call(rbind, column_info)        # Convert list to matrix
colnames(column_info) <- c("Donor", "CellType")   # Assign column names
column_info <- as.data.frame(column_info)

# Define the custom order for cell types
column_info$CellType <- factor(column_info$CellType, levels = custom_celltype_order)

# Define colors for cell types and donors
celltype_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(custom_celltype_order)),
  custom_celltype_order
)

donor_colors <- setNames(
  RColorBrewer::brewer.pal(n = length(unique(column_info$Donor)), "Set1"),
  unique(column_info$Donor)
)

# Create annotation for the columns
col_annotation <- HeatmapAnnotation(
  CellType = column_info$CellType,
  Donor = column_info$Donor,
  col = list(
    CellType = celltype_colors,
    Donor = donor_colors
  ),
  annotation_name_side = "left"
)

## aggregate donors

# Reshape the data into matrix format
heatmap_data <- proportion_data  %>% dplyr::filter(Level3 %in% custom_celltype_order) %>% group_by(Level3) %>% summarise(across(everything(), ~ mean(.))) %>%
#%>% dplyr::select(donor, Level3, all_of(selected_marker))
  pivot_longer(
    cols = -c(Level3),    # Keep donor and Level3 as identifiers
    names_to = "ProteinMarker",  # Name for protein marker column
    values_to = "PositiveProp"          # Name for PositiveProp column
  ) %>%
  pivot_wider(
    names_from = c(Level3), # Use donor and Level3 as columns
    values_from = PositiveProp            # Values are PositiveProps
  )

# Extract matrix data
protein_markers <- heatmap_data$ProteinMarker[-1]     # Save the protein markers
heatmap_matrix <- as.matrix(heatmap_data[-1, -1])   # Convert all but the first column to a matrix
rownames(heatmap_matrix) <- protein_markers       # Assign rownames from protein_markers

# Extract donor and cell type information for the columns
column_info <- strsplit(colnames(heatmap_matrix), "_")
column_info <- do.call(rbind, column_info)        # Convert list to matrix
colnames(column_info) <- c("CellType")   # Assign column names
column_info <- as.data.frame(column_info)

# Define the custom order for cell types
column_info$CellType <- factor(column_info$CellType, levels = custom_celltype_order)

# Define colors for cell types and donors
celltype_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(custom_celltype_order)),
  custom_celltype_order
)

# Create annotation for the columns
col_annotation <- HeatmapAnnotation(
  CellType = column_info$CellType,
  col = list(
    CellType = celltype_colors),
  annotation_name_side = "left"
)

short_list_pro = c("CD44", "CD162", "CD11a", "CD305", "HLA.DR.DP.DQ", "CD33", 
"CD18", "CD62L", "CD45RA", "CD71", "CD43", "CD164", "CD82", "CD117", 
"CD205", "CD133", "CD36", "CD63", "CD41", "CD35", "CD54", "FR.b", 
"CD45", "CD32", "CD11c", "CD101", "CD55", "CD58", "CD116", "CD64", 
"CD354", "CD172a", "CD11b", "CD93", "CD192", "CD38", "CD9", "CD47", 
"CD81", "CD10", "CD72", "CD19", "CD24", "CD2", "CD7", "CD5", 
"CD52", "CD123", "CD84", "CD49f", "CD226", "CD102", "CLEC1B", 
"CD49b", "CD42b", "CD61", "CD62P", "CD151", "CD29", "HLA.ABC", 
"CD13", "CD73", "CD106", "CD304", "CD90", "CD271", "CD27", "CD45RB", 
"CD8", "CD26", "CD4", "CD45RO", "CD127", "CD34", "CD14", "CD235a", "CD16", "CD163")

pdf(paste0(fig_path, "titrated_heatmap_allcelltypeLevel3_allmarker_aggreDonor_short_list.pdf"), width = 18, height = 16)
# Plot the heatmap
Heatmap(
  heatmap_matrix[short_list_pro, ],
  name = "PositiveProp",                # Legend title
  cluster_rows = TRUE,           # Cluster rows (protein markers)
  cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
  top_annotation = col_annotation, # Add column annotations
  column_order = order(column_info$CellType), # Order columns by custom cell type
  column_title_rot = 90,
  gap = unit(0, "cm"),
  column_gap = unit(0, "mm"),       # Remove the white gap between columns
  column_split = column_info$CellType, # Split columns by cell type
  row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
  column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
#   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
  col = colorRamp2(c(min(heatmap_matrix, na.rm = TRUE), 
                     median(heatmap_matrix, na.rm = TRUE), 
                     max(heatmap_matrix, na.rm = TRUE)),
                   c("purple", "black", "yellow"))  # Color gradient
)
dev.off()

## ====
# Create the heatmap object for the totalVI data
heatmap_obj <- Heatmap(
  heatmap_matrix[short_list_pro,],
  name = "PositiveProp",                # Legend title
  cluster_rows = TRUE,           # Cluster rows (protein markers)
  cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
  top_annotation = col_annotation, # Add column annotations
  column_order = order(column_info$CellType), # Order columns by custom cell type
  column_title_rot = 90,
  gap = unit(0, "cm"),
  column_gap = unit(0, "mm"),       # Remove the white gap between columns
  column_split = column_info$CellType, # Split columns by cell type
  row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
  column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
#   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
  col = colorRamp2(c(min(heatmap_matrix, na.rm = TRUE), 
                     median(heatmap_matrix, na.rm = TRUE), 
                     max(heatmap_matrix, na.rm = TRUE)),
                   c("purple", "black", "yellow"))  # Color gradient
)

# Draw the heatmap to get the clustering information
heatmap_drawn <- draw(heatmap_obj)
# Extract the clustered row order
clustered_row_order <- row_order(heatmap_drawn)
# Get the row names in the clustered order
clustered_row_names <- rownames(heatmap_matrix[short_list_pro, ])[clustered_row_order]
# Output the clustered row names
clustered_row_names



## function!!
level3_heatmap = function(input_data, selected_marker, clustered_row_names, fig_path, file_suffix, low_t = -2, mid_t = 0, high_t = 2){

  # Reshape the data into matrix format
  heatmap_data <- input_data %>% dplyr::select(-donor) %>%
  dplyr::group_by(Level3M) %>% 
  dplyr::summarise(across(everything(), ~ mean(.))) %>% 
  dplyr::select(Level3M, all_of(selected_marker)) %>%
    pivot_longer(
      cols = -c(Level3M),    # Keep donor and Level3M as identifiers
      names_to = "ProteinMarker",  # Name for protein marker column
      values_to = "PositiveProp"          # Name for PositiveProp column
    ) %>%
    pivot_wider(
      names_from = c(Level3M), # Use donor and Level3M as columns
      values_from = PositiveProp            # Values are PositiveProps
    )


  # Extract matrix data
  protein_markers <- heatmap_data$ProteinMarker     # Save the protein markers
  heatmap_matrix <- as.matrix(heatmap_data[, -1])   # Convert all but the first column to a matrix
  rownames(heatmap_matrix) <- protein_markers       # Assign rownames from protein_markers

  # Extract donor and cell type information for the columns
  column_info <- strsplit(colnames(heatmap_matrix), "_")
  column_info <- do.call(rbind, column_info)        # Convert list to matrix
  colnames(column_info) <- c("CellType")   # Assign column names
  column_info <- as.data.frame(column_info)

  # Define the custom order for cell types
  column_info$CellType <- factor(column_info$CellType, levels = custom_celltype_order)

  # Define colors for cell types and donors
  celltype_colors <- setNames(
    colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(custom_celltype_order)),
    custom_celltype_order
  )

  # Create annotation for the columns
  col_annotation <- HeatmapAnnotation(
    CellType = column_info$CellType,
    col = list(
      CellType = celltype_colors),
    annotation_name_side = "left"
  )
  # heatmap_matrix %>% round(1) %>% write.csv(paste0(out_path, "/titrated_heatmap_allcelltypeLevel3_allmarker_nonscaled_", file_suffix, ".csv"),  row.names = TRUE)

  pdf(paste0(fig_path, "titrated_heatmap_allcelltypeLevel3_allmarker_aggreDonor_", file_suffix, "_nonscaled.pdf"), width = 18, height = 16)
  # Plot the heatmap
  print(Heatmap(
    heatmap_matrix[clustered_row_names, ],
    name = "PositiveProp",                # Legend title
    cluster_rows = TRUE,           # Cluster rows (protein markers)
    cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
    top_annotation = col_annotation, # Add column annotations
    column_order = order(column_info$CellType), # Order columns by custom cell type
    column_title_rot = 90,
    gap = unit(0, "cm"),
    column_gap = unit(0, "mm"),       # Remove the white gap between columns
    column_split = column_info$CellType, # Split columns by cell type
    row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
    column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
  #   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
    col = colorRamp2(c(min(heatmap_matrix, na.rm = TRUE), 
                      median(heatmap_matrix, na.rm = TRUE), 
                      max(heatmap_matrix, na.rm = TRUE)),
                    c("purple", "black", "yellow"))  # Color gradient
  ))
  print(Heatmap(
    heatmap_matrix[clustered_row_names, ],
    name = "PositiveProp",                # Legend title
    cluster_rows = FALSE,           # Cluster rows (protein markers)
    cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
    top_annotation = col_annotation, # Add column annotations
    column_order = order(column_info$CellType), # Order columns by custom cell type
    column_title_rot = 90,
    gap = unit(0, "cm"),
    column_gap = unit(0, "mm"),       # Remove the white gap between columns
    column_split = column_info$CellType, # Split columns by cell type
    row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
    column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
  #   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
    col = colorRamp2(c(0, 3, 6),
                    c("purple", "black", "yellow"))  # Color gradient
  ))
  print(Heatmap(
    heatmap_matrix[clustered_row_names, ],
    name = "PositiveProp",                # Legend title
    cluster_rows = TRUE,           # Cluster rows (protein markers)
    cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
    top_annotation = col_annotation, # Add column annotations
    column_order = order(column_info$CellType), # Order columns by custom cell type
    column_title_rot = 90,
    gap = unit(0, "cm"),
    column_gap = unit(0, "mm"),       # Remove the white gap between columns
    column_split = column_info$CellType, # Split columns by cell type
    row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
    column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
  #   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
    col = colorRamp2(c(0, 3, 6),
                    c("purple", "black", "yellow"))  # Color gradient
  ))
  dev.off()
  # Function to scale rows to [-1, 1]
  # scale_to_range <- function(x) {
  #   (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * 6 - 3
  # }
  scale_to_range <- function(y) {
    x = log2(y/mean(y, na.rm = TRUE))
    return(x)
    # (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * 6 - 3
  }
  # Apply scaling function to each row
  scaled_heatmap_matrix <- t(apply(heatmap_matrix[clustered_row_names, ], 1, scale_to_range))
  # scaled_heatmap_matrix %>% round(1) %>% write.csv(paste0(out_path, "/titrated_heatmap_allcelltypeLevel3_allmarker_log2FC_byMean_", file_suffix, ".csv"),  row.names = TRUE)

  pdf(paste0(fig_path, "titrated_heatmap_allcelltypeLevel3_allmarker_aggreDonor_", file_suffix, "_log2FC_byMean.pdf"), width = 18, height = 16)
  # Plot the heatmap
  print(Heatmap(
    scaled_heatmap_matrix[clustered_row_names, ],
    name = "PositiveProp",                # Legend title
    cluster_rows = FALSE,           # Cluster rows (protein markers)
    cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
    top_annotation = col_annotation, # Add column annotations
    column_order = order(column_info$CellType), # Order columns by custom cell type
    column_title_rot = 90,
    gap = unit(0, "cm"),
    column_gap = unit(0, "mm"),       # Remove the white gap between columns
    column_split = column_info$CellType, # Split columns by cell type
    row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
    column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
  #   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
    col = colorRamp2(c(-2, 0, 2),
                    c("purple", "black", "yellow"))  # Color gradient
  ))
  print(Heatmap(
    scaled_heatmap_matrix[clustered_row_names, ],
    name = "PositiveProp",                # Legend title
    cluster_rows = TRUE,           # Cluster rows (protein markers)
    cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
    top_annotation = col_annotation, # Add column annotations
    column_order = order(column_info$CellType), # Order columns by custom cell type
    column_title_rot = 90,
    gap = unit(0, "cm"),
    column_gap = unit(0, "mm"),       # Remove the white gap between columns
    column_split = column_info$CellType, # Split columns by cell type
    row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
    column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
  #   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
    col = colorRamp2(c(-2, 0, 2),
                    c("purple", "black", "yellow"))  # Color gradient
  ))
  dev.off()

  scale_to_range <- function(y) {
    x = log2(y/median(y, na.rm = TRUE))
    return(x)
    # (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * 6 - 3
  }
  # Apply scaling function to each row
  scaled_heatmap_matrix <- t(apply(heatmap_matrix[clustered_row_names, ], 1, scale_to_range))

  pdf(paste0(fig_path, "titrated_heatmap_allcelltypeLevel3_allmarker_aggreDonor_", file_suffix, "_log2FC_byMedian.pdf"), width = 18, height = 16)
  # Plot the heatmap
  print(Heatmap(
    scaled_heatmap_matrix[clustered_row_names, ],
    name = "PositiveProp",                # Legend title
    cluster_rows = FALSE,           # Cluster rows (protein markers)
    cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
    top_annotation = col_annotation, # Add column annotations
    column_order = order(column_info$CellType), # Order columns by custom cell type
    column_title_rot = 90,
    gap = unit(0, "cm"),
    column_gap = unit(0, "mm"),       # Remove the white gap between columns
    column_split = column_info$CellType, # Split columns by cell type
    row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
    column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
  #   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
    col = colorRamp2(c(-2, 0, 2),
                    c("purple", "black", "yellow"))  # Color gradient
  ))
  print(Heatmap(
    scaled_heatmap_matrix[clustered_row_names, ],
    name = "PositiveProp",                # Legend title
    cluster_rows = TRUE,           # Cluster rows (protein markers)
    cluster_columns = FALSE,        # Cluster columns (donor and cell type combinations)
    top_annotation = col_annotation, # Add column annotations
    column_order = order(column_info$CellType), # Order columns by custom cell type
    column_title_rot = 90,
    gap = unit(0, "cm"),
    column_gap = unit(0, "mm"),       # Remove the white gap between columns
    column_split = column_info$CellType, # Split columns by cell type
    row_names_gp = gpar(fontsize = 16),   # Adjust font size for row names
    column_names_gp = gpar(fontsize = 0), # gpar(fontsize = 8), # Adjust font size for column names
  #   clustering_method_rows = "ward.D2", # Use Ward clustering for rows
    col = colorRamp2(c(-2, 0, 2),
                    c("purple", "black", "yellow"))  # Color gradient
  ))
  dev.off()

}

totalvi_data = totalVI_input[colnames(seurat_obj), ] %>% data.frame
colnames(totalvi_data) = rownames(adt_obj)
totalvi_data$donor = adt_obj@meta.data$donor
totalvi_data$Level3M = adt_obj@meta.data$Level3M
totalvi_data$Level3M[which(totalvi_data$Level3M == "BMCP")] = "BMCP-1"
totalvi_data$Level3M[which(totalvi_data$Level3M == "Eosinophil-early")] = "BMCP-2"

# unique(totalvi_data$Level3M)[which(!(totalvi_data$Level3M %>% unique %in% custom_celltype_order))]

# custom_celltype_order[which(!(custom_celltype_order %in% unique(totalvi_data$Level3M)))]

selected_marker = rownames(adt_obj)[!grepl("Isotype", rownames(adt_obj))]

level3_heatmap(totalvi_data, clustered_row_names, clustered_row_names, fig_path, "totalVI_shortlist")


clr_data = seurat_obj@assays$ADT@scale.data %>% t %>% data.frame
colnames(clr_data) = rownames(adt_obj)
clr_data$donor = adt_obj@meta.data$donor
clr_data$Level3M = adt_obj@meta.data$Level3M
clr_data$Level3M[which(clr_data$Level3M == "BMCP")] = "BMCP-1"
clr_data$Level3M[which(clr_data$Level3M == "Eosinophil-early")] = "BMCP-2"
level3_heatmap(clr_data, clustered_row_names, clustered_row_names, fig_path, "CLR_adjust")


harmony_data = out %>% data.frame
harmony_data$donor = adt_feature$donor
harmony_data$Level3M = adt_feature$`Level 3`
level3_heatmap(harmony_data, clustered_row_names, clustered_row_names, fig_path, "Harmony_adjust")

res_norm_data = res_norm %>% data.frame
res_norm_data$donor = adt_feature$donor
res_norm_data$Level3M = adt_feature$`Level 3`
level3_heatmap(res_norm_data, clustered_row_names, clustered_row_names, fig_path, "res_norm")
