#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(umap)
library(viridis)


method <- args[1] #"Arcsinh_b5"
run_name = args[2] #"publicData_CITEseq"
out_path = "./manuscript/results/public13Dataset_CITEseq" 
# adt_data = readRDS(file = paste0(out_path, "/RDS/adt_data_RawCount_", run_name, ".rds"))
# adt_feature = readRDS(file = paste0("./results/publicData_CITEseq/RDS/adt_feature_", run_name, ".rds"))
adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", run_name, ".rds"))

ans = readRDS(file = paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds"))
ans[is.na(ans)] = 0
ans_pca = prcomp(ans)$x 
print("Start UMAP")
ans_umap = umap(ans_pca[, 1:min(ncol(ans), 50)])$layout %>% data.frame
saveRDS(list(ans_pca, ans_umap), file = paste0(out_path, "/RDS/adt_", method, "_", run_name, "_pca_umap.rds"))

# tmp = readRDS(file = paste0(out_path, "/RDS/adt_", method, "_", run_name, "_pca_umap.rds"))
# ans_umap = tmp[[2]]


print("Start plotting")
pdf(paste0(out_path, "/Figures/adt_", method, "_", run_name, "_umap.pdf"), width = 15, height = 15)
reindex1 = which(adt_feature$cell_type_l1 == "undefined")
reindex2 = which(adt_feature$cell_type_l1 != "undefined")
reindex = c(reindex1, reindex2)

print(plot_umap(ans_umap[reindex, ], adt_feature[reindex, ], point_size = 0.3, parameter_list = list(
target_feature = "sample", 
color_design = colorRampPalette(brewer.pal(8, "Spectral"))(length(adt_feature$sample %>% unique)),
method_label = method
)) #+ scale_color_brewer(palette = "Dark2")
)

print(plot_umap(ans_umap[reindex, ], adt_feature[reindex, ], point_size = 0.3, parameter_list = list(
target_feature = "study_name", 
color_design = colorRampPalette(brewer.pal(8, "Dark2"))(length(adt_feature$study_name %>% unique)),
method_label = method
)))

print(plot_umap(ans_umap[reindex, ], adt_feature[reindex, ], point_size = 0.3, parameter_list = list(
target_feature = "sample_status", 
color_design = colorRampPalette(brewer.pal(8, "Set1"))(length(adt_feature$sample_status %>% unique)),
method_label = method
)))

adt_feature$cell_type_l1 = factor(adt_feature$cell_type_l1, levels = c("B", "CD4 T", "CD8 T", "monocytes", "NK", "DCs", "undefined"))
print(plot_umap(ans_umap[reindex, ], adt_feature[reindex, ], point_size = 0.3, parameter_list = list(
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
print(plot_umap(ans_umap[reindex, ], adt_feature[reindex, ], point_size = 0.3, parameter_list = list(
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
 
