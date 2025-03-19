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
library(rcartocolor)
library(pheatmap)
library(ADTnorm)

## =====================
## study name and paths
## =====================
run_name = "COVID19"
master_path = "./"
out_path = paste0(master_path, "manuscript/results/", run_name)
fig_path = paste0(out_path, "/Figures/", run_name)

adt_data = readRDS(file = paste0(out_path, "/RDS/adt_data_", run_name, ".rds"))
# rna_data = readRDS(file = paste0(out_path, "/RDS/rna_data_", run_name, ".rds")) 
adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", run_name, ".rds"))
rm_ind = readRDS(file = paste0(out_path, "/RDS/rm_ind1_ind2_", run_name, ".rds"))

adt_feature$CZI_initial_clustering = factor(adt_feature$CZI_initial_clustering, levels = c("CD4", "Treg", "CD8", "gdT", "MAIT", "NK_16hi",  "NK_56hi", "B_cell", "Plasmablast",  "CD14", "CD16", "DCs", "pDC", "HSC",  "Platelets", "RBC", "Lymph_prolif",  "Mono_prolif")
, labels = c("CD4 T", "Treg", "CD8 T", "gamma delta T", "MAIT", "NK CD56-",  "NK CD56+", "B", "Plasmablast",  "Monocytes CD14+", "Monocytes CD16+", "cDC", "pDC", "HSPC",  "Platelets", "RBC", "Prolif lymph",  "Prolif Mono"))

## ===========
## auto-gating
## ===========
cell_x_adt = arcsinh_transform(cell_x_adt = adt_data) 
cell_x_feature = adt_feature %>% data.frame
rownames(cell_x_adt) = paste0("Cell", 1:nrow(cell_x_adt))
rownames(cell_x_feature) = paste0("Cell", 1:nrow(cell_x_feature))
cell_x_feature$sample = cell_x_feature$site_sample


pheno_matrix <- matrix("-", nrow(cell_x_adt), ncol(cell_x_adt))
colnames(pheno_matrix) <- colnames(cell_x_adt)
rownames(pheno_matrix) <- rownames(cell_x_adt)
for(adt_marker_select in colnames(cell_x_adt)){
    if(file.exists(paste0(out_path, "/manualTuning/RDS/peak_valley_locations_", adt_marker_select, "_", run_name, ".rds"))){
        peak_valley_density = readRDS(paste0(out_path, "/manualTuning/RDS/peak_valley_locations_", adt_marker_select, "_", run_name, ".rds"))
    
    }else{
        if(file.exists(paste0(out_path, "/RDS/peak_valley_locations_", adt_marker_select, "_", run_name, ".rds"))){
            peak_valley_density = readRDS(paste0(out_path, "/RDS/peak_valley_locations_", adt_marker_select, "_", run_name, ".rds"))
       
        }else{
            print(adt_marker_select)
        }
    }

   
    peak_info = peak_valley_density$peak_landmark_list
    valley_info = peak_valley_density$valley_landmark_list %>% data.frame
    if(ncol(valley_info) > 1){
        if(sum(is.na(valley_info[, ncol(valley_info)])) == nrow(valley_info)){
            valley_info = valley_info[, -ncol(valley_info), drop = F]
        }
    }
    for(each_sample in unique(cell_x_feature$sample)){
        cell_ind <- which(cell_x_feature$sample == each_sample)
        if((ncol(valley_info) > 1) && !is.na(valley_info[each_sample, ncol(valley_info)])){
            mid_cell_ind <- which(cell_x_adt[cell_ind, adt_marker_select] >= valley_info[each_sample, 1] & cell_x_adt[cell_ind, adt_marker_select] <= valley_info[each_sample, ncol(valley_info)])
            pheno_matrix[cell_ind, adt_marker_select][mid_cell_ind] <- "mid"

            pos_cell_ind <- which(cell_x_adt[cell_ind, adt_marker_select] > valley_info[each_sample, ncol(valley_info)])
            pheno_matrix[cell_ind, adt_marker_select][pos_cell_ind] <- "+"
        }else{
            pos_cell_ind <- which(cell_x_adt[cell_ind, adt_marker_select] > valley_info[each_sample, 1])
            pheno_matrix[cell_ind, adt_marker_select][pos_cell_ind] <- "+"
        }
        
        
    }
}
pheno_matrix <- pheno_matrix %>% data.frame
pheno_matrix %>% head
saveRDS(pheno_matrix, file = paste0(out_path, "/RDS/pheno_matrix_", run_name, ".rds"))

pheno_matrix = readRDS(file = paste0(out_path, "/RDS/pheno_matrix_", run_name, ".rds"))

select_label = c("B", "CD4 T", "Treg", "CD8 T", "MAIT", "gamma delta T", "Monocytes CD14+", "Monocytes CD16+", "NK CD56-", "NK CD56+", "cDC", "pDC", "HSPC",  "Platelets",  "Plasmablast", "RBC")

cell_x_feature$cell_type_select = factor(cell_x_feature$CZI_initial_clustering, levels = c("CD4 T", "Treg", "CD8 T", "gamma delta T", "MAIT", "NK CD56-",  "NK CD56+", "B", "Plasmablast",  "Monocytes CD14+", "Monocytes CD16+", "cDC", "pDC", "HSPC",  "Platelets", "RBC", "Prolif lymph",  "Prolif Mono"), labels = c("CD4 T", "Treg", "CD8 T", "gamma delta T", "MAIT", "NK CD56-",  "NK CD56+", "B", "Plasmablast",  "Monocytes CD14+", "Monocytes CD16+", "cDC", "pDC", "HSPC",  "Platelets", "RBC", "Other Cells",  "Other Cells"))


auto_gating = list()
## for select cell types matching original cell types in the paper
auto_gating[["CD4 T"]] = pheno_matrix %>% dplyr::filter(CD2 == "+", CD3 == "+", CD5 == "+", CD19 == "-", CD20 == "-", CD4 == "+", CD8 == "-", CD14 == "-", CD16 =="-", CD371 == "-", CD56 == "-", HLADR == "-", CD38 == "-", CD11c == "-", CD25 == "-", CD235ab == "-", CD123 == "-") %>% rownames
auto_gating[["Treg"]] = pheno_matrix %>% dplyr::filter(CD3 == "+", CD5 == "+", CD19 == "-", CD20 == "-",  CD4 != "-", CD8 == "-", CD25 == "+", CD127 == "-", CD14 == "-", CD16 =="-", CD371 == "-", CD56 == "-", HLADR == "-", CD235ab == "-", CD123 == "-") %>% rownames

auto_gating[["CD8 T"]] = pheno_matrix %>% dplyr::filter(CD3 == "+", CD5 == "+", CD19 == "-", CD20 == "-", CD4 == "-", CD8 == "+", CD14 == "-", CD16 =="-", CD371 == "-", CD56 == "-", HLADR == "-", CD38 == "-", CD11c == "-", CD235ab == "-", CD123 == "-") %>% rownames
auto_gating[["MAIT"]] = pheno_matrix %>% dplyr::filter(CD3 == "+", CD5 == "+", CD19 == "-", CD20 == "-", CD4 == "-", CD8 == "-", CD161 == "+", CD26 == "+", KLRG1  == "+", CD14 == "-", CD16 =="-", CD371 == "-", CD56 == "-", HLADR == "-", CD38 == "-", CD11c == "-", CD235ab == "-", CD123 == "-") %>% rownames
auto_gating[["gamma delta T"]] = pheno_matrix %>% dplyr::filter(CD3 == "+", CD5 == "+", CD19 == "-", CD20 == "-", CD4 == "-", CD8 == "-", TCRVdelta2 == "+", CD14 == "-", CD16 =="-", CD371 == "-", CD56 == "-", HLADR == "-", CD38 == "-", CD11c == "-", CD235ab == "-", CD123 == "-") %>% rownames

auto_gating[["NK"]] = c(
    pheno_matrix %>% dplyr::filter(CD3 == "-", CD19 == "-", CD20 == "-", CD14 == "-",  HLADR == "-", CD56 == "+", CD335 == "+") %>% rownames, 
    pheno_matrix %>% dplyr::filter(CD3 == "-", CD19 == "-", CD20 == "-", CD14 == "-",  HLADR == "-", CD16 == "+") %>% rownames) %>% unique
auto_gating[["NK CD56+"]] = pheno_matrix %>% dplyr::filter(CD2 == "+", CD3 == "-", CD5 == "-", CD4 == "-", CD8 == "-", CD19 == "-", CD20 == "-", CD14 == "-",  HLADR == "-", CD371 == "-",  CD56 == "+", CD335 == "+") %>% rownames
auto_gating[["NK CD56-"]] = pheno_matrix %>% dplyr::filter(CD3 == "-", CD19 == "-", CD20 == "-", CD14 == "-",  HLADR == "-", CD16 == "+", CD56 == "-", CD371 == "-", CD38 == "+") %>% rownames

auto_gating[["B"]] = pheno_matrix %>% dplyr::filter(CD3 == "-", CD19 == "+", CD20 == "+", CD38 == "-") %>% rownames
auto_gating[["Plasmablast"]] = pheno_matrix %>% dplyr::filter(CD3 == "-", CD4 == "-", CD8 == "-", CD5 == "-", CD19 == "+", CD38 == "+", CD14 == "-", CD16 =="-", CD371 == "-", CD56 == "-",  CD11c == "-", CD235ab == "-", CD123 == "-") %>% rownames

auto_gating[["Monocytes"]] = c(
    pheno_matrix %>% dplyr::filter(CD3 == "-", CD19 == "-", CD20 == "-", CD14 == "+", CD16 == "-", CD371 == "+") %>% rownames, 
    pheno_matrix %>% dplyr::filter(CD3 == "-", CD19 == "-", CD20 == "-", CD14 == "-", CD16 == "+", CD371 == "+") %>% rownames) %>% unique
auto_gating[["Monocytes CD14+"]] = pheno_matrix %>% dplyr::filter(CD2 == "-", CD3 == "-", CD5 == "-", CD19 == "-", CD20 == "-", CD14 == "+", CD16 == "-", CD371 == "+", CD38 == "+", CD235ab == "-", CD123 == "-") %>% rownames
auto_gating[["Monocytes CD16+"]] = pheno_matrix %>% dplyr::filter(CD2 == "-", CD3 == "-", CD5 == "-", CD19 == "-", CD20 == "-", CD14 == "-", CD16 == "+", CD371 == "+", HLADR == "+", CD88 == "+", CD235ab == "-", CD123 == "-") %>% rownames ##HLADR CD88

auto_gating[["DC"]] = c(
    pheno_matrix %>% dplyr::filter(CD3 == "-", CD19 == "-", CD20 == "-", CD14 == "-", CD16 == "-", CD56 == "-", HLADR == "+", CD11c == "+") %>% rownames,
    pheno_matrix %>% dplyr::filter(CD3 == "-", CD19 == "-", CD20 == "-", CD14 == "-", CD16 == "-", CD56 == "-", HLADR == "+", CD123 == "+") %>% rownames) %>% unique
auto_gating[["cDC"]] = pheno_matrix %>% dplyr::filter(CD3 == "-", CD5 == "-", CD8 == "-", CD19 == "-", CD20 == "-", CD14 == "-", CD16 == "-", CD56 == "-", HLADR == "+", CD11c == "+", CD371 == "+", CD34 == "-", CD235ab == "-", CD38 == "-", CD88 == "-", CD123 == "-", CD161 == "-", CD304 == "-") %>% rownames 
auto_gating[["pDC"]] = pheno_matrix %>% dplyr::filter(CD3 == "-", CD8 == "-", CD19 == "-", CD20 == "-", CD14 == "-", CD16 == "-", CD56 == "-", HLADR == "+", CD123 == "+", CD304 == "+") %>% rownames


auto_gating[["HSPC"]] = pheno_matrix %>% dplyr::filter(CD2 == "-", CD3 == "-", CD5 == "-", CD4 == "-", CD8 == "-", CD19 == "-", CD20 == "-", CD14 == "-", CD16 == "-", CD56 == "-", CD34 == "+", HLADR == "+") %>% rownames

auto_gating[["Platelets"]] = pheno_matrix %>% dplyr::filter(CD3 == "-", CD5 == "-", CD4 == "-", CD8 == "-", CD19 == "-", CD20 == "-", CD14 == "-", CD16 == "-", CD56 == "-", HLADR == "-", CD34 == "-", CD41 == "+", CD49f == "+",  CD371 == "-",  CD235ab == "-", CD123 == "-", CD38 == "-", CD11c == "-") %>% rownames

auto_gating[["RBC"]] = pheno_matrix %>% dplyr::filter(CD3 == "-", CD19 == "-", CD20 == "-", CD14 == "-", CD16 == "-", CD56 == "-", HLADR == "-", CD34 == "-",  CD371 == "-", CD235ab == "+", CD123 == "-", CD38 == "-", CD11c == "-") %>% rownames

manual_gating <- list()
for(cT in select_label){
    manual_gating[[cT]] <- cell_x_feature %>% dplyr::filter(cell_type_select == cT) %>% rownames
}
cell_x_feature$auto_gating_overlap = "undefined"
cell_x_feature$auto_gating = "undefined"

acc_each_sample <- c()
for(cT in select_label){
    print(cT)
    print(auto_gating[[cT]] %>% length)
    overlap_ind <- which(auto_gating[[cT]] %in% manual_gating[[cT]])
    overlap_cell_name <- auto_gating[[cT]][overlap_ind]
    cell_x_feature$auto_gating_overlap[which(rownames(cell_x_feature) %in% overlap_cell_name)] = cT
    cell_x_feature$auto_gating[which(rownames(cell_x_feature) %in% auto_gating[[cT]])] = cT
    
    acc_each_sample_tmp <- left_join(
        cell_x_feature[auto_gating[[cT]], c("batch", "sample")] %>% group_by(batch, sample) %>% summarize(freq = n()),
        cell_x_feature[overlap_cell_name, c("batch", "sample")] %>% group_by(batch, sample) %>% summarize(freq = n()),
        by = c("batch", "sample")
    ) %>% mutate(accuracy = freq.y/freq.x * 100, cT = cT, cell_type_label = "Broad Labeling")
    
    if(cT %in% c("CD4 T", "CD8 T", "B", "NK CD56-", "Monocytes CD14+")){
        acc_each_sample_tmp = acc_each_sample_tmp %>% dplyr::filter(freq.y > 50)
    }else{
        acc_each_sample_tmp = acc_each_sample_tmp %>% dplyr::filter(freq.y > 5)
    }
    acc_each_sample = rbind(acc_each_sample, acc_each_sample_tmp)

}

acc_each_sample_df = cell_x_feature %>% 
dplyr::filter(cell_type_select %in% select_label) %>% 
group_by(sample, cell_type_select) %>% 
summarize(cell_num_manual_label = n()) %>% 
mutate(cT = cell_type_select) %>% 
dplyr::select(sample, cT, cell_num_manual_label) %>% 
data.frame %>% 
left_join(acc_each_sample, ., by = c("sample", "cT")) %>% data.frame

acc_each_sample_df$cT = factor(acc_each_sample_df$cT, levels = select_label)
cell_x_feature$auto_gating_select = factor(cell_x_feature$auto_gating_select, levels = c(select_label, "undefined"))


file_list = adt_feature$batch %>% levels

pdf(paste0(out_path, "/Figures/auto_gating_consistency.pdf"), width = 11, height = 11)
acc_each_sample_df %>% dplyr::filter(cell_num_manual_label > 0) %>%  
ggplot(aes(x = cT, y = accuracy)) +
geom_boxplot(outlier.shape = NA, position = position_dodge2(preserve = "single")) +
geom_jitter(size = 2, position = position_dodge(width = 0.4), aes(color = batch)) +
xlab("") +
ylab("Auto-gating Consistency with RNA-based Annotation") +
theme_bw(base_size = 20) +
scale_color_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))) +
rotate_x_text(angle = 45) +
theme(legend.position = "top") +
rremove("legend.title") 
dev.off()

## umap with auto-gating label
method = "ADTnorm_sample_manualtuning"
# filename <- paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds")
# res_norm = readRDS(file = filename)

citeNorm = readRDS(file = paste0(out_path, "/RDS/citeNorm_obj_halfdata_", run_name, "_", method, "_rna30_adt30_knn50_mostseparable.rds"))

ans_umap = Embeddings(citeNorm, reduction = "wnn.umap") %>% data.frame
colnames(ans_umap) = c("X1", "X2")
set.seed(20240215)
cell_ind = sample(1:nrow(adt_data), round(nrow(adt_data)/2), replace = FALSE) %>% sort

cell_x_feature$auto_gating_label = factor(cell_x_feature$auto_gating_overlap, levels = c(select_label, "undefined"))

reindex1 = which(cell_x_feature[cell_ind, ]$auto_gating_label == "undefined")
reindex2 = which(cell_x_feature[cell_ind, ]$auto_gating_label != "undefined")
reindex = c(reindex1, reindex2)

pdf(paste0(out_path, "/Figures/auto_gating_wnn_umap.pdf"), width = 11, height = 11)
print(plot_umap_extend_raster(ans_umap[reindex, ], cell_x_feature[cell_ind, ][reindex, ], point_size = 1, parameter_list = list(
target_feature = "auto_gating_label", 
color_break = c("CD4 T", "Treg", "CD8 T", "gamma delta T", "MAIT", "NK CD56-",  "NK CD56+", "B", "Plasmablast",  "Monocytes CD14+", "Monocytes CD16+", "cDC", "pDC", "HSPC",  "Platelets", "RBC", "undefined"),
    color_design = c("#AADBDF", "#E599E8", "#ABDFAE", "#CC5C5D", "#EAA6AC", "#E4E7A8",  "#385F9E", "#D97E4A", "#4CB45E",  "#ABA9D9", "#F8BD98", "#67A6E4", "#D9DFA9", "#C487EB",  "#F6B6CE", "#E43324", "grey"),
method_label = "Auto-gating"
)))

dev.off()


method = "ADTnorm_sample_manualtuning"
res_norm = readRDS(file = paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds"))
npc = 15
nneighbors = 20
mindist = 0.001
## SVD version
svd_matrix = svd(res_norm[cell_ind, ], nu = npc)$u
umap_svd = uwot::umap(svd_matrix, n_neighbors = nneighbors, min_dist = mindist, n_components = 2)
ans_umap = umap_svd %>% data.frame
colnames(ans_umap) = c("X1", "X2")

cell_x_feature$auto_gating_label = factor(cell_x_feature$auto_gating_overlap, levels = c(select_label, "undefined"))

reindex1 = which(cell_x_feature[cell_ind, ]$auto_gating_label == "undefined")
reindex2 = which(cell_x_feature[cell_ind, ]$auto_gating_label != "undefined")
reindex = c(reindex1, reindex2)

pdf(paste0(out_path, "/Figures/auto_gating_adt_umap.pdf"), width = 11, height = 11)
print(plot_umap_extend_raster(ans_umap[reindex, ], cell_x_feature[cell_ind, ][reindex, ], point_size = 1, parameter_list = list(
target_feature = "auto_gating_label", 
color_break = c("CD4 T", "Treg", "CD8 T", "gamma delta T", "MAIT", "NK CD56-",  "NK CD56+", "B", "Plasmablast",  "Monocytes CD14+", "Monocytes CD16+", "cDC", "pDC", "HSPC",  "Platelets", "RBC", "undefined"),
    color_design = c("#AADBDF", "#E599E8", "#ABDFAE", "#CC5C5D", "#EAA6AC", "#E4E7A8",  "#385F9E", "#D97E4A", "#4CB45E",  "#ABA9D9", "#F8BD98", "#67A6E4", "#D9DFA9", "#C487EB",  "#F6B6CE", "#E43324", "grey"),
method_label = "Auto-gating"
)))

dev.off()

