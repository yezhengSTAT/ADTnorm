library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(httpgd)
library(ggpubr)
library(ADTnorm)

run_name = "public13Dataset_CITEseq"

master_path = "./"
in_path = "./publicData_CITEseq/data/"
out_path = paste0(master_path, "manuscript/results/", run_name)
fcs_path = paste0(out_path, "/FCS/")
fig_path = paste0(out_path, "/Figures/")

file_list <- c(
    "10X_pbmc_10k", "10X_pbmc_1k", "10X_pbmc_5k_v3", "10X_pbmc_5k_nextgem", "10X_malt_10k", 
    "stuart_2019", "granja_2019_pbmc", "granja_2019_bmmc", "hao_2020",
    "kotliarov_2020", "witkowski_2020", "triana_2021", "buus_2021_T"   
)


marker_list = c("CD3", "CD4", "CD8",  "CD14", "CD19", "CD25", "CD45RA", "CD56", "CD127")
adt_data_full = readRDS(file = paste0(out_path, "/RDS/adt_data_full_RawCount_", run_name, ".rds"))
adt_data = readRDS(file = paste0(out_path, "/RDS/adt_data_RawCount_", run_name, ".rds"))
adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", run_name, ".rds"))

cell_x_adt = arcsinh_transform(cell_x_adt = adt_data) 
cell_x_feature = adt_feature
rownames(cell_x_adt) = paste0("Cell", 1:nrow(cell_x_adt))
rownames(cell_x_feature) = paste0("Cell", 1:nrow(cell_x_feature))

pheno_matrix <- matrix("-", nrow(cell_x_adt), ncol(cell_x_adt))
colnames(pheno_matrix) <- colnames(cell_x_adt)
rownames(pheno_matrix) <- rownames(cell_x_adt)

for(adt_marker_select in colnames(cell_x_adt)){
    peak_valley_density = readRDS(paste0(out_path, "/RDS/peak_valley_locations_", adt_marker_select, "_ADTnorm_sample_manual_keepZero.rds"))
    
    peak_info = peak_valley_density$peak_landmark_list
    valley_info = peak_valley_density$valley_landmark_list
    
    for(each_sample in unique(cell_x_feature$sample)){
        cell_ind <- which(cell_x_feature$sample == each_sample)
        if(ncol(valley_info) > 1){
            mid_cell_ind <- which(cell_x_adt[cell_ind, adt_marker_select] >= valley_info[each_sample, 1] & cell_x_adt[cell_ind, adt_marker_select] <= valley_info[each_sample, ncol(valley_info)])
            pheno_matrix[cell_ind, adt_marker_select][mid_cell_ind] <- "mid"
        }
        pos_cell_ind <- which(cell_x_adt[cell_ind, adt_marker_select] > valley_info[each_sample, ncol(valley_info)])
        pheno_matrix[cell_ind, adt_marker_select][pos_cell_ind] <- "+"
        
    }
}
pheno_matrix <- pheno_matrix %>% data.frame
pheno_matrix %>% head

broad_label = c("B", "CD4 T", "CD8 T", "monocytes", "NK", "DCs")
refined_label = c( 'naive CD4','memory CD4', 'Treg', 'naive CD8', 'memory CD8')

auto_gating = list()
## broad cell types
auto_gating[["B"]] = pheno_matrix %>% filter(CD3 == "-", CD19 == "+") %>% rownames
auto_gating[["CD4 T"]] = pheno_matrix %>% filter(CD3 == "+", CD19 == "-", CD4 == "+", CD8 == "-") %>% rownames
auto_gating[["CD8 T"]] = pheno_matrix %>% filter(CD3 == "+", CD19 == "-", CD4 == "-", CD8 == "+") %>% rownames
auto_gating[["monocytes"]] = pheno_matrix %>% filter(CD3 == "-", CD19 == "-", CD14 == "+") %>% rownames
auto_gating[["NK"]] = pheno_matrix %>% filter(CD3 == "-", CD19 == "-", CD14 == "-", CD56 == "+") %>% rownames
auto_gating[["DCs"]] = pheno_matrix %>% filter(CD3 == "-", CD19 == "-", CD14 == "-", CD56 == "-") %>% rownames

## refined cell types
auto_gating[["naive CD4"]] = pheno_matrix %>% filter(CD3 == "+", CD19 == "-", CD4 == "+", CD8 == "-", CD45RA == "+") %>% rownames
auto_gating[["memory CD4"]] = pheno_matrix %>% filter(CD3 == "+", CD19 == "-", CD4 == "+", CD8 == "-", CD45RA == "-") %>% rownames
auto_gating[["Treg"]] = pheno_matrix %>% filter(CD3 == "+", CD19 == "-", CD4 == "+", CD8 == "-", CD25 == "+", CD127 == "-") %>% rownames
auto_gating[["naive CD8"]] = pheno_matrix %>% filter(CD3 == "+", CD19 == "-", CD4 == "-", CD8 == "+", CD45RA == "+") %>% rownames
auto_gating[["memory CD8"]] = pheno_matrix %>% filter(CD3 == "+", CD19 == "-", CD4 == "-", CD8 == "+", CD45RA == "-") %>% rownames

## manual gating label as gold standard
manual_gating <- list()
for(cT in broad_label){
    manual_gating[[cT]] <- cell_x_feature %>% dplyr::filter(cell_type_l1 == cT) %>% rownames

}
for(cT in refined_label){
    manual_gating[[cT]] <- cell_x_feature %>% dplyr::filter(cell_type_l2 == cT) %>% rownames
}

## accuracy check
cell_x_feature$auto_gating_l1 = "undefined"
cell_x_feature$auto_gating_l2 = "undefined"
cell_x_feature$auto_gating_mixed = "undefined"

acc_each_sample <- c()
for(cT in broad_label){
    overlap_ind <- which(auto_gating[[cT]] %in% manual_gating[[cT]])
    overlap_cell_name <- auto_gating[[cT]][overlap_ind]
    cell_x_feature$auto_gating_l1[which(rownames(cell_x_feature) %in% overlap_cell_name)] = cT
    cell_x_feature$auto_gating_mixed[which(rownames(cell_x_feature) %in% overlap_cell_name)] = cT
    
    acc_each_sample <- left_join(
        cell_x_feature[auto_gating[[cT]], c("study_name", "sample")] %>% group_by(study_name, sample) %>% summarize(freq = n()),
        cell_x_feature[overlap_cell_name, c("study_name", "sample")] %>% group_by(study_name, sample) %>% summarize(freq = n()),
        by = c("study_name", "sample")
    ) %>% filter(freq.x >= 50) %>% mutate(accuracy = freq.y/freq.x * 100, cT = cT, cell_type_label = "Broad Labeling") %>% rbind(acc_each_sample, .)

}
for(cT in refined_label){
    overlap_ind <- which(auto_gating[[cT]] %in% manual_gating[[cT]])
    overlap_cell_name <- auto_gating[[cT]][overlap_ind]
    cell_x_feature$auto_gating_l2[which(rownames(cell_x_feature) %in% overlap_cell_name)] = cT
    cell_x_feature$auto_gating_mixed[which(rownames(cell_x_feature) %in% overlap_cell_name)] = cT

    acc_each_sample <- left_join(
        cell_x_feature[auto_gating[[cT]], c("study_name", "sample")] %>% group_by(study_name, sample) %>% summarize(freq = n()),
        cell_x_feature[overlap_cell_name, c("study_name", "sample")] %>% group_by(study_name, sample) %>% summarize(freq = n()),
        by = c("study_name", "sample")
    ) %>% filter(freq.x >= 50) %>% mutate(accuracy = freq.y/freq.x * 100, cT = cT, cell_type_label = "Refined Labeling") %>% rbind(acc_each_sample, .)

}

acc_each_sample_df = rbind(
    cell_x_feature %>% group_by(sample, study_name, cell_type_l1) %>% summarize(cell_num_manual_label = n()) %>% mutate(cT = cell_type_l1) %>% select(study_name, sample, cT, cell_num_manual_label) %>% data.frame,
    cell_x_feature %>% group_by(sample, study_name, cell_type_l2) %>% summarize(cell_num_manual_label = n()) %>% mutate(cT = cell_type_l2) %>% select(study_name, sample, cT, cell_num_manual_label) %>% data.frame
) %>% left_join(acc_each_sample, ., by = c("study_name", "sample", "cT")) %>% data.frame

acc_each_sample_df$cT = factor(acc_each_sample_df$cT, levels = factor(c(broad_label, refined_label)))
cell_x_feature$auto_gating_l1 = factor(cell_x_feature$auto_gating_l1, levels = c("B", "CD4 T", "CD8 T", "monocytes", "NK", "DCs", "undefined"))
cell_x_feature$auto_gating_l2 = factor(cell_x_feature$auto_gating_l2, levels = c(
    "naive CD4", "memory CD4", "Treg", 
    "naive CD8", "memory CD8", "undefined"))
cell_x_feature$auto_gating_mixed[which(cell_x_feature$auto_gating_mixed %in% c("CD4 T", "CD8 T"))] = "undefined"

cell_x_feature$auto_gating_mixed = factor(
    cell_x_feature$auto_gating_mixed, 
    levels = c("B", "naive CD4", "memory CD4", "Treg", "naive CD8", "memory CD8", "monocytes", "NK", "DCs", "undefined"))

acc_each_sample_df %>% filter(cell_num_manual_label > 100) %>%  
ggplot(aes(x = cT, y = accuracy, color = study_name, fill = study_name)) +
geom_boxplot(lwd=1.5, outlier.shape = NA, position = position_dodge2(preserve = "single")) +
facet_grid(~cell_type_label, scale = "free", space = "free") +
xlab("") +
ylab("Auto-gating Accuracy (100%)") +
theme_bw(base_size = 20) +
scale_color_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))) +
scale_fill_manual(breaks = file_list, values = alpha(colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list)), 0.7)) +
rotate_x_text(angle = 0) +
theme(legend.position = "top") +
rremove("legend.title") +
ggsave(filename = paste0(out_path, "/Figures/auto_gating_phenotype_accuracy_breakdowncelltype_breakdownstudy.pdf"), width = 18, height = 10)


acc_each_sample_df %>% filter(cell_num_manual_label > 100) %>% group_by(cT, cell_type_label, study_name) %>% summarize(median_accuracy = median(accuracy), mean_accuracy = mean(accuracy)) %>% 
ggplot(aes(x = cT, y = median_accuracy)) +
geom_boxplot(outlier.shape = NA) +
geom_jitter(aes(color = study_name), width = 0.2, size = 8, shape = 19) +
facet_grid(~cell_type_label, scale = "free", space = "free") +
xlab("") +
ylab("Auto-gating Accuracy (100%)") +
theme_bw(base_size = 23) +
scale_color_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))) +
# scale_fill_manual(breaks = file_list, values = alpha(colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list)), 0.7)) +
rotate_x_text(angle = 0) +
theme(legend.position = "top") +
rremove("legend.title") 
ggsave(filename = paste0(out_path, "/Figures/auto_gating_phenotype_accuracy_breakdowncelltype_breakdownstudy_aggregateboxplot.pdf"), width = 20, height = 12)


color_breaks = c("B", "CD4 T", "CD8 T", "monocytes", "NK", "DCs", "naive CD4", "memory CD4", "Treg", "naive CD8", "memory CD8")
color_values = c("#B00000", "#3F918B", "#896191", "#FF980A", "#226fa7", "#F781BF", "#9ef954", "#3F918B", "#03500e", "#896191", "#350154")

acc_each_sample_df %>% filter(cell_num_manual_label > 100) %>%  
ggplot(aes(x = cT, y = accuracy, color = cT, fill = cT)) +
# geom_bar(stat = "identity", position = "dodge") +
geom_boxplot(outlier.shape = NA, position = position_dodge2(preserve = "single")) +
geom_jitter() +
facet_grid(~cell_type_label, scale = "free", space = "free") +
xlab("") +
ylab("Auto-gating Accuracy (100%)") +
theme_bw(base_size = 20) +
scale_color_manual(breaks = color_breaks, values = color_values) +
scale_fill_manual(breaks = color_breaks, values = alpha(color_values, 0.7)) +
rotate_x_text(angle = 0) +
rremove("legend") +
ggsave(filename = paste0(out_path, "/Figures/auto_gating_phenotype_accuracy_breakdowncelltype.pdf"), width = 16, height = 9)

acc_each_sample_df %>% filter(cell_num_manual_label > 100) %>%  
ggplot(aes(x = study_name, y = accuracy, color = study_name, fill = study_name)) +
# geom_bar(stat = "identity", position = "dodge") +
geom_boxplot(outlier.shape = NA, position = position_dodge2(preserve = "single")) +
geom_jitter() +
xlab("") +
ylab("Auto-gating Accuracy (100%)") +
theme_bw(base_size = 20) +
scale_color_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))) +
scale_fill_manual(breaks = file_list, values = alpha(colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list)), 0.7)) +
rotate_x_text(angle = 20) +
rremove("legend") +
ggsave(filename = paste0(out_path, "/Figures/auto_gating_phenotype_accuracy_breakdownstudy_name.pdf"), width = 16, height = 9)


## umap with auto-gating label
method = "ADTnorm_sample_manual_keepZero"
tmp = readRDS(file = paste0(out_path, "/RDS/adt_", method, "_", run_name, "_pca_umap_nneighbors20_mindist0.001_cosine_uwot.rds"))
ans_umap = tmp[[2]]

pdf(paste0(out_path, "/Figures/auto_gating_umap.pdf"), width = 15, height = 15)
reindex1 = which(cell_x_feature$auto_gating_mixed == "undefined")
reindex2 = which(cell_x_feature$auto_gating_mixed != "undefined")
reindex = c(reindex1, reindex2)

print(plot_umap_raster(ans_umap[reindex, ], cell_x_feature[reindex, ], point_size = 0.3, parameter_list = list(
target_feature = "auto_gating_mixed",
color_design = c("#B00000",  "#9ef954", "#3F918B", "#03500e",
    "#896191", "#350154", "#FF980A", "#226fa7", "#F781BF",
    "grey"),
#color_design = c(colorRampPalette(brewer.pal(8, "Dark2"))(14), "grey"),
color_break = c("B", "naive CD4", "memory CD4", "Treg", 
    "naive CD8", "memory CD8", "monocytes", "NK", "DCs", 
     "undefined"),
method_label = "Auto-gating"
)))

cell_x_feature$auto_gating_mixed = factor(
    cell_x_feature$auto_gating_mixed, 
    levels = c("B", "naive CD4", "memory CD4", "Treg", "naive CD8", "memory CD8", "monocytes", "NK", "DCs", "undefined"), 
    labels = c("B (CD19+)", "naive CD4 (CD3+CD4+CD45RA+)", "memory CD4 (CD3+CD4+CD45RA-)", "Treg (CD3+CD4+CD25+CD127+)", "naive CD8 (CD3+CD8+CD45RA+)", "memory CD8 (CD3+CD8+CD45RA-)", "monocytes (CD3-CD19-CD14+)", "NK (CD3-CD19-CD14-CD56+)", "DCs (CD3-CD19-CD14-CD56-)", "undefined"))
print(plot_umap_raster(ans_umap[reindex, ], cell_x_feature[reindex, ], point_size = 0.3, parameter_list = list(
target_feature = "auto_gating_mixed",
color_design = c("#B00000",  "#9ef954", "#3F918B", "#03500e",
    "#896191", "#350154", "#FF980A", "#226fa7", "#F781BF",
    "grey"),
#color_design = c(colorRampPalette(brewer.pal(8, "Dark2"))(14), "grey"),
color_break = c("B (CD19+)", "naive CD4 (CD3+CD4+CD45RA+)", "memory CD4 (CD3+CD4+CD45RA-)", "Treg (CD3+CD4+CD25+CD127+)", "naive CD8 (CD3+CD8+CD45RA+)", "memory CD8 (CD3+CD8+CD45RA-)", "monocytes (CD3-CD19-CD14+)", "NK (CD3-CD19-CD14-CD56+)", "DCs (CD3-CD19-CD14-CD56-)", "undefined"),
method_label = "Auto-gating"
)))

reindex1 = which(cell_x_feature$auto_gating_l1 == "undefined")
reindex2 = which(cell_x_feature$auto_gating_l1 != "undefined")
reindex = c(reindex1, reindex2)

print(plot_umap_raster(ans_umap[reindex, ], cell_x_feature[reindex, ], point_size = 0.3, parameter_list = list(
target_feature = "auto_gating_l1", 
color_design = c("#B00000", "#3F918B", "#896191", "#FF980A", "#226fa7", "#F781BF", "grey"), 
#color_design = c(colorRampPalette(brewer.pal(8, "Set1"))(6), "grey"),
color_break = c("B", "CD4 T", "CD8 T", "monocytes", "NK", "DCs", "undefined"),
method_label = "Auto-gating"
)))


reindex1 = which(cell_x_feature$auto_gating_l2 == "undefined")
reindex2 = which(cell_x_feature$auto_gating_l2 != "undefined")
reindex = c(reindex1, reindex2)

print(plot_umap_raster(ans_umap[reindex, ], cell_x_feature[reindex, ], point_size = 0.3, parameter_list = list(
target_feature = "auto_gating_l2",
color_design = c(
    "#9ef954", "#3F918B", "#03500e",
    "#896191", "#350154", "grey"),
#color_design = c(colorRampPalette(brewer.pal(8, "Dark2"))(14), "grey"),
color_break = c(
    "naive CD4", "memory CD4", "Treg", 
    "naive CD8", "memory CD8", "undefined"),
method_label = "Auto-gating"
)))

dev.off()