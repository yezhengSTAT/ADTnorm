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
library(gghalves)
library(ADTnorm)
# p_load(flowStats, flowCore, FlowSOM, ncdfFlow, flowViz, pdfCluster, cluster)

## =====================
## study name and paths
## =====================
run_name = "hema_progenitor_titration"
master_path = "./"
out_path = paste0(master_path, "manuscript/results/", run_name)
# fcs_path = paste0(out_path, "/FCS/")
fig_path = paste0(out_path, "/Figures_", run_name, "/")

## =================
## load titration data
## =================
# in_path = "./data/Zhang2024/"
# annotation = fread(paste0(in_path, "/GEO/GSE245108_Level3R-Titrated-Titration-Annotation.txt"), header = TRUE) %>% mutate(filename = gsub(".*\\.", "", UID), barcode = gsub("\\..*", "", UID))
# cd271_obj = Read10X_h5(paste0(in_path, "/GEO/GSE245108_CD271-1_filtered_feature_bc_matrix.h5"))
# protein_list = rownames(cd271_obj$`Antibody Capture`)[1:275]

# adt_data = c()
# adt_feature = c()
# for(item_name in c("CD271-1", "CD271-2", "CD271-3", "CD271-4", "CD34-1", "CD34-2", "CD34-3", "CD34-4", "TNC-1-2-1", "TNC-1-2-2", "TNC-3-4-1", "TNC-3-4-2")){
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
# adt_feature$donor = factor(adt_feature$filename, levels = c("CD271-1", "CD271-2", "CD271-3", "CD271-4", "CD34-1", "CD34-2", "CD34-3", "CD34-4", "TNC-1-2-1", "TNC-1-2-2", "TNC-3-4-1", "TNC-3-4-2"), labels = c(paste0("Donor", c(1:4, 1:4,1:4))))
# adt_feature$subset = factor(adt_feature$filename, levels = c("CD271-1", "CD271-2", "CD271-3", "CD271-4", "CD34-1", "CD34-2", "CD34-3", "CD34-4", "TNC-1-2-1", "TNC-1-2-2", "TNC-3-4-1", "TNC-3-4-2"), labels = c(rep(c("CD271", "CD34", "BMNC"), each = 4)))
# adt_feature$sample = paste0(adt_feature$donor, "_", adt_feature$subset, "_", adt_feature$Concentration)
# adt_feature$batch = adt_feature$sample

# cell_x_adt = data.frame(t(adt_data))
# cell_x_feature = data.frame(adt_feature)
# colnames(cell_x_adt) = colnames(cell_x_adt) %>% gsub("HuMs\\.", "", .) %>% gsub("Hu\\.", "", .) %>% gsub("Isotype_", "Isotype", .) %>% gsub("_.*", "", .)

# rm_cell = which(adt_feature$Concentration == "Titrated" | adt_feature$Concentration == "1.5")
# cell_x_adt = cell_x_adt[-rm_cell, ]
# cell_x_feature = cell_x_feature[-rm_cell, ]
# saveRDS(cell_x_adt, paste0(out_path, "/RDS/cell_x_adt_", run_name, ".rds"))
# saveRDS(cell_x_feature, paste0(out_path, "/RDS/cell_x_feature_", run_name, ".rds"))

adt_data = readRDS(paste0(out_path, "/RDS/cell_x_adt_", run_name, ".rds"))
adt_feature = readRDS(paste0(out_path, "/RDS/cell_x_feature_", run_name, ".rds"))

## ================
## QC and ADTnorm
## ================

## check the unique value number for each protein marker for each sample in adt_feature
uniq_num = cbind(adt_data, adt_feature) %>% dplyr::group_by(sample) %>% summarise_all(n_distinct) %>% as.data.frame
uniq_num[, 2:276] %>% colMeans %>% sort
uniq_num[, 2:276] %>% apply(., 2, median) %>% sort
uniq_num[, 2:276] %>% apply(., 2, max) %>% sort

## only process adt markers that at least have 10 counts per sample per concentration
adt_marker_filter = uniq_num[, 2:276] %>% apply(., 2, max) %>% data.frame %>% dplyr::filter(. >=10) %>% rownames

# adt_feature$batch = factor(adt_feature$Concentration, levels = c("0.25", "0.5", "1", "2", "4"))
adt_feature$sample = factor(
    adt_feature$sample, 
    levels = paste(
        rep(rep(c("Donor1", "Donor2", "Donor3", "Donor4"), each = 5), 3),
        rep(c("CD271", "CD34", "BMNC"), each = 20),
        rep(c("0.25", "0.5", "1", "2", "4"), 12), sep = "_")
)
adt_feature$batch = adt_feature$sample

res_norm = ADTnorm(customize_landmark = TRUE, # FALSE
save_fig = FALSE,
save_landmark = FALSE,
cell_x_adt = adt_data,
cell_x_feature = adt_feature,
save_outpath = out_path,
study_name = run_name,
trimodal_marker = NULL, #c("CD4"),
bw_smallest_bi = 1.1,
bw_smallest_tri = 1.1,
bw_smallest_adjustments = list(CD3 = 0.8, CD4 = 1.1, CD8 = 0.8),
shoulder_valley_slope = -1,
marker_to_process = "IgM",
exclude_zeroes = FALSE,
bimodal_marker = colnames(adt_data)[colnames(adt_data) != "CD328"],
positive_peak = NULL,
quantile_clip = 1,
peak_type = "mode",
multi_sample_per_batch = FALSE,
shoulder_valley = TRUE,
valley_density_adjust = 3,
landmark_align_type = "negPeak_valley_posPeak",
midpoint_type = "valley",
neg_candidate_thres = asinh(2/5 + 1),
lower_peak_thres = 0.01,
brewer_palettes = "Set1",
detect_outlier_valley = FALSE,
target_landmark_location = NULL,
clean_adt_name = FALSE,
override_landmark = NULL,
verbose = TRUE)

saveRDS(res_norm, file = paste0(out_path, "/RDS/adt_data_norm_ADTnorm_sample_manual_", run_name, ".rds"))

## =============
## UMAP
## =============
cell_select = which(adt_feature$Concentration == "1")
adt_data_select = adt_data[cell_select, ]
adt_feature_select = adt_feature[cell_select, ]

adt_obj = CreateSeuratObject(counts = adt_data_select, data = adt_data_select)
adt_obj = SetAssayData(adt_obj, slot = "scale.data", new.data = t(res_norm))
adt_obj = FindVariableFeatures(adt_obj, nfeature = 2000)
adt_obj = AddMetaData(adt_obj, metadata = adt_feature)
adt_obj = RunPCA(adt_obj)
adt_obj = FindNeighbors(adt_obj, dims = 1:50, reduction = "pca")
adt_obj = FindClusters(adt_obj, resolution = 2.5, algorithm = 2)
adt_obj = RunUMAP(adt_obj, dims = 1:50, reduction = "pca", min.dist = 0.01, n.neighbors = 20)

DimPlot(adt_obj, reduction = "umap", label = TRUE, group.by = "seurat_clusters", cols = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(26))
DimPlot(adt_obj, reduction = "umap", label = TRUE, group.by = "type")
DimPlot(adt_obj, reduction = "umap", label = TRUE, group.by = "batch")


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
peak_num_summary = c()
peak_sep_summary = c()

for(adt_marker_select in colnames(cell_x_adt)){ ##
    
    peak_valley_density = readRDS(paste0(out_path, "/RDS/peak_valley_locations_", adt_marker_select, "_", run_name, ".rds"))
    
    peak_info = peak_valley_density$peak_landmark_list
    valley_info = peak_valley_density$valley_landmark_list

    for(each_sample in unique(adt_feature$sample)){
        batch_info = each_sample
        each_peak_info = peak_info[each_sample, ]
        peak_num = sum(is.na(each_peak_info) == FALSE)
        peak_num_summary = data.frame(
            peak_num = peak_num, 
            batch = batch_info, 
            sample = each_sample, adt_marker = adt_marker_select) %>% 
            rbind(peak_num_summary, .)
        if(peak_num == 1){
            peak_sep_summary = data.frame(
                sep_power = one_peak_stain_quality(adt_feature, cell_x_adt, adt_marker_select, each_sample, peak_info, valley_info, peak_num = peak_num),
                adt_marker = adt_marker_select, sample = each_sample, batch = batch_info, peak_num = paste0("# of peak: ", peak_num)) %>% rbind(peak_sep_summary, .)
        }
        if(peak_num == 2){
            peak_sep_summary = data.frame(
                sep_power = two_peak_stain_quality(adt_feature, cell_x_adt, adt_marker_select, each_sample, peak_info, valley_info, peak_num = peak_num),
                adt_marker = adt_marker_select, sample = each_sample, batch = batch_info, peak_num = paste0("# of peak: ", peak_num)) %>% rbind(peak_sep_summary, .)
        }
        if(peak_num > 2){
            peak_sep_summary = data.frame(
                sep_power = multi_peak_stain_quality(adt_feature, cell_x_adt, adt_marker_select, each_sample, peak_info, valley_info, peak_num = peak_num),
                adt_marker = adt_marker_select, sample = each_sample, batch = batch_info, peak_num = paste0("# of peak: ", peak_num)) %>% rbind(peak_sep_summary, .)
        }
    }
}

peak_sep_summary$Donor = factor(gsub("_.*", "", peak_sep_summary$sample), levels = paste0("Donor", c(1:4)))
peak_sep_summary$Capture = factor(gsub(".*_(.*)_.*", "\\1", peak_sep_summary$sample), levels = c("CD271", "CD34", "BMNC"))
peak_sep_summary$Concentration = factor(gsub(".*_", "", peak_sep_summary$sample), levels = c("0.25", "0.5", "1", "2", "4"))

saveRDS(peak_sep_summary, file = paste0(out_path, "/RDS/peak_sep_summary_", run_name, ".rds"))

peak_sep_summary = readRDS(paste0(out_path, "/RDS/peak_sep_summary_", run_name, ".rds"))
peak_sep_summary$Capture = factor(peak_sep_summary$Capture, levels = c("CD34","CD271","BMNC"), labels = c("CD34high", "CD34+CD271+", "BMNC"))
head(peak_sep_summary)


library(gghalves)

p_allTitration = peak_sep_summary %>% 
ggplot(aes(x = Capture, y = sep_power, fill = Concentration)) +
geom_violin(position = position_dodge(width = 0.9)) +
geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9, seed = 1), size = 0.5) +
theme_bw(base_size = 25) +
scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(5)) +
scale_color_brewer(palette = "Set1") +
ylab("Stain Quality Score") +
xlab("") +
theme(legend.position = "top") +
# rotate_x_text(angle = 45) +
scale_y_continuous(trans = "log10") +
geom_hline(yintercept = 2, linewidth =2, linetype = "dashed", color = "grey30")

p_allTitration
ggsave(paste0(fig_path, "titration_stain_quality_score_acrossCapture_violin.pdf"), width = 8, height = 9)


p_filterTitration = peak_sep_summary %>% dplyr::filter(sep_power > 2) %>%
  ggplot(aes(x = Capture, y = sep_power, fill = Concentration)) +
  geom_boxplot() +
  theme_bw(base_size = 25) +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(5)) +
  ylab("Stain Quality Score") +
  xlab("") +
  theme(legend.position = "top") +
  # rotate_x_text(angle = 45) +
  scale_y_continuous(trans = "log10")

p_filterTitration
ggsave(paste0(fig_path, "titration_stain_quality_score_acrossCapture_filter_boxplot.pdf"), width = 8, height = 9)

adt_marker_select = peak_sep_summary %>% dplyr::filter(sep_power > 2) %$% adt_marker %>% table %>% sort %>% data.frame %>% dplyr::filter(Freq > 1) %$% .

## ================
## density plot of 
#  [1] IgM     CD235ab CD94    CD49a   CD46    CD86    CD39    CD49d   HLA.DR 
# [10] CD244   CD99    CLEC12A
## ================
require(magrittr)
require(ggplot2)
require(RColorBrewer)
require(tidyr)
require(ggridges)
require(cytotidyr)

plot_adt_density <- function(cell_x_adt, cell_x_feature, adt_marker_list, unit, parameter_list = NULL){
    
    if ("method_label" %in% names(parameter_list)) {
        method_label = parameter_list[["method_label"]]
    } else {
        method_label <- ""
    }
    if ("bw" %in% names(parameter_list)) {
        bw = parameter_list$bw
    } else {
        bw = 0.1
    }
    if(unit == "study"){
      tmpProfile <- cell_x_adt %>% 
        dplyr::select(all_of(adt_marker_list)) %>%
        data.frame() %>%
        gather(key = "ADT", value = "counts") %>%
        mutate(
            sample = rep(cell_x_feature$study_name, length(adt_marker_list)),
            batch = rep(cell_x_feature$study_name, length(adt_marker_list)),
            disease = rep(cell_x_feature$sample_status, length(adt_marker_list))
        ) 
    }else if(unit == "sample"){
      tmpProfile <- cell_x_adt %>% 
        dplyr::select(all_of(adt_marker_list)) %>%
        data.frame() %>%
        gather(key = "ADT", value = "counts") %>%
        mutate(
            sample = rep(cell_x_feature$sample, length(adt_marker_list)),
            batch = rep(cell_x_feature$batch, length(adt_marker_list)),
            disease = rep(cell_x_feature$subset, length(adt_marker_list))
        ) 
    }
    fillColor <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(unique(tmpProfile$batch)))
    resPlot <- ggplot(tmpProfile, aes(x = counts, y = sample, fill = batch)) +
        geom_density_ridges(bandwidth = bw, scale = 2.5) +
        facet_grid( ~ factor(ADT, levels = adt_marker_list), scales = "free") + #, space = "free"
        theme_bw(base_size = 40) +
        xlab(method_label) +
        ylab("") +
        ggpubr::rotate_x_text(angle = 90) +
        ggpubr::rremove("legend") +
        scale_fill_manual(values = fillColor) +
        ggpubr::rremove("legend.title")
        #  +
        # coord_cartesian(xlim = c(0, 7.5))

    if(("arcsinhTransform" %in% names(parameter_list)) && (parameter_list[["arcsinhTransform"]] == TRUE)){
        
        resPlot = resPlot + 
        scale_x_continuous(trans = cytotidyr:::asinh_trans(cofactor = 5), breaks = c(0, 1, 3, 5, 10, 100, 1000, 5000))
            #scale_x_continuous(trans = 'asinh', breaks = c(0, 10, 100, 1000, 5000)) +
           
    }

    return(resPlot)
}

adt_feature$batch = adt_feature$subset

igm1 = plot_adt_density(
cell_x_adt = cell_x_adt, 
cell_x_feature = adt_feature, 
adt_marker_list = "IgM", 
unit = "sample",
parameter_list = list(method_label = "ADTnorm Recommended", arcsinhTransform = TRUE, bw = 0.02)
) + scale_fill_manual(breaks = levels(adt_feature$batch), values = colorRampPalette(RColorBrewer::brewer.pal(3, "Dark2"))(length(unique(adt_feature$batch)))) +  coord_cartesian(xlim = c(0.5, 5)) 

igm2 = plot_adt_density(
cell_x_adt = res_norm, 
cell_x_feature = adt_feature, 
adt_marker_list = "IgM", 
unit = "sample",
parameter_list = list(method_label = "ADTnorm Recommended", arcsinhTransform = TRUE, bw = 0.02)
) + scale_fill_manual(breaks = levels(adt_feature$batch), values = colorRampPalette(RColorBrewer::brewer.pal(3, "Dark2"))(length(unique(adt_feature$batch)))) +  coord_cartesian(xlim = c(0.5, 4)) 

ggarrange(igm1, igm2, ncol = 2, nrow = 1)
ggsave(paste0(fig_path, "density_plot_igm.pdf"), width = 20, height = 30)

den = print(plot_adt_density(
cell_x_adt = cell_x_adt, 
cell_x_feature = adt_feature, 
adt_marker_list = adt_marker_select[which(!(adt_marker_select %in% colnames(adt_data_final)))][2:12], 
unit = "sample",
parameter_list = list(method_label = "ADTnorm Recommended", arcsinhTransform = TRUE, bw = 0.0083)
) + scale_fill_manual(breaks = levels(adt_feature$batch), values = colorRampPalette(RColorBrewer::brewer.pal(3, "Dark2"))(length(unique(adt_feature$batch)))) +  coord_cartesian(xlim = c(0.75, 1.25))
)

den
ggsave(paste0(fig_path, "density_plot_adt_markers_not_in_paper.pdf"), width = 35, height = 30)



saveRDS(list(igm1, igm2, p_allTitration, p_filterTitration, den), file = paste0(out_path, "/RDS/figs_obj_igm1_igm2_p_allTitration_p_filterTitration_", run_name, ".rds"))



















## read in h5ad file
adata = readH5AD(paste0(in_path, "/CITE-Seq_Titration_Datasets/adata_combined_rna_adt_annotated-titration.h5ad"))
adata %>% rownames %>% tail(400)
## 274 Protein markers
colData(adata) %>% head
colData(adata)[, c("Capture", "Sample")] %>% unique

> colData(adata)[, c("Concentration", "Sample")] %>% table

# Concentration BF32-CD271 BF32-CD34 BF32-TNC BM32-CD34 BM32-TNC WF32-CD271
#          0.25       4269      5213     1898      6826     1628       4839
#          0.5        3921      5136     1975      7099     1521       5189
#          1.0        3981      5063     1654      6686     1648       5062
#          1.5        3135      3929        0      5370        0       3924
#          2.0        3879      4862     1961      6469     1489       4989
#          4.0        3843      5424     1894      7197     1621       3869
#              Sample
# Concentration WF32-CD34 WF32-TNC WM22-CD271 WM22-CD34 WM22-TNC
#          0.25      7007     1585       4345      8755     1304
#          0.5       6840     1533       4438      8691     1186
#          1.0       6722     1528       4738      8306     1375
#          1.5       5066        0       3379      6386        0
#          2.0       6735     1501       4181      8506     1464
#          4.0       7114     1127       4327      9104     1739
## That is why they dropped the 1.5x concentration.
colData(adata)$Donor_Capture = paste0(colData(adata)$Donor, "_", colData(adata)$Capture)
colData(adata)[, c("Concentration", "Donor_Capture")] %>% table


adt_index = 36602:36876
cell_index = which(adata$Concentration != 1.5)
## totalVI count
adt_data = adata@assays@data$X[adt_index, ]



adata_raw = readH5AD(paste0(in_path, "/CITE-Seq_Titration_Datasets/adata_combined_rna_titration-raw-counts.h5ad"))
colData(adata_raw)[, c("Concentration", "Sample")] %>% table

aml_obj = Read10X_h5(paste0(in_path, "/GEO/GSE245108_AML-7_2_filtered_feature_bc_matrix.h5"))

# Read specific data from the file
# Replace "group/dataset_name" with the actual path to the dataset
data <- h5read(h5_file, "group/dataset_name")


