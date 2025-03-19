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
library(ADTnorm)
# p_load(flowStats, flowCore, FlowSOM, ncdfFlow, flowViz, pdfCluster, cluster)

## =====================
## study name and paths
## =====================
run_name = "titration_Felix2022"
master_path = "./"
out_path = paste0(master_path, "manuscript/results/", run_name)
# fcs_path = paste0(out_path, "/FCS/")
fig_path = paste0(out_path, "/Figures/", run_name)

## ======================================
## Clean and organize titration datasets
## ======================================
## specific data to be analyzed
cellranger_matrix = Read10X(data.dir = paste0(out_path, "/Data/"))
# sce_obj = CreateSeuratObject(counts = cellranger_matrix$`Gene Expression`)
# sce_obj[['ADT']] = CreateAssayObject(counts = cellranger_matrix$`Antibody Capture`)

## load data directly
adt_data_full = cellranger_matrix$`Antibody Capture` %>% t %>% as.matrix
colnames(adt_data_full) = colnames(adt_data_full) %>% gsub("\\.1", "", .) %>% gsub("_Recombinant", "", .) %>% gsub("_", "", .) %>% gsub("-", "", .) %>% gsub("\\.", "", .)

marker_res = fread(paste0(out_path, "/Data/Supp5.csv"), header = TRUE)
optim_res = marker_res[, 1:2]
colnames(optim_res) = c("marker", "optim")
optim_res$marker = optim_res$marker %>% gsub("\\.", "", .) %>% gsub("-", "", .) %>% gsub("\\s", "", .)

remove_marker = colnames(adt_data_full)[which(!(colnames(adt_data_full) %in% optim_res$marker))]
# optim_res$marker[which(!(optim_res$marker %in% colnames(adt_data_full)))]
remove_marker = c(remove_marker, optim_res$marker[optim_res$optim == "not detectable"])

adt_feature = fread(paste0(out_path, "/Data/metadata.csv")) %>% data.frame
file_list = c("0.04X", "0.2X", "1X", "2X")
adt_feature$sample = factor(adt_feature$Calls, levels = file_list)
adt_feature$batch = factor(adt_feature$Calls, levels = file_list)
cellid = adt_feature$V1 %>% gsub(".*_", "", .) %>% gsub("\\.1", "-1", .) 

adt_data = adt_data_full[cellid, which(!(colnames(adt_data_full) %in% remove_marker))]
marker_list = colnames(adt_data) %>% sort #c("CD1c", "CD1d", "CD2", "CD3", "CD4", "CD5", "CD7", "CD8", "CD10", "CD11a", "CD11b", "CD11c", "CD14", "CD15", "CD16", "CD18", "CD19",  
# "CD20", "CD21", "CD22", "CD23", "CD24", "CD25", "CD26", "CD27", "CD28", "CD29", "CD30", "CD31", "CD32", "CD33", "CD35", "CD36",  "CD38", "CD39", "CD40",
# "CD41", "CD44", "CD45", "CD45RA", "CD45RO", "CD47", "CD49b", "CD49d", "CD49f", "CD52", "CD54", "CD56", "CD57", "CD58", "CD62L", "CD62P", "CD64", 
# "CD66ace", "CD66b", "CD69",  "CD70", "CD73", "CD79b", "CD81", "CD82", "CD85j", "CD86", "CD88", "CD94", "CD95", "CD96", "CD98", "CD99",  
# "CD101", "CD107a",  "CD122", "CD123", "CD127", "CD137", "CD141", "CD150", "CD154", "CD155", "CD158", "CD158b", "CD158e1", "CD158f",  
# "CD161", "CD163", "CD183", "CD185", "CD194", "CD195", "CD196", "CD197", "CD223", "CD224", "CD226",  "CD244", "CD254", "CD267", "CD268",  "CD272", 
# "CD278", "CD279",  "CD305", "CD307d", "CD314", "CD319", "CD328", "CD335", "CD337", "CD360", "CD366", 
# "CLEC12A", "CX3CR1", "GARP", "HLADR", "IgA", "IgD", "IgGFc", "IgM", "Iglightchainkappa", "Iglightchainlambda", "KLRG1", "TCRab", "TIGIT", "integrinB7") # colnames(adt_data) %>% sort

head(adt_data)
head(adt_feature)
saveRDS(adt_data, file = paste0(out_path, "/RDS/cell_x_adt_", run_name, ".rds"))
saveRDS(adt_feature, file = paste0(out_path, "/RDS/cell_x_feature_", run_name, ".rds"))

## Baseline Methods
arcsinh_param <- list(a = 1, b = 1 / 5, c = 0)
normalized_method <- list(
  Arcsinh = list(f = arcsinh_basic_transform, start_adt_method = "RawCount", param = arcsinh_param),
  CLR = list(f = clr_transform, start_adt_method = "RawCount", param = arcsinh_param),
  CLR2 = list(f = clr2_transform, start_adt_method = "RawCount", param = arcsinh_param),
  logCPM = list(f = log1pcpm_transform, start_adt_method = "RawCount"),
  ADTnorm_sample_default = list(
    f = ADTnorm_transform, start_adt_method = "RawCount", 
    param = list(
      unit = "sample",
      save_outpath = out_path,
      study_name = run_name,
      trimodal_marker = c("CD4"),
      lower_peak_thres = 0.001,
      bw_smallest_bi = 0.8,
      bw_smallest_tri = 0.8,
      bw_smallest_adjustments = list(CD3 = 0.8, CD4 = 0.8, CD8 = 0.8),
      shoulder_valley_slope = -1
    )
  ),
  Harmony_Arcsinh_sample = list(f = harmony_transform, start_adt_method = "Arcsinh", param = list(batch_column = "sample")),
  Harmony_CLR_sample = list(f = harmony_transform, start_adt_method = "CLR", param = list(batch_column = "sample"))
)

## run normalization
baseline_methods <- c("Arcsinh", "CLR", "CLR2", "logCPM", "ADTnorm_sample_default", "ADTnorm_sample_manual", "Harmony_Arcsinh_sample", "Harmony_CLR_sample") 

for(method in baseline_methods) {
  print(method)
  filename <- paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds")
  norm_f <- normalized_method[[method]][["f"]]
  start_adt_method <- normalized_method[[method]][["start_adt_method"]]
  method_param <- normalized_method[[method]][["param"]]

  if(start_adt_method == "RawCount"){
    cell_x_adt <- adt_data
  }else{
    cell_x_adt <- readRDS(file = paste0(out_path, "/RDS/adt_data_norm_", start_adt_method, "_", run_name, ".rds"))
  }
  if(stringr::str_detect(method, "CytofRUV")){
    cell_x_adt_colnames = colnames(cell_x_adt)
    colnames(cell_x_adt) = paste0("tmpName", 1:ncol(cell_x_adt))
    setwd("~")
  }
  
  ans <- do.call(
    norm_f,
    list(cell_x_adt = cell_x_adt, cell_x_feature = adt_feature, parameter_list = method_param)
  )
  if(stringr::str_detect(method, "CytofRUV")){
    colnames(ans) = cell_x_adt_colnames
  }
   
  saveRDS(ans, file = filename)

}

plot_adt_density <- function(cell_x_adt, cell_x_feature, adt_marker_list, unit, parameter_list = NULL, ncol = 6){
    
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
      tmpProfile <- cell_x_adt %>% data.frame %>%
        dplyr::select(all_of(adt_marker_list)) %>%
        data.frame() %>%
        gather(key = "ADT", value = "counts") %>%
        mutate(
            sample = rep(cell_x_feature$sample, length(adt_marker_list)),
            batch = rep(cell_x_feature$batch, length(adt_marker_list))
        ) 
    }else if(unit == "sample"){
      tmpProfile <- cell_x_adt %>% data.frame %>% 
        dplyr::select(all_of(adt_marker_list)) %>%
        data.frame() %>%
        gather(key = "ADT", value = "counts") %>%
        mutate(
            sample = rep(cell_x_feature$sample, length(adt_marker_list)),
            batch = rep(cell_x_feature$batch, length(adt_marker_list))
        ) 
    }
    fillColor <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(unique(tmpProfile$batch)))
    resPlot <- ggplot(tmpProfile, aes(x = counts, y = sample, fill = batch)) +
        geom_density_ridges(bandwidth = bw, scale = 2.5) +
        facet_wrap( ~ factor(ADT, levels = adt_marker_list), scales = "free", ncol = ncol) + #, space = "free"
        theme_bw(base_size = 50) +
        xlab(method_label) +
        ylab("") +
        ggpubr::rotate_x_text(angle = 90) +
        scale_fill_manual(values = fillColor) +
        rremove("legend.title") +
        theme(
        legend.position = "top",
        axis.text.y = element_blank(), # Remove y-axis labels
        panel.spacing = unit(1, "lines") # Reduce the gap between panels, adjust the unit value as needed
        )

    if(("arcsinhTransform" %in% names(parameter_list)) && (parameter_list[["arcsinhTransform"]] == TRUE)){
        
        resPlot = resPlot + 
        scale_x_continuous(trans = cytotidyr:::asinh_trans(cofactor = 5), breaks = c(0, 10, 100, 1000, 5000)) 
    }

    return(resPlot)
}

bw_list = list("Arcsinh" = 0.08, "CLR" = 0.1, "CLR2" = 0.1, "logCPM" = 0.08,  "ADTnorm_sample_default" = 0.1, "ADTnorm_sample_manual" = 0.2, "Harmony_Arcsinh_sample" = 0.1, "Harmony_CLR_sample" = 0.1)
trans_list = list("Arcsinh" = FALSE, "CLR" = FALSE, "CLR2" = FALSE, "logCPM" = TRUE, "ADTnorm_sample_default" = FALSE, "ADTnorm_sample_manual" = FALSE, "Harmony_Arcsinh_sample" = FALSE, "Harmony_CLR_sample" = FALSE)
baseline_methods = c("Arcsinh", "CLR", "CLR2", "logCPM",  "ADTnorm_sample_default", "ADTnorm_sample_manual", "Harmony_Arcsinh_sample", "Harmony_CLR_sample") 

for(method in baseline_methods){
  print(method)
  filename <- paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds")
  ans = readRDS(file = filename)

  fig_height = 150
  pdf(paste0(out_path, "/Figures/", method, "_sample_", run_name, "_adt_density.pdf"), width = 80, height = fig_height)
  print(plot_adt_density(
  cell_x_adt = ans %>% data.frame, 
  cell_x_feature = adt_feature, 
  adt_marker_list = marker_list, 
  unit = "sample",
  parameter_list = list(method_label = method, arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
  ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))))
  dev.off()

  
}

## summary marker that failed and worked
method = "ADTnorm_sample_manual"
filename <- paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds")
ans = readRDS(file = filename)

marker_category = "mixed"
marker_list_fail = c("CD2", "CD7", "CD26", "CD27", "CD28", "CD52", "CD82", "CD127", "CD122", "CD224", "CD244", "CD305", "CD328", "TCRab")
marker_list_mixed = c("CD11b", "CD31", "CD39", "CD40", "CD54", "CD73", "CD94",  "CD101", "CD158b", "CD194", "CD335",  "KLRG1")
marker_list_success = c("CD3", "CD4", "CD5", "CD8", "CD11c", "CD16", "CD19", "CD21", "CD32", "CD33", "CD35", "CD36", "CD45RA", "CD45RO", "CD49f", "CD56", "CD58", "CD62L", "CD69", "CD95", "CD158", "CD158e1", "CD226", "CD328", "HLADR", "IgD", "IgM", "IgGFc", "Iglightchainkappa", "Iglightchainlambda")

fig_height = 45
pdf(paste0(out_path, "/Figures/", method, "_sample_", run_name, "_adt_density_", marker_category, ".pdf"), width = 80, height = fig_height)
print(plot_adt_density(
cell_x_adt = ans %>% data.frame, 
cell_x_feature = adt_feature, 
adt_marker_list = marker_list_mixed, 
unit = "sample",
parameter_list = list(method_label = method, arcsinhTransform = trans_list[[method]], bw = bw_list[[method]]),
ncol = 4
) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))))
dev.off()

# ## add refined annotation
method = "ADTnorm_sample_manual"
filename <- paste0(out_path, "/RDS/adt_data_norm_ADTnorm_sample_manual_", run_name, ".rds")
res_norm = readRDS(file = filename)
adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_withRefinedAnnotation_ADTnorm_sample_manual_titration_Felix2022.rds"))

## UMAP
generate_umap_titration = function(ans_umap, adt_feature, method, run_name, npc){
  pdf(paste0(out_path, "/Figures/adt_", method, "_", run_name, "_umap_raster_cosine_uwot", npc, ".pdf"), width = 7, height = 7)

  print(plot_umap(ans_umap, adt_feature, point_size = 1, parameter_list = list(
  target_feature = "sample", 
  color_design = colorRampPalette(brewer.pal(8, "Dark2"))(length(adt_feature$sample %>% unique)),
  method_label = method
  )) #+ scale_color_brewer(palette = "Dark2")
  )
  print(plot_umap(ans_umap, adt_feature, point_size = 1, parameter_list = list(
  target_feature = "annotation", 
  color_design = c("#E41A1C", 
    "#4C71A4", "#47A265", "#934c93", 
    "#CB6651", "#FFAF13", 
    "#E8D430", "#B05B3A", "#F781BF", "grey"),
    #color_design = c(colorRampPalette(brewer.pal(8, "Dark2"))(14), "grey"),
  color_break = c("B",
    "Naive CD4", "Memory CD4", "Tissue-resident Memory CD4", 
    "Naive CD8", "Memory CD8", 
    "Monocytes", "CD16+ NK", "CD16- NK", "undefined"),
    method_label = method
  )))
  adt_feature$type = factor(adt_feature$type, levels = c("B", "CD4T", "CD8T", "CM", "NK"), labels = c("B", "CD4 T", "CD8 T", "Monocytes", "NK"))
  print(plot_umap(ans_umap, adt_feature, point_size = 1, parameter_list = list(
  target_feature = "type", 
  color_design = c("#B00000", "#3F918B", "#896191", "#FF980A", "#226fa7", "#F781BF", "grey"), 
  color_break = c("B", "CD4 T", "CD8 T", "Monocytes", "NK", "DCs", "undefined"),
  method_label = method
  )))

  dev.off()
}

set.seed(20240129)
method = "Harmony_Arcsinh_sample"
filename <- paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds")
res_norm = readRDS(file = filename)
index = which(!(colnames(adt_data) %in% c("CD7", "CD2", "CD26", "CD27", "CD28", "CD52", "CD82", "CD127", "CD305", "TCRab")))
for(npc in c(10, 15, 30, 50)){  
    ans_pca = prcomp(res_norm[, index])$x 
    ans_umap = uwot::umap(ans_pca[, 1:min(ncol(res_norm), npc)], metric = "cosine", seed = "20240129") %>% data.frame
    # n_neighbors = 20, min_dist = 0.001,
    generate_umap_titration(ans_umap, adt_feature, method, run_name, npc = npc)  
}

## ================
## staining quality
## ================
cell_x_adt = arcsinh_transform(cell_x_adt = adt_data) 
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


peak_sep_summary$batch = factor(peak_sep_summary$batch, levels = file_list)
peak_sep_summary$sample = factor(peak_sep_summary$sample, levels = file_list)
peak_sep_summary$adt_marker = factor(peak_sep_summary$adt_marker, levels = marker_list)
peak_sep_summary$peak_num[which(peak_sep_summary$adt_marker == "CD49f")] = "# of peak: 2"

peak2marker = peak_sep_summary$adt_marker[which(peak_sep_summary$peak_num == "# of peak: 2")] %>% table
peak2marker = which(peak2marker > 1) %>% names
list4 = c("CD36")
list2 = c("CD16", "CD56", "CD62L", "CD33", "CD45RA", "CD49f",  "CD158", "CD158e1", "HLADR", "IgD", "IgM")
listhigh = c("CD26", "CD45RO", "CD127", "CD226", "CD305", "IgGFc")
listsame = c("CD11b", "CD21", "CD35", "CD45RO", "CD58", "CD95", "CD226", "IgGFc")

peak_sep_select = peak_sep_summary %>% dplyr::filter(adt_marker %in% peak2marker)
peak_sep_select$category = "1X"
peak_sep_select$category[which(peak_sep_select$adt_marker %in% list4)] = "0.04X"
peak_sep_select$category[which(peak_sep_select$adt_marker %in% list2)] = "0.2X"
peak_sep_select$category[which(peak_sep_select$adt_marker %in% listhigh)] = "2X+"
peak_sep_select$category[which(peak_sep_select$adt_marker %in% listsame)] = "Same"

peak_sep_select$decision = "notSelected"
peak_sep_select$decision[which(peak_sep_select$adt_marker %in% list4 & peak_sep_select$sample == "0.04X")] = "Select"
peak_sep_select$decision[which(peak_sep_select$adt_marker %in% list2 & peak_sep_select$sample == "0.2X")] = "Select"
peak_sep_select$decision[which(peak_sep_select$category == "1X" & peak_sep_select$sample == "1X")] = "Select"
peak_sep_select$decision[which(peak_sep_select$category == "2X+" & peak_sep_select$sample == "2X")] = "Select"
# peak_sep_select$decision[which(peak_sep_select$category == "Same")] = "Select"

peak_sep_select$batch = factor(peak_sep_select$batch, levels = file_list)

pdf(paste0(out_path, "/Figures/marker_stain_quality_score_concentration_recommendation.pdf"), width = 35, height = 15)
peak_sep_select %>%
  ggplot(aes(x = adt_marker, y = sep_power, fill = batch)) +
  geom_bar(stat = "identity", position = position_dodge(), aes(alpha = decision)) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "navy") +
  facet_grid(~category, scales = "free_x", space = "free") +
  theme_bw(base_size = 25) +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(peak_sep_select$batch)))) +
  scale_alpha_discrete(range = c(0.45, 1)) +  # Adjust as necessary for your data
  ylab("Separation Power") +
  xlab("") +
  theme(legend.position = "top") +
  rotate_x_text(angle = 45)
dev.off()

pdf(paste0(out_path, "/Figures/marker_staining_peak_distance_sd_deep_auc_withMulti_withOne_eachsample_PeakNum.pdf"), width = 45, height = 15)
peak_sep_summary %>% ggplot(aes(x = adt_marker, y = sep_power, color = batch, fill = batch)) +
geom_bar(stat = "identity", position = position_dodge2(width = 0.7, preserve = "single")) +
facet_grid(~peak_num, scales = "free_x", space = "free") +
theme_bw(base_size = 25) +
scale_color_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))) +
scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))) +
ylab("Separation Power") +
xlab("") +
theme(legend.position = "top") +
rotate_x_text(angle = 45)
dev.off()

## marker that failed at low concentrations: marker_list_fail
pdf(paste0(out_path, "/Figures/marker_stain_quality_score_failed.pdf"), width = 15, height = 8)

peak_sep_select %>% dplyr::filter(adt_marker %in% marker_list_fail) %>% 
  ggplot(aes(x = adt_marker, y = sep_power, fill = batch)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "navy") +
  facet_grid(~category, scales = "free_x", space = "free") +
  theme_bw(base_size = 25) +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(peak_sep_select$batch)))) +
  scale_alpha_discrete(range = c(0.45, 1)) +  # Adjust as necessary for your data
  ylab("Separation Power") +
  xlab("") +
  theme(legend.position = "top") +
  rotate_x_text(angle = 45)
  
dev.off()

## ===========
## auto-gating
## ===========
cell_x_adt = arcsinh_transform(cell_x_adt = adt_data) 
cell_x_feature = adt_feature
rownames(cell_x_adt) = paste0("Cell", 1:nrow(cell_x_adt))
rownames(cell_x_feature) = paste0("Cell", 1:nrow(cell_x_feature))

pheno_matrix <- matrix("-", nrow(cell_x_adt), ncol(cell_x_adt))
colnames(pheno_matrix) <- colnames(cell_x_adt)
rownames(pheno_matrix) <- rownames(cell_x_adt)

for(adt_marker_select in colnames(cell_x_adt)){
    peak_valley_density = readRDS(paste0(out_path, "/RDS/peak_valley_locations_", adt_marker_select, "_", run_name, ".rds"))
   
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

broad_label = c("B", "CD4 T", "CD8 T", "Monocytes", "NK")
refined_label = c("Naive CD4", "Memory CD4", "Tissue-resident Memory CD4", "Naive CD8", "Memory CD8", "CD16+ NK", "CD16- NK")

adt_feature$cell_type_l1 = factor(adt_feature$type, levels = c("B", "CD4T", "CD8T", "CM", "NK"), labels = broad_label)
adt_feature$cell_type_l2 = factor(adt_feature$annotation, levels = c("B", "Naive CD4", "Memory CD4", "Tissue-resident Memory CD4", "Naive CD8", "Memory CD8", "Monocytes", "CD16+ NK", "CD16- NK"))
cell_x_feature = adt_feature
rownames(cell_x_feature) = paste0("Cell", rownames(cell_x_feature))

FeaturePlot(adt_obj, order = TRUE, features = c("CD16"), cols = c("yellow", "black", "purple"))

auto_gating = list()
## broad cell types
auto_gating[["B"]] = pheno_matrix %>% dplyr::filter(CD3 == "-", CD19 == "+") %>% rownames
auto_gating[["CD4 T"]] = pheno_matrix %>% dplyr::filter(CD3 == "+", CD19 == "-", CD4 == "+", CD8 == "-") %>% rownames
auto_gating[["CD8 T"]] = pheno_matrix %>% dplyr::filter(CD3 == "+", CD19 == "-", CD4 == "-", CD8 == "+") %>% rownames
auto_gating[["Monocytes"]] = c(pheno_matrix %>% dplyr::filter(CD3 == "-", CD19 == "-", CD14 == "+", CD16 == "-") %>% rownames, pheno_matrix %>% dplyr::filter(CD3 == "-", CD19 == "-", CD11c == "+", CD16 == "-") %>% rownames) %>% unique
auto_gating[["NK"]] = c(pheno_matrix %>% dplyr::filter(CD3 == "-", CD19 == "-", CD14 == "-", CD56 == "+") %>% rownames, pheno_matrix %>% dplyr::filter(CD3 == "-", CD19 == "-", CD14 == "-", CD16 == "+") %>% rownames) %>% unique

## refined cell types
auto_gating[["Naive CD4"]] = pheno_matrix %>% dplyr::filter(CD3 == "+", CD19 == "-", CD4 == "+", CD8 == "-", CD45RA == "+", CD45RO == "-") %>% rownames
auto_gating[["Memory CD4"]] = pheno_matrix %>% dplyr::filter(CD3 == "+", CD19 == "-", CD4 == "+", CD8 == "-", CD45RA == "-", CD45RO == "+", CD49f == "-") %>% rownames
auto_gating[["Tissue-resident Memory CD4"]] = pheno_matrix %>% dplyr::filter(CD3 == "+", CD19 == "-", CD4 == "+", CD8 == "-", CD45RA == "-", CD45RO == "+", CD49f == "+") %>% rownames
auto_gating[["Naive CD8"]] = pheno_matrix %>% dplyr::filter(CD3 == "+", CD19 == "-", CD4 == "-", CD8 == "+", CD45RA == "+", CD45RO == "-") %>% rownames
auto_gating[["Memory CD8"]] = pheno_matrix %>% dplyr::filter(CD3 == "+", CD19 == "-", CD4 == "-", CD8 == "+", CD45RA == "-", CD45RO == "+") %>% rownames
auto_gating[["CD16+ NK"]] = pheno_matrix %>% dplyr::filter(CD3 == "-", CD19 == "-", CD14 == "-", CD56 == "+", HLADR == "-", CD16 == "+") %>% rownames
auto_gating[["CD16- NK"]] = pheno_matrix %>% dplyr::filter(CD3 == "-", CD19 == "-", CD14 == "-", CD56 == "+", HLADR == "-", CD16 == "-") %>% rownames

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
        cell_x_feature[auto_gating[[cT]], c("batch", "sample")] %>% group_by(batch, sample) %>% summarize(freq = n()),
        cell_x_feature[overlap_cell_name, c("batch", "sample")] %>% group_by(batch, sample) %>% summarize(freq = n()),
        by = c("batch", "sample")
    ) %>% dplyr::filter(freq.x >= 5, freq.y >=10) %>% mutate(accuracy = freq.y/freq.x * 100, cT = cT, cell_type_label = "Broad Labeling") %>% rbind(acc_each_sample, .)

}
for(cT in refined_label){
    overlap_ind <- which(auto_gating[[cT]] %in% manual_gating[[cT]])
    overlap_cell_name <- auto_gating[[cT]][overlap_ind]
    cell_x_feature$auto_gating_l2[which(rownames(cell_x_feature) %in% overlap_cell_name)] = cT
    cell_x_feature$auto_gating_mixed[which(rownames(cell_x_feature) %in% overlap_cell_name)] = cT

    acc_each_sample <- left_join(
        cell_x_feature[auto_gating[[cT]], c("batch", "sample")] %>% group_by(batch, sample) %>% summarize(freq = n()),
        cell_x_feature[overlap_cell_name, c("batch", "sample")] %>% group_by(batch, sample) %>% summarize(freq = n()),
        by = c("batch", "sample")
    ) %>% dplyr::filter(freq.x >= 5, freq.y >=10) %>% mutate(accuracy = freq.y/freq.x * 100, cT = cT, cell_type_label = "Refined Labeling") %>% rbind(acc_each_sample, .)

}

acc_each_sample_df = rbind(
    cell_x_feature %>% dplyr::filter(cell_type_l1 %in% broad_label) %>% group_by(sample, cell_type_l1) %>% summarize(cell_num_manual_label = n()) %>% mutate(cT = cell_type_l1) %>% dplyr::select(sample, cT, cell_num_manual_label) %>% data.frame,
    cell_x_feature %>% dplyr::filter(cell_type_l2 %in% refined_label) %>% group_by(sample, cell_type_l2) %>% summarize(cell_num_manual_label = n()) %>% mutate(cT = cell_type_l2) %>% dplyr::select(sample, cT, cell_num_manual_label) %>% data.frame
) %>% left_join(acc_each_sample, ., by = c("sample", "cT")) %>% data.frame

acc_each_sample_df$cT = factor(acc_each_sample_df$cT, levels = c("B", "CD4 T", "CD8 T", "Monocytes", "NK", "Naive CD4", "Memory CD4", 
"Tissue-resident Memory CD4", "Naive CD8", "Memory CD8", "CD16+ NK", "CD16- NK"))
cell_x_feature$auto_gating_l1 = factor(cell_x_feature$auto_gating_l1, levels = c(broad_label, "undefined"))
cell_x_feature$auto_gating_l2 = factor(cell_x_feature$auto_gating_l2, levels = c(refined_label, "undefined"))

acc_each_sample_df %>% dplyr::filter(cell_num_manual_label > 10) %>%  
ggplot(aes(x = cT, y = accuracy)) +
geom_boxplot(outlier.shape = NA, position = position_dodge2(preserve = "single")) +
geom_jitter(size = 5, position = position_dodge(width = 0.4), aes(color = batch)) +
facet_grid(~cell_type_label, scale = "free", space = "free") +
xlab("") +
ylab("Auto-gating Accuracy (100%)") +
theme_bw(base_size = 20) +
scale_color_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))) +
rotate_x_text(angle = 45) +
theme(legend.position = "top") +
rremove("legend.title") +
ggsave(filename = paste0(out_path, "/Figures/auto_gating_phenotype_accuracy.pdf"), width = 10, height = 10)


## umap with auto-gating label
method = "ADTnorm_sample_manual"
filename <- paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds")
res_norm = readRDS(file = filename)
index = which(!(colnames(adt_data) %in% c("CD7", "CD2", "CD26", "CD27", "CD28", "CD52", "CD82", "CD127", "CD305", "TCRab")))
for(npc in c(10, 15, 30, 50)){  
    ans_pca = prcomp(res_norm[, index])$x 
    ans_umap = uwot::umap(ans_pca[, 1:min(ncol(res_norm), npc)], metric = "cosine", seed = "20240129") %>% data.frame
    pdf(paste0(out_path, "/Figures/umap_auto_gating_", method, "_", run_name, "_npc_", npc, ".pdf"), width = 10, height = 10)
        ## broad label
        reindex1 = which(cell_x_feature$auto_gating_l1 == "undefined")
        reindex2 = which(cell_x_feature$auto_gating_l1 != "undefined")
        reindex = c(reindex1, reindex2)
        cell_x_feature$auto_gating_l1_label = factor(cell_x_feature$auto_gating_l1, levels = c(broad_label, "undefined"), labels = c("B (CD19+)", "CD4 T (CD3+CD4+CD8-)", "CD8 T (CD3+CD4-CD8+)", "Monocytes (CD3-CD19-CD14+/CD11c+)", "NK (CD3-CD19-CD14-CD56+)", "undefined"))
        print(plot_umap(ans_umap[reindex, ], cell_x_feature[reindex, ], point_size = 1, parameter_list = list(
        target_feature = "auto_gating_l1_label", 
        color_design = c("#B00000", "#3F918B", "#896191", "#FF980A", "#226fa7", "grey"), 
        color_break = c("B (CD19+)", "CD4 T (CD3+CD4+CD8-)", "CD8 T (CD3+CD4-CD8+)", "Monocytes (CD3-CD19-CD14+/CD11c+)", "NK (CD3-CD19-CD14-CD56+)", "undefined"),
        method_label = "Auto-gating"
        )))

        ## refined label
        reindex1 = which(cell_x_feature$auto_gating_mixed == "undefined")
        reindex2 = which(cell_x_feature$auto_gating_mixed != "undefined" & !(cell_x_feature$auto_gating_mixed %in% c("NK", "CD4 T", "CD8 T")))
        reindex = c(reindex1, reindex2)
        cell_x_feature$auto_gating_mixed_label = factor(cell_x_feature$auto_gating_mixed, levels = c("B",
            "Naive CD4", "Memory CD4", "Tissue-resident Memory CD4", 
            "Naive CD8", "Memory CD8", 
            "Monocytes", "CD16+ NK", "CD16- NK", "undefined"), 
            labels = c("B (CD19+)",
            "Naive CD4 (CD4+CD45RA+CD45RO-)", "Memory CD4 (CD4+CD45RA-CD45RO+CD49f-)", "Tissue-resident Memory CD4 (CD4+CD45RA-CD45RO+CD49f+)", 
            "Naive CD8 (CD8+CD45RA+CD45RO-)", "Memory CD8 (CD8+CD45RA-CD45RO+)", 
            "Monocytes (CD3-CD19-CD14+/CD11c+)", "CD16+ NK (CD56+CD16+)", "CD16- NK (CD56+CD16-)", "undefined"))
        print(plot_umap(ans_umap[reindex, ], cell_x_feature[reindex, ], point_size = 1, parameter_list = list(
        target_feature = "auto_gating_mixed_label",
        color_design = c("#E41A1C", 
        "#4C71A4", "#47A265", "#934c93", 
        "#CB6651", "#FFAF13", 
        "#E8D430", "#B05B3A", "#F781BF", "grey"),
        #color_design = c(colorRampPalette(brewer.pal(8, "Dark2"))(14), "grey"),
        color_break = c("B (CD19+)",
            "Naive CD4 (CD4+CD45RA+CD45RO-)", "Memory CD4 (CD4+CD45RA-CD45RO+CD49f-)", "Tissue-resident Memory CD4 (CD4+CD45RA-CD45RO+CD49f+)", 
            "Naive CD8 (CD8+CD45RA+CD45RO-)", "Memory CD8 (CD8+CD45RA-CD45RO+)", 
            "Monocytes (CD3-CD19-CD14+/CD11c+)", "CD16+ NK (CD56+CD16+)", "CD16- NK (CD56+CD16-)", "undefined"),
        method_label = "Auto-gating"
        ))+ guides(color = guide_legend(nrow = 4))) 

    dev.off()

}



  