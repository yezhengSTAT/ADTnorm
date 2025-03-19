library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(httpgd)
library(ggpubr)
library(splines)

run_name = "public13Dataset_CITEseq"

master_path = "./"
in_path = "./publicData_CITEseq/data/"
out_path = paste0(master_path, "manuscript/results/", run_name)
fcs_path = paste0(out_path, "/FCS/")
fig_path = paste0(out_path, "/Figures/")

marker_list = c("CD3", "CD4", "CD8",  "CD14", "CD19", "CD25", "CD45RA", "CD56", "CD127")
adt_data_full = readRDS(file = paste0(out_path, "/RDS/adt_data_full_RawCount_", run_name, ".rds"))
adt_data = readRDS(file = paste0(out_path, "/RDS/adt_data_RawCount_", run_name, ".rds"))
adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", run_name, ".rds"))

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
    y_valley = y[x_valley]
    auc_scaler = area_under_curve(density_res, threshold = valley_info[each_sample, 1], peak_num = peak_num)
    deep_scaler = (1 + y_peakR - y_valley) * (1 + auc_scaler)
    return(deep_scaler)
}

two_peak_stain_quality = function(cell_x_feature, cell_x_adt, adt_marker_select, each_sample, peak_info, valley_info, peak_num = 2){

    ind = which(cell_x_feature$sample == each_sample)
    mean_diff = (peak_info[each_sample, ncol(peak_info)] - peak_info[each_sample, 1])
    adt_tmp = cell_x_adt[ind, adt_marker_select]
    within_peak_sd = sqrt((sum((adt_tmp[which(adt_tmp < valley_info[each_sample, 1])] - peak_info[each_sample, 1])^2) + sum((adt_tmp[which(adt_tmp > valley_info[each_sample, 1])] - peak_info[each_sample, ncol(peak_info)])^2))/length(adt_tmp))
    
    deep_scaler = valley_deep_scaler(cell_x_feature, cell_x_adt, adt_marker_select, each_sample, peak_info, valley_info, peak_num = peak_num)
    # return(mean_diff/within_peak_sd)
    return((mean_diff*deep_scaler)/within_peak_sd)
}

multi_peak_stain_quality = function(cell_x_feature, cell_x_adt, adt_marker_select, each_sample, peak_info, valley_info, peak_num = 3){

    ind = which(cell_x_feature$sample == each_sample)
    mean_diff = (peak_info[each_sample, ncol(peak_info)] - peak_info[each_sample, 1])
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
    
    peak_valley_density = readRDS(paste0(out_path, "/RDS/peak_valley_locations_", adt_marker_select, "_ADTnorm_sample_manual_keepZero.rds"))
    
    peak_info = peak_valley_density$peak_landmark_list
    valley_info = peak_valley_density$valley_landmark_list

    for(each_sample in unique(adt_feature$sample)){
        batch_info = unique(adt_feature$study_name[which(adt_feature$sample == each_sample)]) %>% as.vector
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

file_list <- c(
    "10X_pbmc_10k", "10X_pbmc_1k", "10X_pbmc_5k_v3", "10X_pbmc_5k_nextgem", "10X_malt_10k", 
    "stuart_2019", "granja_2019_pbmc", "granja_2019_bmmc", "hao_2020",
    "kotliarov_2020", "witkowski_2020", "triana_2021", "buus_2021_T"   
)

peak_sep_summary$batch = factor(peak_sep_summary$batch, levels = file_list)
peak_sep_summary$sample = factor(peak_sep_summary$sample, levels = unique(cell_x_feature$sample))
