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
library(ADTnorm)


## =====================
## study name and paths
## =====================
run_name = "public13Dataset_CITEseq"

master_path = "./"
in_path = "./publicData_CITEseq/data/"
out_path = paste0(master_path, "manuscript/results/", run_name)
fcs_path = paste0(out_path, "/FCS/")
fig_path = paste0(out_path, "/Figures/")

## ======================================
## Clean and organize 13 public datasets
## ======================================
## specific data to be analyzed
file_list <- c(
    "10X_pbmc_10k", "10X_pbmc_1k", "10X_pbmc_5k_v3", "10X_pbmc_5k_nextgem", "10X_malt_10k", 
    "stuart_2019", "granja_2019_pbmc", "granja_2019_bmmc", "hao_2020",
    "kotliarov_2020", "witkowski_2020", "triana_2021", "buus_2021_T"   
)

## load data directly
marker_list = c("CD3", "CD4", "CD8",  "CD14", "CD19", "CD25", "CD45RA", "CD56", "CD127")


## ======================================
## Create T cell only data 
## ======================================
adt_data_full = readRDS(file = paste0(out_path, "/RDS/adt_data_full_RawCount_", run_name, ".rds"))
adt_data = readRDS(file = paste0(out_path, "/RDS/adt_data_RawCount_", run_name, ".rds"))
adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", run_name, ".rds"))
adt_feature_unique = adt_feature %>% select(sample, batch, sample_status, study_name) %>% unique
adt_feature_unique = adt_feature_unique[match(adt_feature$sample %>% levels, adt_feature_unique$sample),]

gex_data = readRDS(file = './results/publicData_CITEseq/RDS/rna_data_common_RawCount_public13Dataset_CITEseq.rds')

# cell_toremove = which(adt_feature$study_name %in% c("10X_pbmc_10k", "10X_pbmc_1k", "10X_pbmc_5k_v3", "10X_pbmc_5k_nextgem", "stuart_2019") & !(adt_feature$cell_type_l1 %in% c("CD4 T", "CD8 T")))
# cell_toremove = which(adt_feature$study_name %in% c("hao_2020", "triana_2021") & !(adt_feature$cell_type_l1 %in% c("CD4 T", "CD8 T")))
cell_toremove = c(
    which(adt_feature$study_name %in% c("hao_2020") & !(adt_feature$cell_type_l1 %in% c("CD4 T", "CD8 T"))),
    which(adt_feature$study_name %in% c("triana_2021") & !(adt_feature$cell_type_l1 %in% c("CD8 T")))
)
adt_data_t = adt_data[-cell_toremove, ]
adt_feature_t = adt_feature[-cell_toremove, ]
gex_data_t = gex_data[-cell_toremove, ]

run_name = "public13Dataset_CITEseq_Tcelllargeprop_cd8T"
out_path = paste0(master_path, "manuscript/results/", run_name)
fcs_path = paste0(out_path, "/FCS/")
fig_path = paste0(out_path, "/Figures/")
saveRDS(adt_data_t, file = paste0(out_path, "/RDS/adt_data_RawCount_", run_name, ".rds"))
saveRDS(adt_feature_t, file = paste0(out_path, "/RDS/adt_feature_", run_name, ".rds"))
saveRDS(gex_data_t, file = paste0(out_path, "/RDS/rna_data_RawCount_", run_name, ".rds"))


## ==================================
## Normalization Methods Processing
## ==================================
run_name = "public13Dataset_CITEseq_Tcelllargeprop_cd8T" #"public13Dataset_CITEseq_Tcelllargeprop" ## "public13Dataset_CITEseq_Tcelllargeprop_cd8T"
out_path = paste0(master_path, "manuscript/results/", run_name)
fcs_path = paste0(out_path, "/FCS/")
fig_path = paste0(out_path, "/Figures/")

adt_data = readRDS(file = paste0(out_path, "/RDS/adt_data_RawCount_", run_name, ".rds"))
adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", run_name, ".rds"))
adt_feature_unique = adt_feature %>% select(sample, batch, sample_status, study_name) %>% unique
adt_feature_unique = adt_feature_unique[match(adt_feature$sample %>% levels, adt_feature_unique$sample),]

positive_peak_study_Tcelllargeprop_cd8T = list(
  ADT = c(
    "CD3", "CD3", "CD3", "CD8"
  ),
  sample = c(
    "hao_2020", "triana_2021", "buus_2021_T", "triana_2021"
  )
)
sample_tmp = adt_feature %>% dplyr::filter(study_name %in% c("hao_2020", "triana_2021", "buus_2021_T")) %$% sample %>% unique %>% as.vector
sample_tmp2 = adt_feature %>% dplyr::filter(study_name == "triana_2021") %$% sample %>% unique %>% as.vector
positive_peak_sample_Tcelllargeprop_cd8T = list(
  ADT = c(
    rep("CD3", length(sample_tmp)), rep("CD8", length(sample_tmp2))
  ),
  sample = c(
    sample_tmp, sample_tmp2
  )
)
## Baseline Methods
arcsinh_param <- list(a = 1, b = 1 / 5, c = 0)
normalized_method <- list(
  Arcsinh = list(f = arcsinh_basic_transform, start_adt_method = "RawCount", param = arcsinh_param),
  CLR = list(f = clr_transform, start_adt_method = "RawCount", param = arcsinh_param),
  logCPM = list(f = log1pcpm_transform, start_adt_method = "RawCount"),
  Arcsinh_CLR = list(f = arcsinh_clr_transform, start_adt_method = "Arcsinh", param = arcsinh_param),
  CytofRUV_study = list(
    f = cytofruv_transform, 
    start_adt_method = "RawCount", 
    param = list(
      tmp_path = paste0(out_path, "/xlsx/CytofRUV_study/"), 
      fcs_path = fcs_path,
      batch = adt_feature_unique$study_name,
      condition = adt_feature_unique$study_name,
      patient_id = adt_feature_unique$study_name,
      clusters_nb = length(unique(adt_feature$cell_type_l1)), ## use broad clustering label number as the k mean clustering number
      start_adt_method = "RawCount")
  ),
  CytofRUV_sample = list(
    f = cytofruv_transform, 
    start_adt_method = "RawCount", 
    param = list(
      tmp_path = paste0(out_path, "/xlsx/CytofRUV_sample/"), 
      fcs_path = fcs_path,
      batch = adt_feature_unique$sample,
      condition = adt_feature_unique$sample_status,
      patient_id = adt_feature_unique$sample,
      clusters_nb = length(unique(adt_feature$cell_type_l1)), ## use broad clustering label number as the k mean clustering number
      start_adt_method = "RawCount")
  ),
  fastMNN_study = list(f = fastMNN_transform, start_adt_method = "RawCount", param = list(norm_unit = "study")),
  fastMNN_sample = list(f = fastMNN_transform, start_adt_method = "RawCount", param = list(norm_unit = "sample")),
  DSB = list(f = dsb_transform, start_adt_method = "RawCount", param = list(use_isotype_control = FALSE)),
  ADTnorm_sample = list(
    f = ADTnorm_transform, start_adt_method = "RawCount", 
    param = list(
      unit = "sample",
      save_outpath = out_path,
      study_name = run_name,
      trimodal_marker = c("CD4", "CD45RA"),
      positive_peak = positive_peak_sample_Tcelllargeprop_cd8T,
      lower_peak_thres = 0.05
    )
  ),
  ADTnorm_study = list(
    f = ADTnorm_transform, start_adt_method = "RawCount", 
    param = list(
      unit = "study",
      save_outpath = out_path,
      study_name = run_name,
      trimodal_marker = c("CD4", "CD45RA"),
      positive_peak = positive_peak_study_Tcelllargeprop_cd8T,
      lower_peak_thres = 0.05
    )
  )
  # Harmony_RawCount_study = list(f = harmony_transform, start_adt_method = "RawCount", param = list(batch_column = "study_name")), 
  # Harmony_RawCount_sample = list(f = harmony_transform, start_adt_method = "RawCount", param = list(batch_column = "sample")),
  # Harmony_Arcsinh_study = list(f = harmony_transform, start_adt_method = "Arcsinh", param = list(batch_column = "study_name")),
  # Harmony_Arcsinh_sample = list(f = harmony_transform, start_adt_method = "Arcsinh", param = list(batch_column = "sample")),
  # Harmony_CLR_study = list(f = harmony_transform, start_adt_method = "CLR", param = list(batch_column = "study_name")),
  # Harmony_CLR_sample = list(f = harmony_transform, start_adt_method = "CLR", param = list(batch_column = "sample")),
  # Harmony_logCPM_study = list(f = harmony_transform, start_adt_method = "logCPM", param = list(batch_column = "study_name")),
  # Harmony_logCPM_sample = list(f = harmony_transform, start_adt_method = "logCPM", param = list(batch_column = "sample"))
)

## run normalization
baseline_methods <- c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "CytofRUV_study", "CytofRUV_sample", "fastMNN_study", "fastMNN_sample", "DSB", "Harmony_RawCount_study", "Harmony_RawCount_sample", "Harmony_Arcsinh_study", "Harmony_Arcsinh_sample", "Harmony_CLR_study", "Harmony_CLR_sample", "Harmony_logCPM_study", "Harmony_logCPM_sample", "ADTnorm_sample", "ADTnorm_study", "totalVI_sample_GPU", "totalVI_sample_CPU", "totalVI_study_GPU", "totalVI_study_CPU", "sciPENN_sample_GPU", "sciPENN_sample_CPU", "sciPENN_study_GPU", "sciPENN_study_CPU",  "ADTnorm_study_manual_keepZero", "ADTnorm_sample_manual_keepZero") ## totalVI , sciPENN, scAR. decontPro


for(method in baseline_methods[c(1:9)]) {
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


## ======================================
## 2. plot adt expression density profile
## ======================================
## each panel is an ADT marker
## each track is a sample
## color by batch
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
            batch = rep(cell_x_feature$study_name, length(adt_marker_list)),
            disease = rep(cell_x_feature$sample_status, length(adt_marker_list))
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
        scale_x_continuous(trans = cytotidyr:::asinh_trans(cofactor = 5), breaks = c(0, 10, 100, 1000, 5000))
            #scale_x_continuous(trans = 'asinh', breaks = c(0, 10, 100, 1000, 5000)) +
            coord_cartesian(xlim = c(-1, 10000)) 
    }

    return(resPlot)
}
bw_list = list("Arcsinh" = 0.05, "CLR" = 0.05, "logCPM" = 0.01, "Arcsinh_CLR" = 0.003, "Harmony_RawCount_sample" = 0.1, "Harmony_RawCount_study" = 0.1, "Harmony_Arcsinh_sample" = 0.05, "Harmony_Arcsinh_study" = 0.05, "Harmony_CLR_sample" = 0.05, "Harmony_CLR_study" = 0.05,  "Harmony_logCPM_sample" = 0.05, "Harmony_logCPM_study" = 0.05,"CytofRUV_study" = 0.05, "CytofRUV_sample" = 0.05, "fastMNN_study" = 0.05,  "fastMNN_sample" = 0.05,"DSB" = 0.2, "ADTnorm_sample" = 0.1, "ADTnorm_study" = 0.2, "totalVI_sample_GPU" = 0.1, "totalVI_sample_CPU" = 0.1, "totalVI_study_GPU" = 0.1, "totalVI_study_CPU" = 0.1, "sciPENN_sample_GPU" = 0.05, "sciPENN_sample_CPU" = 0.05, "sciPENN_study_GPU" = 0.05, "sciPENN_study_CPU" = 0.05, "ADTnorm_study_manual_keepZero" = 0.2, "ADTnorm_sample_manual_keepZero" = 0.15,
"Harmony_ADTnorm_study" = 0.05, "Harmony_ADTnorm_sample" = 0.05)

trans_list = list("Arcsinh" = FALSE, "CLR" = FALSE, "logCPM" = TRUE, "Arcsinh_CLR" = FALSE, "Harmony_RawCount_sample" = TRUE, "Harmony_RawCount_study" = TRUE, "Harmony_Arcsinh_sample" = FALSE, "Harmony_Arcsinh_study" = FALSE, "Harmony_CLR_sample" = FALSE, "Harmony_CLR_study" = FALSE, "Harmony_logCPM_sample" = TRUE, "Harmony_logCPM_study" = TRUE, "CytofRUV_study" = FALSE, "CytofRUV_sample" = FALSE, "fastMNN_study" = FALSE, "fastMNN_sample" = FALSE, "DSB" = FALSE, "ADTnorm_sample" = FALSE, "ADTnorm_study" = FALSE, "totalVI_sample_GPU" = TRUE, "totalVI_sample_CPU" = TRUE, "totalVI_study_GPU" = TRUE, "totalVI_study_CPU" = TRUE, "sciPENN_sample_GPU" = FALSE, "sciPENN_sample_CPU" = FALSE, "sciPENN_study_GPU" = FALSE, "sciPENN_study_CPU" = FALSE, "ADTnorm_study_manual_keepZero" = FALSE, "ADTnorm_sample_manual_keepZero" = FALSE,
"Harmony_ADTnorm_study" = FALSE, "Harmony_ADTnorm_sample" = FALSE)

for(method in baseline_methods[1:17]){
  print(method)
  filename <- paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds")
  ans = readRDS(file = filename)

  if(grepl("study", method, ignore.case = TRUE)){
    
    fig_height = 16
    pdf(paste0(out_path, "/Figures/", method, "_", run_name, "_adt_density.pdf"), width = 60, height = fig_height)
    print(plot_adt_density(
    cell_x_adt = ans %>% data.frame, 
    cell_x_feature = adt_feature, 
    adt_marker_list = marker_list, 
    unit = "study",
    parameter_list = list(method_label = method, arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
    ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))))
    dev.off()

  }else if(grepl("sample", method, ignore.case = TRUE)){

    fig_height = 50
    pdf(paste0(out_path, "/Figures/", method, "_", run_name, "_adt_density.pdf"), width = 40, height = fig_height)
    print(plot_adt_density(
    cell_x_adt = ans %>% data.frame, 
    cell_x_feature = adt_feature, 
    adt_marker_list = marker_list, 
    unit = "sample",
    parameter_list = list(method_label = method, arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
    ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))))
    dev.off()

  }else{
    fig_height = 16
    pdf(paste0(out_path, "/Figures/", method, "_study_", run_name, "_adt_density.pdf"), width = 60, height = fig_height)
    print(plot_adt_density(
    cell_x_adt = ans %>% data.frame, 
    cell_x_feature = adt_feature, 
    adt_marker_list = marker_list, 
    unit = "study",
    parameter_list = list(method_label = method, arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
    ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))))
    dev.off()

    fig_height = 50
    pdf(paste0(out_path, "/Figures/", method, "_sample_", run_name, "_adt_density.pdf"), width = 40, height = fig_height)
    print(plot_adt_density(
    cell_x_adt = ans %>% data.frame, 
    cell_x_feature = adt_feature, 
    adt_marker_list = marker_list, 
    unit = "sample",
    parameter_list = list(method_label = method, arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
    ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))))
    dev.off()

  }
  
}

 
## totalVI density plot
dataset_name_list = run_name
unit_list = c("study", "sample")
study_name_levels = adt_feature$study_name %>% levels
sample_levels = adt_feature$sample %>% levels

for(dataset_name in dataset_name_list){
  print(dataset_name)
  for(unit in unit_list){
    print(unit)
    method = paste0("totalVI_", unit, "_GPU")

    if(grepl("study", unit, ignore.case = TRUE)){
      cell_x_adt = read.csv(paste0(out_path, "/CSV/cell_x_adt_totalVI_", unit, "_name_GPU_", dataset_name, ".csv"))
      cell_x_feature = read.csv(paste0(out_path, "/CSV/cell_x_feature_totalVI_", unit, "_name_GPU_", dataset_name, ".csv"))
      cell_x_feature$study_name = factor(cell_x_feature$study_name, levels = study_name_levels)
      cell_x_feature$sample = factor(cell_x_feature$sample, levels = sample_levels)

      saveRDS(cell_x_adt, paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds"))
      saveRDS(cell_x_feature, paste0(out_path, "/RDS/adt_feature_", method, "_", run_name, ".rds"))

      fig_height = 16
      pdf(paste0(out_path, "/Figures/totalVI_", unit, "_GPU_", dataset_name, "_adt_density.pdf"), width = 40, height = fig_height)
      print(plot_adt_density(
      cell_x_adt = cell_x_adt %>% data.frame, 
      cell_x_feature = cell_x_feature, 
      adt_marker_list = marker_list, 
      unit = "study",
      parameter_list = list(method_label = paste0("totalVI_", unit), arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
      ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))))
      dev.off()

    } 
    if(grepl("sample", unit, ignore.case = TRUE)){
      cell_x_adt = read.csv(paste0(out_path, "/CSV/cell_x_adt_totalVI_", unit, "_GPU_", dataset_name, ".csv"))
      cell_x_feature = read.csv(paste0(out_path, "/CSV/cell_x_feature_totalVI_", unit, "_GPU_", dataset_name, ".csv"))
      cell_x_feature$study_name = factor(cell_x_feature$study_name, levels = study_name_levels)
      cell_x_feature$sample = factor(cell_x_feature$sample, levels = sample_levels)
      
      saveRDS(cell_x_adt, paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds"))
      saveRDS(cell_x_feature, paste0(out_path, "/RDS/adt_feature_", method, "_", run_name, ".rds"))

      fig_height = 50
      pdf(paste0(out_path, "/Figures/totalVI_", unit, "_GPU_", dataset_name, "_adt_density.pdf"), width = 40, height = fig_height)
      print(plot_adt_density(
      cell_x_adt = cell_x_adt %>% data.frame, 
      cell_x_feature = cell_x_feature, 
      adt_marker_list = marker_list, 
      unit = "sample",
      parameter_list = list(method_label = paste0("totalVI_", unit), arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
      ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))))
      dev.off()
    }
  }
}



## sciPENN density plot
dataset_name_list = run_name
unit_list = c("study", "sample")
study_name_levels = adt_feature$study_name %>% levels
sample_levels = adt_feature$sample %>% levels

for(dataset_name in dataset_name_list){
  print(dataset_name)
  for(unit in unit_list){
    print(unit)
    method = paste0("sciPENN_", unit, "_GPU")

    if(grepl("study", unit, ignore.case = TRUE)){
      cell_x_adt = read.csv(paste0(out_path, "/CSV/sciPENN_Imputed_Protein_", dataset_name, "_", unit, "_name.csv"))
      cell_x_adt$rownames = NULL
      cell_x_feature = read.csv(paste0(out_path, "/CSV/sciPENN_Embedding_DataFeature_", dataset_name, "_", unit, "_name.csv"))
      cell_x_feature$rownames = NULL
      cell_x_feature$study_name = factor(cell_x_feature$study_name, levels = study_name_levels)
      cell_x_feature$sample = factor(cell_x_feature$sample, levels = sample_levels)

      saveRDS(cell_x_adt, paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds"))
      saveRDS(cell_x_feature, paste0(out_path, "/RDS/adt_feature_", method, "_", run_name, ".rds"))

      fig_height = 16
      pdf(paste0(out_path, "/Figures/sciPENN_",  unit, "_", dataset_name, "_adt_density.pdf"), width = 40, height = fig_height)
      print(plot_adt_density(
      cell_x_adt = cell_x_adt %>% data.frame, 
      cell_x_feature = cell_x_feature, 
      adt_marker_list = marker_list, 
      unit = "study",
      parameter_list = list(method_label = paste0("sciPENN_", unit), arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
      ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))))
      dev.off()

    } 
    if(grepl("sample", unit, ignore.case = TRUE)){
      
      cell_x_adt = read.csv(paste0(out_path, "/CSV/sciPENN_Imputed_Protein_", dataset_name, "_", unit, ".csv"))
      cell_x_adt$rownames = NULL
      cell_x_feature = read.csv(paste0(out_path, "/CSV/sciPENN_Embedding_DataFeature_", dataset_name, "_", unit, ".csv"))
      cell_x_feature$rownames = NULL
      cell_x_feature$study_name = factor(cell_x_feature$study_name, levels = study_name_levels)
      cell_x_feature$sample = factor(cell_x_feature$sample, levels = sample_levels)
      
      saveRDS(cell_x_adt, paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds"))
      saveRDS(cell_x_feature, paste0(out_path, "/RDS/adt_feature_", method, "_", run_name, ".rds"))

      fig_height = 50
      pdf(paste0(out_path, "/Figures/sciPENN_", unit, "_", dataset_name,  "_adt_density.pdf"), width = 40, height = fig_height)
      print(plot_adt_density(
      cell_x_adt = cell_x_adt %>% data.frame, 
      cell_x_feature = cell_x_feature, 
      adt_marker_list = marker_list, 
      unit = "sample",
      parameter_list = list(method_label = paste0("sciPENN_", unit), arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
      ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))))
      dev.off()

    }
  }
}

## ================================================================
## 6. Abnormal/Arbitrary manipulation of normalized protein counts
## ================================================================
select_methods <- c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount_study", "Harmony_RawCount_sample", "Harmony_Arcsinh_study", "Harmony_Arcsinh_sample", "Harmony_CLR_study", "Harmony_CLR_sample", "Harmony_logCPM_study", "Harmony_logCPM_sample", "fastMNN_study", "fastMNN_sample",  "CytofRUV_study", "CytofRUV_sample", "DSB","totalVI_sample_GPU", "totalVI_study_GPU", "sciPENN_sample_GPU", "sciPENN_study_GPU", "ADTnorm_sample", "ADTnorm_study", "ADTnorm_study_manual_keepZero", "ADTnorm_sample_manual_keepZero") 
mess_select_methods =  c("Arcsinh", "CLR", "Harmony_CLR_study", "fastMNN_study", "CytofRUV_study", "DSB", "totalVI_study_GPU", "sciPENN_study_GPU", "ADTnorm_study_manual_keepZero")

for(data_name in file_list){
  for(protein_marker in c("CD3", "CD4", "CD8", "CD19", "CD25", "CD56", "CD45RA")){
    protein_sym = sym(protein_marker)
    mess_plot = list()
      for(method in mess_select_methods){
        if(method %in% c("totalVI_sample_GPU", "totalVI_study_GPU", "sciPENN_sample_GPU", "sciPENN_study_GPU")){
          adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", method, "_", run_name, ".rds"))
        }else{
            adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", run_name, ".rds"))
        }
        ans = readRDS(file = paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds"))
        cell_index = which(adt_feature$study_name == data_name)
        mess_df =  data.frame(ans)[cell_index, ] %>% cbind(adt_feature[cell_index, ])
        mess_df$cell_type_l1 = factor(mess_df$cell_type_l1, levels = c("B", "CD4 T", "CD8 T", "monocytes", "NK", "DCs", "undefined") )
        
        mess_plot[[method]] = mess_df %>%
        ggplot(aes(x = cell_type_l1, y = !!protein_sym, color = cell_type_l1)) +
        geom_jitter(aes(x = cell_type_l1, y = !!protein_sym, color = cell_type_l1), data = mess_df %>% slice_sample(n = min(10000, nrow(mess_df))), size = 0.5, alpha = 0.5) +
        geom_violin(color = "black", alpha = 0.3) +
        geom_boxplot(color = "black", alpha = 0.3) +
        theme_bw(base_size = 25) +
        scale_color_manual(
          breaks = c("B", "CD4 T", "CD8 T", "monocytes", "NK", "DCs", "undefined"), 
          values = c(colorRampPalette(brewer.pal(8, "Set1"))(6), "grey")) +
        ggtitle(method) +
        rremove("legend") +
        xlab("")+
        rotate_x_text(angle = 20)
      }

      ggarrange(mess_plot[[mess_select_methods[1]]], mess_plot[[mess_select_methods[2]]], mess_plot[[mess_select_methods[3]]], mess_plot[[mess_select_methods[4]]], mess_plot[[mess_select_methods[5]]], mess_plot[[mess_select_methods[6]]], mess_plot[[mess_select_methods[7]]] + coord_cartesian(ylim = c(0, 2000)) + scale_y_continuous(trans = "log1p"), mess_plot[[mess_select_methods[8]]], mess_plot[[mess_select_methods[9]]], nrow = 3, ncol = 3)

      ggsave(file = paste0(out_path, "/Figures/", run_name, "_", data_name, "_messyPlot_", protein_marker, ".pdf"), width = 15, height = 20)
  }
}

library(grid)
raw_ans = readRDS(file = paste0(out_path, "/RDS/adt_data_RawCount_", run_name, ".rds"))
for(data_name in file_list){
  for(protein_marker in c("CD3", "CD4", "CD8", "CD19", "CD25", "CD56", "CD45RA")){
    protein_sym = sym(protein_marker)
    order_plot = list()
    for(method in mess_select_methods){
      if(method %in% c("totalVI_sample_GPU", "totalVI_study_GPU", "sciPENN_sample_GPU", "sciPENN_study_GPU")){
        adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", method, "_", run_name, ".rds"))
      }else{
          adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", run_name, ".rds"))
      }
      ans = readRDS(file = paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds"))
      cell_index = which(adt_feature$study_name == data_name)
      order_df =  data.frame(
        RawCount = data.frame(raw_ans)[cell_index, protein_marker] %>% rank,
        NormalizedCount = data.frame(ans)[cell_index, protein_marker] %>% rank
      ) %>% cbind(adt_feature[cell_index, ])
      order_df$cell_type_l1 = factor(order_df$cell_type_l1, levels = c("B", "CD4 T", "CD8 T", "monocytes", "NK", "DCs", "undefined") )

      order_plot[[method]] = order_df %>% ggplot(aes(x = RawCount, y = NormalizedCount, color = cell_type_l1)) +
      geom_point(alpha = 0.5) +
      theme_bw(base_size = 25) +
      scale_color_manual(
        breaks = c("B", "CD4 T", "CD8 T", "monocytes", "NK", "DCs", "undefined"), 
        values = c(colorRampPalette(brewer.pal(8, "Set1"))(6), "grey")) +
      ggtitle(method %>% gsub("_study.*", "", .)) +
      rremove("legend.title") +
      labs(x = NULL, y=NULL) + 
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank())

    }
    combined_plot = ggarrange(order_plot[[mess_select_methods[1]]], order_plot[[mess_select_methods[2]]], order_plot[[mess_select_methods[3]]], order_plot[[mess_select_methods[4]]], order_plot[[mess_select_methods[5]]], order_plot[[mess_select_methods[6]]], order_plot[[mess_select_methods[7]]], order_plot[[mess_select_methods[8]]], order_plot[[mess_select_methods[9]]], nrow = 3, ncol = 3, common.legend = TRUE, legend = "top")

    pdf(file = paste0(out_path, "/Figures/RankPlot_", run_name, "_", data_name, "_", protein_marker, ".pdf"), width = 15, height = 15)
    print(annotate_figure(combined_plot, 
                    bottom = text_grob("Cell Rank Based on Raw Counts", size = 30),
                    left = text_grob("Cell Rank Based on Normalized Counts", size = 30, rot = 90)))
    grid.text(protein_marker, x = 0.05, y = 0.98, just = c("left", "top"), gp = gpar(fontsize = 30))
    dev.off()
  }
}

## revision 20241203
fig_path = "./manuscript/results/revision/Figures/" #paste0(out_path, "/Figures/")

mess_select_methods =  c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount_study", "Harmony_RawCount_sample", "Harmony_Arcsinh_study", "Harmony_Arcsinh_sample", "Harmony_CLR_study", "Harmony_CLR_sample", "Harmony_logCPM_study", "Harmony_logCPM_sample", "fastMNN_study", "fastMNN_sample",  "CytofRUV_study", "CytofRUV_sample", "DSB","totalVI_sample_GPU", "totalVI_study_GPU", "sciPENN_sample_GPU", "sciPENN_study_GPU", "ADTnorm_sample", "ADTnorm_study", "ADTnorm_study_manual_keepZero")

data_name = "10X_malt_10k"
protein_marker = "CD19"
protein_sym = sym(protein_marker)
mess_plot = list()
for(method in mess_select_methods){
  if(method %in% c("totalVI_sample_GPU", "totalVI_study_GPU", "sciPENN_sample_GPU", "sciPENN_study_GPU")){
    adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", method, "_", run_name, ".rds"))
  }else{
      adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", run_name, ".rds"))
  }
  ans = readRDS(file = paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds"))
  cell_index = which(adt_feature$study_name == data_name)
  mess_df =  data.frame(ans)[cell_index, ] %>% cbind(adt_feature[cell_index, ])
  mess_df$cell_type_l1 = factor(mess_df$cell_type_l1, levels = c("B", "CD4 T", "CD8 T", "monocytes", "NK", "DCs", "undefined") )
  
  mess_plot[[method]] = mess_df %>%
  ggplot(aes(x = cell_type_l1, y = !!protein_sym, color = cell_type_l1)) +
  geom_jitter(aes(x = cell_type_l1, y = !!protein_sym, color = cell_type_l1), data = mess_df %>% slice_sample(n = min(10000, nrow(mess_df))), size = 0.5, alpha = 0.5) +
  geom_violin(color = "black", alpha = 0.3) +
  geom_boxplot(color = "black", alpha = 0.3) +
  theme_bw(base_size = 25) +
  scale_color_manual(
    breaks = c("B", "CD4 T", "CD8 T", "monocytes", "NK", "DCs", "undefined"), 
    values = c(colorRampPalette(brewer.pal(8, "Set1"))(6), "grey")) +
  ggtitle(method) +
  rremove("legend") +
  xlab("")+
  rotate_x_text(angle = 20)
}

pdf(paste0(fig_path, "messyPlot_scaling_severe.pdf"), width = 18, height = 11)
print(ggarrange(mess_plot[["Arcsinh"]], mess_plot[["CLR"]], mess_plot[["logCPM"]], mess_plot[["Arcsinh_CLR"]], nrow = 1, ncol = 4))
dev.off()

pdf(paste0(fig_path, "messyPlot_harmony_severe.pdf"), width = 18, height = 11)
print(ggarrange(mess_plot[["Harmony_RawCount_sample"]] + ylim(-1000, 1000), mess_plot[["Harmony_Arcsinh_sample"]], mess_plot[["Harmony_CLR_study"]], mess_plot[["Harmony_logCPM_sample"]], nrow = 1, ncol = 4))
dev.off()
