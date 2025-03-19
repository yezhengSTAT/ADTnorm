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
library(shapr)
library(fda)
library(pryr)
library(flowStats)
library(flowCore)
library(ncdfFlow)
library(flowViz)
library(pdfCluster)
library(cluster)
library(FlowSOM)
library(CytofRUV)
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

## clean data
clean_adt_name <- function(adt_name) {
    adt_rename <- adt_name %>%
        gsub("CD19-CAR", "CARCD19", .) %>%
        gsub("_PROT", "", .) %>%
        gsub("_TotalSeqB", "", .) %>%
        gsub("_control", "", .) %>%
        gsub("-GA.*", "", .) %>%
        gsub("-GC.*", "", .) %>%
        gsub("-GT.*", "", .) %>%
        gsub("-GG.*", "", .) %>%
        gsub("-CA.*", "", .) %>%
        gsub("-CC.*", "", .) %>%
        gsub("-CT.*", "", .) %>%
        gsub("-CG.*", "", .) %>%
        gsub("-AA.*", "", .) %>%
        gsub("-AC.*", "", .) %>%
        gsub("-AT.*", "", .) %>%
        gsub("-AG.*", "", .) %>%
        gsub("-TA.*", "", .) %>%
        gsub("-TC.*", "", .) %>%
        gsub("-TG.*", "", .) %>%
        gsub("-TT.*", "", .) %>%
        gsub("_.*", "", .) %>%
        gsub("ADT-", "", .) 

    if (!("CD8" %in% adt_rename)) {
        if ("CD8a" %in% adt_rename | "CD8A" %in% adt_rename) {
            ind <- which(adt_rename == "CD8a" | adt_rename == "CD8A")
            adt_rename[ind] <- "CD8"
        } else if ("CD8b" %in% adt_rename | "CD8B" %in% adt_rename) {
            ind <- which(adt_rename == "CD8B" | adt_rename == "CD8b")
            adt_rename[ind] <- "CD8"
        }
    }
    adt_rename <- replace(adt_rename, adt_rename == "MouseIgG1kappaisotype", "IgG1")
    adt_rename <- replace(adt_rename, adt_rename == "MouseIgG2akappaisotype", "IgG2a")
    adt_rename <- replace(adt_rename, adt_rename == "Mouse IgG2bkIsotype", "IgG2b")
    adt_rename <- replace(adt_rename, adt_rename == "RatIgG2bkIsotype", "IgG2b-Rat")

    adt_rename <- replace(adt_rename, adt_rename == "Mouse-IgG1", "IgG1")
    adt_rename <- replace(adt_rename, adt_rename == "Mouse-IgG2a", "IgG2a")
    adt_rename <- replace(adt_rename, adt_rename == "Mouse-IgG2b", "IgG2b")
    adt_rename <- replace(adt_rename, adt_rename == "Rat-IgG2b", "IgG2b-Rat")

    return(adt_rename)
}

clean_public_data = function(file_name, in_path){

    ## 10X_pbmc_10k
    ## one sample, manual gating cell types
    if(file_name == "10X_pbmc_10k"){
      data = readRDS(paste0(in_path, file_name, ".rds"))
      colData(data)$sample = paste0(file_name, "_sample1")
      colData(data)$batch = paste0(file_name, "_batch1")
      colData(data)$sample_status = "Healthy"
      colData(data)$study_name = file_name
    }

    ## "10X_pbmc_1k"
    ## one sample, manual gating cell types
    if(file_name == "10X_pbmc_1k"){
      data = readRDS(paste0(in_path, file_name, ".rds"))
      colData(data)$sample = paste0(file_name, "_sample1")
      colData(data)$batch = paste0(file_name, "_batch1")
      colData(data)$sample_status = "Healthy"
      colData(data)$study_name = file_name
    }

    ## 10X_pbmc_5k_v3
    ## one sample, manual gating cell types
    if(file_name == "10X_pbmc_5k_v3"){
      data = readRDS(paste0(in_path, file_name, ".rds"))
      colData(data)$sample = paste0(file_name, "_sample1")
      colData(data)$batch = paste0(file_name, "_batch1")
      colData(data)$sample_status = "Healthy"
      colData(data)$study_name = file_name
    }

    ## 10X_pbmc_5k_nextgem
    ## one sample, manual gating cell types
    if(file_name == "10X_pbmc_5k_nextgem"){
      data = readRDS(paste0(in_path, file_name, ".rds"))
      colData(data)$sample = paste0(file_name, "_sample1")
      colData(data)$batch = paste0(file_name, "_batch1")
      colData(data)$sample_status = "Healthy"
      colData(data)$study_name = file_name
    }

    ## 10X_malt_10k
    ## one sample, manual gating cell types
    if(file_name == "10X_malt_10k"){
      data = readRDS(paste0(in_path, file_name, ".rds"))
      colData(data)$sample = paste0(file_name, "_sample1")
      colData(data)$batch = paste0(file_name, "_batch1")
      colData(data)$sample_status = "MALTtumor"
      colData(data)$study_name = file_name
    }

    ## stuart_2019
    ## one sample, hashtag by 55 tags
    if(file_name == "stuart_2019"){
      data = readRDS(paste0(in_path, file_name, ".rds"))
      colData(data)$sample = paste0(file_name, "_sample1")
      colData(data)$batch = paste0(file_name, "_batch1")
      colData(data)$sample_status = "Healthy"
      colData(data)$study_name = file_name
      rownames(altExp(data)) = rownames(altExp(data)) %>% gsub("-.*", "", .)
    }

    ## granja_2019_pbmc
    ## 1 donor: 2 replicates
    if(file_name == "granja_2019_pbmc"){
      data = readRDS(paste0(in_path, file_name, ".rds"))
      colData(data)$sample = paste0(file_name, "_", colData(data)$donor)
      colData(data)$batch = paste0(file_name, "_batch1")
      colData(data)$sample_status = "Healthy"
      colData(data)$study_name = file_name
    }

    ## granja_2019_bmmc
    ## 1 donor: 2 replicates
    if(file_name == "granja_2019_bmmc"){
      data = readRDS(paste0(in_path, file_name, ".rds"))
      colData(data)$sample = paste0(file_name, "_", colData(data)$donor)
      colData(data)$batch = paste0(file_name, "_batch1")
      colData(data)$sample_status = "Healthy"
      colData(data)$study_name = file_name
    }

    ## hao_2020
    ## 8 volunteers x 3 time points
    ## 13 lanes
    if(file_name == "hao_2020"){
      data = readRDS(paste0(in_path, file_name, ".rds"))
      colData(data)$sample = paste0(file_name, "_", colData(data)$orig.ident) ## "P2_7", "P1_7", "P4_3", "P3_7", "P4_7", "P3_3", "P1_3", "P3_0", "P1_0", "P2_3", "P4_0", "P2_0", "P5_7", "P7_0", "P6_3", "P8_0", "P7_7", "P6_0", "P8_7", "P7_3", "P8_3", "P6_7", "P5_0", "P5_3"
      colData(data)$batch = paste0(file_name, "_", colData(data)$lane) ## "L1", "L2", "L3", "L4", "L5", "E2L1", "E2L2", "E2L3", "E2L4", "E2L5", "E2L6", "E2L7", "E2L8
      colData(data)$sample_status = "HIV Vaccine" 
      colData(data)$study_name = file_name
      rm_index = which(stringr::str_detect(rownames(altExp(data)), "-2"))
      altExp(data) = altExp(data)[-rm_index, ]
      rm_index = which(rownames(altExp(data)) == "CD8a")
      altExp(data) = altExp(data)[-rm_index, ]
      rownames(altExp(data)) = rownames(altExp(data)) %>% gsub("-1", "", .)
    }

    ## kotliarov_2020
    ## only have d0 data
    ## 20 patients that classified into low and high antibody responders
    if(file_name == "kotliarov_2020"){
      data = readRDS(paste0(in_path, file_name, ".rds"))
      colData(data)$sample = paste0(file_name, "_", colData(data)$sampleid)  ## "256" "273" "200" "233" "245" "212" "277" "207" "237" "261" "209" "205" "268" "234" "250" "215" "279" "229" "201" "236"
      colData(data)$batch = paste0(file_name, "_batch", colData(data)$batch) ## "1" "2"
      colData(data)$sample_status = "Lupus" 
      colData(data)$study_name = file_name
    }

    ## witkowski_2020
    ## 4 patients x 3 time points
    if(file_name == "witkowski_2020"){
      data = readRDS(paste0(in_path, file_name, ".rds"))
      colData(data)$sample = paste0(file_name, "_", colData(data)$sample_name %>% gsub("_CITE", "", .))
      colData(data)$batch = paste0(file_name, "_", colData(data)$library)
      colData(data)$sample_status = "B-ALL" 
      colData(data)$study_name = file_name
      rownames(altExp(data)) <- rownames(altExp(data)) %>% clean_adt_name() %>% gsub("HLA-DR", "HLADR", .) %>% gsub("HLA\\.DR", "HLADR", .) %>% gsub("PD-1", "PD1", .) %>% gsub(" ", "", .)  %>% gsub("CD8a", "CD8", .) %>% gsub("\\s*\\([^\\)]+\\)", "", .) %>% gsub("CD279", "PD1", .) %>% gsub("CD274", "PDL1", .)
      rownames(altExp(data)) = rownames(altExp(data)) %>% gsub("-.*", "", .)
    }

    ## only use "triana_2021.BM_All"
    ## refer to supptable 3
    ## 6 healthy donor: 3 young 3old; 3 AML patient
    if(file_name == "triana_2021"){
      data = readRDS(paste0(in_path, file_name, ".rds"))
      colData(data)$sample = paste0(file_name, "_", colData(data)$Batch)  ## "AML1"  "Q59"   "AML3"  "Aged2" "Aged1" "Aged3" "BM2"   "BM3"   "BM1" 
      colData(data)$batch = paste0(file_name, "_", colData(data)$Batch)
      colData(data)$sample_status = "Healthy"
      colData(data)$sample_status[which(colData(data)$Batch %in% c("AML1", "AML3", "Q59"))] = "AML" 
      colData(data)$study_name = file_name
      rownames(altExp(data)) = rownames(altExp(data)) %>% gsub("-AB", "", .)
    }

    ## buus 2021
    ## only keep the T cells filtered by CD3 positive
    if(file_name == "buus_2021_T"){
        data = readRDS(paste0(in_path, file_name, ".rds"))
        colData(data)$sample = paste0(file_name, "_", colData(data)$group)
        colData(data)$batch = paste0(file_name, "_", colData(data)$group)
        colData(data)$sample_status = "Healthy"
        colData(data)$sample_status[which(colData(data)$tissue == "Lung")] = "LungTumor" 
        colData(data)$study_name = file_name
    }

    ## Unify the surface marker name
    rownames(altExp(data)) <- rownames(altExp(data)) %>% clean_adt_name() %>% gsub("HLA-DR", "HLADR", .) %>% gsub("HLA\\.DR", "HLADR", .) %>% gsub("PD-1", "PD1", .) %>% gsub(" ", "", .)  %>% gsub("CD8a", "CD8", .) %>% gsub("\\s*\\([^\\)]+\\)", "", .) %>% gsub("CD279", "PD1", .) %>% gsub("CD274", "PDL1", .) %>% gsub("\\.", "", .) %>% gsub("-", "", .)

    rownames(altExp(data)) = rownames(altExp(data)) %>% toupper()

    colData(data)$ADTseqDepth <- colSums(counts(altExp(data)))
    return(data)
}

## read in data into list
## clean data regarding samples, batches, cell type, adt marker name
## data should have colData: sample, batch
adt_list <- c()
adt_data_dim <- 0
for (file in file_list) {
  print(file)
  data = clean_public_data(file, in_path)
  print(rownames(altExp(data)))
  adt_list <- rbind(adt_list, data.frame(adt = rownames(altExp(data)), file = file))
  adt_data_dim = adt_data_dim + ncol(data)
}

adt_full = adt_list$adt %>% unique %>% sort
adt_data_full = matrix(NA, ncol = adt_list$adt %>% unique %>% length, nrow = adt_data_dim)
colnames(adt_data_full) = adt_full
adt_feature = c()
cell_index = 1
for (file in file_list) {
  print(file)
  data = clean_public_data(file, in_path)
  adt_feature = rbind(adt_feature, colData(data) %>% data.frame %>% dplyr:: select(sample, batch, sample_status, study_name, ADTseqDepth, cell_type_l1, cell_type_l2))

  adt_marker_exist = adt_full[which(adt_full %in% rownames(altExp(data)))]
  adt_data_full[cell_index:(cell_index + ncol(data) - 1), adt_marker_exist] = t(counts(altExp(data)))[, adt_marker_exist] %>% as.matrix
  cell_index = cell_index + ncol(data)
   
}

## common ADT marker
marker_list = c("CD3", "CD4", "CD8",  "CD14", "CD19", "CD25", "CD45RA", "CD56", "CD127")

## cell_x_adt raw count matrix and cell_x_feature matrix
adt_data = adt_data_full[, marker_list] %>% data.frame
adt_feature = adt_feature %>% data.frame
adt_feature$study_name = factor(adt_feature$study_name, levels = file_list) ## setting levels to match the sample order
adt_feature$sample = factor(adt_feature$sample, levels = unique(adt_feature$sample)) ## setting levels to match the sample order
adt_feature_unique <- adt_feature %>% dplyr::select(sample, batch, sample_status, study_name) %>% unique
adt_feature_unique <- adt_feature_unique[match(adt_feature$sample %>% levels, adt_feature_unique$sample),]

saveRDS(adt_data_full, file = paste0(out_path, "/RDS/adt_data_full_RawCount_", run_name, ".rds"))
saveRDS(adt_data, file = paste0(out_path, "/RDS/adt_data_RawCount_", run_name, ".rds"))
saveRDS(adt_feature, file = paste0(out_path, "/RDS/adt_feature_", run_name, ".rds"))


## load data directly
marker_list = c("CD3", "CD4", "CD8",  "CD14", "CD19", "CD25", "CD45RA", "CD56", "CD127")
adt_data_full = readRDS(file = paste0(out_path, "/RDS/adt_data_full_RawCount_", run_name, ".rds"))
adt_data = readRDS(file = paste0(out_path, "/RDS/adt_data_RawCount_", run_name, ".rds"))
adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", run_name, ".rds"))
adt_feature_unique = adt_feature %>% dplyr::select(sample, batch, sample_status, study_name) %>% unique
adt_feature_unique = adt_feature_unique[match(adt_feature$sample %>% levels, adt_feature_unique$sample),]
baseline_mem = mem_used() ## 1.56GB~1.73GB for reading the data

## ==================================
## Normalization Methods Processing
## ==================================
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
      positive_peak = list(ADT = "CD3", sample = "buus_2021_T_PBMC_50ul_1_1000k"),
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
      positive_peak = list(ADT = "CD3", sample = "buus_2021_T"),
      lower_peak_thres = 0.05
    )
  ),
  # Harmony_RawCount_study = list(f = harmony_transform, start_adt_method = "RawCount", param = list(batch_column = "study_name")), 
  # Harmony_RawCount_sample = list(f = harmony_transform, start_adt_method = "RawCount", param = list(batch_column = "sample")),
  # Harmony_Arcsinh_study = list(f = harmony_transform, start_adt_method = "Arcsinh", param = list(batch_column = "study_name")),
  # Harmony_Arcsinh_sample = list(f = harmony_transform, start_adt_method = "Arcsinh", param = list(batch_column = "sample")),
  # Harmony_CLR_study = list(f = harmony_transform, start_adt_method = "CLR", param = list(batch_column = "study_name")),
  # Harmony_CLR_sample = list(f = harmony_transform, start_adt_method = "CLR", param = list(batch_column = "sample")),
  # Harmony_logCPM_study = list(f = harmony_transform, start_adt_method = "logCPM", param = list(batch_column = "study_name")),
  # Harmony_logCPM_sample = list(f = harmony_transform, start_adt_method = "logCPM", param = list(batch_column = "sample")),
  decontPro = list(f = decontPro_transform, start_adt_method = "RawCount", param = NULL)
)

## run normalization
baseline_methods <- c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "CytofRUV_study", "CytofRUV_sample", "fastMNN_study", "fastMNN_sample", "DSB", "Harmony_RawCount_study", "Harmony_RawCount_sample", "Harmony_Arcsinh_study", "Harmony_Arcsinh_sample", "Harmony_CLR_study", "Harmony_CLR_sample", "Harmony_logCPM_study", "Harmony_logCPM_sample", "ADTnorm_sample", "ADTnorm_study", "totalVI_sample_GPU", "totalVI_sample_CPU", "totalVI_study_GPU", "totalVI_study_CPU", "sciPENN_sample_GPU", "sciPENN_sample_CPU", "sciPENN_study_GPU", "sciPENN_study_CPU",  "ADTnorm_study_manual_keepZero", "ADTnorm_sample_manual_keepZero", "decontPro") ## scAR

## ====================================
## 1. run time and memory performance
## ====================================
run_time_memory <- c()
for(method in baseline_methods[c(30)]) {
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
  
  ## start of measuring the runtime and memory
  start_time <- Sys.time()
  mem_before <- mem_used()
  ans <- do.call(
    norm_f,
    list(cell_x_adt = cell_x_adt, cell_x_feature = adt_feature, parameter_list = method_param)
  )
  if(stringr::str_detect(method, "CytofRUV")){
    colnames(ans) = cell_x_adt_colnames
  }
  mem_after <- mem_used()
  mem_change <- mem_after - mem_before
  end_time <- Sys.time()
  run_time <- as.numeric(end_time - start_time, units = "secs")
  
  run_time_memory <- data.frame(run_time_seconds = run_time %>% as.numeric, memory_after_MB = (mem_after %>% as.numeric)/1000000, memory_change_MB = (mem_change %>% as.numeric)/1000000, method = method, run_name = run_name) %>% rbind(run_time_memory, .)
   
  saveRDS(ans, file = filename)

}
saveRDS(run_time_memory, file = paste0(out_path, "/RDS/decontPro_run_time_memory_", run_name, ".rds"))

# saveRDS(run_time_memory, file = paste0(out_path, "/RDS/baselineMethods_run_time_memory_", run_name, ".rds"))
# saveRDS(run_time_memory, file = paste0(out_path, "/RDS/ADTnorm_run_time_memory_", run_name, ".rds"))


# run_time_memory_harmony = rbind(
#   readRDS(file = paste0(out_path, "/RDS/Harmony2_run_time_memory_", run_name, ".rds")),
#   readRDS(file = paste0(out_path, "/RDS/Harmony_RawCount_run_time_memory_", run_name, ".rds"))
# )
# run_time_memory = rbind(run_time_memory, run_time_memory_harmony)
# run_time_memory$method = factor(run_time_memory$method, levels = baseline_methods)
# saveRDS(run_time_memory, file = paste0(out_path, "/RDS/baselineMethods_run_time_memory_", run_name, ".rds"))

# run_time_memory_adtnorm = run_time_memory
# run_time_memory = rbind(run_time_memory, run_time_memory_adtnorm)
# saveRDS(run_time_memory, file = paste0(out_path, "/RDS/AllMethods_run_time_memory_", run_name, ".rds"))

run_time_memory = readRDS(file = paste0(out_path, "/RDS/AllMethods_run_time_memory_", run_name, ".rds"))
run_time_memory$memory_after_MB[10:17] = run_time_memory$memory_after_MB[10:17] + as.numeric(baseline_mem)/1000000
run_time_memory$method = factor(run_time_memory$method, levels = baseline_methods)
run_time_memory$memory_after_GB = run_time_memory$memory_after_MB/1024 %>% as.numeric

## add runtime and memory for totalVI
study_mem = fread("./manuscript/results/public13Dataset_CITEseq/CSV/public13Dataset_CITEseq_totalVI_memory_usage.log")
study_mem$V3[2]
sample_mem = fread("./manuscript/results/public13Dataset_CITEseq/CSV/public13Dataset_CITEseq_totalVI_memory_usage_sample.log")
sample_mem$V3[2]
study_runtime = fread("./manuscript/results/public13Dataset_CITEseq/CSV/public13Dataset_CITEseq_totalVI_runtime.csv", header = FALSE, sep = ":")
study_runtime_sum = study_runtime[, 3] %>% gsub("}.*", "", .) %>% as.numeric + study_runtime[, 5] %>% gsub("}.*", "", .) %>% as.numeric
sample_runtime = fread("./manuscript/results/public13Dataset_CITEseq/CSV/public13Dataset_CITEseq_totalVI_runtime_sample.csv", header = FALSE, sep = ":")
sample_runtime_sum = sample_runtime[, 3] %>% gsub("}.*", "", .) %>% as.numeric + sample_runtime[, 5] %>% gsub("}.*", "", .) %>% as.numeric

run_time_memory = rbind(run_time_memory, data.frame(run_time_seconds = study_runtime_sum, memory_after_GB = study_mem$V3[2], memory_after_MB = study_mem$V3[2] * 1024, memory_change_MB = 0, method = "totalVI_study_CPU", run_name = run_name))
run_time_memory = rbind(run_time_memory, data.frame(run_time_seconds = sample_runtime_sum, memory_after_GB = sample_mem$V3[2], memory_after_MB = sample_mem$V3[2] * 1024, memory_change_MB = 0, method = "totalVI_sample_CPU", run_name = run_name))

study_mem = fread("./manuscript/results/public13Dataset_CITEseq/CSV/public13Dataset_CITEseq_totalVI_memory_usage_study_name_GPU.log")
study_mem$V3[2]
sample_mem = fread("./manuscript/results/public13Dataset_CITEseq/CSV/public13Dataset_CITEseq_totalVI_memory_usage_sample_GPU.log")
sample_mem$V3[2]
study_runtime = fread("./manuscript/results/public13Dataset_CITEseq/CSV/public13Dataset_CITEseq_totalVI_runtime_study_name_GPU.csv", header = FALSE, sep = ":")
study_runtime_sum = study_runtime[, 3] %>% gsub("}.*", "", .) %>% as.numeric + study_runtime[, 5] %>% gsub("}.*", "", .) %>% as.numeric
sample_runtime = fread("./manuscript/results/public13Dataset_CITEseq/CSV/public13Dataset_CITEseq_totalVI_runtime_sample_GPU.csv", header = FALSE, sep = ":")
sample_runtime_sum = sample_runtime[, 3] %>% gsub("}.*", "", .) %>% as.numeric + sample_runtime[, 5] %>% gsub("}.*", "", .) %>% as.numeric

run_time_memory = rbind(run_time_memory, data.frame(run_time_seconds = study_runtime_sum, memory_after_GB = study_mem$V3[2], memory_after_MB = study_mem$V3[2] * 1024, memory_change_MB = 0, method = "totalVI_study_GPU", run_name = run_name))
run_time_memory = rbind(run_time_memory, data.frame(run_time_seconds = sample_runtime_sum, memory_after_GB = sample_mem$V3[2], memory_after_MB = sample_mem$V3[2] * 1024, memory_change_MB = 0, method = "totalVI_sample_GPU", run_name = run_name))

## ad runtime and memory for sciPENN
study_mem = fread("./manuscript/results/public13Dataset_CITEseq/CSV/sciPENN_memory_usage_public13Dataset_CITEseq_study_name_GPU.csv")
study_mem$V3[2]
sample_mem = fread("./manuscript/results/public13Dataset_CITEseq/CSV/sciPENN_memory_usage_public13Dataset_CITEseq_sample_GPU.csv")
sample_mem$V3[2]

study_runtime = fread("./manuscript/results/public13Dataset_CITEseq/CSV/sciPENN_runtime_public13Dataset_CITEseq_study_name_GPU.csv", header = FALSE, sep = ":")
study_runtime_sum = study_runtime %>% as.numeric
sample_runtime = fread("./manuscript/results/public13Dataset_CITEseq/CSV/sciPENN_runtime_public13Dataset_CITEseq_sample_GPU.csv", header = FALSE, sep = ":")
sample_runtime_sum = sample_runtime %>% as.numeric

run_time_memory = rbind(run_time_memory, data.frame(run_time_seconds = study_runtime_sum, memory_after_GB = study_mem$V3[2], memory_after_MB = study_mem$V3[2] * 1024, memory_change_MB = 0, method = "sciPENN_study_GPU", run_name = run_name))
run_time_memory = rbind(run_time_memory, data.frame(run_time_seconds = sample_runtime_sum, memory_after_GB = sample_mem$V3[2], memory_after_MB = sample_mem$V3[2] * 1024, memory_change_MB = 0, method = "sciPENN_sample_GPU", run_name = run_name))

study_mem = fread("./manuscript/results/public13Dataset_CITEseq/CSV/sciPENN_memory_usage_public13Dataset_CITEseq_study_name_CPU.csv")
study_mem$V3[2]
sample_mem = fread("./manuscript/results/public13Dataset_CITEseq/CSV/sciPENN_memory_usage_public13Dataset_CITEseq_sample_CPU.csv")
sample_mem$V3[2]

study_runtime = fread("./manuscript/results/public13Dataset_CITEseq/CSV/sciPENN_runtime_public13Dataset_CITEseq_study_name_CPU.csv", header = FALSE, sep = ":")
study_runtime_sum = study_runtime %>% as.numeric
sample_runtime = fread("./manuscript/results/public13Dataset_CITEseq/CSV/sciPENN_runtime_public13Dataset_CITEseq_sample_CPU.csv", header = FALSE, sep = ":")
sample_runtime_sum = sample_runtime %>% as.numeric

run_time_memory = rbind(run_time_memory, data.frame(run_time_seconds = study_runtime_sum, memory_after_GB = study_mem$V3[2], memory_after_MB = study_mem$V3[2] * 1024, memory_change_MB = 0, method = "sciPENN_study_CPU", run_name = run_name))
run_time_memory = rbind(run_time_memory, data.frame(run_time_seconds = sample_runtime_sum, memory_after_GB = sample_mem$V3[2], memory_after_MB = sample_mem$V3[2] * 1024, memory_change_MB = 0, method = "sciPENN_sample_CPU", run_name = run_name))

run_time_memory

saveRDS(run_time_memory, file = paste0(out_path, "/RDS/AllMethods_totalVI_sciPENN_run_time_memory_", run_name, ".rds"))

## add decontPro memory and runtime
run_time_memory = readRDS(file = paste0(out_path, "/RDS/AllMethods_totalVI_sciPENN_run_time_memory_", run_name, ".rds"))

run_time_memory_decontPro = readRDS(file = paste0(out_path, "/RDS/decontPro_run_time_memory_", run_name, ".rds"))

run_time_memory_decontPro$memory_after_GB = run_time_memory_decontPro$memory_after_MB/1024 %>% as.numeric
run_time_memory = rbind(run_time_memory, run_time_memory_decontPro)
saveRDS(run_time_memory, file = paste0(out_path, "/RDS/AllMethods_totalVI_sciPENN_decontPro_run_time_memory_", run_name, ".rds"))

filtered_data1 <- run_time_memory %>% dplyr::filter(!grepl("study", method))
filtered_data1$method = filtered_data1$method %>% gsub("_sample", "", .)
filtered_data1$method = factor(filtered_data1$method, levels = c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount", "Harmony_Arcsinh", "Harmony_CLR", "Harmony_logCPM",  "fastMNN", "CytofRUV", 
"DSB", "decontPro", "totalVI_GPU", "totalVI_CPU", "sciPENN_GPU", "sciPENN_CPU", "ADTnorm"))

fillColor <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(unique(filtered_data1$method))) %>% rev

filtered_data1 %>% ggplot(aes(x = method, y = run_time_seconds, label = round(run_time_seconds, 1), fill = method)) +
geom_bar(stat = "identity") +
geom_text(vjust = -0.1, size = 7) +
theme_bw(base_size = 30) +
scale_fill_manual(values = fillColor) +
xlab("") +
ylab("Runtime (seconds)") +
rremove("legend") +
rotate_x_text(angle = 90) +
scale_y_break(c(8000, 36000), scales = "free") +
ggtitle("Unit: Sample")

ggsave(paste0(fig_path, "run_time_sample.pdf"), width = 13, height = 15)

filtered_data1 %>% ggplot(aes(x = method, y = round(memory_after_GB, 1), label = round(memory_after_GB, 1), fill = method)) +
geom_bar(stat = "identity") +
geom_text(vjust = -0.1, size = 7) +
theme_bw(base_size = 30) +
scale_fill_manual(values = fillColor) +
xlab("") +
ylab("Memory (GB)") +
rremove("legend") +
rotate_x_text(angle = 90) +
scale_y_continuous(trans = "log1p", breaks = c(0, round(as.numeric(baseline_mem)/(1000000*1024), 2), 30, 200, 400, 600)) + 
ggtitle("Unit: Sample") +
geom_hline(yintercept = 1.6, linetype = "dashed", color = "grey20", size = 1)
# geom_hline(yintercept = round(as.numeric(baseline_mem)/(1000000*1024), 1), linetype = "dashed", color = "grey20", size = 1)
ggsave(paste0(fig_path, "run_memory_sample.pdf"), width = 13, height = 13)

filtered_data2 <- run_time_memory %>% dplyr::filter(!grepl("sample", method))
filtered_data2$method = filtered_data2$method %>% gsub("_study", "", .)
filtered_data2$method = factor(filtered_data2$method, levels = c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount", "Harmony_Arcsinh", "Harmony_CLR", "Harmony_logCPM",  "fastMNN", "CytofRUV", 
"DSB", "decontPro", "totalVI_GPU", "totalVI_CPU", "sciPENN_GPU", "sciPENN_CPU", "ADTnorm"))


filtered_data2 %>% ggplot(aes(x = method, y = run_time_seconds, label = round(run_time_seconds, 1), fill = method)) +
geom_bar(stat = "identity") +
geom_text(vjust = -0.1, size = 7) +
theme_bw(base_size = 30) +
scale_fill_manual(values = fillColor) +
xlab("") +
ylab("Runtime (seconds)") +
rremove("legend") +
rotate_x_text(angle = 90)+
scale_y_break(c(4500, 40000), scales = "free") +
ggtitle("Unit: Study")

ggsave(paste0(fig_path, "run_time_study.pdf"), width = 13, height = 15)

filtered_data2 %>% ggplot(aes(x = method, y = round(memory_after_GB, 1), label = round(memory_after_GB, 1), fill = method)) +
geom_bar(stat = "identity") +
geom_text(vjust = -0.1, size = 7) +
theme_bw(base_size = 30) +
scale_fill_manual(values = fillColor) +
xlab("") +
ylab("Memory (GB)") +
rremove("legend") +
rotate_x_text(angle = 90)+
ggtitle("Unit: Study") +
scale_y_continuous(trans = "log1p", breaks = c(0, round(as.numeric(baseline_mem)/(1000000*1024), 2), 30, 200, 400, 600)) + 
geom_hline(yintercept = 1.6, linetype = "dashed", color = "grey20", size = 1)
# geom_hline(yintercept = as.numeric(baseline_mem)/(1000000*1024), linetype = "dashed", color = "grey20", size = 1)
ggsave(paste0(fig_path, "run_memory_study.pdf"), width = 13, height = 13)



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
"Harmony_ADTnorm_study" = 0.05, "Harmony_ADTnorm_sample" = 0.05, "decontPro" = 0.2)

trans_list = list("Arcsinh" = FALSE, "CLR" = FALSE, "logCPM" = TRUE, "Arcsinh_CLR" = FALSE, "Harmony_RawCount_sample" = TRUE, "Harmony_RawCount_study" = TRUE, "Harmony_Arcsinh_sample" = FALSE, "Harmony_Arcsinh_study" = FALSE, "Harmony_CLR_sample" = FALSE, "Harmony_CLR_study" = FALSE, "Harmony_logCPM_sample" = TRUE, "Harmony_logCPM_study" = TRUE, "CytofRUV_study" = FALSE, "CytofRUV_sample" = FALSE, "fastMNN_study" = FALSE, "fastMNN_sample" = FALSE, "DSB" = FALSE, "ADTnorm_sample" = FALSE, "ADTnorm_study" = FALSE, "totalVI_sample_GPU" = TRUE, "totalVI_sample_CPU" = TRUE, "totalVI_study_GPU" = TRUE, "totalVI_study_CPU" = TRUE, "sciPENN_sample_GPU" = FALSE, "sciPENN_sample_CPU" = FALSE, "sciPENN_study_GPU" = FALSE, "sciPENN_study_CPU" = FALSE, "ADTnorm_study_manual_keepZero" = FALSE, "ADTnorm_sample_manual_keepZero" = FALSE, "Harmony_ADTnorm_study" = FALSE, "Harmony_ADTnorm_sample" = FALSE, "decontPro" = FALSE)

for(method in baseline_methods[30]){
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

## decontPro density Plot
method = "decontPro"  
filename <- paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds")
ans = readRDS(file = filename)
fig_height = 16
pdf(paste0(out_path, "/Figures/", method, "_study_", run_name, "_adt_density.pdf"), width = 60, height = fig_height)
print(plot_adt_density(
cell_x_adt = ans %>% data.frame, 
cell_x_feature = adt_feature, 
adt_marker_list = marker_list, 
unit = "study",
parameter_list = list(method_label = method, arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))) +
scale_x_continuous(trans = "log1p", breaks = c(1, 5, 10^seq(1, 4)))
)
dev.off()

fig_height = 50
pdf(paste0(out_path, "/Figures/", method, "_sample_", run_name, "_adt_density.pdf"), width = 40, height = fig_height)
print(plot_adt_density(
cell_x_adt = ans %>% data.frame, 
cell_x_feature = adt_feature, 
adt_marker_list = marker_list, 
unit = "sample",
parameter_list = list(method_label = method, arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list)))+
scale_x_continuous(trans = "log1p", breaks = c(1, 5, 10^seq(1, 4))))
dev.off()

## totalVI density plot
dataset_name_list = c("public12Dataset_CITEseq_notriana", "public13Dataset_CITEseq")
unit_list = c("study", "sample")
study_name_levels = adt_feature$study_name %>% levels
sample_levels = adt_feature$sample %>% levels

for(dataset_name in dataset_name_list[1]){
  print(dataset_name)
  for(unit in unit_list){
    print(unit)
    method = paste0("totalVI_", unit, "_GPU")

    if(grepl("study", unit, ignore.case = TRUE)){
      cell_x_adt = read.csv(paste0(out_path, "/CSV/cell_x_adt_totalVI_", dataset_name, "_", unit, "_name.csv"))
      cell_x_feature = read.csv(paste0(out_path, "/CSV/cell_x_feature_totalVI_", dataset_name, "_", unit, "_name.csv"))
      cell_x_feature$study_name = factor(cell_x_feature$study_name, levels = study_name_levels)
      cell_x_feature$sample = factor(cell_x_feature$sample, levels = sample_levels)

      saveRDS(cell_x_adt, paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds"))
      saveRDS(cell_x_feature, paste0(out_path, "/RDS/adt_feature_", method, "_", run_name, ".rds"))

      # fig_height = 16
      # pdf(paste0(out_path, "/Figures/totalVI_", dataset_name, "_", unit, "_adt_density.pdf"), width = 40, height = fig_height)
      # print(plot_adt_density(
      # cell_x_adt = cell_x_adt %>% data.frame, 
      # cell_x_feature = cell_x_feature, 
      # adt_marker_list = marker_list, 
      # unit = "study",
      # parameter_list = list(method_label = paste0("totalVI_", unit), arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
      # ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))))
      # dev.off()

    } 
    if(grepl("sample", unit, ignore.case = TRUE)){
      cell_x_adt = read.csv(paste0(out_path, "/CSV/cell_x_adt_totalVI_", dataset_name, "_", unit, ".csv"))
      cell_x_feature = read.csv(paste0(out_path, "/CSV/cell_x_feature_totalVI_", dataset_name, "_", unit, ".csv"))
      cell_x_feature$study_name = factor(cell_x_feature$study_name, levels = study_name_levels)
      cell_x_feature$sample = factor(cell_x_feature$sample, levels = sample_levels)
      
      saveRDS(cell_x_adt, paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds"))
      saveRDS(cell_x_feature, paste0(out_path, "/RDS/adt_feature_", method, "_", run_name, ".rds"))

      # fig_height = 50
      # pdf(paste0(out_path, "/Figures/totalVI_", dataset_name, "_", unit, "_adt_density.pdf"), width = 40, height = fig_height)
      # print(plot_adt_density(
      # cell_x_adt = cell_x_adt %>% data.frame, 
      # cell_x_feature = cell_x_feature, 
      # adt_marker_list = marker_list, 
      # unit = "sample",
      # parameter_list = list(method_label = paste0("totalVI_", unit), arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
      # ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))))
      # dev.off()
    }
  }
}


for(dataset_name in dataset_name_list[2]){
  print(dataset_name)
  for(unit in unit_list){
    print(unit)
    method = paste0("totalVI_", unit, "_GPU_test_withoutcelltypecolumns") #"_include_protein_background"
    if(grepl("study", unit, ignore.case = TRUE)){
      cell_x_adt = read.csv(paste0(out_path, "/CSV/cell_x_adt_totalVI_", dataset_name, "_", unit, "_name", suffix, ".csv"))
      cell_x_feature = read.csv(paste0(out_path, "/CSV/cell_x_feature_totalVI_", dataset_name, "_", unit, "_name", suffix, ".csv"))
      cell_x_feature$study_name = factor(cell_x_feature$study_name, levels = study_name_levels)
      cell_x_feature$sample = factor(cell_x_feature$sample, levels = sample_levels)
      
      saveRDS(cell_x_adt, paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds"))
      saveRDS(cell_x_feature, paste0(out_path, "/RDS/adt_feature_", method, "_", run_name, ".rds"))

      # fig_height = 16
      # pdf(paste0(out_path, "/Figures/totalVI_", dataset_name, "_", unit, "_include_protein_background_adt_density.pdf"), width = 40, height = fig_height)
      # print(plot_adt_density(
      # cell_x_adt = cell_x_adt %>% data.frame, 
      # cell_x_feature = cell_x_feature, 
      # adt_marker_list = marker_list, 
      # unit = "study",
      # parameter_list = list(method_label = paste0("totalVI_", unit), arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
      # ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))))
      # dev.off()

    } 
    if(grepl("sample", unit, ignore.case = TRUE)){
      cell_x_adt = read.csv(paste0(out_path, "/CSV/cell_x_adt_totalVI_", dataset_name, "_", unit, suffix, ".csv"))
      cell_x_feature = read.csv(paste0(out_path, "/CSV/cell_x_feature_totalVI_", dataset_name, "_", unit, suffix, ".csv"))
      cell_x_feature$study_name = factor(cell_x_feature$study_name, levels = study_name_levels)
      cell_x_feature$sample = factor(cell_x_feature$sample, levels = sample_levels)
      
      saveRDS(cell_x_adt, paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds"))
      saveRDS(cell_x_feature, paste0(out_path, "/RDS/adt_feature_", method, "_", run_name, ".rds"))

      # fig_height = 50
      # pdf(paste0(out_path, "/Figures/totalVI_", dataset_name, "_", unit, "_include_protein_background_adt_density.pdf"), width = 40, height = fig_height)
      # print(plot_adt_density(
      # cell_x_adt = cell_x_adt %>% data.frame, 
      # cell_x_feature = cell_x_feature, 
      # adt_marker_list = marker_list, 
      # unit = "sample",
      # parameter_list = list(method_label = paste0("totalVI_", unit), arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
      # ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))))
      # dev.off()

    }
  }
}

## sciPENN density plot
dataset_name_list = c("public12Dataset_CITEseq_notriana", "public13Dataset_CITEseq")
unit_list = c("study", "sample")
study_name_levels = adt_feature$study_name %>% levels
sample_levels = adt_feature$sample %>% levels

for(dataset_name in dataset_name_list[1]){
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

      # fig_height = 16
      # pdf(paste0(out_path, "/Figures/sciPENN_", dataset_name, "_", unit, "_adt_density.pdf"), width = 40, height = fig_height)
      # print(plot_adt_density(
      # cell_x_adt = cell_x_adt %>% data.frame, 
      # cell_x_feature = cell_x_feature, 
      # adt_marker_list = marker_list, 
      # unit = "study",
      # parameter_list = list(method_label = paste0("sciPENN_", unit), arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
      # ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))))
      # dev.off()

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

      # fig_height = 50
      # pdf(paste0(out_path, "/Figures/sciPENN_", dataset_name, "_", unit, "_adt_density.pdf"), width = 40, height = fig_height)
      # print(plot_adt_density(
      # cell_x_adt = cell_x_adt %>% data.frame, 
      # cell_x_feature = cell_x_feature, 
      # adt_marker_list = marker_list, 
      # unit = "sample",
      # parameter_list = list(method_label = paste0("sciPENN_", unit), arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
      # ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))))
      # dev.off()

    }
  }
}

## =============================================================================
## 3. umap 
## =============================================================================
## transfer csv file into rds file using previous density plot generation codes
## generate the umap by publicData_CITEseq_umap.R


## =========================
## 4. Summary of Evaluation
## =========================
## run evaluations on cell type and batches 

select_methods <- c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount_study", "Harmony_RawCount_sample", "Harmony_Arcsinh_study", "Harmony_Arcsinh_sample", "Harmony_CLR_study", "Harmony_CLR_sample", "Harmony_logCPM_study", "Harmony_logCPM_sample", "fastMNN_study", "fastMNN_sample",  "CytofRUV_study", "CytofRUV_sample", "DSB", "decontPro", "totalVI_sample_GPU", "totalVI_study_GPU", "sciPENN_sample_GPU", "sciPENN_study_GPU", "ADTnorm_sample", "ADTnorm_study", "ADTnorm_study_manual_keepZero", "ADTnorm_sample_manual_keepZero") #,  "ADTnorm_study_manual_keepZero", "ADTnorm_sample_manual_keepZero"

ari_summary_reso0.2 = c()
for(method in select_methods){
    
  tmp = readRDS(paste0(out_path, "/RDS/adt_", method, "_", run_name, "_ari_resolution0.2.rds"))
  ari_summary_reso0.2 = data.frame(
    ari_sample = tmp[[1]],
    ari_site = tmp[[2]],
    ari_refineCT = tmp[[3]],
    ari_broadCT = tmp[[4]],
    method = method
  ) %>% rbind(ari_summary_reso0.2, .)
  
}
ari_summary_reso0.15 = c()
for(method in select_methods){
    
  tmp = readRDS(paste0(out_path, "/RDS/adt_", method, "_", run_name, "_ari_resolution0.15.rds"))
  ari_summary_reso0.15 = data.frame(
    ari_sample = tmp[[1]],
    ari_site = tmp[[2]],
    ari_refineCT = tmp[[3]],
    ari_broadCT = tmp[[4]],
    method = method
  ) %>% rbind(ari_summary_reso0.15, .)
  
}

ari_all_reso0.2 = ari_summary_reso0.2 %>% data.frame
ari_all_reso0.2$method = factor(ari_all_reso0.2$method, levels = select_methods)
# Calculate mean and standard deviation for each group
summary_df_reso0.2 <- ari_all_reso0.2 %>%
  group_by(method) %>%
  summarise(
    mean_ari_sample = mean(ari_sample),
    sd_ari_sample = sd(ari_sample),
    mean_ari_refineCT = mean(ari_refineCT),
    sd_ari_refineCT = sd(ari_refineCT),
    mean_ari_site = mean(ari_site),
    sd_ari_site = sd(ari_site),
    mean_ari_broadCT = mean(ari_broadCT),
    sd_ari_broadCT = sd(ari_broadCT)
  )

ari_all_reso0.15 = ari_summary_reso0.15 %>% data.frame
ari_all_reso0.15$method = factor(ari_all_reso0.15$method, levels = select_methods)
# Calculate mean and standard deviation for each group
summary_df_reso0.15 <- ari_all_reso0.15 %>%
  group_by(method) %>%
  summarise(
    mean_ari_sample = mean(ari_sample),
    sd_ari_sample = sd(ari_sample),
    mean_ari_refineCT = mean(ari_refineCT),
    sd_ari_refineCT = sd(ari_refineCT),
    mean_ari_site = mean(ari_site),
    sd_ari_site = sd(ari_site),
    mean_ari_broadCT = mean(ari_broadCT),
    sd_ari_broadCT = sd(ari_broadCT)
  )


pdf(paste0(out_path, "/Figures/Evaluation_metrics_ARI_13clusters_withDecontPro.pdf"), width = 15, height = 11)

rename_methods = c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony(RawCount)", "Harmony(Arcsinh)", "Harmony(CLR)", "Harmony(logCPM)", "fastMNN", "CytofRUV", "DSB", "decontPro", "totalVI", "sciPENN", "ADTnorm(default)", "ADTnorm(customized)")
summary_df_study = summary_df_reso0.15 %>% dplyr::filter(method %in% c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount_study", "Harmony_Arcsinh_study", "Harmony_CLR_study", "Harmony_logCPM_study", "fastMNN_study", "CytofRUV_study", "DSB",  "decontPro", "totalVI_study_GPU", "sciPENN_study_GPU", "ADTnorm_study", "ADTnorm_study_manual_keepZero")) %>% data.frame
summary_df_study$method = factor(summary_df_study$method, levels = c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount_study", "Harmony_Arcsinh_study", "Harmony_CLR_study", "Harmony_logCPM_study", "fastMNN_study", "CytofRUV_study", "DSB",  "decontPro", "totalVI_study_GPU", "sciPENN_study_GPU", "ADTnorm_study", "ADTnorm_study_manual_keepZero"), labels = rename_methods)

summary_df_sample = summary_df_reso0.2 %>% dplyr::filter(method %in% c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount_sample", "Harmony_Arcsinh_sample", "Harmony_CLR_sample", "Harmony_logCPM_sample", "fastMNN_sample", "CytofRUV_sample", "DSB",  "decontPro", "totalVI_sample_GPU", "sciPENN_sample_GPU", "ADTnorm_sample_manual_keepZero", "ADTnorm_sample")) %>% data.frame
summary_df_sample$method = factor(summary_df_sample$method, levels = c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount_sample", "Harmony_Arcsinh_sample", "Harmony_CLR_sample", "Harmony_logCPM_sample", "fastMNN_sample", "CytofRUV_sample", "DSB",  "decontPro", "totalVI_sample_GPU", "sciPENN_sample_GPU", "ADTnorm_sample_manual_keepZero", "ADTnorm_sample"), labels = rename_methods)

fillColor <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(rename_methods)) %>% rev

p1 = summary_df_study %>% ggplot(aes(x = mean_ari_site, y = mean_ari_broadCT, label = method, color = method)) +
  geom_point(aes(color = method), size = 4) + # Dots for the mean values
  geom_errorbar(aes(ymin = mean_ari_broadCT - sd_ari_broadCT, ymax = mean_ari_broadCT + sd_ari_broadCT), 
                width = 0.001, color = 'grey') + # Vertical error bars
  geom_errorbarh(aes(xmin = mean_ari_site - sd_ari_site, xmax = mean_ari_site + sd_ari_site), 
                 height = 0.01, color = 'grey') + # Horizontal error bars
  xlab("ARI of Batch Effect (Study)") +
  ylab("ARI of Cell Type Separation (Broad)") + 
  geom_text_repel(size = 5) +
  theme_bw(base_size = 30) +
  scale_color_manual(values = fillColor) +
  rremove("legend") 

p2 = summary_df_sample %>% ggplot(aes(x = mean_ari_sample, y = mean_ari_refineCT, label = method, color = method)) +
  geom_point(aes(color = method), size = 4) + # Dots for the mean values
  geom_errorbar(aes(ymin = mean_ari_refineCT - sd_ari_refineCT, ymax = mean_ari_refineCT + sd_ari_refineCT), 
                width = 0.001, color = 'grey') + # Vertical error bars
  geom_errorbarh(aes(xmin = mean_ari_sample - sd_ari_sample, xmax = mean_ari_sample + sd_ari_sample), 
                 height = 0.01, color = 'grey') + # Horizontal error bars
  xlab("ARI of Batch Effect (Sample)") +
  ylab("ARI of Cell Type Separation (Refined)") + 
  geom_text_repel(size = 5) +
  theme_bw(base_size = 30) +
  scale_color_manual(values = fillColor) +
  rremove("legend")

ggarrange(p1, p2, ncol = 2, nrow = 1)


dev.off()

seed = 20220524
generate_index = function(seed){
    set.seed(seed)
    subsample_index = c()
    for(study_each in unique(adt_feature$study_name)){
        study_index = which(adt_feature$study_name == study_each)
        if(length(study_index) < 5000){
            subsample_index = c(subsample_index, study_index)
        }else(
            subsample_index = c(subsample_index, sample(study_index, 5000, replace = FALSE))
        )
    }
    print(length(subsample_index))

    return(subsample_index)
}
subsample_index = generate_index(seed)
ct2_label = adt_feature$cell_type_l2[subsample_index] %>% as.character
cell_select = which(ct2_label != "undefined")
ct2_sublabel = ct2_label[cell_select]

ct1_label = adt_feature$cell_type_l1[subsample_index] %>% as.character
cell_select = which(ct1_label != "undefined")
ct1_sublabel = ct1_label[cell_select]

si_dist = c()
for(method in select_methods){
  tmp = readRDS(paste0(out_path, "/RDS/adt_", method, "_", run_name, "_umap_si.rds"))
  si_dist = data.frame(si_score = tmp[[1]][1, ], method = method, eval = "sample", label = adt_feature$sample[subsample_index]) %>% rbind(si_dist, .) 
  si_dist = data.frame(si_score = tmp[[2]][1, ], method = method, eval = "study", label = adt_feature$study_name[subsample_index]) %>% rbind(si_dist, .) 
  si_dist = data.frame(si_score = tmp[[3]][1, ], method = method, eval = "celltype_refined", label = ct2_sublabel) %>% rbind(si_dist, .) 
  si_dist = data.frame(si_score = tmp[[4]][1, ], method = method, eval = "celltype_broad", label = ct1_sublabel) %>% rbind(si_dist, .) 
}

si_dist %>% dplyr::filter(eval == "celltype_refined") %>% ggplot(aes(x = label, y = si_score, fill = label)) +
geom_boxplot() +
facet_wrap(~method, scales = "free") +
theme_bw(base_size = 30) +
rremove("legend")

si_dist %>% dplyr::filter(eval == "celltype_refined") %>% ggplot(aes(x = method, y = si_score, fill = method)) +
geom_boxplot() +
# facet_wrap(~label) +
theme_bw(base_size = 30) +
rremove("legend") +
rotate_x_text(angle = 45)

si_norm = function(si_score, min_val = -1){

  max_val = 1
  si_norm = (si_score - min_val)/(max_val - min_val)
  return(si_norm)

}

si_summary = c()
for(method in select_methods){
    
  tmp = readRDS(paste0(out_path, "/RDS/adt_", method, "_", run_name, "_umap_si.rds"))
  si_summary = data.frame(
    si_sample = rowMedians(tmp[[1]]) %>% as.vector %>% si_norm(., min_val = -0.7),
    si_site = rowMedians(tmp[[2]]) %>% as.vector %>% si_norm(., min_val = -0.4),
    si_refineCT = rowMedians(tmp[[3]]) %>% as.vector %>% si_norm(),
    si_broadCT = rowMedians(tmp[[4]]) %>% as.vector %>% si_norm(),
    method = method
  ) %>% rbind(si_summary, .)
  
}

si_all = si_summary %>% data.frame
si_all$method = factor(si_all$method, levels = select_methods)

# Calculate mean and standard deviation for each group
summary_df <- si_all %>%
  group_by(method) %>%
  summarise(
    mean_si_sample = mean(si_sample),
    sd_si_sample = sd(si_sample),
    mean_si_refineCT = mean(si_refineCT),
    sd_si_refineCT = sd(si_refineCT),
    mean_si_site = mean(si_site),
    sd_si_site = sd(si_site),
    mean_si_broadCT = mean(si_broadCT),
    sd_si_broadCT = sd(si_broadCT)
  )


pdf(paste0(out_path, "/Figures/Evaluation_metrics_Si_median_withDecontPro.pdf"), width = 16, height = 11)

rename_methods = c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony(RawCount)", "Harmony(Arcsinh)", "Harmony(CLR)", "Harmony(logCPM)", "fastMNN", "CytofRUV", "DSB", "decontPro", "totalVI", "sciPENN", "ADTnorm(default)", "ADTnorm(customized)")
summary_df_study = summary_df %>% dplyr::filter(method %in% c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount_study", "Harmony_Arcsinh_study", "Harmony_CLR_study", "Harmony_logCPM_study", "fastMNN_study", "CytofRUV_study", "DSB", "decontPro", "totalVI_study_GPU", "sciPENN_study_GPU", "ADTnorm_study", "ADTnorm_study_manual_keepZero")) %>% data.frame
summary_df_study$method = factor(summary_df_study$method, levels = c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount_study", "Harmony_Arcsinh_study", "Harmony_CLR_study", "Harmony_logCPM_study", "fastMNN_study", "CytofRUV_study", "DSB", "decontPro", "totalVI_study_GPU", "sciPENN_study_GPU", "ADTnorm_study", "ADTnorm_study_manual_keepZero"), labels = rename_methods)

summary_df_sample = summary_df %>% dplyr::filter(method %in% c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount_sample", "Harmony_Arcsinh_sample", "Harmony_CLR_sample", "Harmony_logCPM_sample", "fastMNN_sample", "CytofRUV_sample", "DSB", "decontPro", "totalVI_sample_GPU", "sciPENN_sample_GPU", "ADTnorm_study", "ADTnorm_sample_manual_keepZero")) %>% data.frame
summary_df_sample$method = factor(summary_df_sample$method, levels = c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount_sample", "Harmony_Arcsinh_sample", "Harmony_CLR_sample", "Harmony_logCPM_sample", "fastMNN_sample", "CytofRUV_sample", "DSB", "decontPro", "totalVI_sample_GPU", "sciPENN_sample_GPU", "ADTnorm_sample_manual_keepZero", "ADTnorm_study"), labels = rename_methods)

fillColor <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(rename_methods)) %>% rev

p1 = summary_df_study %>% ggplot(aes(x = mean_si_site, y = mean_si_broadCT, label = method, color = method)) +
  geom_point(aes(color = method), size = 4) + # Dots for the mean values
  geom_errorbar(aes(ymin = mean_si_broadCT - sd_si_broadCT, ymax = mean_si_broadCT + sd_si_broadCT), 
                width = 0.001, color = 'grey') + # Vertical error bars
  geom_errorbarh(aes(xmin = mean_si_site - sd_si_site, xmax = mean_si_site + sd_si_site), 
                 height = 0.01, color = 'grey') + # Horizontal error bars
  xlab("Silhouette of Batch Effect (Study)") +
  ylab("Silhouette of Cell Type Separation (Broad)") + 
  geom_text_repel(size = 5) +
  theme_bw(base_size = 30) +
  scale_color_manual(values = fillColor) +
  rremove("legend")

p2 = summary_df_sample %>% ggplot(aes(x = mean_si_sample, y = mean_si_refineCT, label = method, color = method)) +
  geom_point(aes(color = method), size = 4) + # Dots for the mean values
  geom_errorbar(aes(ymin = mean_si_refineCT - sd_si_refineCT, ymax = mean_si_refineCT + sd_si_refineCT), 
                width = 0.001, color = 'grey') + # Vertical error bars
  geom_errorbarh(aes(xmin = mean_si_sample - sd_si_sample, xmax = mean_si_sample + sd_si_sample), 
                 height = 0.01, color = 'grey') + # Horizontal error bars
  xlab("Silhouette of Batch Effect (Sample)") +
  ylab("Silhouette of Cell Type Separation (Refined)") + 
  geom_text_repel(size = 5) +
  theme_bw(base_size = 30) +
  scale_color_manual(values = fillColor) +
  rremove("legend")

ggarrange(p1, p2, ncol = 2, nrow = 1)

dev.off()

## LISI
seed = 20220524
generate_index = function(seed){
    set.seed(seed)
    subsample_index = c()
    for(study_each in unique(adt_feature$study_name)){
        study_index = which(adt_feature$study_name == study_each)
        if(length(study_index) < 5000){
            subsample_index = c(subsample_index, study_index)
        }else(
            subsample_index = c(subsample_index, sample(study_index, 5000, replace = FALSE))
        )
    }
    print(length(subsample_index))

    return(subsample_index)
}
subsample_index = generate_index(seed)
ct2_label = adt_feature$cell_type_l2[subsample_index] %>% as.character
cell_select = which(ct2_label != "undefined")
ct2_sublabel = ct2_label[cell_select]

ct1_label = adt_feature$cell_type_l1[subsample_index] %>% as.character
cell_select = which(ct1_label != "undefined")
ct1_sublabel = ct1_label[cell_select]

lisi_dist = c()
for(method in select_methods){ #
  tmp = readRDS(paste0(out_path, "/RDS/adt_", method, "_", run_name, "_pca_lisi.rds"))
  lisi_dist = data.frame(lisi_score = tmp[[1]][1, ], method = method, eval = "sample", label = adt_feature$sample[subsample_index]) %>% rbind(lisi_dist, .) 
  lisi_dist = data.frame(lisi_score = tmp[[2]][1, ], method = method, eval = "study", label = adt_feature$study_name[subsample_index]) %>% rbind(lisi_dist, .) 
  lisi_dist = data.frame(lisi_score = tmp[[3]][1, ], method = method, eval = "celltype_refined", label = adt_feature$cell_type_l2[subsample_index]) %>% rbind(lisi_dist, .) 
  lisi_dist = data.frame(lisi_score = tmp[[4]][1, ], method = method, eval = "celltype_broad", label = adt_feature$cell_type_l1[subsample_index]) %>% rbind(lisi_dist, .) 
  lisi_dist = data.frame(lisi_score = tmp[[5]][1, ], method = method, eval = "celltype_refined_noundefined", label = ct2_sublabel) %>% rbind(lisi_dist, .) 
  lisi_dist = data.frame(lisi_score = tmp[[6]][1, ], method = method, eval = "celltype_broad_noundefined", label = ct1_sublabel) %>% rbind(lisi_dist, .) 
}


library(ggridges)
lisi_dist$method = factor(lisi_dist$method, levels = select_methods)

pdf(paste0(out_path, "/Figures/Evaluation_metrics_lisi_distribution_withDecontPro.pdf"), width = 15, height = 11)

lisi_dist %>% dplyr::filter(eval %in% c("sample", "study")) %>% ggplot(aes(x = lisi_score, fill= method, y = method)) +
geom_density_ridges() +
facet_grid(~eval, scales = "free_x") +
theme_bw(base_size = 30) +
scale_fill_manual(breaks = methods, values = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(select_methods))) +
# coord_cartesian(xlim = c(0.5, 3)) +
# scale_x_continuous(breaks = 0:8) +
xlab("LISI of Batch Effect")+
rremove("legend") +
ylab("")

lisi_dist %>% dplyr::filter(eval %in% c("celltype_refined", "celltype_broad")) %>% ggplot(aes(x = lisi_score, fill= method, y = method)) +
geom_density_ridges() +
facet_grid(~eval, scales = "free_x") +
theme_bw(base_size = 30) +
scale_fill_manual(breaks = methods, values = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(select_methods))) +
# coord_cartesian(xlim = c(0.5, 3)) +
# scale_x_continuous(breaks = 0:8) +
xlab("LISI of Cell Type Separation")+
rremove("legend") +
ylab("") +
xlim(0, 3)


lisi_dist %>% dplyr::filter(eval %in% c("celltype_refined_noundefined", "celltype_broad_noundefined")) %>% ggplot(aes(x = lisi_score, fill= method, y = method)) +
geom_density_ridges() +
facet_grid(~eval, scales = "free_x") +
theme_bw(base_size = 30) +
scale_fill_manual(breaks = methods, values = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(select_methods))) +
# coord_cartesian(xlim = c(0.5, 3)) +
# scale_x_continuous(breaks = 0:8) +
xlab("LISI of Cell Type Separation")+
rremove("legend") +
ylab("") +
xlim(0, 3)
dev.off()


get_mode = function(in_matrix, adjust = 3){
  mode_list = c()
  for(i in 1:nrow(in_matrix)){
    den = density(in_matrix[i,], adjust = adjust)
    mode_list = c(mode_list, den$x[which.max(den$y)])
  }
  return(mode_list)
}



lisi_summary = c()
for(method in select_methods){
    
  tmp = readRDS(paste0(out_path, "/RDS/adt_", method, "_", run_name, "_pca_lisi.rds"))
  lisi_summary = data.frame(
    lisi_sample = tmp[[1]] %>% get_mode(adjust = 3), #rowMedians(tmp[[1]]) %>% as.vector,
    lisi_site = tmp[[2]] %>% get_mode(adjust = 3), #rowMedians(tmp[[2]]) %>% as.vector,
    lisi_refineCT = tmp[[3]] %>% get_mode(adjust = 3), #rowMedians(tmp[[3]]) %>% as.vector,
    lisi_broadCT = tmp[[4]] %>% get_mode(adjust = 3), #rowMedians(tmp[[4]]) %>% as.vector,
    method = method
  ) %>% rbind(lisi_summary, .)
  
}

lisi_all = lisi_summary %>% data.frame
lisi_all$method = factor(lisi_all$method, levels = select_methods)

# Calculate mean and standard deviation for each group
summary_df <- lisi_all %>%
  group_by(method) %>%
  summarise(
    mean_lisi_sample = mean(lisi_sample),
    sd_lisi_sample = sd(lisi_sample),
    mean_lisi_refineCT = mean(lisi_refineCT),
    sd_lisi_refineCT = sd(lisi_refineCT),
    mean_lisi_site = mean(lisi_site),
    sd_lisi_site = sd(lisi_site),
    mean_lisi_broadCT = mean(lisi_broadCT),
    sd_lisi_broadCT = sd(lisi_broadCT)
  )


pdf(paste0(out_path, "/Figures/Evaluation_metrics_lisi_peakmode.pdf"), width = 15, height = 11)

rename_methods = c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony(RawCount)", "Harmony(Arcsinh)", "Harmony(CLR)", "Harmony(logCPM)", "fastMNN", "CytofRUV", "DSB", "decontPro", "totalVI", "sciPENN", "ADTnorm(default)", "ADTnorm(customized)")
summary_df_study = summary_df %>% dplyr::filter(method %in% c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount_study", "Harmony_Arcsinh_study", "Harmony_CLR_study", "Harmony_logCPM_study", "fastMNN_study", "CytofRUV_study", "DSB", "decontPro", "totalVI_study_GPU", "sciPENN_study_GPU", "ADTnorm_study", "ADTnorm_study_manual_keepZero")) %>% data.frame
summary_df_study$method = factor(summary_df_study$method, levels = c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount_study", "Harmony_Arcsinh_study", "Harmony_CLR_study", "Harmony_logCPM_study", "fastMNN_study", "CytofRUV_study", "DSB", "decontPro", "totalVI_study_GPU", "sciPENN_study_GPU", "ADTnorm_study", "ADTnorm_study_manual_keepZero"), labels = rename_methods)

summary_df_sample = summary_df %>% dplyr::filter(method %in% c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount_sample", "Harmony_Arcsinh_sample", "Harmony_CLR_sample", "Harmony_logCPM_sample", "fastMNN_sample", "CytofRUV_sample", "DSB", "decontPro","totalVI_sample_GPU", "sciPENN_sample_GPU", "ADTnorm_sample", "ADTnorm_sample_manual_keepZero")) %>% data.frame
summary_df_sample$method = factor(summary_df_sample$method, levels = c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount_sample", "Harmony_Arcsinh_sample", "Harmony_CLR_sample", "Harmony_logCPM_sample", "fastMNN_sample", "CytofRUV_sample", "DSB", "decontPro", "totalVI_sample_GPU", "sciPENN_sample_GPU",  "ADTnorm_sample_manual_keepZero", "ADTnorm_sample"), labels = rename_methods)

fillColor <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(rename_methods)) %>% rev

p1 = summary_df_study %>% ggplot(aes(x = mean_lisi_site, y = mean_lisi_broadCT, label = method, color = method)) +
  geom_point(aes(color = method), size = 4) + # Dots for the mean values
  geom_errorbar(aes(ymin = mean_lisi_broadCT - sd_lisi_broadCT, ymax = mean_lisi_broadCT + sd_lisi_broadCT), 
                width = 0.001, color = 'grey') + # Vertical error bars
  geom_errorbarh(aes(xmin = mean_lisi_site - sd_lisi_site, xmax = mean_lisi_site + sd_lisi_site), 
                 height = 0.001, color = 'grey') + # Horizontal error bars
  xlab("LISI of Batch Effect (Study)") +
  ylab("LISI of Cell Type Separation (Broad)") + 
  geom_text_repel(size = 5) +
  theme_bw(base_size = 30) +
  scale_color_manual(values = fillColor) +
  rremove("legend")

p2 = summary_df_sample %>% ggplot(aes(x = mean_lisi_sample, y = mean_lisi_refineCT, label = method, color = method)) +
  geom_point(aes(color = method), size = 4) + # Dots for the mean values
  geom_errorbar(aes(ymin = mean_lisi_refineCT - sd_lisi_refineCT, ymax = mean_lisi_refineCT + sd_lisi_refineCT), 
                width = 0.001, color = 'grey') + # Vertical error bars
  geom_errorbarh(aes(xmin = mean_lisi_sample - sd_lisi_sample, xmax = mean_lisi_sample + sd_lisi_sample), 
                 height = 0.001, color = 'grey') + # Horizontal error bars
  xlab("LISI of Batch Effect (Sample)") +
  ylab("LISI of Cell Type Separation (Refined)") + 
  geom_text_repel(size = 5) +
  theme_bw(base_size = 30) +
  scale_color_manual(values = fillColor) +
  rremove("legend")
ggarrange(p1, p2, ncol = 2, nrow = 1)

dev.off()




## ====================================================================================
## 5. missing marker density plot for: CLR, logCPM, Arcsinh, ADTnorm, sciPENN, totalVI
## ====================================================================================

# run_name <- "publicData_CITEseq_12buusT"
# adt_data = adt_data_full %>% data.frame

# methods <- c("Arcsinh_b5", "CLR", "logCPM", "Arcsinh_b5_CLR", "HarmonyRawCount_sample",  "HarmonyRawCount_study", "HarmonyArcsinh_b5_sample",  "HarmonyArcsinh_b5_study", "HarmonyCLR_sample", "HarmonyCLR_study", "HarmonylogCPM_sample",  "HarmonylogCPM_study","CytofRUV_study", "CytofRUV_sample", "fastMNN_study_log", "fastMNN_sample_log", "DSB", "ADTnorm") #names(normalized_method)
# for (method in methods[c(1:2)]) {
#   print(method)
#   filename <- paste0(out_path, "/RDS/adt_data_norm_", method, "_allMarker_", study_name, ".rds")
#   norm_f <- normalized_method[[method]][["f"]]
#   start_adt_method <- normalized_method[[method]][["start_adt_method"]]
#   method_param <- normalized_method[[method]][["param"]]

#   if(start_adt_method == "RawCount"){
#     cell_x_adt <- adt_data
#   }else{
#     cell_x_adt <- readRDS(file = paste0(out_path, "/RDS/adt_data_norm_", start_adt_method, "_", study_name, ".rds"))
#   }
#   if(stringr::str_detect(method, "CytofRUV")){
#     cell_x_adt_colnames = colnames(cell_x_adt)
#     colnames(cell_x_adt) = paste0("tmpName", 1:ncol(cell_x_adt))
#     setwd("~")
#   }
  
#   ans <- do.call(
#     norm_f,
#     list(cell_x_adt = cell_x_adt, cell_x_feature = adt_feature, parameter_list = method_param)
#   )
#   if(stringr::str_detect(method, "CytofRUV")){
#     colnames(ans) = cell_x_adt_colnames
#   }

#   saveRDS(ans, file = filename)

# }

# marker_list = c("CD5", "CD11B", "CD11C", "CD15", "CD16", "CD27", "CD45RO",  "CD62L", "CD197", "PD1")

# bw_list = list("Arcsinh_b5" = 0.05, "CLR" = 0.05, "logCPM" = 0.1, "Arcsinh_b5_CLR" = 0.003, "HarmonyRawCount_sample" = 0.05, "HarmonyRawCount_study" = 0.05, "HarmonyArcsinh_b5_sample" = 0.05, "HarmonyArcsinh_b5_study" = 0.05, "HarmonyCLR_sample" = 0.03, "HarmonyCLR_study" = 0.03,  "HarmonylogCPM_sample" = 0.05, "HarmonylogCPM_study" = 0.05,"CytofRUV_study" = 0.05, "CytofRUV_sample" = 0.05, "fastMNN_study_log" = 0.005,  "fastMNN_sample_log" = 0.005,"DSB" = 0.1, "ADTnorm" = 0.1)
# trans_list = list("Arcsinh_b5" = FALSE, "CLR" = FALSE, "logCPM" = TRUE, "Arcsinh_b5_CLR" = FALSE, "HarmonyRawCount_sample" = TRUE, "HarmonyRawCount_study" = TRUE, "HarmonyArcsinh_b5_sample" = FALSE, "HarmonyArcsinh_b5_study" = FALSE, "HarmonyCLR_sample" = FALSE, "HarmonyCLR_study" = FALSE, "HarmonylogCPM_sample" = TRUE, "HarmonylogCPM_study" = TRUE, "CytofRUV_study" = FALSE, "CytofRUV_sample" = FALSE, "fastMNN_study_log" = FALSE, "fastMNN_sample_log" = FALSE, "DSB" = FALSE, "ADTnorm" = FALSE)
# for(method in methods[c(18)]){
#   filename <- paste0(out_path, "/RDS/adt_data_norm_", method, "_allMarker_", study_name, ".rds")
#   ans = readRDS(file = filename)

#   pdf(paste0(out_path, "/figures/", method, "_missingMarker_adt_density.pdf"), width = 100, height = 20)
#   print(plot_adt_density(
#   cell_x_adt = ans %>% data.frame, 
#   cell_x_feature = adt_feature, 
#   adt_marker_list = marker_list, 
#   parameter_list = list(method_label = method, arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
#   ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list))))

#   dev.off()

# }


## ================================================================
## 6. Abnormal/Arbitrary manipulation of normalized protein counts
## ================================================================
select_methods <- c("Arcsinh", "CLR", "logCPM", "Arcsinh_CLR", "Harmony_RawCount_study", "Harmony_RawCount_sample", "Harmony_Arcsinh_study", "Harmony_Arcsinh_sample", "Harmony_CLR_study", "Harmony_CLR_sample", "Harmony_logCPM_study", "Harmony_logCPM_sample", "fastMNN_study", "fastMNN_sample",  "CytofRUV_study", "CytofRUV_sample", "DSB","totalVI_sample_GPU", "totalVI_study_GPU", "sciPENN_sample_GPU", "sciPENN_study_GPU", "ADTnorm_sample", "ADTnorm_study", "ADTnorm_study_manual_keepZero", "ADTnorm_sample_manual_keepZero") 
mess_select_methods =  c("Arcsinh", "CLR", "Harmony_CLR_study", "fastMNN_study", "CytofRUV_study", "DSB", "totalVI_study_GPU", "sciPENN_study_GPU", "ADTnorm_study_manual_keepZero")

data_name = "10X_malt_10k"

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

## Try on more T cell only datasets

### Manual tuning for shiny App demonstration
library(shiny)
method = "ADTnorm_study"
adt_feature$sample = factor(adt_feature$study_name, levels = file_list)
adt_feature$batch = factor(adt_feature$study_name, levels = file_list)
# override_landmark = NULL
# all_landmarks = load_landmarks(paste0(out_path, "/manualTuning/"))
# override_landmark = all_landmarks

# marker_to_process = "CD305"
marker_to_process = "CD3"
res_norm = ADTnorm(
    marker_to_process = marker_to_process,
    target_landmark_location = NULL, #c(1, 3),
    cell_x_adt = adt_data, 
    cell_x_feature = adt_feature, 
    customize_landmark = TRUE,
    multi_sample_per_batch = FALSE,
    save_fig = FALSE,
    save_landmark = FALSE,
    save_outpath = paste0(out_path, "/test/"),
    study_name = run_name,
    trimodal_marker = c("CD4", "CD45RA"),
    bw_smallest_bi = 1.1,
    bw_smallest_tri = 1.1,
    bw_smallest_adjustments = list(CD3 = 0.8, CD4 = 0.8, CD8 = 0.8),
    shoulder_valley_slope = -1,
    exclude_zeroes = FALSE,
    bimodal_marker = NULL,
    positive_peak = NULL,
    quantile_clip = 1,
    peak_type = "mode",
    shoulder_valley = TRUE,
    valley_density_adjust = 3,
    landmark_align_type = "negPeak_valley_posPeak",
    midpoint_type = "valley",
    neg_candidate_thres = asinh(5/5 + 1),
    lower_peak_thres = 0.0001,
    brewer_palettes = "Set1",
    detect_outlier_valley = FALSE,
    clean_adt_name = FALSE,
    override_landmark = NULL,
    verbose = TRUE)

