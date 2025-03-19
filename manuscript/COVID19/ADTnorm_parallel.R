#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

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
library(ADTnorm)
# p_load(flowStats, flowCore, FlowSOM, ncdfFlow, flowViz, pdfCluster, cluster)

## =====================
## study name and paths
## =====================
run_name = "COVID19"
master_path = "./"
out_path = paste0(master_path, "manuscript/results/", run_name)
# fcs_path = paste0(out_path, "/FCS/")
fig_path = paste0(out_path, "/Figures/", run_name)

adt_data = readRDS(file = paste0(out_path, "/RDS/adt_data_", run_name, ".rds"))
adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", run_name, ".rds")) %>% data.frame

arcsinh_param <- list(a = 1, b = 1 / 5, c = 0)

method = "ADTnorm_sample_manual"
adt_feature$sample = factor(adt_feature$site_sample)
res_norm = ADTnorm(
    cell_x_adt = adt_data, 
    cell_x_feature = adt_feature, 
    marker_to_process = args[1], 
    customize_landmark = FALSE,
    save_fig = TRUE,
    save_landmark = TRUE,
    save_outpath = out_path,
    study_name = run_name,
    trimodal_marker = c("CD4", "CD45RA", "CD49f", "CD62L"),
    bw_smallest_bi = 1.1,
    bw_smallest_tri = 1.1,
    bw_smallest_adjustments = list(CD3 = 1.1, CD4 = 1.1, CD8 = 1.1),
    shoulder_valley_slope = -1,
    exclude_zeroes = FALSE,
    bimodal_marker = NULL,
    positive_peak = NULL,
    quantile_clip = 1,
    peak_type = "mode",
    multi_sample_per_batch = TRUE,
    shoulder_valley = TRUE,
    valley_density_adjust = 3,
    landmark_align_type = "negPeak_valley_posPeak",
    midpoint_type = "valley",
    neg_candidate_thres = asinh(3/5 + 1),
    lower_peak_thres = 0.005,
    brewer_palettes = "Set1",
    detect_outlier_valley = FALSE,
    target_landmark_location = NULL,
    clean_adt_name = FALSE,
    override_landmark = NULL,
    verbose = TRUE)


saveRDS(res_norm, file = paste0(out_path, "/RDS/adt_data_norm_", method, "_", args[1], "_", run_name, ".rds"))

