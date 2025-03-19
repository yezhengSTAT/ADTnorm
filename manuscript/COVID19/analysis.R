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
# p_load(flowStats, flowCore, FlowSOM, ncdfFlow, flowViz, pdfCluster, cluster)

## =====================
## study name and paths
## =====================
run_name = "COVID19"
master_path = "./"
out_path = paste0(master_path, "manuscript/results/", run_name)
# fcs_path = paste0(out_path, "/FCS/")
fig_path = paste0(out_path, "/Figures/", run_name)

## ======================================
## Clean and organize titration datasets
## ======================================
sce = readH5AD("./data/covid_portal_210320_with_raw.h5ad")
adt_data = sce@assays@data@listData$raw[24738:24929, ]
adt_czi_obj = readRDS("./data/10_1038_s41591_021_01329_2_protein.rds")

rownames(adt_data) = rownames(adt_czi_obj) %>% gsub("\\(.*", "", .) %>% gsub("\\.", "", .) %>% gsub("_", "", .) %>% gsub("-", "", .) %>% gsub("/", "", .) %>% gsub(" ", "", .) %>% gsub(",", "", .) %>% gsub("α", "alpha", .) %>% gsub("κ", "kappa", .) %>% gsub("λ", "lambda", .) %>% gsub("β", "beta", .) %>% gsub("γ", "gamma", .) %>% gsub("δ", "delta", .) 
colnames(adt_data) = colnames(sce)
adt_data = adt_data %>% t() 

rna_data = sce@assays@data@listData$raw[1:24737, ]
rownames(rna_data) = rownames(sce)[1:24737]
colnames(rna_data) = colnames(sce)

adt_czi_obj_metadata = adt_czi_obj@meta.data
colnames(adt_czi_obj_metadata) = paste0("CZI_", colnames(adt_czi_obj_metadata))
adt_feature = cbind(colData(sce), adt_czi_obj_metadata)
adt_feature$sample = factor(adt_feature$sample_id)
adt_feature$batch = factor(adt_feature$Site)
adt_feature$site_sample = paste0(adt_feature$Site, "_", adt_feature$sample)
adt_feature$adt_depth = rowSums(adt_data)
adt_feature$adt_nonzero = rowSums(adt_data > 0)
adt_feature$adt_nonzero_prop = adt_feature$adt_nonzero / ncol(adt_data) * 100

rm_ind1 = which(adt_feature$adt_depth <= 500) ## according to the density plot
rm_ind2 = which(adt_feature$adt_nonzero_prop <= 70) ## according to the density plot
adt_data = adt_data[-union(rm_ind1, rm_ind2), ] 
adt_feature = adt_feature[-union(rm_ind1, rm_ind2), ]
rna_data = rna_data[-union(rm_ind1, rm_ind2), ]
saveRDS(union(rm_ind1, rm_ind2), file = paste0(out_path, "/RDS/rm_ind1_ind2_", run_name, ".rds"))
saveRDS(adt_data, file = paste0(out_path, "/RDS/adt_data_", run_name, ".rds"))
saveRDS(rna_data, file = paste0(out_path, "/RDS/rna_data_", run_name, ".rds"))
saveRDS(adt_feature, file = paste0(out_path, "/RDS/adt_feature_", run_name, ".rds"))

adt_dsb = sce@assays@data@listData$X[(24929 - 191):24929, -rm_ind] %>% t()
colnames(adt_dsb) = colnames(adt_data)
rownames(adt_dsb) = colnames(sce)[-rm_ind]
saveRDS(adt_dsb, file = paste0(out_path, "/RDS/adt_data_norm_paperOriginal_", run_name, ".rds"))
adt_dsb = readRDS(file = paste0(out_path, "/RDS/adt_data_norm_paperOriginal_", run_name, ".rds"))

## fastMNN in original paper
# adt = adt_dsb %>% t()
# adt1 <- adt[, adt_feature$Site=="Cambridge"]
# adt2 <- adt[, adt_feature$Site=="Ncl"]
# adt3 <- adt[, adt_feature$Site=="Sanger"]
# param <- MulticoreParam(workers=4)
# print("Starting MNN")
# set.seed(300)
# mnncor <- batchelor::fastMNN(adt2,adt1,adt3,
# 		     BPPARAM=param, # commented this out for singularity
# 		     k=20,
# 		     d=50,
# 		     merge.order=c(1,2,3),
# 		     BSPARAM=IrlbaParam(deferred=TRUE),
# 		     assay.type="counts",
# 		     #                      BNPARAM=AnnoyParam(),
# 		     cos.norm=TRUE, # prob necesarry
# 		     correct.all=TRUE)
# print("Done MNN")

# adt_fastmnn <- assay(mnncor,"reconstructed")
# adt_fastmnn <- adt_fastmnn[, colnames(adt)] %>% t()
# saveRDS(mnncor, file = paste0(out_path, "/RDS/mnncor_", run_name, ".rds"))
# saveRDS(adt_fastmnn, file = paste0(out_path, "/RDS/adt_data_norm_fastMNN_", run_name, ".rds"))


adt_data = readRDS(file = paste0(out_path, "/RDS/adt_data_", run_name, ".rds"))
# rna_data = readRDS(file = paste0(out_path, "/RDS/rna_data_", run_name, ".rds")) 
adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", run_name, ".rds"))
rm_ind = readRDS(file = paste0(out_path, "/RDS/rm_ind1_ind2_", run_name, ".rds"))

adt_feature$CZI_initial_clustering = factor(adt_feature$CZI_initial_clustering, levels = c("CD4", "Treg", "CD8", "gdT", "MAIT", "NK_16hi",  "NK_56hi", "B_cell", "Plasmablast",  "CD14", "CD16", "DCs", "pDC", "HSC",  "Platelets", "RBC", "Lymph_prolif",  "Mono_prolif")
, labels = c("CD4 T", "Treg", "CD8 T", "gamma delta T", "MAIT", "NK CD56-",  "NK 56+", "B", "Plasmablast",  "Monocytes CD14+", "Monocytes CD16+", "cDC", "pDC", "HSPC",  "Platelets", "RBC", "Prolif lymph",  "Prolif Mono"))

## =================
## RNA UMAP
## =================
sce = readH5AD("./data/covid_portal_210320_with_raw.h5ad")
rna_umap = sce@int_colData@listData$reducedDims@listData$X_umap %>% data.frame
rna_umap
colnames(rna_umap) = c("X1", "X2")
saveRDS(rna_umap, file = paste0(out_path, "/RDS/rna_umap_", run_name, ".rds"))

pdf(paste0(out_path, "/Figures/rna_umap_raster.pdf"), width = 13, height = 13)
print(plot_umap_extend_raster(rna_umap[-rm_ind, ], adt_feature, point_size = 0.3, parameter_list = list(
    target_feature = "full_clustering", 
    color_design = c(brewer.pal(8, "Set1"), carto_pal(11, "Bold"), brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"), brewer.pal(12, "Set3")),
    method_label = method
))
)

print(plot_umap_extend_raster(rna_umap[-rm_ind, ], adt_feature, point_size = 0.3, parameter_list = list(
    target_feature = "CZI_initial_clustering", 
    color_break = c("CD4 T", "Treg", "CD8 T", "gamma delta T", "MAIT", "NK CD56-",  "NK 56+", "B", "Plasmablast",  "Monocytes CD14+", "Monocytes CD16+", "cDC", "pDC", "HSPC",  "Platelets", "RBC", "Prolif lymph",  "Prolif Mono"),
    color_design = c("#AADBDF", "#E599E8", "#ABDFAE", "#CC5C5D", "#EAA6AC", "#E4E7A8",  "#385F9E", "#D97E4A", "#4CB45E",  "#ABA9D9", "#F8BD98", "#67A6E4", "#D9DFA9", "#C487EB",  "#F6B6CE", "#E43324", "#6FC6CA",  "#5F9CC1"),
    method_label = method
))
)

print(plot_umap_extend_raster(rna_umap[-rm_ind, ], adt_feature, point_size = 0.3, parameter_list = list(
    target_feature = "CZI_cell_type", 
    color_design = c(brewer.pal(8, "Set1"), carto_pal(11, "Bold"), brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"), brewer.pal(12, "Set3")),
    method_label = method
))
)
dev.off()

## ======================
## subsample the samples
## ======================
## select samples per condition per site
# paste0(adt_feature$Site, "_", adt_feature$sample)
sample_select = c(
  "Cambridge_BGCV14_CV0940", "Cambridge_BGCV10_CV0939", "Ncl_MH8919282", "Ncl_MH8919332", ## healthy, no healthy sample from sanger
  "Sanger_AP11", "Sanger_AP4", "Ncl_MH9143426", "Ncl_MH8919326", "Cambridge_BGCV04_CV0100", "Cambridge_BGCV08_CV0073", # mild
  "Sanger_AP6", "Sanger_AP3", "Cambridge_BGCV06_CV0234", "Cambridge_BGCV06_CV0037", "Ncl_MH9143270", "Ncl_MH9179824", # Moderate
  "Sanger_AP5", "Sanger_AP8", "Cambridge_BGCV12_CV0178", "Cambridge_BGCV06_CV0178", "Ncl_MH9143320", "Ncl_MH9143274", # severe
  "Sanger_AP2", "Sanger_AP12", "Cambridge_BGCV06_CV0201", "Cambridge_BGCV03_CV0200", "Ncl_MH8919328", "Ncl_newcastle21" # critical
)

sample_select_small = c(
  "Cambridge_BGCV14_CV0940", "Ncl_MH8919282", ## healthy, no healthy sample from sanger
  "Sanger_AP11",  "Ncl_MH9143426", "Cambridge_BGCV04_CV0100",  # mild
  "Sanger_AP6",  "Cambridge_BGCV06_CV0234",  "Ncl_MH9143270",  # Moderate
  "Sanger_AP5", "Cambridge_BGCV12_CV0178", "Ncl_MH9143320", # severe
  "Sanger_AP2", "Cambridge_BGCV06_CV0201", "Ncl_MH8919328" # critical
)

marker_list = c("CD1a", "CD1c", "CD1d", "CD2", "CD3", "CD4", "CD5", "CD7", "CD8", "CD10", "CD11a", "CD11b", "CD11c", "CD14", "CD15", "CD16", "CD18","CD19", "CD20", "CD21", "CD22", "CD23", "CD24", "CD25", "CD26", "CD27", "CD28", "CD29", "CD30", "CD31", "CD32", "CD33", "CD34", "CD35", "CD36", "CD38", "CD39", "CD40", "CD41", "CD44", "CD45", "CD45RA", "CD45RO", "CD47", "CD49b", "CD49d", "CD49f", "CD52", "CD54", "CD56", "CD57", "CD58", "CD62L", "CD62P", "CD64", "CD66ace", "CD66b", "CD69", "CD70", "CD71", "CD73", "CD79b", "CD80", "CD81", "CD82", "CD83", "CD85j", "CD86", "CD88", "CD90", "CD94", "CD95", "CD96", "CD98", "CD99", "CD101", "CD103", "CD106", "CD107a", "CD112", "CD117", "CD122", "CD123", "CD124", "CD127", "CD133","CD137", "CD137L", "CD138", "CD141", "CD144", "CD146", "CD150", "CD152", "CD154", "CD155", "CD158", "CD158b", "CD158e", "CD158f1", "CD161", "CD163", "CD169", "CD178", "CD183", "CD184", "CD185", "CD193", "CD194", "CD195", "CD196", "CD197", "CD204", "CD206", "CD207", "CD209", "CD223", "CD224", "CD226", "CD235ab", "CD244", "CD252", "CD254", "CD257", "CD258", "CD267", "CD268", "CD269", "CD272", "CD273", "CD274", "CD275", "CD278", "CD279", "CD294", "CD303", "CD304", "CD305", "CD307d", "CD307e", "CD309", "CD314", "CD319", "CD324", "CD326", "CD328", "CD335", "CD336", "CD337", "CD357", "CD360", "CD366", "CD370", "CD371", "B7H4", "cMet", "CX3CR1", "DR3", "FcεRIalpha", "GARP", "HLAA2", "HLAABC", "HLADR", "HLAF", "IgA", "IgD", "IgG1kappa", "IgG2akappa", "IgG2bkappa", "IgG2bkappaIsotypeCtrl", "IgGFc", "Iglightchainkappa", "Iglightchainlambda", "IgM", "Integrinbeta7", "KLRG1", "LOX1", "Mac2", "NLRP2", "PhosphoTau", "Podocalyxin", "Podoplanin", "TCR", "TCRValpha24Jalpha18", "TCRValpha72", "TCRVbeta131", "TCRVgamma9", "TCRVdelta2", "TCRgammadelta", "TIGIT", "TSLPR", "XCR1")

## Baseline Methods
arcsinh_param <- list(a = 1, b = 1 / 5, c = 0)
normalized_method <- list(
  Arcsinh = list(f = arcsinh_basic_transform, start_adt_method = "RawCount", param = arcsinh_param),
  CLR = list(f = clr_transform, start_adt_method = "RawCount", param = arcsinh_param),
  logCPM = list(f = log1pcpm_transform, start_adt_method = "RawCount"),
  ADTnorm_sample_default = list(
    f = ADTnorm_transform, start_adt_method = "RawCount", 
    param = list(
      unit = "sample",
      save_outpath = out_path,
      study_name = run_name,
      trimodal_marker = c("CD4", "CD45RA", "CD49f", "CD62L"),
      lower_peak_thres = 0.01,
      bw_smallest_bi = 1.1,
      bw_smallest_tri = 0.8,
      bw_smallest_adjustments = list(CD3 = 0.8, CD4 = 0.8, CD8 = 0.8),
      shoulder_valley_slope = -1,
      neg_candidate_thres = asinh(2/5 + 1)
    )
  ),
  Harmony_Arcsinh_sample = list(f = harmony_transform, start_adt_method = "Arcsinh", param = list(batch_column = "sample")),
  Harmony_CLR_sample = list(f = harmony_transform, start_adt_method = "CLR", param = list(batch_column = "sample"))
)

## run normalization
baseline_methods <- c("Arcsinh", "CLR", "logCPM", "Harmony_Arcsinh_sample", "Harmony_CLR_sample", "ADTnorm_sample_default", "paperOriginal", "fastMNN") # "ADTnorm_sample_manual", , 
file_list = adt_feature$batch %>% levels

for(method in baseline_methods[1:3]) {
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

plot_adt_density <- function(cell_x_adt, cell_x_feature, adt_marker_list, unit, parameter_list = NULL, ncol = 8){
    
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
bw_list = list("Arcsinh" = 0.08, "CLR" = 0.1, "CLR2" = 0.1, "logCPM" = 0.05,  "ADTnorm_sample_default" = 0.1, "ADTnorm_sample_manual" = 0.2, "Harmony_Arcsinh_sample" = 0.1, "Harmony_CLR_sample" = 0.1, "ADTnorm_sample_default" = 0.1, "ADTnorm_sample_manualtuning" = 0.1, "paperOriginal" = 0.1, "fastMNN" = 0.01)
trans_list = list("Arcsinh" = FALSE, "CLR" = FALSE, "CLR2" = FALSE, "logCPM" = TRUE, "ADTnorm_sample_default" = FALSE, "ADTnorm_sample_manual" = FALSE, "Harmony_Arcsinh_sample" = FALSE, "Harmony_CLR_sample" = FALSE, "ADTnorm_sample_default" = FALSE, "ADTnorm_sample_manualtuning" = FALSE, "paperOriginal" = FALSE, "fastMNN" = FALSE)

marker_list1 = marker_list[1:64]
marker_list2 = marker_list[65:128]
marker_list3 = marker_list[129:192]
show_sample = which(adt_feature$site_sample %in% sample_select_small)
for(method in baseline_methods[8]){
  print(method)
  filename <- paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds")
  ans = readRDS(file = filename)[show_sample, ]

  fig_height = 150
  pdf(paste0(out_path, "/Figures/", method, "_sample_", run_name, "_adt_density.pdf"), width = 100, height = fig_height)
  print(plot_adt_density(
  cell_x_adt = ans %>% data.frame, 
  cell_x_feature = adt_feature[show_sample, ], 
  adt_marker_list = marker_list1, 
  unit = "sample",
  parameter_list = list(method_label = method, arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
  ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list) + 1)))

  print(plot_adt_density(
  cell_x_adt = ans %>% data.frame, 
  cell_x_feature = adt_feature[show_sample, ], 
  adt_marker_list = marker_list2, 
  unit = "sample",
  parameter_list = list(method_label = method, arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
  ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list) + 1)))
  
  print(plot_adt_density(
  cell_x_adt = ans %>% data.frame, 
  cell_x_feature = adt_feature[show_sample, ], 
  adt_marker_list = marker_list3, 
  unit = "sample",
  parameter_list = list(method_label = method, arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
  ) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list) + 1)))
  
  dev.off()
}

## ADTnorm - run through ADTnorm_parallel.R for parallel running
## Concatenate the results from ADTnorm_parallel.R
res_norm_full = as.data.table(adt_data)
method = "ADTnorm_sample_default"
for(adt_marker_each in marker_list){
    if(file.exists(paste0(out_path, "/RDS/adt_data_norm_ADTnorm_sample_manual_", adt_marker_each, "_COVID19.rds"))){
        tmp = readRDS(paste0(out_path, "/RDS/adt_data_norm_ADTnorm_sample_manual_", adt_marker_each, "_COVID19.rds")) %>% as.vector
        set(res_norm_full, j = adt_marker_each, value = tmp)       
    }else{
        print(adt_marker_each)
    }

}
saveRDS(res_norm, file = paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds"))

method = "ADTnorm_sample_default"
res_norm = readRDS(file = paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds"))

for(adt_marker_each in marker_list){
    if(file.exists(paste0(out_path, "/RDS/adt_data_norm_ADTnorm_sample_manualtuning_", adt_marker_each, "_COVID19.rds"))){
        tmp = readRDS(paste0(out_path, "/RDS/adt_data_norm_ADTnorm_sample_manualtuning_", adt_marker_each, "_COVID19.rds")) %>% as.vector
        res_norm[[adt_marker_each]] = tmp    
    }else{
        print(adt_marker_each)
    }

}
method = "ADTnorm_sample_manualtuning"
# saveRDS(res_norm, file = paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds"))

res_norm = readRDS(file = paste0(out_path, "/RDS/adt_data_norm_", method, "_", run_name, ".rds")) %>% data.frame


fig_height = 150
pdf(paste0(out_path, "/Figures/", method, "_sample_", run_name, "_adt_density.pdf"), width = 100, height = fig_height)
print(plot_adt_density(
cell_x_adt = res_norm[show_sample, ] %>% data.frame, 
cell_x_feature = adt_feature[show_sample, ], 
adt_marker_list = marker_list1, 
unit = "sample",
parameter_list = list(method_label = method, arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list) + 1)))

print(plot_adt_density(
cell_x_adt = res_norm[show_sample, ] %>% data.frame, 
cell_x_feature = adt_feature[show_sample, ], 
adt_marker_list = marker_list2, 
unit = "sample",
parameter_list = list(method_label = method, arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list) + 1)))

print(plot_adt_density(
cell_x_adt = res_norm[show_sample, ] %>% data.frame, 
cell_x_feature = adt_feature[show_sample, ], 
adt_marker_list = marker_list3, 
unit = "sample",
parameter_list = list(method_label = method, arcsinhTransform = trans_list[[method]], bw = bw_list[[method]])
) + scale_fill_manual(breaks = file_list, values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(file_list) + 1)))

dev.off()


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
    # within_peak_sd = sqrt((sum((adt_tmp[which(adt_tmp < valley_info[each_sample, 1])] - peak_info[each_sample, 1])^2) + sum((adt_tmp[which(adt_tmp > valley_info[each_sample, 1])] - median_pos)^2))/length(adt_tmp))    
    # within_peak_sd = sqrt(sum((adt_tmp[which(adt_tmp < valley_info[each_sample, 1])] - peak_info[each_sample, 1])^2)/length(which(adt_tmp < valley_info[each_sample, 1])))    
    if(sum(adt_tmp == arcsinh_transform(0))/length(adt_tmp) > 0.98){
        within_peak_sd = 1
    }else{
        within_peak_sd = sqrt(sd(adt_tmp))
    }
    deep_scaler = valley_deep_scaler(cell_x_feature, cell_x_adt, adt_marker_select, each_sample, peak_info, valley_info, peak_num = peak_num)
    # return(mean_diff/within_peak_sd)
    return((mean_diff*deep_scaler)/within_peak_sd)
 

}

get_positive_cell_prop = function(cell_x_feature, cell_x_adt, adt_marker_select, each_sample, valley_info){
    ind = which(cell_x_feature$sample == each_sample)
    adt_tmp = cell_x_adt[ind, adt_marker_select]
    valley = valley_info[each_sample, 1]

    prop = sum(adt_tmp >= valley)/length(adt_tmp) * 100
    return(prop)
}

## use valley and peak
adt_feature$sample = adt_feature$site_sample

peak_valley_locations_list = list()
for(adt_marker_select in colnames(cell_x_adt)){ ##
    peak_valley_density = readRDS(paste0(out_path, "/manualTuning/RDS_afterTune/peak_valley_locations_", adt_marker_select, "_", run_name, ".rds"))
    peak_valley_locations_list[[adt_marker_select]] = peak_valley_density
}

peak_num_summary = c()
peak_sep_summary = c()
positive_prop_summary = c()
for(adt_marker_select in colnames(cell_x_adt)){ ##
    
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
    valley_info = peak_valley_density$valley_landmark_list

    for(each_sample in unique(adt_feature$sample)){
        batch_info = each_sample
        each_peak_info = peak_info[each_sample, ]
        peak_num = sum(is.na(each_peak_info) == FALSE)
        positive_prop_summary = data.frame(
            positive_prop = get_positive_cell_prop(cell_x_feature = adt_feature, cell_x_adt, adt_marker_select, each_sample, valley_info), 
            batch = batch_info, 
            sample = each_sample, adt_marker = adt_marker_select) %>% 
            rbind(positive_prop_summary, .)

    }
}
method = "ADTnorm_sample_manualtuning"
# saveRDS(list(peak_sep_summary, peak_num_summary, positive_prop_summary), file = paste0(out_path, "/RDS/peak_stain_quality_", run_name, "_", method, "_peak_sep_num_prop.rds"))
tmp = readRDS(file = paste0(out_path, "/RDS/peak_stain_quality_", run_name, "_", method, "_peak_sep_num_prop.rds"))
peak_sep_summary = tmp[[1]]
peak_num_summary = tmp[[2]]
positive_prop_summary = tmp[[3]]

peak_sep_summary$batch = peak_sep_summary$batch %>% gsub("_.*", "", .)
peak_sep_summary$batch = factor(peak_sep_summary$batch, levels = file_list)
peak_sep_summary$adt_marker = factor(peak_sep_summary$adt_marker, levels = marker_list)

peak_num_summary$batch = peak_num_summary$batch %>% gsub("_.*", "", .)
peak_num_summary$batch = factor(peak_num_summary$batch, levels = file_list)
peak_num_summary$adt_marker = factor(peak_num_summary$adt_marker, levels = marker_list)

positive_prop_summary$batch = positive_prop_summary$batch %>% gsub("_.*", "", .)
positive_prop_summary$batch = factor(positive_prop_summary$batch, levels = file_list)
positive_prop_summary$adt_marker = factor(positive_prop_summary$adt_marker, levels = marker_list)


peak2marker = peak_sep_summary$adt_marker[which(peak_sep_summary$peak_num != "# of peak: 1")] %>% table
peak2marker = which(peak2marker > 0) %>% names

peak2markerselect = c("CD1d", "CD2", "CD3", "CD4", "CD5", "CD7", "CD8", "CD11b", 
"CD11c", "CD14", "CD16", "CD19", "CD20", "CD21", 
"CD26", "CD27", "CD28", "CD31", "CD32", "CD33", "CD35", "CD36", 
"CD38", "CD39", "CD40", "CD41", "CD45RA", "CD45RO",  
"CD52", "CD54", "CD56", "CD58", "CD62L", "CD62P", "CD64", 
"CD71", "CD73", "CD88", "CD94", "CD95", 
"CD99", "CD101", "CD127", 
"CD158e", "CD169", "CD235ab", "CD244", 
"CD268", "CD269", "CD305", "CD328", "CD371", "HLADR")

positive_prop_summary$adt_marker %>% unique %>% length

pos_marker = positive_prop_summary %>% dplyr::filter(positive_prop > 30) %$% adt_marker %>% unique

## to save data to share with Daniel
cell_x_adt_posNeg = matrix("neg", nrow = nrow(cell_x_adt), ncol = ncol(cell_x_adt))
colnames(cell_x_adt_posNeg) = colnames(cell_x_adt)
for(adt_marker_select in colnames(cell_x_adt)){ ##
    
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
    valley_info = peak_valley_density$valley_landmark_list

    for(each_sample in unique(adt_feature$sample)){
        batch_info = each_sample
        each_peak_info = peak_info[each_sample, ]
        ind = which(adt_feature$sample == each_sample)
        adt_tmp = cell_x_adt[ind, adt_marker_select]
        valley = valley_info[each_sample, 1]

        cell_x_adt_posNeg[ind, adt_marker_select][which(adt_tmp >= valley)] = "pos"
    }
}
rownames(cell_x_adt_posNeg) = rownames(cell_x_adt)
saveRDS(list(cell_x_adt_posNeg, adt_feature), file = paste0(out_path, "/RDS/cell_x_adt_posNeg_cell_x_feature_", run_name, ".rds"))


## compare centers
pdf(paste0(out_path, "/Figures/stain_quality_compare_across_center.pdf"), width = 9, height = 13)
peak_sep_summary %>% dplyr::filter(adt_marker %in% peak2markerselect) %>% group_by(adt_marker, batch) %>% summarize(mean_sep_power = mean(sep_power)) %>% 
ggplot(aes(x = batch, y = mean_sep_power, fill = batch)) +
geom_boxplot(outlier.shape = NA) +
geom_jitter(width = 0.2) +
theme_bw(base_size = 25) +
scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(peak_sep_summary$batch)) + 1)) +
xlab("") +
rremove("legend.title") +
rremove("legend") +
ylab("Average Separation Power per Protein Marker") +
ggtitle("Protein Markers with Multiple Peaks")

peak_num_summary %>% dplyr::filter(adt_marker %in% peak2markerselect) %>% group_by(adt_marker, batch) %>% summarize(mean_peak_num = mean(peak_num)) %>% 
ggplot(aes(x = batch, y = mean_peak_num, fill = batch)) +
geom_boxplot(outlier.shape = NA) +
geom_jitter(width = 0.2) +
theme_bw(base_size = 25) +
scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(peak_sep_summary$batch)) + 1)) +
xlab("") +
rremove("legend.title") +
rremove("legend") +
ylab("Average Peak Number per Protein Marker") +
ggtitle("Protein Markers with Multiple Peaks")

positive_prop_summary %>% dplyr::filter(adt_marker %in% peak2markerselect) %>% group_by(adt_marker, batch) %>% summarize(mean_pos_prop = mean(positive_prop)) %>% 
ggplot(aes(x = batch, y = mean_pos_prop, fill = batch)) +
geom_boxplot(outlier.shape = NA) +
geom_jitter(width = 0.2) +
theme_bw(base_size = 25) +
scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(peak_sep_summary$batch)) + 1)) +
xlab("") +
rremove("legend.title") +
rremove("legend") +
ylab("Average Positive Proportion per Protein Marker") +
ggtitle("Protein Markers with Multiple Peaks")


left_join(
    peak_sep_summary %>% dplyr::filter(adt_marker %in% peak2markerselect) %>% group_by(adt_marker, batch)%>% summarize(mean_sep_power = mean(sep_power)),
    peak_num_summary %>% dplyr::filter(adt_marker %in% peak2markerselect) %>% group_by(adt_marker, batch) %>% summarize(mean_peak_num = mean(peak_num)), by = c("adt_marker", "batch")
) %>% ggplot(aes(x = mean_peak_num, y = mean_sep_power, color = batch)) +
geom_point() +
geom_smooth(method = "lm", se = FALSE) +
theme_bw(base_size = 25) +
scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(peak_sep_summary$batch)) + 1)) +
xlab("Average Peak Number per Protein Marker") +
ylab("Average Separation Power per Protein Marker")
dev.off()


pdf(paste0(out_path, "/Figures/stain_quality_acrossMarkers_across_center.pdf"), width = 25, height = 10)

adt_marker_sorted = peak_sep_summary %>% dplyr::filter(adt_marker %in% peak2marker) %>% group_by(adt_marker)%>% summarize(mean_sep_power = mean(sep_power)) %>% arrange(mean_sep_power) %$% adt_marker

peak_sep_summary %>% dplyr::filter(adt_marker %in% peak2marker) %>% group_by(adt_marker, batch)%>% summarize(mean_sep_power = mean(sep_power)) %>% 
ggplot(aes(x = factor(adt_marker, adt_marker_sorted), y = mean_sep_power)) +
geom_boxplot(outlier.shape = NA) +
geom_point(size = 2, aes(color = batch)) +
theme_bw(base_size = 25) +
xlab("") +
ylab("Average Separation Power") +
rotate_x_text(angle = 90) + 
rremove("legend.title") +
scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(peak_sep_summary$batch)) + 1))

adt_marker_sorted = peak_num_summary %>% dplyr::filter(adt_marker %in% peak2marker) %>% group_by(adt_marker)%>% summarize(mean_peak_num = mean(peak_num)) %>% arrange(mean_peak_num) %$% adt_marker

peak_num_summary %>% dplyr::filter(adt_marker %in% peak2marker) %>% group_by(adt_marker, batch)%>% summarize(mean_peak_num = mean(peak_num)) %>% 
ggplot(aes(x = factor(adt_marker, adt_marker_sorted), y = mean_peak_num)) +
geom_boxplot(outlier.shape = NA) +
geom_point(size = 2, aes(color = batch)) +
theme_bw(base_size = 25) +
xlab("") +
ylab("Average Peak Number") +
rotate_x_text(angle = 90) + 
rremove("legend.title") +
scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(peak_sep_summary$batch)) + 1))

adt_marker_sorted = positive_prop_summary %>% dplyr::filter(adt_marker %in% peak2marker) %>% group_by(adt_marker) %>% summarize(mean_pos_prop = mean(positive_prop)) %>% arrange(mean_pos_prop) %$% adt_marker

positive_prop_summary %>% dplyr::filter(adt_marker %in% peak2marker) %>% group_by(adt_marker, batch)%>% summarize(mean_pos_prop = mean(positive_prop)) %>% 
ggplot(aes(x = factor(adt_marker, adt_marker_sorted), y = mean_pos_prop)) +
geom_boxplot(outlier.shape = NA) +
geom_point(size = 2, aes(color = batch)) +
theme_bw(base_size = 25) +
xlab("") +
ylab("Average Positive Proportion") +
rotate_x_text(angle = 90) + 
rremove("legend.title") +
scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(peak_sep_summary$batch)) + 1))

dev.off()

## patient comparison
target_status = c("Healthy", "Asymptomatic", "Mild", "Moderate", "Severe", "Critical")

pdf(paste0(out_path, "/Figures/sample_number_distribution_across_status_center.pdf"), width = 8, height = 8)
adt_feature %>% data.frame %>% dplyr::select(Status_on_day_collection_summary, batch, sample) %>% unique %>% group_by(Status_on_day_collection_summary, batch) %>% summarize(freq = n()) %>% 
ggplot(aes(x = factor(Status_on_day_collection_summary, levels = c(target_status, "LPS_90mins", "LPS_10hours", "Non_covid")), y = freq, fill = batch)) +
geom_bar(stat = "identity", position = position_dodge()) +
theme_bw(base_size = 25) +
scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(4)) +
xlab("") +
ylab("Sample Number") +
rremove("legend.title") +
theme(legend.position = "top") +
rotate_x_text(angle = 45)
dev.off()

peak_sep_summary_info = left_join(peak_sep_summary, data.frame(adt_feature)[, c("Collection_Day", "Sex", "Status", "Swab_result", "Status_on_day_collection_summary", "Worst_Clinical_Status", "sample")] %>% unique, by = "sample")
peak_num_summary_info = left_join(peak_num_summary, data.frame(adt_feature)[, c("Collection_Day", "Sex", "Status", "Swab_result", "Status_on_day_collection_summary", "Worst_Clinical_Status", "sample")] %>% unique, by = "sample")
positive_prop_summary_info = left_join(positive_prop_summary, data.frame(adt_feature)[, c("Collection_Day", "Sex", "Status", "Swab_result", "Status_on_day_collection_summary", "Worst_Clinical_Status", "sample")] %>% unique, by = "sample")


pdf(paste0(out_path, "/Figures/stain_quality_acrossStatus.pdf"), width = 30, height = 10)
adt_marker_sorted = peak_sep_summary %>% dplyr::filter(adt_marker %in% peak2marker) %>% group_by(adt_marker)%>% summarize(mean_sep_power = mean(sep_power)) %>% arrange(mean_sep_power) %$% adt_marker

peak_sep_summary_info %>% dplyr::filter(Status_on_day_collection_summary %in% target_status) %>% dplyr::filter(adt_marker %in% peak2markerselect) %>% 
ggplot(aes(x = factor(adt_marker, levels = adt_marker_sorted), y = sep_power, fill = factor(Status_on_day_collection_summary, levels = target_status))) +
geom_boxplot(outlier.shape = NA) +
# geom_jitter(width = 0.2) +
theme_bw(base_size = 25) +
scale_fill_brewer(palette = "Set1") +
xlab("") +
rremove("legend.title") +
scale_y_continuous(trans = "log2") +
rremove("legend.title") +
theme(legend.position = "top") +
rotate_x_text(angle = 90) +
ylab("Stain Quality (Separation Power)")

adt_marker_sorted = peak_num_summary %>% dplyr::filter(adt_marker %in% peak2marker) %>% group_by(adt_marker)%>% summarize(mean_peak_num = mean(peak_num)) %>% arrange(mean_peak_num) %$% adt_marker

peak_num_summary_info %>% dplyr::filter(Status_on_day_collection_summary %in% target_status) %>% dplyr::filter(adt_marker %in% peak2markerselect) %>% 
ggplot(aes(x = factor(adt_marker, levels = adt_marker_sorted), y = peak_num, fill = factor(Status_on_day_collection_summary, levels = target_status))) +
geom_boxplot(outlier.shape = NA) +
# geom_jitter(width = 0.2) +
theme_bw(base_size = 25) +
scale_fill_brewer(palette = "Set1") +
xlab("") +
rremove("legend.title") +
scale_y_continuous(trans = "log2") +
rremove("legend.title") +
theme(legend.position = "top") +
rotate_x_text(angle = 90) +
ylab("Peak Number")

adt_marker_sorted = positive_prop_summary %>% dplyr::filter(adt_marker %in% peak2marker) %>% group_by(adt_marker)%>% summarize(mean_pos_prop = mean(positive_prop)) %>% arrange(mean_pos_prop) %$% adt_marker

positive_prop_summary_info %>% dplyr::filter(Status_on_day_collection_summary %in% target_status) %>% dplyr::filter(adt_marker %in% peak2markerselect) %>% 
ggplot(aes(x = factor(adt_marker, levels = adt_marker_sorted), y = positive_prop, fill = factor(Status_on_day_collection_summary, levels = target_status))) +
geom_boxplot(outlier.shape = NA) +
# geom_jitter(width = 0.2) +
theme_bw(base_size = 25) +
scale_fill_brewer(palette = "Set1") +
xlab("") +
rremove("legend.title") +
scale_y_continuous(trans = "log2") +
rremove("legend.title") +
theme(legend.position = "top") +
rotate_x_text(angle = 90) +
ylab("Positive Cells Proportion")
dev.off()

## heatmap of positive proportion
summary_data <- positive_prop_summary_info %>% dplyr::filter(Status_on_day_collection_summary %in% target_status) %>%
  group_by(adt_marker, Status_on_day_collection_summary) %>%
  summarise(mean_positive_prop = mean(positive_prop, na.rm = TRUE)) %>%
  ungroup() %>% dplyr::filter(adt_marker %in% peak2markerselect)

heatmap_data <- summary_data %>%
  pivot_wider(names_from = adt_marker, values_from = mean_positive_prop)
heatmap_data = heatmap_data[match(target_status, heatmap_data$Status_on_day_collection_summary), ]
# Prepare data for pheatmap (remove the first column which contains row names)
row_names <- heatmap_data$Status_on_day_collection_summary
heatmap_matrix <- as.matrix(heatmap_data[, -1])
rownames(heatmap_matrix) = row_names
# Generate the heatmap
pheatmap(heatmap_matrix,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "none",
         filename = paste0(out_path, "/Figures/mean_positive_prop_acrossStatus_heatmap_noscaled.pdf"),
         width = 16, 
         height = 3.5,
         border_color = NA  # Hide cell borders for a cleaner look
)
pheatmap(heatmap_matrix,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "column",
        #  annotation_row = annotation_row,
         filename = paste0(out_path, "/Figures/mean_positive_prop_acrossStatus_heatmap_scaledcolumn.pdf"),
         width = 16, 
         height = 3.5,
         border_color = NA  # Hide cell borders for a cleaner look
)

## heatmap removing the NCL data
summary_data <- positive_prop_summary_info %>% dplyr::filter(batch != "Ncl") %>%
  dplyr::filter(Status_on_day_collection_summary %in% target_status) %>%
  group_by(adt_marker, Status_on_day_collection_summary) %>%
  summarise(mean_positive_prop = mean(positive_prop, na.rm = TRUE)) %>%
  ungroup() %>% dplyr::filter(adt_marker %in% peak2markerselect)

heatmap_data <- summary_data %>%
  pivot_wider(names_from = adt_marker, values_from = mean_positive_prop)
heatmap_data = heatmap_data[match(target_status, heatmap_data$Status_on_day_collection_summary), ]
# Prepare data for pheatmap (remove the first column which contains row names)
row_names <- heatmap_data$Status_on_day_collection_summary
heatmap_matrix <- as.matrix(heatmap_data[, -1])
rownames(heatmap_matrix) = row_names
# Generate the heatmap
pheatmap(heatmap_matrix,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "none",
        #  annotation_row = annotation_row,
         filename = paste0(out_path, "/Figures/mean_positive_prop_acrossStatus_heatmap_noscaled_removeNCL.pdf"),
         width = 16, 
         height = 3.5,
         border_color = NA  # Hide cell borders for a cleaner look
)
pheatmap(heatmap_matrix,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "column",
        #  annotation_row = annotation_row,
         filename = paste0(out_path, "/Figures/mean_positive_prop_acrossStatus_heatmap_scaledcolumn_removeNCL.pdf"),
         width = 16, 
         height = 3.5,
         border_color = NA  # Hide cell borders for a cleaner look
)

## using all marker
summary_data <- positive_prop_summary_info %>% 
  dplyr::filter(Status_on_day_collection_summary %in% target_status) %>%
  group_by(adt_marker, Status_on_day_collection_summary) %>%
  summarise(mean_positive_prop = mean(positive_prop, na.rm = TRUE)) %>%
  ungroup() #%>% dplyr::filter(adt_marker %in% peak2markerselect)

heatmap_data <- summary_data %>%
  pivot_wider(names_from = adt_marker, values_from = mean_positive_prop)
heatmap_data = heatmap_data[match(target_status, heatmap_data$Status_on_day_collection_summary), ]
# Prepare data for pheatmap (remove the first column which contains row names)
row_names <- heatmap_data$Status_on_day_collection_summary
heatmap_matrix <- as.matrix(heatmap_data[, -1])
rownames(heatmap_matrix) = row_names

annotation_col = data.frame(peakType = c(rep("TwoPeaks", length(peak2marker)), rep("OnePeak", length(marker_list) - length(peak2marker))))
rownames(annotation_col) = c(peak2marker, marker_list[which(!(marker_list %in% peak2marker))])
annotation_color = list(peakType = c(OnePeak = "#3b5998", TwoPeaks = "#e8702a"))

# Generate the heatmap
pheatmap(heatmap_matrix,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "none",
         annotation_col = annotation_col,
         annotation_colors = annotation_color,
        #  annotation_row = annotation_row,
         filename = paste0(out_path, "/Figures/mean_positive_prop_acrossStatus_heatmap_noscaled_allmarker.pdf"),
         width = 23, 
         height = 3.5,
         border_color = NA  # Hide cell borders for a cleaner look
)
pheatmap(heatmap_matrix,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "column",
         annotation_col = annotation_col,
         annotation_colors = annotation_color,
        #  annotation_row = annotation_row,
         filename = paste0(out_path, "/Figures/mean_positive_prop_acrossStatus_heatmap_scaledcolumn_allmarker.pdf"),
         width = 23, 
         height = 3.5,
         border_color = NA  # Hide cell borders for a cleaner look
)


summary_data <- positive_prop_summary_info %>% dplyr::filter(batch != "Ncl") %>%
  dplyr::filter(Status_on_day_collection_summary %in% target_status) %>%
  group_by(adt_marker, Status_on_day_collection_summary) %>%
  summarise(mean_positive_prop = mean(positive_prop, na.rm = TRUE)) %>%
  ungroup() #%>% dplyr::filter(adt_marker %in% peak2markerselect)

heatmap_data <- summary_data %>%
  pivot_wider(names_from = adt_marker, values_from = mean_positive_prop)
heatmap_data = heatmap_data[match(target_status, heatmap_data$Status_on_day_collection_summary), ]
# Prepare data for pheatmap (remove the first column which contains row names)
row_names <- heatmap_data$Status_on_day_collection_summary
heatmap_matrix <- as.matrix(heatmap_data[, -1])
rownames(heatmap_matrix) = row_names

annotation_col = data.frame(peakType = c(rep("TwoPeaks", length(peak2marker)), rep("OnePeak", length(marker_list) - length(peak2marker))))
rownames(annotation_col) = c(peak2marker, marker_list[which(!(marker_list %in% peak2marker))])
annotation_color = list(peakType = c(OnePeak = "#3b5998", TwoPeaks = "#e8702a"))

# Generate the heatmap
pheatmap(heatmap_matrix,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "none",
         annotation_col = annotation_col,
         annotation_colors = annotation_color,
        #  annotation_row = annotation_row,
         filename = paste0(out_path, "/Figures/mean_positive_prop_acrossStatus_heatmap_noscaled_removeNCL_allmarker.pdf"),
         width = 23, 
         height = 3.5,
         border_color = NA  # Hide cell borders for a cleaner look
)
pheatmap(heatmap_matrix,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "column",
         annotation_col = annotation_col,
         annotation_colors = annotation_color,
        #  annotation_row = annotation_row,
         filename = paste0(out_path, "/Figures/mean_positive_prop_acrossStatus_heatmap_scaledcolumn_removeNCL_allmarker.pdf"),
         width = 23, 
         height = 3.5,
         border_color = NA  # Hide cell borders for a cleaner look
)

## ==========
## Mono DE
## ==========
get_positive_cell_prop_cT = function(cell_x_feature, cell_x_adt, adt_marker_select, each_sample, valley_info, cell_type_select){
    ind = which(cell_x_feature$sample == each_sample & cell_x_feature$full_clustering == cell_type_select)
    adt_tmp = cell_x_adt[ind, adt_marker_select]
    valley = valley_info[each_sample, 1]
    total_cell_num = length(which(cell_x_feature$sample == each_sample))
    
    pos_num = sum(adt_tmp >= valley)
    cT_cell_num = length(adt_tmp)
    prop = pos_num/total_cell_num * 100

  
    return(list(pos_num, cT_cell_num, prop, total_cell_num))
}

get_positive_cell_prop_cT = function(cell_x_feature, cell_x_adt, adt_marker_select, each_sample, valley_info, cell_type_select){
    ind = which(cell_x_feature$sample == each_sample & cell_x_feature$full_clustering == cell_type_select)
    adt_tmp = cell_x_adt[ind, adt_marker_select]
    valley = valley_info[each_sample, 1]
    # total_cell_num = length(which(cell_x_feature$sample == each_sample))
    total_cell_num = length(ind)

    pos_num = sum(adt_tmp >= valley)
    cT_cell_num = length(adt_tmp)
    prop = pos_num/total_cell_num * 100

  
    return(list(pos_num, cT_cell_num, prop, total_cell_num))
}

## use valley and peak
adt_feature$sample = adt_feature$site_sample

positive_prop_perCT_summary = c()
for(adt_marker_select in colnames(cell_x_adt)){ ## colnames(cell_x_adt)
    print(adt_marker_select)
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
    valley_info = peak_valley_density$valley_landmark_list
    for(cell_type_select in c("CD83_CD14_mono", "CD14_mono", "CD16_mono", "C1_CD16_mono")){
        for(each_sample in unique(adt_feature$sample)){
            batch_info = each_sample
            tmp = get_positive_cell_prop_cT(cell_x_feature = adt_feature, cell_x_adt, adt_marker_select, each_sample, valley_info, cell_type_select)

            positive_prop_perCT_summary = data.frame(
                pos_num = tmp[[1]],
                cT_cell_num = tmp[[2]],
                positive_prop = tmp[[3]],
                total_cell_num = tmp[[4]],
                batch = batch_info, 
                sample = each_sample, 
                adt_marker = adt_marker_select,
                cell_type_select = cell_type_select) %>% 
                rbind(positive_prop_perCT_summary, .)
        }
    }
}


method = "ADTnorm_sample_manualtuning"
saveRDS(positive_prop_perCT_summary, file = paste0(out_path, "/RDS/peak_stain_quality_", run_name, "_", method, "_positive_prop_perCT_summary_update.rds"))

positive_prop_perCT_summary = readRDS(file = paste0(out_path, "/RDS/peak_stain_quality_", run_name, "_", method, "_positive_prop_perCT_summary.rds"))

positive_prop_perCT_summary$batch = positive_prop_perCT_summary$batch %>% gsub("_.*", "", .)
positive_prop_perCT_summary$batch = factor(positive_prop_perCT_summary$batch, levels = file_list)
positive_prop_perCT_summary$adt_marker = factor(positive_prop_perCT_summary$adt_marker, levels = marker_list)

positive_prop_perCT_summary$cT_cell_prop = positive_prop_perCT_summary$cT_cell_num/positive_prop_perCT_summary$total_cell_num * 100

positive_prop_perCT_summary_info = left_join(positive_prop_perCT_summary, data.frame(adt_feature)[, c("Collection_Day", "Sex", "Status", "Swab_result", "Status_on_day_collection_summary", "Worst_Clinical_Status", "sample")] %>% unique, by = "sample")

target_status = c("Healthy", "Asymptomatic", "Mild", "Moderate", "Severe", "Critical")
filter_prop_perCT = positive_prop_perCT_summary_info %>% 
dplyr::filter(Status_on_day_collection_summary %in% target_status) %>% 
dplyr::filter(!(adt_marker %in% c("IgG1kappa", "IgG2akappa", "IgG2bkappa", "IgG2bkappaIsotypeCtrl", "FcεRIalpha", "IgGFc", "Iglightchainkappa", "Iglightchainlambda"))) %>% 
dplyr::filter(cT_cell_num > 70) %>%
dplyr::filter(!(batch == "Cambridge" & adt_marker == "CD31")) 

filter_prop_perCT$Status_on_day_collection_summary = factor(filter_prop_perCT$Status_on_day_collection_summary, levels = target_status)

wtest_perCT_summary = c()
for(cell_type_each in c("CD83_CD14_mono", "CD14_mono", "CD16_mono", "C1_CD16_mono")){
    print(cell_type_each)
    ## to select adt_marker
    if(cell_type_each %in% c("CD83_CD14_mono", "CD14_mono", "C1_CD16_mono")){
    select_adtmarker = filter_prop_perCT %>% 
    dplyr::filter(cell_type_select == cell_type_each) %>%
    dplyr::filter(Status_on_day_collection_summary %in% target_status) %>% 
    dplyr::filter(!(adt_marker %in% c("IgG1kappa", "IgG2akappa", "IgG2bkappa", "IgG2bkappaIsotypeCtrl", "FcεRIalpha", "IgGFc", "Iglightchainkappa", "Iglightchainlambda"))) %>% 
    dplyr::filter(positive_prop > 10) %$% adt_marker %>% unique()
    }else{
    select_adtmarker = filter_prop_perCT %>% 
    dplyr::filter(cell_type_select == cell_type_each) %>%
    dplyr::filter(Status_on_day_collection_summary %in% target_status) %>% 
    dplyr::filter(!(adt_marker %in% c("IgG1kappa", "IgG2akappa", "IgG2bkappa", "IgG2bkappaIsotypeCtrl", "FcεRIalpha", "IgGFc", "Iglightchainkappa", "Iglightchainlambda"))) %>% 
    dplyr::filter(positive_prop > 5) %$% adt_marker %>% unique()
    }


    for(marker_each in select_adtmarker){
        prop_tmp = filter_prop_perCT %>% dplyr::filter(adt_marker == marker_each, cell_type_select == cell_type_each)
        hd_ind = which(prop_tmp$Status_on_day_collection_summary %in% c("Healthy"))
        pt_ind = which(prop_tmp$Status_on_day_collection_summary %in% c("Asymptomatic", "Mild", "Moderate", "Severe", "Critical"))
        if(length(hd_ind) > 2 && length(pt_ind) > 2){
            wt_p = wilcox.test(prop_tmp$positive_prop[hd_ind], prop_tmp$positive_prop[pt_ind], alternative = "two.sided")$p.value
            tt_p = t.test(prop_tmp$positive_prop[hd_ind], prop_tmp$positive_prop[pt_ind], alternative = "two.sided")$p.value
            test_logfc = log2(mean(prop_tmp$positive_prop[pt_ind], na.rm = TRUE)/mean(prop_tmp$positive_prop[hd_ind], na.rm = TRUE))
            if(mean(prop_tmp$positive_prop[hd_ind], na.rm = TRUE) == 0){
                test_logfc = 0
            }
            wtest_perCT_summary = data.frame(adt_marker = marker_each, wt_p = wt_p, tt_p = tt_p, logfc = test_logfc, mean_hd = mean(prop_tmp$positive_prop[hd_ind], na.rm = TRUE), mean_pt = mean(prop_tmp$positive_prop[pt_ind], na.rm = TRUE), cell_type_select = cell_type_each) %>% rbind(wtest_perCT_summary, .)
        }

    }
}

wtest_perCT_summary$wt_color = "nonsign"
wtest_perCT_summary$wt_color[which(wtest_perCT_summary$wt_p < 0.01 & wtest_perCT_summary$logfc > log2(1.5))] = "up"
wtest_perCT_summary$wt_color[which(wtest_perCT_summary$wt_p < 0.01 & wtest_perCT_summary$logfc < -log2(1.5))] = "down"
wtest_perCT_summary$tt_color = "nonsign"
wtest_perCT_summary$tt_color[which(wtest_perCT_summary$tt_p < 0.01 & wtest_perCT_summary$logfc > log2(1.5))] = "up"
wtest_perCT_summary$tt_color[which(wtest_perCT_summary$tt_p < 0.01 & wtest_perCT_summary$logfc < -log2(1.5))] = "down"

p1 = wtest_perCT_summary %>% dplyr::filter(cell_type_select == "CD14_mono") %>% ggplot(aes(x = logfc, y = -log10(tt_p), label = adt_marker, color = tt_color)) +
geom_point() +
geom_label_repel(box.padding = 0.1, max.overlaps = 15) +
theme_bw(base_size = 20) +
xlab("log2(Fold Change)") +
ylab("-log10(P Value)") +
ggtitle("CD14 Monocytes") +
geom_vline(xintercept = -log2(1.5), linetype = "dashed", color = "grey50", size = 1) +
geom_vline(xintercept = log2(1.5), linetype = "dashed", color = "grey50", size = 1) +
geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey50", size = 1) +
scale_color_manual(values = c("grey", "#A93226", "#2471A3"), breaks = c("nonsign", "up", "down")) +
rremove("legend")


p2 = wtest_perCT_summary %>% dplyr::filter(cell_type_select == "CD16_mono") %>% ggplot(aes(x = logfc, y = -log10(tt_p), label = adt_marker, color = tt_color)) +
geom_point() +
geom_label_repel(box.padding = 0.1, max.overlaps = 15) +
theme_bw(base_size = 20) +
xlab("log2(Fold Change)") +
ylab("-log10(P Value)") +
ggtitle("CD16 Monocytes") +
geom_vline(xintercept = -log2(1.5), linetype = "dashed", color = "grey50", size = 1) +
geom_vline(xintercept = log2(1.5), linetype = "dashed", color = "grey50", size = 1) +
geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey50", size = 1) +
scale_color_manual(values = c("grey", "#A93226", "#2471A3"), breaks = c("nonsign", "up", "down")) +
rremove("legend")


p3 = wtest_perCT_summary %>% dplyr::filter(cell_type_select == "CD83_CD14_mono") %>% ggplot(aes(x = logfc, y = -log10(tt_p), label = adt_marker, color = tt_color)) +
geom_point() +
geom_label_repel(box.padding = 0.1, max.overlaps = 15) +
theme_bw(base_size = 20) +
xlab("log2(Fold Change)") +
ylab("-log10(P Value)") +
ggtitle("CD83+ CD14 Monocytes") +
geom_vline(xintercept = -log2(1.5), linetype = "dashed", color = "grey50", size = 1) +
geom_vline(xintercept = log2(1.5), linetype = "dashed", color = "grey50", size = 1) +
geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey50", size = 1) +
scale_color_manual(values = c("grey", "#A93226", "#2471A3"), breaks = c("nonsign", "up", "down")) +
rremove("legend")

pdf(paste0(out_path, "/Figures/Covid19_mono_posprop_volcano.pdf"), width = 15, height = 8)
ggarrange(p1, p2, p3, ncol = 3, nrow = 1)
dev.off()

## heatmap of mean value of adt marker and positive proportion
## compare between healthy and covid patients

target_markers = c("CD38", "CD64", "CD169")
target_marker_summary = c()
for(cell_type_each in c("CD83_CD14_mono", "CD14_mono", "CD16_mono")){
    
    cell_hd = which(adt_feature$full_clustering == cell_type_each & adt_feature$Status_on_day_collection_summary == "Healthy")
    cell_pt = which(adt_feature$full_clustering == cell_type_each & adt_feature$Status_on_day_collection_summary %in% c("Asymptomatic", "Mild", "Moderate", "Severe", "Critical"))
    filter_prop_perCT_ind = which(filter_prop_perCT$cell_type_select == cell_type_each)

    target_marker_summary = cbind(res_norm[cell_hd, target_markers], adt_feature[cell_hd, "sample", drop = FALSE]) %>% 
    data.frame %>% dplyr::group_by(sample) %>% 
    summarize(CD38 = mean(CD38), CD64 = mean(CD64), CD169 = mean(CD169)) %>% 
    tidyr::pivot_longer(cols = c(CD38, CD64, CD169), names_to = "adt_marker", values_to = "mean_value") %>% 
    mutate(type = "Healthy", cell_type = cell_type_each) %>%     
    left_join(., filter_prop_perCT[filter_prop_perCT_ind, c("sample", "adt_marker", "positive_prop")], by = c("sample", "adt_marker")) %>% 
    dplyr::filter(!(is.na(positive_prop))) %>%
    rbind(target_marker_summary, .)

    target_marker_summary = cbind(res_norm[cell_pt, target_markers], adt_feature[cell_pt,  "sample", drop = FALSE]) %>%
    data.frame %>% dplyr::group_by(sample) %>%
    summarize(CD38 = mean(CD38), CD64 = mean(CD64), CD169 = mean(CD169)) %>%
    tidyr::pivot_longer(cols = c(CD38, CD64, CD169), names_to = "adt_marker", values_to = "mean_value") %>%
    mutate(type = "COVID-19", cell_type = cell_type_each) %>% 
    left_join(., filter_prop_perCT[filter_prop_perCT_ind, c("sample", "adt_marker", "positive_prop")], by = c("sample", "adt_marker")) %>% 
    dplyr::filter(!(is.na(positive_prop))) %>%
    rbind(target_marker_summary, .)
    
}
target_marker_summary$adt_marker = factor(target_marker_summary$adt_marker, levels = target_markers)

b1 = target_marker_summary %>% dplyr::filter(!is.na(positive_prop)) %>% 
ggplot(aes(x = cell_type, y = mean_value, fill = type)) +
geom_boxplot() +
facet_grid(~adt_marker) +
theme_bw(base_size = 20) +
xlab("") +
ylab("Normalized ADT Counts") +
scale_fill_brewer(palette = "Set1") +
rremove("legend.title") +
theme(legend.position = "top") +
rotate_x_text(angle = 20)

b2 = target_marker_summary %>% dplyr::filter(!is.na(positive_prop)) %>% 
ggplot(aes(x = cell_type, y = positive_prop, fill = type)) +
geom_boxplot() +
facet_grid(~adt_marker) +
theme_bw(base_size = 20) +
xlab("") +
ylab("ADT Marker Positive Proportion") +
scale_fill_brewer(palette = "Set1") +
rremove("legend.title") +
theme(legend.position = "top") +
rotate_x_text(angle = 20)

pdf(paste0(out_path, "/Figures/Covid19_mono_posprop_boxplot_update.pdf"), width = 11, height = 16)
ggarrange(b1, b2, ncol = 1, nrow = 2, common.legend = TRUE)
dev.off()

pdf(paste0(out_path, "/Figures/Covid19_mono_posprop_bubble_update.pdf"), width = 15, height = 7)
target_marker_summary %>% group_by(type, cell_type, adt_marker) %>% summarize(NormalizedCount = mean(mean_value), PositiveProportion = mean(positive_prop)) %>% ungroup %>% group_by(adt_marker) %>%
  mutate(z_score = scale(NormalizedCount, center = TRUE, scale = TRUE)) %>%
  ungroup() %>%
 ggplot(aes(x = type, y = adt_marker, color = z_score, size = PositiveProportion)) +
geom_point() +
theme_bw(base_size = 25) +
scale_color_viridis(option = "magma", begin = 0.9, end = 0) +
facet_grid(~cell_type) +
xlab("") +
ylab("") +
scale_size_continuous(range = c(1, 40))
dev.off()

