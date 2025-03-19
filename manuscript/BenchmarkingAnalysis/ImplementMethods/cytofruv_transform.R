## CytofRUV normalization
require(writexl)
require(CATALYST)

cytofruv_transform <- function(cell_x_adt = NULL, cell_x_feature = NULL, parameter_list = NULL) {

    ## read parameters
    tmp_path <- parameter_list$tmp_path
    fcs_path <- parameter_list$fcs_path
    condition <- parameter_list$condition
    patient_id <- parameter_list$patient_id
    batch <- parameter_list$batch
    clusters_nb <- parameter_list$clusters_nb
    start_adt_method <- parameter_list$start_adt_method

    ## change column name to ensure all end with number
    ori_col_name <- colnames(cell_x_adt)
    colnames(cell_x_adt) <- paste0(colnames(cell_x_adt), "0")

    for (each_sample in cell_x_feature$sample %>% levels()) {
        if (!dir.exists(paste0(fcs_path, "/", start_adt_method, "_for_CytofRUV"))) {
            dir.create(paste0(fcs_path, "/", start_adt_method, "_for_CytofRUV"))
        }
        fcs_file_name <- paste0(fcs_path, "/", start_adt_method, "_for_CytofRUV/", each_sample, ".fcs")
        if (!file.exists(fcs_file_name)) {
            sample_ind <- which(cell_x_feature$sample == each_sample)
            fcs_count <- cell_x_adt[sample_ind, ] %>% as.matrix()
            fcs <- flowFrame(fcs_count)
            write.FCS(fcs, filename = fcs_file_name)
        }
    }
    md <- data.frame(
        file_name = paste0(levels(cell_x_feature$sample), ".fcs"),
        sample_id = levels(cell_x_feature$sample),
        condition = condition,
        patient_id = patient_id,
        batch = batch
    )

    panel <- data.frame(
        fcs_colname = colnames(cell_x_adt), 
        antigen = colnames(cell_x_adt), 
        marker_class = "type"
    )
    write_xlsx(x = md, path = paste0(fcs_path, "/", start_adt_method, "_for_CytofRUV/", "/Metadata.xlsx"))
    write_xlsx(x = panel, path = paste0(fcs_path, "/", start_adt_method, "_for_CytofRUV/", "/Panel.xlsx"))
    
    cytof_data <- load_data(paste0(fcs_path, "/", start_adt_method, "_for_CytofRUV/"), metadata_filename = "Metadata.xlsx", panel_filename = "Panel.xlsx")
    cytof_data$daf <- cluster_data(
        cytof_data$daf, 
        seed = 12345, 
        markers_to_use = cytof_data$daf %>% rownames(), 
        clusters_nb = clusters_nb)
    cytof_data$lineage_markers <- cytof_data$daf %>% rownames()

    dir_name_norm_data <- "CytofRUV"
    raw_data <- data.frame(
        sample = cytof_data$daf$sample_id,
        cluster = cluster_ids(cytof_data$daf), t(SummarizedExperiment::assay(cytof_data$daf, "exprs"))
    )
    rep_samples <- levels(cell_x_feature$sample) %>% list() # list(file_list)
    cluster_list_rep_samples <- list(seq(1, clusters_nb))
    k_value <- 5
    seed <- 1234

    normalise_data(
        data = cytof_data,
        raw_data = raw_data,
        rep_samples = rep_samples,
        norm_clusters = cluster_list_rep_samples,
        k = k_value,
        num_clusters = clusters_nb,
        wd_data = fcs_path,
        dir_norm_data = dir_name_norm_data
    )
    out <- c()
    for (each_sample in levels(cell_x_feature$sample)) {
        tmpNorm <- read.FCS(filename = paste0(fcs_path, dir_name_norm_data, "/Norm_RUVIII_k", k_value, "_", each_sample, ".fcs"))

        out <- exprs(tmpNorm) %>%
            round(3) %>%
            data.frame() %>%
            rbind(out, .)
    }
    colnames(out) <- ori_col_name
    return(out)


}
