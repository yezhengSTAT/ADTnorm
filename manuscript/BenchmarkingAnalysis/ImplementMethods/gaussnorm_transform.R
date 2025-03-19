## gaussNorm transformation
require(flowStats)

gaussnorm_transform <- function(cell_x_adt = NULL, cell_x_feature = NULL, parameter_list = NULL) {

    ## read in parameters
    parameter_list_name <- names(parameter_list)
    if ("fcs_path" %in% parameter_list_name) {
        fcs_path <- parameter_list$fcs_path
    } else {
        return("fcs_path is not provided in the parameter_list!")
    }

    if ("start_adt_method" %in% parameter_list_name) {
        start_adt_method <- parameter_list$start_adt_method
    } else {
        return("start_adt_method is not provided in the paramter_list!")
    }

    ## write out fcs files
    for (each_sample in cell_x_feature$sample %>% levels()) {
        if (!dir.exists(paste0(fcs_path, "/", start_adt_method))) {
            dir.create(paste0(fcs_path, "/", start_adt_method))
        }
        fcs_file_name <- paste0(fcs_path, "/", start_adt_method, "/", each_sample, ".fcs")
        if (!file.exists(fcs_file_name)) {
            sample_ind <- which(cell_x_feature$sample == each_sample)
            fcs_count <- cell_x_adt[sample_ind, ] %>% as.matrix()
            fcs <- flowFrame(fcs_count)
            write.FCS(fcs, filename = fcs_file_name)
        }
    }

    ## gaussnorm
    file_fcs <- read.ncdfFlowSet(paste0(fcs_path, "/", start_adt_method, "/", levels(cell_x_feature$sample), ".fcs"))
    adt_gaussnorm <- gaussNorm(file_fcs, colnames(cell_x_adt))$flowset

    out <- c()
    for (i in 1:length(adt_gaussnorm)) {
        out <- adt_gaussnorm[[i]] %>%
            exprs() %>%
            rbind(out, .)
    }

    return(out)
}
