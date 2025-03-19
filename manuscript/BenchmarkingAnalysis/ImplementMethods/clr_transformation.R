## CLR transformation
require(dplyr)
require(Seurat)

clr_transform <- function(cell_x_adt = NULL, cell_x_feature = NULL, parameter_list = NULL){

    ## output
    out <- cell_x_adt %>%
        NormalizeData(normalization.method = "CLR", scale.factor = 1000000, margin = 1)
    return(out)
}