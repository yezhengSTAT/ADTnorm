## Harmony remove batch effect
require(harmony)
require(dplyr)

harmony_transform <- function(cell_x_adt = NULL, cell_x_feature = NULL, parameter_list = NULL) {
    out <- HarmonyMatrix(
        data_mat = cell_x_adt,
        meta_data = cell_x_feature$batch,
        do_pca = FALSE
    )
    return(out)
}
