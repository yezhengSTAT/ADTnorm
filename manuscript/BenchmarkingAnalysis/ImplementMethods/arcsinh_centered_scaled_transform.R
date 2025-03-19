## Arcsinh centered scaled transformation
require(dplyr)
require(Seurat)

arcsinh_centered_scaled_transform <- function(cell_x_adt = NULL, cell_x_feature = NULL, parameter_list = NULL){
    ## parameters
    a = 1
    b = 1/5
    c = 0
    if(!is.null(parameter_list)){
        if("a" %in% names(parameter_list)){
            a = parameter_list[["a"]]
        }
        if("b" %in% names(parameter_list)){
            b = parameter_list[["b"]]
        }
        if("c" %in% names(parameter_list)){
            c = parameter_list[["c"]]
        }
    }
    ## transformation
    asinhTrans <- arcsinhTransform(transformationId = "ln-transformation", a = a, b = b, c = c)

    ## output
    out <- ((cell_x_adt %>%
        asinhTrans() %>%
        NormalizeData(normalization.method = "CLR", scale.factor = 1000000, margin = 1) %>%
        exp()) - 1) * 3
    return(out)
}