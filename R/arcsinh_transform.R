#' Arcsine transformation.
#'
#' This function arcsine transform the input cell_x_adt matrix with co-factor 5. The definition of this function is x_new <- asinh(a + b * x) + c)
#' @param cell_x_adt Matrix where rows are cells and columns are ADT markers.
#' @param parameter_list Parameter list for a: positive double that corresponds to a shift about 0; b: positive double that corresponds to a scale factor; c: positive double. By default a = 1, b = 1/5 and c = 0.
#' @export
#' @examples
#' \dontrun{
#' arcsinh_transform(cell_x_adt)
#' }

arcsinh_transform = function(cell_x_adt = NULL, parameter_list = NULL){
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
    asinhTrans = flowCore::arcsinhTransform(transformationId = "ln-transformation", a = a, b = b, c = c)

    ## output
    out = asinhTrans(cell_x_adt)
    return(out)
}
