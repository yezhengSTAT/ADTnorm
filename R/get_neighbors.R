#' Identify the valley outliers and impute by valley by closet neighbor samples on the graph.
#'
#' This function identify the valley that tend to be outliers compared to other valley locations and try to find the closest samples that have similar density distribution to imput the valley. If no neighbor sample is detected, the valley will remain as original.
#' @param target_sample The target sample whose valley need to be imputed. Find the neighbor samples whose density distribution is close to the target sample.
#' @param adt_marker_select The marker whose valley need to be imputed. Fine the neighbor samples whose density distribution is close to the target sample of the same ADT marker.
#' @param cell_x_adt Matrix of ADT raw counts in cells (rows) by ADT markers (columns) format.
#' @param cell_x_feature Matrix of cells (rows) by cell features (columns) such as cell type, sample, and batch related information.
#' @param nearest_neighbor_n Number of top nearest neighbor samples to detect.
#' @param nearest_neighbor_threshold Threshold to call neighbor samples.
#' @export
#' @examples
#' \dontrun{
#' get_neighbors(target_sample, adt_marker_select, cell_x_adt, cell_x_feature)
#' }
# require(EMDomics)
# require(dplyr)
get_neighbors <- function(target_sample, adt_marker_select, cell_x_adt, cell_x_feature, nearest_neighbor_n = 3, nearest_neighbor_threshold = NULL){

    knn_res <- c()
    target_cell_ind <- which(cell_x_feature$sample == target_sample)
    sample_list <- setdiff(cell_x_feature$sample %>% unique, target_sample)
    for(sample in sample_list){
         cell_ind <- which(cell_x_feature$sample == sample)
         exp_data <- c(cell_x_adt[target_cell_ind, adt_marker_select], cell_x_adt[cell_ind, adt_marker_select])
         names(exp_data) <- rownames(cell_x_adt)[c(target_cell_ind, cell_ind)]
         labels <- c(rep("target", length(target_cell_ind)), rep("sample", length(cell_ind)))
         names(labels) <- rownames(cell_x_adt)[c(target_cell_ind, cell_ind)]

         knn_res <- c(knn_res, EMDomics::calculate_emd_gene(exp_data, labels, names(exp_data)))
    }
    names(knn_res) <- sample_list
    if(is.null(nearest_neighbor_threshold)){
        return(knn_res %>% sort %>% utils::head(nearest_neighbor_n) %>% names)
    }else{
        return(knn_res[knn_res <= nearest_neighbor_threshold] %>% sort %>% utils::head(nearest_neighbor_n) %>% names)
    }


}
