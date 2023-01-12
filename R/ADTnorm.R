#' ADTnorm normalization to remove the technical variations across samples for each ADT marker.
#'
#' This function removes the technical variations such as batch effect, sequencing depth biases, antibody selection difference and antibody concentration differences, etc. The normalized samples are ready for integration across studies.
#' @param cell_x_adt Matrix of ADT raw counts in cells (rows) by ADT markers (columns) format.
#' @param cell_x_feature Matrix of cells (rows) by cell features (columns) such as cell type, sample, and batch-related information.
#' @param save_outpath The path to save the results.
#' @param study_name Name of this run.
#' @param marker_to_process Markers to normalize. Leaving empty to process all the ADT markers in cell_x_adt matrix.
#' @param exclude_zeroes Whether to allow events with zero counts on a particular marker to be altered by landmark registration. Recommend TRUE if zeroes in the data represent dropout (likely for large ADT panels, big datasets, or undersequenced data)
#' @param bimodal_marker Specify ADT markers that are likely to have two peaks based on researchers' prior knowledge or preliminary observation of particular data to be processed. Leaving it as default, ADTnorm will try to find the bimodal peak in all markers that are not listed in `trimodal_marker.`
#' @param trimodal_marker Index of the ADT markers that tend to have three peaks based on researchers' prior knowledge (e.g., CD4) or preliminary observation on particular data to be processed.
#' @param positive_peak A list variable containing a vector of ADT marker(s) and a corresponding vector of sample name(s) in matching order to specify that the uni-peak detected should be aligned to positive peaks. For example, for samples that only contain T cells. The only CD3 peak should be aligned to the positive peaks of other samples.
#' @param bw_smallest_bi The smallest bandwidth parameter value for bi-modal peaks. Recommend 1.1.
#' @param bw_smallest_tri The smallest bandwidth parameter value for tri-modal peaks. Recommend the same value for CD4, such as 0.5.
#' @param bw_smallest_adjustments A named list of floats, with names matching marker names, specifying the smallest bandwidth parameter value. Recommend 0.5 or 0.8.
#' @param input_raw_counts Whether to (TRUE), by default, expect raw counts to be inputted and arcsin transformation to be performed, or (FALSE) bypass the arcsin transformation. Recommend trying log-normalization.
#' @param quantile_clip Implement an upper quantile clipping (useful to avoid warping function errors caused by few events of very high expression). If greater than 0, clips expression of events in that upper quantile.
#' @param peak_type The type of peak to be detected. Select from "midpoint" for setting the peak landmark to the midpoint of the peak region being detected or "mode" for setting the peak landmark to the mode location of the peak. "midpoint" can be generally more robust across samples and less impacted by the bandwidth. "mode" can be more accurate in determining the peak location if the bandwidth is generally ideal for the target marker.
#' @param multi_sample_per_batch Set it to TRUE to discard the positive peak that only appear in one sample per batch (sample number is >=3 per batch).
#' @param shoulder_valley Indicator to specify whether a shoulder valley is expected in case of the heavy right tail where the population of cells should be considered as a positive population.
#' @param shoulder_valley_slope The slope on the ADT marker density distribution to call shoulder valley.
#' @param valley_density_adjust Parameter for `density` function: bandwidth used is adjust*bw. This makes it easy to specify values like ‘half the default’ bandwidth.
#' @param landmark_align_type Algin the peak and valleys using one of the "negPeak", "negPeak_valley", "negPeak_valley_posPeak", and "valley" alignment modes.
#' @param midpoint_type Fill in the missing first valley by the midpoint of two positive peaks ("midpoint") or impute by other valleys ("valley").
#' @param neg_candidate_thres The upper bound for the negative peak. Users can refer to their IgG samples to obtain the minimal upper bound of the IgG sample peak. It can be one of the values of asinh(4/5+1), asinh(6/5+1), or asinh(8/5+1) if the right 95% quantile of IgG samples is large.
#' @param lower_peak_thres The minimal ADT marker densipty height of calling it a real peak. Set it to 0.01 to avoid a suspicious positive peak. Set it to 0.001 or smaller to include some small but tend to be real positive peaks, especially for markers like CD19.
#' @param brewer_palettes Set the color scheme of color brewer.
#' @param save_intermediate_rds Save the rds file for the intermediate objects.
#' @param save_intermediate_fig Save the density plot figure for checking the peak and valley location detection.
#' @param detect_outlier_valley Detect outlier valley and impute by the neighbor samples.
#' @param target_landmark_location Align the landmarks to a fixed location or, by default, align to the mean across samples for each landmark. The default value is NULL. Setting it to "fixed" will align the negative peak to 1 and the right-most positive peak to 3. Users can also assign a two-element vector indicating the location of the negative and most positive peaks to be aligned.
#' @param clean_adt_name Clean the ADT marker name
#' @param verbose Set the verbosity of the function.
#' @examples
#' \dontrun{
#' ADTnorm(
#'   cell_x_adt = cell_x_adt,
#'   cell_x_feature = cell_x_feature,
#'   save_outpath = save_outpath,
#'   study_name = study_name,
#'   marker_to_process = c("CD3", "CD4", "CD8")
#'  )
#' }
#' @export
#' @import dplyr ggplot2
ADTnorm = function(cell_x_adt = NULL, cell_x_feature = NULL, save_outpath = NULL, study_name = "ADTnorm", marker_to_process = NULL, exclude_zeroes = FALSE, bimodal_marker = NULL, trimodal_marker = NULL, positive_peak = NULL, bw_smallest_bi = 1.1, bw_smallest_tri = 0.8, bw_smallest_adjustments=list(), input_raw_counts = TRUE, quantile_clip=0, peak_type = "midpoint", multi_sample_per_batch = FALSE, shoulder_valley = FALSE, shoulder_valley_slope = -0.5, valley_density_adjust = 3, landmark_align_type = "negPeak_valley_posPeak", midpoint_type = "valley", neg_candidate_thres = asinh(8/5 + 1), lower_peak_thres = 0.001, brewer_palettes = "Set1", save_intermediate_rds = FALSE, save_intermediate_fig = TRUE, detect_outlier_valley = FALSE, target_landmark_location = NULL, clean_adt_name = FALSE, manual_peak_overrides=list(), verbose=FALSE){
    ## input parameter checking
    if(is.null(cell_x_adt)){
        stop("Please provide ADT raw count matrix in the cell as row and adt marker as column format.")
    }
    if(nrow(cell_x_adt) < ncol(cell_x_adt)){
        warning("Please check if the ADT raw count matrix has cell as row and adt marker as column.")
    }
    if(is.null(cell_x_feature)){
        stop("Please provide cell meta information including the sample and/or batch.")
    }
    if(clean_adt_name){
        colnames(cell_x_adt) = colnames(cell_x_adt) %>% clean_adt_name
        marker_to_process = marker_to_process %>% clean_adt_name
    }
    
    if(sum(!(bimodal_marker %in% colnames(cell_x_adt))) > 0){
        stop("Please provide consistent bimodal marker name as the column name of input ADT count matrix.")
    }
    if(save_intermediate_rds && is.null(save_outpath)){
        stop("Please provide the save_outpath to save the intermediate results in rds.")
    }
    if(save_intermediate_fig && is.null(save_outpath)){
        stop("Please provide the save_outpath to save the intermediate figures in pdf.")
    }
    if(!is.null(target_landmark_location)){
        if(target_landmark_location == "fixed"){
            print("Will align negative peak to 1 and right-most positive peak to 3.")
            target_landmark_location = c(1, 3)
        }else{
            if(length(target_landmark_location) == 2 && target_landmark_location[1] < target_landmark_location[2]){
                print(paste0("Will align negative peak to", target_landmark_location[1], " and right-most positive peak to ", target_landmark_location[2]))
            }else{
                stop("Please provide two elements vector to target_landmark_location where the first element is smaller!")
            }
        }
    }
    if(!(peak_type %in% c("mode", "midpoint"))){
        stop("Please specify the peak type to be either 'mode' or 'midpoint'.")
    }
    
    if(exclude_zeroes){
        na_mask = is.na(cell_x_adt)
        # Set all cell_x_adt to NA to avoid transformation during arcsin (if used)
        cell_x_adt[cell_x_adt==0] <- NA
    }
    
    ## preprocess the input data
    if(input_raw_counts){
        ## Check that cell_x_adt is integers only
        matrix = as.matrix(cell_x_adt)
        as_int = as.matrix(cell_x_adt)
        mode(as_int) <- "integer"
        if(!all(as_int[!is.na(as_int)] == matrix[!is.na(matrix)]) | !all(matrix[!is.na(matrix)]>=0)){
            stop("When using input_raw_counts, please only input positive integers for expression values.")
        }
        
        cell_x_adt = arcsinh_transform(cell_x_adt = cell_x_adt) ## Arcsinh transformation
    }
    all_marker_name = colnames(cell_x_adt) ## save the original marker name
    if(!is.factor(cell_x_feature$sample)){
        cell_x_feature$sample = factor(cell_x_feature$sample, levels = unique(cell_x_feature$sample))
    }
    ## get the index of bimodal marker
    if(is.null(bimodal_marker)){
        bimodal_marker = setdiff(all_marker_name, trimodal_marker)
    }
    bimodal_marker_index = which(all_marker_name %in% bimodal_marker)
    if(!is.null(trimodal_marker)){
        trimodal_marker_index = which(all_marker_name %in% trimodal_marker)
    }else{
        trimodal_marker_index = NULL
    }
    ## get the index of marker whose peak should be aligned to positive peak
    if(!is.null(positive_peak)){
        positive_peak_marker_index = c()
        for(adt_each in positive_peak[["ADT"]]){
            positive_peak_marker_index = c(positive_peak_marker_index, which(all_marker_name %in% adt_each))
        }
        positive_peak[["ADT_index"]] = positive_peak_marker_index
        
    }
    
    ## get the index of important lineage marker
    # if(is.null(cd3_index) && "CD3" %in% all_marker_name){
    #     cd3_index = which(all_marker_name == "CD3")
    # }
    # if(is.null(cd4_index) && "CD4" %in% all_marker_name){
    #     cd4_index = which(all_marker_name == "CD4")
    # }
    # if(is.null(cd8_index) && "CD8" %in% all_marker_name){
    #     cd8_index = which(all_marker_name == "CD8")
    # }
    
    colnames(cell_x_adt) = paste0("tmpName", 1:ncol(cell_x_adt)) ## replace the marker name by temp simple name to avoid symbol change by r
    
    
    ## get the ADT marker index that need normalization
    if(is.null(marker_to_process)){
        print(paste0("ADTnorm will process all the ADT markers from the ADT matrix:", paste(all_marker_name,  collapse = ", ")))
        adt_marker_index_list = 1:ncol(cell_x_adt)
        cell_x_adt_norm = cell_x_adt ## to record normalization results
    }else{
        print(paste0("ADTnorm will process the following ADT markers as provided:", paste(marker_to_process, collapse = ", ")))
        adt_marker_index_list = which(all_marker_name %in% marker_to_process)
        cell_x_adt_norm = cell_x_adt[, adt_marker_index_list] %>% t %>% t ## to record normalization results
        colnames(cell_x_adt_norm) = paste0("tmpName", 1:ncol(cell_x_adt))[adt_marker_index_list]
        if(length(adt_marker_index_list) == 0){
            print("Please provide consistent ADT marker name as the input ADT expression count matrix column name.")
        }
    }
    
    ## ADTnorm each marker
    for(adt_marker_index in adt_marker_index_list){ ## process each ADT marker respectively -- pending parallel running setup
        adt_marker_select = colnames(cell_x_adt)[adt_marker_index]
        adt_marker_select_name = all_marker_name[adt_marker_index]
        print(adt_marker_select_name)
        if(exclude_zeroes){
            exclude_zeroes_mask <- !(is.na(cell_x_adt[adt_marker_select]))
            cell_x_adt_total <- cell_x_adt
            cell_x_adt <- subset(cell_x_adt,exclude_zeroes_mask)
            cell_x_adt_norm_total <- cell_x_adt_norm
            cell_x_adt_norm <- subset(cell_x_adt_norm,exclude_zeroes_mask)
            cell_x_feature_total <- cell_x_feature
            cell_x_feature <- subset(cell_x_feature,exclude_zeroes_mask)
            if(verbose){
                print('Zeroes excluded...')
            }
        }
        if(quantile_clip>0){
            if(verbose){
                print('Performing quantile clip')
            }
            quant = quantile(cell_x_adt[[adt_marker_select]], 1 - quantile_clip)
            cell_x_adt[cell_x_adt[adt_marker_select]>quant,adt_marker_select] <- quant
        }
        ## smallest bw for density curve
        ## get the peak mode location
        
        bwFac_smallest = bw_smallest_bi
        if(adt_marker_index %in% trimodal_marker_index)
            bwFac_smallest = bw_smallest_tri
        if(adt_marker_select_name %in% names(bw_smallest_adjustments))
            bwFac_smallest = bw_smallest_adjustments[[adt_marker_select_name]]
        
        if(verbose)
            print(paste("Using bandwidth of",bwFac_smallest))
        
        if(peak_type == "mode"){
            peak_mode_res <- get_peak_mode(cell_x_adt, cell_x_feature, adt_marker_select, adt_marker_index, bwFac_smallest, bimodal_marker_index, trimodal_marker_index, positive_peak, neg_candidate_thres = neg_candidate_thres, lower_peak_thres = lower_peak_thres)
        } else if(peak_type == "midpoint"){
            peak_mode_res <- get_peak_midpoint(cell_x_adt, cell_x_feature, adt_marker_select, adt_marker_index, bwFac_smallest, bimodal_marker_index, trimodal_marker_index, positive_peak, neg_candidate_thres = neg_candidate_thres, lower_peak_thres = lower_peak_thres)
        }
        if(adt_marker_select_name %in% names(manual_peak_overrides)){
            peak_overrides <- manual_peak_overrides[[adt_marker_select_name]]
            peak_mode_res
            for(sample in rownames(peak_mode_res)){
                if(sample %in% names(peak_overrides)){
                    len_needed = length(peak_mode_res[sample,])
                    len_provided = length(peak_overrides[[sample]])
                    overrides <- peak_overrides[[sample]][0:len_needed]
                    ddata <- cell_x_adt[(cell_x_feature$sample == sample),c(adt_marker_select)] 
                    if ((max(overrides,na.rm = TRUE)>max(ddata,na.rm = TRUE)) | (min(overrides,na.rm = TRUE) < min(ddata,na.rm = TRUE)) | (is.unsorted(overrides,na.rm=TRUE,TRUE))){
                        stop(paste0("Peak overrides must be within the bounds of the data (",min(ddata,na.rm = TRUE),', ',max(ddata,na.rm = TRUE),") and must be strictly increasing.",
                                    "\n You submitted: ",list(overrides)," for marker: ", adt_marker_select_name))
                    }
                    peak_mode_res[sample,] <- overrides
                }
            }
            
        }
        if(verbose){
            print(peak_mode_res)
        }
        
        ## get the valley location
        peak_valley_list <- get_valley_location(cell_x_adt, cell_x_feature, adt_marker_select, peak_mode_res, shoulder_valley, positive_peak, multi_sample_per_batch, adjust = valley_density_adjust, min_fc = 20, shoulder_valley_slope = shoulder_valley_slope, neg_candidate_thres = neg_candidate_thres)
        
        valley_location_res <- peak_valley_list$valley_location_list
        peak_mode_res <- peak_valley_list$peak_landmark_list
        
        if(detect_outlier_valley){ ## detect if valley is outlier and impute by neighbor samples if found
            valley_location_res <- detect_impute_outlier_valley(valley_location_res, adt_marker_select, cell_x_adt, cell_x_feature, scale = 3, method = "MAD", nearest_neighbor_n = 3, nearest_neighbor_threshold = 0.75)
        }
        
        ## density plot for peak and valley location checking
        if(input_raw_counts){
            method_label <- "Arcsinh Transformation"
            figure_label <- "ArcsinhTransformation"
        }else{
            method_label <- "Input Transformation"
            figure_label <- "InputTransformation"
        }
        density_plot <- plot_adt_density_with_peak_valley_each(cell_x_adt[, adt_marker_select], cell_x_feature, peak_landmark_list = peak_mode_res, valley_landmark_list = valley_location_res, brewer_palettes = brewer_palettes, parameter_list = list(
            bw = 0.1,
            method_label = method_label
        ))
        
        ## save all three objects
        peak_valley <- list(
            peak_landmark_list = peak_mode_res,
            valley_landmark_list = valley_location_res
        )
        if(save_intermediate_rds){
            if(!dir.exists(paste0(save_outpath, "/RDS"))){
                dir.create(paste0(save_outpath, "/RDS"), recursive = TRUE)
            }
            
            saveRDS(peak_valley, file = paste0(save_outpath, "/RDS/peak_valley_raw_", adt_marker_select_name, "_", study_name, ".rds"), compress = FALSE)
            saveRDS(density_plot, file = paste0(save_outpath, "/RDS/density_raw_", adt_marker_select_name, "_", study_name, ".rds"), compress = FALSE)
            
        }
        if(save_intermediate_fig){
            if(!dir.exists(paste0(save_outpath, "/figures"))){
                dir.create(paste0(save_outpath, "/figures"), recursive = TRUE)
            }
            
            grDevices::pdf(paste0(save_outpath, "/figures/ArcsinhTransform_", adt_marker_select_name, "_", study_name, ".pdf"), width = 11, height = ceiling(length(levels(cell_x_feature$sample)) * 0.4))
            print(density_plot)
            grDevices::dev.off()
        }
        
        
        landmark_matrix <- landmark_fill_na(
            peak_landmark_list = peak_mode_res,
            valley_landmark_list = valley_location_res,
            landmark_align_type = landmark_align_type,
            midpoint_type = midpoint_type,
            neg_candidate_thres = neg_candidate_thres
        )
        if(!is.null(target_landmark_location)){ ## specify the locations where negative peak, valley, and positive peak should align to
            if(landmark_align_type == "negPeak_valley_posPeak"){
                if(ncol(peak_mode_res) == 1){
                    target_landmark = c(target_landmark_location[1], round(mean(target_landmark_location), 1))
                }else if(ncol(peak_mode_res) == 2){
                    target_landmark = c(target_landmark_location[1], round(mean(target_landmark_location), 1), target_landmark_location[2])
                }else if(ncol(peak_mode_res) > 2){
                    target_landmark = c(target_landmark_location[1], (target_landmark_location[1] + round(mean(target_landmark_location), 1))/2, target_landmark_location[2])
                }
            }else if(landmark_align_type == "negPeak_valley"){
                target_landmark = c(target_landmark_location[1], round(mean(target_landmark_location), 1))
            }else if(landmark_align_type == "negPeak"){
                target_landmark = c(target_landmark_location[1])
            }else if(landmark_align_type == "valley"){
                target_landmark = c(round(mean(target_landmark_location), 1))
            }else {
                stop("Please provide one of the landmark_align_type from: negPeak, negPeak_valley, negPeak_valley_posPeak, valley")
            }
            
        }else{
            target_landmark = NULL
        }
        peak_alignment_res = peak_alignment(cell_x_adt[, adt_marker_select], cell_x_feature, landmark_matrix, target_landmark = target_landmark)
        cell_x_adt_norm[, adt_marker_select] = peak_alignment_res[[1]]
        
        if(ncol(peak_alignment_res[[2]]) == 2){
            peak_mode_norm_res = peak_alignment_res[[2]][, 1,drop=FALSE]
            valley_location_norm_res = peak_alignment_res[[2]][, 2,drop=FALSE]
        }else if(ncol(peak_alignment_res[[2]]) == 3){
            peak_mode_norm_res = peak_alignment_res[[2]][, c(1, 3), drop=FALSE]
            valley_location_norm_res = peak_alignment_res[[2]][, 2,drop=FALSE]
        }else if(ncol(peak_alignment_res[[2]]) == 5){
            peak_mode_norm_res = peak_alignment_res[[2]][, c(1, 3, 5),drop=FALSE]
            valley_location_norm_res = peak_alignment_res[[2]][, c(2, 4),drop=FALSE]
        }
        density_norm_plot <- plot_adt_density_with_peak_valley_each(cell_x_adt_norm[, adt_marker_select], cell_x_feature, peak_landmark_list = peak_mode_norm_res, valley_landmark_list = valley_location_norm_res, brewer_palettes = brewer_palettes, parameter_list = list(
            bw = 0.2,
            method_label = "ADTnorm"
        ))
        if(save_intermediate_rds){
            saveRDS(density_norm_plot, file = paste0(save_outpath, "/RDS/density_ADTnorm_", adt_marker_select_name, "_", study_name, ".rds"), compress = FALSE)
        }
        if(save_intermediate_fig){
            grDevices::pdf(paste0(save_outpath, "/figures/ADTnorm_", adt_marker_select_name, "_", study_name, ".pdf"), width = 11, height = ceiling(length(levels(cell_x_feature$sample)) * 0.4))
            print(density_norm_plot)
            grDevices::dev.off()
        }
        if(exclude_zeroes){ # Insert zeroes back into their proper places
            cell_x_adt_total[exclude_zeroes_mask,] <- cell_x_adt
            cell_x_adt_norm_total[exclude_zeroes_mask,] <- cell_x_adt_norm
            cell_x_adt <- cell_x_adt_total
            cell_x_adt_norm <- cell_x_adt_norm_total
            cell_x_feature <- cell_x_feature_total
        }
    }
    
    if(exclude_zeroes){
        # Set all cell_x_adt NAs back to 0
        cell_x_adt_norm[is.na(cell_x_adt_norm)] <- 0
        cell_x_adt_norm[na_mask] <- NA
    }
    colnames(cell_x_adt_norm) = all_marker_name[adt_marker_index_list]
    return(cell_x_adt_norm)
    
}
