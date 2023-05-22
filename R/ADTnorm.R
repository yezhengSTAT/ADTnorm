#' ADTnorm normalization to remove the technical variations across samples for each ADT marker.
#'
#' This function removes the technical variations such as batch effect, sequencing depth biases, antibody selection difference and antibody concentration differences, etc. The normalized samples are ready for integration across studies.
#' @param cell_x_adt Matrix of ADT raw counts in cells (rows) by ADT markers (columns) format. By default, ADTnorm expects raw counts to be provided and arcsin transformation to be performed by ADTnorm internally. If ADTnorm detects that the input count matrix is non-integer, it will skipped the arcsine transformation and user need to tune the parameters to fit their input transformation.
#' @param cell_x_feature Matrix of cells (rows) by cell features (columns) such as sample, batch or other cell-type related information. Please ensure that cell_x_feature matrix at least contains sample and batch columns with the exact "sample" and "batch" column name. Also one sample should only has one unique batch label to avoid confusion in the final normalization results visualization.
#' @param save_outpath The path to save the results.
#' @param study_name Name of this run.
#' @param marker_to_process Markers to normalize. Leaving empty to process all the ADT markers in cell_x_adt matrix.
#' @param exclude_zeroes Indicator to consider zeros as NA, i.e., missing values. Recommend TRUE if zeroes in the data represent dropout, likely for large ADT panels, big datasets, or undersequenced data. Default is FALSE.
#' @param bimodal_marker Specify ADT markers that are likely to have two peaks based on researchers' prior knowledge or preliminary observation of particular data to be processed. Leaving it as default, ADTnorm will try to find the bimodal peak in all markers that are not listed in `trimodal_marker.`
#' @param trimodal_marker Index of the ADT markers that tend to have three peaks based on researchers' prior knowledge (e.g., CD4) or preliminary observation on particular data to be processed.
#' @param positive_peak A list variable containing a vector of ADT marker(s) and a corresponding vector of sample name(s) in matching order to specify that the uni-peak detected should be aligned to positive peaks. For example, for samples that only contain T cells. The only CD3 peak should be aligned to the positive peaks of other samples.
#' @param bw_smallest_bi The smallest bandwidth parameter value for bi-modal peaks. Recommend 1.1.
#' @param bw_smallest_tri The smallest bandwidth parameter value for tri-modal peaks. Recommend the same value for CD4, such as 0.5.
#' @param bw_smallest_adjustments A named list of floats, with names matching marker names, specifying the smallest bandwidth parameter value. The default value is bw_smallest_adjustments = list(CD3 = 0.8, CD4 = 0.8, CD8 = 0.8). Recommend 0.5 or 0.8 for common multi-modal marker.
#' @param quantile_clip Implement an upper quantile clipping to avoid warping function errors caused by outlier measurement of extremely high expression. Provide the quantile threshold to remove outlier points above such qunatile. Default is 1, meaning no filtering. 0.99 means 99th quantile and points above 99th quantile will be discard.
#' @param peak_type The type of peak to be detected. Select from "midpoint" for setting the peak landmark to the midpoint of the peak region being detected or "mode" for setting the peak landmark to the mode location of the peak. "midpoint" can be generally more robust across samples and less impacted by the bandwidth. "mode" can be more accurate in determining the peak location if the bandwidth is generally ideal for the target marker.
#' @param multi_sample_per_batch Set it to TRUE to discard the positive peak that only appear in one sample per batch (sample number is >=3 per batch).
#' @param shoulder_valley Indicator to specify whether a shoulder valley is expected in case of the heavy right tail where the population of cells should be considered as a positive population. Default is TRUE.
#' @param shoulder_valley_slope The slope on the ADT marker density distribution to call shoulder valley. Default is -0.5
#' @param valley_density_adjust Parameter for `density` function: bandwidth used is adjust*bw. This makes it easy to specify values like ‘half the default’ bandwidth. Default is 3.
#' @param landmark_align_type Algin the peak and valleys using one of the "negPeak", "negPeak_valley", "negPeak_valley_posPeak", and "valley" alignment modes. Default is "negPeak_valley_posPeak".
#' @param midpoint_type Fill in the missing first valley by the midpoint of two positive peaks ("midpoint") or impute by other valleys ("valley"). Default is valley.
#' @param neg_candidate_thres The upper bound for the negative peak. Users can refer to their IgG samples to obtain the minimal upper bound of the IgG sample peak. It can be one of the values of asinh(4/5+1), asinh(6/5+1), or asinh(8/5+1) if the right 95% quantile of IgG samples is large. Default is asinh(8/5+1) for raw count input. This filtering will be disabled if input is not raw count data. 
#' @param lower_peak_thres The minimal ADT marker density height of calling it a real peak. Set it to 0.01 to avoid a suspicious positive peak. Set it to 0.001 or smaller to include some small but tend to be real positive peaks, especially for markers like CD19. Default is 0.001.
#' @param brewer_palettes Set the color scheme of color brewer. Default is "Set1".
#' @param save_intermediate_rds Save the rds file for the intermediate objects. Default is FALSE
#' @param save_intermediate_fig Save the density plot figure for checking the peak and valley location detection. Default is TRUE.
#' @param detect_outlier_valley Detect outlier valley and impute by the neighbor samples. Outlier detection methods, choose from "MAD" (Median Absolute Deviation) or "IQR" (InterQuartile Range). Recommend trying "MAD" first if needed. Default is FALSE.
#' @param target_landmark_location Align the landmarks to a fixed location or, by default, align to the mean across samples for each landmark. The default value is NULL. Setting it to "fixed" will align the negative peak to 1 and the right-most positive peak to 5. Users can also assign a two-element vector indicating the location of the negative and most positive peaks to be aligned.
#' @param clean_adt_name Clean the ADT marker name. Default is FALSE.
#' @param customized_landmark By setting it to be TRUE, ADTnorm will trigger the interactive landmark tuning function and pop out a shiny application for user's manual setting of the peaks and valleys location. We recommend using this function after initial rounds of ADTnorm normalization with a few parameters tuning attempt. It is better to narrow down a few ADT markers that do need manual tuning and provide the list to marker_to_process as the interactive function will pop out for every marker being processed. Default is FALSE.
#' @param verbose Set the verbosity of the function. Default is FALSE.
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
#' @import dplyr ggplot2 shiny
ADTnorm = function(cell_x_adt = NULL, cell_x_feature = NULL, save_outpath = NULL, study_name = "ADTnorm", marker_to_process = NULL, exclude_zeroes = FALSE, bimodal_marker = NULL, trimodal_marker = NULL, positive_peak = NULL, bw_smallest_bi = 1.1, bw_smallest_tri = 0.8, bw_smallest_adjustments = list(CD3 = 0.8, CD4 = 0.8, CD8 = 0.8), quantile_clip = 1, peak_type = "midpoint", multi_sample_per_batch = FALSE, shoulder_valley = TRUE, shoulder_valley_slope = -0.5, valley_density_adjust = 3, landmark_align_type = "negPeak_valley_posPeak", midpoint_type = "valley", neg_candidate_thres = NULL, lower_peak_thres = 0.001, brewer_palettes = "Set1", save_intermediate_rds = FALSE, save_intermediate_fig = TRUE, detect_outlier_valley = FALSE, target_landmark_location = NULL, clean_adt_name = FALSE, customized_landmark = FALSE, verbose = FALSE){
    
    ## =========================
    ## INPUT PARAMETER CHECKING
    ## =========================
    if(is.null(cell_x_adt)){ ## if ADT count matrix is provided
        stop("Please provide ADT raw count matrix in the cell as row and adt marker as column format.")
    }else{
        matrix = as.matrix(cell_x_adt)
        as_int = as.matrix(cell_x_adt)
        mode(as_int) = "integer"
        if(any(as_int[!is.na(as_int)] != matrix[!is.na(matrix)])){ ## if any input count is not integer
            print("Note: Input ADT count matrix values are not all integer hence the input is NOT considered as raw count data. ADTnorm will assume transformation or scaling has been done and will skip the arcsine transformation......")
            arcsine_transform_flag = FALSE
            if(is.null(neg_candidate_thres)){
                print("neg_candidate_thres is not set hence the suspecious negative peak candidate filtering function will be diabled. Come back to set it if several negative peaks are mistakenly detected due to the discrete values around 0.")
                neg_candidate_thres = min(matrix) ## disable the negative peak candidates filtering
            }
           

        }else{ ## all input counts are integer, check if counts are nonnegative
            if(any(matrix[!is.na(matrix)] < 0)){
                stop("ADTnorm detects integer input ADT count matrix and hence assume input is raw count. Please check why there is non-negative value in the raw count matrix.")
            }
            arcsine_transform_flag = TRUE

            if(is.null(neg_candidate_thres)){
                neg_candidate_thres = asinh(8/5 + 1) ## give the default value for neg_candidate_thres for raw count input
            }
            
        }
    }
    if(nrow(cell_x_adt) < ncol(cell_x_adt)){ ## ensure the input matrix in the right format: cell as row and adt marker as column
        stop("Please check if the ADT raw count matrix has cell as row and adt marker as column.")
    }
    if(is.null(cell_x_feature)){ ## check if feature information is provided that include sample and/or batch
        stop("Please provide cell meta information including the sample and/or batch.")
    }else{
        if(!("sample" %in% colnames(cell_x_feature))){ ## ensure that feature matrix has a column indicating the sample label
            stop("Please provide the label for each sample and name the column 'sample'. sample is the smallest unit that you want to group the cells. It can be a donor or a donor under one condition or within one batch run. Please note that one sample should only has one unique batch label to avoid confusion in the normalization.")
        }else{
            print("Reminder: Please note that one sample should only has one unique batch label to avoid confusion in the normalization and visualization......")
        }
        if(!("batch" %in% colnames(cell_x_feature))){ ## check 
            stop("Please provide the label for each batch run and name the column 'batch'. batch refers to the experimental run or scenarios where unwanted technical noise was introduced. Please note that one sample should only has one unique batch label to avoid confusion in the normalization.")
        }
        if(nrow(cell_x_feature) != nrow(cell_x_adt)){
            stop("Please check if the cell number matches between cell_x_adt and cell_x_feature matrices.")
        }
    }
    if(clean_adt_name){ ## if user wants to clean the adt name
        colnames(cell_x_adt) = colnames(cell_x_adt) %>% clean_adt_name
        marker_to_process = marker_to_process %>% clean_adt_name
    }
    if(any(!(marker_to_process %in% colnames(cell_x_adt)))){ ## check if all the markers required normalization exist in the cell_x_adt
        print("Warning: Please double check the marker_to_process value or the column name of cell_x_adt. The following markers are not provided in the cell_x_adt input matrix:")
        not_exist_marker = marker_to_process[which(!(marker_to_process %in% colnames(cell_x_adt)))]
        print(paste(not_exist_marker, collapse = ", "))
        print("The remaining existing markers will be processed. Stop the job and rerun if you want to modify marker_to_process......")
    }
    if(sum(!(bimodal_marker %in% colnames(cell_x_adt))) > 0){ ## check if adt marker provided are not among the input adt count matrix columns.
        stop("Please provide consistent bimodal marker name as the column name of input ADT count matrix.")
    }
    if(sum(!(trimodal_marker %in% colnames(cell_x_adt))) > 0){ ## check if adt marker provided are not among the input adt count matrix columns.
        stop("Please provide consistent trimodal marker name as the column name of input ADT count matrix.")
    }
    if(save_intermediate_rds && is.null(save_outpath)){ ## check the output path to save intermediate data
        stop("Please provide the save_outpath to save the intermediate results in rds.")
    }
    if(save_intermediate_fig && is.null(save_outpath)){ ## check the output path to save intermediate figures
        stop("Please provide the save_outpath to save the intermediate figures in pdf.")
    }
    if(!is.null(target_landmark_location)){ ## Check if user wants to align to a fixed location. 
        if(target_landmark_location == "fixed"){ ## Currently default fixed alignemnt locations are set to 1 and 5 for better visualization. 
            print("Will align negative peak to 1 and right-most positive peak to 5.")
            target_landmark_location = c(1, 5)
        }else{ ## If user provide the fixed alignment location, align to the user-set locations to align the negative and the right-most positive peak.
            if(length(target_landmark_location) == 2 && target_landmark_location[1] < target_landmark_location[2]){
                print(paste0("Will align negative peak to", target_landmark_location[1], " and right-most positive peak to ", target_landmark_location[2]))
            }else{
                stop("Please provide two elements vector to target_landmark_location where the first element is smaller!")
            }
        }
    }
    if(!(peak_type %in% c("mode", "midpoint"))){ ## To define the peak landmark. We can use the midpoint of peak region (more robust) or the mode of the peak (more accurate)
        stop("Please specify the peak type to be either 'mode' or 'midpoint'.")
    }
    
    ## ==============================
    ## DATA PROCESSING AND CLEANNING
    ## ==============================
    ## Remove ADT marker if it only has zero value across all the cells
    col_sums = colSums(cell_x_adt, na.rm = TRUE)
    if (any(col_sums == 0)){
        message("Markers with zero counts will be ignored")
        cell_x_adt = cell_x_adt[, col_sums > 0, drop = FALSE] 
    }

    ## Whether to consider zero as missing value and set to NA.
    if(exclude_zeroes){
        na_mask = is.na(cell_x_adt)
        # Set all cell_x_adt to NA to avoid transformation during arcsin (if used)
        cell_x_adt[cell_x_adt == 0] = NA

    }

    ## Arcsinh transformation with default shift factor 1 and scale factor 1/5
    if(arcsine_transform_flag){
        cell_x_adt = arcsinh_transform(cell_x_adt = cell_x_adt) 
    }

    ## save the original marker name
    all_marker_name = colnames(cell_x_adt)
    
    ## factorize the sample label
    if(!is.factor(cell_x_feature$sample)){
        cell_x_feature$sample = factor(cell_x_feature$sample, levels = unique(cell_x_feature$sample))
    }

    ## get the index of bimodal markers and trimodal markers
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
    
    ## replace the marker name by temp simple name to avoid symbol change by r
    colnames(cell_x_adt) = paste0("tmpName", 1:ncol(cell_x_adt)) 
    
    
    ## get the ADT marker index that need normalization
    if(is.null(marker_to_process)){
        print(paste0("ADTnorm will process all the ADT markers from the ADT matrix: ", paste(all_marker_name,  collapse = ", ")))
        adt_marker_index_list = 1:ncol(cell_x_adt)
        cell_x_adt_norm = cell_x_adt ## to record normalization results, other markers that are not listed for normalization will remain the same as input value.
    }else{
        adt_marker_index_list = which(all_marker_name %in% marker_to_process)
        if(length(adt_marker_index_list) == 0){ ## may not need it anymore 
            stop("Please provide consistent ADT marker name as the input ADT expression count matrix column name.")
        }else{
            print(paste0("Note: ADTnorm will process the following ADT markers as provided: ", paste(all_marker_name[adt_marker_index_list], collapse = ", ")))
            cell_x_adt_norm = cell_x_adt[, adt_marker_index_list, drop=FALSE] #%>% t %>% t ## to record normalization results
            colnames(cell_x_adt_norm) = paste0("tmpName", 1:ncol(cell_x_adt))[adt_marker_index_list]
        }
    }
    
    ## ADTnorm each marker
    for(adt_marker_index in adt_marker_index_list){ ## process each ADT marker respectively -- pending parallel running setup
        adt_marker_select = colnames(cell_x_adt)[adt_marker_index] ## temp name
        adt_marker_select_name = all_marker_name[adt_marker_index] ## original name
        
        if(verbose){
            print(paste0("Normalizing marker: ", adt_marker_select_name))
        }


        if(quantile_clip < 1){ 
            quant = quantile(cell_x_adt[[adt_marker_select]], quantile_clip, na.rm = TRUE)
            if(verbose){
                print(paste0("Performing quantile clip to remove outliers beyond the claimed quantile ", quantile_clip, ":", quant, "......"))
            }
            
            cell_x_adt[[adt_marker_select]][which(cell_x_adt[[adt_marker_select]] > quant)] <- NA
        }
        
        ## smallest bw for density curve
        bwFac_smallest = bw_smallest_bi
        if(adt_marker_index %in% trimodal_marker_index){
            bwFac_smallest = bw_smallest_tri
        }   
        if(adt_marker_select_name %in% names(bw_smallest_adjustments)){
            bwFac_smallest = bw_smallest_adjustments[[adt_marker_select_name]]
        }
        
        if(verbose){
            print(paste0("Using bandwidth for constructing density distribution: ", bwFac_smallest))
        }

        ## get the peak mode location    
        if(peak_type == "mode"){
            peak_mode_res = get_peak_mode(cell_x_adt, cell_x_feature, adt_marker_select, adt_marker_index, bwFac_smallest, bimodal_marker_index, trimodal_marker_index, positive_peak, neg_candidate_thres = neg_candidate_thres, lower_peak_thres = lower_peak_thres, arcsine_transform_flag = arcsine_transform_flag)
        } else if(peak_type == "midpoint"){
            peak_mode_res = get_peak_midpoint(cell_x_adt, cell_x_feature, adt_marker_select, adt_marker_index, bwFac_smallest, bimodal_marker_index, trimodal_marker_index, positive_peak, neg_candidate_thres = neg_candidate_thres, lower_peak_thres = lower_peak_thres, arcsine_transform_flag = arcsine_transform_flag)
        }

        ## get the valley location
        peak_valley_list = get_valley_location(cell_x_adt, cell_x_feature, adt_marker_select, peak_mode_res, shoulder_valley, positive_peak, multi_sample_per_batch, adjust = valley_density_adjust, min_fc = 20, shoulder_valley_slope = shoulder_valley_slope, neg_candidate_thres = neg_candidate_thres, arcsine_transform_flag = arcsine_transform_flag)
        
        valley_location_res = peak_valley_list$valley_location_list
        peak_mode_res = peak_valley_list$peak_landmark_list
        
        if(detect_outlier_valley != FALSE){ ## detect if valley is outlier and impute by neighbor samples if found
            valley_location_res = detect_impute_outlier_valley(valley_location_res, adt_marker_select, cell_x_adt, cell_x_feature, scale = 3, method = detect_outlier_valley, nearest_neighbor_n = 3, nearest_neighbor_threshold = 0.75)
        }
        ## Manual overwrite the peak location - peak and valley 
        if(customized_landmark){
            cell_x_adt_sample = data.frame(adt = cell_x_adt[, adt_marker_select], sample = cell_x_feature[, "sample"], batch = cell_x_feature[, "batch"])
            num_landmark = ncol(peak_mode_res) + ncol(valley_location_res)

            landmark_pos = matrix(nrow = nrow(peak_mode_res), ncol = num_landmark)
            landmark_pos[, seq(1, num_landmark, 2)] = peak_mode_res
            landmark_pos[, seq(2, num_landmark, 2)] = valley_location_res
            rownames(landmark_pos) = rownames(peak_mode_res)
            colnames(landmark_pos)[seq(1, num_landmark, 2)] = paste0("peak", 1:ncol(peak_mode_res))
            colnames(landmark_pos)[seq(2, num_landmark, 2)] = paste0("valley", 1:ncol(valley_location_res))

            landmark_pos_customized = get_customize_landmark(cell_x_adt_sample, landmark_pos, bw = 0.2, adt_marker_select_name = adt_marker_select_name)
            peak_mode_res = landmark_pos_customized[, seq(1, num_landmark, 2), drop = FALSE]
            valley_location_res = landmark_pos_customized[, seq(2, num_landmark, 2), drop = FALSE]

            peak_mode_res_na <- apply(peak_mode_res, 2, function(x) all(is.na(x)))
            if(any(peak_mode_res_na)){
                peak_mode_res_filter <- peak_mode_res[, -which(peak_mode_res_na), drop = FALSE]
                peak_mode_res <- peak_mode_res_filter
            }
            if(verbose){
                print("Customized landmark positions:")
                print(landmark_pos_customized)
            }
        }   

        
        ## density plot for peak and valley location checking
        if(arcsine_transform_flag){
            run_label = "Arcsinh Transformation"
        }else{
            run_label = "Input Transformation"
        }
        density_plot = plot_adt_density_with_peak_valley_each(cell_x_adt[, adt_marker_select], cell_x_feature, peak_landmark_list = peak_mode_res, valley_landmark_list = valley_location_res, brewer_palettes = brewer_palettes, parameter_list = list(
            bw = 0.1,
            run_label = run_label
        ))
        
        ## save peak and valley objects
        peak_valley = list(
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
        
        
        landmark_matrix = landmark_fill_na(
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
            peak_mode_norm_res = peak_alignment_res[[2]][, 1, drop = FALSE]
            valley_location_norm_res = peak_alignment_res[[2]][, 2, drop = FALSE]
        }else if(ncol(peak_alignment_res[[2]]) == 3){
            peak_mode_norm_res = peak_alignment_res[[2]][, c(1, 3), drop = FALSE]
            valley_location_norm_res = peak_alignment_res[[2]][, 2, drop = FALSE]
        }else if(ncol(peak_alignment_res[[2]]) == 5){
            peak_mode_norm_res = peak_alignment_res[[2]][, c(1, 3, 5), drop = FALSE]
            valley_location_norm_res = peak_alignment_res[[2]][, c(2, 4), drop = FALSE]
        }
        density_norm_plot = plot_adt_density_with_peak_valley_each(
            cell_x_adt_norm[, adt_marker_select], cell_x_feature, peak_landmark_list = peak_mode_norm_res, valley_landmark_list = valley_location_norm_res, brewer_palettes = brewer_palettes, 
            parameter_list = list(
                bw = 0.2,
                run_label = "ADTnorm"
            )
        )
        if(save_intermediate_rds){
            saveRDS(density_norm_plot, file = paste0(save_outpath, "/RDS/density_ADTnorm_", adt_marker_select_name, "_", study_name, ".rds"), compress = FALSE)
        }
        if(save_intermediate_fig){
            grDevices::pdf(paste0(save_outpath, "/figures/ADTnorm_", adt_marker_select_name, "_", study_name, ".pdf"), width = 11, height = ceiling(length(levels(cell_x_feature$sample)) * 0.4))
            print(density_norm_plot)
            grDevices::dev.off()
        }

    }
    colnames(cell_x_adt_norm) = all_marker_name[adt_marker_index_list]
    return(cell_x_adt_norm)
  
    
}

