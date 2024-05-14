#' Get the peak landmarks location using the peak region midpoint
#'
#' This function detects the peak landmark locations for each sample per ADT markers based on the midpoint of the peak region. Using the peak midpoint instead of the peak mode can be more stable across samples and less affected by the bandwidth.
#' @param cell_x_adt Matrix of ADT raw counts in cells (rows) by ADT markers (columns) format.
#' @param cell_x_feature Matrix of cells (rows) by cell features (columns) such as cell type, sample, and batch-related information.
#' @param adt_marker_select Markers to normalize. Leaving empty to process all the ADT markers in the cell_x_adt matrix.
#' @param adt_marker_index Index of the ADT markers that will be normalized. Leaving empty if adt_marker_select is not specified.
#' @param bwFac_smallest The smallest bandwidth parameter value. Recommend 1.1 for general bi-modal ADT markers except for CD3, CD4, and CD8.
#' @param bimodal_marker_index Index of the ADT markers that tend to have two peaks based on researchers' prior knowledge or preliminary observation on particular data to be processed.
#' @param trimodal_marker_index Index of the ADT markers that tend to have three peaks based on researchers' prior knowledge (e.g., CD4) or preliminary observation on particular data to be processed.
#' @param positive_peak A list variable containing a vector of ADT marker(s) and a corresponding vector of sample name(s) in matching order to specify that the uni-peak detected should be aligned to positive peaks. For example, for samples that only contain T cells. The only CD3 peak should be aligned with the positive peaks of other samples.
#' @param neg_candidate_thres The upper bound for the negative peak. Users can refer to their IgG samples to obtain the minimal upper bound of the IgG sample peak. It can be one of the values of asinh(4/5+1), asinh(6/5+1), or asinh(8/5+1) if the right 95% quantile of IgG samples is large.
#' @param lower_peak_thres The minimal ADT marker density height to call it a real peak. Set it to 0.01 to avoid a suspicious positive peak. Set it to 0.001 or smaller to include some small but tend to be real positive peaks, especially for markers like CD19.
#' @param arcsinh_transform_flag The flag indicates if the input is raw count and arcsinh transformation is implemented.
#' @export
#' @examples
#' \dontrun{
#' get_peak_midpoint(
#'   cell_x_adt = cell_x_adt,
#'   cell_x_feature = cell_x_feature,
#'   adt_marker_select = "CD3",
#'   neg_candidate_thres = asinh(6/5+1)
#' )
#' }
# require(flowStats)
# require(dplyr)
get_peak_midpoint = function(cell_x_adt = NULL, cell_x_feature = NULL, adt_marker_select = NULL, adt_marker_index = NULL, bwFac_smallest = 1.1, bimodal_marker_index = NULL, trimodal_marker_index = NULL, positive_peak = NULL, neg_candidate_thres = asinh(10/5 + 1), lower_peak_thres = 0.001, arcsinh_transform_flag = TRUE) {

    if (length(adt_marker_select) > 1){
        stop("adt_marker_select should be a single marker name")
    }
    
    ## set parameters
    bwFac = 1.2
    border = 0.01
    peakNr = NULL
    densities = NULL
    n = 201
    indices = FALSE

    ## get peak mode for each sample of this processing ADT marker
    sample_name_list = levels(cell_x_feature$sample) ## user provided or auto-detected
    ## get ADT value range with a slightly extension
    from = min(cell_x_adt[, adt_marker_select], na.rm = TRUE) - diff(range(cell_x_adt[, adt_marker_select], na.rm = TRUE)) * 0.15
    to = max(cell_x_adt[, adt_marker_select], na.rm = TRUE) + diff(range(cell_x_adt[, adt_marker_select], na.rm = TRUE)) * 0.15

    peak_num = 0
    peak_mode = list()
    peak_region = list()

    for(sample_name in sample_name_list){
        ## extract the ADT counts for this sample
        cell_ind_tmp = which(cell_x_feature$sample == sample_name)
        cell_notNA = which(!is.na(cell_x_adt[cell_ind_tmp, adt_marker_select]))
        cell_ind = cell_ind_tmp[cell_notNA]
        if(length(cell_ind) > 0){
            # fcs_count = cell_x_adt[cell_ind, ] %>% as.matrix()
            # fcs = flowCore::flowFrame(fcs_count)
            fcs_count = cell_x_adt[cell_ind, adt_marker_select, drop = FALSE] %>% as.matrix
            colnames(fcs_count) = adt_marker_select
            fcs = flowCore::flowFrame(fcs_count)

            unique_value = length(unique(cell_x_adt[cell_ind, adt_marker_select]))
            if(unique_value == 1){
                    ## only one value for this marker
                    peak_mode[[sample_name]] = NA
                    peak_region[[sample_name]] = matrix(NA, ncol = 2, nrow = 1)
                    print(paste0("Warning: ", sample_name, " only has a single unique value!"))
            } else {
                ## get the proportion of cells near zero to diagnoise the negative peak enrichment
                if(arcsinh_transform_flag){
                    zero_prop = sum(cell_x_adt[cell_ind, adt_marker_select] < 2) / length(cell_x_adt[cell_ind, adt_marker_select])
                }else{
                    zero_prop = sum(cell_x_adt[cell_ind, adt_marker_select] < neg_candidate_thres) / length(cell_x_adt[cell_ind, adt_marker_select])
                }
                
                # zero_prop_list[[sample_name]] = zero_prop ## zero proportion
                adt_expression = cell_x_adt[cell_ind, adt_marker_select] ## adt value for this marker and this sample

                ## if most are around 0 and there are very few unique value: add random small number
                if(zero_prop > 0.95){
                    if(length(unique(adt_expression)) < 50){
                        adt_expression = adt_expression + stats::rnorm(length(adt_expression), mean = 0, sd = 0.05)
                        # exprs(dat[[sample_name]])[, adt_marker_select] = adt_expression
                    }
                }
                fres1 = flowCore::filter(fcs, flowStats::curv1Filter(adt_marker_select, bwFac = max(2, bwFac_smallest)))
                fres2 = flowCore::filter(fcs, flowStats::curv1Filter(adt_marker_select, bwFac = max(3, bwFac_smallest)))
                fres3 = flowCore::filter(fcs, flowStats::curv1Filter(adt_marker_select, bwFac = max(3.1, bwFac_smallest)))

                ## different bandwidth w.r.t the zero proportion.
                if (zero_prop > 0.5) {
                    fres = fres3
                } else if (zero_prop > 0.3) {
                    fres = fres2
                } else {
                    fres = fres1
                }

                ## processing CD4
                if(adt_marker_index %in% trimodal_marker_index){ #if (!is.null(cd4_index) && adt_marker_index == cd4_index) {
                    fres = flowCore::filter(fcs, flowStats::curv1Filter(adt_marker_select, bwFac = bwFac_smallest))
                    peak_info = flowStats::curvPeaks(
                        x = fres,
                        dat = adt_expression,
                        borderQuant = border,
                        from = from,
                        to = to
                    )
                    if(length(peak_info$midpoint) != 3){ ## if not obtain 3 peaks, better to use a larger bw
                        fres = flowCore::filter(fcs, flowStats::curv1Filter(adt_marker_select, bwFac = bwFac_smallest + 0.5))
                        peak_info = flowStats::curvPeaks(
                        x = fres,
                        dat = adt_expression,
                        borderQuant = border,
                        from = from,
                        to = to
                        )
                    }
                    peak_ind = peak_info$peaks[, "y"] > lower_peak_thres
                    res = peak_info$midpoint[peak_ind]
                    res_region = peak_info$regions[peak_ind, ]
                # } else if(adt_marker_index %in% trimodal_marker_index){ ## trimodal marker
                #     fres = flowCore::filter(fcs, flowStats::curv1Filter(adt_marker_select, bwFac = bwFac_smallest))
                #     peak_info = flowStats::curvPeaks(
                #         x = fres,
                #         dat = adt_expression,
                #         borderQuant = border,
                #         from = from,
                #         to = to
                #     )
                #     if(length(peak_info$midpoint) != 3){ ## if not obtain 3 peaks, better to use a larger bw
                #         fres = flowCore::filter(fcs, flowStats::curv1Filter(adt_marker_select, bwFac = bwFac_smallest + 0.5))
                #         peak_info = flowStats::curvPeaks(
                #         x = fres,
                #         dat = adt_expression,
                #         borderQuant = border,
                #         from = from,
                #         to = to
                #         )
                #     }
                #     res = peak_info$midpoint
                #     res_region = peak_info$regions
                } else { ## other marker
                    peak_info = flowStats::curvPeaks(
                        x = fres,
                        dat = adt_expression,
                        borderQuant = border,
                        from = from,
                        to = to
                    )
                    ## if no peak is detected
                    if(is.na(peak_info$midpoints[1])){
                        ## try using the smallest bw
                        fres0 = flowCore::filter(fcs, flowStats::curv1Filter(adt_marker_select, bwFac = bwFac_smallest))
                        peak_info = flowStats::curvPeaks(
                            x = fres0,
                            dat = adt_expression,
                            borderQuant = 0,
                            from = from,
                            to = to
                        )
                        ## if still no peak detected
                        if(is.na(peak_info$midpoints[1])){
                            adt_expression = adt_expression + stats::rnorm(length(adt_expression), mean = 0, sd = 0.05)
                            fcs_count = adt_expression %>% t %>% t %>% as.matrix()
                            colnames(fcs_count) = adt_marker_select
                            fcs = flowCore::flowFrame(fcs_count)
                            ## update fres
                            fres1 = flowCore::filter(fcs, flowStats::curv1Filter(adt_marker_select, bwFac = 2))
                            fres2 = flowCore::filter(fcs, flowStats::curv1Filter(adt_marker_select, bwFac = 3))
                            fres3 = flowCore::filter(fcs, flowStats::curv1Filter(adt_marker_select, bwFac = 3.1))
                            fres = fres3
                            peak_info = flowStats::curvPeaks(
                                x = fres,
                                dat = adt_expression,
                                borderQuant = border,
                                from = from,
                                to = to
                            )
                        }
                    }

                    ## User defined the marker that is known to usually have multiple peaks (n = 2)
                    if (adt_marker_index %in% bimodal_marker_index) {
                        if (length(peak_info$midpoint) == 2) {
                            ## 2 peaks, perfect!
                            res = peak_info$midpoint
                            res_region = peak_info$regions
                        } else if (length(peak_info$midpoint) > 2) {
                            ## more than 2 peaks, consider filtering out very low density peaks
                            peak_ind = peak_info$peaks[, "y"] > lower_peak_thres
                            res = peak_info$midpoint[peak_ind]
                            res_region = peak_info$regions[peak_ind, ]

                        } else if (zero_prop > 0.3 && length(peak_info$midpoint) < 2) {
                            ## less than 2 peaks and zero proportion is larger than 0.3, use finer bandwidth:fres1 instead of fres2
                            peak_info = flowStats::curvPeaks(
                                x = fres1,
                                dat = adt_expression,
                                borderQuant = 0,
                                from = from,
                                to = to
                            )
                            if (length(peak_info$midpoint) == 2) {
                                ## peak number ==2 output results.
                                y_sum = peak_info$peaks[, 'y'] %>% sum
                                res = peak_info$midpoint
                                res_region = peak_info$regions

                                fres0 = flowCore::filter(fcs, flowStats::curv1Filter(adt_marker_select, bwFac = bwFac_smallest))

                                peak_info = flowStats::curvPeaks(
                                    x = fres0,
                                    dat = adt_expression,
                                    borderQuant = 0,
                                    from = from,
                                    to = to
                                )
                                if((length(peak_info$midpoint) >= 2) && (sum(peak_info$peaks[, 'y']) > y_sum) && (peak_info$midpoint[2] - peak_info$midpoint[1] > 0.3)){
                                    ## if using the smallest bw obtain better peak mode, switch from fres1 to fres0 results
                                    res = peak_info$midpoint[peak_info$peaks[, "y"] > lower_peak_thres]
                                    res_region = peak_info$regions[peak_info$peaks[, "y"] > lower_peak_thres, ]
                                }


                            } else if (length(peak_info$midpoint) > 2) {
                                ## using new bandwidth, too many peaks, consider filtering out very low density peaks
                                res = peak_info$midpoint[peak_info$peaks[, "y"] > lower_peak_thres]
                                res_region = peak_info$regions[peak_info$peaks[, "y"] > lower_peak_thres, ]
                            } else if (length(peak_info$midpoint) < 2) {
                                ## try with smallest bw to get more peak modes
                                fres0 = flowCore::filter(fcs, flowStats::curv1Filter(adt_marker_select, bwFac = bwFac_smallest))
                                peak_info = flowStats::curvPeaks(
                                    x = fres0,
                                    dat = adt_expression,
                                    borderQuant = 0,
                                    from = from,
                                    to = to
                                )
                                if(any(is.na(peak_info$midpoint)) || (length(peak_info$midpoint) >= 2 && peak_info$midpoint[2] - peak_info$midpoint[1] < 0.5)){

                                    ## smallest bw may lead to NA midpoint or peaks that are too close due to discrete value
                                    peak_info = flowStats::curvPeaks(
                                        x = fres1,
                                        dat =  adt_expression,
                                        borderQuant = 0,
                                        from = from,
                                        to = to
                                    )
                                    res = peak_info$midpoint
                                    res_region = peak_info$regions
                                }else if (length(peak_info$midpoint) >= 2) {
                                    ## using new bandwidth, too many peaks, consider filtering out very low density peaks
                                    res = peak_info$midpoint[peak_info$peaks[, "y"] > lower_peak_thres]
                                    res_region = peak_info$regions[peak_info$peaks[, "y"] > lower_peak_thres, ]
                                }else{
                                    ## still one peak left
                                    res = peak_info$midpoint
                                    res_region = peak_info$regions
                                }
                            }
                        } else if (zero_prop <= 0.3 && length(peak_info$midpoint) < 2) {
                            ## less than 2 peaks and small zero proportion, user finer bandwidth: fres0 instead of fres1
                            fres0 = flowCore::filter(fcs, flowStats::curv1Filter(adt_marker_select, bwFac = bwFac_smallest)) ## 1.5
                            peak_info = flowStats::curvPeaks(
                                x = fres0,
                                dat =  adt_expression,
                                borderQuant = 0,
                                from = from,
                                to = to
                            )
                            if(any(is.na(peak_info$midpoint)) || (length(peak_info$midpoint) >= 2 && peak_info$midpoint[2] - peak_info$midpoint[1] < 0.5)){
                                    peak_info = flowStats::curvPeaks(
                                        x = fres1,
                                        dat =  adt_expression,
                                        borderQuant = 0,
                                        from = from,
                                        to = to
                                    )
                                    res = peak_info$midpoint
                                    res_region = peak_info$regions
                            }else if (length(peak_info$midpoint) <= 2) {
                                res = peak_info$midpoint
                                res_region = peak_info$regions
                            } else if (length(peak_info$midpoint) > 2) {
                                res = peak_info$midpoint[peak_info$peaks[, "y"] > lower_peak_thres]
                                res_region = peak_info$regions[peak_info$peaks[, "y"] > lower_peak_thres, ]
                            }
                        } else {
                            ## no other cases?
                            res = peak_info$midpoint[peak_info$peaks[, "y"] > lower_peak_thres]
                            res_region = peak_info$regions[peak_info$peaks[, "y"] > lower_peak_thres, ]
                        }
                    } else {
                        ## not in user defined marker list, can have 1 peaks. Filter very low density peaks
                        res = peak_info$midpoint[peak_info$peaks[, "y"] > lower_peak_thres]
                        res_region = peak_info$regions[peak_info$peaks[, "y"] > lower_peak_thres, ]
                    }
                } ## end of other marker processing


                ## all the multiple peaks are around 0
                if (length(res) > 1 && zero_prop <= 0.3 && (sum(res < neg_candidate_thres) == length(res))) {
                    ## use broader bandwidth to merge multiple peaks around 0. Use fres2 instead fres1
                    peak_infoTmp = flowStats::curvPeaks(
                        x = fres2,
                        dat =  adt_expression,
                        borderQuant = border,
                        from = from,
                        to = to
                    )
                    # peak_infoTmp$midpoints = peak_infoTmp$peaks[, "x"]

                    if ((adt_marker_select %in% bimodal_marker_index) && (length(peak_infoTmp$midpoints) == 2)) {
                        ## if user define this marker to have 2 peaks.
                        resTmp = peak_infoTmp$midpoints
                        res_regionTmp = peak_infoTmp$regions
                    } else {
                        resTmp = peak_infoTmp$midpoints[peak_infoTmp$peaks[, "y"] > lower_peak_thres]
                        res_regionTmp = peak_infoTmp$regions[peak_infoTmp$peaks[, "y"] > lower_peak_thres, ]
                    }

                    indTmp = which(!is.na(resTmp))
                    resTmp = resTmp[indTmp]
                    if (is.null(nrow(res_regionTmp))) {
                        res_regionTmp = res_regionTmp %>%
                            as.matrix() %>%
                            t()
                    }
                    res_regionTmp = res_regionTmp[indTmp, ]
                    if (length(resTmp) > 1 && (sum(resTmp < 2) < length(resTmp))) {
                        res = resTmp
                        res_region = res_regionTmp
                    }
                }

                ## remove small negative peak around 0
                if (length(res) > 1 && zero_prop < 0.3 && (sum(res < neg_candidate_thres) < length(res))) {
                    check_index = which(res < neg_candidate_thres)
                    for(check_item in check_index){
                        if(peak_info$peaks[check_item, "y"] < 1 && peak_info$peak[max(check_index)+1, "y"]/peak_info$peak[check_item, "y"] > 5){
                            res = res[-1]
                            res_region = res_region[-1, ]
                        }
                    }
                    # if (peak_info$peaks[1, "x"] < 0.9 && peak_info$peaks[1, "y"] < 1 && peak_info$peaks[2, "x"] > 2 && peak_info$peaks[2, "y"] / peak_info$peaks[1, "y"] > 5) {
                    #     res = res[-1]
                    #     res_region = res_region[-1, ]
                    # }
                }

                ## all the peaks around 2 and zero proportion very large. Highly likely to have only one peak.
                if (length(res) > 1 && zero_prop > 0.5 && (sum(res < neg_candidate_thres) == length(res))) {
                    res = peak_info$midpoint[which(peak_info$peaks[, "y"] == max(peak_info$peak[, "y"]))]
                    res_region = peak_info$regions[which(peak_info$peaks[, "y"] == max(peak_info$peak[, "y"])), ]
                }

                ## record peak mode and peak regions
                peak_num = max(peak_num, length(res))
                peak_mode[[sample_name]] = res
                peak_region[[sample_name]] = matrix(NA, ncol = 2, nrow = length(res))
                peak_region[[sample_name]][1:length(res), ] = res_region

            }

        }else{
            peak_mode[[sample_name]] = NA
            peak_region[[sample_name]] = NA

        }

    } ## end of for loop for sample_name in sample_name_list

    ## initiate landmark to record peak mode location
    landmark = matrix(NA, ncol = peak_num, nrow = length(sample_name_list))
    landmarkRegion = list()
    for (i in seq_len(peak_num)) {
        landmarkRegion[[i]] = matrix(NA, ncol = 2, nrow = length(sample_name_list))
        rownames(landmarkRegion[[i]]) = sample_name_list
    }
    rownames(landmark) = sample_name_list
    for (i in names(peak_mode)) { ## go through samples
        if (!is.na(peak_mode[[i]][1])) { ## if the peak modes are detected
            peak_modeNum = length(peak_mode[[i]])

            if (peak_modeNum == 1) {

                ## check if the only peak should be positive peak
                pos_marker_index = which(paste0("tmpName", positive_peak$ADT_index) == adt_marker_select)
                pos_sample_index = which(positive_peak$sample == i)

                if (length(intersect(pos_marker_index, pos_sample_index)) > 0) {
                    landmark[i, min(2, peak_num)] = peak_mode[[i]]
                    landmarkRegion[[min(2, peak_num)]][i, ] = peak_region[[i]]
                } else {
                    landmark[i, 1] = peak_mode[[i]]
                    landmarkRegion[[1]][i, ] = peak_region[[i]]
                }
            } else if (peak_modeNum == 2) {
                landmark[i, c(1, max(2, peak_num))] = peak_mode[[i]]

                landmarkRegion[[1]][i, ] = peak_region[[i]][1, ]
                landmarkRegion[[max(2, peak_num)]][i, ] = peak_region[[i]][2, ]
            } else if (peak_modeNum == 3) {
                landmark[i, c(1, 2, max(3, peak_num))] = peak_mode[[i]]
                landmarkRegion[[1]][i, ] = peak_region[[i]][1, ]
                landmarkRegion[[2]][i, ] = peak_region[[i]][2, ]
                landmarkRegion[[max(3, peak_num)]][i, ] = peak_region[[i]][3, ]
            } else {
                landmark[i, 1:peak_modeNum] = peak_mode[[i]]
                for (k in 1:peak_modeNum) {
                    landmarkRegion[[k]][i, ] = peak_region[[i]][k, ]
                }
            }
        }
    }

    ## if all the peaks are within 1 - highly likely that there is only one negative peak
    if (max(landmark[!is.na(landmark)]) < neg_candidate_thres) {
        landmark_new = matrix(NA, ncol = 1, nrow = nrow(landmark))
        landmarkAllMedian = stats::median(landmark[!is.na(landmark)])
        for (i in 1:nrow(landmark)) {
            landmark_nonNA = landmark[i, !is.na(landmark[i, ])]
            if (length(landmark_nonNA) > 0) {
                landmark_new[i, ] = landmark[i, which.min(abs(landmark_nonNA - landmarkAllMedian))]
            } else {
                landmark_new[i, ] = NA
            }
        }
        landmark = landmark_new
    }
    rownames(landmark) = sample_name_list


    return(landmark)
}

