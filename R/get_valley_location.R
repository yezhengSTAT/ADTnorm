#' Get the valley landmark locations.
#'
#' This function detect the valley locations either between every two peak landmarks or cut at the right heavy tails. If specified positive uni-peak, the valley location will be set at the left side of the uni-peak.
#' @param cell_x_adt Matrix of ADT raw counts in cells (rows) by ADT markers (columns) format.
#' @param cell_x_feature Matrix of cells (rows) by cell features (columns) such as cell type, sample, and batch related information.
#' @param adt_marker_select Markers to normalize. Leaving empty to process all the ADT markers in cell_x_adt matrix.
#' @param peak_mode_res The peak landmark results coming out of `get_peak_mode` or `get_peak_midpoint` function.
#' @param shoulder_valley Indictor to specify whether a shoulder valley is expected in case of heavy right tail where the population of cells should be considered as positive population.
#' @param positive_peak A list variable containing a vector of ADT marker(s) and a corresponding vector of sample name(s) in matching order to specify that the uni-peak detected should be aligned to positive peaks. For example, for samples that only contain T cells. The only CD3 peak should be aligned to positive peaks of other samples.
#' @param multi_sample_per_batch Set it to TRUE to discard the positive peak that only appear in one sample per batch (sample number is >=3 per batch).
#' @param adjust Parameter for `density` function: bandwidth used is actually adjust*bw. This makes it easy to specify values like ‘half the default’ bandwidth.
#' @param min_fc Mimimal fold change between the highest peak density height and candidate valley density height. Default is 20.
#' @param shoulder_valley_slope The slope on the ADT marker density distribution to call shoulder valley.
#' @param neg_candidate_thres The upper bound for the negative peak. Users can refer to their IgG samples to obtain the minimal upper bound of the IgG sample peak. It can be one of the values of asinh(4/5+1), asinh(6/5+1), or asinh(8/5+1) if the right 95% quantile of IgG samples are large.
#' @param lower_peak_thres The minimal ADT marker density height to call it a real peak. Set it to 0.01 to avoid suspecious positive peak. Set it to 0.001 or smaller to include some small but tend to be real positive peaks, especially for markers like CD19.
#' @param arcsine_transform_flag The flag indicating if input is raw count and arcsine transformation is implemented.
#' @examples
#' \dontrun{
#' get_valley_location(cell_x_adt, cell_x_feature, peak_mode_res)
#' }
#' @export
#' @importFrom magrittr %$%
get_valley_location = function(cell_x_adt = NULL, cell_x_feature = NULL, adt_marker_select = NULL, peak_mode_res = NULL, shoulder_valley = TRUE, positive_peak = NULL, multi_sample_per_batch = FALSE, adjust = 1.5, min_fc = 20, shoulder_valley_slope = -1, lower_peak_thres = 0.01, neg_candidate_thres = asinh(10/5 + 1), arcsine_transform_flag = TRUE) {

    peak_landmark_list = peak_mode_res

    ## tag sample if it is the only sample within each batch that have different number of peaks.
    # batch_num = cell_x_feature$batch %>% unique %>% length
    tag_row = vector("list", length = nrow(peak_landmark_list))
    names(tag_row) = rownames(peak_landmark_list)

    if(multi_sample_per_batch && ncol(peak_landmark_list) > 1){ ## multiple sample per batch remove the positive peak that only appear in one sample
        for(batch_each in unique(cell_x_feature$batch)){
            sample_index = cell_x_feature %>% dplyr::filter(batch == batch_each) %$% sample %>% unique
            ## if there is only one second peak, remove it
            if(length(sample_index) > 2){
                col_index = which(colSums(!is.na(peak_landmark_list[sample_index, ])) == 1) ## which col has only 1 peak per batch
                # common_peak_col = which(colSums(!is.na(peak_landmark_list[sample_index, ])) != 1) #ncol(peak_landmark_list[sample_index, ]) - length(col_index)
                for(c in col_index){ ## across those columns that only have 1 peak per batch
                    row_index = sample_index[which(!is.na(peak_landmark_list[sample_index, c]))] ## get the sample name
                    tag_row[[row_index]] = c(tag_row[[row_index]], c) #common_peak_n #non-zero value, record the peak number that is common across samples in this batch
                }
            }
        }
    }

    valley_location_list = matrix(NA, nrow = nrow(peak_landmark_list), ncol = max(1, ncol(peak_landmark_list) - 1))
    rownames(valley_location_list) = cell_x_feature$sample %>% levels()
    for (sample_name in cell_x_feature$sample %>% levels()) {
        # cell_ind = which(cell_x_feature$sample == sample_name)
        cell_ind_tmp = which(cell_x_feature$sample == sample_name)
        cell_notNA = which(!is.na(cell_x_adt[cell_ind_tmp, adt_marker_select]))
        cell_ind = cell_ind_tmp[cell_notNA]
        if(length(cell_ind) > 0){
            peak_landmark = peak_landmark_list[sample_name, ]
            if(arcsine_transform_flag){
                zero_prop = sum(cell_x_adt[cell_ind, adt_marker_select] < 2) / length(cell_x_adt[cell_ind, adt_marker_select])
            }else{
                zero_prop = sum(cell_x_adt[cell_ind, adt_marker_select] < neg_candidate_thres) / length(cell_x_adt[cell_ind, adt_marker_select])
            }

            ## check if user define single peak to be positive peak
            pos_marker_index = which(paste0("tmpName", positive_peak$ADT_index) == adt_marker_select)
            pos_sample_index = which(positive_peak$sample == sample_name)
            if (length(intersect(pos_marker_index, pos_sample_index)) > 0) {
                lower_valley = TRUE
            } else {
                lower_valley = FALSE
            }


            density_res = stats::density(
                cell_x_adt[which(cell_x_feature$sample == sample_name), adt_marker_select],
                adjust = adjust, na.rm = TRUE
            )
            x = density_res$x
            y = density_res$y

            sign_diff = sign(diff(y))
            diff_sign_diff = diff(sign_diff)
            peak = which(diff_sign_diff == -2) + 1
            valley = which(diff_sign_diff == 2) + 1

            x_valley = x[valley]
            y_valley = y[valley]
            x_peak = x[peak]
            real_peak = peak_landmark[!is.na(peak_landmark)] # peak
            np = length(real_peak)

            ## if this sample is tagged, choose the highest peak
            if (!is.null(tag_row[[sample_name]])) {## there are peaks need to be removed
                peak_landmark_y = c()
                for (peak_landmark_each in peak_landmark) {
                    peak_landmark_y = c(peak_landmark_y, y[which.min(abs(x - peak_landmark_each))]) ## get the density value for the peak location
                }
                common_peak_n = length(peak_landmark_y) - length(tag_row[[sample_name]])
                real_peak = real_peak[sort(order(peak_landmark_y, decreasing = T)[1:common_peak_n])] ## get the highest peak
                np = length(real_peak) ## update number of peak i.e, np
                peak_landmark_list[sample_name, ] = NA
                peak_landmark_list[sample_name, (1:ncol(peak_landmark_list))[-tag_row[[sample_name]]]] = real_peak ## remove the outlier peak
            }


            if (np > 1) { ## two peaks or more
                real_valley = c()
                for (i in 1:(np - 1)) {
                    tmp_valley = x_valley[(x_valley > real_peak[i]) & (x_valley < real_peak[i + 1])]
                    tmp_real_valley = tmp_valley[which.min(y_valley[(x_valley > real_peak[i]) & (x_valley < real_peak[i + 1])])]
                    if(length(tmp_real_valley) == 0){ ## no minimal point between two peak, ensure tmp_real_valley is not empty
                        tmp_real_valley = (real_peak[i] + real_peak[i+1]) / 2
                    }
                    if(i == 1 && shoulder_valley){
                        shoulder_cand_index = which(diff(y)/diff(x) > shoulder_valley_slope)
                        first_peak_index = (which(x > max(x_peak[1], real_peak[1])) %>% min) + 50
                        x_shoulder = x[shoulder_cand_index[shoulder_cand_index > first_peak_index][1]]
                        real_valley = c(real_valley, min(x_shoulder, tmp_real_valley, na.rm = T))
                    }else{
                        if(length(tmp_valley) == 0){ ## there is no local minimum that fall within two peaks
                            real_valley = c(real_valley, (real_peak[i] + real_peak[i+1]) / 2) ## use the midpoint of two peak location --- maybe shoulder peak instead?
                        }else{
                            real_valley = c(real_valley, tmp_real_valley) #c(real_valley, tmp_valley[which.min(y_valley[(x_valley > real_peak[i]) & (x_valley < real_peak[i + 1])])]) ## local minimum based on the density
                        }

                    }
                }
            } else if(np == 1) { ## one peak
                if (lower_valley == FALSE) { ## one peak is negative peak
                    real_valley =  x_valley[x_valley > real_peak[1] + 0.1][1]  #x[which(y < max(y) / min_fc)[which(y < max(y) / min_fc) > max(which(y == max(y)), which(x > real_peak[1]) %>% min())] %>% min()]
                    if(shoulder_valley){
                        ## if one peak & go for shoulder threshold
                        shoulder_cand_index = which(diff(y)/diff(x) > shoulder_valley_slope)
                        first_peak_index = (which(x > max(x_peak[1], real_peak[1])) %>% min) + 50
                        x_shoulder = x[shoulder_cand_index[shoulder_cand_index > first_peak_index][1]]
                        real_valley = min(x_shoulder, real_valley, na.rm = T)

                    }else{
                        ## check if no valley is detected due to shoulder peak
                        # if(length(y_valley[x_valley > real_peak[1]]) == 0 || (y_valley[x_valley >  real_peak[1] + 0.1][1] < 0.05)){
                        #     ## if one peak consider the shoulder point
                        #     shoulder_cand_index = which(diff(y)/diff(x) > shoulder_valley_slope)
                        #     first_peak_index = (which(x > max(x_peak[1], real_peak[1])) %>% min) + 50
                        #     x_shoulder = x[shoulder_cand_index[shoulder_cand_index > first_peak_index][1]]
                        #     real_valley = min(x_shoulder, real_valley, na.rm = T)
                        # }
                        if (length(y_valley[!is.na(x_valley) & x_valley > real_peak[1]]) == 0 || (y_valley[!is.na(x_valley) & x_valley > real_peak[1] + 0.1][1] < 0.05)) {
                            shoulder_cand_index = which(diff(y)/diff(x) > shoulder_valley_slope)
                            first_peak_index = (which(x > max(x_peak[1], real_peak[1])) %>% min) + 50
                            x_shoulder = x[shoulder_cand_index[shoulder_cand_index > first_peak_index][1]]
                            real_valley = min(x_shoulder, real_valley, na.rm = TRUE)
                        }
                    }
                    if (zero_prop > 0.8) {
                        real_valley = max(neg_candidate_thres, real_valley)
                    }
                } else { ## one peak is positive peak
                    real_valley = x[which((y < max(y) / min_fc) | (y < lower_peak_thres))[which((y < max(y) / min_fc) | (y < lower_peak_thres))< min(which(y == max(y)), which(x < real_peak[1]) %>% max())] %>% max()]
                    # peak_landmark_list[sample_name, ] = asinh(0/5 + 1)
                    # peak_landmark_list[sample_name, ncol(peak_landmark_list)] = real_peak[1]
                    # real_valley = min(real_valley, min(cell_x_adt[cell_ind, adt_marker_select])) #min(real_valley, asinh(0/5 + 1))
                }
            } else { ## no peak
                real_valley = NA
            }
            ## check if no valley is detected due to shoulder peak
            if (length(real_valley) == 0 || all(is.na(real_valley))) {
                if(length(real_peak) >= 2){ ## midpoint of two peak
                    real_valley = (real_peak[1] + real_peak[2]) / 2
                }else{## shoulder valley
                    shoulder_cand_index = which(diff(y)/diff(x) > shoulder_valley_slope)
                    first_peak_index = (which(x > max(x_peak[1], real_peak[1])) %>% min) + 50
                    x_shoulder = x[shoulder_cand_index[shoulder_cand_index > first_peak_index][1]]
                    if(!is.na(x_shoulder)){
                        real_valley = min(x_shoulder, real_valley, na.rm = T)
                    }
                }

            }
            valley_location_list[sample_name, 1:length(real_valley)] = real_valley
        }


    }

    ## update peak_landmark_list if needed
    if(any(colSums(!is.na(peak_landmark_list)) == 0)){ ## remove column that only has NA values
        col_index = which(colSums(!is.na(peak_landmark_list)) == 0)
        peak_landmark_update = matrix(NA, nrow(peak_landmark_list), ncol(peak_landmark_list) - length(col_index))
        i = 1
        for (c in 1:ncol(peak_landmark_list)) {
            if (!(c %in% col_index)) {
            peak_landmark_update[, i] = peak_landmark_list[, c]
            i = i + 1
            }
        }
        rownames(peak_landmark_update) = rownames(peak_landmark_list)
        peak_landmark_list = peak_landmark_update
    }
    return(list(valley_location_list = valley_location_list, peak_landmark_list = peak_landmark_list))
}
