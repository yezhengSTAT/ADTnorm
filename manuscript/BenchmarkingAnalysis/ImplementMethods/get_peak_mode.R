## get peak location
require(flowStats)
require(dplyr)

get_peak_mode <- function(cell_x_adt = NULL, cell_x_feature = NULL, adt_marker_select = NULL, parameter_list = NULL) {

    ## get parameters
    user_define_marker <- NULL ## known marker that tend to have bimodal
    user_define_peak <- NULL ## known peak that tend to be positive if only one peak detected
    bwFac <- 1.2
    border <- 0.01
    peakNr <- NULL
    densities <- NULL
    n <- 201
    indices <- FALSE
    lower_peak_thres <- 0.001

    start_adt_method <- parameter_list$start_adt_method
    fcs_path <- parameter_list$fcs_path
    parameter_list_name <- names(parameter_list)

    if ("user_define_marker" %in% parameter_list_name) {
        user_define_marker <- parameter_list$user_define_marker ## known marker that tend to have bimodal
    }
    if ("user_define_peak" %in% parameter_list_name) {
        user_define_peak <- parameter_list$user_define_peak ## known peak that tend to be positive if only one peak detected
    }
    if ("bwFac_smallest" %in% parameter_list_name) {
        bwFac_smallest <- parameter_list$bwFac_smallest
    } else {
        bwFac_smallest <- 1.5
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

    dat <- read.ncdfFlowSet(paste0(fcs_path, "/", start_adt_method, "/", levels(cell_x_feature$sample), ".fcs"))

    ranges <- fsApply(dat, range)
    from <- min(sapply(ranges, function(z) z[1, adt_marker_select] - diff(z[, adt_marker_select]) * 0.15), na.rm = TRUE)
    to <- max(sapply(ranges, function(z) z[2, adt_marker_select] + diff(z[, adt_marker_select]) * 0.15), na.rm = TRUE)

    peak_num <- 0
    peak_mode <- list()
    peak_region <- list()
    zero_prop_list <- list()

    for (sample_name in sampleNames(dat)) {
        try_out <- tryCatch(filter(dat[sample_name], curv1Filter(adt_marker_select, bwFac = 2)), error = function(e) {
            c(1, 2)
        })
        if (length(try_out) == 1) {
            fres1 <- filter(dat[sample_name], curv1Filter(adt_marker_select, bwFac = 2))
            fres2 <- filter(dat[sample_name], curv1Filter(adt_marker_select, bwFac = 3))
            fres3 <- filter(dat[sample_name], curv1Filter(adt_marker_select, bwFac = 3.1))

            cell_ind <- which(cell_x_feature$sample == gsub(".fcs", "", sample_name))
            zero_prop <- sum(cell_x_adt[cell_ind, adt_marker_select] < 2) / length(cell_x_adt[cell_ind, adt_marker_select])
            zero_prop_list[[gsub(".fcs", "", sample_name)]] <- zero_prop

            ## different bandwidth w.r.t the zero proportion.
            if (zero_prop > 0.5) {
                fres <- fres3
            } else if (zero_prop > 0.3) {
                fres <- fres2
            } else {
                fres <- fres1
            }

            if (adt_marker_select == "CD4") {
                fres <- filter(dat[sample_name], curv1Filter(adt_marker_select, bwFac = bwFac_smallest))
            }
            peak_info <- flowStats:::curvPeaks(
                x = fres[[sample_name]],
                dat = exprs(dat[[sample_name]])[, adt_marker_select],
                borderQuant = border,
                from = from,
                to = to
            )
            # peak_info$midpoints = peak_info$peaks[, "x"]

            ## User defined the marker that is known to usually have multiple peaks (n = 2)
            if (adt_marker_select %in% user_define_marker) {
                if (length(peak_info$midpoints) == 2) {
                    ## 2 peaks, perfect!
                    res <- peak_info$midpoints
                    res_region <- peak_info$regions
                } else if (length(peak_info$midpoints) > 2) {
                    ## more than 2 peaks, consider filtering out very low density peaks
                    peak_ind <- peak_info$peaks[, "y"] > lower_peak_thres
                    res <- peak_info$midpoints[peak_ind]
                    res_region <- peak_info$regions[peak_ind, ]
                } else if (zero_prop > 0.3 && length(peak_info$midpoints) < 2) {
                    ## less than 2 peaks and zero proportion is larger than 0.3, use finer bandwidth:fres1 instead of fres2
                    peak_info <- flowStats:::curvPeaks(
                        x = fres1[[sample_name]],
                        dat = exprs(dat[[sample_name]])[, adt_marker_select],
                        borderQuant = 0,
                        from = from,
                        to = to
                    )
                    # peak_info$midpoints = peak_info$peaks[, "x"]
                    if (length(peak_info$midpoints) <= 2) {
                        ## peak number <=2 output results.
                        res <- peak_info$midpoints
                        res_region <- peak_info$regions
                    } else if (length(peak_info$midpoints > 2)) {
                        ## using new bandwidth, too many peaks, consider filtering out very low density peaks
                        res <- peak_info$midpoints[peak_info$peaks[, "y"] > lower_peak_thres]
                        res_region <- peak_info$regions[peak_info$peaks[, "y"] > lower_peak_thres, ]
                    }
                } else if (zero_prop <= 0.3 && length(peak_info$midpoints) < 2) {
                    ## less than 2 peaks and small zero proportion, user finer bandwidth: fres0 instead of fres1
                    fres0 <- filter(dat[sample_name], curv1Filter(adt_marker_select, bwFac = bwFac_smallest)) ## 1.5
                    peak_info <- flowStats:::curvPeaks(
                        x = fres0[[sample_name]],
                        dat = exprs(dat[[sample_name]])[, adt_marker_select],
                        borderQuant = 0,
                        from = from,
                        to = to
                    )
                    # peak_info$midpoints = peak_info$peaks[, "x"]
                    if (length(peak_info$midpoints) <= 2) {
                        res <- peak_info$midpoints
                        res_region <- peak_info$regions
                    } else if (length(peak_info$midpoints > 2)) {
                        res <- peak_info$midpoints[peak_info$peaks[, "y"] > lower_peak_thres]
                        res_region <- peak_info$regions[peak_info$peaks[, "y"] > lower_peak_thres, ]
                    }
                } else {
                    ## no other cases?
                    res <- peak_info$midpoints[peak_info$peaks[, "y"] > lower_peak_thres]
                    res_region <- peak_info$regions[peak_info$peaks[, "y"] > lower_peak_thres, ]
                }
            } else {
                ## not in user defined marker list, can have 1 peaks. Filter very low density peaks
                res <- peak_info$midpoints[peak_info$peaks[, "y"] > lower_peak_thres]
                res_region <- peak_info$regions[peak_info$peaks[, "y"] > lower_peak_thres, ]
            }


            ## all the multiple peaks are around 0
            if (length(res) > 1 && zero_prop <= 0.3 && (sum(res < 2) == length(res))) {
                ## use broader bandwidth to merge multiple peaks around 0. Use fres2 instead fres1
                peak_infoTmp <- flowStats:::curvPeaks(
                    x = fres2[[sample_name]],
                    dat = exprs(dat[[sample_name]])[, adt_marker_select],
                    borderQuant = border,
                    from = from,
                    to = to
                )
                # peak_infoTmp$midpoints = peak_infoTmp$peaks[, "x"]

                if (adt_marker_select %in% user_define_marker) {
                    ## if user define this marker to have 2 peaks.
                    if (length(peak_infoTmp$midpoints) == 2) {
                        resTmp <- peak_infoTmp$midpoints
                        res_regionTmp <- peak_infoTmp$regions
                    } else {
                        resTmp <- peak_infoTmp$midpoints[peak_infoTmp$peaks[, "y"] > lower_peak_thres]
                        res_regionTmp <- peak_infoTmp$regions[peak_infoTmp$peaks[, "y"] > lower_peak_thres, ]
                    }
                } else {
                    resTmp <- peak_infoTmp$midpoints[peak_infoTmp$peaks[, "y"] > lower_peak_thres]
                    res_regionTmp <- peak_infoTmp$regions[peak_infoTmp$peaks[, "y"] > lower_peak_thres, ]
                }

                indTmp <- which(!is.na(resTmp))
                resTmp <- resTmp[indTmp]
                if (is.null(nrow(res_regionTmp))) {
                    res_regionTmp <- res_regionTmp %>%
                        as.matrix() %>%
                        t()
                }
                res_regionTmp <- res_regionTmp[indTmp, ]
                if (length(resTmp) > 1 && (sum(resTmp < 2) < length(resTmp))) {
                    res <- resTmp
                    res_region <- res_regionTmp
                }
            }

            ## remove small negative peak around 0
            if (length(res) > 1 && zero_prop < 0.3 && (sum(res < 2) < length(res))) {
                if (peak_info$peaks[1, "x"] < 0.9 && peak_info$peaks[1, "y"] < 1 && peak_info$peaks[2, "x"] > 2 && peak_info$peaks[2, "y"] / peak_info$peaks[1, "y"] > 5) {
                    res <- res[-1]
                    res_region <- res_region[-1, ]
                }
            }

            ## all the peaks around 2 and zero proportion very large. Highly likely to have only one peak.
            if (length(res) > 1 && zero_prop > 0.5 && (sum(res < 2) == length(res))) {
                res <- peak_info$midpoints[which(peak_info$peaks[, "y"] == max(peak_info$peak[, "y"]))]
                res_region <- peak_info$regions[which(peak_info$peaks[, "y"] == max(peak_info$peak[, "y"])), ]
            }



            peak_num <- max(peak_num, length(res))
            peak_mode[[sample_name]] <- res
            peak_region[[sample_name]] <- matrix(NA, ncol = 2, nrow = length(res))
            peak_region[[sample_name]][1:length(res), ] <- res_region
        } else {
            ## only one value for this marker
            peak_mode[[sample_name]] <- NA
            peak_region[[sample_name]] <- matrix(NA, ncol = 2, nrow = 1)
            print(paste0(i, "-Single Value!"))
        }
    }

    landmark <- matrix(NA, ncol = peak_num, nrow = length(peak_mode))
    landmarkRegion <- list()
    for (i in 1:peak_num) {
        landmarkRegion[[i]] <- matrix(NA, ncol = 2, nrow = length(peak_mode))
    }
    for (i in 1:length(peak_mode)) {
        if (!is.na(peak_mode[[i]][1])) {
            peak_modeNum <- length(peak_mode[[i]])
            if (peak_modeNum == 1) {
                if (paste0(names(peak_mode)[i], ">", adt_marker_select) %in% user_define_peak) {
                    landmark[i, min(2, peak_num)] <- peak_mode[[i]]
                    landmarkRegion[[min(2, peak_num)]][i, ] <- peak_region[[i]]
                } else {
                    landmark[i, 1] <- peak_mode[[i]]
                    landmarkRegion[[1]][i, ] <- peak_region[[i]]
                }
            } else if (peak_modeNum == 2) {
                landmark[i, c(1, max(2, peak_num))] <- peak_mode[[i]]

                landmarkRegion[[1]][i, ] <- peak_region[[i]][1, ]
                landmarkRegion[[max(2, peak_num)]][i, ] <- peak_region[[i]][2, ]
            } else if (peak_modeNum == 3) {
                landmark[i, c(1, 2, max(3, peak_num))] <- peak_mode[[i]]
                landmarkRegion[[1]][i, ] <- peak_region[[i]][1, ]
                landmarkRegion[[2]][i, ] <- peak_region[[i]][2, ]
                landmarkRegion[[max(3, peak_num)]][i, ] <- peak_region[[i]][3, ]
            } else {
                landmark[i, 1:peak_modeNum] <- peak_mode[[i]]
                for (k in 1:peak_modeNum) {
                    landmarkRegion[[k]][i, ] <- peak_region[[i]][k, ]
                }
            }
        }
    }

    ## if all the peaks are within 1 - highly likely that there is only one negative peak
    if (max(landmark[!is.na(landmark)]) < 2) {
        landmark_new <- matrix(NA, ncol = 1, nrow = nrow(landmark))
        landmarkAllMedian <- median(landmark[!is.na(landmark)])
        for (i in 1:nrow(landmark)) {
            landmark_nonNA <- landmark[i, !is.na(landmark[i, ])]
            if (length(landmark_nonNA) > 0) {
                landmark_new[i, ] <- landmark[i, which.min(abs(landmark_nonNA - landmarkAllMedian))]
            } else {
                landmark_new[i, ] <- NA
            }
        }
        landmark <- landmark_new
    }

    rownames(landmark) <- levels(cell_x_feature$sample)
    return(list(peak_landmark_list = landmark, zero_prop_list = zero_prop_list))
}
