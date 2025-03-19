## get valley location
require(dplyr)

get_valley_location <- function(cell_x_adt = NULL, cell_x_feature = NULL, adt_marker_select = NULL, parameter_list = NULL, ....) {
    adjust <- parameter_list$adjust
    peak_landmark_list <- parameter_list$peak_mode_res$peak_landmark_list
    zero_prop_list <- parameter_list$peak_mode_res$zero_prop_list
    user_define_peak <- parameter_list$user_define_peak

    parameter_list_name <- names(parameter_list)

    if ("min_fc" %in% parameter_list_name) {
        min_fc <- parameter_list$min_fc
    }
    valley_location_list <- matrix(NA, nrow = nrow(peak_landmark_list), ncol = max(1, ncol(peak_landmark_list) - 1))
    rownames(valley_location_list) <- cell_x_feature$sample %>% levels()
    for (sample_name in cell_x_feature$sample %>% levels()) {
        peak_landmark <- peak_landmark_list[sample_name, ]
        zero_prop <- zero_prop_list[[sample_name]]

        ## check if user define single peak to be positive peak
        if (paste0(sample_name, ">", adt_marker_select) %in% user_define_peak) {
            lower_valley <- TRUE
        } else {
            lower_valley <- FALSE
        }

        density_res <- density(
            cell_x_adt[which(cell_x_feature$sample == sample_name), adt_marker_select],
            adjust = adjust
        )
        x <- density_res$x
        y <- density_res$y

        sign_diff <- sign(diff(y))
        diff_sign_diff <- diff(sign_diff)
        # peak <- which(diff_sign_diff == -2) + 1
        valley <- which(diff_sign_diff == 2) + 1

        # yp <- y[peak]
        x_valley <- x[valley]
        y_valley <- y[valley]
        #   stdy <- std(yp)
        real_peak <- peak_landmark[!is.na(peak_landmark)] # peak
        np <- length(real_peak)


        if (np > 1) {
            real_valley <- c()
            for (i in 1:(np - 1)) {
                tmp_valley <- x_valley[(x_valley > real_peak[i]) & (x_valley < real_peak[i + 1])]
                real_valley <- c(real_valley, tmp_valley[which.min(y_valley[(x_valley > real_peak[i]) & (x_valley < real_peak[i + 1])])])
            }
        } else {
            if (lower_valley == FALSE) {
                real_valley <- x[which(y < max(y) / min_fc)[which(y < max(y) / min_fc) > max(which(y == max(y)), which(x > real_peak[1]) %>% min())] %>% min()] # x_valley[x_valley > real_peak[1]][1]
                if (zero_prop > 0.8) {
                    real_valley <- max(2, real_valley)
                }
            } else {
                real_valley <- x[which(y < max(y) / min_fc)[which(y < max(y) / min_fc) < min(which(y == max(y)), which(x < real_peak[1]) %>% max())] %>% max()]
            }
        }
        ## check if no valley is detected due to shoulder peak
        if (length(real_valley) == 0) {
            real_valley <- (real_peak[1] + real_peak[2]) / 2
        }

        valley_location_list[sample_name, 1:length(real_valley)] <- real_valley
    }


    ## rownames(valley_location_list) <- cell_x_feature$sample %>% levels()
    return(valley_location_list)
}
