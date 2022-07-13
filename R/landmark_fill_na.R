#' Merge the peak and valley landmarks locations and fill in NA if landmark is not detected.
#'
#' This function detect the valley locations either between every two peak landmarks or cut at the right heavy tails. If specified positive uni-peak, the valley location will be set at the left side of the uni-peak. 
#' @param peak_landmark_list Matrix of peak landmark detection results. Rows are samples and column(s) are the peak locations.
#' @param valley_landmark_list Matrix of valley landmark detection results. Rows are samples and column(s) are the valley locations.
#' @param landmark_align_type Algin the peak and valleys using one of the "negPeak", "negPeak_valley", "negPeak_valley_posPeak", and "valley" alignment modes.
#' @param midpoint_type Fill in the missing first valley by the midpoint of two positive peaks ("midpoint") or impute by other valley ("valley").
#' @param neg_candidate_thres The upper bound for the negative peak. Users can refer to their IgG samples to obtain the minimal upper bound of the IgG sample peak. It can be one of the values of asinh(4/5+1), asinh(6/5+1), or asinh(8/5+1) if the right 95% quantile of IgG samples are large. 
#' @export
#' @examples
#' landmark_fill_na(peak_landmark_list, valley_landmark_list, landmark_align_type = "negPeak_valley_posPeak")
 
landmark_fill_na <-  function(peak_landmark_list = NULL, valley_landmark_list = NULL, landmark_align_type = NULL, midpoint_type = "valley", neg_candidate_thres = asinh(10/5 + 1)){
    if(!landmark_align_type %in% c("negPeak", "negPeak_valley", "negPeak_valley_posPeak", "valley")){
      return("Please provide one of the landmark_align_type from: negPeak, negPeak_valley, negPeak_valley_posPeak, valley")
    }
    if(landmark_align_type == "valley"){
      ## only use the first valley to align
      landmark_matrix <- valley_landmark_list[, 1] %>% t %>% t
      landmark_matrix[is.na(landmark_matrix), 1] <- neg_candidate_thres #2 ## fill in na by background level 2 after arcinsh_b5_a1 transformation
    }else{
      ## involve negative peaks in landmark alignment
      if(ncol(peak_landmark_list) == 1){
        ## only have negative peaks
        landmark_matrix <- cbind(
          peak_landmark_list[, 1],
          valley_landmark_list[, 1]
        )
        ## fill in na
        landmark_matrix[is.na(landmark_matrix[, 1]), 1] <- median(landmark_matrix[!is.na(landmark_matrix[, 1]), 1])
        alter_pos = median(landmark_matrix[!is.na(landmark_matrix[, 2]), 2])
        for(tmpIndex in which(is.na(landmark_matrix[, 2]))){
          if(neg_candidate_thres > landmark_matrix[tmpIndex, 1]){
            landmark_matrix[tmpIndex, 2] <- neg_candidate_thres
          }else{
            landmark_matrix[tmpIndex, 2] <- alter_pos
          }
        }
        # landmark_matrix[is.na(landmark_matrix[, 2]), 2] <- neg_candidate_thres
      }else{
        ## have positive peaks
        if(ncol(peak_landmark_list) > 2 && midpoint_type == "midpoint"){
          landmark_matrix <- cbind(
            peak_landmark_list[, 1],
            rowMeans(valley_landmark_list, na.rm = TRUE),
            peak_landmark_list[, ncol(peak_landmark_list)]
          )
        }else{
          landmark_matrix <- cbind(
            peak_landmark_list[, 1],
            valley_landmark_list[, 1],
            peak_landmark_list[, ncol(peak_landmark_list)]
          )
        }

        ## fill in na
        ## fill in valley first
        landmark_matrix[is.na(landmark_matrix[, 2]), 2] <- neg_candidate_thres
        ## due to user_define_peak where unique peak is deemed positive.
        ## fill in either 0.5 or half of the first valley
        landmark_matrix[is.na(landmark_matrix[, 1]), 1] <- landmark_matrix[is.na(landmark_matrix[, 1]), 2]/2 #0.5 
        ## fill in the last positive peak: add on the valley using the median distance from the last positive peak to the first valley
        landmark_matrix[is.na(landmark_matrix[, 3]), 3] <- landmark_matrix[is.na(landmark_matrix[, 3]), 2] + median(landmark_matrix[!is.na(landmark_matrix[, 3]), 3] - landmark_matrix[!is.na(landmark_matrix[, 3]), 2])

      }
    }

    ## only provide negative peak location
    if (landmark_align_type == "negPeak") {
      return(landmark_matrix[, 1] %>% t %>% t)
    }

    ## only provide negative peak and valley location
    if (landmark_align_type == "negPeak_valley") {
      return(landmark_matrix[, 1:2])
    }

    ## provide negative peak, first valley and last postiive peak location
    return(landmark_matrix)
}
