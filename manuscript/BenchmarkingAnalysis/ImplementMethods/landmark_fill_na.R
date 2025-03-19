## merge peak and valley and fill in NA

landmark_fill_na <-  function(peak_landmark_list = NULL, valley_landmark_list = NULL, landmark_align_type = NULL){
    if(!landmark_align_type %in% c("negPeak", "negPeak_valley", "negPeak_valley_posPeak", "valley")){
      return("Please provide one of the landmark_align_type from: negPeak, negPeak_valley, negPeak_valley_posPeak, valley")
    }
    if(landmark_align_type == "valley"){
      ## only use the first valley to align
      landmark_matrix <- valley_landmark_list[, 1] %>% t %>% t
      landmark_matrix[is.na(landmark_matrix), 1] <- 2 ## fill in na by background level 2 after arcinsh_b5_a1 transformation
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
        landmark_matrix[is.na(landmark_matrix[, 2]), 2] <- 2
      }else{
        ## have positive peaks
        landmark_matrix <- cbind(
          peak_landmark_list[, 1],
          valley_landmark_list[, 1],
          peak_landmark_list[, ncol(peak_landmark_list)]
        )

        ## fill in na
        ## fill in valley first
        landmark_matrix[is.na(landmark_matrix[, 2]), 2] <- 2
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
