#' Align the peak and valley landmarks by the warpset function
#'
#' This function monotonously transforms the ADT marker counts to align the landmarks detected in previous steps. By aligning the landmarks, ADTnorm removes the batch effect and allows integration across batches/studies.
#' @param cell_x_adt Matrix of ADT raw counts in cells (rows) by ADT markers (columns) format.
#' @param cell_x_feature Matrix of cells (rows) by cell features (columns) such as cell type, sample, and batch-related information.
#' @param landmark_matrix Matrix of peak and valley landmarks after filling in NA using the `landmark_fill_na` function.
#' @param target_landmark Leave it as NULL to align the landmark to the mean location across samples. Denote it by a vector of the same length as the column number of the landmark to align the negative peak, valley, and positive peak(s) to the specified fixed location.
#' @export
#' @examples
#' \dontrun{
#' peak_alignment(cell_x_adt, cell_x_feature, landmark_matrix)
#' }
# require(dplyr)
# require(flowStats)
# require(fda)
peak_alignment = function(cell_x_adt, cell_x_feature = NULL, landmark_matrix = NULL, target_landmark = NULL, neg_candidate_thres = asinh(8/5 + 1)) {
  ## get parameters
  grouping = NULL
  monwrd = TRUE
  subsample = NULL
  peakNr = NULL
  clipRange = 0.01
  nbreaks = 11 #11
  bwFac = 2
  warpFuns = FALSE
  chunksinze = 10
  newNcFile = NULL
  z = NULL
  nb = 1001

  exp_data = cell_x_adt
  cell_x_adt_norm = cell_x_adt
  samples = levels(cell_x_feature$sample) ## sampleNames(exp_data)

  ## set up fda parameters
  extend = 0.15
  from = min(c(min(cell_x_adt, na.rm = TRUE),  target_landmark[1], min(landmark_matrix))) - diff(range(cell_x_adt, na.rm = TRUE)) * extend
  to = max(c(max(cell_x_adt, na.rm = TRUE), target_landmark[length(target_landmark)], max(landmark_matrix))) + diff(range(cell_x_adt, na.rm = TRUE)) * extend
  
  lower_bound = min(cell_x_adt, na.rm = TRUE) - diff(range(cell_x_adt, na.rm = TRUE)) * extend
  upper_bound = max(cell_x_adt, na.rm = TRUE) + diff(range(cell_x_adt, na.rm = TRUE)) * extend
  
  wbasis = fda::create.bspline.basis(
    rangeval = c(from, to),
    norder = 4, breaks = seq(from, to, len = nbreaks)
  )
  Wfd0 = fda::fd(matrix(0, wbasis$nbasis, 1), wbasis)
  WfdPar = fda::fdPar(Wfd0, 1, 1e-5)


  density_y = c()
  for(sample in samples){
    cell_ind_tmp = which(cell_x_feature$sample == sample)
    cell_ind = cell_ind_tmp[which(!is.na(cell_x_adt[cell_ind_tmp]))]
    if(length(cell_ind) > 0){
      density_y = cbind(density_y, stats::density(cell_x_adt[cell_ind], from = from, to = to, n = nb, na.rm = TRUE)$y)
    }else{
      density_y = cbind(density_y, rep(NA, nb))
    }

  }
  colnames(density_y) = samples

  arg_vals = seq(from, to, len = nb)
  fdobj = fda::Data2fd(arg_vals, density_y, wbasis)

  if (ncol(landmark_matrix) == 1) { ## only one peak no valley: offset
    offsets = landmark_matrix - stats::median(landmark_matrix, na.rm = TRUE)
    names(offsets) = samples
    funs = funsBack = vector("list", nrow(landmark_matrix))
    names(funs) = samples
    names(funsBack) = samples
    for (j in seq_along(funs)) {
      funs[[samples[j]]] = function(x) x - z
      e1 = new.env(hash = TRUE)
      e1$z = offsets[samples[j]]
      environment(funs[[samples[j]]]) = e1
      funsBack[[samples[j]]] = function(x) x + z
      e2 = new.env(hash = TRUE)
      e2$z = offsets[samples[j]]
      environment(funsBack[[samples[j]]]) = e2
    }
  } else { ## more than one landmark: warping
    ## if any valley is beyond the upper bound of range, replace by the upper bound
    if(any(landmark_matrix[, 2] > upper_bound)){
      landmark_matrix[which(landmark_matrix[, 2] > upper_bound), 2] = max(cell_x_adt, na.rm = TRUE) 
      print(paste0("Warning: some valley landmarks are larger the upper bound of the range. They are replaced by the maximum value of cell_x_adt. Please consider reduce 'neg_candidate_thres' value. The default value for 'neg_candidate_thres' is asinh(8/5 + 1) and the current value for 'neg_candidate_thres' is ", neg_candidate_thres, "."))
    }
    if(any(landmark_matrix[, 1] < lower_bound)){
      landmark_matrix[which(landmark_matrix[, 1] < lower_bound), 1] = min(cell_x_adt, na.rm = TRUE) 
      print("Warning: some valley landmarks are smaller than the lower bound of the range. They are replaced by the minimum value of cell_x_adt.")
    }
    args = list("unregfd" = fdobj, "fdobj"=fdobj, "ximarks"=landmark_matrix, "WfdPar"=WfdPar, "monwrd"=monwrd)
    if(!is.null(target_landmark)){
      args[['x0marks']] = target_landmark
    }else{
      args[['x0marks']] = colMeans(landmark_matrix, na.rm = TRUE)
    }
    args_run = args[intersect(names(formals(fda::landmarkreg)), names(args))]
    regDens = do.call(fda::landmarkreg, args_run, quote = TRUE)

    # if(is.null(target_landmark)){
    #   regDens = fda::landmarkreg(fdobj, landmark_matrix, WfdPar = WfdPar, monwrd = monwrd)
    # }else{
    #   regDens = fda::landmarkreg(fdobj, landmark_matrix, x0marks = target_landmark, WfdPar = WfdPar, monwrd = monwrd)
    # }
    # if(is.null(target_landmark)){
    #   regDens = fda::landmarkreg(fdobj, landmark_matrix, WfdPar = WfdPar)
    # }else{
    #   regDens = fda::landmarkreg(fdobj, landmark_matrix, x0marks = target_landmark, WfdPar = WfdPar)
    # }
    warpfdobj = regDens$warpfd
    warpedX = fda::eval.fd(warpfdobj, arg_vals)
    warpedX[1, ] = utils::head(arg_vals, 1)
    warpedX[nrow(warpedX), ] = utils::tail(arg_vals, 1)
    funs = apply(warpedX, 2, stats::approxfun, arg_vals)
    funsBack = apply(warpedX, 2, function(a, b) stats::approxfun(b, a), arg_vals)
  }


  names(funs) = names(funsBack) = samples

  warped_landmark_matrix = landmark_matrix
  leftBoard = rightBoard = vector("list", length(funs))
  newRange = c(Inf, -Inf)

  ## transform the raw data using the warping functions
  for (i in seq_along(funs)) {
    # cell_index = which(cell_x_feature$sample == samples[i])
    cell_ind_tmp = which(cell_x_feature$sample == samples[i])
    cell_index = cell_ind_tmp[which(!is.na(cell_x_adt[cell_ind_tmp]))]
    if(length(cell_index) > 0){
      thisDat = t(t(cell_x_adt[cell_index]))
      newDat = as.matrix(funs[[i]](thisDat))
      newDat[is.na(newDat)] = thisDat[is.na(newDat)]
      cell_x_adt_norm[cell_index] = newDat
      warped_landmark_matrix[i, ] = funs[[i]](landmark_matrix[i, ])
    }else{
      warped_landmark_matrix[i, ] = NA
    }

  }
  return(list(cell_x_adt_norm = cell_x_adt_norm, landmark_matrix_norm = warped_landmark_matrix))
}
