#' Align the peak and valley landmarks by the warpset function
#'
#' This function detect the valley locations either between every two peak landmarks or cut at the right heavy tails. If specified positive uni-peak, the valley location will be set at the left side of the uni-peak.
#' @param cell_x_adt Matrix of ADT raw counts in cells (rows) by ADT markers (columns) format.
#' @param cell_x_feature Matrix of cells (rows) by cell features (columns) such as cell type, sample, and batch related information.
#' @param landmark_matrix Matrix of peak and valley landmarks after fill in NA using `landmark_fill_na` function.
#' @param target_landmark Leave it as NULL to align the landmark to the mean location across samples. Denote it by a vector of the same length of the column number of landmark to align the negative peak, valley and positive peak(s) to the specified fixed location.
#' @export
#' @examples
#' \dontrun{
#' peak_alignment(cell_x_adt, cell_x_feature, landmark_matrix)
#' }
# require(dplyr)
# require(flowStats)
# require(fda)
peak_alignment = function(cell_x_adt, cell_x_feature = NULL, landmark_matrix = NULL, target_landmark = NULL) {
  ## get parameters
  grouping = NULL
  monwrd = TRUE
  subsample = NULL
  peakNr = NULL
  clipRange = 0.01
  nbreaks = 11
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
  from = min(cell_x_adt, na.rm = TRUE) - diff(range(cell_x_adt, na.rm = TRUE)) * extend
  to = max(cell_x_adt, na.rm = TRUE) + diff(range(cell_x_adt, na.rm = TRUE)) * extend
  wbasis = create.bspline.basis(
    rangeval = c(from, to),
    norder = 4, breaks = seq(from, to, len = nbreaks)
  )
  Wfd0 = fd(matrix(0, wbasis$nbasis, 1), wbasis)
  WfdPar = fdPar(Wfd0, 1, 1e-4)


  density_y = c()
  for(sample in samples){
    cell_ind_tmp = which(cell_x_feature$sample == sample)
    cell_ind = cell_ind_tmp[which(!is.na(cell_x_adt[cell_ind_tmp]))]
    if(length(cell_ind) > 0){
      density_y = cbind(density_y, density(cell_x_adt[cell_ind], from = from, to = to, n = nb, na.rm = TRUE)$y)
    }else{
      density_y = cbind(density_y, rep(NA, nb))
    }

  }
  colnames(density_y) = samples

  arg_vals = seq(from, to, len = nb)
  fdobj = Data2fd(arg_vals, density_y, wbasis)

  if (ncol(landmark_matrix) == 1) { ## only one peak no valley: offset
    offsets = landmark_matrix - median(landmark_matrix)
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
    warpedX = eval.fd(warpfdobj, arg_vals)
    warpedX[1, ] = head(arg_vals, 1)
    warpedX[nrow(warpedX), ] = tail(arg_vals, 1)
    funs = apply(warpedX, 2, approxfun, arg_vals)
    funsBack = apply(warpedX, 2, function(a, b) approxfun(b, a), arg_vals)
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
