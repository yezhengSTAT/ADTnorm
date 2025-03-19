## warpset peak valley alignment
require(dplyr)
require(flowStats)
warpset_transform <- function(landmark_matrix = NULL, cell_x_feature = NULL, adt_marker_select = NULL, parameter_list = NULL) {
  ## get parameters
  grouping <- NULL
  monwrd <- TRUE
  subsample <- NULL
  peakNr <- NULL
  clipRange <- 0.01
  nbreaks <- 11
  bwFac <- 2
  warpFuns <- FALSE
  chunksinze <- 10
  newNcFile <- NULL
  z <- NULL
  nb <- 1001
  eps <- .Machine$double.eps
  target <- NULL


  start_adt_method <- parameter_list$start_adt_method
  fcs_path <- parameter_list$fcs_path
  
  ## read in original adt expression matrix
  exp_data <- read.ncdfFlowSet(paste0(fcs_path, "/", start_adt_method, "/", levels(cell_x_feature$sample), ".fcs")) %>% as.flowSet()

  ## check if target marker is available
  adt_marker_flag <- adt_marker_select %in% colnames(exp_data)
  if (!all(adt_marker_flag)) {
    stop(
      "Invalid stain(s) not matching the flowSet:\n    ",
      paste(adt_marker_select[!adt_marker_flag], collapse = ", ")
    )
  }

  ## expression range and sample name list
  ranges <- fsApply(exp_data, range)
  samples <- sampleNames(exp_data)

  ## set up fda parameters
  extend <- 0.15
  from <- min(sapply(ranges, function(z) z[1, adt_marker_select] - diff(z[, adt_marker_select]) * extend), na.rm = TRUE)
  to <- max(sapply(ranges, function(z) z[2, adt_marker_select] + diff(z[, adt_marker_select]) * extend), na.rm = TRUE)
  wbasis <- create.bspline.basis(
    rangeval = c(from, to),
    norder = 4, breaks = seq(from, to, len = nbreaks)
  )
  Wfd0 <- fd(matrix(0, wbasis$nbasis, 1), wbasis)
  WfdPar <- fdPar(Wfd0, 1, 1e-4)

  density_y <- t(fsApply(exp_data, function(exp_data_each) { # t(fsApply(thisX, function(z){
    exp_data_each_range <- range(exp_data_each)[, adt_marker_select]
    exp_data_each <- exprs(exp_data_each)
    exp_data_each <- exp_data_each[exp_data_each[, adt_marker_select] > exp_data_each_range[1] + eps & exp_data_each[, adt_marker_select] < exp_data_each_range[2] - eps, adt_marker_select]
    density(exp_data_each, from = from, to = to, n = nb, na.rm = TRUE)$y
  }))
  arg_vals <- seq(from, to, len = nb)
  fdobj <- Data2fd(arg_vals, density_y, wbasis)

  if (ncol(landmark_matrix) == 1) { ## only one peak: offset
    if (is.null(target)) {
      offsets <- landmark_matrix - median(landmark_matrix)
      names(offsets) <- sampleNames(exp_data)
    } else {
      offsets <- landmark_matrix - landmark_matrix[sampleNames(exp_data) %in% target]
      names(offsets) <- sampleNames(exp_data)
    }
    funs <- funsBack <- vector("list", length(landmark_matrix))
    names(funs) <- samples
    names(funsBack) <- samples
    for (j in seq_along(funs)) {
      funs[[samples[[j]]]] <- function(x) x - z
      e1 <- new.env(hash = TRUE)
      e1$z <- offsets[samples[[j]]]
      environment(funs[[samples[[j]]]]) <- e1
      funsBack[[samples[[j]]]] <- function(x) x + z
      e2 <- new.env(hash = TRUE)
      e2$z <- offsets[samples[[j]]]
      environment(funsBack[[samples[[j]]]]) <- e2
    }
  } else { ## multiple peaks: warping
    if (is.null(target)) {
      capture.output(regDens <- landmarkreg(fdobj, landmark_matrix, WfdPar = WfdPar, monwrd = monwrd))
    } else {
      capture.output(regDens <- landmarkreg(fdobj, landmark_matrix, x0marks = apply(landmark_matrix, 2, jitter)[rownames(landmark_matrix) %in% target, ], WfdPar = WfdPar, monwrd = monwrd))
    }
    warpfdobj <- regDens$warpfd
    warpedX <- eval.fd(warpfdobj, arg_vals)
    warpedX[1, ] <- head(arg_vals, 1)
    warpedX[nrow(warpedX), ] <- tail(arg_vals, 1)
    ## compute warping functions
    ## funs <-  apply(warpedX, 2, function(y) approxfun(arg_vals, y))
    funs <- apply(warpedX, 2, approxfun, arg_vals)
    funsBack <- apply(warpedX, 2, function(a, b) approxfun(b, a), arg_vals)
  }


  names(funs) <- names(funsBack) <- samples # sampleNames(thisX)

  warped_landmark_matrix <- landmark_matrix
  leftBoard <- rightBoard <- vector("list", length(funs))
  # chunkleftBoard<-chunkrightBoard<-rep(list(length(funs)),max(1:length(funs)%/%chunksize)+1)
  newRange <- c(Inf, -Inf)

  ## transform the raw data using the warping functions
  for (i in seq_along(funs)) {
    thisDat <- exprs(exp_data[[i]][, adt_marker_select])
    lb <- thisDat < ranges[[i]][1, adt_marker_select] + eps
    lb[is.na(lb)] <- TRUE
    leftBoard[[i]] <- lb
    rb <- thisDat > ranges[[i]][2, adt_marker_select] - eps
    rb[is.na(rb)] <- TRUE
    rightBoard[[i]] <- rb
    # Include ALL data, none of this thresholding crap at borders.
    sel <- leftBoard[[i]] | rightBoard[[i]]
    # sel<-rep(FALSE,length(thisDat))
    # browser();
    newDat <- as.matrix(funs[[i]](thisDat[!sel, ]))
    newDat[is.na(newDat)] <- thisDat[!sel, ][is.na(newDat)]
    exprs(exp_data[[i]])[!sel, adt_marker_select] <- newDat
    warped_landmark_matrix[i, ] <- funs[[i]](landmark_matrix[i, ])
    newRange[1] <- min(newRange[1], min(exprs(exp_data[[i]])[, adt_marker_select], na.rm = TRUE))
    newRange[2] <- max(newRange[2], max(exprs(exp_data[[i]])[, adt_marker_select], na.rm = TRUE))
  }
  ## make sure that edge envents are set to the extreme values
  ## of the warped data range and update the parameters slot
  ## accordingly
  for (i in seq_along(funs)) {
    minSel <- leftBoard[[i]]
    maxSel <- rightBoard[[i]]
    exprs(exp_data[[i]])[minSel, adt_marker_select] <- as.matrix(rep(
      newRange[1],
      sum(minSel, na.rm = TRUE)
    ),
    ncol = 1
    )
    exprs(exp_data[[i]])[maxSel, adt_marker_select] <- as.matrix(rep(
      newRange[2],
      sum(maxSel, na.rm = TRUE)
    ),
    ncol = 1
    )
    ip <- match(adt_marker_select, pData(parameters(exp_data[[i]]))$name)
    tmp <- parameters(exp_data[[i]])
    oldRanges <- unlist(range(exp_data[[i]])[, adt_marker_select])
    pData(tmp)[ip, c("minRange", "maxRange")] <- c(
      min(oldRanges[1], newRange[1]),
      max(oldRanges[2], newRange[2])
    )
    exp_data[[i]]@parameters <- tmp
  }


  # exp_data
  exp_data_aligned <- as(exp_data, "flowSet")
  phenoData(exp_data_aligned) <- phenoData(exp_data)
  exp_data_aligned <- exp_data_aligned[sampleNames(exp_data)]

  ## convert into matrix
  cell_x_adt_aligned <- c()
  for(sample_name in paste0(levels(cell_x_feature$sample), ".fcs")){
    cell_x_adt_aligned <- rbind(cell_x_adt_aligned, exp_data_aligned[[sample_name]] %>% exprs)
  }
  return(cell_x_adt_aligned %>% data.frame)
}
