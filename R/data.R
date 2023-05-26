#' A matrix of raw count for the cell by ADT markers
#'
#' A dataset containing 422682 cells and 9 ADT markers for the CITE-seq raw measurement of 13 publicly available CITE-seq datasets.
#'
#' @format A data frame with 422682 rows and 9 variables:
#' \describe{
#'   \item{CD3}{CD3 ADT marker raw count across each cell}
#'   \item{CD4}{CD4 ADT marker raw count across each cell}
#'   \item{CD8}{CD8 ADT marker raw count across each cell}
#'   \item{CD14}{CD14 ADT marker raw count across each cell}
#'   \item{CD19}{CD19 ADT marker raw count across each cell}
#'   \item{CD25}{CD25 ADT marker raw count across each cell}
#'   \item{CD45RA}{CD45RA ADT marker raw count across each cell}
#'   \item{CD56}{CD56 ADT marker raw count across each cell}
#'   \item{CD127}{CD127 ADT marker raw count across each cell}
#' }
#' @source See detailed description in the manuscript
"cell_x_adt"

#' A matrix of raw count for the cell by features
#'
#' A dataset containing 422682 cells and 7 feature categories for the CITE-seq raw measurement of 13 publicly available CITE-seq datasets.
#'
#' @format A data frame with 422682 rows and 7 variables:
#' \describe{
#'   \item{sample}{Sample name. In this demo data, the sample name is the same as the study_name, assuming that one study is one sample.}
#'   \item{batch}{Batch ID. In this demo data, the batch ID is the same as the study_name.}
#'   \item{sample_status}{Sample status, i.e., Healthy, MALTtumor, HIV Vaccine, Lupus, B-ALL, AML.}
#'   \item{study_name}{Name of the data set/study.}
#'   \item{ADTseqDepth}{Total UMI per cell.}
#'   \item{cell_type_l1}{Broad level of cell type annotation using manual gating.}
#'   \item{cell_type_l2}{Fine level of cell type annotation using manual gating.}
#' }
#' @source See detailed description in the manuscript
"cell_x_feature"