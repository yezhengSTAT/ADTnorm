% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peak_alignment.R
\name{peak_alignment}
\alias{peak_alignment}
\title{Align the peak and valley landmarks by the warpset function}
\usage{
peak_alignment(
  cell_x_adt,
  cell_x_feature = NULL,
  landmark_matrix = NULL,
  target_landmark = NULL
)
}
\arguments{
\item{cell_x_adt}{Matrix of ADT raw counts in cells (rows) by ADT markers (columns) format.}

\item{cell_x_feature}{Matrix of cells (rows) by cell features (columns) such as cell type, sample, and batch-related information.}

\item{landmark_matrix}{Matrix of peak and valley landmarks after filling in NA using the \code{landmark_fill_na} function.}

\item{target_landmark}{Leave it as NULL to align the landmark to the mean location across samples. Denote it by a vector of the same length as the column number of the landmark to align the negative peak, valley, and positive peak(s) to the specified fixed location.}
}
\description{
This function monotonously transforms the ADT marker counts to align the landmarks detected in previous steps. By aligning the landmarks, ADTnorm removes the batch effect and allows integration across batches/studies.
}
\examples{
\dontrun{
peak_alignment(cell_x_adt, cell_x_feature, landmark_matrix)
}
}
