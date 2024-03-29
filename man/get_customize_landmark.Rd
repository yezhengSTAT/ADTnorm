% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_customize_landmark.R
\name{get_customize_landmark}
\alias{get_customize_landmark}
\title{Prompt Shiny browser to manually customize peaks and valleys locations.}
\usage{
get_customize_landmark(
  cell_x_adt_sample,
  landmark_pos,
  bw,
  adt_marker_select_name,
  brewer_palettes = "Set1"
)
}
\arguments{
\item{cell_x_adt_sample}{Matrix of ADT counts of the selected marker, with columns of sample and batch information for each row of cells.}

\item{landmark_pos}{Matrix of landmark location including peaks and valleys.}

\item{bw}{Bandwidth for the density plot.}

\item{adt_marker_select_name}{The ADT marker needed to be manually processed to set the landmarks.}

\item{brewer_palettes}{Set the color scheme of the color brewer. The default is "Set1".}
}
\description{
This function will launch a shiny app allowing the user to set the location of peaks and valleys manually. The function will output the landmark positions that the user has set.
}
\examples{
\dontrun{
get_customize_landmark(
cell_x_adt_sample, 
landmark_pos, 
bw, 
adt_marker_select_name, 
brewer_palettes = "Set1"
)}
}
