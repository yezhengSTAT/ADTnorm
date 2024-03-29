% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cell_x_feature}
\alias{cell_x_feature}
\title{A matrix of raw count for the cell by features}
\format{
A data frame with 422682 rows and 7 variables:
\describe{
\item{sample}{Sample name. In this demo data, the sample name is the same as the study_name, assuming that one study is one sample.}
\item{batch}{Batch ID. In this demo data, the batch ID is the same as the study_name.}
\item{sample_status}{Sample status, i.e., Healthy, MALTtumor, HIV Vaccine, Lupus, B-ALL, AML.}
\item{study_name}{Name of the data set/study.}
\item{ADTseqDepth}{Total UMI per cell.}
\item{cell_type_l1}{Broad level of cell type annotation using manual gating.}
\item{cell_type_l2}{Fine level of cell type annotation using manual gating.}
}
}
\source{
See detailed description in the manuscript
}
\usage{
cell_x_feature
}
\description{
A dataset containing 422682 cells and 7 feature categories for the CITE-seq raw measurement of 13 publicly available CITE-seq datasets.
}
\keyword{datasets}
