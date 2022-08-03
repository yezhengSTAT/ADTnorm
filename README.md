# ADTnorm

<!-- badges: start -->
[![R-CMD-check](https://github.com/yezhengSTAT/ADTnorm/workflows/R-CMD-check/badge.svg)](https://github.com/yezhengSTAT/ADTnorm/actions)
<!-- badges: end -->

## What is `ADTnorm`?

CITE-seq technology enables the direct measurement of protein expression, known as antibody-derived tags (ADT), in addition to RNA expression. The increase in the copy number of protein molecules leads to a more robust detection of protein features compared to RNA, providing a deep definition of cell types. However, due to added discrepancies of antibodies, such as the different types or concentrations of IgG antibodies, the batch effects of the ADT component of CITE-seq can dominate over biological variations, especially for the across-study integration. We present ADTnorm as a normalization and integration method designed explicitly for the ADT counts of CITE-seq data. Benchmarking with existing scaling and normalization methods, ADTnorm achieves a fast and accurate matching of the negative and positive peaks of the ADT counts across samples, efficiently removing technical variations across batches. Further quantitative evaluations confirm that ADTnorm achieves the best cell-type separation while maintaining the minimal batch effect. Therefore, ADTnorm facilitates the scalable ADT count integration of massive public CITE-seq datasets with distinguished experimental designs, which are essential for creating a corpus of well-annotated single-cell data with deep and standardized annotations.

Manuscript: [Zheng et al. Robust Normalization and Integration of Single-cell Protein Expression across CITE-seq Datasets. BioRxiv. 2022](https://www.biorxiv.org/content/10.1101/2022.04.29.489989v1)


## ADT Normaliztion Pipeline

<img src="./man/figures/pipeline_202208.png" alt="ADTnorm" width="600px">


## Installation

``` R
# install.packages("remotes")
remotes::install_github("yezhengSTAT/ADTnorm", build_vignettes = FALSE)
```

## Input Data

Demo data sets are ```cell_x_adt``` and ```cell_x_feature``` where ```cell_x_adt``` contains a matrix of raw count for the cell by ADT markers. 

- ```cell_x_adt``` is a data frame with 422682 cells (row) and 9 ADT markers (column): CD3, CD4, CD8, CD14, CD19, CD25, CD45RA, CD56, CD127.

- ```cell_x_feature``` is a data frame with 422682 cells (row) and 7  feature variables (column):
  
    - sample: Sample name.
    
    - batch: Batch ID. In this demo data, the batch ID is the same as the study_name.
    
    - sample_status: Sample status, i.e., Healthy, MALTtumor, HIV Vaccine, Lupus, B-ALL, AML.

    - study_name: Name of the data set/study.

    - ADTseqDepth: Total UMI per cell.
    
    - cell_type_l1: Broad level of cell type annotation using manual gating.
    
    - cell_type_l2: Fine level of cell type annotation using manual gating.


``` R
data(cell_x_adt)
data(cell_x_feature) 
```

## Usage

``` R
library(ADTnorm)

save_outpath <- "/path/to/output/location"
study_name <- "ADTnorm_demoRun"

cell_x_adt_norm <- ADTnorm(
  cell_x_adt = cell_x_adt, 
  cell_x_feature = cell_x_feature,
  save_outpath = save_outpath, 
  study_name = study_name, 
  marker_to_process = c("CD3", "CD4", "CD8", "CD45RA"), 
  trimodal_marker = c("CD4", "CD45RA"), 
  positive_peak = list(ADT = "CD3", sample = "buus_2021_T"),
  save_intermediate = TRUE
)
```

Basic parameters introduction. The full parameter explanation for the ```ADTnorm``` function can be found at [Reference - ADTnorm](yezhengstat.github.io/ADTnorm/reference/adtnorm.html). 

```
cell_x_adt:         Matrix of ADT raw counts in cells (rows) by ADT markers (columns) format.

cell_x_feature:     Matrix of cells (rows) by cell features (columns) such as cell type, sample, and batch related information.

save_outpath:       The path to save the results.

study_name:         Name of this run.

marker_to_process:  Markers to normalize. Leaving empty to process all the ADT markers in cell_x_adt matrix.

bimodal_marker:     Specify ADT markers that tend to have two peaks based on researchers' prior knowledge or preliminary observation on particular data to be processed.

trimodal_marker:    Index of the ADT markers that tend to have three peaks based on researchers' prior knowledge (e.g. CD4) or preliminary observation on particular data to be processed.

positive_peak:      A list variable containing a vector of ADT marker(s) and a corresponding vector of sample name(s) in matching order to specify that the uni-peak detected should be aligned to positive peaks. For example, for samples that only contain T cells. The only CD3 peak should be aligned to positive peaks of other samples.

save_intermediate:  Save the rds file for the intermediate objects and generate the density plot figure for checking the peak and valley location detection.
```

For more detailed and typical parameter tuning examples, please visit [tutorial website](yezhengstat.github.io/ADTnorm/articles/ADTnorm-tutorial.html). We will illustrate using the demo data.


## Results

```ADTnorm``` function will generate a matrix of rows of the same number as input ```cell_x_adt``` row number and columns are ADT markers specificed in ```marker_to_process```. The value in the matrix is normalized value by ADTnorm.  In the `save_outpath` specified by the users, there will be two subfolders, `figures` and `RDS`, containing the intermediate object and density plot of detected peak and valley landmarks before and after ADTnorm. Those figures can be used to check if further parameter tuning is needed for certain ADT markers.  

### Raw Counts 

<img src="./man/figures/RawCount.png" alt="RawCount" width="700px">

### ADTnorm Counts

<img src="./man/figures/ADTnorm.png" alt="Normalization" width="700px">


## Contact

Email: yzheng23@fredhutch.org
