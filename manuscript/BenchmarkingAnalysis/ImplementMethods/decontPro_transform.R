## decontPro transformation
require(decontX)
require(Seurat)
require(dplyr)

decontPro_transform <- function(cell_x_adt = NULL, cell_x_feature = NULL, parameter_list = NULL){

  if(!is.null(parameter_list)){
    delta_sd = parameter_list$delta_sd,
    background_sd = parameter_list$background_sd
  }else{
    ## give a default value
    delta_sd = 2e-5
    background_sd = 2e-6
  }
  adt_seurat <- CreateSeuratObject(t(cell_x_adt), assay = "ADT")
  adt_seurat <- NormalizeData(adt_seurat, normalization.method = "CLR", margin = 2) %>%
  ScaleData(assay = "ADT") %>%
  RunPCA(assay = "ADT", features = rownames(adt_seurat), npcs = min(ncol(cell_x_adt), 10), reduction.name = "pca_adt") 
  
  adt_seurat = FindNeighbors(adt_seurat, dims = 1:min(ncol(cell_x_adt)-1, 10), assay = "ADT", reduction = "pca_adt")  %>%
  FindClusters(resolution = 0.5)

  clusters <- as.integer(Idents(adt_seurat))
  
  rm_cell_index = which(rowSums(cell_x_adt) == 0)
  out <- decontPro(t(cell_x_adt[-rm_cell_index, ]),
                 clusters[-rm_cell_index],
                 delta_sd = delta_sd,
                 background_sd = background_sd)

  decontaminated_counts <- out$decontaminated_counts
  norm_counts <- matrix(0, nrow(decontaminated_counts), nrow(cell_x_adt))
  norm_counts[, -rm_cell_index] <- decontaminated_counts

  return(t(norm_counts))
}
