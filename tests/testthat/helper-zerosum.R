library(ADTnorm)

make_test_data_zerosum <- function(n_cells = 100) {
    withr::with_seed(42, {
        markers <- c("CD3", "CD4", "ZERO", "CD19", "CD56")
        counts <- suppressWarnings(matrix(
            rnbinom(n_cells * length(markers), size = 10, mu = 50),
            nrow = n_cells, ncol = length(markers)
        ))
        colnames(counts) <- markers
        counts[, "ZERO"] <- 0L
        list(
            cell_x_adt = as.data.frame(counts),
            cell_x_feature = data.frame(
                cell = paste0("cell", seq_len(n_cells)),
                sample = rep(c("S1", "S2"), each = n_cells / 2)
            ),
            markers = markers
        )
    })
}

make_test_data_allzero <- function(n_cells = 100) {
    withr::with_seed(42, {
        markers <- c("CD3", "CD4", "CD19")
        counts <- matrix(
            0L, nrow = n_cells, ncol = length(markers),
            dimnames = list(NULL, markers)
        )
        list(
            cell_x_adt = as.data.frame(counts),
            cell_x_feature = data.frame(
                cell = paste0("cell", seq_len(n_cells)),
                sample = rep(c("S1", "S2"), each = n_cells / 2)
            )
        )
    })
}
