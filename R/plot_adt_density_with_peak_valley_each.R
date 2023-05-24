#' Plot the expression density profile for ONE ADT marker with identifies peak and valley locations
#'
#' This function plots adt expression density profile with identifies peak and valley locations for only one ADT marker. Each track is a sample. Color by batch
#' @param adt_count Matrix of ADT raw counts in cells (rows) by one target ADT marker (column) format.
#' @param cell_x_feature Matrix of cells (rows) by cell features (columns) such as cell type, sample, and batch related information.
#' @param peak_landmark_list Matrix of peak landmark locations with rows being samples and columns being the peaks.
#' @param valley_landmark_list Matrix of valley landmark locations with rows being samples and columns being the valleys.
#' @param brewer_palettes Set the color scheme of color brewer.
#' @param parameter_list Users can specify: "run_label" to give name for this run; "bw" to adjust the band width of the density plot.
#' @export
#' @examples
#' \dontrun{
#' plot_adt_density_with_peak_valley_each(
#'   cell_x_adt,
#'   cell_x_feature,
#'   peak_landmark_list = peak_mode_norm_res,
#'   valley_landmark_list = valley_location_norm_res,
#'   brewer_palettes = "Set1",
#'   parameter_list = list(bw = 0.1, run_label = "ADTnorm")
#' )
#' }
# require(ggplot2)
# require(RColorBrewer)
# require(tidyr)
# require(ggridges)
# require(ggpubr)
plot_adt_density_with_peak_valley_each = function(adt_count, cell_x_feature, peak_landmark_list, valley_landmark_list, brewer_palettes = "Set1", parameter_list = NULL) {
    if (is.null(parameter_list)) {
        return("parameter_list is NULL!")
    }
    parameter_list_name = names(parameter_list)

    run_label = ""
    bw = 1
    if (!is.null(parameter_list)) {
        if ("run_label" %in% parameter_list_name) {
            run_label = parameter_list[["run_label"]]
        }
        if("bw" %in% parameter_list_name){
            bw = parameter_list$bw
        }
    }

    # peak_landmark_list = parameter_list$peak_landmark_list
    # valley_landmark_list = parameter_list$valley_landmark_list
    # brewer_palettes = parameter_list$brewer_palettes
  
    # If there is no batch, add a dummy variable
    if (! "batch" %in% colnames(cell_x_feature)){ cell_x_feature$batch <- 1 }
      
    tmpProfile = data.frame(counts = adt_count) %>%
        mutate(
            sample = rep(cell_x_feature$sample, 1),
            batch = rep(cell_x_feature$batch, 1)
        ) #%>% dplyr::filter(!is.na(counts))

    peak_location = list()
    valley_location = list()
    for (i in 1:ncol(peak_landmark_list)) {
        peak_location[[i]] = data.frame(
            # ADT = adt_marker_select,
            sample = cell_x_feature$sample %>% levels(),
            peakx = peak_landmark_list[, i],
            peaky = 0.5,
            peaks = 1:length(levels(cell_x_feature$sample))
        )
        if (i <= ncol(valley_landmark_list)) {
            valley_location[[i]] = data.frame(
                # ADT = adt_marker_select,
                sample = cell_x_feature$sample %>% levels(),
                peakx = valley_landmark_list[, i],
                peaky = 0.5,
                peaks = 1:length(levels(cell_x_feature$sample))
            )
        }
    }
    fillColor = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, brewer_palettes))(length(unique(tmpProfile$batch)))
    resPlot = ggplot(tmpProfile, aes_string(x = "counts", y = "sample")) +
        ggridges::geom_density_ridges(aes(fill = factor(batch)), bandwidth = bw) +
        geom_segment(data = peak_location[[1]], aes(x = peakx, xend = peakx, y = peaks, yend = peaky + peaks), size = 1) +
        geom_segment(data = valley_location[[1]], aes(x = peakx, xend = peakx, y = peaks, yend = peaky + peaks), size = 1, color = "grey") +
        # facet_wrap(~ factor(ADT), scales = "free_x") +
        theme_bw(base_size = 20) +
        xlab(run_label) +
        ylab("") +
        ggpubr::rotate_x_text(angle = 90) +
        ggpubr::rremove("legend") +
        scale_fill_manual(values = fillColor) +
        ggpubr::rremove("legend.title")


    if (ncol(peak_landmark_list) >= 2) {
        resPlot = resPlot +
            geom_segment(data = peak_location[[2]], aes(x = peakx, xend = peakx, y = peaks, yend = peaky + peaks), size = 1)
    }

    if (ncol(peak_landmark_list) >= 3) {
        resPlot = resPlot +
            geom_segment(data = peak_location[[3]], aes(x = peakx, xend = peakx, y = peaks, yend = peaky + peaks), size = 1) +
            geom_segment(data = valley_location[[2]], aes(x = peakx, xend = peakx, y = peaks, yend = peaky + peaks), size = 1, color = "grey")
    }
    if (ncol(peak_landmark_list) >= 4) {
        resPlot = resPlot +
            geom_segment(data = peak_location[[4]], aes(x = peakx, xend = peakx, y = peaks, yend = peaky + peaks), size = 1) +
            geom_segment(data = valley_location[[3]], aes(x = peakx, xend = peakx, y = peaks, yend = peaky + peaks), size = 1, color = "grey")
    }

    return(resPlot)
}
