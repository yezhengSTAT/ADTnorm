#' Clean ADT marker name
#'
#' This function enables the general cleaning of ADT marker names. Regardless, users should try cleaning and unifying their ADT marker name first. Please also ensure there is no "/" in the ADT name, such as "TCRγ/δ".
#' @param adt_name A vector of ADT marker name.
#' @export
#' @examples
#' \dontrun{
#' clean_adt_name(colnames(cell_x_adt))
#' }


clean_adt_name <- function(adt_name) {
    adt_rename <- adt_name %>%
        gsub("CD19-CAR", "CARCD19", .) %>%
        gsub("_PROT", "", .) %>%
        gsub("_TotalSeqB", "", .) %>%
        gsub("_control", "", .) %>%
        gsub("-GA.*", "", .) %>%
        gsub("-GC.*", "", .) %>%
        gsub("-GT.*", "", .) %>%
        gsub("-GG.*", "", .) %>%
        gsub("-CA.*", "", .) %>%
        gsub("-CC.*", "", .) %>%
        gsub("-CT.*", "", .) %>%
        gsub("-CG.*", "", .) %>%
        gsub("-AA.*", "", .) %>%
        gsub("-AC.*", "", .) %>%
        gsub("-AT.*", "", .) %>%
        gsub("-AG.*", "", .) %>%
        gsub("-TA.*", "", .) %>%
        gsub("-TC.*", "", .) %>%
        gsub("-TG.*", "", .) %>%
        gsub("-TT.*", "", .) %>%
        gsub("_.*", "", .) %>%
        gsub("ADT-", "", .)

    if (!("CD8" %in% adt_rename)) {
        if ("CD8a" %in% adt_rename | "CD8A" %in% adt_rename) {
            ind <- which(adt_rename == "CD8a" | adt_rename == "CD8A")
            adt_rename[ind] <- "CD8"
        } else if ("CD8b" %in% adt_rename | "CD8B" %in% adt_rename) {
            ind <- which(adt_rename == "CD8B" | adt_rename == "CD8b")
            adt_rename[ind] <- "CD8"
        }
    }
    adt_rename <- replace(adt_rename, adt_rename == "MouseIgG1kappaisotype", "IgG1")
    adt_rename <- replace(adt_rename, adt_rename == "MouseIgG2akappaisotype", "IgG2a")
    adt_rename <- replace(adt_rename, adt_rename == "Mouse IgG2bkIsotype", "IgG2b")
    adt_rename <- replace(adt_rename, adt_rename == "RatIgG2bkIsotype", "IgG2b-Rat")

    adt_rename <- replace(adt_rename, adt_rename == "Mouse-IgG1", "IgG1")
    adt_rename <- replace(adt_rename, adt_rename == "Mouse-IgG2a", "IgG2a")
    adt_rename <- replace(adt_rename, adt_rename == "Mouse-IgG2b", "IgG2b")
    adt_rename <- replace(adt_rename, adt_rename == "Rat-IgG2b", "IgG2b-Rat")

    adt_rename <- adt_rename %>% gsub("HLA-DR", "HLADR", .) %>% gsub("HLA\\.DR", "HLADR", .) %>% gsub("PD-1", "PD1", .) %>% gsub(" ", "", .) %>% gsub("\\s*\\([^\\)]+\\)", "", .) %>% gsub("CD279", "PD1", .) %>% gsub("CD274", "PDL1", .) %>% gsub("\\.", "", .) %>% gsub("-", "", .) %>% gsub("/", "",.)

    return(adt_rename)
}
