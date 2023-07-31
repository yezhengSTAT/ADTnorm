#' Load saved landmarks into a list of matrices, as required for input into "override_landmark".
#'
#' @param dir Directory where landmark location .rds files are stored
#' @param append_rds Whether to append "/RDS" to the provided directory path. The default is TRUE, which aligns with ADTnorm's default save location.
#' @export
#' @examples
#' \dontrun{
#' override_landmark = load_landmarks(dir = save_outpath)
#' }

load_landmarks = function(dir, append_rds = TRUE){
    if(append_rds)
        dir <- paste(dir,'RDS',sep='/')
    
    override_landmark <- list()
    for(file in list.files(dir,pattern='^peak_valley_locations_.*.rds',full.names=TRUE)){
        file_items <- tail(unlist(strsplit(file,'/')),1)
        file_items <- unlist(strsplit(file_items,'_'))
        marker <- paste(file_items[4:(length(file_items)-2)],sep='_')
        override_landmark[[marker]] <- readRDS(file)
    }
    return(override_landmark)
}