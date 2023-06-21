#' Extract the metadata from a study
#'
#' @param study_id The cbioportal study ID
#'
#' @return Tibble with the study metadata
#' @export
release_metadata <- function(study_id) {
  cbio <- cBioPortal("www.cbioportal.org")
  metadata <- clinicalData(cbio, study_id)
  return(metadata)
}
