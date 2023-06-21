#' Extract the metadata from a study
#'
#' @param study_id The cbioportal study ID
#'
#' @return Tibble with the study metadata
#' @export
#'
#' @examples
#' release_metadata(study_id = "brca_cptac_2020")
release_metadata <- function(study_id) {
  cbio <- cBioPortal()
  metadata <- clinicalData(cbio, study_id)
  return(metadata)
}
