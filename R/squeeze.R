#' Extract studies from cbioportal with the specifcied type of data and minimum number of samples
#'
#' @param n Minimum number of samples
#' @param ... Data types
#'
#' @return Tibble with all studies
#' @export
squeeze <- function(n, ...){

  all_studies <-read_csv("../inst/extdata/all_studies_info.csv")

  # Find studies with specified data
  all_studies_info <-
    all_studies |>
    dplyr::select(studyId, name, cancerTypeId, referenceGenome, importDate,...)

  # Filter according to minimum sample size
  candidate_studies <-
    all_studies_info |>
    filter(if_all(6:last_col(), ~ . > n))

  return(candidate_studies)

}
