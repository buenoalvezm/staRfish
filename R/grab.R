
#' Grab studies from cbioportal with the specifcied type of data
#'
#' @param n Minimum number of samples required
#' @param ...
#'
#' @return Table of candidate studies
#' @export
#'
#' @examples
#' grab(n = 20, massSpectrometrySampleCount, mrnaRnaSeqSampleCount)
grab <- function(n, ...){

  # Retrieve all studies
  all_studies <-
    cbioportalR::available_studies("www.cbioportal.org")

  # Find studies with specified data
  all_studies_info <-
    map_df(all_studies$studyId, function(study) {
      get_study_info(study) |>
        as_tibble() |>
        select(studyId, name, cancerTypeId, referenceGenome, importDate,...)
    })

  # Find studies with both transcriptomics and proteomics data
  candidate_studies <-
    all_studies_info |>
    filter(if_all(6:last_col(), ~ . > n))

  return(candidate_studies)

}


