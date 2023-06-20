
#' Grab studies from cbioportal with the specifcied type of data
#'
#' @param n Minimum number of samples required
#' @param ...
#'
#' @return Table of candidate studies
#' @export
#'
#' @examples
#' grab()
grab <- function(n, ...){

  # Find all available studies
  all_studies <-
    cbioportalR::available_studies("www.cbioportal.org")

  # Get info for all studies
  all_studies_info <-
    map_df(all_studies$studyId, function(study) {
      get_study_info(study_id = study, base_url = "wwww.cbioportal.org") |>
        as_tibble()
    })

  # Save results
  write_csv(all_studies_info, "data/all_studies_info.csv")

}







