
#' Grab studies from cbioportal with the specifcied type of data
#'
#' @param n Minimum number of samples required
#' @param ...
#'
#' @return Table of candidate studies
#' @export
grab <- function(){

  # Find all available studies
  all_studies <-
    cbioportalR::available_studies("www.cbioportal.org")

  # Get info for all studies
  all_studies_info <-
    map_df(all_studies$studyId, function(study) {
      get_study_info(study_id = study, base_url = "www.cbioportal.org") |>
        as_tibble()
    })

  # Save results
  write_csv(all_studies_info, "inst/extdata/all_studies_info.csv")

}







