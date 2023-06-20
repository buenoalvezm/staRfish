#' Title Extract the molecular data from a study
#'
#' @param study_id The cbioportal study ID
#' @param molecular_type The type of data to extract
#'
#' @return Tibble with the selected molecular data
#' @export
#'
#' @examples
#' release_data(study_id = "coad_silu_2022", molecular_type = "mrnaRnaSeqSampleCount")
#' release_data(study_id = "brca_cptac_2020", molecular_type = "mrnaRnaSeqV2SampleCount")
#' release_data(study_id = "brca_cptac_2020", molecular_type = "massSpectrometrySampleCount")
release_data <- function(study_id, molecular_type) {

  cbio <- cBioPortal()

  # Extract all protein-coding gene ids
  all_genes <-
    geneTable(cbio, pageSize = 90000) |>
    filter(type == "protein-coding")

  # Define molecular_type id
  if(molecular_type == "mrnaRnaSeqSampleCount") {
    molecular_profile_id <- paste0(study_id, "_rna_seq_mrna")
  } else if (molecular_type == "mrnaRnaSeqV2SampleCount") {
    molecular_profile_id <- paste0(study_id, "_rna_seq_v2_mrna")
  } else if (molecular_type == "massSpectrometrySampleCount") {
    molecular_profile_id <- paste0(study_id, "_protein_quantification")
  } else {
    print(paste("Molecular type", molecular_type, "not allowed. Supported molecular types are mrnaRnaSeqSampleCount, mrnaRnaSeqV2SampleCount and massSpectrometrySampleCount."))
  }

  # Extract data for all genes
  molecular_data <-
    getDataByGenes(
      cbio,
      studyId = study_id,
      genes = all_genes$entrezGeneId,
      by =  c("entrezGeneId", "hugoGeneSymbol"),
      molecularProfileId = molecular_profile_id)

  return(molecular_data[[1]])
}
