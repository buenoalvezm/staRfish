#' Extract the molecular data from a study
#'
#' @param study_id The cbioportal study ID
#' @param molecular_type The type of data to extract
#'
#' @return Tibble with the selected molecular data
#' @export
release_data <- function(study_id, molecular_type) {

  cbio <- cBioPortal("www.cbioportal.org")

  # Extract all protein-coding gene ids
  all_genes <-
    geneTable(cbio, pageSize = 90000) |>
    filter(type == "protein-coding")

  profiles_list <- molecularProfiles(api = cbio, studyId = study_id)[["molecularProfileId"]]

  if(any(grepl("v2",profiles_list)) &  startsWith(molecular_type, prefix = "mrna")){
    molecular_profile_id <- paste0(study_id,"_rna_seq_v2_mrna")
  } else if (startsWith(molecular_type, prefix = "mrna")){
    molecular_profile_id <- paste0(study_id,"_rna_seq_mrna")
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

  final_data <-
    molecular_data[[1]] |>
    dplyr::select(gene=hugoGeneSymbol,pID=patientId,value=value) |>
    pivot_wider(id_cols=gene,names_from=pID,values_from=value)


  return(final_data)
}
