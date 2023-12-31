% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/digest.R
\name{create_rna_prot_correlation}
\alias{create_rna_prot_correlation}
\title{create_rna_prot_correlation
Creates a table with correlation between RNA and protein levels for each gene found in both the RNA and protein data.}
\usage{
create_rna_prot_correlation(rna, protein)
}
\arguments{
\item{rna}{a tibble containing mRNA expression values for genes across patient samples}

\item{protein}{a tibble containing protein expression values for genes across patient samples}
}
\value{
a single tibble containing the correlation of RNA and proteins across all shared genes.
}
\description{
create_rna_prot_correlation
Creates a table with correlation between RNA and protein levels for each gene found in both the RNA and protein data.
}
