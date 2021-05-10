#' @title Genomic locations of cytoband labels
#' @description A dataset containing the chr, start and end position for cytobands
#' according to hg38.
#' @format A data frame with 863 rows and 6 variables:
#' \describe{
#'     \item{chrom}{chromosome}
#'     \item{chromStart}{start position}
#'     \item{chromEnd}{end position}
#'     \item{name}{cytoband name}
#'     \item{gieStain}{color}
#'     \item{color}{HEX color}
#' }
#' @source \url{https://genome.ucsc.edu/cgi-bin/hgTables}
"cytoband_data"

#' @title Names of 2018 TCGA studies from cBioPortal
#' @description A dataset containing the names and studyIds of the 2018 TCGA studies from cBioPortal.
#' @format A data frame with 32 rows and 2 variables:
#' \describe{
#'     \item{Cancer}{Name of diagnosis and sample size}
#'     \item{studyId}{studyId that can be used in the cBioPortalData R package}
#'     }
#' @source \url{https://github.com/waldronlab/cBioPortalData} See data-raw folder.
"cbio_studies"

#' @title Data from 2018 TCGA studies from cBioPortal
#' @description A dataset containing the study name and aggregated gene level copy number data
#' @format A data frame with 14944 rows and 6 variables:
#' \describe{
#'     \item{hugoGeneSymbol}{hugo gene symbol}
#'     \item{Gain}{proportion of cohort with gain in this gene}
#'     \item{Amplification}{proportion of cohort with amplification in this gene}
#'     \item{ShallowDeletion}{proportion of cohort with shallow deletion in this gene}
#'     \item{DeepDeletion}{proportion of cohort with deep deletion in this gene}
#'     \item{study_name}{cancer type and sample size}
#'     }
#' @source \url{https://github.com/waldronlab/cBioPortalData} See data-raw folder.
"all_tcga2018_data"

#' @title Probe data for vignette example
#' @docType data
#' @usage data(probe_data)
#' @description A dataset containing simulated probe data as sample input for launchCNViz
#' @format A data frame with 2006 rows and 6 variables:
#' \describe{
#'     \item{chr}{chromosome}
#'     \item{start}{start location}
#'     \item{end}{end location}
#'     \item{gene}{gene name}
#'     \item{log2}{log2 copy number ratio}
#'     \item{weight}{weight given to log2 value}
#'     }
#' @source Center for Personalized Diagnositics at the University of Pennsylvania
"probe_data"

#' @title Gene data for vignette example
#' @docType data
#' @usage data(gene_data)
#' @description A dataset containing simulated gene data as sample input for launchCNViz
#' @format A dataframe with 112 rows and 6 variables
#' \describe{
#'     \item{chr}{chromosome}
#'     \item{start}{start location}
#'     \item{end}{end location}
#'     \item{gene}{gene name}
#'     \item{log2}{log2 copy number ratio}
#'     \item{weight}{weight given to log2 value}
#'     \item{loh}{loss of heterozygosity}
#'     }
#' @source Center for Personalized Diagnositics at the University of Pennsylvania
"gene_data"

#' @title Segment data for vignette example
#' @docType data
#' @usage data(segment_data)
#' @description A dataset containing simulated segment data as sample input for launchCNViz
#' @format A dataframe with 101 rows and 5 variables
#' \describe{
#'     \item{chr}{chromosome}
#'     \item{start}{start location}
#'     \item{end}{end location}
#'     \item{log2}{log2 copy number ratio}
#'     \item{loh}{loss of heterozygosity}
#'     }
#' @source Center for Personalized Diagnositics at the University of Pennsylvania
"segment_data"

#' @title Variant data for vignette example
#' @docType data
#' @usage data(variant_data)
#' @description A dataset containing simulated SNV and indel data as sample input for launchCNViz
#' @format A dataframe with 119 rows and 4 variables
#' \describe{
#'     \item{gene}{gene name}
#'     \item{mutation_id}{string with information about snv}
#'     \item{depth}{read depth}
#'     \item{start}{starting location}
#'     }
#' @source Center for Personalized Diagnositics at the University of Pennsylvania
"variant_data"

#' @title Metadata for vignette example
#' @docType data
#' @usage data(meta_data)
#' @description A dataset containing simulated metadata as sample input for launchCNViz
#' @format A dataframe with 1 rows and 2 variables
#' \describe{
#'     \item{purity}{sample purity}
#'     \item{ploidy}{tumor ploidy}
#'     }
#' @source Center for Personalized Diagnositics at the University of Pennsylvania
"meta_data"
