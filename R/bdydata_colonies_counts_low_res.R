#' Bird Dynamic low resolution marine bird counts
#'
#' This table gives accurate counts of breeding pairs of marine birds along the french coastline, between 2009 and 2021,
#' based on data aggregated by the french Groupement d'Intéret Scientifique Oiseaux Marins (GISOM).
#' It is the primary data source used by the model to compute projections of future population dynamics.
#' /!\ CAUTION /!\ Because of the sensitive nature of these data, colonies have been clustered to lower the spatial resolution,
#' thus each line of the table represents a group of colonies (within 5km of each other + isolated colonies within 20km), with the centroid given as coordinates.
#' This may slightly alter the results of the population models. If you need to use the unaltered dataset, you can either contact the GISOM, of use the dedicated online
#' application at https://shiny.cefe.cnrs.fr/Shiny_Bird_dynamic.
#'
#' Includes the following columns: \itemize{
#'      \item group: unique identifier for the group of colonies
#'      \item lat: latitude of the group centroid
#'      \item lon: longitude of the group centroid
#'      \item seafront: on which seafront is the group located (either "atlantic" or "mediterranean)
#'      \item species_latin: latin name of the species
#'      \item species_fr: french vernacular name of the species
#'      \item species_en: english vernacular name of the species
#'      \item year: year at which the count was made
#'      \item count: number of breeding couples recorded for the given colony group, species, and year
#'      \item colony: unique identifier, redundant whith column group but kept for consistency
#'      \item regroup: empty column, kept for consistency
#' }
#'
#' @docType data
#'
#' @usage bdydata_colonies_counts_low_res
#'
#' @format Dataframe (.rda)
"bdydata_colonies_counts_low_res"
