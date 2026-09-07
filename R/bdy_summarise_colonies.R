#' Summarizes a table of bird counts into a table of colonies
#'
#' @param countData data frame of bird counts (at species level), where each line gives the number of birds counted for a given colony, species, and year. Must have at least the following columns:
#'                    \itemize{
#'                    \item 'colony': character, name of the colony
#'                    \item 'lat': numeric, latitude of the colony
#'                    \item 'lon': numeric, longitude of the colony
#'                    \item 'seafront': character, name of the seafront on which the colony is located (for instance 'Atlantic', 'Mediterranean'...)
#'                    }
#'
#' @returns Sf table of colony information ('colony', 'lat', 'lon', 'seafront', 'colony_code', and 'geometry'), where each line is a different colony and the following columns:
#' @export
#'
#'
bdy_summarise_colonies <- function(countData){

  colonies00 <- dplyr::select(.data = countData,colony,lat, lon,seafront) %>% unique

  colonies00$colony_code <- paste0("colo_", sprintf("%03d", 1:nrow(colonies00)))
  colonies00 <- colonies00[which(!is.na(colonies00$lat)),]

  colonies00 <- st_as_sf(colonies00, coords = c("lon", "lat"), crs = 4326, agr = "constant",remove=F) %>% st_transform(2154)

  return(colonies00)
}
