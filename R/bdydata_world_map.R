#' Polygons of the world
#'
#' Map of the world, to be used in functions [bdy_get_cost_raster()] and [bdy_calculate_sea_area()]. Be aware that at the moment the package use Lambert93 projection so can only be used for European regions (mainly designed for France).
#'
#' @docType data
#'
#' @usage bdydata_world_map
#'
#' @source Taken from https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-0-countries/
#'
#' @format sf object (.rda)
"bdydata_world_map"
