#' Prepare polygon of country shorelines
#'
#' This polygon of shorelines in the area covered by colonies selected is used in maps and in geospatial analyses (e.g., in [bdy_calculate_sea_area()]). The map is based on the global map included in the package [bdydata_world_map].
#'
#' @param colonies sf object including colonies to cover by the polygon
#' @param buffer buffer size (in km) around the colonies to include in the polygon
#'
#' @returns sf object of countries shorline in EPSG 2154
#' @export
#'
#' @examples
#' \dontrun{
#' data_pts <- st_as_sf(data.frame(lon=c(-2.3, -5.1, -1.9), lat=c(46.7, 48.5, 49.8)), coords = c("lon","lat"), remove = FALSE, crs=st_crs(4326))
#' countryShape <- bdy_prepare_countryShape(data_pts, buffer=100)
#' ggplot(countryShape)+geom_sf()
#' }
#'
bdy_prepare_countryShape <- function(colonies, buffer){

  # Prepare countryShape
  crop_extent <- colonies %>% st_bbox() %>% st_as_sfc() %>% st_buffer(., buffer*1000) %>% st_transform(., 4326)

  countryShape <- bdydata_world_map %>%
    st_crop(., crop_extent) %>%
    st_transform(., st_crs(2154))

  return(countryShape)
}
