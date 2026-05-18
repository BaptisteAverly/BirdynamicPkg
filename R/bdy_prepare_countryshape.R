
#' Title
#'
#' @param geom 
#' @param buffer 
#'
#' @returns
#' @export
#'
#' @examples
#' 
bdy_prepare_countryShape <- function(geom,buffer){
  
  # Prepare countryShape
  crop_extent <- geom %>% st_bbox() %>% st_as_sfc() %>% st_buffer(., buffer*1000) %>% st_transform(., 4326)
  
  countryShape <- bdydata_world_map %>%
    st_crop(., crop_extent) %>%
    st_transform(., st_crs(2154))
  
  return(countryShape)
}
