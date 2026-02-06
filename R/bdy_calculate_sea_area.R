#' Calculate proportion of sea around colony
#'
#' Calculates the proportion of sea area around the colonies of interest (the ratio between marine surface and terrestrial surface within the foraging range of the species of interest).
#' Used by [bdy_apportionning()] to compute mortality weights between colonies.
#'
#' @param shapeFile shapefile of the area of interest
#' @param max_foraging_range numeric, max foraging range of the species of interest in kilometers
#'  @param colonies sf object with each row being a colony where the species of interest is present. Must have at least the following columns:
#'                \itemize{
#'                \item 'geometry': coordinates of the colony
#'                \item 'code_colonie': character, unique identifier for the colony.
#'                }
#'
#' @returns named numeric vector indicating the proportion of sea area for each colony provided
#' @export

bdy_calculate_sea_area <- function(shapeFile,max_foraging_range,colonies){

  # on rassemble les regions pour avoir la france en entier dans le mm shape
  system.time( shape <- st_union(shapeFile) )

  # en L93
  system.time(shape_L93 <- st_transform(shape, 2154))

  load(paste0("3.BV_Data_processed/processed_data_per_species/clean_data_",esp,".RData"))

  # buffer de X km autour de tes points
  buf <- st_buffer(colonies, set_units(max_foraging_range_km, "km"))
  buf$id <- 1:nrow(buf)
  buf$buffer_area <- st_area(buf)

  # aire de "terre = france" dans chaque buffer qui recoupe la "terre"
  system.time(
    intersect_buf <- st_intersection(st_as_sf(shape_L93), buf) %>%
      mutate(intersect_area = st_area(.))   # create new column with shape area
  )

  ### Pourcentage #####

  # on r?cup?re les info des chevauchements des buffers avec la terre
  st_geometry(buf) <- NULL # on retire la g?om?trie
  st_geometry(intersect_buf) <- NULL # on retire la g?om?trie
  intersect_buf <- intersect_buf[,c("id","intersect_area")]
  buf_ok <- merge(buf, intersect_buf, by = "id", all.x = T) # on merge

  # on remplace NA dans colonne "intersect_buf" par 0
  buf_ok$intersect_area[is.na(buf_ok$intersect_area)] <- 0

  # on calcule le pourcentage de terre
  buf_ok$pct_terre <- as.numeric(buf_ok$intersect_area)/as.numeric(buf_ok$buffer_area) # entre 0 et 1
  buf_ok$pct_terre

  # pourcentage de mer
  sea_area <- (1 - buf_ok$pct_terre)
  names(sea_area) <- colonies$code_colonie

  return(sea_area)
}
