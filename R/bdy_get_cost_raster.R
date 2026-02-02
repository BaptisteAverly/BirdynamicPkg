#' Title
#'
#' @param world_map
#' @param region
#' @param colonies_geom
#' @param parcs_geom
#' @param N_buffer
#' @param S_buffer
#' @param W_buffer
#' @param E_buffer
#' @param pixel_size
#'
#' @returns
#' @export
#'
#' @examples

bdy_get_cost_raster <- function(world_map,region=c("Northern Europe","Southern Europe","Western Europe"),
                                colonies_geom,parcs_geom,N_buffer=5000,S_buffer=5000,W_buffer=5000,E_buffer=5000,pixel_size=1000){


  countries <- world_map[world_map$SUBREGION %in% region,]

  # Transform in L93 projection
  countries_L93 <- st_transform(countries, crs = 2154)

  # Take a subset of the raster based on the extent of colonies and parcs locations
  ext_col <- rbind(colonies_geom, parcs_geom) %>% st_bbox

  #st_bbox(colonies_L93)

  ext_col <- extent(c(
    ext_col$xmin-W_buffer,
    ext_col$xmax+E_buffer,
    ext_col$ymin-S_buffer,
    ext_col$ymax+N_buffer))

  # create a raster with this extent which will serve as the cost surface
  rast_cost <- raster(ext_col)

  # pixel size (m) : 1000 m.
  res(rast_cost) <- pixel_size

  # transfer the coordinate system to the raster
  crs(rast_cost) <- CRS("+init=epsg:2154") %>% suppressWarnings()

  # Now you can add data to the cells in your raster to mark the ones that fall within your polygon.
  (rast_cost[as_Spatial(countries_L93),] <- 10000) %>% suppressWarnings()

  ## Create cost surface : "land" = high cost, sea = low cost
  #rast_cost[rast_cost == 1] <- 10000
  rast_cost[is.na(rast_cost)] <- 1

  ## Produce transition matrices, and correct because 8 directions
  system.time( trCost <- transition(1/rast_cost, mean, directions=8) )
  system.time( trCost <- geoCorrection(trCost, type="c") )

  return(list("transition_layer"=trCost,"cost_raster"=rast_cost))
}
