#' Get cost transition matrix
#'
#' Produces the transition matrix and cost raster used to calculate shortest path distances for birds with marine exclusive livestyles.
#'
#' @param world_map sf object of the world with every country. Can be downloaded here:
#'                  https://github.com/nvkelso/natural-earth-vector/blob/master/10m_cultural/ne_10m_admin_0_countries.shp
#' @param region regions of the world where the analysis should take place. Subsetting on the regions of interest can speed up the calculation.
#'                Can be any or a combination of the following: \cr
#'                "All","South-Eastern Asia","South America","Western Asia","Southern Asia","Eastern Asia",
#'                "Eastern Africa","Western Europe","Northern Africa","Central America","Middle Africa",
#'                "Eastern Europe","Southern Africa","Caribbean","Central Asia","Northern Europe",
#'                "Southern Europe","Western Africa","Northern America","Melanesia","Antarctica",
#'                "Australia and New Zealand","Polynesia","Micronesia"
#' @param colonies_geom vector of geometries indicating the coordinates of all the bird colonies of interest.
#'                      This is used to get the extent of the raster to calculate.
#' @param parcs_geom vector of geometries indicating the coordinates of all the wind farms of interest.
#'                    This is used to get the extent of the raster to calculate.
#' @param N_buffer north buffer to increase the extent of the raster beyond the northernmost colony or parc (numeric, in kilometers)
#' @param S_buffer south buffer to increase the extent of the raster beyond the southernmost colony or parc (numeric, in kilometers)
#' @param W_buffer west buffer to increase the extent of the raster beyond the westernmost colony or parc (numeric, in kilometers)
#' @param E_buffer east buffer to increase the extent of the raster beyond the easternmost colony or parc (numeric, in kilometers)
#' @param pixel_size numeric, size in meter of the cells represented by one pixel of the raster. Lower numbers give higher resolution.
#'
#' @returns List of 2 elements:
#'          \itemize{
#'          \item 'transition_matrix': transitionLayer object where sea is represented as low cost and land is represented as high cost
#'          \item 'cost_raster': raster object where sea is represented as low cost and land is represented as high cost
#'          }
#' @export
#'

bdy_get_cost_raster <- function(world_map,region=c("Northern Europe","Southern Europe","Western Europe"),
                                colonies_geom,parcs_geom,N_buffer=10,S_buffer=10,W_buffer=10,E_buffer=10,pixel_size=10){

  if(region == "All"){
    countreis <- world_map
  }else{
    countries <- world_map[world_map$SUBREGION %in% region,]
  }

  # Transform in L93 projection
  countries_L93 <- st_transform(countries, crs = 2154)

  # Take a subset of the raster based on the extent of colonies and parcs locations
  ext_col <- rbind(colonies_geom, parcs_geom) %>% st_bbox

  #st_bbox(colonies_L93)

  ext_col <- extent(c(
    ext_col$xmin-(W_buffer*1000),
    ext_col$xmax+(E_buffer*1000),
    ext_col$ymin-(S_buffer*1000),
    ext_col$ymax+(N_buffer*1000)))

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

  return(list("transition_matrix"=trCost,"cost_raster"=rast_cost))
}
