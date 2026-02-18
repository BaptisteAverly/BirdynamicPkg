#' Calculates distance parc-colonies
#'
#' Calculates the shortest path distances between marine bird colonies and wind farms
#'
#' @param colonies table where each row is a colony of interest.
#'                Must contain one column named "code_colonie" with unique identifiers for each colony
#'                and one spatial "geometry" column containing the coordinates of the colony centroids
#' @param parcs table where each row is a wind famr of interest.
#'                Must contain one column named "NAME" with unique identifiers for each wind farm
#'                and one spatial "geometry" column containing the coordinates of the wind farm centroids
#' @param costMatrix cost raster (Transition object from package gdistance) representing the land areas as costly and marine areas as cost free.
#'                   Used to calculate shortest path distances for strictly marine birds which do not fly over the land.
#' @param doShpa boolean. If set to false, only the euclidean distance is calculated which makes the computation faster, especially if there are many wind farms.
#' @param progress R shiny Progress object used to show the calculation progress to the user. Only for use within a shiny application.
#'
#' @returns List containing 2 matrices with identical format (colonies as rows and wind farms as columns): \cr
#'          'eucl_dist' contains the euclidian distances and 'shpa_dist' contains the shortest path distances (filled with NAs if doShpa=F)
#'
#' @seealso [bdy_get_cost_raster()]
#' @export
#'
bdy_get_distances <- function(colonies,parcs,costMatrix,doShpa=T,progress=NULL){

  ## Calculate Euclidean Distance
  eucl_dist <- (st_distance(x = colonies, y = parcs)) %>% set_units(., km)
  rownames(eucl_dist) <- colonies$code_colonie
  colnames(eucl_dist) <- parcs$NAME

  n_colo <- nrow(colonies)
  n_parc <- length(parcs$NAME)

  ## Calculate shortest path between colonies and wind farms
  goal <- st_coordinates(st_centroid(colonies)) %>% suppressWarnings()
  iMax <- n_parc

  # Object to store shortest path distances between each colony and each parc
  shpa_dist <- eucl_dist
  shpa_dist[] <- NA
  rownames(shpa_dist) <- colonies$code_colonie
  colnames(shpa_dist) <- parcs$NAME

  if(!is.null(progress)){
    progress$set(0.05, detail = paste0("Parc 1  / ",iMax))
  }

  if(doShpa){
    for(i in 1:iMax){

      origin <- st_coordinates(st_centroid(parcs))[i,] %>% t()

      shortPath <- shortestPath(x = costMatrix,
                                origin = origin,
                                goal = goal,
                                output = "SpatialLines") %>% suppressWarnings()

      #crs(shortPath) <- CRS("+init=epsg:2154") %>% suppressWarnings()

      shpa_dist[,i] <- st_length(st_as_sf(shortPath), which = "Euclidean") /1000

      # Manually input distance of 10,000 km if not on the same facade
      shpa_dist[colonies$facade!=parcs$Facade[i],i] <- 10000

      #updating progress
      if(!is.null(progress)){
        if(i < iMax){
          progress$inc(1/iMax, detail = paste0("Parc ", i+1,"  / ",iMax))
        }else{
          progress$inc(1/iMax, detail = paste0("Parc ", iMax,"  / ",iMax))
        }

      }
    } # i
  }
  return(list(eucl_dist=eucl_dist,shpa_dist=shpa_dist))
}
