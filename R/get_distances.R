
#' Calcul des distances entre parc éoliens et colonies
#'
#' @param colonies tableau de colonies
#' @param parcs shapefile des parcs
#' @param costMatrix matrice de couts
#' @param doShpa booleen
#' @param progress uniquement pour une utilisation avec shiny
#'
#' @returns matrice de distance parc-colonies
#' @export
#'
bdy_get_distances <- function(colonies,parcs,costMatrix,doShpa=T,progress=NULL){

  eucl_dist <- (st_distance(x = colonies, y = parcs)) %>% set_units(., km)
  rownames(eucl_dist) <- colonies$code_colonie
  colnames(eucl_dist) <- parcs$NAME

  n_colo <- nrow(colonies)
  n_parc <- length(parcs$NAME)

  ## Calculate Euclidean Distance
  if(!is.null(progress)){
    progress$set(1, detail = paste0("Parc ",n_parc,"  / ",n_parc))
  }

  ## Calculate shortest path between colonies and wind farms
  goal <- st_coordinates(st_centroid(colonies)) %>% suppressWarnings()
  iMax <- n_parc

  # Object to store shortest path distances between each colony and each parc
  shpa_dist <- eucl_dist
  shpa_dist[] <- NA
  rownames(shpa_dist) <- colonies$code_colonie
  colnames(shpa_dist) <- parcs$NAME

  if(doShpa){
    for(i in 1:iMax){

      origin <- st_coordinates(st_centroid(parcs))[i,] %>% t()

      shortPath <- shortestPath(x = costMatrix,
                                origin = origin,
                                goal = goal,
                                output = "SpatialLines") %>% suppressWarnings()

      #crs(shortPath) <- CRS("+init=epsg:2154") %>% suppressWarnings()

      shpa_dist[,i] <- st_length(st_as_sf(shortPath), which = "Euclidean") /1000

      #updating progress
      if(!is.null(progress)){
        progress$inc(1/iMax, detail = paste0("Parc ", i,"  / ",n_parc))
      }
    } # i
  }
  return(list(eucl_dist=eucl_dist,shpa_dist=shpa_dist))
}
