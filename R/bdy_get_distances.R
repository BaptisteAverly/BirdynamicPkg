#' Calculates distance parc-colonies
#'
#' Calculates the shortest path distances between marine bird colonies and wind farms
#'
#' @param colonies table where each row is a colony of interest.
#'                Must contain the columns:
#'                \itemize{
#'                \item 'code_colonie': character, unique identifiers for each colony
#'                \item 'geometry': sf coordinates of the colony centroids
#'                \item 'facade': character, which sea front is the colony located in
#'                }
#' @param parcs table where each row is a wind farm of interest.
#'                Must contain the columns:
#'                \itemize{
#'                \item 'NAME': character, unique identifier for each wind farm
#'                \item 'geometry': sf coordiantes of the wind famr centroid
#'                \item 'Facade': character, which sea front is the colony located in
#'                }
#' @param costMatrix named list containing cost rasters (Transition object from package gdistance) representing the land areas as costly and marine areas as cost free.
#'                   Must contain one object per sea front, with names identical to names in column 'facade' of argument 'colonies' and column 'Facade' of argument 'parcs'.
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
  goals <- st_coordinates(st_centroid(colonies)) %>% suppressWarnings()
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

      #coordinates of the current parc
      origin <- st_coordinates(st_centroid(parcs))[i,] %>% t()

      #which façade is the park in
      facade <- parcs$Facade[i]

      #idx to select only the colonies which are on this façade
      facadeIdx <- which(colonies$facade == facade)

      shortPath <- shortestPath(x = costMatrix[[facade]],
                                origin = origin,
                                goal = goals[facadeIdx,],
                                output = "SpatialLines") %>% suppressWarnings()

      #crs(shortPath) <- CRS("+init=epsg:2154") %>% suppressWarnings()

      #storing the distances for the colonies on the right façade in the table
      shpa_dist[facadeIdx,i] <- st_length(st_as_sf(shortPath), which = "Euclidean") /1000
      #colonies which are not on the same façade as the parc get very high shortest path distances (no influence)
      shpa_dist[which(colonies$facade != facade),i] <- 10000

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
