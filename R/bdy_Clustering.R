#' Clustering
#'
#'Clustering using Mean Shift Algorithm
#'
#' @param coord
#' @param h_km Define the bandwidth (in km) used by the Clustering Algorithm
#' @param dist_regr_isol Define the max distance to regroup an isolated colony to closest cluster
#' @param regroupIsolates
#'
#' @returns
#' @export
#'
#' @examples
bdy_clustering <- function(coord,h_km=5,dist_regr_isol=4*h_km,regroupIsolates=T){

  ## Define range of Lat/Lon values
  rg <-  apply(coord, 2, range)

  ## Convert the bandwith in terms of % of Lat/Lon range extent
  # Latitude: 1 deg = 110.574 km.
  h_Y <- h_km / (diff(rg[,"lat"]) * 110.574)

  # Longitude: 1 deg = 111.320*cos(latitude) km.
  deg2rad <- function(deg) {(deg * pi) / (180)}
  h_X <- h_km / (diff(rg[,"lon"]) * cos(deg2rad(mean(rg[,"lat"]))) * 111.320)

  ## Apply the Mean Shift Algorithm (from LPCM package)
  if(h_km > 0){
    ms_res <- st_coordinates(coord) %>% LPCM::ms(., h = c(h_X, h_Y), plot = F)
    group <- as.factor(ms_res$cluster.label)
  }else{
    group <- as.factor(1:nrow(coord))
  }

  ## Regrouping of isolated cluster (small spatial scale : 15km)
  ## Regroupement des cluster isolés

  if(regroupIsolates){

    old_group <- group

    if(h_km > 0){
      # Identifier les clusters/groupes isolés
      iso_clus <- which(table(group) == 1)

      # Distances entre paires de colonies
      colo_dist <- (st_distance(coord)) %>% set_units(., km)
      diag(colo_dist) <- Inf

      ## Loop over
      if(length(iso_clus) > 0){
        for(j in 1:length(iso_clus)){
          iso_colo_j <- which(old_group == iso_clus[j])  # identify the colony index from that isolated group
          colo_closest <- which.min(colo_dist[iso_colo_j,])         # identify its closest friend-colony

          # Apply its friend group number only if distance is less than XXXX
          if(colo_dist[iso_colo_j, colo_closest] < set_units(dist_regr_isol, km)){
            group[iso_colo_j] <- group[colo_closest]
          }
        } # j
      } # if 1

      ## Renumber groups
      group <- group %>% factor %>% as.numeric %>% as.factor
    }
  }
  return(group)
}
