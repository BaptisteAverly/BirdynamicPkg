#' Calculates weights to distribute bird mortality between colonies and groups of colonies for a given species.
#' for a given colony, these weights are influenced by: (1) Distance between the colony and the wind parcs,
#' (2) The population size of the species of interest for that colony and (3), the proportion of marine surface surrounding the colony.
#'
#' @param max_foraging_range_km numeric, maximum foraging range in kilometers for the species of interest
#' @param colonies data frame with each row being a colony where the species of interest is present. Must have at least the following columns:
#'                'group': numeric, colonies with the same number belong to the same cluster
#'                'avg': numeric, average population size in that colony over the years of interests
#'                'code_colonie': character, unique identifier for the colony.
#' @param sea_area named numeric vector giving for each colony the proportion of marine surface surrounding it.
#'                  Each value should be named with the unique identifier for the corresponding colony.
#' @param tbl_dist matrix giving the distances between colonies (rows) and wind farms (columns), as outputed by bdy_get_distances()
#' @param incl_pop_size boolean, whether population size should influence the weight calculation
#' @param incl_sea_area boolean, whether the proportion of marine surface surrounding the colony should influence the weight calculation
#'
#' @returns list of 4 matrices:
#'          'AW_colo':
#' @export
#'
#' @examples

bdy_apportionning <- function(max_foraging_range_km,colonies,sea_area,tbl_dist,incl_pop_size=T,incl_sea_area=T){
  n_parc = ncol(tbl_dist)
  n_group <- nlevels(colonies$group) #B

  ### APPORTIONNING ###
  ## Get Relative Weights for Apportioning

  ### Create Relative Weights (AW) for Apportioning fatalities

  ## Filter colonies for the given species
  #tbl_dist <- tbl_dist[which(rownames(tbl_dist) %in% colonies$code_colonie), ]

  ## Ensure it is ordered by code_colonie
  tbl_dist <- tbl_dist[order(row.names(tbl_dist)),,drop=F]

  ## Table for Apportioning
  parcNames = colnames(tbl_dist)
  tbl_app <- as.data.frame(tbl_dist)
  colnames(tbl_app) = parcNames

  ## Add info colony size
  tbl_app$size <- colonies$avg[match(row.names(tbl_dist),colonies$code_colonie)]
  #tbl_app$size <- colonies$last[match(row.names(tbl_dist),colonies$code_colonie)]

  ## Add info sea_area
  tbl_app$sea_area <- sea_area[row.names(tbl_dist)]

  ## Avoid negative values for sea area
  tbl_app$sea_area[tbl_app$sea_area < 0] <- min(tbl_app$sea_area[tbl_app$sea_area > 0])

  # check: which(!( (names(sea_area)) %in% ((colonies$code_colonie))) )
  # check: cbind((names(sea_area)), (colonies$code_colonie) )

  ## Define function to get relative weights
  # rel_weight <- function(x) if(sum(x) == 0) 0 else x/sum(x)
  rel_weight <- function(x) x/sum(x)

  ## Relation "risque relatif" et distance (d²)
  #d <- seq(10, 100, length.out = 100)
  #rr <- (1/rel_weight(d^2)) %>% rel_weight
  # plot(x = d, y = rr, type = "l", ylab = "Risque relatif", xlab = "Distance (km)")


  #####################################
  #### COLONY SCALE
  #####################################

  ## Make table of Absolute Weights
  AW_dist <- AW_size <- AW_area <- as.data.frame(matrix(NA, nrow = nrow(tbl_app), ncol = n_parc, dimnames = list(rownames(tbl_app), colnames(tbl_dist))))

  ## Relative weight including or not pop size and sea area
  for(i in 1:n_parc){
    AW_dist[,i] <- as.numeric(1/rel_weight((tbl_app[,i]^2)))
    AW_size[,i] <- rel_weight(tbl_app$size)
    AW_area[,i] <- rel_weight(1/tbl_app$sea_area)
  } # i
  AW_dist %>% colSums
  AW_size %>% colSums
  AW_area %>% colSums

  (AW_dist * AW_size * AW_area) %>% colSums

  ## Exclude distances beyond foraging range(==> AW = 0)
  lim_dist_km <- max_foraging_range_km
  AW_dist[tbl_dist > set_units(lim_dist_km, km)] <- 0
  AW_dist %>% colSums

  ## Total AW
  AW_colo <- AW_dist
  if(incl_pop_size) AW_colo <- AW_colo * AW_size
  if(incl_sea_area) AW_colo <- AW_colo * AW_area
  AW_colo %>% colSums

  ## Rescale AW : RELATIVE WEIGHT
  RW_colo <- sapply(AW_colo, rel_weight)
  for(j in 1:ncol(RW_colo)) RW_colo[,j][is.nan(RW_colo[,j])] <- 0
  rownames(RW_colo) <- rownames(tbl_app)
  RW_colo %>% colSums


  #####################################
  #### CLUSTER/GROUP SCALE
  #####################################
  group_dist <- as.data.frame(matrix(NA, nrow = n_group, ncol = n_parc, dimnames = list(levels(colonies$group), colnames(tbl_dist))))
  group_size <- group_area <- c()

  for(gg in 1:n_group){

    ## Distance moyenne au group pondérée par la taille de colonie
    group_dist[gg,] <-
      (rel_weight(tbl_app[colonies$code_colonie[colonies$group == gg], "size"]) *
         tbl_app[colonies$code_colonie[colonies$group == gg], 1:n_parc,drop=F]
      ) %>%
      colSums

    ## Relative weight Mean pop size of the group
    group_size[gg] <- sum(tbl_app[colonies$code_colonie[colonies$group == gg], "size"])

    ## Mean sea area of the group, weighted by colony size
    group_area[gg] <-
      (rel_weight(tbl_app[colonies$code_colonie[colonies$group == gg], "size"]) *
         tbl_app[colonies$code_colonie[colonies$group == gg], "sea_area"]
      ) %>%
      sum

  } # gg
  rm(gg)

  ## Make table of Absolute Weights
  AW_dist <- AW_size <- AW_area <- as.data.frame(matrix(NA, nrow = n_group, ncol = n_parc, dimnames = list(levels(colonies$group), colnames(tbl_dist))))
  #dim(AW_size)

  for(i in 1:n_parc){
    AW_dist[,i] <- as.numeric(1/rel_weight((group_dist[,i]^2)))
    AW_size[,i] <- rel_weight(group_size)
    AW_area[,i] <- rel_weight(1/group_area)
  } # i

  AW_dist %>% colSums
  AW_size %>% colSums
  AW_area %>% colSums

  (AW_dist * AW_size * AW_area) %>% colSums

  ## Exclude distances beyond foraging range(==> AW = 0)
  lim_dist_km <- max_foraging_range_km
  AW_dist[group_dist > lim_dist_km] <- 0
  AW_dist %>% colSums

  ## Total AW
  AW_group <- AW_dist
  if(incl_pop_size) AW_group <- AW_group * AW_size
  if(incl_sea_area) AW_group <- AW_group * AW_area
  AW_group %>% colSums

  ## Rescale AW : RELATIVE WEIGHT
  RW_group <- sapply(AW_group, rel_weight)
  for(j in 1:ncol(RW_group)) RW_group[,j][is.nan(RW_group[,j])] <- 0
  rownames(RW_group) <- sort(unique(colonies$group))

  ## Ensure proper order by group number
  RW_group <- RW_group[order(as.numeric(rownames(RW_group))),,drop=F]
  RW_group %>% colSums

  return(list("AW_colo"=AW_colo,"RW_colo"=RW_colo,"AW_group"=AW_group,"RW_group"=RW_group))
}
