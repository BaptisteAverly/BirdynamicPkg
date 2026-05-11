

bdy_check_colonies <- function(effecTable){
  
  #-------add checks-------
  
  
  
  #------------------------
  
  colonies00 <- dplyr::select(.data = effecTable,colony,lat, lon,seafront) %>% unique
  
  colonies00$colony_code <- paste0("colo_", sprintf("%03d", 1:nrow(colonies00)))
  colonies00 <- colonies00[which(!is.na(colonies00$lat)),]
  
  colonies00 <- st_as_sf(colonies00, coords = c("lon", "lat"), crs = 4326, agr = "constant",remove=F) %>% st_transform(2154)
  
  return(colonies00)
}
