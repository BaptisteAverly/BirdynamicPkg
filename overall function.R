
setwd(dir=ifelse(file.exists("C:/Users/Victor"),
                 "C:/Users/Victor/Documents/Projects/Birdynamic", # Setwd Victor
                 "D:/Documents/JOBS/Independant/Birdynamic")) # Setwd Baptiste

#these should be included in the package
world_map <- st_read("0.BV_Data_input/world_shpFile/ne_10m_admin_0_countries.shp")
load("0.BV_Data_input/species_parameters.RData")

#these are user input
countData <-  read.csv2("0.BV_Data_input/BD_Effectifs_clean.csv")
load("Test_data_Shiny/clean_data_parcs.RData")
parcs <- parcs_L93
parcs$NAME <- bdy_clean_names(parcs$NAME)
parcs$seafront <- "atlantique"
species = seasons$species_latin[1:4]
countryShape <- st_read("0.BV_Data_input/France_shp/france.shp") #shapefile for france could be included in package as default
#ou sea_area function could be made better to get country of interest from world map
timeRange = c(2009,2021)
mortality <- read.csv(file="Test_data_Shiny/mortality_all/formatted/04_collisionRisk_Cormoran huppé.csv")
mortality <- bdy_check_mortality(mortality,seasons$species_latin)$table #include inside overall function ?
mortality$species_latin <- mortality$espece_latin

bdy_run_analysis <- function(species,countData,countryShape,parcs,timeRange = c(2009,2021),n_iteration=1000,...){

  ##checking count data and adding colonies code
  colonies_all <- bdy_check_colonies(countData)

  ## Add colony_code  to "countData"
  countData$colony_code <- colonies_all$colony_code[match(countData$colony,colonies_all$colony)]

  #calculating cost matrix for each seafront
  cost_matrix <- list()
  for(seafront in unique(colonies_all$seafront)){
    seafrontIdx <- which(colonies_all$seafront==seafront)
    cost_matrix[[seafront]] <- bdy_get_cost_raster(world_map,
                                                   pixel_size=10000, #to speed up testing, remove in final function
                                                   geom=colonies_all$geometry[seafrontIdx])$transition_matrix
  }

  marine_sp <- species[!foraging_range_tbl$terrestrial_habits[match(species,foraging_range_tbl$species_latin)]]

  all_distances <- bdy_get_distances(colonies=colonies_all,
                                     parcs = parcs,
                                     costMatrix=cost_matrix,
                                     doShpa = length(marine_sp)>0)

  for(sp in species){

    count_processed <- bdy_process_count_data(sp=sp,
                                              effec00=countData,
                                              colonies00=colonies_all,first_year=timeRange[1],last_year = timeRange[2])

    # !! prevoir cas ou l'espece d'interet n'est pas dans la liste !!
    max_foraging_range = foraging_range_tbl$max_km[which(foraging_range_tbl$species_latin==sp)]
    terrestrial_habbit = foraging_range_tbl$terrestrial_habits[which(foraging_range_tbl$species_latin==sp)]

    distances <-  all_distances[[match(terrestrial_habbit,c(T,F))]]
    idx <- which(row.names(distances) %in% colonies$colony_code)
    distances <- distances[idx,]

    sea_area <- bdy_calculate_sea_area(countryShape,
                           max_foraging_range = max_foraging_range,
                           count_processed$colonies_sp)

    apportionning <- bdy_apportionning(max_foraging_range_km=max_foraging_range,
                                       colonies=count_processed$colonies_sp,
                                       sea_area=sea_area,
                                       tbl_dist = distances)

    morta_iter_group <- bdy_processing_mortality(
      collision = mortality[which(mortality$species_latin==sp),],
      season=seasons[which(seasons$species_latin==sp),][month.abb],
      n_iteration=n_iteration,
      RW_group = apportionning$RW_group
    )

  }



}
