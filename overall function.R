
setwd(dir=ifelse(file.exists("C:/Users/Victor"),
                 "C:/Users/Victor/Documents/Projects/Birdynamic", # Setwd Victor
                 "D:/Documents/JOBS/Independant/Birdynamic")) # Setwd Baptiste

#these should be included in the package
world_map <- st_read("0.BV_Data_input/world_shpFile/ne_10m_admin_0_countries.shp")
load("0.BV_Data_input/species_parameters.RData")
modelFilePath <- "app/0.Data/model_01.txt"

#these are user input
countData <-  read.csv2("0.BV_Data_input/BD_Effectifs_clean.csv")
load("Test_data_Shiny/clean_data_parcs.RData")
parcs <- parcs_L93
parcs$NAME <- bdy_clean_names(parcs$NAME)
parcs$seafront <- "atlantique"
species = seasons$species_latin[1:4]
countryShape <- st_read("0.BV_Data_input/France_shp/france.shp") #shapefile for france could be included in package as default
#or sea_area function could be made better to get country of interest from world map
timeRange = c(2009,2021)
mortality <- read.csv(file="Test_data_Shiny/mortality_all/formatted/01_collisionRisk_Fou de Bassan.csv")
mortality <- bdy_check_mortality(mortality,seasons$species_latin)$table #include inside overall function ?
mortality$species_latin <- mortality$espece_latin

bdy_run_analysis <- function(species,countData,countryShape,parcs,timeRange = c(2009,2021),
                             n_iteration=1000,ni_noImpact=10000,ni_withImpact=1000,...){

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

  #which species of interest are strictly marine
  marine_sp <- species[!foraging_ranges$terrestrial_habits[match(species,foraging_ranges$species_latin)]]

  #calculating distances
  all_distances <- bdy_get_distances(colonies=colonies_all,
                                     parcs = parcs,
                                     costMatrix=cost_matrix,
                                     doShpa = length(marine_sp)>0)

  #preparing a list to store model outpu for each species
  model_output <- list()

  for(sp in species){

    #processing count data
    count_processed <- bdy_process_count_data(sp=sp,
                                              effec00=countData,
                                              colonies00=colonies_all,first_year=timeRange[1],last_year = timeRange[2])

    # !! prevoir cas ou l'espece d'interet n'est pas dans la liste !!
    foraging_range_sp = foraging_ranges$max_km[which(foraging_ranges$species_latin==sp)]
    terrestrial_habbit = foraging_ranges$terrestrial_habits[which(foraging_ranges$species_latin==sp)]

    #selecting apropriate distance table (shpa or eucl), and only for colonies where species is present
    distances <-  all_distances[[match(terrestrial_habbit,c(T,F))]]
    idx <- which(row.names(distances) %in% count_processed$colonies_sp$colony_code)
    distances <- distances[idx,]

    #calculating sea area for each colony
    sea_area <- bdy_calculate_sea_area(countryShape,
                           max_foraging_range = foraging_range_sp,
                           count_processed$colonies_sp)

    #apportionning
    apportionning <- bdy_apportionning(max_foraging_range_km=foraging_range_sp,
                                       colonies=count_processed$colonies_sp,
                                       sea_area=sea_area,
                                       tbl_dist = distances)

    #distributing mortality
    morta_iter_group <- bdy_processing_mortality(
      collision = mortality[which(mortality$species_latin==sp),],
      season=seasons[which(seasons$species_latin==sp),][month.abb],
      n_iteration=n_iteration,
      RW_group = apportionning$RW_group
    )

    #model without impact
    no_impact_output <- bdy_model_no_impact(count_data=count_processed$group_counts_sp,
                                            PI=count_processed$ppa_yr_sp,
                                            survival=vital_rates[vital_rates$species_latin== sp, "survival"],
                                            fecundity=vital_rates[vital_rates$species_latin == sp, "fecundity"],
                                            propRepro=vital_rates[vital_rates$species_latin == sp, "propRepro"],
                                            modelFile=modelFilePath, nimble=F,
                                            ni=ni_noImpact,
                                            lightResults = T)

    #model with impact
    with_impact_output <- bdy_model_with_impact(n_group=nlevels(count_processed$colonies_sp$group),
                                                ny_data=ncol(count_processed$group_counts_sp),
                                                posterior=no_impact_output,
                                                mortality=morta_iter_group,
                                                ni=ni_withImpact)

    #preparing output
    units(foraging_range_sp)<-"km"
    distanceBIN <- distances<foraging_range_sp
    rownames(distanceBIN)<-rownames(distances)
    names(distanceBIN) <- names(distances)

    model_output[[sp]] <- list(no_impact=with_impact_output$sc0,
                               with_impact=with_impact_output$sc1,
                               colonies=count_processed$colonies_sp,
                               mortality=morta_iter_group,
                               distance=round(distances),
                               distanceBIN=distanceBIN)
  }

  #raw results
  Raw_ResTables <- bdy_raw_res_tables(mod_out=model_output)

  #pretty results
  Pretty_table_national <- bdy_pretty_result_table(Raw_ResTables, type="national")

  #summary figure
  Summary_Figures <- bdy_summary_figure(Raw_ResTables, Pretty_table_national)
}
