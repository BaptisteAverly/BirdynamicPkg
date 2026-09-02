#' Run complete analysis
#'
#' This function can be used to run the complete analysis (from calculating cost raster to getting model outputs). See the package vignette for more information.
#'
#' @param species vector of species names to include in the analysis
#' @param countData table of bird counts
#' @param windfarms sf object indicating the position of windfarms, for instance using [bdydata_windfarm_example]; must be in EPSG 2154
#' @param timeRange minimum and maximum years of bird counts to integrate in trends estimates
#' @param foraging_ranges table of species foraging ranges, for instance using [bdydata_foraging_ranges]
#' @param mortality table of collision, for instance using [bdydata_mortality_example]
#' @param n_iteration Number of iterations for the collision model (parameter used in [bdy_process_mortality()])
#' @param ni_noImpact Number of iterations to use in [bdy_model_no_impact()]
#' @param ni_withImpact Number of iterations to use in [bdy_model_with_impact()]
#'
#' @returns List of model outputs for each species, including:
#'          \itemize{
#'          \item 'no_impact': Bird counts assuming no impact from windfarm; see [bdy_model_no_impact()] for more details
#'          \item 'with_impact': Bird counts assuming no impact from windfarm; see [bdy_model_with_impact()] for more details
#'          \item 'colonies': Table of colonies where the species is present
#'          \item 'mortality': Table of collision per group of colony (see more details in [bdy_process_mortality()])
#'          \item 'distance': Table of distances between colonies and windfarms
#'          \item 'distanceBIN': Binarised table of distances between colonies and windfarms binarised, reporting whether the windfarm is accessible to the individuals from the colony (i.e., the distance is smaller than the maximum foraging distance of the species)
#'          }
#' @export
#'

bdy_run_analysis <- function(species,countData,windfarms,timeRange = c(2009,2021), foraging_ranges, mortality,
                             n_iteration=1000,ni_noImpact=10000,ni_withImpact=1000,...){

  ##checking count data and adding colonies code
  colonies_all <- bdy_summarise_colonies(countData)

  #calculating cost matrix for each seafront
  cost_matrix <- bdy_get_cost_raster(colonies_all)

  #calculating distances
  all_distances <- bdy_get_distances(colonies=colonies_all,
                                     windfarms = windfarms,
                                     costMatrix=cost_matrix,
                                     doShpa = any(!subset(foraging_ranges, species_latin %in% species)$terrestrial_habits))

  #preparing a list to store model output for each species
  model_output <- list()

  for(sp in species){

    # !! prevoir cas ou l'espece d'interet n'est pas dans la liste !!
    foraging_range_sp = foraging_ranges$max_km[which(foraging_ranges$species_latin==sp)]
    terrestrial_habit = foraging_ranges$terrestrial_habits[which(foraging_ranges$species_latin==sp)]

    #processing count data
    count_processed <- bdy_process_count_data(sp=sp,
                                              countData=countData,
                                              colonies=colonies_all,first_year=timeRange[1],last_year = timeRange[2],
                                              max_foraging_range_km=foraging_range_sp)

    #selecting apropriate distance table (shpa or eucl), and only for colonies where species is present
    if(terrestrial_habit){
      distances <-  all_distances$eucl_dist
    }else{
      distances <-  all_distances$shpa_dist
    }
    distances <- distances[rownames(distances) %in% count_processed$colonies_sp$colony_code,]

    #apportionning
    apportionning <- bdy_apportionning(max_foraging_range_km=foraging_range_sp,
                                       colonies=count_processed$colonies_sp,
                                       sea_area=count_processed$sea_area_sp,
                                       tbl_dist = distances)

    #distributing mortality
    morta_iter_group <- bdy_process_mortality(
      collision = mortality[which(mortality$species_latin==sp),],
      season=seasons[which(seasons$species_latin==sp),][month.abb],
      n_iteration=n_iteration,
      RW_group = apportionning$RW_group
    )

    #model without impact
    no_impact_output <- bdy_model_no_impact(group_counts=count_processed$group_counts_sp,
                                            ppa=count_processed$ppa_yr_sp,
                                            survival=vital_rates[vital_rates$species_latin== sp, "survival"],
                                            fecundity=vital_rates[vital_rates$species_latin == sp, "fecundity"],
                                            propRepro=vital_rates[vital_rates$species_latin == sp, "propRepro"],
                                            modelFile=modelFilePath, nimble=T,
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

  return(model_output)
}
