#' Check data conformity
#'
#' Checks that all data required for the analyses are loaded; see Vignette for more details.
#'
#' @param windfarms_L93 sf object indicating the position of windfarms, for instance using [bdydata_parcs_example]; must be in EPSG 2154
#' @param formatted_mortality table of collision from [bdy_check_mortality()], for instance using [bdydata_mortality_example]
#' @param colonies_all table of colonies from [bdy_check_colonies()]
#' @param countData table of bird counts
#' @param vital_rates table of species vital rates, for instance using [bdydata_vital_rates]
#' @param seasons table of species seasonal activities, for instance using [bdydata_seasons]
#' @param foraging_ranges table of species foraging ranges, for instance using [bdydata_foraging_ranges]
#' @param countryShape sf object with the shape of countries shoreline, for instance cropped from [bdydata_world_map]
#'
#' @returns a text message if all data are valid, an error if any problem was detected
#' @export
#'

bdy_check_alldata <- function(windfarms_L93, formatted_mortality, colonies_all, countData, vital_rates, seasons, foraging_ranges, countryShape){

  ### Check column names
  # windfarms_L93
  col_missing_windfarms <- c("NAME", "seafront") %>% .[! . %in% names(windfarms_L93)]
  if(length(col_missing_windfarms)>0){stop(paste0("Missing columns in 'windfarms_L93': ", paste0(col_missing_windfarms, collapse=", ")))}

  # formatted_mortality
  col_missing_morta <- c("species_latin", "parc", "mois", "iteration", "coefficient") %>% .[! . %in% names(formatted_mortality)]
  if(length(col_missing_morta)>0){stop(paste0("Missing columns in 'formatted_mortality': ", paste0(col_missing_morta, collapse=", ")))}

  # colonies_all
  col_missing_colonies <- c("seafront", "colony_code", "colony") %>% .[! . %in% names(colonies_all)]
  if(length(col_missing_colonies)>0){stop(paste0("Missing columns in 'colonies_all': ", paste0(col_missing_colonies, collapse=", ")))}

  # countData
  col_missing_count <- c("species_latin", "seafront", "year", "count_mean", "colony", "colony_code") %>% .[! . %in% names(countData)]
  if(length(col_missing_count)>0){stop(paste0("Missing columns in 'countData': ", paste0(col_missing_count, collapse=", ")))}

  # vital_rates
  col_missing_vital <- c("species_latin", "class", "age", "survival", "fecundity", "propRepro") %>% .[! . %in% names(vital_rates)]
  if(length(col_missing_vital)>0){stop(paste0("Missing columns in 'vital_rates': ", paste0(col_missing_vital, collapse=", ")))}

  # seasons
  col_missing_seasons <- c("species_latin", month.abb) %>% .[! . %in% names(seasons)]
  if(length(col_missing_seasons)>0){stop(paste0("Missing columns in 'seasons': ", paste0(col_missing_seasons, collapse=", ")))}

  # foraging_ranges
  col_missing_foraging <- c("species_latin", "terrestrial_habits", "max_km") %>% .[! . %in% names(foraging_ranges)]
  if(length(col_missing_foraging)>0){stop(paste0("Missing columns in 'foraging_ranges': ", paste0(col_missing_foraging, collapse=", ")))}



  ### Check windfarm data
  # Check names match between windfarms_L93 and formatted_mortality
  Farms_shp <- unique(windfarms_L93$NAME)
  Farms_morta <- unique(formatted_mortality$parc)

  Farms_shp_missing <- unique(Farms_shp[! Farms_shp %in% Farms_morta])
  if(length(Farms_shp_missing)>0){stop(paste0("Some windfarms are mentionned in 'windfarms_L93' but not in 'formatted_mortality': ", paste0(Farms_shp_missing, collapse=", ")))}
  Farms_morta_missing <- unique(Farms_morta[! Farms_morta %in% Farms_shp])
  if(length(Farms_morta_missing)>0){stop(paste0("Some windfarms are mentionned in 'formatted_mortality' but not in 'windfarms_L93': ", paste0(Farms_morta_missing, collapse=", ")))}

  # Check all windfarms have names
  if((TRUE %in% is.na(windfarms_L93$NAME)) | ("" %in% windfarms_L93$NAME)){stop("Some windfarm names are missing from 'windfarms_L93$NAME'")}

  # Check no duplicated names
  Dupli_names <- windfarms_L93$NAME %>% table(.) %>% subset(., .>1) %>% names(.)
  if(length(Dupli_names)>0){stop(paste0("Some windfarm names are duplicated: ", paste0(Dupli_names, collapse=", "), ". Each windfarm should have a different name"))}



  ### Check species names match between formatted_mortality and all species data
  SP_list <- unique(formatted_mortality$species_latin)

  # countData
  SP_missing_count <- SP_list[! SP_list %in% countData$species_latin]
  if(length(SP_missing_count)>0){stop(paste0("Some species are included in 'formatted_mortality' but are missing from 'countData': ", paste0(SP_missing_count, collapse=", ")))}

  # vital_rates
  SP_missing_vital <- SP_list[! SP_list %in% vital_rates$species_latin]
  if(length(SP_missing_vital)>0){stop(paste0("Some species are included in 'formatted_mortality' but are missing from 'vital_rates': ", paste0(SP_missing_vital, collapse=", ")))}

  # seasons
  SP_missing_seasons <- SP_list[! SP_list %in% seasons$species_latin]
  if(length(SP_missing_seasons)>0){stop(paste0("Some species are included in 'formatted_mortality' but are missing from 'seasons': ", paste0(SP_missing_seasons, collapse=", ")))}

  # foraging_ranges
  SP_missing_foraging <- SP_list[! SP_list %in% foraging_ranges$species_latin]
  if(length(SP_missing_foraging)>0){stop(paste0("Some species are included in 'formatted_mortality' but are missing from 'foraging_ranges': ", paste0(SP_missing_foraging, collapse=", ")))}



  ### Check formatted_mortality columns

  # Check that month is numeric, with all values, only 1 to 12
  if(all(is.na(formatted_mortality$mois))){stop("The column 'mois' from 'formatted_mortality' is not numeric.")}
  if(any(!(1:12) %in% formatted_mortality$mois)){stop("The column 'mois' from 'formatted_mortality' does not include all month values (1 to 12 are required).")}
  if(any(!formatted_mortality$mois %in% 1:12)){stop("The column 'mois' from 'formatted_mortality' includes values other than 1:12.")}

  ### check that coefficient is numeric
  if(all(is.na(formatted_mortality$coefficient))){stop("The column 'coefficient' from 'formatted_mortality' is not numeric.")}
  if(any(formatted_mortality$coefficient < 0)){stop("The column 'coefficient' from 'formatted_mortality' has some negative values; they should all be positive.")}

  ### check that iteration is numeric
  if(all(is.na(formatted_mortality$iteration))){stop("The column 'iteration' from 'formatted_mortality' is not numeric.")}



  ### Check colonies




  return("All required datasets are loaded and valid")
}
