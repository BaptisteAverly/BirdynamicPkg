#' Add additional count data
#'
#' Used to add to or replace the default count data included in the package.
#'
#' @param new_count_table data frame containing the count data to be added to the default one. Must contain the following columns:
#'                  \itemize{
#'                  \item 'colony': character, unique identifier for the colony in which the count was done
#'                  \item 'lon': numeric, longitude of the colony (WGS84)
#'                  \item 'lat': numeric, latitude of the colony (WGS84)
#'                  \item 'seafront': character, which sea front is the colony located in
#'                  \item 'species_latin': character, latin name of the species counted. For a list of available species, see bdydata_seasons$species_latin
#'                  \item 'year': numeric, year during which the count was done
#'                  \item 'count': numeric, number of breeding pair counted for a given colony, species, and year
#'                  }
#' @param removeNewDuplicates boolean. If TRUE (default), then potential duplicated rows (same longitude, latitude, species, and count year) will be removed from the count table provided by the user (argument new_count_table).
#'                                     if FALSE, then potential duplicated rows will be removed from the default count table (bdydata_colonies_counts_low_res).
#' @param replaceDefaultTable boolean. If FALSE (default), then provided table (argument new_count_table), will be appended to the default count table (bdydata_colonies_counts_low_res).
#'                                     If TRUE, then provided table will completely remplace the default table.
#'
#' @returns Updated data frame with formatted count data
#' @export
#'
#' @examples
#'
bdy_add_colonies <- function(new_count_table, removeNewDuplicates = T, replaceDefaultTable = F){

  ## check for missing columns

  columns <- c("colony","lon","lat","seafront","species_latin","year","count")

  missingCol <- columns[which(!columns %in% colnames(new_count_table))]

  if(length(missingCol) > 0){
    stop(paste0("Column(s) missing:",paste0("'",missingCol,"'",collapse=", ")))
  }

  ## check that species are legal

  wrongSpeciesBool <- !new_count_table$species_latin %in% bdydata_seasons$species_latin

  if(all(wrongSpeciesBool)){

    stop("No species in the table is currently implemented in the package (or all latin names are misspelled).
         For a list of available species, see bdydata_seasons$species_latin")
  }

  wrongSpecies <- unique(new_count_table$species_latin[which(wrongSpeciesBool)])

  if(length(wrongSpecies) > 0){
    warning(paste0("The following species are currently not implemented in the package and were removed from the table (",
                    length(which(wrongSpeciesBool))," rows):", paste0("'",wrongSpecies,"'",collapse=", "),
                   "
                   For a list of available species, see bdydata_seasons$species_latin"))
  }

  if("species_fr" %in% colnames(new_count_table)){
    new_count_table$species_fr <- bdydata_seasons$species_fr[match(new_count_table$species_latin,bdydata_seasons$species_latin)]
  }

  if("species_en" %in% colnames(new_count_table)){
    new_count_table$species_en <- bdydata_seasons$species_en[match(new_count_table$species_latin,bdydata_seasons$species_latin)]
  }


  ## check that numeric columns are numeric

  numCol <- c("lon","lat","year","count")

  notNum <- numCol[which(!apply(new_count_table[,numCol],2,is.numeric))]

  if(length(notNum) > 0){
    stop(paste0("Following column(s) should be numeric:",paste0("'",notNum,"'",collapse=", ")))
  }

  ## removing unnecessary columns

  extraCol <- which(!colnames(new_count_table) %in% colnames(bdydata_colonies_counts_low_res))
  if(length(extraCol) > 0){
    new_count_table <- new_count_table[,-extraCol]
  }

  ## updating the count table

  if(replaceDefaultTable){

    new_data <- new_count_table

  }else{

    if(removeNewDuplicates){
      new_data <- rbind(bdydata_colonies_counts_low_res,new_count_table)
    }else{
      new_data <- rbind(new_count_table,bdydata_colonies_counts_low_res)
    }
  }

  ## checking for duplicates

  temp <- new_data
  temp$lon <- round(temp$lon,3)
  temp$lat <- round(temp$lat,3)

  dupl <- which(duplicated(temp[,c("lon","lat","species_latin","year")]))

  if(length(dupl) > 0 ){
    newData <- newData[-dupl,]
    warning(paste0(length(dupl)," rows in the final table were duplicates of other rows (same longitude, latitude, species, and count year) and thus were deleted"))
  }

  return(newData)
}
