#' Pretty result table
#'
#' Transform raw result tables to a table ready to be used for reporting
#'
#' @param Raw_ResTables Raw result tables from function [bdy_raw_res_tables()]
#' @param type Character, either "national" or "subpop", indicating what table to transform
#'
#' @returns A clean table with initial population, mortality, relative impact, extinction probability and their confidence intervals
#' @export
#'

bdy_pretty_result_table <- function(Raw_ResTables, type){

  # Select table depending on type (national or subpop)
  if(type=="national"){
    ResTable <- Raw_ResTables$Tableau_National
  } else {
    ResTable <- Raw_ResTables$Tableau_Subpop
  }

  # Format columns
  ResTable <- ResTable %>%
    mutate(
      Population_initiale = paste0(round(Pop_Init_med), " [", round(Pop_Init_2.5), "-", round(Pop_Init_97.5), "]"),
      Nombre_collisions = paste0(round(Mortality_med, 2), " [", round(Mortality_2.5, 2), "-", round(Mortality_97.5, 2), "]"),
      Impact_Relatif = paste0(round(RelImpact_med, 4), " [", round(RelImpact_2.5, 4), "-", round(RelImpact_97.5, 4), "]"),
      Augmentation_Extinction=round(Ext_Relative,4)
    ) %>%
    dplyr::rename(Espece=Species) %>%
    subset(., select=names(.)[grepl("_med", names(.))==F & grepl("_2.5", names(.))==F & grepl("_97.5", names(.))==F & grepl("Ext_", names(.))==F])

  # Add max impact
  if(type=="national"){
    ResTable_MaxIR <- Raw_ResTables$Tableau_Subpop %>% group_by(Species) %>%  slice(which.max(RelImpact_med))
    ResTable_MaxIR$Text <- ifelse(ResTable_MaxIR$RelImpact_med==0, 0, paste0(round(ResTable_MaxIR$RelImpact_med,4), " (Groupe de colonies : ", ResTable_MaxIR$GrColo, ")"))
    ResTable_MaxExt <- Raw_ResTables$Tableau_Subpop %>% group_by(Species) %>%  slice(which.max(Ext_Relative))
    ResTable_MaxExt$Text <- ifelse(ResTable_MaxExt$Ext_Relative==0, 0, paste0(round(ResTable_MaxExt$Ext_Relative,4), " (Groupe de colonies : ", ResTable_MaxExt$GrColo, ")"))
    ResTable$MAX_Impact_relatif <- ResTable_MaxIR$Text[match(ResTable$Espece, ResTable_MaxIR$Species)]
    ResTable$MAX_Augmentation_extinction <- ResTable_MaxExt$Text[match(ResTable$Espece, ResTable_MaxExt$Species)]
  }


  # Return
  return(ResTable)
}
