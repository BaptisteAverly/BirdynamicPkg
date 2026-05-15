#' Windfarm impact
#'
#' Plot impact of each windfarm on each species
#'
#' @param Raw_ResTables Raw result tables from function [bdy_raw_res_tables()]
#' @param mortalities Mortality data formatted with function [bdy_check_mortality()]
#'
#' @returns ggplot object
#'

bdy_plot_windfarm_impact <- function(Raw_ResTables, mortalities){

  moyenne <- mortalities %>%
    dplyr::group_by(species_latin, parc) %>%
    dplyr::summarise(Mean=mean(coefficient, na.rm=T)) %>%
    group_by(species_latin) %>%
    mutate(Prop_morta = Mean/sum(Mean,na.rm=T))

  moyenne$IR <- Raw_ResTables$Tableau_National$RelImpact_med[match(moyenne$species_latin, Raw_ResTables$Tableau_National$Species)]
  moyenne$Impact_parc <- moyenne$IR * moyenne$Prop_morta

  G <- ggplot(moyenne)+
    geom_bar(aes(x=species_latin, y=Impact_parc, fill=parc), stat="identity")+
    scale_fill_brewer(palette="Set3")+
    ylab("Impact relatif par parc (%)")+xlab("")+
    theme_minimal()

  return(G)

}
