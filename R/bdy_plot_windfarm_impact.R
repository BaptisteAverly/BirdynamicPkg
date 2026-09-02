#' Windfarm impact
#'
#' Plot impact of each windfarm on each species
#'
#' @param Raw_ResTables Raw result tables from function [bdy_raw_res_tables()]
#' @param formatted_mortality formatted table of collision, for instance using [bdydata_mortality_example]
#'
#' @returns ggplot object
#'
#' @export

bdy_plot_windfarm_impact <- function(Raw_ResTables, formatted_mortality){

  moyenne <- formatted_mortality %>%
    dplyr::group_by(species_latin, windfarm) %>%
    dplyr::summarise(Mean=mean(coefficient, na.rm=T)) %>%
    group_by(species_latin) %>%
    mutate(Prop_morta = Mean/sum(Mean,na.rm=T))

  moyenne$IR <- Raw_ResTables$Tableau_National$RelImpact_med[match(moyenne$species_latin, Raw_ResTables$Tableau_National$Species)]
  moyenne$Impact_parc <- moyenne$IR * moyenne$Prop_morta

  G <- ggplot(moyenne)+
    geom_bar(aes(x=species_latin, y=Impact_parc, fill=windfarm), stat="identity")+
    scale_fill_brewer(palette="Set3")+
    ylab("Impact relatif par parc (%)")+xlab("")+
    theme_minimal()

  return(G)

}
