#' Create plot with windfarm impact by parc
#'
#' Plot impact of each windfarm on each species
#'
#' @param Raw_ResTables Raw result tables from function [bdy_raw_res_tables()]
#' @param formatted_mortality formatted table of collision, for instance using [bdydata_mortality_example]
#'
#' @returns list with a ggplot object and a table of raw results
#'
#' @export

bdy_plot_windfarm_impact <- function(Raw_ResTables, formatted_mortality){

  # Fix species french / latin depending on app or package
  if(any(Raw_ResTables$Simulated_National$Species %in% formatted_mortality$species_latin)){
    formatted_mortality$Species_col <- formatted_mortality$species_latin
  } else {
    formatted_mortality$Species_col <- formatted_mortality$espece
  }

  # Calculate mean mortality per species and windfarm
  moyenne <- formatted_mortality %>%
    dplyr::group_by(Species_col, windfarm) %>%
    dplyr::summarise(Mean_mortality=round(mean(coefficient, na.rm=T),4)) %>%
    group_by(Species_col) %>%
    mutate(Prop_morta = round(Mean_mortality/sum(Mean_mortality,na.rm=T),4))

  # Split relative impact by mortality proportions
  moyenne$IR <- Raw_ResTables$Tableau_National$RelImpact_95[match(moyenne$Species_col, Raw_ResTables$Tableau_National$Species)]
  moyenne$Impact_windfarm <- round(moyenne$IR * moyenne$Prop_morta, 4)
  moyenne$IR <- NULL

  # Plot
  G <- ggplot(moyenne)+
    geom_bar(aes(x=Species_col, y=Impact_windfarm, fill=windfarm), stat="identity", width=min(c(0.8), 0.1*(nrow(Raw_ResTables$Tableau_National)+1)))+
    scale_fill_brewer(palette="Set3")+
    ylab("Impact relatif par parc (%)")+xlab("")+
    theme_minimal()

  return(list(plot=G, table=moyenne))

}
