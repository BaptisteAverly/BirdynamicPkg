#' Summary trends figure
#'
#' Plots figure of projected national trends
#'
#' @param Raw_ResTables Raw result tables from function [bdy_raw_res_tables()]
#' @param Pretty_table_national Pretty table from function [bdy_pretty_result_table()]
#'
#' @returns A ggplot object with a facet plot per species showing trends with and without impact
#' @export
#'

bdy_summary_figure <- function(Raw_ResTables, Pretty_table_national){

  # Prepare trends plot
  Summary_simulated <- ddply(Raw_ResTables$Simulated_National, .(Year, Species), function(x){data.frame(
    Count_noimpact_mean = mean(x$Sum_noimpact, na.rm=T),
    Count_withimpact_mean = mean(x$Sum_withimpact, na.rm=T),
    Count_noimpact_2.5 = quantile(x$Sum_noimpact, probs=0.025, na.rm=T),
    Count_withimpact_2.5 = quantile(x$Sum_withimpact, probs=0.025, na.rm=T),
    Count_noimpact_97.5 = quantile(x$Sum_noimpact, probs=0.975, na.rm=T),
    Count_withimpact_97.5 = quantile(x$Sum_withimpact, probs=0.975, na.rm=T),
    Text=NA
  )})

  # Add text info (relative impact + extinction)
  for(SP in Pretty_table_national$Espece){
    Summary_simulated$Text[Summary_simulated$Species==SP] <- paste0("<b>Impact relatif : </b>", Pretty_table_national$Impact_Relatif[Pretty_table_national$Espece==SP], "\n<b>Augmentation extinction : </b>", Pretty_table_national$Augmentation_Extinction[Pretty_table_national$Espece==SP])
  }

  # Plot
  G <- ggplot(Summary_simulated)+
    geom_ribbon(aes(x=Year, ymin=Count_noimpact_2.5/Count_noimpact_mean*100, ymax=Count_noimpact_97.5/Count_noimpact_mean*100), fill="#4dac26", alpha=0.2)+
    geom_ribbon(aes(x=Year, ymin=Count_withimpact_2.5/Count_noimpact_mean*100, ymax=0.001+Count_withimpact_97.5/Count_noimpact_mean*100), fill="#c51b7d", alpha=0.2)+
    geom_line(aes(x=Year, y=Count_noimpact_mean/Count_noimpact_mean*100, col="Sans parcs éoliens"), linewidth=1.5)+
    geom_line(aes(x=Year, y=Count_withimpact_mean/Count_noimpact_mean*100, text=Text, col="Avec parcs éoliens"), linetype="dashed", linewidth=1.5)+
    facet_wrap(~Species)+
    scale_color_manual(name="Modèle", breaks=c("Sans parcs éoliens", "Avec parcs éoliens"), values=c("Sans parcs éoliens"="#4dac26", "Avec parcs éoliens"="#c51b7d"))+
    coord_cartesian(ylim=c(0, min(c(150, max(Summary_simulated$Count_withimpact_97.5/Summary_simulated$Count_noimpact_mean*100)))))+
    xlab("Année")+ylab("Proportion de la population sans impact (%)")+
    theme_minimal()

  return(G)
}
