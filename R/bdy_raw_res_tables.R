#' Raw result tables
#'
#' Compute result tables (national and by group of colonies)
#'
#' @param mod_out Output of the [bdy_model_with_impact()] function, completed with colonies, mortality and distance data
#'
#' @returns Three tables: \cr
#'          - Simulated_National: Simulated counts for the 30 projected years at national level (used to plot temporal trends with [bdy_summary_figure()])
#'          - Tableau_Subpop: Summary table of results by colony (with relative impact, increase in probability of extinction, mortality, etc)
#'          - Tableau_National: Summary table of results at national level (with relative impact, increase in probability of extinction, mortality, etc)
#' @export
#'

bdy_raw_res_tables <- function(mod_out){

  ### Merge simulated data of all species
  Simulated_Counts <- data.frame()
  for(SP in 1:length(mod_out)){

    Simulated_SP <- mod_out[[SP]]$no_impact %>%
      reshape2::melt(., value.name="Count_noimpact") %>%
      dplyr::rename(., GrColo=Var1, Year=Var2, Iteration=Var3) %>%
      mutate(Species=names(mod_out)[SP]) %>%
      relocate(., "Species", .before="GrColo")

    Simulated_SP$Count_withimpact <- mod_out[[SP]]$with_impact %>%
      reshape2::melt(., value.name="Count") %>%
      .$Count

    Simulated_Counts <- rbind(Simulated_Counts, Simulated_SP)

  }

  ### Get parameters
  ny_proj <- max(Simulated_Counts$Year) # Number of years


  ### Impact per group of colonies -----

  ## Relative Impact (the last year of simulation)
  Simulated_Counts$PopInit <- subset(Simulated_Counts, Year==1)$Count_noimpact[match(paste0(Simulated_Counts$Species, Simulated_Counts$GrColo, Simulated_Counts$Iteration), paste0(subset(Simulated_Counts, Year==1)$Species, subset(Simulated_Counts, Year==1)$GrColo, subset(Simulated_Counts, Year==1)$Iteration))]
  Relative_Impact <- Simulated_Counts %>%
    subset(., Year==ny_proj) %>%
    mutate(Rel_impact=100*(Count_noimpact-Count_withimpact)/Count_noimpact,
           Trend_no=(Count_noimpact-PopInit)/PopInit,
           Trend_impact=(Count_withimpact-PopInit)/PopInit,
    )

  ## Subpopulation table
  Tableau_Subpop <- Relative_Impact %>%
    dplyr::group_by(Species, GrColo) %>%
    dplyr::summarise(Pop_Init_med=round(median(PopInit, na.rm=T)),
                     Pop_Init_2.5=round(quantile(PopInit, probs=0.025, na.rm=T)),
                     Pop_Init_97.5=round(quantile(PopInit, probs=0.975, na.rm=T)),

                     Trend_no_med=median(Trend_no, na.rm=T),
                     Trend_no_2.5=quantile(Trend_no, probs=0.025, na.rm=T),
                     Trend_no_97.5=quantile(Trend_no, probs=0.975, na.rm=T),

                     Trend_impact_med=median(Trend_impact, na.rm=T),
                     Trend_impact_2.5=quantile(Trend_impact, probs=0.025, na.rm=T),
                     Trend_impact_97.5=quantile(Trend_impact, probs=0.975, na.rm=T),

                     Mortality_med=NA, Mortality_2.5=NA, Mortality_97.5=NA, Parc_proche=NA,

                     RelImpact_med=median(Rel_impact, na.rm=T),
                     RelImpact_2.5=quantile(Rel_impact, probs=0.025, na.rm=T),
                     RelImpact_97.5=quantile(Rel_impact, probs=0.975, na.rm=T),

                     Ext_Noimpact=mean(Count_noimpact==0, na.rm=T),
                     Ext_Withimpact=mean(Count_withimpact==0, na.rm=T),
                     .groups="keep"
    ) %>%
    mutate(Ext_Relative = 100*(Ext_Withimpact-Ext_Noimpact) %>% ifelse(is.na(.), 0, .))



  ## Add mortality and distance to parcs from other table
  for(SP in 1:length(mod_out)){

    # Mortality
    Morta <- as.data.frame(t(apply(mod_out[[SP]]$mortality, 2, quantile, probs = c(0.025, 0.5, 0.975))))

    Tableau_Subpop$Mortality_med[Tableau_Subpop$Species==names(mod_out)[SP]] <- Morta$`50%`
    Tableau_Subpop$Mortality_2.5[Tableau_Subpop$Species==names(mod_out)[SP]] <- Morta$`2.5%`
    Tableau_Subpop$Mortality_97.5[Tableau_Subpop$Species==names(mod_out)[SP]] <- Morta$`97.5%`

    # Distance to parcs
    Dist_qtt <- mod_out[[SP]]$distance %>%
      as.data.frame() %>%
      mutate(
        Dist_min = apply(., 1, min, na.rm=T),
        Parc_min = names(.)[apply(., 1, which.min)]
      )

    Colonies_SP <- mod_out[[SP]]$colonies
    Colonies_SP$Dist_min <- Dist_qtt$Dist_min[match(Colonies_SP$code_colonie, rownames(Dist_qtt))]
    Colonies_SP$Parc_min <- Dist_qtt$Parc_min[match(Colonies_SP$code_colonie, rownames(Dist_qtt))]

    GrColo_SP <- Colonies_SP %>%
      group_by(group) %>%
      summarise(colonies=paste0(unique(code_colonie), collapse=", "),
                Dist_min=min(Dist_min, na.rm=T),
                Parc_min=Parc_min[which.min(Dist_min)]
      ) %>%
      mutate(Parc_proche = paste0(Parc_min, " (", Dist_min, "km)"))

    Tableau_Subpop$Parc_proche[Tableau_Subpop$Species==names(mod_out)[SP]] <- GrColo_SP$Parc_proche[match(Tableau_Subpop$GrColo[Tableau_Subpop$Species==names(mod_out)[SP]], GrColo_SP$group)]

  }


  ### Impact National ------
  Simulated_National <- Simulated_Counts %>%
    dplyr::group_by(Species, Year, Iteration) %>%
    dplyr::summarise(Sum_PopInit=sum(PopInit, na.rm=T),
                     Sum_noimpact=sum(Count_noimpact, na.rm=T),
                     Sum_withimpact=sum(Count_withimpact, na.rm=T),
                     .groups="keep")

  ## Relative Impact (the last year of simulation)
  Relative_Impact_National <- Simulated_National %>%
    subset(., Year==ny_proj) %>%
    mutate(Rel_impact=100*(Sum_noimpact-Sum_withimpact)/Sum_noimpact)


  ## Probability of extinction (National)
  Tableau_National <- Relative_Impact_National %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(
      Pop_Init_med=round(median(Sum_PopInit, na.rm=T)),
      Pop_Init_2.5=round(quantile(Sum_PopInit, probs=0.025, na.rm=T)),
      Pop_Init_97.5=round(quantile(Sum_PopInit, probs=0.975, na.rm=T)),

      Mortality_med=NA, Mortality_2.5=NA, Mortality_97.5=NA,

      RelImpact_med = median(Rel_impact, na.rm=T),
      RelImpact_2.5 = quantile(Rel_impact, probs=0.025, na.rm=T),
      RelImpact_97.5 = quantile(Rel_impact, probs=0.975, na.rm=T),

      Ext_Noimpact=mean(Sum_noimpact==0, na.rm=T),
      Ext_Withimpact=mean(Sum_withimpact==0, na.rm=T),
      .groups="keep"
    ) %>%
    mutate(Ext_Relative = 100*(Ext_Withimpact-Ext_Noimpact) %>% ifelse(is.na(.), 0, .))


  # Add mortality from another table
  for(SP in 1:length(mod_out)){

    Morta <- quantile(apply(mod_out[[SP]]$mortality, 1, sum), probs = c(0.025, 0.5, 0.975))

    Tableau_National$Mortality_med[Tableau_National$Species==names(mod_out)[SP]] <- Morta["50%"]
    Tableau_National$Mortality_2.5[Tableau_National$Species==names(mod_out)[SP]] <- Morta["2.5%"]
    Tableau_National$Mortality_97.5[Tableau_National$Species==names(mod_out)[SP]] <- Morta["97.5%"]

  }

  return(list(Simulated_National=Simulated_National,
              Tableau_Subpop=Tableau_Subpop,
              Tableau_National=Tableau_National))

}
