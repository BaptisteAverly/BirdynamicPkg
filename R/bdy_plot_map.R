#' Windfarm impact
#'
#' Plot impact of each windfarm on each species
#'
#' @param mod_out Raw outputs of Bird dynamic models
#' @param Raw_ResTables Raw result tables from function [bdy_raw_res_tables()]
#' @param parcs sf object of windfarms
#' @param BufferForaging sf object of species buffers of foraging around colonies to be plotted on the interactive map
#'
#' @returns htmlwidget (leaflet)
#'

bdy_plot_map <- function(mod_out, Raw_ResTables, parcs, BufferForaging=NULL){

  ### Create empty map
  Leaf <- leaflet() %>%
    addProviderTiles("Esri.OceanBasemap", options = list(maxZoom=9)) %>%
    fitBounds(-6, 42, 8.6, 51.5) %>%
    addLayersControl(baseGroups=c("OpenStreetMap"), position="bottomleft") %>%
    addScaleBar(position="bottomright") %>%
    addMouseCoordinates()

  ### Add colonies
  for(SP in 1:length(mod_out)){

    # Get if colony in range of a parc
    Dist <- mod_out[[SP]]$distanceBIN %>% data.frame()
    Dist$Parcs_foraging <- apply(Dist, 1, function(x) paste0(names(data.frame(mod_out[[SP]]$distance))[x == T], collapse = ";"))

    # Distance to parcs
    Dist_qtt <- mod_out[[SP]]$distance %>%
      data.frame() %>%
      mutate(
        Dist_min = apply(., 1, min, na.rm=T),
        Parc_min = names(.)[apply(., 1, which.min)]
      )

    # Load colonies
    Colonies_SP <- mod_out[[SP]]$colonies %>%
      st_transform(., st_crs(4326)) %>%
      mutate(Lon=st_coordinates(.)[,1], Lat=st_coordinates(.)[,2]) %>%
      mutate(Species=names(mod_out)[SP])
    Colonies_SP$Popup <- Colonies_SP$colony
    Colonies_SP$Parcs_foraging <- Dist$Parcs_foraging[match(Colonies_SP$colony_code, rownames(Dist))]
    Colonies_SP$Dist_min <- Dist_qtt$Dist_min[match(Colonies_SP$colony_code, rownames(Dist_qtt))]
    Colonies_SP$Parc_min <- Dist_qtt$Parc_min[match(Colonies_SP$colony_code, rownames(Dist_qtt))]

    GrColo_SP_MCP <- Colonies_SP %>%
      st_buffer(., 1) %>%
      group_by(group) %>%
      summarise(geometry=st_union(geometry),
                colonies=paste0(unique(Popup), collapse=", "),
                Dist_min=min(Dist_min, na.rm=T),
                Parc_min=Parc_min[which.min(Dist_min)],
                Parcs_foraging=Parcs_foraging %>% strsplit(., ";") %>% unlist(.) %>% unique(.) %>% sort(.) %>% paste0(., collapse=", ") %>% ifelse(.=="", "aucun", .)
      ) %>%
      st_make_valid() %>%
      st_convex_hull()
    GrColo_SP <- st_centroid(GrColo_SP_MCP)

    # Add Relative Impact and PopInit
    Raw_ResSP <- subset(Raw_ResTables$Tableau_Subpop, Species==names(mod_out)[SP])
    GrColo_SP$RelImpact <- Raw_ResSP$RelImpact_med[match(GrColo_SP$group, Raw_ResSP$GrColo)]
    GrColo_SP$PopInit <- Raw_ResSP$Pop_Init_med[match(GrColo_SP$group, Raw_ResSP$GrColo)]


    # Add popup
    Res_match <- Raw_ResSP[match(GrColo_SP$group, Raw_ResSP$GrColo),]
    GrColo_SP$Popup <- paste0(
      "<h2>", names(mod_out)[SP], "</h2>",
      "<b>Groupe de colonies n°", GrColo_SP$group, "</b><br>",
      "<b>Colonies : </b>", GrColo_SP$colonies, "<br>",
      "<b>Parc le plus proche : </b>", GrColo_SP$Parc_min, " (", GrColo_SP$Dist_min, "km)<br>",
      "<b>Parcs dans zone d'alimentation: </b>", GrColo_SP$Parcs_foraging, "<br><br>",

      "<p style='color:FireBrick'><b>Impact relatif sur les tendances : </b>", round(Res_match$RelImpact_med,2), "% [", round(Res_match$RelImpact_2.5,2), " ; ", round(Res_match$RelImpact_97.5,2), "]<br>",
      "<b>Augmentation de la probabilité d'extinction : </b>", round(Res_match$Ext_Relative,2), "%</p>",

      "<b>Population en 2021 : </b>", Res_match$Pop_Init_med, " [", Res_match$Pop_Init_2.5, " ; ", Res_match$Pop_Init_97.5, "]<br>",
      "<b>Nombre de collisions : </b>", round(Res_match$Mortality_med,2), " [", round(Res_match$Mortality_2.5,2), " ; ", round(Res_match$Mortality_97.5,2), "]<br>", # Verifier que ca correspond bien a "Mortality"

      "<b>Tendance sans parc : </b>", 100*round(Res_match$Trend_no_med,3), "% [", 100*round(Res_match$Trend_no_2.5,3), " ; ", 100*round(Res_match$Trend_no_97.5,3), "]<br>",
      "<b>Tendance avec parc : </b>", 100*round(Res_match$Trend_impact_med,3), "% [", 100*round(Res_match$Trend_impact_2.5,3), " ; ", 100*round(Res_match$Trend_impact_97.5,3), "]<br>",

      "<b>Probabilité d'extinction sans parc : </b>", 100*round(Res_match$Ext_Noimpact,3), "%<br>",
      "<b>Probabilité d'extinction avec parc : </b>", 100*round(Res_match$Ext_Withimpact,3), "%<br>",
      "<i>Les tendances estimées pour certaines espèces sont sujettes à caution, cela n'empêche pas les impacts relatifs d'être robustes"
    )

    # Prepare Buffer
    if(is.null(BufferForaging)==F){
      Buff_toplot <- BufferForaging[BufferForaging$Espece==names(mod_out)[SP],] %>% st_transform(., st_crs(GrColo_SP_MCP))
      if(nrow(Buff_toplot)>0){
        Leaf <- Leaf %>%
          addPolygons(data=Buff_toplot, fill=F, fillOpacity=0, stroke=T, weight=2, dashArray="5, 5", color="#8c2d04", group=names(mod_out)[SP])
      }
    }

    # Prepare color and size
    ColPal <- colorNumeric("Reds", domain=c(0,100))
    GrColo_SP$LeafSize <- 20 * log10(GrColo_SP$PopInit) / max(log10(GrColo_SP$PopInit), na.rm=T)

    # Add colonies to leaflet map
    Leaf <- Leaf %>%
      #addPolygons(data=GrColo_SP_MCP, fillOpacity=0.7, fillColor="#f1b6da", stroke=F, group=names(mod_out)[SP]) %>%
      #addCircleMarkers(lng=Colonies_SP$Lon, lat=Colonies_SP$Lat, color = "#8e0152", fillOpacity=0.7, label = paste0("Colonie : ", Colonies_SP$Popup), group=names(mod_out)[SP], radius=6) %>%
      addCircleMarkers(data=GrColo_SP, radius=GrColo_SP$LeafSize, color=ifelse(GrColo_SP$Parcs_foraging=="aucun", "#7fbc41", "#8e0152"), fill=T, fillColor = ColPal(GrColo_SP$RelImpact), label=paste0("Cliquez pour les résultats du groupe de colonies n°", GrColo_SP$group), fillOpacity=1, popup = GrColo_SP$Popup, group=names(mod_out)[SP])
  }

  ### Add control layer per species
  Leaf <- Leaf %>%
    hideGroup(., group=names(mod_out)) %>%
    addPolygons(data=st_buffer(st_transform(parcs, st_crs(4326)), 1000), fill="black", col="black", fillOpacity=1, label=paste0("Parc éolien : ", parcs$NAME)) %>%
    addLayersControl(baseGroups = c("", names(mod_out)), position="bottomleft", options = layersControlOptions(collapsed = F))%>%
    addLegend(pal=ColPal, values=c(0,100), title="Impact relatif (%)") %>%
    addLegend(colors=HTML("black; width:20px; height:20px; border:3px solid black; border-radius:0%"), opacity=1, labels="Parc éolien") %>%
    addLegend(colors=c(HTML("white; border:5px solid #d01c8b; border-radius:50%"), HTML("white; border:5px solid #b8e186; border-radius:50%")), labels=c("Groupe impacté", "Groupe non impacté"))

  if(is.null(BufferForaging)==F){
    Leaf <- Leaf %>% addLegend(colors=HTML("#white; width:20px; height:20px; border:3px dashed #8c2d04; border-radius:0%"), labels="Zone d'alimentation")
  }

  return(Leaf)
}
