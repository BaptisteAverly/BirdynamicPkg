#' Check windfarms sf object
#'
#' Check the sf object LoadedShape to verify that each windfarm has a different name, corresponding with names from the mortality file.
#'
#' @param LoadedShape sf object (point or polygon) of windfarms
#' @param Names_mortality names of windfarms from the mortality table (a perfect match is required)
#'
#' @returns an HTML text message, in green if the sf object is valid, in red with details otherwise.
#' @export
#'

bdy_check_shapefile <- function(LoadedShape, Names_mortality){

  Warning <- NULL

  ### Check all parcs have names
  if((TRUE %in% is.na(LoadedShape$NAME)) | ("" %in% LoadedShape$NAME)){Warning <- "Certains parcs n'ont pas de nom attribué.<br>"}

  ### Check names from Parc are the same than names in mortality
  Names_shpONLY <- unique(LoadedShape$NAME[! LoadedShape$NAME %in% Names_mortality]) %>% subset(., . !="" & is.na(.)==F)
  Names_mortaONLY <- unique(Names_mortality[! Names_mortality %in% LoadedShape$NAME]) %>% subset(., . !="" & is.na(.)==F)
  if(length(Names_shpONLY)>0){Warning <- paste0(Warning, "Certains noms de parcs qui apparaissent sur la carte ne sont pas présents dans le fichier de mortalité : ", paste0(Names_shpONLY, collapse=", "), ".<br>")}
  if(length(Names_mortaONLY)>0){Warning <- paste0(Warning, "Certains noms de parcs sont présents dans le fichier de mortalité mais pas sur la carte : ", paste0(Names_mortaONLY, collapse=", "), ".<br>")}

  ### Check no duplicated names
  Dupli_names <- LoadedShape$NAME %>% table(.) %>% subset(., .>1) %>% names(.)
  if(any(Dupli_names=="")){
    Dupli_names <- Dupli_names[-which(Dupli_names==(""))]
  }
  if(length(Dupli_names)>0){Warning <- paste0(Warning, "Certains parcs sur la carte interactive ont le meme nom : ", paste0(Dupli_names, collapse=", "))}

  ### HTML warning
  if(is.null(Warning)){
    Warning <- "<p style='color:green'>Les données spatiales semblent conformes.</p>"
  } else {
    Warning <- paste0("<p style='color:red'>", Warning, "</p>")
  }

  ### Check parcs are within France limits (don't show names warning if it's not the case); only tested in Merge ObserveEvent
  if("CoverLimit_France" %in% names(LoadedShape)){
    if(0 %in% LoadedShape$CoverLimit_Europe){Warning <- "<p style='color:red'>Certains parcs sont hors des limites autorisées pour cette application (en rouge sur la carte); ne gardez que des parcs en pleine mer et localises sur les côtes atlantiques francaises (ou des pays adjacents).</p>"}
    if(0 %in% LoadedShape$CoverLimit_France & grepl("color:green", Warning)){Warning <- "<p style='color:purple'>Les données spatiales semblent conformes; notez toutefois que Bird Dynamic ne calcule que l'impact sur les colonies d'oiseaux marins en France alors que certains parcs que vous avez renseignés sont localisés hors de France.</p>"}
  }

  ### Check we have a parc
  Warn_NoParc <- "<p style='color:red'>Les coordonnées d'un parc éolien doivent être fournies pour lancer l'analyse</p>"
  if(is.null(nrow(LoadedShape))){Warning <- Warn_NoParc} else {if(nrow(LoadedShape)==0){Warning <- Warn_NoParc}}

  return(Warning)
}
