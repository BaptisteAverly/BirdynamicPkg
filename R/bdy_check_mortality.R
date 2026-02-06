#' Check mortality table
#'
#' Performs some checks and cleaning on the mortality table provided by a users to make sure it follows the format required by the app.
#'
#' @param Mortality A data frame giving mortality estimates for the species of interest and the wind farms of interest.Correctly formated ones should have the following columns:
#'                  \itemize{
#'                  \item 'espece_latin': character, name of the species in french
#'                  \item 'parc': character, names of the parc for which the mortality due to collisions is estimated
#'                  \item 'month': numeric, month of the year (1 to 12) for which the mortality is estimated
#'                  \item 'iteration': numeric, iteration index of the collision model.
#'                  \item 'coefficient': numeric, estimated mortality coefficient from the collision modèle, for a given combination of species, parc, month and iteration
#'                  }

#' @param speciesList character vector giving the names of all the species which are acceptable entries for the analysis
#'
#' @returns list of 2:
#'          \itemize{
#'          \item warning: string giving warnings, errors, or otherwise a green pass about the formatting of the table
#'          \item table: the mortality data frame with some dditional formatting
#'          }
#' @export
#'

bdy_check_mortality <- function(Mortality,speciesList){

  tableNames <- c("espece_latin","iteration","mois","parc","coefficient")

  Warning <- NULL
  Info <- NULL

  if(nrow(Mortality)==0){
    Warning <- paste0(Warning, "Le tableau chargé est vide.<br>")
  }else{

    ### check that there are no weird characters because of encoding
    error <- 0
    for(i in 1:ncol(Mortality)){
      if(is.character(Mortality[,i])){
        errorUTF <- !validUTF8(Mortality[,i])
        if(any(errorUTF)){
          error <- 1
          Mortality[which(errorUTF),i] <- "ERREUR_NOM"
        }
      }
    }
    colnames(Mortality)[!validUTF8(colnames(Mortality))] <- "ERREUR_NOM"

    if(error){
      Warning <- paste0(Warning,"Présence d'accents ou de caractères spéciaux, veuillez les enlever ou vérifier l'encodage.<br>")
    }

    ### check that all necessary columns are present
    colnames(Mortality) <- tolower(colnames(Mortality))
    missingColumns <- tableNames[which(!tableNames %in% colnames(Mortality))]
    if(length(missingColumns) > 0){
      Warning <- paste0(Warning, "Les colonnes suivantes sont manquantes: ",paste0(missingColumns,collapse = ", "), ".<br>")
    }

    if(is.null(Warning)){ # if all columns are present and no names errors

      loadedSpecies <- unique(Mortality$espece_latin)

      ###check that month is numeric
      Mortality$mois <- as.numeric(Mortality$mois)
      if(all(is.na(Mortality$mois))){
        Warning <- paste0(Warning,"La colonne mois n'est pas numérique", ".<br>")
      }else{
        ###check that each month is present at least once
        if(any(!(1:12) %in% Mortality$mois)){
          Warning <- paste0(Warning,"Certains mois ne sont pas présent dans le tableau", ".<br>")
        }

        ###check that months are between 1 and 12
        if(any(!Mortality$mois %in% 1:12)){
          Warning <- paste0(Warning,"Valeurs aberrantes dans la colonne mois: les indices doivent être compris entre 1 et 12.<br>")
        }
      }

      ### check that iteration is numeric
      #maybe add additional checks depending on what is expected for this
      Mortality$coefficient <- as.numeric(Mortality$coefficient)
      if(all(is.na(Mortality$coefficient))){
        Warning <- paste0(Warning,"La colonne coefficient n'est pas numérique.<br>")
      }

      ### check that species in the table correspond to the list
      wrongSpecies <- loadedSpecies[which(!loadedSpecies %in% speciesList)]
      if(length(wrongSpecies > 0)){
        Warning <- paste0(Warning, "Les noms d'espèces suivants ne correspondent pas aux oiseaux marins répertoriés: ",paste0(wrongSpecies,collapse = ", "), ".<br>")
      }
    }
  }

  ### if the blocking checks were passed, we do some formatting and non-blocking testing
  if(is.null(Warning)){
    Warning <- "<p style='color:green'>Les données de mortalité semblent conformes.</p>"

    Mortality$espece_latin[which(Mortality$espece_latin=="Gulosus aristotelis")] <- "Phalacrocorax aristotelis"
    Mortality$espece_latin[which(Mortality$espece_latin=="Larus melanocephalus")] <- "Ichthyaetus melanocephalus"
    Mortality$espece_latin[which(Mortality$espece_latin=="Larus ridibundus")] <- "Chroicocephalus ridibundus"

    Mortality$espece_latin <- as.factor(Mortality$espece_latin)
    Mortality$mois <- as.factor(Mortality$mois)
    Mortality$parc <- as.factor(Bdy_clean_names(Mortality$parc))
    Mortality$mois <- as.factor(Mortality$mois)
    if("espece" %in% names(Mortality)){Mortality$espece <- as.factor(Mortality$espece)}

    ### check for and remove duplicates
    dupl <- which(duplicated(Mortality[,c("espece_latin","iteration","mois","parc")]))
    if(length(dupl) > 0){
      Mortality <- Mortality[-dupl,]
      Warning <- paste0(Warning,"<p style='color:purple'>Les lignes ",paste(dupl,collapse=", ")," sont des doublons et ont été supprimées.</p>")
    }

    if(anyNA(Mortality$coefficient)){
      Warning <- paste0(Warning,"<p style='color:purple'>Certaines valeures dans la colonne coefficient ne sont pas numériques, création de NAs</p>")
    }
    return(list(warning=Warning,table=Mortality))

  }else {
    Warning <- paste0("<p style='color:red'>", Warning, "</p>")
    return(list(warning=Warning,table=Mortality))
  }
}
