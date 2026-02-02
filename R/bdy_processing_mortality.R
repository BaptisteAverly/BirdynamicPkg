#' Title
#'
#' @param collision
#' @param season
#' @param n_iteration
#' @param parcNames
#' @param RW_group
#'
#' @returns
#' @export
#'
#' @examples

bdy_processing_mortality <- function(collision,season,n_iteration=1000,parcNames,RW_group){

  ### Make the table to store distribution of collision risk for that species (Iter x Parc)
  #morta_distri <- matrix(NA, nrow = max(collision$iteration), ncol = n_parc, dimnames = list(NULL, parcs_L93$NAME))
  morta_distri <- matrix(NA, nrow = n_iteration, ncol = length(parcNames),dimnames=list(NULL, parcNames))
  collision$coefficient <- as.numeric(collision$coefficient)
  ######################################################################-
  ### 2. PARC : Loop over windfarms (Parc)
  for(kk in 1:length(parcNames)){

    sel <- which(collision$parc == parcNames[kk])

    if(length(sel) == 0){} else{
      cr <- collision[sel,]

      ### 3. OPTION : select the option to keep - B no need for options anymore
      #tmp <- aggregate(coefficient ~ Option, data = cr, mean)
      #sel_option <- tmp$Option[which.max(tmp$coefficient)]
      #cr <- cr[cr$Option == sel_option,]

      ## Months when only local individuals are present
      # cr_sel <- cr[cr$Month %in% months_breeding,]
      sel_local <- which(season %in% c("B", "R"))
      cr_sel <- cr[cr$mois %in% sel_local,]
      cr_local <- aggregate(coefficient ~ iteration, data = cr_sel, sum) # sum over the Breeding season
      cr_local <- (cr_local$coefficient * rep(1, n_iteration)) #B - /!\ This gives a warning /!\
      rm(cr_sel) ; rm(sel_local)

      ## Months when both local and migrants are present (X% affects local populations)
      # cr_sel <- cr[cr$Month %in% months_winter,]
      sel_mix <- which(season %in% c("M", "T"))
      if(length(sel_mix) > 0){
        cr_sel <- cr[cr$mois %in% sel_mix,]
        cr_mix <- aggregate(coefficient ~ iteration, data = cr_sel, sum) # sum over the Breeding season
        # Incertitude sur la proportion des collisions concernant les populations locales
        min_PROP_LOCAL <- min((mean(cr_local)/mean(cr_mix$coefficient)), 1, na.rm = TRUE) # si il y a moins de collisions en période hors-repro qu'en repro, on considère que hors période repro les collisions affectent à 100% les individus locaux
        max_PROP_LOCAL <- 1
        cr_mix <- (cr_mix$coefficient * runif(n_iteration, min_PROP_LOCAL, max_PROP_LOCAL)) #B - /!\ This gives a warning /!\
      } else {
        cr_mix <- 0
      } # close if


      ## Draw many (n_iteration) values (= shuffled distribution)
      cr_mix <- sample(cr_mix, size = n_iteration, replace = TRUE)
      cr_local <- sample(cr_local, size = n_iteration, replace = TRUE)

      ### Fill the table (Iter x Parc)
      # morta_distri[ , parcs_L93$NAME[kk]] <- cr_annual$coefficient
      morta_distri[ , parcNames[kk]] <- cr_local + cr_mix

    } # close if
    rm(sel)
  } # kk
  dim(morta_distri)
  head(morta_distri)


  ## Replace remaining NA by ZERO (these are correct !)
  morta_distri[is.na(morta_distri)] <- 0
  dim(morta_distri)
  head(morta_distri)

  ## Draw many (n_iteration) values for each parc (= shuffled distribution)
  morta_iteration_parc <- apply(morta_distri, 2, sample, size = n_iteration, replace = TRUE)
  #head(morta_iteration_parc)
  #dim(morta_iteration_parc)

  ## Apportioning mortalities per cluster or group - dim : Iterations X Groups
  morta_iteration_gp <- morta_iteration_parc %*% t(RW_group)
  return(morta_iteration_gp)
}
