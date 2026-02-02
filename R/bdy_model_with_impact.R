
#' Title
#'
#' @param n_group
#' @param ny_data
#' @param posterior
#' @param mortality
#' @param ni
#' @param ny_proj
#' @param nimble
#' @param progress
#'
#' @returns
#' @export
#'
#' @examples

#ni = number of iteration
bdy_model_with_impact <- function(n_group,ny_data,posterior,mortality,ni=1000,ny_proj=30,nimble=F,progress=NULL){

  ###############################################################################
  #### IMPACT SIMULATIONS                                                   #####
  #### WITH PARAMETER UNCERTAINTY and DEMOGRAPHIC STOCHASTICITY             #####
  ##############################################################################

  ##Select posterior of target parameters
  # Initial subpop size : prospective period
  sel <- which(colnames(posterior) %in% paste0("n_TOT[", 1:n_group, ",", ny_data+1, "]"))
  post_n0 <- posterior[,sel]

  # Annual subpop growth rate  : prospective period
  sel <- which(colnames(posterior) %in% paste0("growth_i_proj[", 1:n_group,"]"))
  post_lam <- posterior[,sel]

  ## Draw values of initial subpop sizes
  #n0_distri <- apply(post_n0, 2, sample, size = ni, replace = TRUE)

  #colMeans(morta_distri)
  #head(mortality)
  #(morta_iter_parc) %>% colMeans


  ## Buil object to store subpop sizes (under both scenarios)
  n_sc0 <- array(NA, dim = c(n_group, ny_proj, ni))
  morta_rate <- array(NA, dim = c(n_group, ni))

  ## Fill initial subpop sizes
  spl <- sample(nrow(post_n0), size = ni, replace = TRUE)
  #n_sc0[,1,] <- round(t(post_n0[spl,]))
  n_sc0[,1,] <- apply( round(t(post_n0[spl,])) , c(1,2), max, 1)
  n_sc1 <- n_sc0

  ## Draw iterations for subpop growth rates
  #### ATTENTION : DO NOT resample 'spl" here ; needs to keep correlation structure ; otherwise UNCERTAINTY is underestimated
  lam <- t(post_lam[spl,])
  # lam[,1]

  ##### Run #####-
  # mu_m <- apply(mortality, 2, median) # without uncertainty on mortality

  ### Loop over iterations
  time_sim <- system.time(

    for(k in 1:ni){

      #if((k %% (ni/100)) == 0) print(paste("simulation :", (k/ni)*100, "%"))

      ### Mortalities (collisions)
      # mu_m <- mortality[k,] # NUMBER of mortalities (independent from Current Pop Size) -- with uncertainty on mortality
      mu_m <- mortality[k,]/n_sc0[,1,k] # mortality RATE -- with uncertainty on mortality
      morta_rate[,k] <- mu_m

      for(t in 2:ny_proj){

        # sc0 : no collision fatalities
        n_sc0[,t,k] <- rpois(n=n_group, n_sc0[,t-1,k]*lam[,k])

        ## save demographic stochasticity
        fac_ds <- 1 + ((n_sc0[,t,k] - (n_sc0[,t-1,k]*lam[,k]))/(n_sc0[,t-1,k]*lam[,k]))
        fac_ds[is.nan(fac_ds)] <- 1 ## for cases where n = 0

        # sc1: with collision fatalities : MORTALITY RATE
        n_tmp <- round((n_sc1[,t-1,k]*lam[,k]) * fac_ds)

        #real_m <- rpois(n=n_group, mu_m) # Apply mortality as a NUMBER independent from Current Pop Size
        real_m <- rpois(n=n_group, (mu_m*n_tmp)) # Apply mortality as a RATE

        for(i in 1:n_group) real_m[i] <- min(real_m[i], n_tmp[i])
        n_sc1[,t,k] <- n_tmp - real_m

      } # t

      if(((k %% (ni/100)) == 0) && !is.null(progress)){
        #if(!is.null(progress)){
        progress$inc(0.01, detail = paste0("Iteration n° ", k,"  / ",ni))
        if(k==ni)Sys.sleep(1)
      }

    } # k (simulation iterations)
  ) # time_sim
  print(time_sim)
  return(list(sc0=n_sc0, sc1=n_sc1))
}
