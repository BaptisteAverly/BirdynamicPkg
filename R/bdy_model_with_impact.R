#' Model with impact
#'
#' Computes the population model with impact from wind farms for the birds species of interest
#'
#' @param n_group numeric, number of groups of colonies for the bird species of interest
#' @param ny_data numeric, number of years for which count data exist for the bird species of interest
#' @param posterior dataframe, estimated output from [bdy_model_no_impact()]  \cr
#'                  columns contain estimated group size for the first year of projection (one column per group of colonies) and estimated annual growth rate (one column per group of colonies).
#' @param mortality Matrix (rows = iterations, columns = groups of colonies) giving the distribution of mortality accross groups of colonies, as outputted by [bdy_process_mortality()]
#' @param ni number of iterations
#' @param ny_proj number of years from the last annual count for which to compute population projections
#' @param progress R shiny Progress object used to show the calculation progress to the user. Only for use within a shiny application.
#'
#' @returns List of two arrays of size (n_group)*(ny_data)*(ni) containing projected population size:
#'          \itemize{
#'          \item sc0:  projections with no impact from wind farms (null model).
#'          \item sc1: projections with impacts from wind farms.
#'          }
#' @export
#'

bdy_model_with_impact <- function(n_group,ny_data,posterior,mortality,ni=1000,ny_proj=30,progress=NULL){

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

  ## Build object to store subpop sizes (under both scenarios)
  n_sc0 <- array(NA, dim = c(n_group, ny_proj, ni))
  morta_rate <- array(NA, dim = c(n_group, ni))

  ## Fill initial subpop sizes
  spl <- sample(nrow(post_n0), size = ni, replace = TRUE)
  #n_sc0[,1,] <- round(t(post_n0[spl,]))
  n_sc0[,1,] <- apply( round(t(post_n0[spl,])) , c(1,2), max, 1)
  n_sc1 <- n_sc0

  ## Draw iterations for subpop growth rates
  lam <- t(post_lam[spl,])

  ##### Run #####-

  ### Loop over iterations
  time_sim <- system.time(

    for(k in 1:ni){

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
  return(list(sc0=n_sc0, sc1=n_sc1))
}
