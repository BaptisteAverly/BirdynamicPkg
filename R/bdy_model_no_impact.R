#' Model with no impact
#'
#' Computes the null population model (without impact from wind farms) for the birds species of interest
#'
#' @param count_data matrix giving annual bird counts for each group of colonies (rows) and each year of interest (columns)
#' @param PI matrix giving the proportion of colonies monitored for each group of colonies (rows) and each year of interest (columns)
#' @param survival numeric vector giving the survival rates for the different age classes of the species of interest
#' @param fecundity numeric vector giving the fecundity rates for the different age classes of the species of interest
#' @param propRepro numeric vector giving the proportion of reproductive individual for the different age classes of the species of interest
#' @param modelFile character, relative path for the text file containing the population model as jags code
#' @param nimble boolean, whether to use nimble of jags to compute the bayesian model
#' @param lightResults boolean, whether the output should only be the estimated output useful for the model with impact, or the full posterior distribution (can be quite heavy)
#' @param ny_proj number of years from the last annual count for which to compute population projections
#' @param na Number of iterations to run in the JAGS adaptive phase
#' @param nb Number of iterations at the beginning of the chain to discard (i.e., the burn-in). Does not include the adaptive phase iterations.
#' @param ni Total number of iterations per chain (including burn-in).
#' @param nc Number of Markov chains to run.
#' @param nt Thinning rate. Must be a positive integer.
#'
#' @returns Data frame with nrow = ni. \cr
#'          If lightResults=T, columns contain estimated group size for the first year of projection (one column per group of colonies) and estimated annual growth rate (one column per group of colonies).\cr
#'          If lightRestuls=F, columns contain the full posterior distribution: different estimated output for each group of colonies, and each year of interest (count years + projection years). \cr
#'          Refer to the documentation of packages 'nimble' or 'jagsUI' for details.
#' @export
#'

bdy_model_no_impact <- function(count_data,PI,survival,fecundity,propRepro,modelFile,nimble=T,lightResults=T,
                                ny_proj=30,na=5000,nb=1000,ni=20000 + nb,nc=3,nt=5){
  # B - from our understanding, not sure that ny_proj really as to be set to 30, since we only output the first year of projection

  ny_data <- ncol(count_data)
  ny_full <- ny_data + ny_proj
  count_data[, (ny_data+(1:ny_proj))] <- NA
  colnames(count_data)[ny_data+(1:ny_proj)] <- paste0("X", as.numeric(gsub("X", "", names(count_data)[ny_data])) + (1:ny_proj))
  PI[, (ny_data+1:ny_proj)] <- 1

  ### VITAL RATES & SAD #####
  # Get SAD factor "g"
  # from vital rates for that species
  g <- sum(bdy_get_pop_vector(nb_pair = 1000, s = survival, f = fecundity, pr = propRepro)[-1])/1000

  ### BAYESIAN ANALYSIS #####
  # Define max.N0 for the model
  max.N0 <- round(sum(count_data[,1]/(PI[,1]), na.rm = TRUE)*1.5)


  jags.data <- list(y = as.matrix(count_data), I = nrow(count_data), T = ny_full, ny_data = ny_data,
                    ny_proj = ny_proj, max.N0 = max.N0, PI=PI, g=g)

  # Parameters monitored
  parameters <- c("mu.lam_0", "b0", "b1", "N", "n", "n_TOT", "N_TOT", "gamma", "sig.y",
                  "growth_i_data", "growth_i_proj", "growth_i_full",
                  "growth_N_data", "growth_N_proj","growth_N_full")

  if(nimble){

    m01 <- bdy_modelCode()

    # Initial values
    y<-as.matrix(count_data)
    inits <- function(){
      list(
        N = round(sum(y[,1], na.rm = TRUE)*(1+runif(1,0,0.5))),
        n = round(y)
      )
    }

    # Call NIMBLE
    outNimble <-
      nimbleMCMC(
        code = m01,
        constants = jags.data,
        inits = inits(),
        monitors = parameters,
        thin = nt,
        niter = ni,
        nburnin = nb,
        nchains = nc,
        WAIC = FALSE
      )  # close nimbleMCMC

    posterior <- rbind(outNimble$chain1,outNimble$chain2,outNimble$chain3)
    nTotCol <- which(colnames(posterior) %in% paste0("n_TOT[", 1:nlevels(colonies$group), ", ", ncol(group_counts)+1, "]"))
    growthCol <- which(colnames(posterior) %in% paste0("growth_i_proj[", 1:nlevels(colonies$group),"]"))

    no_impact_output <- posterior[,c(nTotCol,growthCol)]
    colnames(no_impact_output) <- str_replace(colnames(no_impact_output)," ","")

  }else{
    system.time(
      outJags <- jags(data = jags.data, inits = NULL,
                      parameters.to.save = parameters,
                      model.file = modelFile,
                      n.chains = nc, n.iter = ni,
                      n.burnin = nb, n.adapt = na,
                      n.thin = nt, DIC = FALSE,
                      parallel = FALSE
      )
    )
    posterior <- rbind(outJags$samples[[1]],outJags$samples[[2]],outJags$samples[[3]])
    nTotCol <- which(colnames(posterior) %in% paste0("n_TOT[", 1:nlevels(colonies$group), ",", ncol(group_counts)+1, "]"))
    growthCol <- which(colnames(posterior) %in% paste0("growth_i_proj[", 1:nlevels(colonies$group),"]"))

    no_impact_output <- posterior[,c(nTotCol,growthCol)]
  }
  if(lightResults){
    return(no_impact_output)
  }else{
    return(posterior)
  }
}
