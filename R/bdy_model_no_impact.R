
#' Title
#'
#' @param count_data
#' @param PI
#' @param survival
#' @param fecundity
#' @param propRepro
#' @param modelFile
#' @param nimble
#' @param ny_proj
#' @param na
#' @param nb
#' @param ni
#' @param nc
#' @param nt
#'
#' @returns
#' @export
#'
#' @examples
#'

#count_data = comptages annuels
#PI = Proportion des colonies suivies par groupe, chaque année
#ny_proj = Ajouter des NA's pour les projections démo
bdy_model_no_impact <- function(count_data,PI,survival,fecundity,propRepro,modelFile,nimble=T,
                                ny_proj=30,na=5000,nb=1000,ni=20000 + nb,nc=3,nt=5){
  ny_data <- ncol(count_data)
  ny_full <- ny_data + ny_proj
  count_data[, (ny_data+(1:ny_proj))] <- NA
  colnames(count_data)[ny_data+(1:ny_proj)] <- paste0("X", as.numeric(gsub("X", "", names(count_data)[ny_data])) + (1:ny_proj))
  PI[, (ny_data+1:ny_proj)] <- 1

  ### VITAL RATES & SAD #####
  # Get SAD factor "g"
  # from vital rates for that species
  g <- sum(pop_vector(nb_pair = 1000, s = survival, f = fecundity, pr = propRepro)[-1])/1000

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

    m01 <- Bdy_modelCode()

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
  return(no_impact_output)
}
