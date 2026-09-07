#' Get a population size vector
#'
#' Builds a vector of population size, for each age class, using a Leslie matrix to get the Stable Age Distribution (SAD). Used by [bdy_model_no_impact]
#'
#' @param nb_pair a single number. Number of Pairs of reproductive individuals.
#' @param s numeric vector giving the survival rates for the different age classes of the species of interest, for example from column 'survival' of [bdydata_vital_rates]
#' @param f numeric vector giving the fecundity rates for the different age classes of the species of interest, for example from column 'fecundity' of [bdydata_vital_rates]
#' @param pr numeric vector giving the proportion of reproductive individual for the different age classes of the species of interest, for example from column 'propRepro' of [bdydata_vital_rates]
#' @param type character, either "pre" or "post" breeding
#'
#' @returns numeric vector giving for each age class the estimated total number of individuals
#'
#' @export

bdy_get_pop_vector <- function(nb_pair, s, f, pr, type="post"){

  N00 <- nb_pair*2

  nac <- length(s)
  vr_list <- as.list(c(s,f))
  names(vr_list) <- c(paste0("s", (1:nac)-1), paste0("f", (1:nac)-1))
  vital_rates <- unlist(vr_list)

  # If Pre-Breeding
  if(type == "pre"){
    A <- diag(0, nrow = nac-1)

    A[1,] <- paste0(  paste0("f",1:(nac-1)), "*", "s0"  )

    if((nac-1) < 3){
      A[-1,] <- paste0("s", 1:(nac-1))
    }else{
      diag(A[-1,]) <- paste0("s", 1:(nac-2))
      A[nac-1,nac-1] <- paste0("s", nac-1)
    }

    # If Post-Breeding
  }else{
    A <- diag(0, nrow = nac)

    A[1,] <-  paste0(paste0("s", (1:nac)-1),"*",paste0("f", c((2:nac)-1,nac-1)))
    if(nac < 3){
      A[-1,] <- paste0("s", (1:nac)-1)
    }else{
      diag(A[-1,]) <- paste0("s", (1:(nac-1))-1)
      A[nac,nac] <- paste0("s", nac-1)
    }

  } # enf if

  symbolic <- noquote(A)
  elements <- parse(text=t(A))

  A <- matrix( sapply(elements, eval, vr_list), nrow=nrow(A), byrow=TRUE)

  SAD <- stable.stage(A)
  mature <- which(f != 0) # identify mature age classes

  Ntot <- (N00/sum(SAD*pr))
  N0 <- round(Ntot*SAD)

  return(N0)
}
