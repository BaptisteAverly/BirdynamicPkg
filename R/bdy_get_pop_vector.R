#' Get a population size vector
#'
#' Builds a vector of population size, for each age class, using a Leslie matrix to get the Stable Age Distribution (SAD)
#'
#' @param nb_pair a single number. Number of Pairs of reproductive individuals.
#' @param s a numeric vector of survival probabilities for each age class
#' @param f a numeric vector of fecundity values for each age class
#' @param pr a numeric vector of proportions of reproductive individuals for each age class
#' @param type character, either "pre" or "post" breeding
#'
#' @returns numeric vector giving for each age class the estimated total number of individuals
#'
#' @export

bdy_get_pop_vector <- function(nb_pair, s, f, pr, type="post"){

  N00 <- nb_pair*2

  # Get the LESLIE matrix, SAD and mature age classes
  A <- build_Leslie(s, f, type = "post")
  build_Leslie <- function(s, f, type = "pre"){ elements_Leslie(s=s, f=f, type=type)$A }


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
