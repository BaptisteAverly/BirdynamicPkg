#' Code null population model
#'
#' Nimble code for the null population model (without impact from wind farms). Meant to be called by [bdy_model_no_impact()]
#'
#' @returns The population model as nimble code
#' @export
#'

bdy_model_code <- function(){
  return(
    nimbleCode( {

      ## National Scale
      # Prior on lambda (population growth rate)
      mu.lam_0 ~ dunif(0, 10)

      # Prior on N0
      N[1] ~ dunif(0, max.N0)

      # Process model over time (population growth)
      for (t in 1:(T-1)){
        lambda[t] <- mu.lam_0
        N[t+1] <- N[t] * lambda[t]
      } # t


      ## Colony-groups Scale
      # Priors
      for(i in 1:I){
        b0[i] ~ dunif(0,100)
        b1[i] ~ dunif(0,100)
      } # i

      for(t in 1:T){
        for(i in 1:I){
          w[i,t] <- b0[i] + b1[i]*(t-1)
          gamma[i,t] <- w[i,t]/sum(w[1:I,t])
        } # i
      } # t

      # Process model : repartition of N among colony-groups
      for(t in 1:T){
        # NN[t] <- round(N[t])
        nb[1:I,t] <- gamma[1:I,t]*N[t]
        n[1:I,t] <- nb[1:I,t]
        n_TOT[1:I,t] <- g*n[1:I,t]

        N_TOT[t] <- sum(n_TOT[1:I,t])

      } # t

      ## Observation process
      # Priors
      sig.y ~ dunif(0.1, 500)
      tau.y <- pow(sig.y, -2)

      # Obs model
      for(i in 1:I){
        for(t in 1:T){
          y[i,t] ~ dnorm(n[i,t]*PI[i,t], sd= sig.y)
        } # t
      } # i

      ## Derived parameters : relative impact and growth_i (colony scale)
      for(i in 1:I){
        growth_i_data[i] <- pow((n_TOT[i,ny_data]/n_TOT[i,1]), (1/(ny_data-1)))
        growth_i_proj[i] <- pow((n_TOT[i,T]/n_TOT[i,ny_data]), (1/(ny_proj)))
        growth_i_full[i] <- pow((n_TOT[i,T]/n_TOT[i,1]), (1/(T-1)))

      } # i

      growth_N_data <- pow((N[ny_data]/N[1]), (1/(ny_data-1)))
      growth_N_proj <- pow((N[T]/N[ny_data]), (1/ny_proj))
      growth_N_full <- pow((N[T]/N[1]), (1/(T-1)))

    } ) # End model description)
  )
}
