#' Title
#'
#' @param sp
#' @param effec00
#' @param colonies00
#' @param first_year
#' @param last_year
#'
#' @returns
#' @export
#'
#' @examples
#'
#'
bdy_process_count_data <- function(sp,countData,colonies,first_year,last_year){

  #---check species name-----


  #--------------------------

  #selecting only the count data for the species of interest
  effecSp <- countData[which(countData$species_latin == sp &
                             countData$year >= first_year &
                             countData$year <= last_year), ]

  #adding colony code to count data
  effecSp$colony_code <- colonies$colony_code[match(effecSp$colony,colonies$colony)]

  #this is for consistency with Thierry's code, but unlikely to be used by users
  if(is.null(effecSp$regroup)){
    effecSp$regroup <- ""
  }
  effecSp$rgp <- FALSE
  effecSp$rgp[(effecSp$regroup != "")] <- TRUE

  #selecting only the colonies for the species of interest
  coloniesSp <- colonies[which(colonies$colony_code %in% effecSp$colony_code),]

  ## Add column "effectifs"
  agrEff <- round(tapply(effecSp$count_mean,effecSp$colony_code,mean),1)
  coloniesSp$mean <- agrEff[match(coloniesSp$colony_code,names(agrEff))]

  ## Remove colonies with effectif max = 0
  coloniesSp <- coloniesSp[which(coloniesSp$mean > 0), ]

  ## Sort colonies : north --> south
  coloniesSp <- arrange(coloniesSp, desc(lat))
  # colonies <- arrange(colonies, desc(scale(Lat*1 + Lon*0)))

  effecSp <- effecSp[which(effecSp$colony_code %in% coloniesSp$colony_code), ]

  ##make clusters
  coloniesSp$group <- bdy_clustering(coord=coloniesSp[c("lon","lat")])
  n_group <- nlevels(coloniesSp$group)

  ## Get count data
  ## Données de comptages agrégés : colonie x an

  counts00 <- as.data.frame(tapply(effecSp$count_mean,list(effecSp$colony_code,effecSp$year),mean))
  rgp00 <- as.data.frame(tapply(effecSp$rgp,list(effecSp$colony_code,effecSp$year),function(x)as.numeric(isTRUE(x))))

  ## Define missing year
  u_yr <- (effecSp$year %>% unique %>% sort)
  all_yr <- (first_year:last_year)
  miss_yr <- all_yr[!(all_yr %in% u_yr)]

  rowNames <- rownames(counts00)
  ## Add missing years (all NA's)
  if(length(miss_yr) > 0){
    for(k in 1:length(miss_yr)){
      counts00 <- counts00 %>% add_column(new = NA, .after = paste(miss_yr[k]-1))
      colnames(counts00)[colnames(counts00) == "new"] <- paste(miss_yr[k])

      rgp00 <- rgp00 %>% add_column(new = NA, .after = paste(miss_yr[k]-1))
      colnames(rgp00)[colnames(rgp00) == "new"] <- paste(miss_yr[k])
    } #k
  } # end if
  rownames(counts00) <- rowNames
  rownames(rgp00) <- rowNames

  ## Total count et Moyenne (sur les années) par colonie

  ## Function "get_last survey"
  get_last <- function(x) x[max(which(!is.na(x)))]

  counts01 <- as.matrix(counts00)
  for(j in 1:nrow(counts00)){
    sel <- which(rgp00[j,] == 0) # use only data that are not part of a "regroup" thing
    if(length(sel)>0){
      counts00$tot[j] <- sum(counts01[j,sel])
      counts00$avg[j] <- mean(counts01[j,sel])
      counts00$last[j] <- get_last(counts01[j,sel])
    }else{
      counts00$tot[j] <- sum(counts01[j,], na.rm = TRUE)
      counts00$avg[j] <- mean(counts01[j,], na.rm = TRUE)
      counts00$last[j] <- get_last(counts01[j,])
    }

  } # j
  rm(counts01)

  ## Rename years
  names(counts00)[names(counts00) %in% first_year:last_year] <- paste0("X",first_year:last_year)
  counts00 <- cbind("colony_code"=row.names(counts00),
                    "group"=coloniesSp$group[match(row.names(counts00),coloniesSp$colony_code)],
                    counts00)

  ## Aggregate colonie_counts by groups
  sum2 <- function(x)ifelse(all(is.na(x)), NA, sum(x, na.rm = T))

  group_counts <- as.data.frame(apply(counts00[grep("X",colnames(counts00))],2,
                        FUN = function(x)tapply(x,counts00$group,sum2)))


  ## Proportion of colonies (from a group) surveyed each year
  # as a measure of proportional "count effort"
  ncs_yr <- ppc_yr <- ppa_yr <- group_counts
  ncs_yr[] <- ppc_yr[] <- ppa_yr[] <- 0

  ## Loop over each group
  for(i in 1:n_group){
    gp_i <- levels(coloniesSp$group)[i]

    tbl00 <- as.data.frame(rbind(
      xtabs(as.formula(paste0("X", first_year, "~ colony_code")), data = counts00, addNA = TRUE, subset = group == gp_i)
    ))

    ## Loop over years
    for(yr in (first_year+1):last_year){
      tbl00 <- rbind(tbl00,
                     xtabs(as.formula(paste0("X", yr, "~ colony_code")), data = counts00, addNA = TRUE, subset = group == gp_i)
      )
    } # yr

    ## Proportion of abundance of each colony
    tot_avg <- xtabs(avg ~ colony_code, data = counts00, addNA = TRUE, subset = group == gp_i)
    ppa_col <- tot_avg/sum(tot_avg, na.rm = TRUE)

    ## Year x Colony : surveyed or not ?
    ppa_yr[i,] <- (ifelse(is.na(tbl00), 0, 1) %*% t(t(ppa_col)))

    ## Number and proportion of colonies surveyed each year
    ncs_yr[i,] <- (ifelse(is.na(tbl00), 0, 1) %>% rowSums)
    ppc_yr[i,] <- (ncs_yr[i,]/ncol(tbl00))

  } # i

  rownames(ncs_yr) <- rownames(ppc_yr) <- rownames(ppa_yr) <- paste("group", levels(coloniesSp$group))

  return(list("species"=sp,"group_counts_sp"=group_counts,"ppa_yr_sp"=as.data.frame(ppa_yr),"colonies_sp"=coloniesSp))

}
