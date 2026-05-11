setwd(dir=ifelse(file.exists("C:/Users/Victor"),
                 "C:/Users/Victor/Documents/Projects/Birdynamic/", # Setwd Victor
                 "D:/Documents/JOBS/Independant/Birdynamic/")) # Setwd Baptiste

min_year = 2009      # Remove years before...
last_year = 2021     # Remove years after...
species <- c("Cormoran huppé","Goéland marin")

#effec00 <- read.csv2("0.BV_Data_input/BD_Effectifs_1997-2021.csv")
load("app/0.Data/clean_data_all_species.RData")

#this is to get the data frame as it should be entered by the user
#(and in the same order as Thierry's code for comparison)
effec00 <- arrange(effec01,facade, secteur, desc(Lat + Lon), Colonie, sous_unite)
effec00 <- data.frame("Colonie"=paste(effec00$secteur,effec00$Colonie,effec00$sous_unite,sep=" - "),
                 "Lat"=effec00$Lat,"Lon"=effec00$Lon,"espece"=effec00$espece,"annee"=effec00$an,
                 "EFF_Min"=effec00$EFF_Min,"EFF_Max"=effec00$EFF_Max,"EFF_Moy"=effec00$EFF_Moy)
#rm(colonies00_L93,effec01)

colonies00 <- dplyr::select(.data = effec00,Colonie,Lat, Lon,) %>% unique

#this is to make sure that colonie code is the same as in original code
colo <- paste(colonies00_L93$secteur,colonies00_L93$Colonie,colonies00_L93$sous_unite,sep=" - ")
colonies00 <- data.frame("Colonie"=colo,"Lat"=colonies00$Lat[match(colo,colonies00$Colonie)],"Lon"=colonies00$Lon[match(colo,colonies00$Colonie)])

## Create column "code_colonie"
colonies00$code_colonie <- paste0("colo_", sprintf("%03d", 1:nrow(colonies00)))
colonies00 <- colonies00[which(!is.na(colonies00$Lat)),]

colonies00 <- st_as_sf(colonies00, coords = c("Lon", "Lat"), crs = 4326, agr = "constant",remove=F) %>% st_transform(2154)

## Add code_colonie to "effec00"
effec00$code_colonie <- colonies00$code_colonie[match(effec00$Colonie,colonies00$Colonie)]

##---------------------------
##here we enter species level
#----------------------------

sp <- species[1]
effecSp <- effec00[which(Bdy_clean_names(effec00$espece) == Bdy_clean_names(sp) &
                           effec00$annee >= min_year &
                           effec00$anne <= last_year), ]

coloniesSp <- colonies00[which(colonies00$code_colonie %in% effecSp$code_colonie),]

## Add column "effectifs"
agrEff <- round(tapply(effecSp$EFF_Moy,effecSp$code_colonie,mean),1)
coloniesSp$Eff <- agrEff[match(coloniesSp$code_colonie,names(agrEff))]

## Remove colonies with effectif max = 0
coloniesSp <- coloniesSp[which(coloniesSp$Eff > 0), ]

## Sort colonies : north --> south
coloniesSp <- arrange(coloniesSp, desc(Lat))
# colonies <- arrange(colonies, desc(scale(Lat*1 + Lon*0)))

effecSp <- effecSp[which(effecSp$code_colonie %in% coloniesSp$code_colonie), ]

##make clusters
coloniesSp$group <- Bdy_Clustering(coord=coloniesSp[c("Lon","Lat")])
n_group <- nlevels(coloniesSp$group)

## Get count data
## Données de comptages agrégés : colonie x an

counts00 <- as.data.frame(tapply(effecSp$EFF_Moy,list(effecSp$code_colonie,effecSp$annee),mean))

## Define missing year
u_yr <- (effecSp$annee %>% unique %>% sort)
all_yr <- (min_year:last_year)
miss_yr <- all_yr[!(all_yr %in% u_yr)]

rowNames <- rownames(counts00)
## Add missing years (all NA's)
if(length(miss_yr) > 0){
  for(k in 1:length(miss_yr)){
    counts00 <- counts00 %>% add_column(new = NA, .after = paste(miss_yr[k]-1))
    colnames(counts00)[colnames(counts00) == "new"] <- paste(miss_yr[k])
  } #k
} # end if
rownames(counts00) <- rowNames

## Total count et Moyenne (sur les années) par colonie

## Function "get_last survey"
get_last <- function(x) x[max(which(!is.na(x)))]

counts01 <- as.matrix(counts00)
for(j in 1:nrow(counts00)){
  counts00$tot[j] <- sum(counts01[j,], na.rm = TRUE)
  counts00$avg[j] <- mean(counts01[j,], na.rm = TRUE)
  counts00$last[j] <- get_last(counts01[j,])
} # j
rm(counts01)

## Rename years
names(counts00)[names(counts00) %in% min_year:last_year] <- paste0("X",min_year:last_year)
counts00 <- cbind("code_colonie"=row.names(counts00),
                  "group"=coloniesSp$group[match(row.names(counts00),coloniesSp$code_colonie)],
                  counts00)

## Aggregate colonie_counts by groups
sum2 <- function(x)ifelse(all(is.na(x)), NA, sum(x, na.rm = T))

group_counts <- apply(counts00[grep("X",colnames(counts00))],2,
                     FUN = function(x)tapply(x,counts00$group,sum2))


## Proportion of colonies (from a group) surveyed each year
# as a measure of proportional "count effort"
ncs_yr <- ppc_yr <- ppa_yr <- group_counts
ncs_yr[] <- ppc_yr[] <- ppa_yr[] <- 0

## Loop over each group
for(i in 1:n_group){
  gp_i <- levels(coloniesSp$group)[i]

  tbl00 <- as.data.frame(rbind(
    xtabs(as.formula(paste0("X", min_year, "~ code_colonie")), data = counts00, addNA = TRUE, subset = group == gp_i)
  ))

  ## Loop over years
  for(yr in (min_year+1):last_year){
    tbl00 <- rbind(tbl00,
                   xtabs(as.formula(paste0("X", yr, "~ code_colonie")), data = counts00, addNA = TRUE, subset = group == gp_i)
    )
  } # yr
  rm(yr)

  ## Proportion of abundance of each colony
  tot_avg <- xtabs(avg ~ code_colonie, data = counts00, addNA = TRUE, subset = group == gp_i)
  ppa_col <- tot_avg/sum(tot_avg, na.rm = TRUE)

  ## Year x Colony : surveyed or not ?
  ppa_yr[i,] <- (ifelse(is.na(tbl00), 0, 1) %*% t(t(ppa_col)))

  ## Number and proportion of colonies surveyed each year
  ncs_yr[i,] <- (ifelse(is.na(tbl00), 0, 1) %>% rowSums)
  ppc_yr[i,] <- (ncs_yr[i,]/ncol(tbl00))

  rm(ppa_col) ; rm(tbl00) ; rm(tot_avg)
} # i

rownames(ncs_yr) <- rownames(ppc_yr) <- rownames(ppa_yr) <- paste("group", levels(coloniesSp$group))
head(ppa_yr)
#ncs_yr ; ppc_yr

## Corrected colonie_counts
corrected_counts <- round(group_counts / ppa_yr)

