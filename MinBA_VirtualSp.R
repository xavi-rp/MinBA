
# Script to create virtual species

library(virtualspecies)


#Creating Raster Stack with variables (WorldCLim)

varbles <- paste0("~/Documents/MinBA_models/", "wc5_other/biovars_eur")
rstrs <- list.files(varbles, pattern = c(".bil$"), full.names = T)
rstrs <- c(rstrs, list.files(varbles, pattern = c(".tif$"), full.names = T))

#vrbles <- stack()
for(rst in 1:length(rstrs)){
  temp <- raster::raster(rstrs[rst])
  if(rst == 1){
    vrbles <- raster::stack(temp)
  }else{
    vrbles <- raster::stack(vrbles, temp)
  }
}

vrbles
graphics.off()
vrbles@extent
plot(vrbles[[1]])


?generateSpFromFun
?generateRandomSp


virtSps <- list()  #list of virtual species complete information to be saved
virtSps_pos <- 0
virtSp_occurrences <- as.data.frame(matrix(nrow = 0, ncol = 0))    #data frame with coordinates of occurrences of all virtual species
repetitions <- 10
vrbles_subset <- c(12, 15, 16, 5:7)
nche_brea <- c("narrow", "wide")
nche_brea <- c("narrow")

#for(nbrth in c("narrow", "wide")){
#for(nbrth in c("narrow")){
for(nbrth in nche_brea){
  #bta <-c(0.1, 0.2, 0.3, 0.4, 0.5)
  bta <-c(0.1, 0.175)
  alp <- c(-0.01, -0.05)
  #for(vsp in 1:5) {
  for(bt in bta) {
    #prvl <- (vsp / 10)
    if(bt < 0.25){
      aph = alp[1]
    }else{
      aph = alp[2]
    }
    rp <- 0
    repeat{
      virtSps_pos <- virtSps_pos + 1
      rp <- rp + 1
      print(paste0(virtSps_pos, ": virtSp_", nbrth, "_beta", gsub("\\.", "", as.character(round(bt, 2))), "_alpha", gsub("\\.", "", as.character(round(aph, 2))), "_sp", rp))
      virtSps[[virtSps_pos]] <- generateRandomSp(#raster.stack = vrbles, 
                                                 raster.stack = vrbles[[vrbles_subset]], 
                                                 #approach = "automatic", # if > 6 variables, automatic = PCA
                                                 approach = "pca", # if > 6 variables, automatic = PCA
                                                 rescale = FALSE,  # If TRUE, the final probability of presence is rescaled between 0 and 1
                                                 convert.to.PA = TRUE, 
                                                 #convert.to.PA = FALSE, 
                                                 #relations = c("gaussian", "linear", "logistic", "quadratic"), 
                                                 #relations = c("quadratic"), 
                                                 relations = c("logistic"), 
                                                 rescale.each.response = TRUE, realistic.sp = TRUE,
                                                 #species.type = c("additive", "multiplicative"), 
                                                 species.type = "additive", 
                                                 #niche.breadth = c("any", "narrow", "wide"),
                                                 niche.breadth = nbrth,
                                                 #sample.points = FALSE, nb.points = 10000,
                                                 sample.points = TRUE, nb.points = 20000,
                                                 PA.method = "probability", 
                                                 alpha = aph, 
                                                 #adjust.alpha = TRUE,
                                                 #beta = "random", 
                                                 beta = bt, 
                                                 #species.prevalence = prvl,
                                                 plot = FALSE)
      
      assign("virtSpi", virtSps[[virtSps_pos]])
      #virtSpi
      #plot(virtSpi$probability.of.occurrence)
      #plot(virtSpi$suitab.raster)
      #plot(virtSpi$pa.raster)
      #head(values(virtSpi$pa.raster))
      #head(coordinates(virtSpi$pa.raster))
      virtSpi_occ <- as.data.frame(coordinates(virtSpi$pa.raster))
      virtSpi_occ$spec <- values(virtSpi$pa.raster)
      
      virtSpi_occ <- virtSpi_occ[virtSpi_occ$spec == TRUE, ]
      virtSpi_occ <- virtSpi_occ[!is.na(virtSpi_occ$spec), ]
      virtSpi_occ$spec <- paste0("virtSp_", nbrth, "_beta", gsub("\\.", "", as.character(round(bt, 2))), "_alpha", gsub("\\.", "", as.character(round(aph, 2))), "_sp", rp)
      names(virtSpi_occ) <- c("lon", "lat", "species")
      #head(virtSpi_occ)
      #nrow(virtSpi_occ)
      virtSp_occurrences <- rbind(virtSp_occurrences, virtSpi_occ)
      
      #assign(paste0("virtSp", vsp), virtSpi)
      if(rp == repetitions) break
    }
  }
}


#virtSps <- list(virtSp1,virtSp2,virtSp3,virtSp4,virtSp5)

write.csv(virtSp_occurrences, paste0(dir2save, "/virtSp_occurrences.csv"), row.names = FALSE)
save(virtSps, file = paste0(dir2save, "/virtSp_allRasters.RData"))


head(virtSp_occurrences)
unique(virtSp_occurrences$species)
table(virtSp_occurrences$species)
nrow(virtSp_occurrences)

#for (i in 1:length(unique(virtSp_occurrences$species))){
#  if(as.numeric(virtSps[[i]]$PA.conversion[3]) < -0.02){
#    print(i)
#    print(virtSps[[i]]$PA.conversion)
#    print(unique(virtSp_occurrences$species)[i])
#    
#  }
#}
#el limit es alpha = -0.02

#load(paste0(dir2save, "/virtSp_allRasters.RData"), verbose = T)
rast2rm <- c()
sp2rm <- c()
for (i in 1:length(unique(virtSp_occurrences$species))){
  if(as.numeric(virtSps[[i]]$PA.conversion[3]) < -0.03){
    print(i)
    print(virtSps[[i]]$PA.conversion)
    rast2rm <- c(rast2rm, i)
    sp2rm <- c(sp2rm, unique(virtSp_occurrences$species)[i])
  }
}
#el limit es alpha = -0.02
#20190904:  donen problemes la 3, 11 i 16
# Eliminaré les alpha > -0.03 (la 3 i 16)

virtSps_good <- virtSps[-c(rast2rm)]
virtSp_occurrences_good <- virtSp_occurrences[!virtSp_occurrences$species %in% sp2rm, ]

write.csv(virtSp_occurrences_good, paste0(dir2save, "/virtSp_occurrences_good.csv"), row.names = FALSE)
save(virtSps_good, file = paste0(dir2save, "/virtSp_allRasters_good.RData"))



#graphics.off()
pdf(paste0(dir2save, "/plot_", "virtSp_", nbrth,
           "_Beta", gsub("\\.", "", paste(range(bta), collapse = "-")), 
           "_Alpha", gsub("\\.", "", paste(gsub("\\-", "", alp), collapse = "-")),
           "_pa.raster.pdf"), width = 80, height = 55)
par(mfrow = c(3, 6))
for(i in c(1:(length(unique(virtSp_occurrences_good$species))))){
  plot(virtSps_good[[i]]$pa.raster, main = unique(virtSp_occurrences_good$species)[i], cex.main = 5, legend.width = 7)
}
dev.off()

pdf(paste0(dir2save, "/plot_", "virtSp_", nbrth,
           "_Beta", gsub("\\.", "", paste(range(bta), collapse = "-")), 
           "_Alpha", gsub("\\.", "", paste(gsub("\\-", "", alp), collapse = "-")),
           "_prob_occ.pdf"), width = 80, height = 55)
par(mfrow = c(3, 6))
for(i in c(1:(length(unique(virtSp_occurrences_good$species))))){
  plot(virtSps_good[[i]]$probability.of.occurrence, main = unique(virtSp_occurrences_good$species)[i], cex.main = 5, legend.width = 7)
}
dev.off()



virtSp_occurrences_good <- read.csv(paste0(dir2save, "/virtSp_occurrences_good.csv"), header = TRUE)
load(paste0(dir2save, "/virtSp_allRasters_good.RData"), verbose = TRUE)

