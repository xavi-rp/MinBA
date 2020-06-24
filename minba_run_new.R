
###########################################################################
########                                                       ############
########             Minimum Background Area for SDMs          ############
########                          MinBA                        ############
########                                                       ############
###########################################################################
#
# minba_run_new.R
#
# Created: Winter 2018
#
# Updated: Spring 2020
#
# Created by: Xavier Rotllan-Puig (xavi.rotllan.puig@gmail.com)
#
# Description: The aim of this script is to run the calculations for the 
#              Case Studies of the paper
#
# ------------------------------------------

#### Call settings ####
if(Sys.info()[4] == "MacBook-MacBook-Pro-de-Xavier.local") {
  source("~/Documents/MinBA_github/minba_00settings.r")
}else{
  source("C:\\Users\\rotllxa\\Desktop\\MinBA_2019/minba_00settings.r")
}

#### Downloading Presence Records ####
if(tolower(pres2bdwnld) != "no"){
  if (data_rep == "gbif")  PreSPickR:::GetBIF(credentials = "gbif_credentials.RData", sp_list = "species_gbif.csv", out_name = "sp_records_gbif")
  if (data_rep == "virtual_sp")  source("MinBA_VirtualSp_iber.R")
  if (data_rep == "bioatles")  PreSPickR:::bioatles(sp_list = "species_bioatles.csv", out_name = "sp_records_bioatles")
}

#### Climatic Data ####
#Downloading data from worldclim.org
if(tolower(clim2bdwnld) != "no"){
  if (data_rep == "gbif"){
    bioclim <- getData('worldclim', var='bio', res = resol, path = paste0(wd))
    save(bioclim, file = "wc5/wc5.RData")
  }
  if (data_rep == "bioatles"){
    bioclim_16 <- getData('worldclim', var='bio', res = resol, lon=3, lat=39,
                          path = paste0(wd))  # importing tile 16
    bioclim <- bioclim_16
    save(bioclim, file = "wc0.5/wc05.RData")
    rm(bioclim_16) ; gc()
  }
  
}else{
  if (data_rep == "gbif")  load(paste0(wd, "/wc5/wc5.RData"), verbose = FALSE)
  if (data_rep == "bioatles")  load(paste0(wd, "/wc0.5/wc05.RData"), verbose = FALSE)
}


#Making a mask (if necessary)
mskng <- "no"
if(tolower(mskng) != "no"){
  msk <- bioclim$bio12
  msk <- reclassify(msk, c(0, msk@data@max, 1)) # reclassify to 1
  }


#### Modelling species ####
library(devtools)
install_github("xavi-rp/MinBAR")
library(MinBAR)
#??MinBAR


if (data_rep == "bioatles"){
  occurrences <- read.csv(paste0(wd, "/sp_records_bioatles.csv"), header = TRUE)
  varbles <- paste0(wd, "/wc0.5")
}else if (data_rep == "gbif"){
  occurrences <- read.csv(paste0(wd, "/sp_records_gbif_native.csv"), header = TRUE)
  varbles <- paste0(wd, "/wc5_other/biovars_eur")
}else if (data_rep == "virtual_sp"){
  #occurrences <- read.csv(paste0(wd, "/sp_records_VirtualSpecies_goo.csv"), header = TRUE)
  #varbles <- paste0(wd, "/wc5_other/biovars_eur")
  #varbles <- paste0("~/Documents/MinBA_models/", "wc5_other/biovars_eur")
  
  virtSp_occurrences <- virtSp_occurrences_good
  #virtSp_occurrences$species <- gsub("virtSp_narrow_", "", virtSp_occurrences$species)
  #virtSp_occurrences$species <- gsub("_alpha-001_", " ", virtSp_occurrences$species)
  #virtSp_occurrences$species <- gsub("_alpha-005_", " ", virtSp_occurrences$species)
  #virtSp_occurrences$species <- gsub("_alpha-005_", " ", virtSp_occurrences$species)
  #virtSp_occurrences$species <- gsub("eta018", "18", virtSp_occurrences$species)
  #virtSp_occurrences$species <- gsub("eta025", "25", virtSp_occurrences$species)
  #virtSp_occurrences$species <- gsub("eta", "", virtSp_occurrences$species)
  #virtSp_occurrences$species <- gsub("sp10", "sp0", virtSp_occurrences$species)
  virtSp_occurrences$species <- gsub("_", " ", virtSp_occurrences$species)  
  virtSp_occurrences$species <- gsub(" ", " sp", virtSp_occurrences$species)  
  virtSp_occurrences$species <- gsub("v", "V", virtSp_occurrences$species)  
  virtSp_occurrences$species <- gsub("sp", "s", virtSp_occurrences$species)
  for(i in 1:9){
    virtSp_occurrences$species <- gsub(paste0("s", i, "$"), paste0("s0", i), virtSp_occurrences$species)
  }
  unique(virtSp_occurrences$species)
  
  occurrences <- virtSp_occurrences
  #varbles <- vrbles[[vrbles_subset]]
  varbles <- vrbles
}
#load(paste0(wd, "/wc5_other/biovars_ent.RData"), verbose = T)
#load(paste0(wd, "/wc5_other/biovars_eur.RData"), verbose = T)
#load(paste0(wd, "/wc0.5/wc05.RData"), verbose = T)


MinBAR:::minba(occ = occurrences,
               varbles = varbles,
               wd = wd,
               prj = 4326,
               num_bands = 10, n_rep = 3,
               maxent_tool = "maxnet")
               #maxent_tool = "dismo")



#### Frequences of Best Buffer ####

best2_bnd_2exp <- read.csv(paste0(dir2save, "/rankingBestBuffer.csv"), header = TRUE)
frec_best_NoTime <- as.data.frame(table(best2_bnd_2exp$Best_Buffer_NoTime))
frec_best_WithTime <- as.data.frame(table(best2_bnd_2exp$Best_Buffer_WithTime))

frec_best <- as.data.frame(matrix(seq(1:num_bands), nrow = num_bands, ncol = 1))
frec_best <- merge(frec_best, frec_best_NoTime, by.x = "V1", by.y = "Var1", all = TRUE)
frec_best <- merge(frec_best, frec_best_WithTime, by.x = "V1", by.y = "Var1", all = TRUE)
frec_best[is.na(frec_best)] <- 0
names(frec_best) <- c("BufferNum", "Frec_Best_NoTime", "Frec_Best_WithTime")

perc1 <- round((prop.table(frec_best$Frec_Best_NoTime)*100), 0)
perc2 <- round((prop.table(frec_best$Frec_Best_WithTime)*100), 0)
perc <- paste0(c(perc1, perc2), "%")
maxY <- max(frec_best_NoTime$Freq, frec_best_WithTime$Freq) + 1
posX <- c(1:num_bands, (num_bands + 2):((2 * (num_bands)) + 1)) - 0.3
posY <- (maxY / 10 * 0.2) + c(frec_best$Frec_Best_NoTime, frec_best$Frec_Best_WithTime)
if (maxY <= 5) legY <- -0.6 else legY <- -2.5  #legY <- -1.5
if (maxY <= 5) posY1 <- -0.4 else posY1 <- -1.7  
if (maxY > 5 & maxY <= 9) legY <- - 1.5 
if (maxY > 5 & maxY <= 9) posY1 <- - 1

palte <- colorRampPalette(colors = c("darkgreen", "green", "yellow", "orange", "red", "darkred"))(num_bands)

#pdf(paste0(dir2save, "/BestBuffers.pdf"))
png(paste0(dir2save, "/BestBuffers.png"))
par(xpd = TRUE, mar = par()$mar + c(3.7, 1, 0, 0))
bpl <- barplot(as.matrix(frec_best[, c(2:3)]),
               main = "Best Buffer With and Without Execution Time", 
               ylab = paste0("Frequencies (n = ", nrow(best2_bnd_2exp), ")"),
               ylim = c(0, maxY),
               beside = TRUE, space = c(0, 1),
               col = palte,
               names.arg = c(1:10, 1:10),
               #names.arg = c("Without Execution Time", "With Execution Time"),
               #cex.names = 1, las = 1)
               cex.names = 0.8, las = 1)
lg <- legend((num_bands/2), legY,
             legend = frec_best$BufferNum,
             fill = palte,
             title = "Buffer Number", cex = 1, 
             ncol = (num_bands/2))
txt1 <- text(x = 6, y = posY1, "Without Execution Time", cex = 1.3)
txt1 <- text(x = 17, y = posY1, "With Execution Time", cex = 1.3)
txt <- text(x = posX, y = posY, perc, cex = 0.8, pos = 4, srt = 45)

dev.off()
#


#### Comparing best model with "real" Virtual Species models ####
res_overlap1 <- as.data.frame(matrix(nrow = 0, ncol = 0))

for(s in 1:length(unique(occurrences$species))){
  sp2check <- unique(occurrences$species)[s]
  
  #Virtual species
  load(paste0(dir2save, "/virtSp_allRasters.RData"), verbose = FALSE)
  #graphics.off()
  #plot(virtSps[[1]]$suitab.raster)
  #virtSps[[1]]$suitab.raster
  #range(virtSps[[1]]$suitab.raster@data@values, na.rm = TRUE)
  vs_rstr <- virtSps[[s]]$probability.of.occurrence
  
  
  #Best model
  sp2check_short <- tolower(paste(substr(sp2check, 1, 3), substr(sub(".* ", "", sp2check), 1, 3), sep = "_"))
  dt2exp_mean <- read.csv(paste0(dir2save, "/results_", sp2check_short, "/info_mod_means_", sp2check_short, ".csv"))
  best_buff <- which(dt2exp_mean$rankFinalNoTime == 1)
  
  
  for(i in 1:10){
    bfr <- i
    
    dir_func <- function(dir2save, sp2check_short, bfr){ 
      res <- tryCatch(
        {
          load(paste0(dir2save, "/results_", sp2check_short, "/model_", sp2check_short,"_", bfr, "_", 1, "/model.RData"), verbose = FALSE)
          #modl
        },
        error = function(con){
          message(con)
          return(NULL)
        }
      )
      if(exists("modl")){ return(modl) }else{ return(NULL) }
    } #end of dir_func
    
    
    
    modl <- dir_func(dir2save, sp2check_short, bfr)
    
    
    #library(ENMeval)
    modl_rstr <- ENMeval::maxnet.predictRaster(mod = modl, env = varbles, 
                                               type = c("logistic"), clamp = TRUE)
    #modl_rstr
    #vs_rstr
    #plot(modl_rstr)
    #plot(vs_rstr)
    
    
    #
    #library(devtools)
    ##install_github("danlwarren/ENMTools", force = TRUE)
    #library(ENMTools)
    res_overlap <- ENMTools::raster.overlap(modl_rstr, vs_rstr, verbose = FALSE)
    res_overlap <- as.data.frame(res_overlap)
    res_overlap$species <- sp2check 
    res_overlap$buffer <- bfr 
    if(i == best_buff){
      res_overlap$best <- "Y"
    }else{
      res_overlap$best <- "N"
    }
    
    res_overlap1 <- rbind(res_overlap1, res_overlap)
    
  }
  
}

View(res_overlap1)
write.csv(res_overlap1, file = paste0(dir2save, "/comp_BestModl_VirtSp.csv"), row.names = FALSE)

#res_overlap1[res_overlap1$best == "Y", ]
#res_overlap1[res_overlap1$species == "b18 sp1", ]

#summary(res_overlap1[res_overlap1$best == "Y", c(1:3)])

summ_overlap <- as.data.frame(matrix(ncol = 0, nrow = 0))
summ_overlap <- rbind(summ_overlap, as.data.frame(t(apply(res_overlap1[res_overlap1$best == "Y", c(1:3)], 2, mean))))
row.names(summ_overlap)[1] <- "mean"
summ_overlap <- rbind(summ_overlap, as.data.frame(t(apply(res_overlap1[res_overlap1$best == "Y", c(1:3)], 2, sd))))
row.names(summ_overlap)[2] <- "sd"
#summ_overlap <- rbind(summ_overlap, as.data.frame(t(apply(res_overlap1[res_overlap1$best == "Y", c(1:3)], 2, function(x) round(length(x), 0)))))
#row.names(summ_overlap)[3] <- "n"
summ_overlap
write.csv(summ_overlap, file = paste0(dir2save, "/comp_BestModl_VirtSp_MeanSD.csv"), row.names = TRUE)






#### GBIF Data sets citation ####
setwd("/Users/xavi_rp/Google Drive/MinBA_OldVersions/gbif_data")

dcs <- list.files(path = "/Users/xavi_rp/Google Drive/MinBA/gbif_data", pattern = "*.RData$", full.names = TRUE)
dcs

lst_rfs <- as.data.frame(matrix(nrow = 0, ncol = 4))
for (d in 1:length(dcs)){
  load(dcs[d], verbose = TRUE)
  if (d == 7){
    d3 <- toupper("Cotoneaster tomentosus Lindl.")
  }else{
    d3 <- rqst_02_meta$request$predicate$predicates[[1]]$value
  }
  d1 <- c(d3, "GBIF", rqst_02_meta$doi, rqst_02_meta$modified)
  lst_rfs[d,] <- d1

}
names(lst_rfs) <- c("species", "source", "DOI", "date_downloaded")

lst_rfs

write.csv(lst_rfs, "citation_data.csv", row.names = FALSE)


