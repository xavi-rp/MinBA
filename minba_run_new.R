
###########################################################################
########                                                       ############
########             Minimum Background Area for SDMs          ############
########                          MinBA                        ############
########                                                       ############
###########################################################################
#
# minba.r
#
# Created on: Winter 2018 (under construction)
#
# Created by: Xavier Rotllan-Puig (xavi.rotllan.puig@gmail.com)
#
# Description: The aim of this script is 
# 
# 
# 


# ------------------------------------------
#source("~/Google Drive/MinBA/minba_run.r")

#### Call settings ####
if(Sys.info()[4] == "MacBook-Pro-de-Xavier.local") {
  source("~/Google Drive/MinBA/minba_00settings.r")
}else{
  source("C:\\Users\\rotllxa\\Desktop\\MinBA_2019/minba_00settings.r")
}
#dir2save <- paste0(wd, "/minba_20180430")    # Europe, etc
#dir2save <- paste0(wd, "/minba_20180506")    # Baelarics

#### Downloading Presence Records ####
if(tolower(pres2bdwnld) != "no"){
  if (data_rep == "gbif")  PreSPickR:::GetBIF(credentials = "gbif_credentials.RData", sp_list = "species_gbif.csv", out_name = "sp_records_gbif")
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
}
#if (data_rep == "gbif")  load(paste0(wd, "/wc5/wc5.RData"), verbose = FALSE)
#if (data_rep == "bioatles")  load(paste0(wd, "/wc0.5/wc05.RData"), verbose = FALSE)


#Making a mask (if necessary)
mskng <- "no"
if(tolower(mskng) != "no"){
  msk <- bioclim$bio12
  msk <- reclassify(msk, c(0, msk@data@max, 1)) # reclassify to 1
  }

#### Modelling species ####

#library(devtools)
#install_github("xavi-rp/MinBAR")
#install.packages(c("Hmisc", "rms"))
library(MinBAR)
??MinBAR

#setwd("C:/Users/rotllxa/Desktop/MinBA_2019")
MinBAR:::minba(occ = paste0(wd, "/sp_records_bioatles.csv"), 
               varbles = paste0(wd, "/wc0.5"), 
               prj = "4326", num_bands = 10,
               n_rep = 3, BI_part = NULL, BI_tot = NULL, SD_BI_part = NULL,
               SD_BI_tot = NULL)


#


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
if (maxY <= 5) legY <- -0.4 else legY <- -1.5

palte <- colorRampPalette(colors = c("darkgreen", "green", "yellow", "orange", "red", "darkred"))(num_bands)

#pdf(paste0(dir2save, "/BestBuffers.pdf"))
png(paste0(dir2save, "/BestBuffers.png"))
par(xpd = TRUE, mar = par()$mar + c(3,1,0,0))
bpl <- barplot(as.matrix(frec_best[,c(2:3)]),
               main = "Best Buffer With and Without Execution Time", ylab = paste0("Frequencies (n = ", nrow(best2_bnd_2exp), ")"), 
               ylim = c(0, maxY),
               beside = TRUE, space = c(0, 1),
               col = palte,
               #names.arg = c(1:10, 1:10),
               names.arg = c("Without Execution Time", "With Execution Time"),
               cex.names = 1, las = 1)
lg <- legend((num_bands/2), legY,
             legend = frec_best$BufferNum, 
             fill = palte,
             title = "Buffer Number", cex = 1, ncol = (num_bands/2))
txt <- text(x = posX, y = posY, perc, cex = 0.8, pos = 4, srt = 45)

dev.off()
# 
  

#### Plot with all buffer plots ####



#### BI by groups ####

dt2exp_mean_all <- as.data.frame(matrix(ncol = 0, nrow = 0))

for (sps in specs){
  dt2exp_mean <- read.csv(paste0(dir2save, "/results_", sps, "/info_mod_means_", sps, ".csv"), header = TRUE)[, 1:4]
  dt2exp_mean_all <- rbind(dt2exp_mean_all, dt2exp_mean)
}

BI_SD_bySpecies <- as.data.frame(dt2exp_mean_all %>% group_by(Species) %>% summarise(Mean_BI_part = mean(BoyceIndex_part),
                                                                                     SD_BI_part = sd(BoyceIndex_part),
                                                                                     Mean_BI_tot = mean(BoyceIndex_tot),
                                                                                     SD_BI_tot = sd(BoyceIndex_tot)))
BI_SD_bySpecies
range(BI_SD_bySpecies$Mean_part)
range(BI_SD_bySpecies$SD_part)
range(BI_SD_bySpecies$Mean_tot)
range(BI_SD_bySpecies$SD_tot)



#### Saving Environment Objects ####

save.image(file = "minba_env.RData")

