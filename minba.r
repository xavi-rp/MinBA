
###########################################################################
########                                                       ############
########             Minimum Background Area for SDMs          ############
########                          MinBA                        ############
########                                                       ############
###########################################################################


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
#source("~/Google Drive/MinBA/minba.r")

#### packages ####
options(java.parameters = "-Xmx4g" )
library(sp)
library(rgdal)
library(raster)
#library(graphics)
#library(rgeos)
#dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_162.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
#library(rJava)
library(dismo)
library(ecospat)
library(geosphere)
library(ENMeval)
#library(rnaturalearth)
library(latticeExtra)
# My useful functions
source("https://raw.githubusercontent.com/xavi-rp/xavi_functions/master/xavi_functions.r")


#### Settings ####
wd <- "~/Google Drive/MinBA"
setwd(wd)
#load(paste0(wd, "/.RData")

# Need to download presence data from GBIF?
# If so, provide a csv with the list of species called "species.csv"
pres2bdwnld <- "no"

# Resolution
resol <- 5  # 5 arcmin ~ 4.5 km2
# Need to download climatic data?
clim2bdwnld <- "no"

# Number of bands
num_bands <- 10

# n to calculate Boyce Index average
n_times <- 3


#### Retrieving Presence Records ####
if(tolower(pres2bdwnld) != "no"){
  source(paste0(wd, "/gbif_data.r"))
  }
presences <- read.csv("gbif_data/sp_records.csv", header = TRUE)
colnames(presences)[1:2] <- c("lat", "lon") 
presences <- presences[, c(2,1,3,4)]
#Selecting presences in Europe+Russia
presences1 <- presences[presences$lon >= -10 & presences$lon <= 180, ]
presences1 <- presences1[presences1$lat >= 27 & presences1$lat <= 82, ]
#coordinates(presences) <- c("lon", "lat")  # setting spatial coordinates
#coordinates(presences1) <- c("lon", "lat")  # setting spatial coordinates
#world <- rnaturalearth::countries110
#msk_eur <- world[world$continent == "Europe",] #to know the coordinates
#msk_eur <- world[world$names == "Portugal",] #to know the coordinates
presences <- presences1
rm(presences1)
gc()

#### Climatic Data ####
#Downloading data from worldclim.org
if(tolower(clim2bdwnld) != "no"){
  bioclim <- getData('worldclim', var='bio', res = resol, path = paste0(wd, "/wc20_5m_biovars"))
  save(bioclim, file = "wc20_5m_biovars/wc5.RData")
}
load("wc20_5m_biovars/wc5.RData", verbose = TRUE)

#Making a mask (if necessary)
mskng <- "no"
if(tolower(mskng) != "no"){
  msk <- bioclim$bio12
  msk <- reclassify(msk, c(0, msk@data@max, 1)) # reclassify to 1
  }

#### Modelling per each species ####
#specs <- unique(presences$sp2)
specs <- unique(presences$sp2)[10:13]

for(sps in specs){
  pres <- presences[presences$sp2 %in% sps, ] # selecting for species
  specs_long <- as.character(unique(pres$species)) # complete name of the species
  pres <- pres[, c(1,2,4)]
  coordinates(pres) <- c("lon", "lat")  # setting spatial coordinates
  
  #### Calculating the centre of the population, its most distant point and "bands" ####
  #pop_cent <- c(mean(presences$x, na.rm =TRUE), mean(presences$y, na.rm =TRUE))
  #pop_cent <- as.data.frame(matrix(pop_cent, ncol = 2))
  #names(pop_cent) <- c("x", "y")
  
  geocntr <- as.data.frame(geomean(pres))  #mean location for spherical (longitude/latitude) coordinates that deals with the angularity
  
  pres$dist2centr <- distGeo(pres, geocntr) #in meters
  pres$dist2centr <- pres$dist2centr/1000   #in km
  furthest <- max(pres$dist2centr) 

  # by % of presences equally distributed
  bndwidth <- as.vector(quantile(pres$dist2centr, probs = seq(0, 1, 1/num_bands), names = TRUE))
  bndwidth <- c(bndwidth[2:(length(bndwidth)-1)], furthest) # Not defined by distance, but by % of presences equally distributed
                                                            # This is particularly useful for very discontinuous distributions (e.g. introduced or invasive species),
                                                            # while not affecting more aggregated populations
  # by equally distant bandwidths
  #bndwidth <- as.vector(seq(0, furthest, furthest/num_bands))[-1]
  
  
  #### Making models for each bandwidth ####
  #table with info to be exported
  dt2exp <- as.data.frame(matrix(ncol = 8, nrow = 0))
  selfinfo2exp <- as.data.frame(matrix(ncol = 9, nrow = 0))
  dt2exp_mean <- as.data.frame(matrix(ncol = 4, nrow = 0))
  
  for (bdw in 1:length(bndwidth)) { # for each bandwidth
    x <- 1
    repeat{   # maybe it can be done directly with maxent; if so, we would also have the "average-model"
      t1 <- Sys.time()
      print(paste0("modelling for ", specs_long, " - bandwidth #", bdw, "_", x))
      # set of presences
      pres4model <- pres[pres$dist2centr <= bndwidth[bdw], ] 

      # croping variables to pres4model extent  + 5%
      ext <- pres4model@bbox
      incr <- apply(ext, 1, function(x) x[2] - x[1]) * 0.05
      ext[1, 1] <- pres4model@bbox[1, 1] - incr[1]
      ext[1, 2] <- pres4model@bbox[1, 2] + incr[1]
      ext[2, 1] <- pres4model@bbox[2, 1] - incr[2]
      ext[2, 2] <- pres4model@bbox[2, 2] + incr[2]
      varbles <- stack(crop(bioclim, ext))
      
      # number of background points (see Guevara et al, 2017)
      num_bckgr <- (varbles@ncols * varbles@nrows) * 50/100
      #if(num_bckgr<100) {
      #  pres4model1 <- sample(1:nrow(pres4model), nrow(pres4model)*0.1) 
      #  pres4model <- pres4model[pres4model1,]
      #}

      # sampling presences for calibrating and predicting (70-30%)
      folds <- sample(1:nrow(pres4model), nrow(pres4model)*0.7)  
      samp <- as.numeric(unlist(folds))
      pres4cali <- pres4model[samp, 1]
      pres4test <- pres4model[-samp, 1]
      rm(pres4model); gc()

      # Running maxent from dismo 
      if(!file.exists(paste0(wd,"/results_", sps))) dir.create(paste0(wd,"/results_", sps))
      if(!file.exists(paste0(wd,"/results_", sps, "/model_", sps, "_", bdw, "_", x))) dir.create(paste0(wd,"/results_", sps, "/model_", sps, "_", bdw, "_", x))
      path <- paste0(wd,"/results_", sps,"/model_", sps, "_", bdw, "_", x)
      
      dir_func <- function(varbles, pres4cali, num_bckgr, path){ # to avoid stop modelling if low number of background points or other errors
        res <- tryCatch(
          {
            modl <- maxent(varbles, pres4cali, removeDuplicates = TRUE,
                           nbg=num_bckgr)
            if(exists("modl")) save(modl, file = paste0(path, "/model.RData"))
          },
          error = function(con){
            message(con)
            return(NULL)
          }
        )
        if(exists("modl")){ return(modl) }else{ return(NULL) } 
      } #end of dir_func
      
      modl <- dir_func(varbles, pres4cali, num_bckgr, path)
      if(is.null(modl)){ break }

      #making predictions
      preds <- predict(modl, varbles, filename=paste0(path, "/predictions"), progress='text', overwrite=TRUE)
      #plot(preds1)
      
      #make evaluations
      bg <- randomPoints(varbles, num_bckgr) # background points
      evs <- evaluate(modl, p=pres4test, a=bg, x=varbles)
      save(evs, file = paste0(path, "/evaluations.RData"))

      #Computing Boyce Index
      byce <- ecospat.boyce(fit = preds, obs = pres4test@coords, nclass=0, window.w="default", res=100, PEplot = TRUE)
      byce
      byce$Spearman.cor
      save(byce, file = paste0(path, "/boyce.RData"))
      
      # gathering info to be exported
      t2 <- Sys.time() - t1
      if(attr(t2, "units") == "hours") {t2 <- t2*60; attr(t2, "units") <- "mins"}
      if(attr(t2, "units") == "secs") {t2 <- t2/60; attr(t2, "units") <- "mins"}
      dt2exp_2 <- as.data.frame(matrix(c(specs_long, paste(bdw, x, sep="_"), bndwidth[bdw], nrow(modl@presence), evs@np, num_bckgr, evs@auc, byce$Spearman.cor), 1, 8, byrow = TRUE))
      dt2exp <- rbind(dt2exp, dt2exp_2)  
      selfinfo2exp_2 <- as.data.frame(matrix(c(specs_long, paste(bdw, x, sep="_"), nrow(pres4cali), nrow(modl@presence), nrow(pres4test), evs@np, num_bckgr, evs@na, t2), 1, 9, byrow = TRUE))
      selfinfo2exp <- rbind(selfinfo2exp, selfinfo2exp_2)
      
      if (x == n_times){ break }else{ x <- x +1 }
    } #end of repeat n times
    
    if(is.null(modl)){ print("jumping to next bandwidth"); next }
    
    print(paste0("computing average for ", specs_long, " - bandwidth #", bdw)) 
    dt2exp[,-c(1:6)] <- data.frame(lapply(dt2exp[-c(1:6)], function(x) as.numeric(as.character(x))))
    dt2exp_m <- mean(dt2exp[(nrow(dt2exp)-n_times+1):nrow(dt2exp), ncol(dt2exp)], na.rm = TRUE)
    selfinfo2exp[,-c(1:8)] <- data.frame(lapply(selfinfo2exp[-c(1:8)], function(x) as.numeric(as.character(x))))
    dt2exp_m1 <- mean(selfinfo2exp[(nrow(selfinfo2exp)-n_times+1):nrow(selfinfo2exp), ncol(selfinfo2exp)], na.rm = TRUE)
    dt2exp_mean_2 <- as.data.frame(matrix(c(specs_long, bndwidth[bdw], dt2exp_m, dt2exp_m1), 1, 4, byrow = TRUE))
    dt2exp_mean <- rbind(dt2exp_mean, dt2exp_mean_2)
    rm(modl, preds, bg, evs, byce); gc()

  } # end of for each bandwidth
  
  names(dt2exp_mean) <- c("Species", "Bandwidth", "BoyceIndex", "ExecutionTime")
  names(dt2exp) <- c("Species", "ModelNum", "Bandwidth", "numPresencesCalib", "numPresencesTest", "numBackground", "AUC", "BoyceIndex")
  names(selfinfo2exp) <- c("Species", "ModelNum", "num_pres_calib", "num_pres_calib_used", "num_pres_test", "num_pres_test_used", "num_background", "num_bckgrnd_used", "exec_time")
  write.csv(dt2exp_mean, paste0(wd, "/results_", sps, "/info_mod_means_", sps, ".csv"), row.names = FALSE)
  write.csv(dt2exp, paste0(wd, "/results_", sps, "/info_mod_", sps, ".csv"), row.names = FALSE)
  write.csv(selfinfo2exp, paste0(wd, "/results_", sps, "/selfinfo_mod_", sps, ".csv"), row.names = FALSE)
  
  #### Making a plot ####
  graphics.off()
  dt2exp_mean[,-1] <- data.frame(lapply(dt2exp_mean[-1], function(x) as.numeric(as.character(x))))
  dt2exp_mean[,names(dt2exp_mean) %in% "BoyceIndex"] <- round(dt2exp_mean[,names(dt2exp_mean) %in% "BoyceIndex"], 3)
  pdf(paste0(wd, "/results_", sps, "/boyce_bandwidth_", sps, ".pdf"))
  plt <- xyplot(BoyceIndex ~ Bandwidth, dt2exp_mean,
                #scales = list(y = list(log = 10)),
                type = c("p", "smooth"), 
                #type = c("p", "l"), 
                ylim = c(0.8, 1.05),
                main = paste0("Boyce Index (mean of ", n_times, " models) - ", specs_long),
                ylab = "Boyce Index", xlab = "Bandwidth (km)")
  plt1 <- xyplot(ExecutionTime ~ Bandwidth, dt2exp_mean,
                 type = c("p", "smooth"), 
                 #type = c("p", "l"), 
                 #ylim = c(0.8, 1.05),
                 ylab = "Execution Time (min)" )
  dbl_plt <- doubleYScale(plt, plt1, add.ylab2 = TRUE)
  plot(dbl_plt)
  dev.off()
  
  #### Getting loess (smooth) curve ####
  l_curve <- loess(BoyceIndex ~ Bandwidth, dt2exp_mean, span = 0.6)    #span=0.8 is default, for smoothing
  l_preds <- stats::predict(l_curve)
  plot(l_curve)
  l_preds
  l_slope <- diff(l_preds)
  max_slp <- dt2exp_mean[which.max(l_slope), dt2exp_mean$BoyceIndex]
  
} # end of loop for sps






#plot(bioclim$bio1)
#plot(pres4model, add=T, col = 2)
#plot(pres, add=T)
