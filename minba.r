
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
library(lattice)
library(latticeExtra)
# My useful functions
source("https://raw.githubusercontent.com/xavi-rp/xavi_functions/master/xavi_functions.r")


#### Settings ####
wd <- "~/Google Drive/MinBA"
setwd(wd)
dir2save <- paste0(wd, "/minba_", format(Sys.Date(), format="%Y%m%d"))
if(!file.exists(dir2save)) dir.create(dir2save)


# Need to download presence data from GBIF/Bioatles?
# If != "no", provide a csv with the list of species called "species.csv"
pres2bdwnld <- "no"
data_rep <- "gbif"
data_rep <- "bioatles"

# Resolution
if (data_rep == "gbif") resol <- 5  # 5 arcmin ~ 4.5 km2
if (data_rep == "bioatles") resol <- 0.5  # 30 arcsec ~ 1 km2

# Need to download climatic data?
clim2bdwnld <- "no"

# Number of bands
num_bands <- 10

# n to calculate Boyce Index average
n_times <- 3


#### Retrieving Presence Records ####
if(tolower(pres2bdwnld) != "no"){
  if (data_rep == "gbif")  source(paste0(wd, "/gbif_data.r"))
  if (data_rep == "bioatles")  source(paste0(wd, "/bioatles_data.r"))
}

if (data_rep == "gbif")  presences <- read.csv("gbif_data/sp_records.csv", header = TRUE)
if (data_rep == "bioatles")  presences <- read.csv("bioatles/sp_records.csv", header = TRUE)

colnames(presences)[1:2] <- c("lat", "lon") 
presences <- presences[, c(2,1,3,4)]
#Selecting presences in Europe + Russia + North Africa
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
  if (data_rep == "gbif"){
    bioclim <- getData('worldclim', var='bio', res = resol, path = paste0(wd))
    save(bioclim, file = "wc5/wc5.RData")
  }
  if (data_rep == "bioatles"){
    bioclim_16 <- getData('worldclim', var='bio', res = resol, lon=3, lat=39,
                          path = paste0(wd))  # importing tile 16
    bioclim <- bioclim_16
    save(bioclim, file = "wc0.5/wc05.RData")
  }
}
if (data_rep == "gbif")  load("wc5/wc5.RData", verbose = FALSE)
if (data_rep == "bioatles")  load("wc0.5/wc05.RData", verbose = FALSE)


#Making a mask (if necessary)
mskng <- "no"
if(tolower(mskng) != "no"){
  msk <- bioclim$bio12
  msk <- reclassify(msk, c(0, msk@data@max, 1)) # reclassify to 1
  }

#### Modelling per each species ####
specs <- unique(presences$sp2)
#specs <- unique(presences$sp2)[-c(1:4)]

best2_bnd_2exp <- as.data.frame(matrix(ncol = 0, nrow = 0)) # a table to export rankings of best and 2nd best bandwidth

for(sps in specs){
  pres <- presences[presences$sp2 %in% sps, ] # selecting for species
  specs_long <- as.character(unique(pres$species)) # complete name of the species
  pres <- pres[, c(1,2,4)]
  coordinates(pres) <- c("lon", "lat")  # setting spatial coordinates
  #plot(pres)
  
  #### Calculating the centre of the population, its most distant point and "bands" ####
  #pop_cent <- c(mean(presences$x, na.rm =TRUE), mean(presences$y, na.rm =TRUE))
  #pop_cent <- as.data.frame(matrix(pop_cent, ncol = 2))
  #names(pop_cent) <- c("x", "y")
  
  geocntr <- as.data.frame(geomean(pres))  #mean location for spherical (longitude/latitude) coordinates that deals with the angularity
  #geocntr1 <- geocntr
  #coordinates(geocntr1) <- c("x", "y") 
  #plot(geocntr1, col = "green", add = TRUE)
  
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
  
  #### Croping variables to pres extent  + 5%
  ext1 <- pres@bbox
  incr1 <- apply(ext1, 1, function(x) x[2] - x[1]) * 0.05
  ext1[1, 1] <- pres@bbox[1, 1] - incr1[1]
  ext1[1, 2] <- pres@bbox[1, 2] + incr1[1]
  ext1[2, 1] <- pres@bbox[2, 1] - incr1[2]
  ext1[2, 2] <- pres@bbox[2, 2] + incr1[2]
  varbles1 <- stack(crop(bioclim, ext1))

  # number of background points (see Guevara et al, 2017)
  num_bckgr1 <- (varbles1@ncols * varbles1@nrows) * 50/100 
  
  
  #### Making models for each bandwidth ####
  #table with info to be exported
  dt2exp <- as.data.frame(matrix(ncol = 12, nrow = 0))
  selfinfo2exp <- as.data.frame(matrix(ncol = 9, nrow = 0))
  dt2exp_mean <- as.data.frame(matrix(ncol = 5, nrow = 0))
  
  for (bdw in 1:length(bndwidth)) { # for each bandwidth
    x <- 1
    # set of presences for modeling within the bandwidth
    pres4model <- pres[pres$dist2centr <= bndwidth[bdw], ] 
    
    # croping variables to pres4model extent  + 5%  <-- to fit the model
    ext <- pres4model@bbox
    incr <- apply(ext, 1, function(x) x[2] - x[1]) * 0.05
    ext[1, 1] <- pres4model@bbox[1, 1] - incr[1]
    ext[1, 2] <- pres4model@bbox[1, 2] + incr[1]
    ext[2, 1] <- pres4model@bbox[2, 1] - incr[2]
    ext[2, 2] <- pres4model@bbox[2, 2] + incr[2]
    varbles <- stack(crop(bioclim, ext))
    #dev.off()
    #plot(varbles$bio1)
    #plot(pres, add = TRUE)
    
    # number of background points (see Guevara et al, 2017)
    num_bckgr <- (varbles@ncols * varbles@nrows) * 50/100 
    #if(num_bckgr<100) {
    #  pres4model1 <- sample(1:nrow(pres4model), nrow(pres4model)*0.1) 
    #  pres4model <- pres4model[pres4model1,]
    #}
    
    
    # sampling presences for calibrating and testing (70-30%) within the bandwidth
    folds <- sample(1:nrow(pres4model), nrow(pres4model)*0.7)  
    samp <- as.numeric(unlist(folds))
    pres4cali <- pres4model[samp, 1]
    pres4test <- pres4model[-samp, 1]
    rm(pres4model); gc()
    
    # sampling presences for testing on the whole extent (30% of total presences except those for calibrating)
    pres1 <- pres[-samp, 1]
    folds1 <- sample(1:nrow(pres1), nrow(pres1)*0.3)  
    samp1 <- as.numeric(unlist(folds1))
    pres4test_tot <- pres1[-samp1, 1]

    
    repeat{   # maybe it can be done directly with maxent; if so, we would also have the "average-model"
      t1 <- Sys.time()
      print(paste0("modelling for ", specs_long, " - bandwidth #", bdw, "_", x))

      # Running maxent from dismo 
      if(!file.exists(paste0(dir2save,"/results_", sps))) dir.create(paste0(dir2save,"/results_", sps))
      if(!file.exists(paste0(dir2save,"/results_", sps, "/model_", sps, "_", bdw, "_", x))) dir.create(paste0(dir2save,"/results_", sps, "/model_", sps, "_", bdw, "_", x))
      path <- paste0(dir2save,"/results_", sps,"/model_", sps, "_", bdw, "_", x)
      
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

      #making predictions on the same extent
      preds <- predict(modl, varbles, filename=paste0(path, "/predictions"), progress='text', overwrite=TRUE)

      #make evaluations (on the same extent with 30% to test)
      bg <- randomPoints(varbles, num_bckgr) # background points
      evs <- evaluate(modl, p=pres4test, a=bg, x=varbles)
      save(evs, file = paste0(path, "/evaluations.RData"))
      #load(paste0(path, "/evaluations.RData"), verbose = TRUE)
      
      #Computing Boyce Index (on the same extent with 30% to test)
      byce <- ecospat.boyce(fit = preds, obs = pres4test@coords, nclass=0, window.w="default", res=100, PEplot = TRUE)
      byce$Spearman.cor
      save(byce, file = paste0(path, "/boyce.RData"))
      #load(paste0(path, "/boyce.RData"), verbose = TRUE)
      
      if(bdw != length(bndwidth)){ #except the last bandwith (it makes no sense repeating predictions/evaluations on the same extent)
        
        ## making predictions on the whole species extent
        preds1 <- predict(modl, varbles1, filename=paste0(path, "/predictions_tot"), progress='text', overwrite=TRUE)
        
        #make evaluations
        bg1 <- randomPoints(varbles1, num_bckgr1) # background points
        evs1 <- evaluate(modl, p=pres4test_tot, a=bg1, x=varbles1)
        save(evs1, file = paste0(path, "/evaluations_tot.RData"))
        #load(paste0(path, "/evaluations_tot.RData"), verbose = T)
        
        #Computing Boyce Index
        byce1 <- ecospat.boyce(fit = preds1, obs = pres4test_tot@coords, nclass=0, window.w="default", res=100, PEplot = TRUE)
        byce1$Spearman.cor
        save(byce1, file = paste0(path, "/boyce_tot.RData"))
        #load(paste0(path, "/boyce_tot.RData"), verbose = TRUE)
        
      }else{  # if it is the last bandwidth
        preds1 <- preds
        bg1 <- bg
        evs1 <- evs
        byce1 <- byce
      }
      # gathering info to be exported
      t2 <- Sys.time() - t1
      if(attr(t2, "units") == "hours") {t2 <- t2*60; attr(t2, "units") <- "mins"}
      if(attr(t2, "units") == "secs") {t2 <- t2/60; attr(t2, "units") <- "mins"}
      dt2exp_2 <- as.data.frame(matrix(c(specs_long, paste(bdw, x, sep="_"), bndwidth[bdw], nrow(modl@presence), evs@np, num_bckgr, evs@auc, byce$Spearman.cor, evs1@np, num_bckgr1, evs1@auc, byce1$Spearman.cor), 1, 12, byrow = TRUE))
      dt2exp <- rbind(dt2exp, dt2exp_2)  
      selfinfo2exp_2 <- as.data.frame(matrix(c(specs_long, paste(bdw, x, sep="_"), nrow(pres4cali), nrow(modl@presence), nrow(pres4test), evs@np, num_bckgr, evs@na, t2), 1, 9, byrow = TRUE))
      selfinfo2exp <- rbind(selfinfo2exp, selfinfo2exp_2)
      
      if (x == n_times){ break }else{ x <- x +1 }
    } #end of repeat n times
    
    if(is.null(modl)){ print("jumping to next bandwidth"); next }
    
    print(paste0("computing average for ", specs_long, " - bandwidth #", bdw)) 
    dt2exp[,-c(1:6)] <- data.frame(lapply(dt2exp[-c(1:6)], function(x) as.numeric(as.character(x))))
    dt2exp_m <- mean(dt2exp[(nrow(dt2exp)-n_times+1):nrow(dt2exp), (ncol(dt2exp)-4)], na.rm = TRUE) #mean Boyce partial area
    dt2exp_m2 <- mean(dt2exp[(nrow(dt2exp)-n_times+1):nrow(dt2exp), ncol(dt2exp)], na.rm = TRUE) #mean Boyce whole area
    selfinfo2exp[,-c(1:8)] <- data.frame(lapply(selfinfo2exp[-c(1:8)], function(x) as.numeric(as.character(x))))
    dt2exp_m1 <- mean(selfinfo2exp[(nrow(selfinfo2exp)-n_times+1):nrow(selfinfo2exp), ncol(selfinfo2exp)], na.rm = TRUE)
    dt2exp_mean_2 <- as.data.frame(matrix(c(specs_long, bndwidth[bdw], dt2exp_m, dt2exp_m2, dt2exp_m1), 1, 5, byrow = TRUE))
    dt2exp_mean <- rbind(dt2exp_mean, dt2exp_mean_2)
    dt2exp_mean[, c(2:5)] <- data.frame(lapply(dt2exp_mean[c(2:5)], function(x) as.numeric(as.character(x))))
    
    rm(modl, preds, preds1, bg, evs, byce); gc()

  } # end of for each bandwidth
  
  names(dt2exp_mean) <- c("Species", "Bandwidth", "BoyceIndex_part", "BoyceIndex_tot", "ExecutionTime")
  names(dt2exp) <- c("Species", "ModelNum", "Bandwidth", "numPresencesCalib", "numPresencesTest", "numBackground", "AUC_part", "BoyceIndex", "numPresencesTest_tot", "numBackground_tot", "AUC_tot", "BoyceIndex_tot")
  names(selfinfo2exp) <- c("Species", "ModelNum", "num_pres_calib", "num_pres_calib_used", "num_pres_test", "num_pres_test_used", "num_background", "num_bckgrnd_used", "exec_time")
  
  computing_ranks <- 1
  if (computing_ranks == 1){
    dt2exp_mean$rankBI_part <- rank(-dt2exp_mean$BoyceIndex_part, ties.method = "first")
    dt2exp_mean$rankBI_tot <- rank(-dt2exp_mean$BoyceIndex_tot, ties.method = "first")
    dt2exp_mean$rankTime <- rank(dt2exp_mean$ExecutionTime, ties.method = "first")
    dt2exp_mean$rankFinalNoTime <- rank((dt2exp_mean$rankBI_part + dt2exp_mean$rankBI_tot), ties.method = "first")
    dt2exp_mean$rankFinalWithTime <- rank((dt2exp_mean$rankBI_part + dt2exp_mean$rankBI_tot + dt2exp_mean$rankTime), ties.method = "first")
    
    best2_bnd <- c(sps, 
                   row.names(dt2exp_mean[dt2exp_mean$rankFinalNoTime == 1, ]), 
                   row.names(dt2exp_mean[dt2exp_mean$rankFinalNoTime == 2, ]),
                   row.names(dt2exp_mean[dt2exp_mean$rankFinalWithTime == 1, ]), 
                   row.names(dt2exp_mean[dt2exp_mean$rankFinalWithTime == 2, ]))
    best2_bnd <- as.data.frame(t(best2_bnd))
    names(best2_bnd) <- c("Species", "Best_Bandwidth_NoTime", "SecondBest_bandwidth_NoTime", "Best_Bandwidth_WithTime", "SecondBest_bandwidth_WithTime")
    
    best2_bnd_2exp <- rbind(best2_bnd_2exp, best2_bnd)
    write.csv(best2_bnd_2exp, paste0(dir2save, "/rankingBestBandwidth.csv"), row.names = FALSE)
  }
  
  write.csv(dt2exp_mean, paste0(dir2save, "/results_", sps, "/info_mod_means_", sps, ".csv"), row.names = FALSE)
  write.csv(dt2exp, paste0(dir2save, "/results_", sps, "/info_mod_", sps, ".csv"), row.names = FALSE)
  write.csv(selfinfo2exp, paste0(dir2save, "/results_", sps, "/selfinfo_mod_", sps, ".csv"), row.names = FALSE)
  
  #### Making a plot ####
  graphics.off()
  #dt2exp_mean[,-1] <- data.frame(lapply(dt2exp_mean[-1], function(x) as.numeric(as.character(x))))
  dt2exp_mean[,names(dt2exp_mean) %in% c("BoyceIndex_part", "BoyceIndex_tot")] <- round(dt2exp_mean[,names(dt2exp_mean) %in% c("BoyceIndex_part", "BoyceIndex_tot")], 3)
  pdf(paste0(dir2save, "/results_", sps, "/boyce_bandwidth_", sps, "_part_tot.pdf"))
  plt <- xyplot(BoyceIndex_part ~ Bandwidth, dt2exp_mean,
                #scales = list(y = list(log = 10)),
                type = c("p", "smooth"),
                span = 0.8,
                #type = c("p", "l"), 
                ylim = c(0.8, 1.05),
                col = "blue",
                main = bquote(Boyce~Index~(mean~of~.(n_times)~models)~-~italic(.(specs_long))),
                ylab = "Boyce Index", xlab = "Bandwidth (km)",
                #par.settings = list(par.ylab.text = list(col = "black")),
                #par.settings = simpleTheme(col = 1),
                key=list(#space = "right",
                         x=0.5,y=0.2,
                         lines = list(col=c("blue", "green", "magenta")),
                         text = list(c("Boyce Index Partial","Boyce Index Total", "Execution Time"))
                ))
  plt1 <- xyplot(ExecutionTime ~ Bandwidth, dt2exp_mean,
                 #type = c("p", "smooth"), 
                 type = c("p", "r"), 
                 #ylim = c(0.8, 1.05),
                 ylab = "Execution Time (min)",
                 col = "magenta"
                 #,par.settings = simpleTheme(col = "magenta")
                 )
  dbl_plt <- doubleYScale(plt, plt1, add.ylab2 = TRUE)
  #plot(dbl_plt)
  plt2 <- xyplot(BoyceIndex_tot ~ Bandwidth, dt2exp_mean,
                 type = c("p", "smooth"),
                 span = 0.8,
                 #type = c("p", "l"), 
                 #ylim = c(0.8, 1.05),
                 col = "green")

  plot(dbl_plt + as.layer(plt2))
  dev.off()
  
} # end of loop for sps

#

#### Frequences of Best bandwidth #### 

best2_bnd_2exp <- read.csv(paste0(dir2save, "/rankingBestBandwidth.csv"), header = TRUE)
frec_best_NoTime <- as.data.frame(table(best2_bnd_2exp$Best_Bandwidth_NoTime))
frec_best_WithTime <- as.data.frame(table(best2_bnd_2exp$Best_Bandwidth_WithTime))

frec_best <- as.data.frame(matrix(seq(1:num_bands), nrow = num_bands, ncol = 1))
frec_best <- merge(frec_best, frec_best_NoTime, by.x = "V1", by.y = "Var1", all = TRUE)
frec_best <- merge(frec_best, frec_best_WithTime, by.x = "V1", by.y = "Var1", all = TRUE)
frec_best[is.na(frec_best)] <- 0 
names(frec_best) <- c("BandwidthNum", "Frec_Best_NoTime", "Frec_Best_WithTime")

barchart(Freq ~ Var1, frec_best_NoTime)
barchart(Freq ~ Var1, frec_best_WithTime)

palte <- colorRampPalette(colors = c("darkgreen", "green", "yellow", "orange", "red", "darkred"))(num_bands)

pdf(paste0(dir2save, "/BestBandwidths.pdf"))
par(xpd = TRUE, mar = par()$mar + c(3,1,0,0))
barplot(as.matrix(frec_best[,c(2:3)]), 
        main = "Best Bandwidth With and Without Execution Time", ylab = "Frequencies", 
        ylim = c(0, (max(frec_best_NoTime$Freq, frec_best_WithTime$Freq) + 1)),
        beside = TRUE, space = c(0.1, 1),
        col = palte,
        #names.arg = c(1:10, 1:10),
        names.arg = c("With Execution Time", "Without Execution Time"),
        cex.names = 1, las = 1)
legend((num_bands/2), -1.5, 
       legend = frec_best$BandwidthNum, 
       fill = palte,
       title = "Bandwith Number", cex = 1, ncol = (num_bands/2))

dev.off()
  
  






#### Getting loess (smooth) curve ####
#  l_curve <- loess(BoyceIndex ~ Bandwidth, dt2exp_mean, span = 0.6)    #span=0.8 is default, for smoothing
#  l_preds <- stats::predict(l_curve)
#  plot(l_curve)
#  l_preds
#  l_slope <- diff(l_preds)
#  max_slp <- dt2exp_mean[which.max(l_slope), dt2exp_mean$BoyceIndex]
  # LOESS (SMOOTH)
  
  #DATA
  set.seed(42)
  x = rnorm(20)
  y = rnorm(20)
  x <- dt2exp_mean$Bandwidth
  y <- dt2exp_mean$BoyceIndex_part
  
  #Plot the points
  plot(x, y, type = "p", ylim = c(0.8, 1.2))
  
  #Obtain points for the smooth curve
  temp = loess.smooth(x, y, span = 2/3, family = c("symmetric"), evaluation = 200) 
  temp = loess.smooth(x, y, span = 2/3, family = c("gaussian"), evaluation = 50) 
  temp = loess.smooth(x, y, span = 0.3, family = c("gaussian"), evaluation = 50) 
  
  #Plot smooth curve
  lines(temp$x, temp$y, lwd = 2, col="green")
  
  #Obtain slope of the smooth curve
  slopes <- diff(temp$y)/diff(temp$x)
  slopes <- round(slopes, 5)
  slopes <- ifelse(slopes < 0, NA, slopes)
  min(slopes, na.rm = TRUE)
  order(slopes)
  
  which(slopes == min(slopes, na.rm = TRUE))
  temp$x[15] 
  temp$y[15]
  temp$x[30] 
  temp$y[30]
  slopes[15]
  
  #Add a trend line
  abline(v = temp$x[21], col="pink" )
  abline(v = temp$x[49], col="green" )
  abline(v = temp$x[14], col="yellow" )
  abline(v = temp$x[5], col="red" )
  abline(v = temp$x[30], col="blue" )
  abline(h = temp$y[34], col="blue" )
  abline(lm(y~x))
  
  
  
  
  plot(x, y, type = "p", ylim = c(0.8, 1.35))
  
  y.loess <- loess(y ~ x, span=0.7, data.frame(x=x, y=y))
  y.predict <- predict(y.loess, data.frame(x=x))
  
  lines(x, y.predict, lwd = 2, col="red")
  
  
  peak <- optimize(function(x, model) predict(model, data.frame(x=x)),
                   c(min(x),max(x)),
                   maximum=TRUE,
                   model=y.loess) 
  
  
  points(peak$maximum,peak$objective, pch=FILLED.CIRCLE<-19)
  
  xval <- approx(x = y.predict, y = x, xout = 1)$y    #sembla que és linear
  spline(x = y.predict, y = x, xout = 1) #nonlinear
  
  
  
  
  dev.off()
  
} # end of loop for sps






#plot(bioclim$bio1)
#plot(pres4model, add=T, col = 2)
#plot(pres, add=T)
