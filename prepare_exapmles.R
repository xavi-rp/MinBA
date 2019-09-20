###########################################################################
########                                                       ############
########             Minimum Background Area for SDMs          ############
########                          MinBA                        ############
########                                                       ############
###########################################################################
#
# prepare_examples.R
#
# Created on: Winter 2018
#
# Created by: Xavier Rotllan-Puig (xavi.rotllan.puig@gmail.com)
#
# Description: The aim of this script is to prepare data set for the 
#              examples to be included in the package
#
#
library(rgdal)
library(sp)
library(raster)

e <- extent(-12, 38, 30,  62)
bioscrop <- crop(bios, e)
bioscrop <- readAll(bioscrop)

sprecords <- read.table("/Users/xavi_rp/Google Drive/MinBAR/data/sprecords.csv", header = T, sep = ";")
load("/Users/xavi_rp/Documents/MinBAR/data/sprecords.RData", verbose = T)
nrow(sprecords[sprecords$decimalLongitude > 37, ])
sprecords <- sprecords[sprecords$decimalLatitude > 30, ]
sprecords <- sprecords[which(extract(is.na(bioscrop$bio1), sprecords) == 1), ]
names(sprecords)
sprecords <- sprecords[-711, ]

write.table(sprecords, file = "/Users/xavi_rp/Google Drive/MinBAR/data/sprecords.csv", sep = ";", row.names = FALSE)

names(sprecords)
sp::coordinates(sprecords) <- c("decimalLongitude", "decimalLatitude")
is.projected(sprecords)
proj4string(sprecords) <- sp::CRS(paste0("+init=EPSG:", 4326))
sprecords@proj4string <- sp::CRS(paste0("+init=EPSG:", 4326))

sprecords <- spTransform(sprecords, CRS("+proj=longlat +ellps=GRS80"))

is.projected(sprecords)

??spTransform


geocntr <- as.data.frame(geosphere::geomean(sprecords))
sprecords$dist2centr <- geosphere::distGeo(sprecords, geocntr) #in meters
sprecords$dist2centr/1000

plot(bioscrop$bio1)
plot(sprecords, add =T)
saveRDS(bioscrop, file = 'bios.rds')

length(which(extract(is.na(bioscrop$bio7), sprecords) == 1))

r <- raster(bioscrop$bio1)

bioscrop <- readRDS("/Users/xavi_rp/Google Drive/MinBAR/Data/bios.rds")
bioscrop
save(bioscrop, file = "/Users/xavi_rp/Google Drive/MinBAR/Data/bios.RData")

library(raster)

load("/Users/xavi_rp/Google Drive/MinBAR/Data/bios.RData")
plot(bioscrop)
bioscrop

vrbles <- readRDS("/Users/xavi_rp/Google Drive/MinBAR/Data/bios.rds")
load("/Users/xavi_rp/Documents/MinBAR/Data/bios.RData", verbose = T)

library(MinBAR)
data(sprecords)
sprecords
data(bios)
bioscrop

MinBAR:::minba(occ = sprecords, varbles = bioscrop, wd = tempdir(), prj = 4326, num_bands = 3)

wd = "/Users/xavi_rp/Google Drive/MinBA/"



load("/Users/xavi_rp/Documents/MinBAR/data/sprecords.RData", verbose = T)
load("/Users/xavi_rp/Documents/MinBAR/Data/bios.RData", verbose = T)

MinBAR:::minba(occ = sprecords, varbles = bioscrop,
                  wd = "/Users/xavi_rp/Google Drive/MinBA/",
                  prj = 4326,
                  num_bands = 3, n_rep = 3,
                  maxent_tool = "maxnet")


MinBAR:::minba(occ = sprecords, varbles = bioscrop,
                  wd = "/Users/xavi_rp/Google Drive/MinBA/",
                  prj = 4326,
                  num_bands = 3, n_rep = 3,
                  maxent_tool = "maxnet")



# Preparing biovars for Eurasia and North Africa 
# (al final, millor fer servir "biovars_ent" perquè evito agafar 
#  tot l'est d'Àsia que pràcticament no hi ha presències, assumint que són introduccions)

wd
rstrs <- list.files("~/Google Drive/MinBA/wc5", pattern = c(".bil$"), full.names = T)
rstrs <- c(rstrs, list.files("~/Google Drive/MinBA/wc5", pattern = c(".tif$"), full.names = T))

#vrbles <- stack()
extens <- extent(biovars_eur)

for(rst in 1:length(rstrs)){
  temp <- raster::raster(rstrs[rst])
  biovars_eur1 <- raster::stack(raster::crop(temp, extens))
  writeRaster(biovars_eur1, 
              filename = paste0("~/Google Drive/MinBA/wc5_other/biovars_eur/", sub(".*/", "", rstrs[rst])) 
              )
}

temp <- raster::raster(paste0("~/Google Drive/MinBA/wc5_other/biovars_eur/bio9.bil"))
plot(temp)

plot(vrbles[[2]])
extent(biovars_eur[[2]])
extens <- extent(vrbles[[2]])

biovars_ent <- vrbles
extens@xmax <- 180
extens@xmin <- -14
extens@ymin <- 28
extens@ymax <- 75
extens

biovars_ent <- raster::stack(raster::crop(biovars_ent, extens))
plot(biovars_ent[[1]])

save(biovars_ent, file = "~/Google Drive/MinBA/wc5_other/biovars_ent.RData")


extens <- extent(biovars_eur[[2]])

extens@xmax <- 60
extens

biovars_eur <- raster::stack(raster::crop(biovars_eur, extens))
plot(biovars_eur[[1]])

save(biovars_eur, file = "~/Google Drive/MinBA/wc5_other/biovars_eur.RData")
load( "~/Google Drive/MinBA/wc5_other/biovars_eur.RData", verbose = TRUE)



# Cutting presences for Eurasia and North Africa (using "biovars_eur")

occurrences <- read.csv(paste0("~/Google Drive/MinBA", "/sp_records_gbif.csv"), header = TRUE)
head(occurrences)
unique(occurrences$species)

colnames(occurrences)[1:2] <- c("lon", "lat")
occurrences <- occurrences[, c(2,1,3)]

sp::coordinates(occurrences) <- c("lon", "lat")  # setting spatial coordinates
occurrences@proj4string <- sp::CRS(paste0("+init=EPSG:", 4326))

extens_biovars_eur <- extent(biovars_eur)
occurrences <- raster::crop(occurrences, extens_biovars_eur)
nrow(occurrences)

occ_coord <- cbind(occurrences@coords, occurrences@data)

occurrences <- read.csv(paste0("~/Google Drive/MinBA", "/sp_records_gbif.csv"), header = TRUE)
names(occurrences)
names(occ_coord) <- names(occurrences)
head(occ_coord)
unique(occ_coord$species)

write.csv(occ_coord, paste0("~/Google Drive/MinBA", "/sp_records_gbif_native.csv"), row.names = FALSE)




