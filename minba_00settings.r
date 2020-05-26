###########################################################################
#
# minba_00settings.r
#
# Created on: Winter 2018
#
# Created by: Xavier Rotllan-Puig (xavi.rotllan.puig@gmail.com)
#
# Description: The aim of this script is to gather the settings to run minba.r
#
#
# ------------------------------------------


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
library(lattice)
library(latticeExtra)
#library(ENMeval)
#library(rnaturalearth)
# My useful functions
source("https://raw.githubusercontent.com/xavi-rp/xavi_functions/master/xavi_functions.r")
# To download presences: https://www.researchgate.net/publication/326440673_PreSPickR_Downloading_Species_Presences_Occurrences_From_Public_Repositories
library(devtools)
install_github("xavi-rp/PreSPickR")
library(PreSPickR)

#### Settings ####
if(Sys.info()[4] == "MacBook-MacBook-Pro-de-Xavier.local") {
  wd <- "~/Documents/MinBA_models/"
}else{
  wd <- "C:\\Users\\rotllxa\\Desktop\\MinBA_2019"
}
setwd(wd)
#dir2save <- paste0(wd, "/minba_20190622_maxnet")             # Europe, etc
#dir2save <- paste0(wd, "/minba_20190622_maxnet_balears")     # Baelarics
dir2save <- paste0(wd, "/minba_20190910")                    # virtualspecies
#dir2save <- paste0(wd, "/minba_", format(Sys.Date(), format="%Y%m%d"))
if(!file.exists(dir2save)) dir.create(dir2save)


# Need to download presence data from GBIF/Bioatles?
# If != "no", provide a csv with the list of species called "species.csv"
pres2bdwnld <- "yes"
pres2bdwnld <- "no"
data_rep <- "bioatles"
data_rep <- "gbif"
data_rep <- "virtual_sp"

# Resolution
if (data_rep == "gbif") resol <- 5  # 5 arcmin ~ 4.5 km2
if (data_rep == "virtual_sp") resol <- 5  # 5 arcmin ~ 4.5 km2
if (data_rep == "bioatles") resol <- 0.5  # 30 arcsec ~ 1 km2

# Need to download climatic data?
clim2bdwnld <- "no"
#clim2bdwnld <- "yes"

# Number of buffers
num_bands <- 10

# n to calculate Boyce Index average
n_times <- 3

# Conditions to stop the modelling proces
#   If the three are NULL, all buffers are preocessed
#   If one or the two BIs have a value, this or these are the minimum limit to be reached
#   If SD_BI, this or these are the minimum limit to be reached

#SD_min <- 0.006

BI_tot <- NULL
#BI_tot <- 0.90
BI_part <- NULL
#BI_part <- 0.88

SD_BI_part <- NULL
#SD_BI_part <- 0.04
SD_BI_tot <- NULL




