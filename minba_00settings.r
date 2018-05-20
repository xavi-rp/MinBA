###########################################################################
#
# minba_00settings.r
#
# Created on: Winter 2018 (under construction)
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
library(ENMeval)
#library(rnaturalearth)
library(lattice)
library(latticeExtra)
# My useful functions
source("https://raw.githubusercontent.com/xavi-rp/xavi_functions/master/xavi_functions.r")


#### Settings ####
wd <- "~/Google Drive/MinBA"
setwd(wd)
dir2save <- paste0(wd, "/minba_20180430")    # Europe, etc
dir2save <- paste0(wd, "/minba_20180506")    # Baelarics
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

