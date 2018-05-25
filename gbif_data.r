
###########################################################################
########                                                       ############
########               Downloading data from GBIF              ############
########                                                       ############
########                                                       ############
###########################################################################


# gbif_data.r
#
# Created on: Winter 2018 (under construction)
#
# Created by: Xavier Rotllan-Puig (xavi.rotllan.puig@gmail.com)
#
# Description: The aim of this script is to define the function GetBIF(), which
#              is used to download species occurrences from GBIF (Global Biodiversity
#              Information Facility ), and saves them as a csv data set.
#              It is based on several functions included in the package "rgbif"
#              (Chamberlain, 2017).
#              GetBIF() retrieve your GBIF credentials (user and password) and 
#              automatically checks in a loop until the request of data made to
#              GBIF is ready and starts the download. Finally, it saves the data 
#              in a csv file.
# 
# Inputs:
#       - The location of a csv file with the names of the species to be downloaded
#       - For security reasons, your GBIF credentials (user, password and email) 
#         can be loaded from a RData file (location needs to be given). Otherwise, 
#         they can be passed as arguments
#       
# 
# Outputs:
#       - A csv file with 3 columns (species, decimalLatitude, decimalLongitude)
#         
# References:
#       - Scott Chamberlain (2017). rgbif: Interface to the Global 'Biodiversity' 
#         Information Facility 'API'. R package version 0.9.8. 
#         https://CRAN.R-project.org/package=rgbif


# ------------------------------------------


#### packages ####
library(sp)
library(rgdal)
library(raster)
library(rgbif)
# My useful functions
source("https://raw.githubusercontent.com/xavi-rp/xavi_functions/master/xavi_functions.r")

#### Settings ####
getwd()
wd <- "~/Google Drive/MinBA/gbif_data"
setwd(wd)
#load(paste0(wd, "/.RData")

# Calling GBIF credentials
load("~/Google Drive/GBIF_credentials/gbif_credentials.RData", verbose = TRUE)

# List of species to be downloaded
specs <- as.vector(read.csv("species.csv", header = FALSE)[,1])
specs <- toupper(specs)

# Output name of the data set
out_name <- "sp_records"


#### Downloading Data ####
## Spin up a download request for SEVERAL species data

for (sps in specs){
  print(paste0("Downloading data for ", sps))
  rqst_02 <- occ_download(paste0("scientificName = ", sps), "hasCoordinate = TRUE",
                          type = "and", user = gbif_usr, pwd = gbif_pwrd, email = email)    #prepares the spin up
  # Creates metadata
  rqst_02_meta <- data.frame(status = "INITIAL")
  round <- 1
  while (rqst_02_meta$status != "SUCCEEDED") {
    cat("\r", paste0("round = ", round, " / ", "status = ", rqst_02_meta$status))
    Sys.sleep(60)
    round <- round + 1
    rqst_02_meta <- rqst_02 %>% occ_download_meta
  }
  
  # Start download when meta says "Status: SUCCEEDED"
  dta <- occ_download_get(key = rqst_02_meta$key, path = ".", overwrite = TRUE, curlopts = list(verbose=TRUE))
  
  # saving citation
  citation_02 <- dta %>% gbif_citation
  
  # Saving download info
  save(list = c("rqst_02", "rqst_02_meta", "dta", "citation_02"), file = paste0("download_info_", sps, ".RData"))
  
}


#### Retrieving Data ####
data1 <- data.frame()

for (sps in specs){
  cat(paste0("Reading data for ", sps), "\n")
  load(paste0("download_info_", sps, ".RData"), verbose = TRUE)
  
  # Reading in data
  #load("download_info.RData", verbose = TRUE)
  data02 <- occ_download_import(dta)
  data02 <- data02[!duplicated(data02[,c(133:134)]), ]
  data02 <- data02[, names(data02) %in% c("species", "decimalLatitude", "decimalLongitude")]
  
  data1 <- rbind(data1, data02)
}

data1 <- as.data.frame(data1)  #data set with coordinates and name of species
data1$sp2 <- tolower(paste(substr(data1$species, 1, 3), substr(sub("^\\S+\\s+", '', data1$species), 1, 3), sep = "_"))
#head(data1)
#unq(data1$species, s=F)
#nasbycols(data1)
#table(data1$species)



#### Saving data ####
print("Saving GBIF data as ", paste0(wd, "/", out_name, ".csv"))
write.csv(data1, paste0(wd, "/", out_name, ".csv"), quote = FALSE, row.names = FALSE)
# data1 <- read.csv(paste0(out_name, ".csv"), header = TRUE)
