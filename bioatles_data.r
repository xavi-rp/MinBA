###########################################################################
########                                                       ############
########               Downloading data from GBIF              ############
########                                                       ############
########                                                       ############
###########################################################################


# bioatles_data.r
#
# Created on: Winter-Spring 2018 (under construction)
#
# Created by: Xavier Rotllan-Puig (xavi.rotllan.puig@gmail.com)
#
# Description: The aim of this script is to download species (presences) data from the Bioatles
#              http://bioatles.caib.es  
# 
# 

# ------------------------------------------

getwd()
wd <- "/Users/xavi/Google Drive/MinBA/bioatles"
setwd(wd)

#### packages ####
library(rvest)
library(stringr)
library(purrr)
library(httr)
library(dplyr)
library(sp)

#### Settings ####
# Output name of the data set
out_name <- "sp_records"

#### Reading in data ####
species <- as.vector(read.csv(paste0(wd, "/species.csv"), header = FALSE)$V1)  # species to be downdloaded
page <- 'http://bioatles.caib.es/serproesfront/cuadriculas.do?seccion=distribEspecies' %>% read_html()  # web site of Bioatles


#### Downloading Species Name and BIOATLES code ####
noms <- as.data.frame(page %>% html_nodes('option') %>% html_text())
chks <- which(grepl("^Selecciona" , noms$`page %>% html_nodes("option") %>% html_text()`))
noms1 <- as.data.frame(noms[(chks[length(chks)-1]+1) : (chks[length(chks)]-1), ])

nds1 <- page %>% html_nodes('select#selectEspecie')
nds1[[1]]
x <- list(.name = xml_name(nds1[[1]]))
attrs <- xml_attrs(nds1[[1]])
attrs <- attrs[!grepl("xmlns", names(attrs))]
x <- c(x, attrs)
children <- xml_children(nds1[[1]])
code <- as.data.frame(bind_rows(lapply(xml_attrs(children), function(x) data.frame(as.list(x), stringsAsFactors=FALSE)))$value)

code_name <- cbind(as.data.frame(code[-1, ]), noms1)
names(code_name) <- c("code", "name")
#code_name[code_name$code == 4582, ]


#### Downloading data of species presences ####
#species <- species[7]
data1 <- data.frame()

for (sps in species) {
  print(paste0("Downloading data: ", sps))
  spec2 <- tolower(paste(substr(sps, 1, 3), substr(sub(".* ", "", sps), 1, 3), sep = "_"))
  #spec <- code_name[grepl("halepensis", code_name$name), ] 
  spec <- code_name[code_name$name %in% sps, ] 
  if (nrow(spec) == 0) stop("No data for this species in Bioatles, please check name/spelling")
  spec_code <- as.vector(spec$code)
  
  #page2 <- "http://bioatles.caib.es/serproesfront/registros.do?accion=listarRegistros&codiEspecie=4582&codiFamilia=0&codiGrupo=0" %>% read_html()
  page2 <- paste0("http://bioatles.caib.es/serproesfront/registros.do?accion=listarRegistros&codiEspecie=", spec_code, "&codiFamilia=0&codiGrupo=0") %>% read_html()
  tbl_sp <- page2 %>% html_table(fill = TRUE, header=T)
  pres2export <- tbl_sp[[2]]
  pres2export <- pres2export[-c(1, ncol(pres2export))]

  write.csv(pres2export, paste0(wd, "/pres_bioatles_", spec2, ".csv"))
  
  data02 <- pres2export[, c(7, 1, 2)]
  data02 <- data02[!duplicated(data02), ]
  data02[, 2:3] <- data02[ ,2:3] * 1000
  names(data02) <- c("species", "x", "y")
  
  #### Reproject presences ####
  coordinates(data02) <- c("x", "y")  # setting spatial coordinates
  proj4string(data02) <- CRS("+init=EPSG:23031")  # define projection: European Datum 1950 (31N)
  #summary(data02)
  CRS.new <- CRS("+init=EPSG:4326") # Lat/Long Geographic Coordinates System WGS84
  data02_WGS84 <- spTransform(data02, CRS.new)  #projecting 
  #summary(data02_WGS84)
  
  data02 <- as.data.frame(data02_WGS84@coords)
  data02$species <- sps
  data02$sp2 <- spec2

  data1 <- rbind(data1, data02)
}
data1 <- data1[, c(2, 1, 3, 4)]

#### Saving data ####
print("Saving Bioatles data as ", paste0(wd, "/", out_name, ".csv"))
write.csv(data1, paste0(wd, "/", out_name, ".csv"), quote = FALSE, row.names = FALSE)
# data1 <- read.csv(paste0(out_name, ".csv"), header = TRUE)










