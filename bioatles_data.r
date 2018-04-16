https://www.datacamp.com/community/tutorials/r-web-scraping-rvest



getwd()
setwd("/Users/xavi/Google Drive/MinBA/bioatles")
library(rvest)
library(stringr)
library(purrr)
library(httr)

page <- 'http://bioatles.caib.es/serproesfront/cuadriculas.do?seccion=distribEspecies' %>% read_html()

noms <- as.data.frame(page %>% html_nodes('option') %>% html_text())
View(noms)
chks <- which(grepl("^Selecciona" , noms$`page %>% html_nodes("option") %>% html_text()`))
noms1 <- as.data.frame(noms[(chks[length(chks)-1]+1) : (chks[length(chks)]-1), ])
head(noms1)
nrow(noms1)

kk1 <- page %>% html_nodes('select#selectEspecie')
kk1[[1]]
x <- list(.name = xml_name(kk1[[1]]))
attrs <- xml_attrs(kk1[[1]])
attrs <- attrs[!grepl("xmlns", names(attrs))]
x <- c(x, attrs)
children <- xml_children(kk1[[1]])
library(dplyr)
code <- as.data.frame(bind_rows(lapply(xml_attrs(children), function(x) data.frame(as.list(x), stringsAsFactors=FALSE)))$value)
nrow(code)

#









library(XML)
doc <- htmlParse(page, asText=TRUE)
plain.text <- xpathSApply(doc, "//p", xmlValue)
cat(paste(plain.text, collapse = "\n"))

str(html_text(kk1, trim=F))
html_name(kk1)
html_attr(kk1, "onclick", default = NA_character_)



data.frame(
  Id = kk1[[1]] %>% xml_attrs("select")
)

html_attr(kk1, "select#codiEspecie")

#


pg2 <- "http://bioatles.caib.es/serproesfront/registros.do?accion=listarRegistros&codiEspecie=4582&codiFamilia=0&codiGrupo=0" %>% read_html()
tbl1 <- pg2 %>% html_table(fill = TRUE, header=T)
pres_row <- tbl1[[2]]






page <- 'http://bioatles.caib.es/serproesfront/cuadriculas.do?seccion=distribEspecies' %>% read_html()
sess1<- html_session('http://bioatles.caib.es/serproesfront/cuadriculas.do?seccion=distribEspecies')
sess1

form <- set_values(html_form(page)[[1]], taxonomia = "Chaenorhinum rodriguezii")
form$method

pg1 <- "http://bioatles.caib.es/serproesfront/VisorServlet#" %>% read_html()
sess2<- html_session('http://bioatles.caib.es/serproesfront/VisorServlet#')


kk2<- submit_form(sess1, form)  %>% jump_to()  #, httr::write_disk("chaen.xls",overwrite = TRUE))
kk2<- submit_form(sess1, form)  %>% jump_to(sess2, "http://bioatles.caib.es/serproesfront/VisorServlet#")#, httr::write_disk("chaen.xls",overwrite = TRUE))
html_nodes(kk2, "div.tiendas_resultado_left") 

#station_id <- page %>% html_nodes('select#selectEspecie') %>% html_text() %>% str_trim() %>% unlist()


pg1 <- "http://bioatles.caib.es/serproesfront/VisorServlet#" %>% read_html()
pg1 %>% html_text()
sess2<- html_session('http://bioatles.caib.es/serproesfront/cuadriculas.do?seccion=distribEspecies')
sess2
station_id <- pg1 %>% html_nodes("span") # %>% html_attrs()
html_attrs(station_id)[[2]] <- "Espècie Chaenorhinum rodriguezii"


html %>% 
  html_nodes('time') %>% 
  # The status information is this time a tag attribute
  html_attrs() %>%             
  # Extract the second element
  map(2) %>%                    
  unlist() 













%>% html_text() %>% unlist()

kk1 <-  html_text(kk1, trim=T)

head(kk1)

html_form(page)[[1]]
form <- set_values(html_form(page)[[1]], taxonomia = "Chaenorhinum rodriguezii")
form



a <- unlist(str_split(kk1, '([[:upper:]])'))[-1]
b <- unlist(str_split(kk1, '([[:lower:]])'))
b <- b[b %in% LETTERS]
c <- paste0(b, a)[-1]
head(c)
length(c)
tail(c)
View(c)