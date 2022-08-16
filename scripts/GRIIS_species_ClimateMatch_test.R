library(rgbif)
library(tidyverse)
library(ggplot2)
library(INBOtheme)
library(googledrive)
library(readxl)
library(taxize)
# Credentials #
# in .Renviron save: # gbif_user = "gbif username", gbif_pwd = "gbif password", gbif_email = "gbif email"

gbif_user <- Sys.getenv("gbif_user")
gbif_pwd <- Sys.getenv("gbif_pwd")
gbif_email <- Sys.getenv("gbif_email")

# 1. List of taxonIDs from species in GRIIS checklist Belgium (https://www.gbif.org/dataset/6d9e952f-948c-4483-9807-575348147c7e) ###
setwd("C:/Users/frederique_steen/Documents/Github/priority_species_DVW/")
df <- read_delim("./data/input/GRIIS/dwca-unified-checklist-v1.12/taxon.txt")
df <- df[df$kingdom == "Plantae", ] %>%
  mutate_at("taxonID", str_replace, "https://www.gbif.org/species/", "")
GRIIS_taxonKeys <- df$taxonID
df$speciesname <- str_c(df$genus, " ", df$specificEpithet)
GRIIS_list <- df[, c("taxonID", "scientificName", "speciesname")]

#Filter all native taxa
Flora_native <-
  read_delim("./data/input/Floradatabank/IndigenePlanten_Floradatabank.csv")
Flora_native <- filter(Flora_native, Omschrijving1 != "neofyt")
native_list <- pull(Flora_native, NaamWetenschappelijk1)
native_list <- word(native_list, 1, 2, sep = " ")
native_list <-
  native_list[!native_list %in% c('Crataegus x', 'Scirpus x', 'Taraxacum Wiggers', 'NA')]

native_id <- as.data.frame(
  get_gbifid(
    native_list,
    ask = TRUE,
    messages = TRUE,
    rows = 1,
    rank = "species",
    method = "backbone",
    sciname = NULL,
  )
)

#Split the GBIF download dataset (~10Gb zipped) in 4, as to be able to run script.

native_taxonID <- native_id$ids
GRIIS_list <- GRIIS_list[which(!GRIIS_list$taxonID %in% native_taxonID ),]


#split dataset into 4 based on number of occurrences
n <- 5
GRIIS_split<-split(GRIIS_list$taxonID, sample(1:n, nrow(GRIIS_list), replace=T))
for (i in 1:length(GRIIS_split)){
        name <- paste0("set_",i)
        assign(name, pred_in("taxonKey", GRIIS_split[[i]]))
}
taxonKey_list<-lapply(ls(pattern="set_"), get)
# # initiate download
for (i in 1:length(taxonKey_list)){
                    
                      occ_download(taxonKey_list[[i]],
                      pred("hasCoordinate", TRUE), # Georefereerde gegevens
                      pred_gte("year", 1900), # alle waarnemingen vanaf 2000 (ervan uitgaand dat deze correcter zijn)
                      pred_in("basisOfRecord", c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "UNKNOWN")),
                      pred("hasGeospatialIssue", FALSE),
                      user = gbif_user,
                      pwd = gbif_pwd,
                      email = gbif_email)
}
# # # # download occurrence data from GBIF
 # # # #download will be done directly from GBIF website,need to fill location into zipfile)
 repeat{
   Sys.sleep(time = 5*length(taxonkey_set1))
   test_set1 <- occ_download_meta(set1)
   if(test_set1$status == "SUCCEEDED"){
   rawdata_set1_imported <- occ_download_get(set1,
                                               path = "./data/intermediate/GRIIS_occ",
                                               overwrite = TRUE) %>%
       occ_download_import()
     break
   }
   print(test_set1$status)
 }
} 


# perform climate match
library(readr)
library(googlesheets4)
library(dplyr)
library(rgbif)
library(sp)
library(devtools)
library(trias)

taxonkeys <- GRIIS_list$taxonID
occ_counts <- data.frame() %>% 
  mutate(GBIF_code = NA_integer_,
         n = NA_integer_
  )

for(i in taxonkeys){
  count <- occ_count(taxonKey = i,
                     georeferenced = TRUE,
                     from = 1900)
  
  occ_counts <- occ_counts %>% 
    add_row(GBIF_code = as.integer(i),
            n = count) %>% 
    arrange(desc(n))
}

max(occ_counts$n)/sum(occ_counts$n, na.rm = TRUE)

species_name=c()
for(i in 1:length(occ_counts$GBIF_code)){
  sub <-gbif_name_usage(occ_counts$GBIF_code[i])
  species_name[i] <-sub$scientificName
}
specieslist_1<-cbind(species_name,occ_counts)


limit <- 1500000
x <- 1
sum <- 0
specieslist_2 <- data.frame()

for(y in specieslist_1$GBIF_code){
  temp <- specieslist_1 %>% 
    filter(GBIF_code == y)
  
  sum <- sum + temp$n
  if(sum < limit){
    temp <- temp %>% 
      mutate(group = x)
  }else{
    sum <- temp$n
    x <- x + 1 
    temp <- temp %>% 
      mutate(group = x)
  }
  if(nrow(specieslist_2) == 0){
    specieslist_2 <- temp
  }else{
    specieslist_2 <- rbind(specieslist_2, temp)
  }
}
test <- specieslist_2 %>% 
  group_by(group) %>% 
  summarise(sum(n))


#r climate match}
groups <- unique(specieslist_2$group)
output_raw <- data.frame()
for(a in groups){
  taxonkeys <- specieslist_2 %>% 
    filter(group == a) %>% 
    distinct(GBIF_code)
  
  taxonkeys <- taxonkeys$GBIF_code
  
  zipfile <- "./data/intermediate/GRIIS_occ/0402736-210914110416597.zip"

  cm_output <- climate_match(region = "belgium",
                             taxon_key = taxonkeys,
                             scenario = "all", #as opposed to 2026-2050-A1FI
                             zip_file = zipfile,
                             n_limit = 90,
                             cm_limit = 0.2,
                             maps = FALSE)
  output_raw <- rbind(output_raw, cm_output$cm)
  remove(cm_output)
  gc(reset = TRUE)
}

write_csv(output_raw, 
          "./data/intermediate/GRIIS_cm/cm_output_raw.csv")

#r read previous output, eval=FALSE}
output_raw <- read.csv("./data/output/climate_match/cm_output_raw.csv")





output <- output_raw %>% 
  arrange(acceptedTaxonKey) %>% 
  mutate(perc_climate = round(perc_climate*100, 2)) %>% 
  select(acceptedTaxonKey, acceptedScientificName, Classification, perc_climate, n_totaal) %>% 
  pivot_wider(id_cols = c(acceptedTaxonKey,
                          acceptedScientificName,
                          n_totaal),
              names_from = Classification,
              values_from = perc_climate,
              values_fn = first) %>% 
  select(-Csa) %>% 
  filter(Cfb >= 25 | Cfa >= 25)



#{r checks - species not in cm}
species_not_in_cm <- specieslist_2 %>% 
  filter(!GBIF_code %in% output$acceptedTaxonKey)
write_csv(species_not_in_cm,
          "./data/intermediate/species_not_in_cm.csv")


#{r export results}
write_csv(output, 
          "./data/output/cm_output.csv")