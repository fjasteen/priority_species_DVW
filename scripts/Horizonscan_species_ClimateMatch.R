# Download occurrence data from potentially invasive species from GBIF #
# Libraries #
library(rgbif)
library(tidyverse)
library(ggplot2)
library(INBOtheme)
library(googledrive)
library(readxl)
library(stringr)
library(taxize)
library(readr)
library(googlesheets4)
library(dplyr)
library(rgbif)
library(sp)
library(devtools)
library(trias)

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

# acceptedKeys <- GRIIS_taxonKeys
# taxonkey_set2 <- pred_in("taxonKey", acceptedKeys2)

# initiate download

# set2 <- occ_download(taxonkey_set2,
#                      pred("hasCoordinate", TRUE), # Georefereerde gegevens
#                      #pred_gte("year", 2000), # alle waarnemingen vanaf 2000
#                      user = gbif_user,
#                      pwd = gbif_pwd,
#                      email = gbif_email)

# download data from GBIF
#
# repeat{
#   Sys.sleep(time = 5*length(acceptedKeys2))
#   test_set2 <- occ_download_meta(set2)
#   if(test_set2$status == "SUCCEEDED"){
#     rawdata_set2_imported <- occ_download_get(set2,
#                                               path = "./priorspecies",
#                                               overwrite = TRUE) %>%
#       occ_download_import()
#     break
#   }
#   print(test_set2$status)
# }

# 2. List of taxonIDs from horizonscanningtool CABI #
cabi <-
  read_delim("./data/input/Cabi/Horizon Scanning_20220125154122643.csv")
cabi <- cabi[-c(1:5), ]
colnames(cabi) <- cabi[1, ]
cabi <- cabi[which(cabi$Phylum == 'Spermatophyta' |
                     cabi$Phylum == 'Pteridophyta'), ]

#select with ones with only genus & species name
cabi <- cabi[grep("[ ]", cabi$`Preferred scientific name`),]

# Get GBIF backbone taxon ID from taxonomic names, based on species rank, accepted status, matchtype exact #

library(taxize)
cabi_id <- as.data.frame(
  get_gbifid(
    na.omit(cabi$`Preferred scientific name`),
    ask = TRUE,
    messages = TRUE,
    rows = 1,
    rank = "species",
    method = "backbone",
    sciname = NULL,
  )
)

cabi_taxonID <- cabi_id$ids

# problem: not found: (Senecio squalidus subsp. rupestris) How to solve? Log output and grep not Found, lookup seperately?

#Compare list_taxonID with list from Scheers et al.
Scheers <- read_delim("./data/input/Literature/Lijst_Handel_ScheersK.csv")
Scheers <- Scheers[Scheers$`Taxon category` == "species", ]
Scheers <- Scheers[grep("[ ]{1,}", Scheers$'Horticultural name'),]
Scheers_list <- Scheers$`Horticultural name`

# Get GBIF backbone taxon ID from taxonomic names, based on species rank, accepted status, matchtype exact #
library(taxize)
Scheers_id <- as.data.frame(
  get_gbifid(
    Scheers_list,
    ask = TRUE,
    messages = TRUE,
    rows = 1,
    rank = "species",
    method = "backbone",
    sciname = NULL,
  )
)

Scheers_taxonID <- Scheers_id$ids

# Get GBIF backbone taxon ID from RIPARIAS taxonomic names, based on species rank, accepted status, matchtype exact #
Riparias <- read_delim("./data/input/Riparias/alertlist_RIPARIAS.csv")
Riparias <- Riparias[Riparias$is_plant == TRUE, ]
Riparias_taxonID <- Riparias$gbif_key

# All horizonscan species Scheers et al. + cabi
horizon_taxonID <-
  unique(c(Scheers_taxonID, cabi_taxonID, Riparias_taxonID))

# Remove all potentially indigenous species based on Floradatabank

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

native_taxonID <- native_id$ids


# Duplicate Species in GRIIS list and in horizonscan list?
horizon_taxonKeys <-
  setdiff(horizon_taxonID, GRIIS_taxonKeys) %>% na.omit

horizon_taxonKeys <-
  horizon_taxonKeys[!horizon_taxonKeys %in% native_taxonID]


# taxonkey_set1 <- pred_in("taxonKey", num_occ$horizon_taxonKeys)
#
# # initiate download
#
# set1 <- occ_download(taxonkey_set1,
#                      pred("hasCoordinate", TRUE), # Georefereerde gegevens
#                      pred_gte("year", 1900), # alle waarnemingen vanaf 1900
#                      pred_in("basisOfRecord", c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "UNKNOWN")),
#                      user = gbif_user,
#                      pwd = gbif_pwd,
#                      email = gbif_email)
#
# # download occurrence data from GBIF
# # #download will be done directly from GBIF website,need to fill location into zipfile)
# repeat{
#   Sys.sleep(time = 5*length(acceptedKeys1))
#   test_set1 <- occ_download_meta(set1)
#   if(test_set1$status == "SUCCEEDED"){
#     rawdata_set1_imported <- occ_download_get(set1,
#                                               path = "./data/intermediate/horizonspecies_occ",
#                                               overwrite = TRUE) %>%
#       occ_download_import()
#     break
#   }
#   print(test_set1$status)
# }

# 4. Climate matching horizonscan species

taxonkeys <- horizon_taxonKeys
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
  
  zipfile <- zipfile <- "./data/intermediate/horizonspecies_occ/0193922-210914110416597.zip"

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
          "./data/intermediate/horizon_cm/cm_output_raw.csv")

#r read previous output, eval=FALSE}
output_raw <- read_csv("./data/intermediate/horizon_cm/cm_output_raw.csv")




```{r filter output}
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
```


```{r checks - species not in cm}
species_not_in_cm <- specieslist_2 %>% 
  filter(!GBIF_code %in% output$acceptedTaxonKey)
write_csv(species_not_in_cm,
          "./data/intermediate/species_not_in_cm.csv")
```

```{r export results}
write_csv(output, 
          "./data/output/cm_output.csv")