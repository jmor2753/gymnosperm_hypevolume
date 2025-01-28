rm(list = ls())

setwd("~/Documents/lab/gymnosperm_niches")

library(dplyr)
library(CoordinateCleaner)
library(raster)
library(ape)
library(stringi)
library(hypervolume)
library(tidyr)
library(purrr)
library(rnaturalearth)
library(rcartocolor)
library(ggplot2)
library(patchwork)

#envar <- raster::getData("worldclim", var = "tmin", res = 2.5)
#envar <- raster::getData("worldclim", var = "bio", res = 2.5)

tree <- read.tree("data/stull&al2021_supermatrix_dated.tre")

tax <- read.csv("data/classification.csv", sep = "\t") #WFO 2023
tax <- tax[tax$taxonRank == "species" & tax$taxonomicStatus == "Accepted", ]
tax$scientificName <- stri_replace_all_fixed(tax$scientificName, " ", "_")

## DO NOT RUN! (too long) - skip to line 265
path_envs <- list.files(path = "wc2-5", pattern='.bil$', recursive = TRUE,
                        full.names = TRUE)
envs <- stack(path_envs)

cyc_data <- readr::read_tsv("data/Cycadopsida.txt")
gin_data <- readr::read_tsv("data/Ginkgoopsida.txt")
gne_data <- readr::read_tsv("data/Gnetopsida.txt")
pin_data <- readr::read_tsv("/Volumes/Personal/lab/araucaria_niches/data/Pinopsida.txt")

## outgroups
amb_data <- readr::read_tsv("data/Amborellaceae.txt")
cal_data <- readr::read_tsv("data/Calycanthaceae.txt")
psi_data <- readr::read_tsv("data/Psilotaceae.txt")
mat_data <- readr::read_tsv("data/Matoniaceae.txt")

cyc_clean <- cyc_data %>%
  setNames(tolower(names(.))) %>%
  filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% 
  #filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN", "PRESERVED_SPECIMEN")) %>%
  filter(basisofrecord %in% c("HUMAN_OBSERVATION")) %>%
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% 
  filter(decimallatitude != 0 | decimallongitude != 0) %>%
  #filter(is.na(as.factor(nwords(locality)))) %>%
  cc_cen(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove country centroids within 2km 
  cc_cap(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove capitals centroids within 2km
  cc_inst(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove zoo and herbaria within 2km 
  cc_sea(lon = "decimallongitude", lat = "decimallatitude") %>% # remove from ocean
  distinct(decimallongitude, decimallatitude, specieskey,datasetkey, .keep_all = TRUE) %>%
  glimpse() 

cyc_clean$species <- stri_replace_all_fixed(cyc_clean$species, " ", "_")

gin_clean <- gin_data  %>%
  setNames(tolower(names(.))) %>%
  filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% 
  #filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN", "PRESERVED_SPECIMEN")) %>%
  filter(basisofrecord %in% c("HUMAN_OBSERVATION")) %>%
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% 
  filter(decimallatitude != 0 | decimallongitude != 0) %>%
  #filter(is.na(as.factor(nwords(locality)))) %>%
  cc_cen(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove country centroids within 2km 
  cc_cap(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove capitals centroids within 2km
  cc_inst(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove zoo and herbaria within 2km 
  cc_sea(lon = "decimallongitude", lat = "decimallatitude") %>% # remove from ocean
  distinct(decimallongitude, decimallatitude, specieskey,datasetkey, .keep_all = TRUE) %>%
  glimpse() 

gin_clean$species <- stri_replace_all_fixed(gin_clean$species, " ", "_")

gne_clean <- gne_data  %>%
  setNames(tolower(names(.))) %>%
  filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% 
  #filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN", "PRESERVED_SPECIMEN")) %>%
  filter(basisofrecord %in% c("HUMAN_OBSERVATION")) %>%
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% 
  filter(decimallatitude != 0 | decimallongitude != 0) %>%
  #filter(is.na(as.factor(nwords(locality)))) %>%
  cc_cen(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove country centroids within 2km 
  cc_cap(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove capitals centroids within 2km
  cc_inst(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove zoo and herbaria within 2km 
  cc_sea(lon = "decimallongitude", lat = "decimallatitude") %>% # remove from ocean
  distinct(decimallongitude, decimallatitude, specieskey,datasetkey, .keep_all = TRUE) %>%
  glimpse() 

gne_clean$species <- stri_replace_all_fixed(gne_clean$species, " ", "_")

pin_clean <- pin_data  %>%
  setNames(tolower(names(.))) %>%
  filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% 
  #filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN", "PRESERVED_SPECIMEN")) %>%
  filter(basisofrecord %in% c("HUMAN_OBSERVATION")) %>%
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% 
  filter(decimallatitude != 0 | decimallongitude != 0) %>%
  #filter(is.na(as.factor(nwords(locality)))) %>%
  cc_cen(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove country centroids within 2km 
  cc_cap(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove capitals centroids within 2km
  cc_inst(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove zoo and herbaria within 2km 
  cc_sea(lon = "decimallongitude", lat = "decimallatitude") %>% # remove from ocean
  distinct(decimallongitude, decimallatitude, specieskey,datasetkey, .keep_all = TRUE) %>%
  glimpse() 

pin_clean$species <- stri_replace_all_fixed(pin_clean$species, " ", "_")

amb_clean <- amb_data  %>%
  setNames(tolower(names(.))) %>%
  filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% 
  #filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN", "PRESERVED_SPECIMEN")) %>%
  filter(basisofrecord %in% c("HUMAN_OBSERVATION")) %>%
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% 
  filter(decimallatitude != 0 | decimallongitude != 0) %>%
  #filter(is.na(as.factor(nwords(locality)))) %>%
  cc_cen(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove country centroids within 2km 
  cc_cap(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove capitals centroids within 2km
  cc_inst(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove zoo and herbaria within 2km 
  cc_sea(lon = "decimallongitude", lat = "decimallatitude") %>% # remove from ocean
  distinct(decimallongitude, decimallatitude, specieskey,datasetkey, .keep_all = TRUE) %>%
  glimpse() 

amb_clean$species <- stri_replace_all_fixed(amb_clean$species, " ", "_")

cal_clean <- cal_data  %>%
  setNames(tolower(names(.))) %>%
  filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% 
  #filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN", "PRESERVED_SPECIMEN")) %>%
  filter(basisofrecord %in% c("HUMAN_OBSERVATION")) %>%
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% 
  filter(decimallatitude != 0 | decimallongitude != 0) %>%
  #filter(is.na(as.factor(nwords(locality)))) %>%
  cc_cen(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove country centroids within 2km 
  cc_cap(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove capitals centroids within 2km
  cc_inst(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove zoo and herbaria within 2km 
  cc_sea(lon = "decimallongitude", lat = "decimallatitude") %>% # remove from ocean
  distinct(decimallongitude, decimallatitude, specieskey,datasetkey, .keep_all = TRUE) %>%
  glimpse() 

cal_clean$species <- stri_replace_all_fixed(cal_clean$species, " ", "_")

psi_clean <- psi_data  %>%
  setNames(tolower(names(.))) %>%
  filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% 
  #filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN", "PRESERVED_SPECIMEN")) %>%
  filter(basisofrecord %in% c("HUMAN_OBSERVATION")) %>%
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% 
  filter(decimallatitude != 0 | decimallongitude != 0) %>%
  #filter(is.na(as.factor(nwords(locality)))) %>%
  cc_cen(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove country centroids within 2km 
  cc_cap(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove capitals centroids within 2km
  cc_inst(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove zoo and herbaria within 2km 
  cc_sea(lon = "decimallongitude", lat = "decimallatitude") %>% # remove from ocean
  distinct(decimallongitude, decimallatitude, specieskey,datasetkey, .keep_all = TRUE) %>%
  glimpse() 

psi_clean$species <- stri_replace_all_fixed(psi_clean$species, " ", "_")

mat_clean <- mat_data  %>%
  setNames(tolower(names(.))) %>%
  filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% 
  #filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN", "PRESERVED_SPECIMEN")) %>%
  filter(basisofrecord %in% c("HUMAN_OBSERVATION")) %>%
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% 
  filter(decimallatitude != 0 | decimallongitude != 0) %>%
  #filter(is.na(as.factor(nwords(locality)))) %>%
  cc_cen(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove country centroids within 2km 
  cc_cap(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove capitals centroids within 2km
  cc_inst(buffer = 2000, lon = "decimallongitude", lat = "decimallatitude") %>% # remove zoo and herbaria within 2km 
  cc_sea(lon = "decimallongitude", lat = "decimallatitude") %>% # remove from ocean
  distinct(decimallongitude, decimallatitude, specieskey,datasetkey, .keep_all = TRUE) %>%
  glimpse() 

mat_clean$species <- stri_replace_all_fixed(mat_clean$species, " ", "_")

## creating a SP object to subset env data
cyc_clean_sp <- cyc_clean
gin_clean_sp <- gin_clean
gne_clean_sp <- gne_clean
pin_clean_sp <- pin_clean

amb_clean_sp <- amb_clean
cal_clean_sp <- cal_clean
psi_clean_sp <- psi_clean
mat_clean_sp <- mat_clean

coordinates(cyc_clean_sp) <- ~decimallongitude+decimallatitude
coordinates(gin_clean_sp) <- ~decimallongitude+decimallatitude
coordinates(gne_clean_sp) <- ~decimallongitude+decimallatitude
coordinates(pin_clean_sp) <- ~decimallongitude+decimallatitude

coordinates(amb_clean_sp) <- ~decimallongitude+decimallatitude
coordinates(cal_clean_sp) <- ~decimallongitude+decimallatitude
coordinates(psi_clean_sp) <- ~decimallongitude+decimallatitude
coordinates(mat_clean_sp) <- ~decimallongitude+decimallatitude

crs(cyc_clean_sp) <- crs(envs)
crs(gin_clean_sp) <- crs(envs)
crs(gne_clean_sp) <- crs(envs)
crs(pin_clean_sp) <- crs(envs)

crs(amb_clean_sp) <- crs(envs)
crs(cal_clean_sp) <- crs(envs)
crs(psi_clean_sp) <- crs(envs)
crs(mat_clean_sp) <- crs(envs)

## extracting env data
env_cyc <- raster::extract(envs, cyc_clean_sp, fun = mean, na.rm = T, sp = T)
env_gin <- raster::extract(envs, gin_clean_sp, fun = mean, na.rm = T, sp = T)
env_gne <- raster::extract(envs, gne_clean_sp, fun = mean, na.rm = T, sp = T)
env_pin <- raster::extract(envs, pin_clean_sp, fun = mean, na.rm = T, sp = T)

env_amb <- raster::extract(envs, amb_clean_sp, fun = mean, na.rm = T, sp = T)
env_cal <- raster::extract(envs, cal_clean_sp, fun = mean, na.rm = T, sp = T)
env_psi <- raster::extract(envs, psi_clean_sp, fun = mean, na.rm = T, sp = T)
env_mat <- raster::extract(envs, mat_clean_sp, fun = mean, na.rm = T, sp = T)

env_cyc_df <- data.frame(env_cyc)
env_gin_df <- data.frame(env_gin)
env_gne_df <- data.frame(env_gne)
env_pin_df <- data.frame(env_pin)

env_amb_df <- data.frame(env_amb)
env_cal_df <- data.frame(env_cal)
env_psi_df <- data.frame(env_psi)
env_mat_df <- data.frame(env_mat)

## cleaning the data from extra stuff
env_cyc_df_clean <- env_cyc_df %>% 
  dplyr::select(species, decimallongitude, decimallatitude, bio1:tmin9)

env_gin_df_clean <- env_gin_df %>% 
  dplyr::select(species, decimallongitude, decimallatitude, bio1:tmin9)

env_gne_df_clean <- env_gne_df %>% 
  dplyr::select(species, decimallongitude, decimallatitude, bio1:tmin9)

env_pin_df_clean <- env_pin_df %>% 
  dplyr::select(species, decimallongitude, decimallatitude, bio1:tmin9)

#write.csv(env_cyc_df_clean, "data/Environmental_Data_Cyc.csv")
#write.csv(env_gin_df_clean, "data/Environmental_Data_Gin.csv")
#write.csv(env_gne_df_clean, "data/Environmental_Data_Gne.csv")
#write.csv(env_pin_df_clean, "data/Environmental_Data_Pin.csv")

env_amb_df_clean <- env_amb_df %>% 
  dplyr::select(species, decimallongitude, decimallatitude, bio1:tmin9)

env_cal_df_clean <- env_cal_df %>% 
  dplyr::select(species, decimallongitude, decimallatitude, bio1:tmin9)

env_psi_df_clean <- env_psi_df %>% 
  dplyr::select(species, decimallongitude, decimallatitude, bio1:tmin9)

env_mat_df_clean <- env_mat_df %>% 
  dplyr::select(species, decimallongitude, decimallatitude, bio1:tmin9)

#write.csv(env_amb_df_clean, "data/Environmental_Data_Amb.csv")
#write.csv(env_cal_df_clean, "data/Environmental_Data_Cal.csv")
#write.csv(env_psi_df_clean, "data/Environmental_Data_Psi.csv")
#write.csv(env_mat_df_clean, "data/Environmental_Data_Mat.csv")

env_cyc_df_clean2 <- read.csv("data/Environmental_Data_Cyc.csv",
                              stringsAsFactors = TRUE, header = TRUE) %>%
  dplyr::select(-X) %>% na.omit() %>% droplevels()

env_gin_df_clean2 <- read.csv("data/Environmental_Data_Gin.csv",
                              stringsAsFactors = TRUE, header = TRUE) %>%
  dplyr::select(-X) %>% na.omit() %>% droplevels()

env_gne_df_clean2 <- read.csv("data/Environmental_Data_Gne.csv",
                              stringsAsFactors = TRUE, header = TRUE) %>%
  dplyr::select(-X) %>% na.omit() %>% droplevels()

env_pin_df_clean2 <- read.csv("data/Environmental_Data_Pin.csv",
                              stringsAsFactors = TRUE, header = TRUE) %>%
  dplyr::select(-X) %>% na.omit() %>% droplevels()

env_amb_df_clean2 <- read.csv("data/Environmental_Data_Amb.csv",
                              stringsAsFactors = TRUE, header = TRUE) %>%
  dplyr::select(-X) %>% na.omit() %>% droplevels()

env_cal_df_clean2 <- read.csv("data/Environmental_Data_Cal.csv",
                              stringsAsFactors = TRUE, header = TRUE) %>%
  dplyr::select(-X) %>% na.omit() %>% droplevels()

env_psi_df_clean2 <- read.csv("data/Environmental_Data_Psi.csv",
                              stringsAsFactors = TRUE, header = TRUE) %>%
  dplyr::select(-X) %>% na.omit() %>% droplevels()

env_mat_df_clean2 <- read.csv("data/Environmental_Data_Mat.csv",
                              stringsAsFactors = TRUE, header = TRUE) %>%
  dplyr::select(-X) %>% na.omit() %>% droplevels()

env_df_clean2 <- rbind(env_cyc_df_clean2, env_gin_df_clean2, env_gne_df_clean2,
                       env_pin_df_clean2, env_amb_df_clean2, env_cal_df_clean2, 
                       env_psi_df_clean2, env_mat_df_clean2)

spp_tree <- env_df_clean2$species %in% intersect(env_df_clean2$species,
                                                 tree$tip.label)
env_df_clean3 <- env_df_clean2[spp_tree,] %>% na.omit() %>% droplevels()

## repeating the analyses at genera and family level
gen_tax <- tax$genus
names(gen_tax) <- tax$scientificName

fam_tax <- tax$family
names(fam_tax) <- tax$scientificName

env_df_gen <- env_df_clean3 %>% 
  mutate(genus = ifelse(species %in% names(gen_tax), 
                        gen_tax[as.character(species)],
                        NA)) %>% 
  na.omit() %>% 
  droplevels()

env_df_gen <- env_df_gen[env_df_gen$genus %in% c(names(table(env_df_gen$genus)[table(env_df_gen$genus)>=10])), ]
#save(env_df_gen, file = "data/env_df_gen.RData")

env_df_fam <- env_df_clean3 %>% 
  mutate(family = ifelse(species %in% names(fam_tax), 
                         fam_tax[as.character(species)],
                         NA)) %>% 
  na.omit() %>% 
  droplevels()
#save(env_df_fam, file = "data/env_df_fam.RData")

#load(file = "data/env_df_gen.RData")
#load(file = "data/env_df_fam.RData")

tr_gen <- tree
for (i in 1:length(tr_gen$tip.label)) {
  tr_gen$tip.label[i] <- gen_tax[tr_gen$tip.label[i]]
}
gen <- unique(env_df_gen$genus)
tr_gen <- keep.tip(tr_gen, gen)

tr_fam <- tree
for (i in 1:length(tr_fam$tip.label)) {
  tr_fam$tip.label[i] <- fam_tax[tr_fam$tip.label[i]]
}
fam <- unique(env_df_fam$family)
tr_fam <- keep.tip(tr_fam, fam)

#write.tree(tr_gen, "data/tr_gen.tre")
#write.tree(tr_fam, "data/tr_fam.tre")

## figure 1b

map <- ne_countries(scale = "medium", type = "map_units", returnclass = "sf")

fam_sf <- sf::st_as_sf(env_df_fam %>% filter(!family %in% c("Amborellaceae",
                                                            "Calycanthaceae",
                                                            "Matoniaceae",
                                                            "Psilotaceae")) %>%
                         arrange(family), 
                       coords = c("decimallongitude", "decimallatitude"),
                       crs = "+proj=longlat +datum=WGS84")

col_fam <- carto_pal(12, "Fall")
names(col_fam) <- names(table(fam_sf$family)[order(table(fam_sf$family), decreasing = T)])

ggplot_maps <- ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  geom_sf(data = fam_sf, aes(colour = family), size = 1) + 
  scale_colour_manual(values = col_fam, name = "Family",
                      labels = c("Pinaceae" = "Pinaceae (n = 8)",
                                 "Cupressaceae" = "Cupressaceae (n = 23)",
                                 "Taxaceae" = "Taxaceae (n = 3)",
                                 "Podocarpaceae" = "Podocarpaceae (n = 14)",
                                 "Ginkgoaceae" = "Ginkgoaceae (n = 1)",
                                 "Ephedraceae" = "Ephedraceae (n = 1)",
                                 "Zamiaceae" = "Zamiaceae (n = 9)",
                                 "Araucariaceae" = "Araucariaceae (n = 2)",
                                 "Cycadaceae" = "Cycadaceae (n = 1)",
                                 "Welwitschiaceae" = "Welwitschiaceae (n = 1)",
                                 "Gnetaceae" = "Gnetaceae (n = 1)",
                                 "Sciadopityaceae" = "Sciadopityaceae (n = 1)")) + 
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 10)) +
  guides(col = guide_legend(title.position = "top", title.hjust = 0.5)) 

pdf("figures/Fig1B.pdf")
print(ggplot_maps)  
dev.off()  

######## HPC cluster ########

env_df_gen_abi <- env_df_gen %>%
  filter(genus == "Abies")

env_df_gen_jun <- env_df_gen %>%
  filter(genus == "Juniperus")

env_df_gen_lar <- env_df_gen %>%
  filter(genus == "Larix")

env_df_gen_pic <- env_df_gen %>%
  filter(genus == "Picea")

env_df_gen_pin <- env_df_gen %>%
  filter(genus == "Pinus")

env_df_gen_tax <- env_df_gen %>%
  filter(genus == "Taxus")

env_df_gen <- env_df_gen %>%
  filter(!genus %in% c("Abies", "Juniperus", "Larix", "Picea", "Pinus", "Taxus",
                       "Amborella", "Calycanthus", "Psilotum", "Tmesipteris",
                       "Matonia"))

env_df_gen_abi_50 <- env_df_gen_abi[sample(1:nrow(env_df_gen_abi), 50000), ]
env_df_gen_jun_50 <- env_df_gen_jun[sample(1:nrow(env_df_gen_jun), 50000), ]
env_df_gen_lar_50 <- env_df_gen_lar[sample(1:nrow(env_df_gen_lar), 50000), ]
env_df_gen_pic_50 <- env_df_gen_pic[sample(1:nrow(env_df_gen_pic), 50000), ]
env_df_gen_pin_50 <- env_df_gen_pin[sample(1:nrow(env_df_gen_pin), 50000), ]
env_df_gen_tax_50 <- env_df_gen_tax[sample(1:nrow(env_df_gen_tax), 50000), ]

env_df_gen_50 <- rbind(env_df_gen, env_df_gen_abi_50, env_df_gen_jun_50,
                       env_df_gen_lar_50, env_df_gen_pic_50, env_df_gen_pin_50,
                       env_df_gen_tax_50)

env_df_fam_cup <- env_df_fam %>%
  filter(family == "Cupressaceae")

env_df_fam_pin <- env_df_fam %>%
  filter(family == "Pinaceae")

env_df_fam <- env_df_fam %>%
  filter(!family %in% c("Cupressaceae", "Pinaceae",
                        "Amborellaceae", "Calycanthaceae", "Psilotaceae", 
                        "Matoniaceae"))

env_df_fam_cup_50 <- env_df_fam_cup[sample(1:nrow(env_df_fam_cup), 50000), ]
env_df_fam_pin_50 <- env_df_fam_pin[sample(1:nrow(env_df_fam_pin), 50000), ]

env_df_fam_50 <- rbind(env_df_fam, env_df_fam_cup_50, env_df_fam_pin_50)

calc_pca <- function(data) {
  res <- prcomp(data[, 3:33])
  res2 <- as.data.frame(cbind(as.numeric(data$decimallongitude),
                              as.numeric(data$decimallatitude),
                              as.numeric(res$x[, 1]),
                              as.numeric(res$x[, 2]),
                              as.numeric(res$x[, 3])))
  colnames(res2) <- c("lon", "lat", "PC1", "PC2", "PC3")
  res2
}

env_df_gen_pca_50 <- env_df_gen_50 %>%
  nest(data = -genus) %>%
  mutate(pca = purrr::map(data, function(.x) calc_pca(.x)))

env_df_fam_pca_50 <- env_df_fam_50 %>%
  nest(data = -family) %>% 
  mutate(pca = purrr::map(data, function(.x) calc_pca(.x)))

## calculating means of the PC axes
calc_mean <- function(data) {
  mean1 <- mean(data[, 3])
  mean2 <- mean(data[, 4])
  mean3 <- mean(data[, 5])
  
  centroid <- data.frame(PC1 = mean1, PC2 = mean2, PC3 = mean3)
  names(centroid) <- c("PC1", "PC2", "PC3")
  centroid
}

env_df_gen_pca_50_cent <- env_df_gen_pca_50 %>%
  mutate(centroid = purrr::map(pca, function(.x) calc_mean(.x)))

mean_all <- numeric()
for (i in 1:3) {
  numb <- NULL
  for (j in 1:nrow(env_df_gen_pca_50)) {
    numb <- c(numb, env_df_gen_pca_50$pca[[j]][, 2+i])
  }
  mean_all[i] <- mean(numb)
}

calc_dist <- function(data) {
  distP <- dist(rbind(mean_all, data), method = "euclidean")[1]
  distP
}

env_df_gen_dist <- env_df_gen_pca_50_cent %>%
  mutate(dist_cent = purrr::map(centroid, function(.x) calc_dist(.x)))

save(env_df_gen_dist, file = "data/env_df_gen_dist.RData")

## calculating means of the PC axes
calc_mean <- function(data) {
  mean1 <- mean(data[, 3])
  mean2 <- mean(data[, 4])
  mean3 <- mean(data[, 5])
  
  centroid <- data.frame(PC1 = mean1, PC2 = mean2, PC3 = mean3)
  names(centroid) <- c("PC1", "PC2", "PC3")
  centroid
}

env_df_fam_pca_50_cent <- env_df_fam_pca_50 %>%
  mutate(centroid = purrr::map(pca, function(.x) calc_mean(.x)))

mean_all_f <- numeric()
for (i in 1:3) {
  numb <- NULL
  for (j in 1:nrow(env_df_fam_pca_50)) {
    numb <- c(numb, env_df_fam_pca_50$pca[[j]][, 2+i])
  }
  mean_all_f[i] <- mean(numb)
}

calc_dist <- function(data) {
  distP <- dist(rbind(mean_all_f, data), method = "euclidean")[1]
  distP
}

env_df_fam_dist <- env_df_fam_pca_50_cent %>%
  mutate(dist_cent = purrr::map(centroid, function(.x) calc_dist(.x)))

save(env_df_fam_dist, file = "data/env_df_fam_dist.RData")

## calculating hypervolumes

env_df_gen_hypervolume_50 <- env_df_gen_pca_50 %>%
    mutate(hypervol = purrr::map(pca, function(.x) {
    print("running")
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))
#save(env_df_gen_hypervolume_50, file = "data/cluster/gen/env_df_gen_hypervolume_50.RData")

env_df_fam_hypervolume_50 <- env_df_fam_pca_50 %>%
  mutate(hypervol = purrr::map(pca, function(.x) {
    print("running")
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))
#save(env_df_fam_hypervolume_50, file = "data/cluster/fam/env_df_fam_hypervolume_50.RData")

######## END CLUSTER HPC ########

## outgroups
env_df_gen_out <- env_df_gen %>%
  filter(genus %in% c("Amborella", "Calycanthus", "Psilotum", "Matonia", "Tmesipteris")) %>%
  dplyr::select(-species)

env_df_fam_out <- env_df_fam %>%
  filter(family %in% c("Amborellaceae", "Calycanthaceae", "Psilotaceae", "Matoniaceae")) %>%
  dplyr::select(-species)

env_df_gen_out_pca <- env_df_gen_out %>%
  nest(data = -genus) %>%
  mutate(pca = purrr::map(data, function(.x) calc_pca(.x)))

env_df_fam_out_pca <- env_df_fam_out %>%
  nest(data = -family) %>% 
  mutate(pca = purrr::map(data, function(.x) calc_pca(.x)))

## calculating hypervolumes

env_df_gen_hypervolume_out <- env_df_gen_out_pca %>%
  mutate(hypervol = purrr::map(pca, function(.x) {
    print("running")
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))
#save(env_df_gen_hypervolume_out, file = "data/cluster/gen/env_df_gen_hypervolume_out.RData")

env_df_fam_hypervolume_out <- env_df_fam_out_pca %>%
  mutate(hypervol = purrr::map(pca, function(.x) {
    print("running")
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))
#save(env_df_fam_hypervolume_out, file = "data/cluster/fam/env_df_fam_hypervolume_out.RData")

## figure S1

load(file = "data/cluster/gen/env_df_gen_hypervolume_50.RData")
load(file = "data/cluster/fam/env_df_fam_hypervolume_50.RData")

abi_gen_sf <- sf::st_as_sf(env_df_gen_abi,
                           coords = c("decimallongitude", "decimallatitude"),
                           crs = "+proj=longlat +datum=WGS84")
abi_gen_50_sf <- sf::st_as_sf(env_df_gen_hypervolume_50$data[env_df_gen_hypervolume_50$genus == "Abies"][[1]],
                              coords = c("decimallongitude", "decimallatitude"),
                              crs = "+proj=longlat +datum=WGS84")

jun_gen_sf <- sf::st_as_sf(env_df_gen_jun,
                           coords = c("decimallongitude", "decimallatitude"),
                           crs = "+proj=longlat +datum=WGS84")
jun_gen_50_sf <- sf::st_as_sf(env_df_gen_hypervolume_50$data[env_df_gen_hypervolume_50$genus == "Juniperus"][[1]],
                              coords = c("decimallongitude", "decimallatitude"),
                              crs = "+proj=longlat +datum=WGS84")

lar_gen_sf <- sf::st_as_sf(env_df_gen_lar,
                           coords = c("decimallongitude", "decimallatitude"),
                           crs = "+proj=longlat +datum=WGS84")
lar_gen_50_sf <- sf::st_as_sf(env_df_gen_hypervolume_50$data[env_df_gen_hypervolume_50$genus == "Larix"][[1]],
                              coords = c("decimallongitude", "decimallatitude"),
                              crs = "+proj=longlat +datum=WGS84")

pic_gen_sf <- sf::st_as_sf(env_df_gen_pic,
                           coords = c("decimallongitude", "decimallatitude"),
                           crs = "+proj=longlat +datum=WGS84")
pic_gen_50_sf <- sf::st_as_sf(env_df_gen_hypervolume_50$data[env_df_gen_hypervolume_50$genus == "Picea"][[1]],
                              coords = c("decimallongitude", "decimallatitude"),
                              crs = "+proj=longlat +datum=WGS84")

pin_gen_sf <- sf::st_as_sf(env_df_gen_pin,
                           coords = c("decimallongitude", "decimallatitude"),
                           crs = "+proj=longlat +datum=WGS84")
pin_gen_50_sf <- sf::st_as_sf(env_df_gen_hypervolume_50$data[env_df_gen_hypervolume_50$genus == "Pinus"][[1]],
                              coords = c("decimallongitude", "decimallatitude"),
                              crs = "+proj=longlat +datum=WGS84")

tax_gen_sf <- sf::st_as_sf(env_df_gen_tax,
                           coords = c("decimallongitude", "decimallatitude"),
                           crs = "+proj=longlat +datum=WGS84")
tax_gen_50_sf <- sf::st_as_sf(env_df_gen_hypervolume_50$data[env_df_gen_hypervolume_50$genus == "Taxus"][[1]],
                              coords = c("decimallongitude", "decimallatitude"),
                              crs = "+proj=longlat +datum=WGS84")

cup_fam_sf <- sf::st_as_sf(env_df_fam_cup,
                           coords = c("decimallongitude", "decimallatitude"),
                           crs = "+proj=longlat +datum=WGS84")
cup_fam_50_sf <- sf::st_as_sf(env_df_fam_hypervolume_50$data[env_df_fam_hypervolume_50$family == "Cupressaceae"][[1]],
                           coords = c("decimallongitude", "decimallatitude"),
                           crs = "+proj=longlat +datum=WGS84")

pin_fam_sf <- sf::st_as_sf(env_df_fam_pin,
                           coords = c("decimallongitude", "decimallatitude"),
                           crs = "+proj=longlat +datum=WGS84")
pin_fam_50_sf <- sf::st_as_sf(env_df_fam_hypervolume_50$data[env_df_fam_hypervolume_50$family == "Pinaceae"][[1]],
                              coords = c("decimallongitude", "decimallatitude"),
                              crs = "+proj=longlat +datum=WGS84")

abi_gen_maps <- ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  geom_sf(data = abi_gen_sf, color = col_fam[1],
          size = 1) + 
  theme_void()

abi_gen_50_maps <- ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  geom_sf(data = abi_gen_50_sf, color = col_fam[12],
          size = 1) + 
  theme_void() +
  labs(title = "Abies") + 
  theme(title = element_text(face = "bold", size = 7),
        plot.title = element_text(hjust = 1))

jun_gen_maps <- ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  geom_sf(data = jun_gen_sf, color = col_fam[1],
          size = 1) + 
  theme_void() 

jun_gen_50_maps <- ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  geom_sf(data = jun_gen_50_sf, color = col_fam[12],
          size = 1) + 
  theme_void() +
  labs(title = "Juniperus") + 
  theme(title = element_text(face = "bold", size = 7),
        plot.title = element_text(hjust = 1))

lar_gen_maps <- ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  geom_sf(data = lar_gen_sf, color = col_fam[1],
          size = 1) + 
  theme_void() +
  labs(title = "A") + 
  theme(title = element_text(face = "bold", size = 10))

lar_gen_50_maps <- ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  geom_sf(data = lar_gen_50_sf, color = col_fam[12],
          size = 1) + 
  theme_void() +
  labs(title = "Larix") + 
  theme(title = element_text(face = "bold", size = 7),
        plot.title = element_text(hjust = 1))

pic_gen_maps <- ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  geom_sf(data = pic_gen_sf, color = col_fam[1],
          size = 1) + 
  theme_void() 

pic_gen_50_maps <- ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  geom_sf(data = pic_gen_50_sf, color = col_fam[12],
          size = 1) + 
  theme_void() +
  labs(title = "Picea") + 
  theme(title = element_text(face = "bold", size = 7),
        plot.title = element_text(hjust = 1))

pin_gen_maps <- ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  geom_sf(data = pin_gen_sf, color = col_fam[1],
          size = 1) + 
  theme_void() 

pin_gen_50_maps <- ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  geom_sf(data = pin_gen_50_sf, color = col_fam[12],
          size = 1) + 
  theme_void() +
  labs(title = "Pinus") + 
  theme(title = element_text(face = "bold", size = 7),
        plot.title = element_text(hjust = 1))

tax_gen_maps <- ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  geom_sf(data = tax_gen_sf, color = col_fam[1],
          size = 1) + 
  theme_void() 

tax_gen_50_maps <- ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  geom_sf(data = tax_gen_50_sf, color = col_fam[12],
          size = 1) + 
  theme_void() +
  labs(title = "Taxus") + 
  theme(title = element_text(face = "bold", size = 7),
        plot.title = element_text(hjust = 1))

cup_fam_maps <- ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  geom_sf(data = cup_fam_sf, color = col_fam[1],
          size = 1) + 
  theme_void() +
  labs(title = "B") + 
  theme(title = element_text(face = "bold", size = 10))

cup_fam_50_maps <- ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  geom_sf(data = cup_fam_50_sf, color = col_fam[12],
          size = 1) + 
  theme_void() +
  labs(title = "Cupressaceae") + 
  theme(title = element_text(face = "bold", size = 7),
        plot.title = element_text(hjust = 1))

pin_fam_maps <- ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  geom_sf(data = pin_fam_sf, color = col_fam[1],
          size = 1) + 
  theme_void() 

pin_fam_50_maps <- ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  geom_sf(data = pin_fam_50_sf, color = col_fam[12],
          size = 1) + 
  theme_void() +
  labs(title = "Pinaceae") + 
  theme(title = element_text(face = "bold", size = 7),
        plot.title = element_text(hjust = 1))

all_maps <- list(lar_gen_maps, lar_gen_50_maps,
                 pic_gen_maps, pic_gen_50_maps,
                 tax_gen_maps, tax_gen_50_maps,
                 pin_gen_maps, pin_gen_50_maps,
                 abi_gen_maps, abi_gen_50_maps,
                 jun_gen_maps, jun_gen_50_maps,
                 cup_fam_maps, cup_fam_50_maps, 
                 pin_fam_maps, pin_fam_50_maps)

pdf("figures/FigS1.pdf", height = 15, width = 7)
wrap_plots(all_maps) + plot_layout(ncol = 2)
dev.off()