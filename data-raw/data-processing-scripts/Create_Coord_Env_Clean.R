# Load packages
library(tidyr)
library(dplyr)

# Load habitat type data
Coord_Env_data <- "data-raw/Coords-Enviro.csv"
Coord_Env_data <- read.csv(Coord_Env_data, header = TRUE)

# Create final data set
Coord_Env_data <- Coord_Env_data %>% 
  group_by(city) %>% # Group
  mutate(site = 1:n()) %>% # Number sites sequentially
  filter(city == "Fergus" | city == "Acton") %>% # Only Fergus and Acton
  slice(rep(1:n(), each = 10)) %>% # Add 10 rows per city (i.e. number of plants)
  ungroup() %>% # Ungroup
  group_by(city, site) %>% # Regroup by site
  mutate(plant = 1:n(), # Number plants sequentially within sites
         city_code = ifelse(city == "Fergus", "Fe", "Ac"), # Add city code
         PlantID = paste(city_code, site, plant, sep = ".")) %>% # Add PlantID
  filter(site != 1) %>%
  ungroup() %>% # Ungroup 
  select(PlantID, everything()) %>% #Move PlantID to beginning
  select(PlantID, lat, long, buildings, impervious) %>%
  as.data.frame()

# Create coordinates dataset
Coord_data <- Coord_Env_data %>%
  select(PlantID, long, lat)

# Create environmental dataset
Env_data <- Coord_Env_data %>%
  select(PlantID, buildings, impervious)

# Write cleanes up dataset to CSV
write.csv(Coord_data, "data-clean/Coord-data.csv", row.names = FALSE)
write.csv(Env_data, "data-clean/Env-data.csv", row.names = FALSE)
