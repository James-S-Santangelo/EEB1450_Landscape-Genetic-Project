# Load required packages
library(tibble)
library(dplyr)
library(tidyr)

# Load Microsat data
MicroSat_Data <- "data-raw/Johnson-et-al_8-Cities_MicroSat-Alleles.csv"
MicroSat_Data <- read.csv(MicroSat_Data, header = TRUE)

# Create dataset with only Acton and Fergus
MicroSat_Data <- MicroSat_Data %>%
  filter(City == "Fergus" | City == "Acton") 

# Look at data to confirm correct import
as.tibble(MicroSat_Data)

# Create final dataset
MicroSat_Data <- MicroSat_Data %>% 
  separate(PopID, sep = "[.]", into = c("city_pop", "plant"), remove = TRUE) %>% #Get plant #
  separate(city_pop, sep = 1, into = c("random", "CodePop"), remove = TRUE) %>% # Get Population code (City ID + Pop)
  select(-random) %>% # Remove random letter that was at beginning of ID string
  separate(CodePop, sep = "(?<=[A-Za-z])(?=[0-9])", into = c("city_code", "pop"), remove = TRUE) %>% # Split Population code into City ID and Pop columns
  mutate(PlantID = paste(city_code, pop, plant, sep = ".")) %>% # Concatenate City ID, Pop # and Plant # to form Plant ID
  select(PlantID, everything()) %>% # Move Plant ID column to start of data frame
  filter(pop != 1) %>%
  select(-City, -city_code, -plant, -pop) %>%
  as.data.frame()

# Look at data to confirm correct structure
as.tibble(MicroSat_Data)

# Combine Alleles into single column, separated by colons
MicroSat_Data <- data.frame(MicroSat_Data[1], 
                            mapply( paste, MicroSat_Data[-1][c(T,F)], 
                                    MicroSat_Data[-1][c(F,T)], sep = ":") )

# Write cleanes up dataset to CSV
write.csv(MicroSat_Data, "data-clean/Johnson-et-al_8-Cities_MicroSat-Loci.csv", row.names = FALSE)
