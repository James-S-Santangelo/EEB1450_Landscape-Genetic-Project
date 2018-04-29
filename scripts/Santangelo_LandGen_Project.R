# Load required packages
library(EcoGenetics)
library(adegenet)
library(pegas)
library(poppr)
library(PopGenReport)
library(dplyr)
library(SoDA)
library(hierfstat)
library(vegan)

# Create checkpoint with package versions on date analysis was performed.
# Install packages and dependencies in project root.
# Will return error if R version differs.
# R v.3.4.3 source code can be downloaded from https://cran.rstudio.com/
library(checkpoint)
checkpoint("2018-04-29", project = getwd(),
           checkpointLocation = "./", verbose = TRUE,
           forceInstall = TRUE, forceProject = TRUE)

# Load in all datasets
MicroSat <- read.csv("data-clean/MicroSat-Loci.csv")
Structure <- read.csv("data-clean/Structure-data.csv")
Coord <- read.csv("data-clean/Coord-data.csv")
Enviro <- read.csv("data-clean/Env-data.csv")

# Add rownames and order dataframes
add_rownames <- function(data_frame){
  row.names(data_frame) <- data_frame$PlantID
  data_frame <- data_frame[order(row.names(data_frame)),]
  data_frame <- data_frame %>% 
    select(-PlantID)
  return(data_frame)
}

# Apply function above to dataframes
MicroSat <- add_rownames(MicroSat)
Structure <- add_rownames(Structure)
Coord <- add_rownames(Coord)
Enviro <- add_rownames(Enviro)
  
# Create ecogenetics object
create_ecogen <- function(MicroSat, Structure, Coord, Enviro, City){
  
  # Create list with data frames for ecogen object
  dataFrames <- list(MicroSat = MicroSat, 
                     Structure = Structure, 
                     Coord = Coord, 
                     Enviro = Enviro)
  
  # First two letters of city name will be string used to subset dataframes
  toString(City)
  subsetter = substr(City, 1, 2)
  subsetted_dataFrames <- list() # Empty list to append subsetted dataframes

  # Loop through dataframes and subset by city. Append to list
  for(i in 1:length(dataFrames)){
    rows <- grep(paste0("^", subsetter), rownames(dataFrames[[i]]))
    subsetted <- dataFrames[[i]][rows, ]
    name <- names(dataFrames[i])
    subsetted_dataset_name <- paste0(name, "Sub")
    subsetted_dataFrames[[subsetted_dataset_name]] <- subsetted
  }
  # Create ecogen object using subsetted dataframes
  EcoGen.name <- paste(City, "ecogen", sep = ".")
  EcoGen.name <- ecogen(XY = subsetted_dataFrames$CoordSub, 
                        G = subsetted_dataFrames$MicroSatSub,
                        E = subsetted_dataFrames$EnviroSub, 
                        S = subsetted_dataFrames$StructureSub,
                        G.processed = TRUE, order.G = TRUE, type = "codominant",
                        ploidy = 2, sep = ":", ncod = NULL, missing = "NA",
                        NA.char = "0", poly.level = NULL, rm.empty.ind = TRUE, order.df = TRUE,
                        set.names = NULL, valid.names = FALSE)
  
  # Change coordinates to X, Y in kilometers
  EcoGen.name@XY <- as.data.frame(SoDA::geoXY(EcoGen.name@XY[,'lat'],
                                              EcoGen.name@XY[,'long'],
                                              unit = 1000))
  return(EcoGen.name)
}

# Create Ecogen objects for each city
Fergus.ecogen <- create_ecogen(MicroSat, Structure, Coord, Enviro, "Fergus")
Acton.ecogen <- create_ecogen(MicroSat, Structure, Coord, Enviro, "Acton")

# Create genind object for each city
Fergus.genind <- ecogen2genind(Fergus.ecogen)
Acton.genind <- ecogen2genind(Acton.ecogen)

# Add pop to genind objects
Fergus.genind@pop <- Fergus.genind@strata$pop
Acton.genind@pop <- Acton.genind@strata$pop

# Create Genpop objects. Average X, Y coordinates for each population.
Fergus.genpop <- adegenet::genind2genpop(Fergus.genind, process.other = TRUE)
Acton.genpop <- adegenet::genind2genpop(Acton.genind, process.other = TRUE)

# Inspect and summarize each genind objects
summary(Fergus.genind)
summary(Acton.genind)

# Deviation from HWE for each locus, across populations
round(pegas::hw.test(Acton.genind, B = 1000), digits = 3)
round(pegas::hw.test(Fergus.genind, B = 1000), digits = 3)

# Deviations from HWE for each locus and each population individually
HWE_Loci_Pops <- function(genind_object){
  HWE.test <- data.frame(sapply(seppop(genind_object), 
                                function(ls) pegas::hw.test(ls, B=1000)[,4]))
  HWE.test.MC <- t(data.matrix(HWE.test))
  {cat("Monte Carlo (p-values):", "\n")
    round(HWE.test,3)}
}

HWE_Loci_Pops_Acton <- HWE_Loci_Pops(Acton.genind)
HWE_Loci_Pops_Fergus <- data.frame(HWE_Loci_Pops(Fergus.genind))

# Global LD and LD among marker pairs
LD_Acton <- poppr::ia(Acton.genind, sample = 10)
LD_Fergus <- poppr::ia(Fergus.genind, sample = 1000)
LD_Pair_Acton <- poppr::pair.ia(Acton.genind)
LD_Pair_Fergus <- poppr::pair.ia(Fergus.genind)

create_IBD_plot <- function(genind_object, genpop_object){
  
  # Calculate proportion of shared alleles among populations
  # Returns pairwise Fst matrix
  Dgen <- hierfstat::pairwise.fst(genind_object,
                                  pop = NULL, res.type = c("dist"))
  
  # Calculate geographic distance matrix
  Dgeo <- dist(genpop_object@other$xy[,c("X", "Y")])
  par(mar=c(4,4,0,0))
  
  # Plot Genetic distance against genetic distance 
  dens <- MASS::kde2d(Dgeo, Dgen, n=300) # Estimate 2 dimensional kernal.
  myPal <- colorRampPalette(c("white","blue","gold","orange","red")) # Set color palette for kernal
  plot(Dgeo, Dgen, pch=20, cex=0.5,
       xlab="Geographic Distance", ylab="Genetic Distance") # Create plot
  image(dens, col=transp(myPal(300), 0.7), add=TRUE) # Add kernal density to plot
  abline(lm(Dgen ~ Dgeo)) # Add linear model to plot
  lines(loess.smooth(Dgeo, Dgen), col="red") # Add loess smoother
}

mantel_tests <- function(ecogen_object, genind_object, genpop_object, nclasses){
  Dgen <- hierfstat::pairwise.fst(genind_object,
                                  pop = NULL, res.type = c("dist"))
  Dgeo <- dist(genpop_object@other$xy[,c("X", "Y")])
  
  print(ade4::mantel.randtest(Dgen,Dgeo))
  
  print(eco.cormantel(M = dist(ecogen_object[["A"]]), 
                XY = ecogen_object[["XY"]],  
                nsim = 1000, 
                nclass = nclasses,
                alternative = "less"))
  
}

# Isolation by distance plots
Fergus_IBD_plot <- create_IBD_plot(Fergus.genind, Fergus.genpop)
Acton_IBD_plot <- create_IBD_plot(Acton.genind, Acton.genpop)

# Mantel test and mantel correllelograms
mantel_tests(Fergus.ecogen,Fergus.genind, Fergus.genpop, nclasses = 8)
mantel_tests(Acton.ecogen, Acton.genind, Acton.genpop, nclasses = 5)

# Function to perform RDA
perform_rda <- function(Structure, genpop_object, City){
  
  # Subset Structure dataframe for City and return one observation per population
  toString(City)
  subsetter = substr(City, 1, 2)
  rows <- grep(paste0("^", subsetter), rownames(Structure))
  city.structure <- Structure[rows, ]
  city.habitat <- city.structure %>%
    group_by(pop) %>%
    slice(1)
  
  # Define datasets for RDA
  MicroSats <- makefreq(genpop_object) # Alleles
  Habitats <- city.habitat # Habitat (i.e. urban or rural)
  Geo_coords <- genpop_object@other$xy # Geographic coordinates of population
  
  # Perform RDA for effect of urbanization, conditioned on distance. 
  rda <- vegan::rda(MicroSats ~ Habitat + Condition(Geo_coords$X + Geo_coords$Y), data = Habitats)
  Pvals <- anova.cca(rda, permutations = 10000) # Permute RDA to get P-values
  
  # Return habitat data, rda, and Pvals as list.
  lst <- list(Habitats, rda, Pvals)
  return(lst)
}

# Perform RDAs for each city. 
Acton.RDA <- perform_rda(Structure, Acton.genpop, "Acton")
Fergus.RDA <- perform_rda(Structure, Fergus.genpop, "Fergus")

plot_rda <- function(RDA_object){
  
  # Extract habitat and RDA data from list returned from RDA function
  Habitat_data <- RDA_object[[1]]
  RDA_model <- RDA_object[[2]]
  
  # Color vector for plotting
  colvec <- c("red2", "green4")
  
  # Create plot
  plot(RDA_model, type = "n")
  
  # Add points, colored by habitat
  with(Habitat_data, points(RDA_model, display = "sites", col = colvec[Habitat],
                               scaling = 3, pch = 21, bg = colvec[Habitat]))
  
  # Add legend
  with(Habitat_data, legend("topright", legend = levels(Habitat), bty = "n",
                               col = colvec, pch = 21, pt.bg = colvec))
  
  ordiellipse(RDA_model, Habitat_data$Habitat, col = colvec)
}

# Generate RDA biplots. 
plot_rda(Fergus.RDA)
plot_rda(Acton.RDA)
