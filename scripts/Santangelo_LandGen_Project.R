# Load required packages
library(EcoGenetics)

# library(adegenet)
# library(pegas)
# library(poppr)
# library(dplyr)
# library(hierfstat)

# Load in all datasets
MicroSat <- read.csv("data-clean/Johnson-et-al_8-Cities_MicroSat-Loci.csv")
Habitat <- read.csv("data-clean/Johnson-et-al_8-Cities_Habitat-data.csv")
Coord <- read.csv("data-clean/Johnson-et-al_8-Cities_Coord-data.csv")
Enviro <- read.csv("data-clean/Johnson-et-al_8-Cities_Env-data.csv")


add_rownames <- function(data_frame){
  row.names(data_frame) <- data_frame$PlantID
  data_frame <- data_frame %>% 
    select(-PlantID)
  # data_frame <- data_frame[order(row.names(data_frame)),]
  return(data_frame)
}

MicroSat <- add_rownames(MicroSat)
Habitat <- add_rownames(Habitat)
Coord <- add_rownames(Coord)
Enviro <- add_rownames(Enviro)

MicroSat_NoPop <- MicroSat %>%
  select(-pop)
  
# Create ecogenetics object
cities.ecogen <- ecogen(XY = Coord,
                        G.processed = TRUE, order.G = TRUE, type = "codominant",
                        ploidy = 2, sep = ":", ncod = NULL, missing = "0",
                        NA.char = "0", poly.level = NULL, rm.empty.ind = TRUE, order.df = TRUE,
                        set.names = NULL, valid.names = FALSE)
ecoslot.G(cities.ecogen, type = "codominant", 
          ploidy = 2, sep = ":", NA.char = "0") <- MicroSat_NoPop
ecoslot.E(cities.ecogen) <- Enviro
ecoslot.S(cities.ecogen) <- Habitat

# Create genind object for each city
create_genind <- function(data_frame){
  genind <- df2genind(X = data_frame[,c(6:21)], sep=":",
            ncode = NULL,
            ind.names= data_frame$PlantID,
            loc.names=NULL,
            pop = data_frame$pop,
            NA.char="0",
            ploidy=2,
            type="codom",
            strata=NULL,
            hierarchy=NULL)
  
  return(genind)
}

Fergus.genind <- create_genind(MicroSat_Data[MicroSat_Data$City == "Fergus", ] )
Acton.genind <- create_genind(MicroSat_Data[MicroSat_Data$City == "Acton", ] )

Fergus.genind@other$Enviro <- MicroSat_Data[MicroSat_Data$City == "Fergus", ]$Enviro
Acton.genind@other$Enviro <- MicroSat_Data[MicroSat_Data$City == "Acton", ]$Enviro

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

HWE_Loci_Pops_Acton <- data.frame(HWE_Loci_Pops(Acton.genind))
HWE_Loci_Pops_Acton <- data.frame(HWE_Loci_Pops(Fergus.genind))

# Global LD and LD among marker pairs
LD_Acton <- poppr::ia(Acton.genind, sample = 1000)
LD_Fergus <- poppr::ia(Fergus.genind, sample = 1000)
LD_Pair_Acton <- poppr::pair.ia(Acton.genind)
LD_Pair_Fergus <- poppr::pair.ia(Fergus.genind)

fstat(Acton.genind)
pairwise.fst(Acton.genind, pop=NULL, res.type=c("dist","matrix"))
gstat.randtest(Acton.genind,nsim=99)

Fergus_Clusters <- find.clusters(Fergus.genind, max.n.clust=40)
Fergus_dapc <- dapc(Fergus.genind, Fergus_Clusters$grp)
scatter(Fergus_dapc)

# Genind for both Acton and Fergus combined
Act.Ferg.genind <- create_genind(MicroSat_Data_sub)

Act_Ferg_Clusters <- find.clusters(Act.Ferg.genind, max.n.clust=40)
Act_Ferg_dapc <- dapc(Act.Ferg.genind, Act_Ferg_Clusters$grp, grp = "City")
scatter(Act_Ferg_dapc, posi.da="bottomright", bg="white",
        pch=17:22, cstar=0, scree.pca=TRUE,
        posi.pca="bottomleft")
