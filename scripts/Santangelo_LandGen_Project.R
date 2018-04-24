# Load Cleaned microsat loci data
MicroSat_Data <- read.csv("data-clean/Johnson-et-al_8-Cities_MicroSat-Loci", header = TRUE)

# Load required packages
library(adegenet)
library(pegas)
library(poppr)
library(dplyr)

# Create dataset with only Acton and Fergus
MicroSat_Data_sub <- MicroSat_Data %>%
  filter(City == "Fergus" | City == "Acton") 

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
Fergus.genind <- create_genind(MicroSat_Data_sub[MicroSat_Data_sub$City == "Fergus", ] )
Acton.genind <- create_genind(MicroSat_Data_sub[MicroSat_Data_sub$City == "Acton", ] )

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
