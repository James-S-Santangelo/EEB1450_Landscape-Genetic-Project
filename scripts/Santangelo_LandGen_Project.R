# Load required packages
library(EcoGenetics)
library(adegenet)
library(pegas)
library(poppr)
library(dplyr)
library(SoDA)
library(hierfstat)
library(gstudio)
library(PopGenReport)
library(mmod)

# Load in all datasets
MicroSat <- read.csv("data-clean/Johnson-et-al_8-Cities_MicroSat-Loci.csv")
Structure <- read.csv("data-clean/Johnson-et-al_8-Cities_Structure-data.csv")
Coord <- read.csv("data-clean/Johnson-et-al_8-Cities_Coord-data.csv")
Enviro <- read.csv("data-clean/Johnson-et-al_8-Cities_Env-data.csv")

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
  
  
  dataFrames <- list(MicroSat = MicroSat, 
                     Structure = Structure, 
                     Coord = Coord, 
                     Enviro = Enviro)
  toString(City)
  subsetter = substr(City, 1, 2)
  
  subsetted_dataFrames <- list()

  for(i in 1:length(dataFrames)){
    rows <- grep(paste0("^", subsetter), rownames(dataFrames[[i]]))
    subsetted <- dataFrames[[i]][rows, ]
    name <- names(dataFrames[i])

    subsetted_dataset_name <- paste0(name, "Sub")
    
    subsetted_dataFrames[[subsetted_dataset_name]] <- subsetted
    
  }
  
  EcoGen.name <- paste(City, "ecogen", sep = ".")
  EcoGen.name <- ecogen(XY = subsetted_dataFrames$CoordSub, 
                        G = subsetted_dataFrames$MicroSatSub,
                        E = subsetted_dataFrames$EnviroSub, 
                        S = subsetted_dataFrames$StructureSub,
                        G.processed = TRUE, order.G = TRUE, type = "codominant",
                        ploidy = 2, sep = ":", ncod = NULL, missing = "NA",
                        NA.char = "0", poly.level = NULL, rm.empty.ind = TRUE, order.df = TRUE,
                        set.names = NULL, valid.names = FALSE)
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

# Add pop to genin objects
Fergus.genind@pop <- Fergus.genind@strata$pop
Acton.genind@pop <- Acton.genind@strata$pop

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

HWE_Loci_Pops_Acton <- data.frame(HWE_Loci_Pops(Acton.genind))
HWE_Loci_Pops_Acton <- data.frame(HWE_Loci_Pops(Fergus.genind))

# Global LD and LD among marker pairs
LD_Acton <- poppr::ia(Acton.genind, sample = 1000)
LD_Fergus <- poppr::ia(Fergus.genind, sample = 1000)
LD_Pair_Acton <- poppr::pair.ia(Acton.genind)
LD_Pair_Fergus <- poppr::pair.ia(Fergus.genind)

fstat(Acton.genind, pop = Acton.genind@strata$pop)
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

GD.pop.PairwiseFst.hierfstat <- hierfstat::pairwise.fst(Fergus.genind,
                                                        pop = NULL, res.type = c("dist"))
GD.pop.propShared <- PopGenReport::pairwise.propShared(Fergus.genind)
GD.pop.Nei <- adegenet::dist.genpop(Fergus.genpop, method=1)
GD.pop.Edwards <- adegenet::dist.genpop(Fergus.genpop, method=2)
GD.pop.Reynolds <- adegenet::dist.genpop(Fergus.genpop, method=3)
GD.pop.Rogers <- adegenet::dist.genpop(Fergus.genpop, method=4)
GD.pop.Provesti <- adegenet::dist.genpop(Fergus.genpop, method=5)
GD.pop.Joost <- mmod::pairwise_D(Fergus.genind, linearized = FALSE)
GD.pop.Hedrick <- mmod::pairwise_Gst_Hedrick(Fergus.genind, linearized = FALSE)
GD.pop.NeiGst <- mmod::pairwise_Gst_Nei(Fergus.genind, linearized = FALSE)

GD.pop <- list(pairwiseFst.hierfstat = GD.pop.PairwiseFst.hierfstat,
               propShared.PopGenReport = 1 - GD.pop.propShared,
               Nei.adegenet = GD.pop.Nei,
               Edwards.adegenet = GD.pop.Edwards,
               Reynolds.adegenet = GD.pop.Reynolds,
               Rogers.adegenet = GD.pop.Rogers,
               Provesti.adegenet = GD.pop.Provesti,
               Joost.mmod = GD.pop.Joost,
               Hedrick.mmod = GD.pop.Hedrick,
               Nei.mmod = GD.pop.NeiGst)
round(cor(sapply(GD.pop, function(ls) as.vector(ls))),2)[,1:2]

Dgeo <- dist(Fergus.genpop@other$xy[,c("X", "Y")])
par(mar=c(4,4,0,0))
Dgen <- GD.pop$pairwiseFst.hierfstat
dens <- MASS::kde2d(Dgeo, Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold","orange","red"))
plot(Dgeo, Dgen, pch=20, cex=0.5,
     xlab="Geographic Distance", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(Dgen ~ Dgeo))
lines(loess.smooth(Dgeo, Dgen), col="red")


corm <- eco.cormantel(M = dist(Acton.ecogen[["A"]]), 
                      XY = Acton.ecogen[["XY"]],  
                      nsim = 1000, 
                      nclass = 4,
                      alternative = "less")
corm
eco.plotCorrelog(corm)

GD.pop$propShared.PopGenReport
