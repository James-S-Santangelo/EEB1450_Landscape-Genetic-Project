# Load Cleaned microsat loci data
MicroSat_Data <- read.csv("data-clean/Johnson-et-al_8-Cities_MicroSat-Loci", header = TRUE)

# Load required packages
library(adegenet)
library(pegas)

# Create dataframe for each city. Stored in list
Cities <- split(MicroSat_Data, MicroSat_Data$City)

# Create genind object for each city
for (i in 1:length(Cities)){
  name = sprintf("%s.genind", names(Cities[i]))
  data_frame = Cities[i][[1]]
  assign(name, df2genind(X = data_frame[,c(6:21)], sep=":",
                             ncode = NULL,
                             ind.names= data_frame$PlantID,
                             loc.names=NULL,
                             pop = data_frame$pop,
                             NA.char="0",
                             ploidy=2,
                             type="codom",
                             strata=NULL,
                             hierarchy=NULL))
}

# Inspect and summarize each genind objects


# Deviation from HWE
round(pegas::hw.test(Acton.genind, B = 1000), digits = 3)
round(pegas::hw.test(Brantford.genind, B = 1000), digits = 3)
round(pegas::hw.test(Elmira.genind, B = 1000), digits = 3)
round(pegas::hw.test(Everett.genind, B = 1000), digits = 3)
round(pegas::hw.test(Fergus.genind, B = 1000), digits = 3)
round(pegas::hw.test(Guelph.genind, B = 1000), digits = 3)
round(pegas::hw.test(`Port Hope.genind`, B = 1000), digits = 3)
round(pegas::hw.test(Waterloo.genind, B = 1000), digits = 3)

HWE.test <- data.frame(sapply(seppop(Acton.genind), 
                              function(ls) pegas::hw.test(ls, B=1000)[,4]))
HWE.test.MC <- t(data.matrix(HWE.test))
{cat("Chi-squared test (p-values):", "\n")
  round(HWE.test.MC,3)}

alpha=0.05
Prop.loci.out.of.HWE <- data.frame(MC = apply(HWE.test.MC < alpha, 2, mean))
Prop.loci.out.of.HWE

Prop.pops.out.of.HWE <- data.frame(MC = apply(HWE.test.MC < alpha, 1, mean))
Prop.pops.out.of.HWE

MC.fdr <- matrix(p.adjust(HWE.test.MC, method="fdr"), 
                 nrow = nrow(HWE.test.MC))

Prop.pops.out.of.HWE <- data.frame(MC = apply(HWE.test.MC < alpha, 1, mean),
                                   MC.fdr = apply(MC.fdr < alpha, 1, mean))
Prop.pops.out.of.HWE 
