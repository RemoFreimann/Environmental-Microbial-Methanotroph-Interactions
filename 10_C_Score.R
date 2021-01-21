install.packages("EcoSimR")

library(EcoSimR)
library(dplyr)
library(tibble)
library(vegan)
## the data is loadad by running the pmoA and 16S rRNA workflow scripts first

## mobs
dataCscore.pmoA= decostand((otu_table(bac.rotsee.pmoA)@.Data), "pa")
(dataCscore.pmoA=rownames_to_column(as.data.frame(dataCscore.pmoA), "Species"))

pmoA_Csore_Mod <- cooc_null_model(dataCscore.pmoA, algo="sim9",nReps=50000,burn_in = 1000)

summary(pmoA_Csore_Mod)
plot(pmoA_Csore_Mod,type="hist")
plot(pmoA_Csore_Mod,type="cooc")

## 16s
dataCscore.16s= decostand((otu_table(bac.rotsee.16s)@.Data), "pa")
(dataCscore.16s=rownames_to_column(as.data.frame(dataCscore.16s), "Species"))

bac16s_Csore_Mod <- cooc_null_model(dataCscore.16s, algo="sim9",nReps=10000,burn_in = 5000)

summary(bac16s_Csore_Mod)
plot(bac16s_Csore_Mod,type="hist")
plot(bac16s_Csore_Mod,type="cooc")

################################################################
## End of the script
################################################################