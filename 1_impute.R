#######################################################################################
## Impute missin data in environmental table 
#######################################################################################
library(missMDA)


#######################################################################################
## load data set to be imputed and quickly check its integrity.
## This is a CSV file with Sample IDs as first column and all variables arranged
## as collumns. Only columns with missing (numerical) can be imputed.
#######################################################################################
rotsee.env.table = read.csv("Metafile_to impute.csv" , header = TRUE)
names(rotsee.env.table)
rownames(rotsee.env.table) = rotsee.env.table[,1]
head(rotsee.env.table)

#######################################################################################
## Impute the data and save the new values as a new table for downstream analysis.
## If you dont want to use imputation you can exclude samples 
## with missing values according to the vectors (i.e. column numbers) introduced below
#######################################################################################
# check no of comps to be used for imputation (i.e. continuous variables)
names(rotsee.env.table)
(miss_env_dim = estim_ncpPCA( rotsee.env.table[,c(14,18:47,68,70,71)],
                              ncp.max = 15, scale=TRUE, nbsim = 10000, threshold=1e-6, verbose=TRUE ))
(no_ncp = miss_env_dim$ncp)

## check uncertainty of imputation
## by simulation of PCA
mi= MIPCA(rotsee.env.table[,c(14,18:47,68,70,71)], ncp= no_ncp, nboot=1000, scale=TRUE,
          method="Regularized", method.mi="Bayes", verbose=TRUE, threshold=1e-06)
plot(mi, axes=c(1,2), new.plot=TRUE)
## by overimputation
res.over<-Overimpute(mi)

## If you are happy with the imputed variables quality proceede
## with the actual imputation
rotsee.env.miss = imputePCA( rotsee.env.table[,c(14,18:47,68,70,71)], ncp= no_ncp, scale=TRUE ,
                             threshold=1e-06, method="Regularized")
rotsee.env.miss.table=rotsee.env.miss$completeObs

## bind the metafile (i.e. non imputet columns and imputed columns of the original data set)
rotsee.env.table.imp = cbind(rotsee.env.table[,c(1:13)], rotsee.env.miss.table[,1, drop=FALSE], rotsee.env.table[,c(15,16,17)],
                         rotsee.env.miss.table[,c(2:31)], rotsee.env.table[,c(48:67)], rotsee.env.miss.table[,c(32:34)] )

## set negative values in the metafile to 0, except for 13C-CH4 isotopes
rotsee.env.table.imp[,-which(names(rotsee.env.table.imp)=="X13CH4")][rotsee.env.table.imp[,-which(names(rotsee.env.table.imp)=="X13CH4")] < 0] <- 0



## make the sample IDs the column names
rownames(rotsee.env.table.imp) = rotsee.env.table.imp[,1]
names(rotsee.env.table.imp)


## write the file to folder
write.table(rotsee.env.table.imp, file = "Metafile_imputed.csv", sep = ",", col.names = NA,
            qmethod = "double")

write.table(rotsee.env.table.imp, file = "Metafile_imputed.txt", sep = "\t", col.names = NA,
            qmethod = "double")

## -> this metafile will be used for the whole analysis, impeted profiles were not used for donwstream analysis but rather for specific 
## interpreatations based on correlations with other variables.
## You will have to to remove the specific variables later on (i.e. "NO2", "NO3", "X13CH4", "SO4", "Chla", "PO4", "NH4", "DIC")
## Single impuded variables can be incorporated when not many data is missing.


## transform the csv to a tab delimited txt file for the phyloseq package!! and remove first column
##-------------------------------------------------------------------------

## Amelia can be used as an alternative based on EM algos
## see library(Amelia)

#######################################################################################
## end of script
#######################################################################################
