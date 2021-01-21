################################################################
# load the dataset and packages 
################################################################
install.packages("stringr")

library(ggplot2)
library(phyloseq)
theme_set(theme_bw())
library(reshape2)
library(vegan)
library(fields)
library(RColorBrewer)
library(mvabund)
library(tcltk)  
library(foreach)
library(pheatmap)
library(grDevices)
library(car)
library(leaps)
library(multcomp)
library(relaimpo)
library(MASS)
library(RcmdrMisc)
library(gplots)
library(cooccur)
library(plyr)
library(doParallel)
library(missMDA)
library(corrplot)
#library(glmulti)
library(psadd)
library(ape)
library(Hmisc)
library(corrplot)
library(igraph)


## data ... take accordingly from the pmoa and 16s scripts (i.e. if you want another normalization) or prepare as proposed below

##---------------------------------------------------------------------------------------------------------------------------------  
## pmoA data 
##--------------------------------------------------------------------------------------------------------------------------------- 
otufile     <- "data/P349_ZOTU_sintax.tab"  
mapfile     <- "data/Metafile_imputed.txt" 
treefile    <- "data/P349_ZOTU.tre"
refseqfile  <- "data/P349_ZOTU.fa"

## ---------------------------------------------------------------------------------------
## Create primary S4 data structure from the quime output files to be used with phyloseq
## ---------------------------------------------------------------------------------------
bac.rotsee <- import_qiime(otufilename = otufile, mapfilename = mapfile,
                           treefilename = treefile, refseqfilename = refseqfile)
bac.rotsee.original.pmoA <- bac.rotsee   # save the non-filtered data object

############################################################################################
# order otus according to phylogenetic distances and make a new phyloseq object
############################################################################################
(rotsee.otu.sort <-otu_table(bac.rotsee))
(rotsee.tax.sort <- tax_table(bac.rotsee))
(rotsee.tree.sort <- phy_tree(bac.rotsee))
(rotsee.env.sort <- sample_data(bac.rotsee))
(rotsee.refseq.sort <- refseq(bac.rotsee))

library(ggtree)
d=fortify(rotsee.tree.sort)
dd = subset(d, isTip)
phylogenetic_order_vector=dd$label[order(dd$y, decreasing=TRUE)]
(rotsee.otu.sort <- rotsee.otu.sort[match(phylogenetic_order_vector,rownames(rotsee.otu.sort)), ])
(rotsee.tax.sort <- rotsee.tax.sort[match(phylogenetic_order_vector,rownames(rotsee.tax.sort)), ])

## detach and reload libraries due to some incompatibilities of ggtree and its dependencies
detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}

############################################################################################
## reload the packages as "ggtree" it interferes with "phylo" class
############################################################################################
detachAllPackages()
library(ggplot2)
library(phyloseq)
theme_set(theme_bw())
library(reshape2)
library(vegan)
library(fields)
library(RColorBrewer)
library(mvabund)
library(tcltk)  
library(foreach)
library(pheatmap)
library(grDevices)
library(car)
library(leaps)
library(multcomp)
library(relaimpo)
library(MASS)
library(RcmdrMisc)
library(gplots)
library(cooccur)
library(plyr)
library(doParallel)
library(missMDA)
library(corrplot)
#library(glmulti)
library(psadd)
library(ape)
library(Hmisc)
library(corrplot)
library(igraph)

############################################################################################
## merge to a new phylsoeq object
############################################################################################
(bac.rotsee = merge_phyloseq(rotsee.otu.sort, rotsee.tax.sort, rotsee.env.sort, rotsee.tree.sort , rotsee.refseq.sort))

############################################################################################
## Check and remove OTUs with ZERO reads (here there are mainly OTUs removed that were only present in Magdalenas data set)
############################################################################################
bac.rotsee.prune <- bac.rotsee
bac.rotsee.prune <- prune_taxa(taxa_sums(bac.rotsee.prune) > 0, bac.rotsee.prune)
ntaxa(bac.rotsee.prune) 

############################################################################################
## Filter OTUs that appear less than x times (i.e. read counts) in at least A sample counts
############################################################################################
filt3 <- genefilter_sample(bac.rotsee.prune, filterfun_sample(function(x) x > 50), A=3)   
## more than 50 reads, in at least 3 samples
bac.rotsee.prune.filt3 <- prune_taxa(filt3, bac.rotsee.prune)
ntaxa(bac.rotsee.prune.filt3) 

############################################################################################
##  standardize samples to the median sequencing depth
##(same read counts, i.e. median of all read counts)
############################################################################################
total = median(sample_sums(bac.rotsee.prune.filt3))
standf = function(x, t=total) (t * (x / sum(x)))
bac.rotsee.prune.filt3.seqdepth = transform_sample_counts(bac.rotsee.prune.filt3, standf)
sample_sums(bac.rotsee.prune.filt3.seqdepth)

############################################################################################
##   Create datasets with relative abundance  (samples sum up to 1 each)
##  These Normalized datasets are for analyses that look at relative instead of absolute abundance 
############################################################################################
bac.rotsee.prune.ab <- transform_sample_counts(bac.rotsee.prune, function(x){x/sum(x)})
bac.rotsee.prune.filt3.seqdepth.ab <- transform_sample_counts(bac.rotsee.prune.filt3.seqdepth, function(x){x/sum(x)})
bac.rotsee.prune.filt3.rare.ab <- transform_sample_counts(bac.rotsee.prune.filt3.rare, function(x){x/sum(x)})

############################################################################################
## choose one of the  prefiltered data set for further investigation
## ---------------------------------------------------------------------------------------
## the herein choosen prefilterd data will be used for the rest of the script!!!! 
############################################################################################
bac.rotsee <-  bac.rotsee.prune.filt3.seqdepth # this data set will be converted to be used with vegan funcitonalities
## check sparsity
(pmoA.sparse<-plot_sparsity(bac.rotsee))
## check number of samples where abundance is not null 
taxa_prev(bac.rotsee)
halfsamples= ncol(sample_data(bac.rotsee))*0.5
sparsekeep=names(which(taxa_prev(bac.rotsee)>=halfsamples))
## you might want to use this subset for the pearson/spearman network analysis
pmoA.rotsee.spearman.network.phyloseq =  prune_taxa(sparsekeep, bac.rotsee)
##

## remove non-existent OTUs from the subset
ntaxa(bac.rotsee) 
bac.rotsee <- prune_taxa(taxa_sums(bac.rotsee) > 0, bac.rotsee)
ntaxa(bac.rotsee) 
bac.rotsee.pmoA = bac.rotsee

########################################################################################################################################
# produce OTU, taxanomic and environmental tables for the usage with other packackes 
########################################################################################################################################
## bac.rotsee =  bac.rotsee.prune.filt3.seqdepth = 90 OTU left
options(max.print = 2000000000)

(rotsee.otu.table <- otu_table(bac.rotsee))
(rotsee.tax.table <- tax_table(bac.rotsee))
(rotsee.env.table <- get_variable(bac.rotsee)) #mapfile 
ntaxa(bac.rotsee)
# set negative values in the metafile to 0, except for 13C-CH4 isotopes
rotsee.env.table[,-which(names(rotsee.env.table)=="X13CH4")][rotsee.env.table[,-which(names(rotsee.env.table)=="X13CH4")] < 0] <- 0 
#rotsee.env.table.log = cbind( rotsee.env.table[,c(1:13,15,16,17,37)],
#                              log1p(rotsee.env.table[,-c(1:13,15,16,17, 35,   48:56, 59:67)]))
rotsee.env.table=rotsee.env.table[,-1]
#rotsee.env.table=rotsee.env.table[,-which(names(rotsee.env.table)=="X")]
names(rotsee.env.table)
rotsee.env.table.log = cbind( rotsee.env.table[,c(1:13)],
                              log1p(rotsee.env.table[,c(14),drop=FALSE]),
                              rotsee.env.table[,c(15:17)],
                              log1p(rotsee.env.table[,c(18:34)]),
                              #log1p(abs(rotsee.env.table[,c(35),drop=FALSE])),
                              log1p(rotsee.env.table[,c(35), drop=FALSE] + (1+abs(min(rotsee.env.table$X13CH4)))),
                              log1p(rotsee.env.table[,c(36),drop=FALSE]),
                              rotsee.env.table[,c(37),drop=FALSE],
                              log1p(rotsee.env.table[,c(38:70)]))

head(rotsee.env.table.log)

####################
## remove impuded variables that have no data for one or more seasons
bad.impudation = c("NH4", "PO4", "DIC", "NO2", "NO3", "X13CH4", "Chla", "SO4")
rotsee.env.table      =  rotsee.env.table[,!names(rotsee.env.table) %in% bad.impudation]
rotsee.env.table.log  =  rotsee.env.table.log[,!names(rotsee.env.table.log) %in% bad.impudation]

########################################################################################################################################
# function to produce data sets that can be used for multivariate methods based on the vegan package
########################################################################################################################################
rotsee.otu.table.t<-t(rotsee.otu.table)     # this OTU table can be read by vegan functions, it will plot OTU numbers. See the function below...
rotsee.vegan.otu.table <- rotsee.otu.table.t

#rotsee.asigned.otu.table<- merge(rotsee.otu.table, rotsee.tax.table, by="row.names", all=TRUE) # a merged table

## function which asigns the deepest available Taxonomic rank ("Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species") to a specific OTU

make_vegan_otu_table <- function (otu.table.qiime, tax.table.qiime) {
  
  otu.table.qiime.t<-t(otu.table.qiime)
  otu.table.qiime.to.vegan <- otu.table.qiime.t
  
  
  for (i in c(1:length(colnames(otu.table.qiime.t)))) {
    
    if (!is.na(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Species"][[1]]) &
        toString(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Species"][[1]])!="uncultured_bacterium"){
      colnames(otu.table.qiime.to.vegan)[i]<- tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Species"][[1]]
    }                      # use this bit of code if analysis goes down to Family Level...replace following "if" statement with "else if"
    
    else if (!is.na(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Genus"][[1]]) &
             toString(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Genus"][[1]])!="uncultured_bacterium"){
      colnames(otu.table.qiime.to.vegan)[i]<- tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Genus"][[1]]
    }  
    
    else if (!is.na(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Family"][[1]]) &
             toString(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Family"][[1]])!="uncultured_bacterium"){
      colnames(otu.table.qiime.to.vegan)[i]<- tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Family"][[1]]
    }
    
    else if (!is.na(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Order"][[1]]) &
             toString(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Order"][[1]])!="uncultured_bacterium"){
      colnames(otu.table.qiime.to.vegan)[i]<- tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Order"][[1]]
    }
    
    
    else if (!is.na(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Class"][[1]]) &
             toString(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Class"][[1]])!="uncultured_bacterium"){
      colnames(otu.table.qiime.to.vegan)[i]<- tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Class"][[1]]
    }
    
    
    else if (!is.na(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Phylum"][[1]]) &
             toString(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Phylum"][[1]])!="uncultured_bacterium"){
      colnames(otu.table.qiime.to.vegan)[i]<- tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Phylum"][[1]]
    }
    
    
    # use this bit of code if analysis goes down to Species Level...replace following "if" statement with "else if"
    
    else {colnames(otu.table.qiime.to.vegan)[i]<-"unknown OTU"}
  }
  
  return(otu.table.qiime.to.vegan)
}


rotsee.vegan.otu.table<- make_vegan_otu_table (rotsee.otu.table, rotsee.tax.table) # this is a table with the deepest Phylogenetic information gathered by searching respective databases in QIIME
rotsee.vegan.otu.table.pmoA = rotsee.vegan.otu.table    # will be used for network analysis 

# quickly check if everything went well   
colnames(rotsee.vegan.otu.table)[which(rownames(rotsee.otu.table)=="ZOTU98")]  # should be a RPCS
colnames(rotsee.vegan.otu.table)[which(rownames(rotsee.otu.table)=="ZOTU50")]  # should be a Methylobacter_sp._BB5

## standardize or not (if not use only spearman correlation)
#ClusterHelPearsonpmoA   <-rotsee.vegan.otu.table.pmoA
#ClusterHelPearsonpmoA   <-rotsee.vegan.otu.table.pmoA / rowSums(rotsee.vegan.otu.table.pmoA)
ClusterHelPearsonpmoA   <-decostand(rotsee.vegan.otu.table.pmoA, "hellinger")

(colnames(ClusterHelPearsonpmoA) <- paste(colnames(ClusterHelPearsonpmoA),"_",rownames(rotsee.otu.table), sep=""))
(rotsee.tax.table.pmoA = tax_table(bac.rotsee.pmoA))


##--------------------------------------------------------------------------------------------------------------------------------- 
## 16s data
##--------------------------------------------------------------------------------------------------------------------------------- 
otufile     <- "data/P349_OTU_ab2_utax.tab"          
mapfile     <- "data/Metafile_imputed.txt"          
treefile    <- "data/P349_OTU_ab2.tre"
refseqfile  <- "data/P349_OTU_ab2.fa"

## ---------------------------------------------------------------------------------------
# Create primary S4 data structure from the quime output files
## ---------------------------------------------------------------------------------------
bac.rotsee <- phyloseq:::import_qiime(otufilename = otufile, mapfilename = mapfile,
                                      treefilename = treefile, refseqfilename = refseqfile)
bac.rotsee.original.16s <- bac.rotsee              # save the non-filtered data object

minusvecsites=c("A04")                        # remove as no MOBs are present
bac.rotsee <- subset_samples(bac.rotsee, !X.SampleID == minusvecsites) 

## ---------------------------------------------------------------------------------------
# order otus according to phylogenetic distances and make a new phyloseq object
## ---------------------------------------------------------------------------------------
(rotsee.otu.sort <-otu_table(bac.rotsee))
(rotsee.tax.sort <- tax_table(bac.rotsee))
(rotsee.tree.sort <- phy_tree(bac.rotsee))
(rotsee.env.sort <- sample_data(bac.rotsee))
(rotsee.refseq.sort <- refseq(bac.rotsee))
## order the otu table according to phylogenetic tree
library(ggtree)
d=fortify(rotsee.tree.sort)
dd = subset(d, isTip)
phylogenetic_order_vector=dd$label[order(dd$y, decreasing=TRUE)]
(rotsee.otu.sort <- rotsee.otu.sort[match(phylogenetic_order_vector,rownames(rotsee.otu.sort)), ])
(rotsee.tax.sort <- rotsee.tax.sort[match(phylogenetic_order_vector,rownames(rotsee.tax.sort)), ])
## detach and reload libraries due to some incompatibilities of ggtree and its dependencies
detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}
detachAllPackages()

library(ggplot2)
library(phyloseq)
theme_set(theme_bw())
library(reshape2)
library(vegan)
library(fields)
library(RColorBrewer)
library(mvabund)
library(tcltk)  
library(foreach)
library(pheatmap)
library(grDevices)
library(car)
library(leaps)
library(multcomp)
library(relaimpo)
library(MASS)
library(RcmdrMisc)
library(gplots)
library(cooccur)
library(plyr)
library(doParallel)
library(missMDA)
library(corrplot)
#library(glmulti)
library(psadd)
library(ape)
library(Hmisc)
library(corrplot)
library(igraph)

## merge the 16S rRNA data set
(bac.rotsee = merge_phyloseq(rotsee.otu.sort, rotsee.tax.sort, rotsee.env.sort, rotsee.tree.sort, rotsee.refseq.sort ))

# Removing unwanted taxa (i.e. chroloplast sequence based OTUs)
## ---------------------------------------------------------------------------------------
## !Archaea & Eukary
bac.rotsee.prune <- phyloseq::subset_taxa(bac.rotsee, Kingdom=="Bacteria")
ntaxa(bac.rotsee.prune) 

#save chlorplast data, just in case we need it later
chlpst <- phyloseq::subset_taxa(bac.rotsee, Class %in% "Chloroplast")
ntaxa(chlpst) 

## Remove chloroplasts from Bacterial OTU table
bac.rotsee.prune <- phyloseq::subset_taxa(bac.rotsee.prune, !(Class %in% "Chloroplast"))
ntaxa(bac.rotsee.prune)

## Check and remove OTUs with ZERO reads
bac.rotsee.prune <- phyloseq::prune_taxa(taxa_sums(bac.rotsee.prune) > 0, bac.rotsee.prune)
ntaxa(bac.rotsee.prune) 
#--> Use this data set for analyses requiring unfiltered datatsets, e.g. diversity. This data contains many singletons

## ---------------------------------------------------------------------------------------
#   FILTERING (should be done according to your hypthesis, i.e. if you're looking for rare
#              species also or just want to process the main players)
## ---------------------------------------------------------------------------------------
bac.rotsee.prune.filt3 <- bac.rotsee.prune

## Filter OTUs that appear at least in two different samples
filt3 <- genefilter_sample(bac.rotsee.prune, filterfun_sample(function(x) x > 2), A=1)
bac.rotsee.prune.filt3 <- prune_taxa(filt3, bac.rotsee.prune)
ntaxa(bac.rotsee.prune.filt3)

# filt3 <- genefilter_sample(bac.rotsee.prune, filterfun_sample(function(x) x >= 3), A=2)   
# bac.rotsee.prune.filt3 <- prune_taxa(filt3, bac.rotsee.prune)
# ntaxa(bac.rotsee.prune.filt3)

##  or 
## filter 16s for most abundant and present OTUs i.e. 1% of all reads and at least in 10% of the samples
# (bac.rotsee.prune.filt3 = filter_taxa( bac.rotsee.prune, function(x) sum(x > sample_sums( bac.rotsee.prune)/100*0.01) > (0.1*length(x)), TRUE)) 
#  ntaxa(bac.rotsee.prune.filt3)


## here we quickly check sparsity for spearman/pearson network based analysis (see other scripts later on)
#(bac.sparse<-plot_sparsity(bac.rotsee.prune.filt3))
# taxa_prev(bac.rotsee.prune.filt3)
# halfsamples= ncol(sample_data(bac.rotsee.prune.filt3))*0.5
# sparsekeep=names(which(taxa_prev(bac.rotsee.prune.filt3)>=halfsamples))
# bac.rotsee.spearman.network.phyloseq =  prune_taxa(sparsekeep, bac.rotsee)

## ---------------------------------------------------------------------------------------
## remove MOBs (based on 16S rRNA indetification) 
## ---------------------------------------------------------------------------------------

MOBs = c( "OTU_129",         # Methylosinus          --> sMMO
          "OTU_68",          # Crenothrix            
          "OTU_56",          # Crenothrix            
          "OTU_141",         # Crenothrix            
          "OTU_931",         # Crenothrix
          "OTU_575",         # Crenothrix
          "OTU_433",         # Crenothrix
          "OTU_614",         # Crenothrix
          "OTU_917",         # Crenothrix
          "OTU_319",         # Methylomonas         
          "OTU_2426",        # Methylomonas         
          "OTU_821",         # Methylomonas         --> sMMO
          "OTU_663",         # Methylocaldum        --> sMMO
          "OTU_2276",        # Methylocaldum        --> sMMO
          "OTU_3421",        # Methylocaldum        --> sMMO
          "OTU_3248",        # Methylocaldum        --> sMMO
          "OTU_112",         # Methylacidiphilales  LD19
          "OTU_1686",        # Methylacidiphilales  LD19
          "OTU_566",         # Methylacidiphilales  LD19
          "OTU_1479",        # Methylacidiphilales 
          "OTU_1541"         # Methylacidiphilales 
)

smmo_vector<- c("OTU_129", "OTU_663", "OTU_2276", "OTU_3421", "OTU_3248")     ## vector appointing the sMMO containing OTUs


## produce 16s rRNA data excluding MOBs to be used in combination with the pMOA data 
(allnames = taxa_names(bac.rotsee.prune.filt3))    
(noMOBs <- allnames[! allnames %in% MOBs])
bac.rotsee.prune.noMOBS <- prune_taxa(noMOBs, bac.rotsee.prune.filt3)
bac.rotsee.MOBs = prune_taxa(MOBs, bac.rotsee.prune)

##  produce phyloseq objects with 16s rRNA only MOBs, excluding MOBs or including all 
(rotsee.otu.table <-otu_table(bac.rotsee.prune.noMOBS))
(rotsee.tax.table <- tax_table(bac.rotsee.prune.noMOBS))
(rotsee.tree.table <- phy_tree(bac.rotsee.prune.noMOBS))
(rotsee.env.table <- sample_data(bac.rotsee.prune.noMOBS))
(rotsee.refseq <- refseq(bac.rotsee.prune.noMOBS))
(bac.rotsee.no.MOBs =  merge_phyloseq(rotsee.otu.table, rotsee.tax.table, rotsee.tree.table , rotsee.env.table,  rotsee.refseq) )

## ***************************************************************************************
## Data prefiltering (standardization)
## ***************************************************************************************

## ---------------------------------------------------------------------------------------
## Standardization of sample reads to the median sequencing depth applied to each sample  (readcounds, 100% and 1)
## ---------------------------------------------------------------------------------------
(total = median(sample_sums(bac.rotsee.prune.noMOBS)) )
standf = function(x, t=total) (t * (x / sum(x)))
bac.rotsee.no.MOBs.seqdepth = transform_sample_counts(bac.rotsee.prune.noMOBS, standf)  

## quickly check colsums (i.e. reads)
colSums(otu_table(bac.rotsee.no.MOBs.seqdepth))

## ***************************************************************************************
## ---------------------------------------------------------------------------------------
## choose one of the  prefiltered data set for further investigation
## ---------------------------------------------------------------------------------------
## the herein choosen prefilterd data will be used for the rest of the script!!!! 
## ***************************************************************************************

## remove non-existent OTUs from the newly generated subsets
bac.rotsee.no.MOBs.seqdepth      
ntaxa(bac.rotsee.no.MOBs.seqdepth)  
bac.rotsee.nomobs <- prune_taxa(taxa_sums(bac.rotsee.no.MOBs.seqdepth) > 0, bac.rotsee.no.MOBs.seqdepth)
ntaxa(bac.rotsee.nomobs) 

## ---------------------------------------------------------------------------------------
## again, choose the according dataset for further analysis (i.e. including or excluding MOBs etc.)
## -------------------------------------------------------------------------------------
## check sparsity again
#(bac.sparse<-plot_sparsity(bac.rotsee.nomobs))
# taxa_prev(bac.rotsee.nomobs)
# halfsamples= ncol(sample_data(bac.rotsee.nomobs))*0.5
# sparsekeep=names(which(taxa_prev(bac.rotsee.nomobs)>=halfsamples))
# bac.rotsee.nomobs.spearman.network.phyloseq =  prune_taxa(sparsekeep, bac.rotsee.nomobs)      ## you might want to use this subset for the pearson/spearman network analysis later on

rotsee.otu.table.16s <- otu_table(bac.rotsee.nomobs)
rotsee.tax.table.16s <- tax_table(bac.rotsee.nomobs)

## make a vegan object
rotsee.nomobs.otu.table.t<-t(rotsee.otu.table.16s)

## run the make vegan function
rotsee.vegan.otu.table <- make_vegan_otu_table (rotsee.otu.table.16s, rotsee.tax.table.16s) # this is a table with the deepest Phylogenetic information gathered by searching respective databases in QIIME
rotsee.vegan.otu.table.16s = rotsee.vegan.otu.table    # will be used for network analysis 

# quickly check if everything went well   
colnames(rotsee.vegan.otu.table)[which(rownames(rotsee.otu.table.16s)=="OTU_1686")]  # should be a LD19 or character(0) if MOBs were removed (i.e. if rotsee.nomobs.otu.table.t was used)
colnames(rotsee.vegan.otu.table)[which(rownames(rotsee.otu.table.16s)=="OTU_50")]  # should be a Isosphaeraceae

(colnames(rotsee.vegan.otu.table.16s) <- paste(colnames(rotsee.vegan.otu.table.16s),"_",rownames(rotsee.otu.table.16s), sep=""))

## standardize or not (if not use spearman)
# ClusterHelPearson16s     <-rotsee.vegan.otu.table.16s
# ClusterHelPearson16s    <-rotsee.vegan.otu.table.16s / rowSums(rotsee.vegan.otu.table.16s)
ClusterHelPearson16s    <-decostand(rotsee.vegan.otu.table.16s, "hellinger") 

## remove sites A1 to A5 and D12 as no MOBs based on pmoA are present
ix= which(rownames(ClusterHelPearson16s)%in% c("A01", "A02" ,"A03", "A05", "D12"))
ClusterHelPearson16s = ClusterHelPearson16s[-ix,]

##------------------------------------------------------------------------------
## merge pmoA and 16S rRNA data for the network analysis
##------------------------------------------------------------------------------
# weight between matrices
ClusterHelPearson16s<-ClusterHelPearson16s/(ncol(rotsee.vegan.otu.table.16s)+ncol(ClusterHelPearsonpmoA))*ncol(rotsee.vegan.otu.table.16s)
ClusterHelPearsonpmoA<-ClusterHelPearsonpmoA/(ncol(rotsee.vegan.otu.table.16s)+ncol(ClusterHelPearsonpmoA))*ncol(ClusterHelPearsonpmoA)

## make matrix for cooccurence network
# reorder data frames to match sampling order
ClusterHelPearsonpmoA=ClusterHelPearsonpmoA[match(rownames(ClusterHelPearson16s), rownames(ClusterHelPearsonpmoA)), ]
rcorrMatrix<-as.matrix(cbind(ClusterHelPearsonpmoA  ,  ClusterHelPearson16s))  
rcorrMatrix<-decostand(rcorrMatrix, "hellinger")
(removezerovector = which(colSums(rcorrMatrix)==0))
rcorrMatrix   <- rcorrMatrix[,!colnames(rcorrMatrix) %in% names(removezerovector)]
rcorrMatrix <- as.matrix(rcorrMatrix)

removemissing16sOTUsvector = which(colnames(rotsee.vegan.otu.table.16s) %in% names(removezerovector))
rotsee.tax.table.16s =rotsee.tax.table.16s[-removemissing16sOTUsvector,]

##------------------------------------------------------------------
## cooccurence network analysis function
##------------------------------------------------------------------

library(stringr)

co_occurrence_network<-function(matrix, alpha, p.cutoff){
  
  #correlation analysis based on spearman or pearson's co-efficient
  matrix.dist<-rcorr((matrix),type="spearman") 
  matrix.cor<-matrix.dist$r
  matrix.cor.p<-matrix.dist$P
  
  #Multiple testing correction using Benjamini-Hochberg standard false discovery rate correction ("FDR-BH")
  matrix.cor.p <- p.adjust(matrix.cor.p, method="BH")
  
  #1.Consider positive cooccurence at given coefficient (alpha) and p-value cutoffs
  matrix.cor1<-matrix.cor
  matrix.cor1.p<-matrix.cor.p
  matrix.cor1[which(matrix.cor1 <= alpha)]=0
  matrix.cor1[which(matrix.cor1.p > p.cutoff)]=0
  # drop those rows and columns with sum = 0
  matrix.cor1<-matrix.cor1[which(rowSums(matrix.cor1)!=1),]
  matrix.cor1<-matrix.cor1[,which(colSums(matrix.cor1)!=0)]
  
  
  
  #2.Consider netagive cooccurence at given coefficient (-alpha) and p-value cutoffs
  matrix.cor2<-matrix.cor
  matrix.cor2.p<-matrix.cor.p
  matrix.cor2[which(matrix.cor2 >= (-alpha))]=0
  matrix.cor2[which(matrix.cor2.p > p.cutoff)]=0
  # drop those rows and columns with sum = 0
  matrix.cor2<-matrix.cor2[which(rowSums(matrix.cor2)!=1),]
  #matrix.cor2<-matrix.cor2[,which(colSums(matrix.cor2)!=0)]
  
  
  
  #3.Consider both positive and netagive cooccurence at given coefficient (alpha) and p-value cutoffs
  matrix.cor3<-matrix.cor
  matrix.cor3.p<-matrix.cor.p
  matrix.cor3[which(matrix.cor3>=(-alpha) & matrix.cor3 <= alpha)]=0
  matrix.cor3[which(matrix.cor3.p>p.cutoff)]=0
  # drop those rows and columns with sum = 0
  matrix.cor3<-matrix.cor3[which(rowSums(matrix.cor3)!=1),]
  matrix.cor3<-matrix.cor3[,which(colSums(matrix.cor3)!=0)]
  
  
  
  # generate graphs using igraph
  g1<-graph.adjacency(matrix.cor1,weight=T,mode="undirected")
  g1<-simplify(g1)
  V(g1)$label <- V(g1)$name
  V(g1)$degree <- degree(g1)
  # add names of specific phylogenetic depths to the igraph object
  my.list.g1 <- row.names(matrix.cor1)
  logical.g1 <- colnames(matrix)  %in% my.list.g1
  namevector.gephipmoA <-as.data.frame((rotsee.tax.table.pmoA[,"Genus"]))
  colnames(namevector.gephipmoA)[1]="Group"
  namevector.gephi16s <-as.data.frame((rotsee.tax.table.16s[,"Phylum"]))
  colnames(namevector.gephi16s)[1]="Group"
  gephinamevector <- rbind(namevector.gephipmoA, namevector.gephi16s)
  groupnames.gephi=as.character(gephinamevector[,][logical.g1])
  V(g1)$groupname <- groupnames.gephi
  pmoAinMatrix <- str_count(dimnames(matrix.cor1)[1], "_ZOTU")
  V(g1)$pointsize <- c(rep(50, pmoAinMatrix), rep(10, length(groupnames.gephi) - pmoAinMatrix))    ## change according to groupnames gephi after run till last line
  
  
  g2<-graph.adjacency(abs(matrix.cor2),weight=T,mode="undirected")
  g2<-simplify(g2)
  V(g2)$label <- V(g2)$name
  V(g2)$degree <- degree(g2)
  
  g3<-graph.adjacency(matrix.cor3,weight=T,mode="undirected")
  g3<-simplify(g3)
  V(g3)$label <- V(g3)$name
  V(g3)$degree <- degree(g3)
  
  
  # append the output into results
  result<-list()
  
  result$matrix.cor<-matrix.cor
  result$matrix.cor.p<-matrix.cor.p
  
  result$matrix.cor1<-matrix.cor1
  result$graph1<-g1
  
  result$matrix.cor2<-matrix.cor2
  result$graph2<-g2
  
  result$matrix.cor3<-matrix.cor3
  result$graph3<-g3
  return(result)
  
} 

##########################################################################
## run the function
pattern = co_occurrence_network(rcorrMatrix , alpha=0.7, p.cutoff=0.05)


##2. Creating gml files of network (to be visulized in Gephi or Cytoscape)
write.graph(pattern$graph1,'Pos0.7-NW.gml',format='gml')
write.graph(pattern$graph2,'Neg0.7-NW.gml',format='gml')
write.graph(pattern$graph3,'PosNeg0.7-NW.gml',format='gml')

###3. Calculating network topological properties(for prositive network)
g<-pattern$graph1
c <- cluster_walktrap(g)
# Global toplogical features
modularity(c)
md <- modularity(g, membership(c), weights = NULL)
cc <- transitivity(g, vids = NULL,
                   weights = NULL)
spl <- average.path.length(g, directed=FALSE, unconnected=TRUE)
gd  <- graph.density(g, loops=FALSE)
nd  <- diameter(g, directed = FALSE, unconnected = TRUE, weights = NULL)

node.degree <- degree(g, v = V(g), mode="all")
ad  <- mean(node.degree)

e <- ecount(g)
v <- vcount(g)

e
v

global.topology <- data.frame(e,v,cc,spl,md,gd,nd,ad)
write.csv(global.topology, file="Pos0.6-NW-global.topology.csv")

# Node toplogical features
betweenness.centrality <- betweenness(g, v=V(g), 
                                      directed = FALSE, weights = NULL,
                                      nobigint = TRUE, normalized = FALSE)
closeness.centrality <- closeness(g, vids = V(g),
                                  weights = NULL, normalized = FALSE)
node.transitivity <- transitivity(g, type = c("local"), vids = NULL,
                                  weights = NULL)

node.topology <- data.frame(node.degree, betweenness.centrality, closeness.centrality, node.transitivity)
write.csv(node.topology, file="Pos0.6-NW-node.topology.csv")

# Ggploting node degreee distribution in a log-log plot
degree.df <- data.frame(table(degree=factor(node.degree, levels=seq_len(max(node.degree)))))
degree.df$degree <- as.numeric(as.character(degree.df$degree))

#4. Creating an abundance table for OTUs present in the positive and negative network
my.list1 <- row.names(pattern$matrix.cor1)
my.list2 <- row.names(pattern$matrix.cor2)

logical1 <- colnames(rcorrMatrix)  %in% my.list1
logical2 <- colnames(rcorrMatrix)  %in% my.list2

tab.subset1 <- subset(t(rcorrMatrix),logical1)
tab.subset2 <- subset(t(rcorrMatrix),logical2)

write.table(tab.subset1,'Pos0.6-NW.txt',sep="\t")
write.table(tab.subset2,'Neg0.6-NW.txt',sep="\t")

# Creating table with the names according to taxonomy depth coloring in gephi
namevector.gephipmoA <-as.data.frame((rotsee.tax.table.pmoA[,"Genus"]))
colnames(namevector.gephipmoA)[1]="Group"
namevector.gephi16s <-as.data.frame((rotsee.tax.table.16s[,"Phylum"]))
colnames(namevector.gephi16s)[1]="Group"
gephinamevector <- rbind(namevector.gephipmoA, namevector.gephi16s)
(groupnames.gephi=gephinamevector[,][logical1])
length(groupnames.gephi)


################################################################
## End of the script
################################################################