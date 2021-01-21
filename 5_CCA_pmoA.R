##-----------------------------------------------------------
## Data Import in R
## ---------------------------------------------------------------------------------------
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
#--> Use for analyses requiring unfiltered datatsets, e.g. diversity. This data contains many singletons

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


bac.rotsee <- bac.rotsee.prune.filt3.seqdepth

# bac.rotsee <- subset_samples(bac.rotsee.prune.filt3.seqdepth, Campaign=="A_June_2013")
# 
# bac.rotsee <- subset_samples(bac.rotsee.prune.filt3.seqdepth, Campaign=="B_August_2013")
# 
# bac.rotsee <- subset_samples(bac.rotsee.prune.filt3.seqdepth, Campaign=="C_September_2014")
# 
# bac.rotsee <- subset_samples(bac.rotsee.prune.filt3.seqdepth, Campaign=="D_September_2015")
# 
# bac.rotsee <- subset_samples(bac.rotsee.prune.filt3.seqdepth, Campaign!="D_June_2015") # without June

########################################################################################################################################
# produce OTU, taxanomic and environmental tables for the usage with other packackes 
########################################################################################################################################
## bac.rotsee =  bac.rotsee.prune.filt3.seqdepth = 90 OTU left
options(max.print = 2000000000)

(rotsee.otu.table <- otu_table(bac.rotsee))
(rotsee.tax.table <- tax_table(bac.rotsee))
(rotsee.env.table <- get_variable(bac.rotsee)) #mapfile 
ntaxa(bac.rotsee)
names(rotsee.env.table)


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
names(rotsee.env.table.log)

## remove impuded variables that have no data for one or more seasons
bad.impudation = c("NH4", "PO4", "DIC", "NO2", "NO3", "X13CH4", "Chla", "SO4")
rotsee.env.table      =  rotsee.env.table[,!names(rotsee.env.table) %in% bad.impudation]
rotsee.env.table.log  =  rotsee.env.table.log[,!names(rotsee.env.table.log) %in% bad.impudation]

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

# quickly check if everything went well   
colnames(rotsee.vegan.otu.table)[which(rownames(rotsee.otu.table)=="ZOTU98")]  # should be a RPCS
colnames(rotsee.vegan.otu.table)[which(rownames(rotsee.otu.table)=="ZOTU50")]  # should be a Methylobacter_sp._BB5

## you can use this function to write the produced OTU table into your working directory
# write.csv(rotsee.vegan.otu.table, "vegan_test.csv")

## function to merge otu and phylogenetic name
mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}
library(Biostrings)
refseq_filt3 = refseq(bac.rotsee.prune.filt3)
refseq_filt3@ranges@NAMES <- paste(mgsub(c("\\.", "/","-"),c("","_","_"),(refseq_filt3@ranges@NAMES)),"_",colnames(rotsee.vegan.otu.table), sep="") # replace characters that  are not parsed propperly


## set colors for Genus
genus.numbers <- length(unique(rotsee.tax.table[,"Genus"]))
#col.genus<-(colorRampPalette(brewer.pal(genus.numbers,"Reds"))(genus.numbers))  # choose this or the next color
col.genus.orig<-brewer.pal(genus.numbers, "PuOr")
genus.df<-transform(rotsee.tax.table, id=match(Genus, unique(Genus)))
test=unique(genus.df$Genus)
match(test, sort(unique(genus.df[,"Genus"]))) # copy the vector to the col.genus.orig
col.genus<-col.genus.orig[c( 2, 6, 3, 1, 5, 4)]
## assess size for species points to be plottet
bac.rotsee.prune.filt3.otu = otu_table(bac.rotsee.prune.filt3)
min(taxa_sums(bac.rotsee.prune.filt3.otu))
max(taxa_sums(bac.rotsee.prune.filt3.otu))
species.size.vec <- taxa_sums(bac.rotsee.prune.filt3.otu)
x=species.size.vec
normalized.species.size.vec = sqrt((x-min(x))/(max(x)-min(x))*3+0.5)

ClusterHel           <-decostand(rotsee.vegan.otu.table, "hellinger")     # square root of equally standardized row sums
ClusterNorm          <-decostand(rotsee.vegan.otu.table, method="total")  # standardized to equal row sum

## kick out some variables (i.e. qpcr copies)
names(rotsee.env.table.log)
(PhysicoChemicalLog=(rotsee.env.table.log[,-c(1:18,22,26, 30, 40:59)]))
(PhysicoChemicalLogsub=na.omit(rotsee.env.table.log[,-c(1:18,22,26, 30, 40:59)]))
names(PhysicoChemicalLogsub)

PhysicoChemicalRaw=(rotsee.env.table[,-c(1:18,22,26, 30, 40:59)])
PhysicoChemicalRawsub=na.omit(rotsee.env.table[,-c(1:18,22,26, 30, 40:59)])

ClusterHelsub <- ClusterHel[rownames(PhysicoChemicalLogsub),]

## check optically for some potential colinearity (i.e. with O2 or CH4 gradient)
pairs(PhysicoChemicalLogsub)

library(dplyr)
library(tidyr)
library(stringr)

PhysicoChemicalLogsub %>% gather(-O2, key = "var", value = "value") %>% 
ggplot(aes(x = value, y = O2)) +
  #stat_smooth() +
  geom_smooth(method="lm") +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

PhysicoChemicalLogsub %>%  gather(-CH4, key = "var", value = "value") %>%
  ggplot(aes(x = value, y = CH4)) +
  #stat_smooth() +
  geom_smooth(method="lm") +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

#######################################################################
## remove collinear variables prior to model selection
## stepwise VIF function with preallocated vectors
vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  require(fmsb)
  
  if(class(in_frame) != 'data.frame') in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-vector('list', length = ncol(in_frame))
  names(vif_init) <- names(in_frame)
  var_names <- names(in_frame)
  
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in<-formula(paste(val,' ~ .'))
    vif_init[[val]]<-VIF(lm(form_in,data=in_frame,...))
  }
  vif_max<-max(unlist(vif_init))
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(names(in_frame))
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-vector('list', length = ncol(in_dat))
      names(vif_vals) <- names(in_dat)
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in<-formula(paste(val,' ~ .'))
        vif_add<-VIF(lm(form_in,data=in_dat,...))
        vif_vals[[val]]<-vif_add
      }
      
      max_row<-which.max(vif_vals)[1]
      
      vif_max<-vif_vals[max_row]
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        vif_vals <- do.call('rbind', vif_vals)
        vif_vals
        prmatrix(vif_vals,collab='vif',rowlab=row.names(vif_vals),quote=F)
        cat('\n')
        cat('removed: ', names(vif_max),unlist(vif_max),'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% names(vif_max)]
      
    }
    
    return(names(in_dat))
    
  }
  
}


## Environmental data prep
## remove elements that with VIFs larger than 10 (you can do this also after the AIC further below)

(col.env.vif<- vif_func(in_frame=PhysicoChemicalLogsub,thresh=10,trace=T))
PhysicoChemicalLogsub = PhysicoChemicalLogsub[,c(col.env.vif)]

## you can additionally remove variables that might be corellated or are not of interest for you hypothesis
#remove=c("Mn_DGT", "Mn_diss_nM", "T")
#(PhysicoChemicalLogsub = PhysicoChemicalLogsub[ , !names(PhysicoChemicalLogsub) %in% remove ])

#######################################################################
## use the set for all seasons to select a subset for specific sampling campaigns
# 
# 
# 
# ## June 2013
# ## take the FORMULA subset from all seasons ... run the selection and remove colinear variables
# 
# ## first run this
# remove=c("DIC_mM", "T_C", "Mn_DGT_nM", "TDN_mM", "Mn_diss_nM", "Light_uE..m2s")
# (PhysicoChemicalLogsub = PhysicoChemicalLogsub[ , !names(PhysicoChemicalLogsub) %in% remove ])
# 
# 
# ## then come back to here after selection and proceed with the endmodel
# (Cluster.cca<-cca(FORMULA, PhysicoChemicalLogsub, na.action=na.exclude))   # we use the log transformed env variables
# vif.cca(Cluster.cca) 
# sort(vif.cca(Cluster.cca) )  # remove colinear variables >20, and add or remove terms if whished 
# (Cluster.cca= update(Cluster.cca, . ~ . - NO2._uM   ))
# vif.cca(Cluster.cca)
# sort(vif.cca(Cluster.cca) )
# (Cluster.cca= update(Cluster.cca, . ~ . - SO42._uM   ))
# vif.cca(Cluster.cca)
# sort(vif.cca(Cluster.cca) )
# (Cluster.cca= update(Cluster.cca, . ~ . - Cond_uS.cm   ))
# vif.cca(Cluster.cca)
# sort(vif.cca(Cluster.cca) )
# (Cluster.cca= update(Cluster.cca, . ~ . - PO43._uM   ))
# vif.cca(Cluster.cca)
# sort(vif.cca(Cluster.cca) )
# 
# (PhysicoChemicalLog=(rotsee.env.table.log[,-c(1:18,22,26, 30, 48:67)]))
# (PhysicoChemicalLogsub=na.omit(rotsee.env.table.log[,-c(1:18,22,26, 30, 48:67)]))
# keep=c("pH", "Cu_diss_nM", "O2_uM", "Fe_diss_nM", "Turb_NTU", "NO3._uM")
# (PhysicoChemicalLogsub = PhysicoChemicalLogsub[ , names(PhysicoChemicalLogsub) %in% keep ])
# 
# 
# 
# ## August 2014
# ## take the FORMULA subset from all seasons ... run the selection and remove colinear variables
# 
# ## first run this
# remove=c("DIC_mM", "T_C", "Mn_DGT_nM", "TDN_mM", "Mn_diss_nM", "Light_uE..m2s")
# (PhysicoChemicalLogsub = PhysicoChemicalLogsub[ , !names(PhysicoChemicalLogsub) %in% remove ])
# 
# 
# ## then come back to here after selection and proceed with the endmodel
# (Cluster.cca<-cca(FORMULA, PhysicoChemicalLogsub, na.action=na.exclude))   # we use the log transformed env variables
# vif.cca(Cluster.cca) 
# sort(vif.cca(Cluster.cca) )  # remove colinear variables >20, and add or remove terms if whished 
# (Cluster.cca= update(Cluster.cca, . ~ . - pH    ))
# vif.cca(Cluster.cca)
# sort(vif.cca(Cluster.cca) )
# 
# 
# (PhysicoChemicalLog=(rotsee.env.table.log[,-c(1:18,22,26, 30, 48:67)]))
# (PhysicoChemicalLogsub=na.omit(rotsee.env.table.log[,-c(1:18,22,26, 30, 48:67)]))
# keep=c("CH4_uM", "DOC_mM", "STot_uM", "Cond_uS.cm", "DOC_mM", "Zn_diss_nM", "Cr_tot_nM")
# (PhysicoChemicalLogsub = PhysicoChemicalLogsub[ , names(PhysicoChemicalLogsub) %in% keep ])
# 
# 
#            
# 
# ## September 2014
# ## take the FORMULA subset from all seasons ... run the selection and remove colinear variables
# 
# ## first run this
# remove=c("DIC_mM", "T_C", "Mn_DGT_nM", "TDN_mM", "Mn_diss_nM", "Light_uE..m2s")
# (PhysicoChemicalLogsub = PhysicoChemicalLogsub[ , !names(PhysicoChemicalLogsub) %in% remove ])
# 
# 
# ## then come back to here after selection and proceed with the endmodel
# (Cluster.cca<-cca(FORMULA, PhysicoChemicalLogsub, na.action=na.exclude))   # we use the log transformed env variables
# vif.cca(Cluster.cca) 
# sort(vif.cca(Cluster.cca) )  # remove colinear variables >20, and add or remove terms if whished 
# (Cluster.cca= update(Cluster.cca, . ~ . - O2_uM    ))
# vif.cca(Cluster.cca)
# sort(vif.cca(Cluster.cca) )
# (Cluster.cca= update(Cluster.cca, . ~ . - CH4_uM   ))
# vif.cca(Cluster.cca)
# sort(vif.cca(Cluster.cca) )
# (Cluster.cca= update(Cluster.cca, . ~ . - NH4._uM    ))
# vif.cca(Cluster.cca)
# sort(vif.cca(Cluster.cca) )
# (Cluster.cca= update(Cluster.cca, . ~ . - Cu_DGT_nM   ))
# vif.cca(Cluster.cca)
# sort(vif.cca(Cluster.cca) )
# 
# (PhysicoChemicalLog=(rotsee.env.table.log[,-c(1:18,22,26, 30, 48:67)]))
# (PhysicoChemicalLogsub=na.omit(rotsee.env.table.log[,-c(1:18,22,26, 30, 48:67)]))
# keep=c("SO42._uM", "Fe_diss_nM", "Cr_tot_nM", "Turb_NTU", "Cond_uS.cm")
# (PhysicoChemicalLogsub = PhysicoChemicalLogsub[ , names(PhysicoChemicalLogsub) %in% keep ])
# 
#        
# 
# 
# 
# 
# ## September 2015
# ## take the FORMULA subset from all seasons ... run the selection and remove colinear variables
# 
# ## first run this
# remove=c("DIC_mM", "T_C", "Mn_DGT_nM", "TDN_mM", "Mn_diss_nM", "Light_uE..m2s")
# (PhysicoChemicalLogsub = PhysicoChemicalLogsub[ , !names(PhysicoChemicalLogsub) %in% remove ])
# 
# 
# ## then come back to here after selection and proceed with the endmodel
# (Cluster.cca<-cca(FORMULA, PhysicoChemicalLogsub, na.action=na.exclude))   # we use the log transformed env variables
# vif.cca(Cluster.cca) 
# sort(vif.cca(Cluster.cca) )  # remove colinear variables >20, and add or remove terms if whished 
# 
# 
# (PhysicoChemicalLog=(rotsee.env.table.log[,-c(1:18,22,26, 30, 48:67)]))
# (PhysicoChemicalLogsub=na.omit(rotsee.env.table.log[,-c(1:18,22,26, 30, 48:67)]))
# keep=c("Cond_uS.cm", "CH4_uM", "PO43._uM", "X.13CH4_.")
# (PhysicoChemicalLogsub = PhysicoChemicalLogsub[ , names(PhysicoChemicalLogsub) %in% keep ])
# 
#######################################################################

###################################################################################
## perform model selections on non-colinear variables
###################################################################################
(Cluster.cca<-cca(ClusterHelsub~ ., PhysicoChemicalLogsub, na.action=na.omit)) 

# rerun from here for the single season after vifs are removed

(Cluster.cca.forward <- ordistep (cca(ClusterHelsub~1, PhysicoChemicalLogsub, na.action=na.omit), # ordiR2step must inherit from rda...cca not working
                                  scope= formula(Cluster.cca),
                                  direction="forward",pstep=10000,nperm=999999,  Pin = 0.05, Pout = 0.1))

(Cluster.cca.backward <- ordistep (cca(ClusterHelsub ~., PhysicoChemicalLogsub, na.action=na.omit),
                                   direction="backward",pstep=10000,nperm=999999,  Pin = 0.05, Pout = 0.1))

(Cluster.cca.both <- ordistep (cca(ClusterHelsub~1, PhysicoChemicalLogsub, na.action=na.omit),
                               scope= formula(Cluster.cca),
                               direction="both",pstep=10000,nperm=999999,  Pin = 0.05, Pout = 0.1))



####################################################################################
## the best terms for the final model are assessed
####################################################################################
rownames(Cluster.cca.forward$CCA$biplot)
rownames(Cluster.cca.backward$CCA$biplot)
rownames(Cluster.cca.both$CCA$biplot)
(summarized.selection<- c(rownames(Cluster.cca.forward$CCA$biplot),
                          rownames(Cluster.cca.backward$CCA$biplot),
                          rownames(Cluster.cca.both$CCA$biplot)))

(sortsumsel<- (sort(table(summarized.selection))))
(constrains.freq.table<- as.data.frame (sortsumsel))

## to what ever reason you have to use one of the two FORMULA functions... seems to be a R bug.... just uncomment and run the other if you get an error
## the same is true for the removement of single terms further below

#(FORMULA <-as.formula(paste("ClusterHel ~ ", paste  (rownames((constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=1))])), # here take again complete data set
#                                                     collapse="+"),sep = "")))

(FORMULA <-as.formula(paste("ClusterHel ~ ", paste ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=2))]), # here take again complete data set
                                                    collapse="+"),sep = "")))


###########################################################################################################
#
# Run selected CCA model and tests and apriori fit physico-chemical variables not included in
# the model
#
###########################################################################################################

(Cluster.cca.forwsel <- cca(FORMULA, PhysicoChemicalLog, na.action=na.exclude))
summary(Cluster.cca.forwsel)
anova(Cluster.cca.forwsel,by="terms", permutations = how(nperm=999),
      model = c("direct"), 
      parallel = getOption("mc.cores"), strata = NULL,
      scope = NULL)
# test the variation inflation factors and significance of models, axes and constraints




## function to replace variables from the rda constraints due to colinearity
mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}


#(FORMULA    <-as.formula(paste("ClusterHel ~ ", mgsub(c("\\alpha_.", "\\+pmoA_.", "\\+Depth_DGT_m", "\\+Depth_m"  ),
#                                                      c("",""), paste (rownames (constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=2))]), # here take again complete data set
#                                                     collapse="+")),sep = "")))

#(FORMULA    <-as.formula(paste("ClusterHel ~ ", mgsub(c("\\+DIC_mM"),  c(""),
#                                                      paste ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=1))][]), # here take again complete data set, change [] to -1 if first parameter should be replaced
#                                                             collapse="+")), sep = "")))
(Cluster.cca.forwsel <- cca(FORMULA, PhysicoChemicalLog, na.action=na.exclude))
sort(vif.cca(Cluster.cca.forwsel))


#FORMULA= as.formula(ClusterHel ~ CH4_uM + SO42._uM + DIC_mM + Light_uE..m2s)    #use this one if you want to make a custom selection of variables
#FORMULA = ClusterHel ~ Chla_ug.L + Cond_uS.cm + Zn_diss_nM + CH4_uM

## rerun cca with new formula
(Cluster.cca.forwsel <- cca(FORMULA, PhysicoChemicalLog, na.action=na.exclude))
sort(vif.cca(Cluster.cca.forwsel))

##################----------------------------------------------------------------------------------------------
(Cluster.cca.best.forwsel <- vegan::cca(FORMULA, PhysicoChemicalLog, na.action=na.exclude))   # this is the final model
##################-----------------------------------------------------------------------------------------------
#Cluster.cca.best.forwsel = update(Cluster.cca.best.forwsel, .~. + CH4_uM - Cu_diss_nM)
# (Cluster.cca.best.forwsel <- cca(ClusterHel ~ Cu_DGT_nM + Condition(Depth_m), PhysicoChemicalLog, na.action=na.exclude))   # alternatively make model with specific variables


summary(Cluster.cca.best.forwsel)

anova(Cluster.cca.best.forwsel) # check if total model is singificant


# check for significance of environmental terms (all na data has to be removed!)
#(FORMULAsub    <-as.formula(paste("ClusterHelsub ~ ", mgsub(c("\\+Fe_tot_nM", "\\+Depth_DGT_m", "\\+Light_uE..m2s"  ),
#                                                            c("","",""), paste (rownames (constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=3))]), # here take again complete data set
#                                                                              collapse="+")),sep = "")))

#(FORMULAsub    <-as.formula(paste("ClusterHelsub ~ ", mgsub(c("\\+pmoA_.", "\\+Depth_DGT_m", "\\+Depth_m"  ), c("","",""),
#                                                            paste ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=2))][-1]), # here take again complete data set
#                                                                               collapse="+")),sep = "")))

FORMULAsub=FORMULA

(Cluster.cca.best.forwsel.sub <- cca(FORMULAsub, PhysicoChemicalLogsub, na.action=na.omit))  
anova(Cluster.cca.best.forwsel.sub, by="terms", permutations = how(nperm=999),
      model = c("direct"), 
      parallel = getOption("mc.cores"), strata = NULL,
      scope = NULL)

anova(Cluster.cca.best.forwsel.sub, by="axis", permutations = how(nperm=999),  # check for significant axes
      model = c("direct"),
      strata = NULL,
      cutoff = 1)

# fit the non-used environmental variables in the RDA model
(best.model.env.selection <-c(rownames(Cluster.cca.best.forwsel$CCA$biplot)))
(rest.model.env.selection <-PhysicoChemicalRaw[! c(colnames(PhysicoChemicalRaw)) %in% best.model.env.selection])

(envfit.rest  <- envfit(Cluster.cca.best.forwsel, rest.model.env.selection[,1:length(rest.model.env.selection),drop=FALSE], na.rm=TRUE))

## fit the species in the rda model
#(speciesfit   <- envfit(Cluster.cca.best.forwsel, rotsee.vegan.otu.table))

PerContCon(Cluster.cca.best.forwsel) 

(envfit.apriori  <- envfit(Cluster.cca.best.forwsel, PhysicoChemicalRawsub, na.rm=TRUE, choices=c(1,2))) # you can check with the above output
## which dimensions might be good for a specific model parameter for
## fitting on and plotting of the CCA constraints

##########################################################################################################
# Bi-plot of the CCA model         
##########################################################################################################
my.color1<-c("cyan", "cyan4", "cyan3")


#pdf("CCA pmoA All Seasons.pdf",useDingbats=FALSE)

## Basic Plot
ordiplot<-plot(Cluster.cca.best.forwsel,type="none",choices = c(xx<-1, yy<-2),scaling=3,main="", xlab="",ylab="")     # here you can change the axes to be ploted       


points(ordiplot,"species", cex=normalized.species.size.vec * 2 ,#cex=1.5,
       pch=21, col=col.genus[genus.df$id], bg=alpha(col.genus[genus.df$id], 0.5), lwd=0)
legend("topleft",legend=sort(unique(genus.df[,"Genus"])),pch=c(19),col=alpha(col.genus.orig, 0.99), bty="n")
points(ordiplot,"sites", cex=3, alpha= 2, pch=as.numeric(as.factor(rotsee.env.table$Campaign))+21, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.5))

text(ordiplot$sites[, 1]*1.3, ordiplot$sites[, 2]*1.3,
     labels=unlist(dimnames(ordiplot$sites)[1]),
     col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],cex=0.7)

image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.3,0.7), axis.args = list(cex.axis = 0.7),
           zlim=range(-(rotsee.env.table$Depth_m)) )



## Confidence ellipses for the Locations (SE +-95%)
with(rotsee.env.table, ordiellipse(ordiplot, groups= interaction(Oxidation_Zone), display="sites",
                                   draw="polygon", col=c( "darkblue","green","lightblue","orange"), border= "darkgrey", alpha=20 , label=FALSE))

#(groupz<-sort(unique(rotsee.env.table$Light_Zone)))
#with(rotsee.env.table,ordiellipse(ordiplot, Light_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[1], show.groups=groupz[1],draw="polygon", alpha=30))
#with(rotsee.env.table,ordiellipse(ordiplot, Light_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[2], show.groups=groupz[2],draw="polygon", alpha=30))


## plot the significant environmental variables from the selected canonical model
(curr.plot.multiplier <- vegan:::ordiArrowMul(scores(Cluster.cca.best.forwsel$CCA$biplot[,c(xx,yy)]),
                                              at = c(0,0), fill = 1,display="sites", choices = c(1,2)))# multiplier to fit the current plot ratio
arrows(0, 0, Cluster.cca.best.forwsel$CCA$biplot[, 1]*curr.plot.multiplier*0.7, Cluster.cca.best.forwsel$CCA$biplot[, 2]*curr.plot.multiplier*0.7, lty=1,lwd=1.5,
       length = 0.15, angle = 30,  col="darkgrey")


text(Cluster.cca.best.forwsel$CCA$biplot[, 1]*curr.plot.multiplier*0.9, Cluster.cca.best.forwsel$CCA$biplot[, 2]*curr.plot.multiplier*0.9,
     labels=unlist(dimnames(Cluster.cca.best.forwsel$CCA$biplot)[1]), col="black", cex=0.8)



title(xlab=as.expression(paste("CCA",xx," (",round(c(100/(Cluster.cca.best.forwsel$tot.chi)*(Cluster.cca.best.forwsel$CCA$eig[xx])),2),"%)")),
      ylab=as.expression(paste("CCA",yy," (",round(c(100/(Cluster.cca.best.forwsel$tot.chi)*(Cluster.cca.best.forwsel$CCA$eig[yy])),2),"%)")))


dev.off()

# Optionally: plot the rest of the environmental variables if you wish to....
pval.limit <- which(c(envfit.rest$vectors$pvals)<=1) # set max p-value 
(envfit.vectors       <- scores(envfit.rest,"vectors")[c(pval.limit),1:2])  # these are envfit vectors scaled to each one another (i.e. multipied by square root of column r2)
(env.loading.factor   <- sqrt((envfit.rest$vectors$r)[c(pval.limit)]))      # loading factor for the envfit scores to scale them according to their correlations
(curr.plot.multiplier <- vegan:::ordiArrowMul(scores(envfit.rest,"vectors")[c(pval.limit),1:2],
                                              at = c(0,0), fill = 1,display="sites", choices = c(xx,yy)))# multiplier to fit the current plot ratio
env.arrows.x  <- as.data.frame(scores(envfit.rest,"vectors"))[c(pval.limit),1]
env.arrows.y  <- as.data.frame(scores(envfit.rest,"vectors"))[c(pval.limit),2]
env.names     <-rownames(as.data.frame(scores(envfit.rest,"vectors")))[c(pval.limit)]

arrows(0, 0, env.arrows.x*curr.plot.multiplier*0.55, env.arrows.y*curr.plot.multiplier*0.55, lty=1, col="lightgrey",
       length = 0.15, angle = 30)
text(env.arrows.x*curr.plot.multiplier*0.72, env.arrows.y*curr.plot.multiplier*0.72,
     labels=env.names, col="lightgrey", cex=0.7)    

## perconcont function 
PerContCon<- function(RDAmodel) {                # RDAmodel is a vegan ordination output)     
  
  (coefRDA<-intersetcor(RDAmodel))
  writeLines("Weighted Canonical Coefficients\n")
  print(coefRDA)
  writeLines("\n")
  writeLines("-----------------------------------------------------------------------------------------\n")
  
  envsel<-rownames(attributes(RDAmodel$terminfo$terms)$factors)
  rdasel<-unlist(attributes(RDAmodel$CCA$eig))
  mylist1<-vector("list",length(rdasel))
  mylist2<-vector("list",length(rdasel))
  mylist3<-vector("list",length(envsel))
  mylist4<-vector("list",length(envsel))
  varexp<- 100/RDAmodel$tot.chi *RDAmodel$CCA$tot.chi 
  
  for(i in 1:length(rdasel)){
    
    contri1<-100/RDAmodel$CCA$tot.chi * RDAmodel$CCA$eig[i] / (sum(abs(coefRDA[,i]))) * coefRDA[,i]
    mylist1[i]<-list(contri1)
    
    contri2<-100/RDAmodel$tot.chi * RDAmodel$CCA$eig[i] / (sum(abs(coefRDA[,i]))) * coefRDA[,i]
    mylist2[i]<-list(contri2)
    
    
    
    unname(contri3<-100/(sum(abs(coefRDA[,]))) * (sum(abs(coefRDA[i,]))))
    mylist3[i]<-list(contri3)
    
    unname(contri4<-100/RDAmodel$tot.chi * sum(RDAmodel$CCA$eig[]) / (sum(abs(coefRDA[,]))) * (sum(abs(coefRDA[i,]))))
    mylist4[i]<-list(contri4)
    
    
    
    
  }
  
  # constraint realm output
  #------------------------
  VariationList1<-data.frame(matrix(unlist(mylist1), nrow=length(rdasel), byrow=T))
  #VariationList<-Reduce(rbind, mylist))
  rownames(VariationList1)<-unlist(attributes(RDAmodel$CCA$eig))
  colnames(VariationList1)<-rownames(attributes(RDAmodel$terminfo$terms)$factors)
  totalExplained1<-rowSums(abs(VariationList1))
  
  writeLines("% explained by environmental variables (within constraint realm )\n")
  print(VariationList1)
  writeLines("\n")
  writeLines("total % explained by all environmental variables (within constraint realm)\n")
  print(totalExplained1)
  writeLines("\n")
  writeLines("\n")
  
  VariationList3<-data.frame(matrix(unlist(mylist3), nrow=length(rdasel), byrow=T))
  #VariationList2<-data.frame(Reduce(rbind, mylist2))
  rownames(VariationList3)<-unlist(rownames(attributes(RDAmodel$terminfo$terms)$factors))
  colnames(VariationList3)<-"% total explained"
  VariationList3.splitbyaxes<- sum(VariationList3[,])/totalExplained1*VariationList1
  
  writeLines("total % explained by environmental variable (within constraint realm)\n")
  print(VariationList3)
  writeLines("\n")
  writeLines("% explained by environmental variable  split and scaled to 100% by single canonical axes(within constraint realm)\n")
  print(t(VariationList3.splitbyaxes))
  writeLines("\n")
  writeLines("\n")
  writeLines("-----------------------------------------------------------------------------------------\n")
  
  #total output
  #------------------------
  VariationList2<-data.frame(matrix(unlist(mylist2), nrow=length(rdasel), byrow=T))
  #VariationList2<-data.frame(Reduce(rbind, mylist2))
  rownames(VariationList2)<-unlist(attributes(RDAmodel$CCA$eig))
  colnames(VariationList2)<-rownames(attributes(RDAmodel$terminfo$terms)$factors)
  totalExplained2<-rowSums(abs(VariationList2))
  
  writeLines("% explained by environmental variables (percentages on specific axes)\n")
  print(VariationList2)
  writeLines("\n")
  writeLines("total % explained by all environmental variables\n")
  print(totalExplained2)
  writeLines("\n")
  writeLines("\n")
  
  VariationList4<-data.frame(matrix(unlist(mylist4), nrow=length(rdasel), byrow=T))
  #VariationList2<-data.frame(Reduce(rbind, mylist2))
  rownames(VariationList4)<-unlist(rownames(attributes(RDAmodel$terminfo$terms)$factors))
  colnames(VariationList4)<-"% total explained"
  totalExplained3<-rowSums(abs(VariationList4))
  VariationList4.splitbyaxes<-sum(VariationList4[,])/totalExplained2*VariationList2
  
  writeLines("total % explained by environmental variables\n")
  print(VariationList4)
  writeLines("\n")
  writeLines(paste("% explained by environmental variable  split and scaled to single canonical axes:\n",
                   round(sum(VariationList4[,]),4),"% of variation are covered\n"))
  print(t(VariationList4.splitbyaxes))
  writeLines("\n")
  writeLines("\n")
  writeLines("-----------------------------------------------------------------------------------------\n")
  
  #total model power
  #-------------------------
  writeLines("% explained variation in complete model\n")
  print(varexp)
  writeLines("\n")
  
  
  
}



################################################################
## End of the script
################################################################