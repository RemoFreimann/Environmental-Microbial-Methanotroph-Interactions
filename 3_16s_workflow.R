## =======================================================================================

## Bacteria Community Analysis (16S rRNA gene based) in R after QIIME Pipeline 
## last tested in R version 4.0.0 (2020-04-24) -- "Arbor Day"

## =======================================================================================

## ---------------------------------------------------------------------------------------
# unlock lines constraint for working folder if necessary 
## ---------------------------------------------------------------------------------------
.trPaths <- paste(paste(Sys.getenv('APPDATA'), '\\Tinn-R\\tmp\\',
            sep=''), c('', 'search.txt', 'objects.txt', 'file.r', 'selection.r',
            'block.r', 'lines.r'), sep='')
## ---------------------------------------------------------------------------------------
# install potentially missing packages
## ---------------------------------------------------------------------------------------
install.packages("devtools")
library(devtools)
list.of.packages <- c("vegan","phyloseq", "reshape2", "ggplot2", "tcltk", "fields", "RColorBrewer", "mvabund",
                      "doParallel","foreach", "pheatmap", "grDevices", "relaimpo", "multcomp", "car", "MASS",
                      "leaps","RcmdrMisc", "gplots","cooccur", "ggtree", "missMDA", "corrplot")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=TRUE)
install_github("cpauvert/psadd")
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(c("multtest", "phyloseq"))

BiocManager::install('phyloseq')

## ---------------------------------------------------------------------------------------
# load packages and set working directory (some might be redundant, check also for error messages)
## ---------------------------------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(phyloseq)
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
library(biomformat)
library(missMDA)
library(corrplot)
library(psadd)
library(ggtree)


## ---------------------------------------------------------------------------------------
## set bw theme
## ---------------------------------------------------------------------------------------
theme_set(theme_bw()) 

## ---------------------------------------------------------------------------------------
## Data Import in R
## ---------------------------------------------------------------------------------------

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

## check imported data
sample_data(bac.rotsee)                   # sample variables
head(get_variable(bac.rotsee))            # variable overview
rank_names(bac.rotsee)                    # Available Taxonomic ranks: [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
sort(sample_sums(bac.rotsee))             #Reads per sample
taxa_names(bac.rotsee)
ntaxa(bac.rotsee)
tax_table(bac.rotsee)


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
 library(reshape2)
 library(phyloseq)
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
 library(biomformat)
 library(missMDA)
 library(corrplot)
 library(psadd)
 #library(ggtree)           ## interferes with "phylo" class, so dont load it anymore as its not needed below
 
 
 ## merge the 16S rRNA data set
 (bac.rotsee = merge_phyloseq(rotsee.otu.sort, rotsee.tax.sort, rotsee.env.sort, rotsee.tree.sort, rotsee.refseq.sort ))

#--------------------------------------------------------------------------------------------------------------------
 rotsee.env.table <- sample_data(bac.rotsee) 
## set negative values in the metafile to 0, except for 13C-CH4 isotopes
 rotsee.env.table[,-which(names(rotsee.env.table)=="X13CH4")][rotsee.env.table[,-which(names(rotsee.env.table)=="X13CH4")] < 0] <- 0
 #-------------------------------------------------------------------------------------------------------------------- 

## ---------------------------------------------------------------------------------------
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

## Filter OTUs that appear at least in three different samples
filt3 <- genefilter_sample(bac.rotsee.prune, filterfun_sample(function(x) x > 3), A=1)
bac.rotsee.prune.filt3 <- prune_taxa(filt3, bac.rotsee.prune)
ntaxa(bac.rotsee.prune.filt3) 

## here we quickly check sparsity for spearman/pearson network based analysis (see other scripts later on)
(bac.sparse<-plot_sparsity(bac.rotsee.prune.filt3))
taxa_prev(bac.rotsee.prune.filt3)
halfsamples= ncol(sample_data(bac.rotsee.prune.filt3))*0.5
sparsekeep=names(which(taxa_prev(bac.rotsee.prune.filt3)>=halfsamples))
bac.rotsee.spearman.network.phyloseq =  prune_taxa(sparsekeep, bac.rotsee)

## ---------------------------------------------------------------------------------------
## build phyloseq object for MOBs (based on 16S rRNA indetification) only and merge with prefitlered total community data set
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


## produce 16s rRNA data excluding MOBs to be used in combination with the pMOA data (this is the data that is used in the Network analysis (DIABLO))
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

(rotsee.MOB.otu.table <-otu_table(bac.rotsee.MOBs))
(rotsee.MOB.tax.table <- tax_table(bac.rotsee.MOBs))
(rotsee.MOB.tree.table <- phy_tree(bac.rotsee.MOBs))
(rotsee.MOB.env.table <- sample_data(bac.rotsee.MOBs))
(rotsee.MOB.refseq <- refseq(bac.rotsee.MOBs))

bac.rotsee.and.MOBs.otu <- merge_phyloseq(rotsee.otu.table, rotsee.MOB.otu.table )
bac.rotsee.and.MOBs.tax <- merge_phyloseq(rotsee.tax.table, rotsee.MOB.tax.table )
bac.rotsee.and.MOBs.env <- merge_phyloseq(rotsee.env.table, rotsee.MOB.env.table )
bac.rotsee.and.MOBs.refseq <- merge_phyloseq(rotsee.refseq, rotsee.MOB.refseq )



(bac.rotsee.no.MOBs =  merge_phyloseq(rotsee.otu.table, rotsee.tax.table, rotsee.tree.table , rotsee.env.table,  rotsee.refseq) )
(bac.rotsee.MOBs =  merge_phyloseq(rotsee.MOB.otu.table, rotsee.MOB.tax.table, rotsee.MOB.tree.table , rotsee.MOB.env.table,  rotsee.MOB.refseq) )
(bac.rotsee.and.MOBs = merge_phyloseq(bac.rotsee.and.MOBs.otu, bac.rotsee.and.MOBs.tax, bac.rotsee.and.MOBs.env , bac.rotsee.and.MOBs.refseq) )


## CHOOSE ONE OF THE THREE DATASETS FOR FURTHER PROCESSING ACCORDING TO WHAT YOU WANT TO DO!
#(bac.rotsee.prune.filt3 <- bac.rotsee.and.MOBs)   
#(bac.rotsee.prune.filt3 <- bac.rotsee.no.MOBs)     # is used for the procruste NMDS rotation later on and the network analysis, so no redundancy of 16S rRNA gene based MOBs will occure
#(bac.rotsee.prune.filt3 <- bac.rotsee.MOBs) 

## rarefy curve to check your library coverage by samples
 rarefyOTU <- t(otu_table(bac.rotsee.and.MOBs.otu))
 
 S <- specnumber(rarefyOTU) # observed number of species
 (raremax <- min(rowSums(rarefyOTU)))
 Srare <- rarefy(rarefyOTU, raremax)
 plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
 abline(0, 1)
 rarecurve(rarefyOTU, step = 20, sample = raremax, col = "blue", cex = 0.6)
 
## ***************************************************************************************
## ---------------------------------------------------------------------------------------
##
## Data prefiltering (standardization)
##
## ---------------------------------------------------------------------------------------
## ***************************************************************************************

## ---------------------------------------------------------------------------------------
## check reads distribution
## ---------------------------------------------------------------------------------------

# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(bac.rotsee.and.MOBs))

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "blue", binwidth = 2000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

## ---------------------------------------------------------------------------------------
## Standardization of sample reads to the median sequencing depth applied to each sample  (readcounds, 100% and 1)
## ---------------------------------------------------------------------------------------

total = median(sample_sums(bac.rotsee.and.MOBs))
standf = function(x, t=total) (t * (x / sum(x)))
bac.rotsee.and.MOBs.seqdepth = transform_sample_counts(bac.rotsee.and.MOBs, standf)
standf = function(x) (100 * x / sum(x))
bac.rotsee.and.MOBs.seqdepth_ab_100 = transform_sample_counts(bac.rotsee.and.MOBs.seqdepth, standf)
bac.rotsee.and.MOBs.seqdepth.ab <- transform_sample_counts(bac.rotsee.and.MOBs.seqdepth, function(x){x/sum(x)})

total = median(sample_sums(bac.rotsee.no.MOBs)) 
standf = function(x, t=total) (t * (x / sum(x)))
bac.rotsee.no.MOBs.seqdepth = transform_sample_counts(bac.rotsee.no.MOBs, standf)  
standf = function(x) (100 * x / sum(x))
bac.rotsee.no.MOBs.seqdepth_ab_100 = transform_sample_counts(bac.rotsee.no.MOBs.seqdepth, standf)
bac.rotsee.no.MOBs.seqdepth.ab <- transform_sample_counts(bac.rotsee.no.MOBs.seqdepth, function(x){x/sum(x)})

total = median(sample_sums(bac.rotsee.MOBs))
standf = function(x, t=total) (t * (x / sum(x)))
bac.rotsee.MOBs.seqdepth = transform_sample_counts(bac.rotsee.MOBs, standf) 
standf = function(x) (100 * x / sum(x))
bac.rotsee.MOBs_seqdepth_ab_100 = transform_sample_counts(bac.rotsee.MOBs.seqdepth, standf) 
bac.rotsee.MOBs.seqdepth.ab <- transform_sample_counts(bac.rotsee.MOBs.seqdepth, function(x){x/sum(x)})

## quickly check colsums (i.e. reads)
colSums(otu_table(bac.rotsee.no.MOBs.seqdepth))
colSums(otu_table(bac.rotsee.no.MOBs.seqdepth_ab_100))
colSums(otu_table(bac.rotsee.no.MOBs.seqdepth.ab))

## ---------------------------------------------------------------------------------------
## Rarefy reads for analysis like UNIFRAC
## ---------------------------------------------------------------------------------------
bac.rotsee.and.MOBs.rare <- rarefy_even_depth(bac.rotsee.and.MOBs, sample.size = min(sample_sums(bac.rotsee.and.MOBs)),
                                                rngseed = 33, replace = FALSE, trimOTUs = TRUE)
bac.rotsee.and.MOBs.rare
sample_sums(bac.rotsee.and.MOBs.rare) 


## ***************************************************************************************
## ---------------------------------------------------------------------------------------
## choose one of the  prefiltered data set for further investigation
## ---------------------------------------------------------------------------------------
## the herein choosen prefilterd data will be used for the rest of the script!!!! 
## ***************************************************************************************

## remove non-existent OTUs from the newly generated subsets
bac.rotsee.all <-  bac.rotsee.and.MOBs.seqdepth         
ntaxa(bac.rotsee.all)  
bac.rotsee.all <- prune_taxa(taxa_sums(bac.rotsee.all) > 0, bac.rotsee.all)
ntaxa(bac.rotsee.all) 

bac.rotsee.nomobs <-  bac.rotsee.no.MOBs.seqdepth      
ntaxa(bac.rotsee.nomobs)  
bac.rotsee.nomobs <- prune_taxa(taxa_sums(bac.rotsee.nomobs) > 0, bac.rotsee.nomobs)
ntaxa(bac.rotsee.nomobs) 

bac.rotsee.mobs <-  bac.rotsee.MOBs.seqdepth
ntaxa(bac.rotsee.mobs)
bac.rotsee.mobs <- prune_taxa(taxa_sums(bac.rotsee.mobs) > 0, bac.rotsee.mobs)
ntaxa(bac.rotsee.mobs) 

## ---------------------------------------------------------------------------------------
## again, choose the according dataset for further analysis (i.e. including or excluding MOBs etc.)
## --------------------------------------------------------------------------------------

## check sparsity again
(bac.sparse<-plot_sparsity(bac.rotsee.nomobs))
taxa_prev(bac.rotsee.nomobs)
halfsamples= ncol(sample_data(bac.rotsee.nomobs))*0.5
sparsekeep=names(which(taxa_prev(bac.rotsee.nomobs)>=halfsamples))

bac.rotsee.nomobs.spearman.network.phyloseq =  prune_taxa(sparsekeep, bac.rotsee.nomobs)      ## you might want to use this subset for the pearson/spearman network analysis later on


## ---------------------------------------------------------------------------------------
## Barplots for main phyla
## ---------------------------------------------------------------------------------------



## FUNCTIONS and DEFINITIONS for the plotting in ggplot2:

## Defining a theme (graphic settings) for barplots
theme_HB <- function(base_size = 12, base_family = "") {
  # Starts with theme_bw and then modify some parts
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text         = element_text(size = rel(1.3), color= "black"),
      axis.title        = element_text(size = rel(1.5)),
      legend.text       = element_text(size = rel(1.0)),
      legend.title      = element_text(face = "bold", size = rel(1.0)),
      plot.title        = element_text(face = "bold", size = rel(1.8), hjust=0.5)
    )
}

##Function for ordered barplots (often easier to read)
##From the web: https://github.com/joey711/phyloseq/issues/442
##style- adjusted by Helmut
plot_ordered_bar<-function (physeq, x = "Sample", 
                            y = "Abundance", 
                            fill = NULL, 
                            leg_size = 0.5,
                            title = NULL) {
  require(ggplot2)
  require(phyloseq)
  require(plyr)
  require(grid)
  bb <- psmelt(physeq)
  
  samp_names <- aggregate(bb$Abundance, by=list(bb$Sample), FUN=sum)[,1]
  .e <- environment()
  bb[,fill]<- factor(bb[,fill], rev(sort(unique(bb[,fill])))) #fill to genus
  
  
  bb<- bb[order(bb[,fill]),] # genus to fill
  p = ggplot(bb, aes_string(x = x, y = y, 
                            fill = fill), 
             environment = .e, ordered = FALSE)
  
  
  p = p +geom_bar(stat = "identity", 
                  position = "stack", 
                  color = "black") 
  
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0)) 
  
  p = p + guides(fill = guide_legend(override.aes = list(colour = NULL), reverse=TRUE)) + 
    theme(legend.key = element_rect(colour = "black")) 
  
  p = p + theme(legend.key.size = unit(leg_size, "cm"))

  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}




## plot of phyla when at least 1% abundant in any sample (raw reads, not standardized)

bac.rotsee.test = filter_taxa(bac.rotsee.no.MOBs, function(x) sum(x > sample_sums(bac.rotsee.no.MOBs)/100*1) > (0.1*length(x)), TRUE) # 1% of all reads and at least in 10% of the samples
standf = function(x) (100 * x / sum(x))
bac.rotsee.test_ab_100 = transform_sample_counts(bac.rotsee.test, standf) 
#bac.rotsee.test_ab_100=  tax_glom(bac.rotsee.test_ab_100, "Phylum")
Top12 = names(sort(taxa_sums(bac.rotsee.test_ab_100), TRUE)[])  


#pdf("main phylum 16s_relativ.pdf" ,useDingbats=FALSE)  
title="relative abundance of main phylum"
p = plot_ordered_bar(bac.rotsee.test_ab_100, x="Depth_m", fill = "Phylum", title=title)
p + coord_flip()  +  
  scale_x_reverse(breaks=c(1:16), limits=c(15,1)) +
  ylim(0,101) +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_line(color="grey"),
        #panel.grid.minor = element_line(color="grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        plot.title=element_text(size=15)) + 
  facet_wrap(~Campaign)
#dev.off()






########################################################################################################################################
# this lines are used to run the rest of the manuscript for different subsets (i.e. Campaign or a specific depth zone or all together)
########################################################################################################################################

## choose one of the lines

 bac.rotsee <- bac.rotsee.all
 bac.rotsee.16s <- bac.rotsee.nomobs      # All Seasons, reload these two lines before switching to a specific campaign
 bac.rotsee.mobs
  
# 
# bac.rotsee <- subset_samples(bac.rotsee, Campaign=="A_June_2013")
# bac.rotsee.mobs <- subset_samples(bac.rotsee.mobs, Campaign=="A_June_2013")
# 
# bac.rotsee <- subset_samples(bac.rotsee, Campaign=="B_August_2013")
# bac.rotsee.mobs <- subset_samples(bac.rotsee.mobs, Campaign=="B_August_2013")
# 
# bac.rotsee <- subset_samples(bac.rotsee, Campaign=="C_September_2014")
# bac.rotsee.mobs <- subset_samples(bac.rotsee.mobs, Campaign=="C_September_2014")
# 
# bac.rotsee <- subset_samples(bac.rotsee, Campaign=="D_September_2015")
# bac.rotsee.mobs <- subset_samples(bac.rotsee.mobs, Campaign=="D_September_2015")
# 
# bac.rotsee <- subset_samples(bac.rotsee, Campaign!="D_June_2015")
# bac.rotsee.mobs <- subset_samples(bac.rotsee.mobs, Campaign!="D_June_2015")# without June



## remove non-existent OTUs from the subset
# ntaxa(bac.rotsee) 
# bac.rotsee <- prune_taxa(taxa_sums(bac.rotsee) > 0, bac.rotsee)
# ntaxa(bac.rotsee) 
# ntaxa(bac.rotsee.mobs) 
# bac.rotsee.mobs <- prune_taxa(taxa_sums(bac.rotsee.mobs) > 0, bac.rotsee.mobs)
# (nmobs<-ntaxa(bac.rotsee.mobs) )

########################################################################################################################################
# produce OTU, taxanomic and environmental tables for the usage with other packages (i.e. Vegan functionalities)
########################################################################################################################################

options(max.print = 2000000000)

 ## All 16S rRNA OTUs tables
(rotsee.otu.table <- otu_table(bac.rotsee.all))
(rotsee.tax.table <- tax_table(bac.rotsee.all))
(rotsee.env.table <- get_variable(bac.rotsee.all)) #mapfile 
ntaxa(bac.rotsee)
# set negative values in the metafile to 0, except for 13C-CH4 isotopes
rotsee.env.table[,-which(names(rotsee.env.table)=="X13CH4")][rotsee.env.table[,-which(names(rotsee.env.table)=="X13CH4")] < 0] <- 0 
names(rotsee.env.table)
## make the sample IDs the column names
rotsee.env.table=rotsee.env.table[,-1] ## remove first column
rotsee.env.table.log = cbind( rotsee.env.table[,c(1:13)],
                              log1p(rotsee.env.table[,c(14),drop=FALSE]),
                              rotsee.env.table[,c(15:17)],
                              log1p(rotsee.env.table[,c(18:34)]),
                              #log1p(abs(rotsee.env.table[,c(35),drop=FALSE])),
                              log1p(rotsee.env.table[,c(35), drop=FALSE] + (1+abs(min(rotsee.env.table$X13CH4)))),
                              log1p(rotsee.env.table[,c(36),drop=FALSE]),
                              rotsee.env.table[,c(37),drop=FALSE],
                              log1p(rotsee.env.table[,c(38:70)]))


## define vectors for smmo and pmmo mobs for plot coloring further below
(MOBs_reordered<- MOBs)
(MOBs_reordered<- na.omit(match(MOBs, c(rownames(rotsee.otu.table@.Data)))))
(mobs_vector_all_otus = rownames(rotsee.otu.table@.Data)[MOBs_reordered])
sort(match(smmo_vector,mobs_vector_all_otus))
nmobs
## choose the vector acording to the month
(pMMO_vector = c(1,rep(2,11),1,1,1,1,rep(2,5)))                         # All 
# (pMMO_vector = c(1,rep(2,11),1,1,1,1,rep(2,5)))                         # Without June
# (pMMO_vector = c(1,rep(2,10),1,1,rep(2,6)))                             # June      2013
# (pMMO_vector = c(1,rep(2,9),1,rep(2,3)))                                # August    2014
# (pMMO_vector = c(1,rep(2,11),1,1,1,rep(2,4)))                           # September 2014
#(pMMO_vector = c(1,rep(2,9),1,1,1,1,1,1,rep(2,4)))                       # September 2015



## only 16S rRNA MOBs tables
(rotsee.mobs.otu.table <- otu_table(bac.rotsee.mobs))
(rotsee.mobs.tax.table <- tax_table(bac.rotsee.mobs))
(rotsee.mobs.env.table <- get_variable(bac.rotsee.mobs))
# set negative values in the metafile to 0, except for 13C-CH4 isotopes
rotsee.mobs.env.table[,-which(names(rotsee.mobs.env.table)=="X13CH4")][rotsee.mobs.env.table[,-which(names(rotsee.mobs.env.table)=="X13CH4")] < 0] <- 0 
rotsee.mobs.env.table.log=rotsee.env.table.log


(mobs_vector_all_otus = rownames(rotsee.mobs.otu.table@.Data))
sort(match(smmo_vector,rownames(rotsee.mobs.otu.table@.Data)))



## choose the vector acording to the month

(pMMO_vector_mobs = c(1,rep(2,11),1,1,1,1,rep(2,5)))                         # All
# (pMMO_vector_mobs = c(1,rep(2,11),1,1,1,1,rep(2,5)))                         # without June
# (pMMO_vector_mobs = c(1,rep(2,10),1,1,rep(2,6)))                             # June 2013
# (pMMO_vector_mobs = c(1,rep(2,9),1,rep(2,3)))                                # August 2013
# (pMMO_vector_mobs = c(1,rep(2,11),1,1,1,rep(2,4)))                           # September 2014
#(pMMO_vector_mobs = c(rep(2,1),1,1,1,1,rep(2,8),1,1,1,rep(2,4)))              # September 2015




## only 16S rRNA without MOBs
(rotsee.nomobs.otu.table <- otu_table(bac.rotsee.nomobs))
(rotsee.nomobs.tax.table <- tax_table(bac.rotsee.nomobs))
(rotsee.nomobs.env.table <- get_variable(bac.rotsee.nomobs))
# set negative values in the metafile to 0, except for 13C-CH4 isotopes
rotsee.nomobs.env.table[,-which(names(rotsee.nomobs.env.table)=="X13CH4")][rotsee.nomobs.env.table[,-which(names(rotsee.nomobs.env.table)=="X13CH4")] < 0] <- 0 
rotsee.nomobs.env.table.log=rotsee.env.table.log




(mobs_vector_all_otus = rownames(rotsee.mobs.otu.table@.Data))
sort(match(smmo_vector,rownames(rotsee.mobs.otu.table@.Data)))


## you can view the tables in a new window:
#            utils::View(rotsee.tax.table)

## you can write the tables into csv files for further analysis in excel:
#            write.csv(rotsee.otu.table, file ="rotsee_otu_table.csv")




####################################################################################################################################
## choose the env table with only the variables acutally measured during a specific campaign (only do this if you do not
## want to use the imputed variable table)

## run this four lines and choose the campaign specific lines below
# (rotsee.env.table <- get_variable(bac.rotsee))
# names(rotsee.env.table)
# rotsee.env.table=rotsee.env.table[,-1] ## remove first column
# 
# (rotsee.mobs.env.table <- get_variable(bac.rotsee.mobs))
# names(rotsee.mobs.env.table)
# rotsee.env.table=rotsee.env.table[,-1] ## remove first column
# 
# rotsee.env.table[,-which(names(rotsee.env.table)=="X13CH4")][rotsee.env.table[,-which(names(rotsee.env.table)=="X13CH4")] < 0] <- 0
# rotsee.mobs.env.table[,-which(names(rotsee.mobs.env.table)=="X13CH4")][rotsee.mobs.env.table[,-which(names(rotsee.mobs.env.table)=="X13CH4")] < 0] <- 0
# (rotsee.env.table.sort=rotsee.env.table)
# (rotsee.mobs.env.table.sort=rotsee.mobs.env.table)
# ####
# 
# names(rotsee.env.table.sort)
#  (rotsee.env.table=cbind( rotsee.env.table.sort[,c(1:13,15,16,17,37)],(rotsee.env.table.sort[,c(14,18:36,38:58)])))    # All
#  rotsee.env.table.log=cbind( rotsee.env.table.sort[,c(1:13,15,16,17,37)],log1p(rotsee.env.table.sort[,c(14,18:36,38:58)]))
#  envvector<-c(11,12,14,15,16,17,18,19,20,21,
#               22,23,24,25,26,27,29,32,34,35,37,38,39,40,
#               48,49:56,57,58)
#  head(rotsee.env.table[,envvector])
#  rotsee.mobs.env.table=rotsee.env.table
#  rotsee.mobs.env.table.log=rotsee.env.table.log
# 
# 
# rotsee.env.table=cbind(rotsee.env.table.sort[,c(1:13,15,16,17,37)],(rotsee.env.table.sort[,c(14,18:34,36:40,42,44:58)]))    # without June
# rotsee.env.table.log=cbind( rotsee.env.table.sort[,c(1:13,15,16,17,37)],log1p(rotsee.env.table.sort[,c(14,18:36,38:58)]))
# envvector<-c(17:38,40:43,45:55)
# rotsee.env.table[,envvector]
# rotsee.mobs.env.table=rotsee.env.table
# rotsee.mobs.env.table.log=rotsee.env.table.log
# 
# 
# rotsee.env.table=cbind(rotsee.env.table.sort[,c(1:13,15,16,17,37)],(rotsee.env.table.sort[,c(14,18:34,36:40,42,44:58)]))   # June 2013
# rotsee.env.table.log=cbind( rotsee.env.table.sort[,c(1:13,15,16,17,37)],log1p(rotsee.env.table.sort[,c(14,18:36,38:58)]))
# envvector<-c(17:38,40:43,45:55)
# rotsee.env.table[,envvector]
# rotsee.mobs.env.table=rotsee.env.table
# rotsee.mobs.env.table.log=rotsee.env.table.log
# 
#   
# rotsee.env.table=cbind(rotsee.env.table.sort[,c(1:13,15,16,17,37)],(rotsee.env.table.sort[,c(14,18:34,36:40,42,44:58)]))    # August 2013
# rotsee.env.table.log=cbind( rotsee.env.table.sort[,c(1:13,15,16,17,37)],log1p(rotsee.env.table.sort[,c(14,18:36,38:58)]))
# envvector<-c(17:38,40:43,45:55)
# rotsee.env.table[,envvector]
# rotsee.mobs.env.table=rotsee.env.table
# rotsee.mobs.env.table.log=rotsee.env.table.log
# 
# 
# rotsee.env.table=cbind(rotsee.env.table.sort[,c(1:13,15,16,17,37)],(rotsee.env.table.sort[,c(14,18:34,36:40,42,44:58)]))   # September 2014
# rotsee.env.table.log=cbind( rotsee.env.table.sort[,c(1:13,15,16,17,37)],log1p(rotsee.env.table.sort[,c(14,18:36,38:58)]))
# (envvector<-c(17:47))
# (rotsee.env.table[,envvector])
# (rotsee.mobs.env.table=rotsee.env.table)
# (rotsee.mobs.env.table.log=rotsee.env.table.log)
# 
# 
# rotsee.env.table=cbind(rotsee.env.table.sort[,c(1:13,15,16,17,37)],(rotsee.env.table.sort[,c(14,18:34,36:40,42,44:58)]))    # September 2015
# rotsee.env.table.log=cbind( rotsee.env.table.sort[,c(1:13,15,16,17,37)],log1p(rotsee.env.table.sort[,c(14,18:36,38:58)]))
# envvector<-c(17:27,32:40,42,44:53)
# rotsee.env.table[,envvector]
# rotsee.mobs.env.table=rotsee.env.table
# rotsee.mobs.env.table.log=rotsee.env.table.log
# 

# run this lines 
#rotsee.mobs.env.table=rotsee.env.table
#rotsee.mobs.env.table.log=rotsee.env.table.log



########################################################################################################################################
# plot diversity measures
# use the prefiltered absolut readcount for this (i.e. bac.rotsee.prune.filt3 when incorporating 16S rRNA based MOBs or  bac.rotsee.prune.noMOBS  if not)
########################################################################################################################################

richness.16s = estimate_richness(bac.rotsee.prune.filt3, split = TRUE, measures = c("Observed",  "chao", "Shannon", "Simpson", "InvSimpson", "Fisher") )
richness.16s = cbind(richness.16s, rotsee.env.table[,c("Depth_m", "Campaign")])
library(reshape2)
plot_richness(bac.rotsee.prune.filt3, x="Depth_m", color = "Campaign" ,
              measures=c("Observed",  "chao", "Shannon", "Simpson", "InvSimpson", "Fisher")) 


depthrichplot =function (physeq, x = "samples", color = NULL, shape = NULL, 
          title = NULL, scales = "free_y", nrow = 1, shsi = NULL, measures = NULL, 
          sortby = NULL) 
{
  erDF = estimate_richness(physeq, split = TRUE, measures = measures)
  measures = colnames(erDF)
  ses = colnames(erDF)[grep("^se\\.", colnames(erDF))]
  measures = measures[!measures %in% ses]
  if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
    DF <- data.frame(erDF, sample_data(physeq))
  }
  else {
    DF <- data.frame(erDF)
  }
  if (!"samples" %in% colnames(DF)) {
    DF$samples <- sample_names(physeq)
  }
  if (!is.null(x)) {
    if (x %in% c("sample", "samples", "sample_names", "sample.names")) {
      x <- "samples"
    }
  }
  else {
    x <- "samples"
  }
  mdf = reshape2::melt(DF, measure.vars = measures)
  mdf$se <- NA_integer_
  if (length(ses) > 0) {
    selabs = ses
    names(selabs) <- substr(selabs, 4, 100)
    substr(names(selabs), 1, 1) <- toupper(substr(names(selabs), 
                                                  1, 1))
    mdf$wse <- sapply(as.character(mdf$variable), function(i, 
                                                           selabs) {
      selabs[i]
    }, selabs)
    for (i in 1:nrow(mdf)) {
      if (!is.na(mdf[i, "wse"])) {
        mdf[i, "se"] <- mdf[i, (mdf[i, "wse"])]
      }
    }
    mdf <- mdf[, -which(colnames(mdf) %in% c(selabs, "wse"))]
  }
  if (!is.null(measures)) {
    if (any(measures %in% as.character(mdf$variable))) {
      mdf <- mdf[as.character(mdf$variable) %in% measures, 
                 ]
    }
    else {
      warning("Argument to `measures` not supported. All alpha-diversity measures (should be) included in plot.")
    }
  }
  if (!is.null(shsi)) {
    warning("shsi no longer supported option in plot_richness. Please use `measures` instead")
  }
  if (!is.null(sortby)) {
    if (!all(sortby %in% levels(mdf$variable))) {
      warning("`sortby` argument not among `measures`. Ignored.")
    }
    if (!is.discrete(mdf[, x])) {
      warning("`sortby` argument provided, but `x` not a discrete variable. `sortby` is ignored.")
    }
    if (all(sortby %in% levels(mdf$variable)) & is.discrete(mdf[, 
                                                                x])) {
      wh.sortby = which(mdf$variable %in% sortby)
      mdf[, x] <- factor(mdf[, x], levels = names(sort(tapply(X = mdf[wh.sortby, 
                                                                      "value"], INDEX = mdf[wh.sortby, x], mean, na.rm = TRUE, 
                                                              simplify = TRUE))))
    }
  }
  richness_map = aes_string(x = x, y = "value", colour = color, 
                            shape = shape)
  p = ggplot(mdf, richness_map) + geom_point(na.rm = TRUE)
  if (any(!is.na(mdf[, "se"]))) {
    p = p + geom_errorbar(aes(ymax = value + se, ymin = value - 
                                se), width = 0.1)
  }
  p = p + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, 
                                           hjust = 0))
  p = p + ylab("Alpha Diversity Measure")
  p = p + coord_flip() + scale_x_reverse() 
  p=  p + geom_line()
  p = p + facet_wrap(~variable, nrow = nrow, scales = "free")
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}


#pdf("alpha diversity 16s.pdf", useDingbats=FALSE, width=18, height=10)
depthrichplot(bac.rotsee.prune.filt3, x="Depth_m", color = "Campaign" ,
              measures=c("Observed",  "chao", "Shannon", "Simpson", "InvSimpson", "Fisher"))
dev.off()


########################################################################################################################################
## Check multivariate Beta-Diversity between sample groupings
########################################################################################################################################

## Calculate Bray-Curtis distances between samples
dis <- vegdist(t(otu_table(bac.rotsee.prune.filt3)))
## Calculate multivariate dispersions of different samples grouped by factors
envbetadisp_oxi= get_variable(bac.rotsee.prune.filt3)$Oxidation_Zone
envbetadisp_campaign=get_variable(bac.rotsee.prune.filt3)$Campaign

## choose the the factor to test (below you can test the oxidation zones only)
envbetadisp_oxi_campaign=interaction(envbetadisp_oxi, envbetadisp_campaign)   # oxidation zone * campaign

## compare betadiversity of oxyc zones
mod <- betadisper(dis, group=envbetadisp_oxi_campaign,
                  type="median", bias.adjust=TRUE)
mod

## Perform test
anova(mod)
permutest(mod, pairwise = TRUE, permutations = 99)
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)
which(mod.HSD$group[,4]<=0.05)

## Plot BeataDispersion PCoA
plot(mod, ellipse = TRUE, hull = FALSE) # 1 sd data ellipse
plot(mod, ellipse = TRUE, hull = FALSE, conf = 0.90) # 90% data ellipse

## Draw a boxplot of the distances to centroid for each group
boxplot(mod)

## some additional boxplots between zones and seasons
betadistances=mod$distances
betagroups=mod$group
beta_anova=data.frame(betadistances, (betagroups))
beta_anova$betadistances_logit=car::logit(betadistances, percents=FALSE)


barorder=c("Oxic.A_June_2013", "Oxycline.A_June_2013" , "Oxidation_Zone.A_June_2013",
           "Oxic.B_August_2013", "Oxycline.B_August_2013", "Oxidation_Zone.B_August_2013", "Anoxic.B_August_2013",
           "Oxic.C_September_2014", "Oxycline.C_September_2014", "Oxidation_Zone.C_September_2014", "Anoxic.C_September_2014",
           "Oxic.D_September_2015",  "Oxidation_Zone.D_September_2015", "Anoxic.D_September_2015")

#Anova
options(contrasts=c(factor="contr.sum",ordered="contr.poly"))      ## SPSS type III ... marginal 
(Anova(Tukeylm<-lm(betadistances_logit~X.betagroups., data=(beta_anova),
                   contrasts=list(ind=contr.sum)),
       type=3, singular.ok=TRUE))
(TukeComp<-summary(glht(Tukeylm, linfct= mcp(X.betagroups.="Tukey"), p.adjust="none")))
mod.cld <- cld(TukeComp)

(plot.posthoc.annot <-(mod.cld$mcletters$Letters))
plot.posthoc.annot= data.frame(plot.posthoc.annot)
plot.posthoc.annot$id=rownames(plot.posthoc.annot)
plot.posthoc.annot=transform(plot.posthoc.annot, id=match(id, barorder))
plot.posthoc.annot=arrange(plot.posthoc.annot, id)
## boxplot of betadisper between zones and seasons
p <- ggplot(beta_anova, aes(x = X.betagroups., y = betadistances)) + 
  geom_boxplot( aes(colour=X.betagroups.),
               alpha =0.8,
               cex=1,
               width=0.7) +
  
  scale_x_discrete(name="", limits = barorder,
                   labels=c("Oxic June13", "Oxycline June13", "Oxidation Zone June13",
                            "Oxic August13", "Oxycline August13", "Oxidation Zone August13", "Anoxic August13",
                            "Oxic September14", "Oxycline September14", "Oxidation Zone September14", "Anoxic September14",
                            "Oxic September15", "Oxidation Zone September15", "Anoxic September15")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = c(3.5, 7.5, 11.5)) +
  scale_y_continuous(name="Betadiversity", expand=c(0,0), breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6) , 
                                limits=c(0, 0.62) ) +
  annotate("text", label=paste(plot.posthoc.annot[,1]), y= 0.05 , x= c(1:14), cex=4)
p


## compare betadiversity between oxyc zones over all
envbetadisp_oxi_campaign=interaction(envbetadisp_oxi)      # oxidation Zone

## compare betadiversity of oxyc zones
mod <- betadisper(dis, group=envbetadisp_oxi_campaign,
                  type="median", bias.adjust=TRUE)
mod

## Perform test
anova(mod)
permutest(mod, pairwise = TRUE, permutations = 99)
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)
which(mod.HSD$group[,4]<=0.05)

## Plot BeataDispersion PCoA
plot(mod, ellipse = TRUE, hull = FALSE) # 1 sd data ellipse
plot(mod, ellipse = TRUE, hull = FALSE, conf = 0.90) # 90% data ellipse

## Draw a boxplot of the distances to centroid for each group
boxplot(mod)
## Beatadiver boxplots between zones
## some additional boxplots
betadistances=mod$distances
betagroups=mod$group
beta_anova=data.frame(betadistances, (betagroups))
beta_anova$betadistances_logit=car::logit(betadistances, percents=FALSE)

barorder=c("Oxic", "Oxycline" , "Oxidation_Zone", "Anoxic")


#Anova
options(contrasts=c(factor="contr.sum",ordered="contr.poly"))      ## SPSS type III ... marginal 
(Anova(Tukeylm<-lm(betadistances_logit~X.betagroups., data=(beta_anova),
                   contrasts=list(ind=contr.sum)),
       type=3, singular.ok=TRUE))
(TukeComp<-summary(glht(Tukeylm, linfct= mcp(X.betagroups.="Tukey"), p.adjust="holmes")))
mod.cld <- cld(TukeComp)


(plot.posthoc.annot <-(mod.cld$mcletters$Letters))
plot.posthoc.annot= data.frame(plot.posthoc.annot)
plot.posthoc.annot$id=rownames(plot.posthoc.annot)
plot.posthoc.annot=transform(plot.posthoc.annot, id=match(id, barorder))
plot.posthoc.annot=arrange(plot.posthoc.annot, id)


p <- ggplot(beta_anova, aes(x = X.betagroups., y = betadistances)) + 
  geom_boxplot( aes(colour=X.betagroups.),
                alpha =0.8,
                cex=1,
                width=0.7) +
  
  scale_x_discrete(name="", limits = barorder,
                   labels=c("Oxic", "Oxycline" , "Oxidation_Zone", "Anoxic")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +
  scale_y_continuous(name="Betadiversity", expand=c(0,0), breaks = c(0.2,0.3,0.4,0.5,0.6) , 
                     limits=c(0.19, 0.62) ) +
  annotate("text", label=paste(plot.posthoc.annot[,1]), y= 0.6 , x= c(1:4), cex=4)
p





########################################################################################################################################
# function to produce data sets that can be used for multivariate methods based on  the vegan package
########################################################################################################################################
rotsee.otu.table.t<-t(rotsee.otu.table)     # this OTU table can be read by vegan functions, it will plot OTU numbers. See the function below...
rotsee.nomobs.otu.table.t<-t(rotsee.nomobs.otu.table)
rotsee.mobs.table.t<-t(rotsee.mobs.otu.table)

#rotsee.asigned.otu.table<- merge(rotsee.otu.table, rotsee.tax.table, by="row.names", all=TRUE) # a merged table

## function which asigns the deepest available Taxonomic rank ("Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species") to a specific OTU

make_vegan_otu_table <- function (otu.table.qiime, tax.table.qiime) {
                                  
                                    otu.table.qiime.t<-t(otu.table.qiime)
                                    otu.table.qiime.to.vegan <- otu.table.qiime.t
                                    
                                    for (i in c(1:length(colnames(otu.table.qiime.t)))) {
                                      
                                       if (!is.na(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Species"][[1]])){
                                         colnames(otu.table.qiime.to.vegan)[i]<- tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Species"][[1]]
                                       }                        # use this bit of code if analysis goes down to Species Level...replace following "if" statement with "else if"
                                      
                                      else if (!is.na(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Genus"][[1]])){
                                            colnames(otu.table.qiime.to.vegan)[i]<- tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Genus"][[1]]
                                      }
                                      else if (!is.na(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Family"][[1]])){
                                        colnames(otu.table.qiime.to.vegan)[i]<- tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Family"][[1]]
                                      }
                                      else if (!is.na(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Order"][[1]])){
                                        colnames(otu.table.qiime.to.vegan)[i]<- tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Order"][[1]]
                                      }
                                      else if (!is.na(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Class"][[1]])){
                                        colnames(otu.table.qiime.to.vegan)[i]<- tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Class"][[1]]
                                      }
                                      else if (!is.na(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Phylum"][[1]])){
                                        colnames(otu.table.qiime.to.vegan)[i]<- tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Phylum"][[1]]
                                      }
                                      else if (!is.na(tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Kingdom"][[1]])){
                                        colnames(otu.table.qiime.to.vegan)[i]<- tax.table.qiime[ which(rownames(tax.table.qiime) == colnames(otu.table.qiime.t)[i]),"Kingdom"][[1]]
                                      }
                                      else {colnames(otu.table.qiime.to.vegan)[i]<-"unknown OTU"}
                                    }
                                    
                                    return(otu.table.qiime.to.vegan)
                                }
                                      
## run the function
rotsee.vegan.otu.table<- make_vegan_otu_table (rotsee.otu.table, rotsee.tax.table) # this is a table with the deepest Phylogenetic information gathered by searching respective databases in QIIME
rotsee.vegan.otu.table.16s = rotsee.vegan.otu.table    ## will be used in other scripts

rotsee.vegan.nomobs.otu.table<- make_vegan_otu_table (rotsee.nomobs.otu.table, rotsee.nomobs.tax.table)
rotsee.vegan.mobs.otu.table<- make_vegan_otu_table (rotsee.mobs.otu.table, rotsee.mobs.tax.table)


# quickly check if everything went well   
colnames(rotsee.vegan.otu.table)[which(rownames(rotsee.otu.table)=="OTU_1686")]  # should be a LD19 or character(0) if MOBs were removed (i.e. if rotsee.nomobs.otu.table.t was used)
colnames(rotsee.vegan.nomobs.otu.table)[which(rownames(rotsee.nomobs.otu.table)=="OTU_1686")]
colnames(rotsee.vegan.otu.table)[which(rownames(rotsee.otu.table)=="OTU_50")]  # should be a Isosphaeraceae

## you can use this function to write the produced OTU table into your working directory
# write.csv(rotsee.vegan.otu.table, "vegan_test.csv")


## remove A1 to A5 for the procruste test in the procruste analysis work sheet!
#rotsee.vegan.otu.table=rotsee.vegan.otu.table[-c(which(rownames(rotsee.vegan.otu.table)== c("A01", "A02" ,"A03", "A05", "D12"))),]
#rotsee.vegan.otu.table=rotsee.vegan.otu.table[-c(which(rownames(rotsee.vegan.otu.table)== c("A01", "A02" ,"A03", "A05", "D12"))),]



## function to replace variables
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
## write the fasta files of the MOBs only with otu number and deepest phylogenetic information (can be used for blasting and crosscomparison with the pMOA based phylogeny)
library(Biostrings)
refseq_mobs = refseq(bac.rotsee.MOBs)
refseq_mobs@ranges@NAMES <- paste(mgsub(c("\\.", "/","-"),c("","_","_"),(refseq_mobs@ranges@NAMES)),"_",colnames(rotsee.vegan.mobs.otu.table), sep="") # replace characters that  are not parsed propperly
writeXStringSet(refseq_mobs,file="seq_MOBs.fas",format="fasta")  

## data is now prepared for some multivariate analysis with vegan...


## this can be used for the phylocom analysis
# comm.16s =otu_table(bac.rotsee)
# (newtaxanames <- paste(colnames(rotsee.vegan.otu.table),"_",rownames(comm.16s), sep=""))
# bac.rotsee.phytree=bac.rotsee
# taxa_names(bac.rotsee.phytree) <- newtaxanames
# phy.tree.16s =phy_tree(bac.rotsee.phytree)
# write.tree(phy.tree.16s, file = "testtree16s")


########################################################################################################################################
########################################################################################################################################
## "in-depth" NMDS analysis
########################################################################################################################################
########################################################################################################################################

########################################################################################################################################
## NMDS
set.seed(77)
rotsee.mds.allOTUs.16s <- rotsee.mds.allOTUs <- metaMDS(rotsee.vegan.otu.table, distance="bray", trace=TRUE, plot=TRUE, k=2)
## some analytics of the quality of the NMDS 
stressplot(rotsee.mds.allOTUs)
(goodpoint=goodness(rotsee.mds.allOTUs)) #sum of squared values is equal to squared stress

## you can rotate the MDS according to represent the depth succession on the first axis (or any other paramter...)
rotsee.mds.allOTUs <- with(rotsee.env.table, MDSrotate(rotsee.mds.allOTUs, Depth_m))
plot(rotsee.mds.allOTUs)

########################################################################################################################################
## define some colors for the plot
col.depth<-(colorRampPalette(brewer.pal(5,"Blues"))(1000)) # depth scale
col.fact1<-(colorRampPalette(brewer.pal(5,"Greens"))(1000))
col.fact2<-(colorRampPalette(brewer.pal(5,"Reds"))(1000))
YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
col.fact3<-(colorRampPalette(YlOrBr, space="Lab"))

(genus.numbers <- length(unique(rotsee.tax.table[,"Phylum"])))
#col.genus<-(colorRampPalette(brewer.pal(genus.numbers,"Reds"))(genus.numbers))  # choose this or the next color
col.genus<-(colorRampPalette(brewer.pal(5,"Greens"))(genus.numbers))
#col.genus<-colorRampPalette(c("brown", "red", "yellow","blue","green"))(genus.numbers)
#col.genus<-distinctColorPalette(genus.numbers)

# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# col.genus=sample(col_vector, genus.numbers)

#col.genus<-colorRampPalette(c("red", "orange","cyan", "darkblue"))(genus.numbers) # or choose manyd different colors
genus.df<-transform(rotsee.tax.table, id=match(Phylum, unique(Phylum)))

(genus.numbers.mobs <- length(unique(rotsee.mobs.tax.table[,"Class"])))
col.genus.mobs<-colorRampPalette(c("red", "orange"))(genus.numbers.mobs) # or choose manyd different colors
genus.df.mobs<-transform(rotsee.mobs.tax.table, id=match(Class, unique(Class)))

# vector for species point sizes according to read counts standardized to seqdepth
bac.rotsee.otu.sub = otu_table(bac.rotsee)
min(taxa_sums(bac.rotsee.otu.sub))
max(taxa_sums(bac.rotsee.otu.sub))
species.size.vec <- taxa_sums(bac.rotsee.otu.sub)
x=species.size.vec
normalized.species.size.vec = (x-min(x))/(max(x)-min(x))*3+0.5

## subset for the procrust analysis with pmoA ... just run this to get rid of 5 sites not included in the pmoA sequencing
## otherwise skip this line

########################################################################################################################################
## basic nmds plot... will be used for further analysis.... use this code block whenever you want to reset and add stuff from below
#----------------------------------------------------------------------------------------------------------------------------------------
#pdf("NMDS_16sAllSeasonsMDS.pdf" ,useDingbats=FALSE)      ## if you use this command, there will be a pdf produced of anything which gets printed to device below...until you type "dev.off()"

par(mar=c(4.1, 4.1, 4.1, 7.1), xpd=TRUE)
 resetPar <- function() {
 dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
 }
 #par(resetPar())   # you can run this function / line to reset par values to default
plot(rotsee.mds.allOTUs, type="n", xlim=c(-1.7,1.6), ylim=c(-1.8, 1.6))
points(rotsee.mds.allOTUs, display ="species", cex=0.8, pch=21, col="grey", bg=alpha("lightgrey", 0.8), lwd=0)
#points(rotsee.mds.allOTUs, display ="species", cex=1,#normalized.species.size.vec * 1.2 ,#cex=1.5,
#       pch=21, col=alpha(col.genus[genus.df$id], 1), bg=alpha(col.genus[genus.df$id], 0.5), lwd=0)
#legend("topright", inset=c(-0.22,0),c(unique(rotsee.tax.table[,"Phylum"])),pch=c(19),col=alpha(unique(col.genus[genus.df$id]), 0.7), bty="n", cex=0.6)
#points(rotsee.mds.allOTUs, display ="sites", cex=c(goodpoint)*400, alpha= 2, pch=as.numeric(rotsee.env.table$Campaign)+20, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
#       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.5))
## make pch not round A viereggli, B sind diamonds, C sind drüüeggli, D sind inversed drüüeglli
points(rotsee.mds.allOTUs, display ="sites", cex=3.5, alpha= 2, pch=as.numeric(as.factor(rotsee.env.table$Campaign))+21, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.8))
text(rotsee.mds.allOTUs$points[, 1]+0.05, rotsee.mds.allOTUs$points[, 2]+0.05,
     labels=unlist(dimnames(rotsee.mds.allOTUs$points)[1]),
     col="darkblue",#col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
     cex=0.7)
legend(-1.2,1.4,legend=as.expression(paste("Stress:\n",
                                          round(c(rotsee.mds.allOTUs$stress)*100,2),
                                          "%")),cex=0.8,bty="n")
scaleplot=image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.145,0.33),axis.args = list(cex.axis = 0.7), 
           zlim=range(-(rotsee.env.table$Depth_m)) )
########################################################################################################################################
#----------------------------------------------------------------------------------------------------------------------------------------
## next bolcks are adding different informations to the basic plot (i.e. identification of 16S rRNA based MOBs etc.)
#----------------------------------------------------------------------------------------------------------------------------------------
## plot MOB OTU's on the total community plot with the color according to if smmo is potentially present in the respective OTU
#----------------------------------------------------------------------------------------------------------------------------------------
pmmo_color=c("red", "violet")
(mobmatch = is.na(match(MOBs, rownames(rotsee.otu.table))))
MOBs_present<-MOBs[mobmatch=="FALSE"]
MOBsOTUsnumbers=c()
for (i in c(1:length(MOBs_present))){
  
  MOBsOTUsnumbers[i]<- which(colnames(t(rotsee.otu.table)) == MOBs_present[i])
}
MOBsOTUsnumbers
points(rotsee.mds.allOTUs$species[MOBs_reordered,1], rotsee.mds.allOTUs$species[MOBs_reordered,2],
       cex=2, pch=24, col="orange", bg=alpha(pmmo_color[pMMO_vector], 0.2))
text(rotsee.mds.allOTUs$species[MOBsOTUsnumbers,1]-0.05, rotsee.mds.allOTUs$species[MOBsOTUsnumbers,2]-0.05,
     colnames(rotsee.vegan.otu.table)[MOBsOTUsnumbers], cex=0.7, col=pmmo_color[pMMO_vector])
#dev.off()

## or assign with OTU
     text(rotsee.mds.allOTUs$species[MOBsOTUsnumbers,1], rotsee.mds.allOTUs$species[MOBsOTUsnumbers,2]+0.05,
     colnames(rotsee.otu.table.t)[MOBsOTUsnumbers], col=pmmo_color[pMMO_vector], cex = 0.7)
     
########################################################################################################################################
## plot OTU's on the total community plot that are highly abundant
#----------------------------------------------------------------------------------------------------------------------------------------   
plotabundat<- genefilter_sample(bac.rotsee, filterfun_sample(function(nrotu) nrotu > 1000), A=10) # larger 10000 reads at least in 10 samples
(abundantvect = as.numeric(which(plotabundat)))
         points(rotsee.mds.allOTUs$species[abundantvect,1], rotsee.mds.allOTUs$species[abundantvect,2],
          cex=2, pch=24, col="Gold", bg=alpha("Gold", 0.5))
          text(rotsee.mds.allOTUs$species[abundantvect ,1], rotsee.mds.allOTUs$species[abundantvect ,2],
          colnames(rotsee.vegan.otu.table)[abundantvect ], cex=0.7, col="black")
        
          ## or assign with OTU
          text(rotsee.mds.allOTUs$species[abundantvect,1], rotsee.mds.allOTUs$species[abundantvect,2]+0.1,
          colnames(rotsee.otu.table.t)[abundantvect], col="darkgrey", cex = 0.7)
     
########################################################################################################################################     
## additional graphical elements such as environmental fittings or GAMs
#----------------------------------------------------------------------------------------------------------------------------------------     
     
## environmenal fitting (this is a linear fitting on the plot, it is usually better to use the GAM (Generalized Additive Model)
## fitting from further below)-
(mds.allOTUenvfit <- envfit(rotsee.mds.allOTUs, rotsee.env.table[,-c(1:10,12)], na.rm=TRUE, permutations=1000))
     #sink(file = "nvfit on all 16s.txt", append = FALSE, type = c("output"),
     #     split = FALSE)
     #mds.allOTUenvfit 
     #sink()
     which(mds.allOTUenvfit$vectors$r >= 0.6)
     plot(mds.allOTUenvfit, choices = c(1,2), axis = FALSE, p.max = 0.05, col = "lightgrey",  add = TRUE,
     cex=0.7)

## most of the variables are significantly correlated with the OTU structure. This is due to the "congruent" structuring throughout the depth profile
## also note that half of the samples were kicked out as there are missing values(influencing the varibables that actually have no missing values as well).
#---------------------------------------------------------------------------------------------------------------------------------------- 
## OTU fitting, here you can aposteriori fitt OTUs to check out how well they're representated on the NMDS
#(mds.allOTUfit <- envfit(rotsee.mds.allOTUs, rotsee.vegan.otu.table[,1:10], na.rm=TRUE, permutations=100))
#plot(mds.allOTUfit, choices = c(1,2), axis = FALSE, p.max = 0.01, col = "lightgrey",  add = TRUE, cex=0.8)
mds.allOTUfit <- envfit(rotsee.mds.allOTUs, rotsee.vegan.otu.table[,1:10], na.rm=TRUE, permutations=100)
plot(mds.allOTUfit, choices = c(1,2), axis = FALSE, p.max = 0.01, col = "lightgrey",  add = TRUE, cex=0.8)

#---------------------------------------------------------------------------------------------------------------------------------------- 
## Confidence ellipses
with(rotsee.env.table, ordiellipse(rotsee.mds.allOTUs, groups= interaction(Campaign), display="sites",
                                draw="polygon", col=c( "red","darkorange", "orange"), border= "lightgrey", alpha=20 ))

with(rotsee.env.table, ordiellipse(rotsee.mds.allOTUs, groups= interaction(Oxidation_Zone), display="sites",
                                   draw="polygon", col=c( "darkblue","green","lightblue","orange"), border= "lightgrey", alpha=20 , label=FALSE))

with(rotsee.env.table, ordiellipse(rotsee.mds.allOTUs, groups= interaction(Oxidation_Zone), display="sites",
                                   draw="polygon", col=c( "darkblue","green","lightblue","orange"), border= "lightgrey", alpha=20 , label=FALSE, 
                                   kind="se", conf=0.95))

#---------------------------------------------------------------------------------------------------------------------------------------- 
## GAMS for environmental variables can be plotted on the ordination
gam.PhysicoChemRaw<-ordisurf(as.formula(paste("rotsee.mds.allOTUs ~ ",paste(colnames(rotsee.env.table[parame<-"Depth_m"])))),rotsee.env.table,
                             choices=c(1,2), knots=10,                                         
                             family="gaussian", isotropic= TRUE,
                             col="black",scaling=3,
                             add = TRUE, display = "sites",
                             nlevels = 15, labcex = 0.8,
                             bubble = TRUE, cex = 1, select = TRUE, method = "REML",
                             gamma = 1, plot = TRUE, lwd=0.5)
(summary.gam.PhysicoChemRaw <- summary(gam.PhysicoChemRaw) ) # this summary shows how well a variable correlates with the OTU structure

legend(0.5,-1.2,legend=as.expression(paste("Explained deviance of\n",parame,":\n",
                                           round(c(summary.gam.PhysicoChemRaw$dev.expl)*100,2),
                                           "%,   ","P=",
                                           round(c(summary.gam.PhysicoChemRaw$s.pv[[1]]),3))),cex=0.8,bty="n")

#---------------------------------------------------------------------------------------------------------------------------------------- 
## plot GAMs for specific depth range, this can be interesting, to within which depth range actual correlations occure 

lowerlimit= 9.5
upperlimit= 2

# above lower limit [m]
rotsee.mds.allOTUs.8<-scores(rotsee.mds.allOTUs)[which(rotsee.env.table$Depth_m <= lowerlimit ),c(1,2)]
rotsee.env.table.8<-rotsee.env.table[which(rotsee.env.table$Depth_m <= lowerlimit),]

gam.PhysicoChemRaw<-ordisurf(rotsee.mds.allOTUs.8[,c(1,2)],  rotsee.env.table.8[,parame<-"mmoX_copies.L"],
                             choices=c(1,2), knots=10,                                         
                             family="gaussian", isotropic= TRUE,
                             col="orange",scaling=3,
                             add = TRUE, display = "sites",
                             nlevels = 10, labcex = 0.8,
                             bubble = TRUE, cex = 1, select = TRUE, method = "REML",
                             gamma = 1, plot = TRUE, lwd=0.5)
(summary.gam.PhysicoChemRaw <- summary(gam.PhysicoChemRaw) )

legend(-1.5,-0.7,legend=as.expression(paste("Explained deviance of\n",parame,"\n lower water column:\n",
                                           round(c(summary.gam.PhysicoChemRaw$dev.expl)*100,2),
                                           "%,   ","P=",
                                           round(c(summary.gam.PhysicoChemRaw$s.pv[[1]]),3))),cex=0.8,bty="n")


# between upper and lower limit
rotsee.mds.allOTUs.8<-scores(rotsee.mds.allOTUs)[which(rotsee.env.table$Depth_m <= lowerlimit & rotsee.env.table$Depth_m > upperlimit),c(1,2)]
rotsee.env.table.8<-rotsee.env.table[which(rotsee.env.table$Depth_m <= lowerlimit & rotsee.env.table$Depth_m > upperlimit),]

gam.PhysicoChemRaw<-ordisurf(rotsee.mds.allOTUs.8[,c(1,2)],  rotsee.env.table.8[,parame<-"Cu_part_nM"],
                             choices=c(1,2), knots=10,                                         
                             family="gaussian", isotropic= TRUE,
                             col="red",scaling=3,
                             add = TRUE, display = "sites",
                             nlevels = 10, labcex = 0.8,
                             bubble = TRUE, cex = 1, select = TRUE, method = "REML",
                             gamma = 1, plot = TRUE, lwd=0.5)
(summary.gam.PhysicoChemRaw <- summary(gam.PhysicoChemRaw) )

legend(-1.3,-1.2,legend=as.expression(paste("Explained deviance of\n",parame,"\n between 3m and 8m:\n",
                                            round(c(summary.gam.PhysicoChemRaw$dev.expl)*100,2),
                                            "%,   ","P=",
                                            round(c(summary.gam.PhysicoChemRaw$s.pv[[1]]),3))),cex=0.8,bty="n")

#---------------------------------------------------------------------------------------------------------------------------------------- 
## plot bar for depth with locater
#position<-locator(n=1)
#position$x= 1.5
#$y= -1.3
#colorbar.plot(position$x, position$y, strip=c(1:1000), strip.width = 0.05, strip.length = 0.5, zrange=
#              adj.x = 0.5, adj.y = 0.5, col = rev(col.depth), 
#              horizontal = FALSE)
#text(position$x,position$y,"depth")


#---------------------------------------------------------------------------------------------------------------------------------------- 
## make pdf for gams of all environmental variables (will be safed in working directory)

# here choose the env variables without NAs for specific campaign 
# normally use the envvector defined for the specific campaign before
head(rotsee.env.table[,-c(1:13,15:17)])

gamvector=c(1:length(rotsee.env.table))
gamvector=gamvector[-c(1:13,15:17, 48:67 )]  # remove factors

sink(file = "GAMs on 16s.txt", append = FALSE, type = c("output"),
     split = FALSE)

for (i in c(gamvector)){
  
  
  pdf(paste(colnames(rotsee.env.table[i]), " GAM fitting",".pdf", sep=""),useDingbats=FALSE) #, width=4.5, height=4)
  
  par(mar=c(4.1, 4.1, 4.1, 7.1), xpd=TRUE)
  plot(rotsee.mds.allOTUs, type="n", xlim=c(-1.7,1.6), ylim=c(-1.8, 1.6))
  
  points(rotsee.mds.allOTUs, display ="species", cex=0.5, pch=21, col="orange", bg=alpha("orange", 0.01), lwd=0)
  #points(rotsee.mds.allOTUs, display ="species", cex=normalized.species.size.vec * 1.3 ,#cex=1.5,
  #       pch=21, col=alpha(col.genus[genus.df$id], 0.5), bg=alpha(col.genus[genus.df$id], 0.3), lwd=0)
  #legend("topright", inset=c(-0.22,0),c(unique(rotsee.tax.table[,"Phylum"])),pch=c(19),col=alpha(unique(col.genus[genus.df$id]), 0.5), bty="n", cex=0.6)
  
  #points(rotsee.mds.allOTUs, display ="sites", cex=c(goodpoint)*400, alpha= 2, pch=as.numeric(rotsee.env.table$Campaign)+20, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
  #       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.5))
  
  points(rotsee.mds.allOTUs, display ="sites", cex=4.5, alpha= 2, pch=as.numeric(as.factor(rotsee.env.table$Campaign))+21, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
         bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.6))
  text(rotsee.mds.allOTUs$points[, 1]+0.05, rotsee.mds.allOTUs$points[, 2]+0.05,
       labels=unlist(dimnames(rotsee.mds.allOTUs$points)[1]),
       col="darkblue",#col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
       cex=0.7)
  
  legend(-1.2,1.4,legend=as.expression(paste("Stress:\n",
                                             round(c(rotsee.mds.allOTUs$stress)*100,2),
                                             "%")),cex=0.8,bty="n")
  
 image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.145,0.33),axis.args = list(cex.axis = 0.7), 
                       zlim=range(-(rotsee.env.table$Depth_m)) )
  
  gam.PhysicoChemRaw<-ordisurf(as.formula(paste("rotsee.mds.allOTUs ~ ",paste(colnames(rotsee.env.table[i])))),rotsee.env.table,
                               choices=c(1,2), knots=10,                                         
                               family="gaussian", isotropic= TRUE,
                               col="darkgrey",scaling=3,
                               add = TRUE, display = "sites",
                               nlevels = 15, labcex = 0.8,
                               bubble = TRUE, cex = 1, select = TRUE, method = "REML",
                               gamma = 1, plot = TRUE)
  print(colnames(rotsee.env.table[i]))
  print(summary.gam.PhysicoChemRaw <- summary(gam.PhysicoChemRaw) )
  
  legend(-1,-1.3,legend=as.expression(paste("Explained deviance of\n ", colnames(rotsee.env.table[i]),
                                           round(c(summary.gam.PhysicoChemRaw$dev.expl)*100,2),
                                           "%,   ","P=",
                                           round(c(summary.gam.PhysicoChemRaw$s.pv[[1]]),3))),cex=0.8,bty="n")
  
    dev.off()
}

 sink() 
  
#---------------------------------------------------------------------------------------------------------------------------------------- 
## here you can interactively identify potential bacterial interactors


par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
plot(rotsee.mds.allOTUs, type="n", xlim=c(-1.5,1.6))
points(rotsee.mds.allOTUs, display ="species", cex=0.1, pch=21, col="orange", bg=alpha("orange", 0.01), lwd=0)
#points(rotsee.mds.allOTUs, display ="species", cex=normalized.species.size.vec * 1.3 ,#cex=1.5,
#       pch=21, col=alpha(col.genus[genus.df$id], 0.5), bg=alpha(col.genus[genus.df$id], 0.3), lwd=0)
#legend("topright", inset=c(-0.22,0),c(unique(rotsee.tax.table[,"Phylum"])),pch=c(19),col=alpha(unique(col.genus[genus.df$id]), 0.5), bty="n", cex=0.6)

points(rotsee.mds.allOTUs, display ="sites", cex=3.5, alpha= 2, pch=as.numeric(as.factor(rotsee.env.table$Oxidation_Zone))+21, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.5))

#points(rotsee.mds.allOTUs, display ="sites", cex=4.5, alpha= 2, pch=as.numeric(as.factor(rotsee.env.table$Campaign))+21, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
#       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.6))
# text(rotsee.mds.allOTUs$points[, 1]+0.05, rotsee.mds.allOTUs$points[, 2]+0.05,
#      labels=unlist(dimnames(rotsee.mds.allOTUs$points)[1]),
#      col="darkblue",#col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
#      cex=0.7)
#legend(-1.2,1.4,legend=as.expression(paste("Stress:\n",
#                                           round(c(rotsee.mds.allOTUs$stress)*100,2),
#                                           "%")),cex=0.8,bty="n")

scaleplot=image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.145,0.33),axis.args = list(cex.axis = 0.7), 
                     zlim=range(-(rotsee.env.table$Depth_m)) )

identify(rotsee.mds.allOTUs$species, labels=rownames(rotsee.mds.allOTUs$species), cex=0.6, col="black" )
identify(rotsee.mds.allOTUs$species, labels=rownames(rotsee.otu.table), cex=0.6, col="black" )
identify(rotsee.mds.allOTUs$species, labels=(rotsee.tax.table@.Data[,"Class"]), cex=0.6, col="grey" )

#---------------------------------------------------------------------------------------------------------------------------------------- 
## add points of species found in the OTU list that migth be interesting
(bacterianumber <- which(rownames(rotsee.mds.allOTUs$species)==c("Magnetospirillum")))  # you can identify the species you're looking for 
(bacterianumber <- which(colnames(rotsee.otu.table.t)==("OTU_1686")))                   # or a specific OTU number
points(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
       labels=rownames(rotsee.mds.allOTUs$species)[bacterianumber], cex=2, pch=25,col="black", bg=alpha("yellow", alpha=0.4))
text(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
     labels=rownames(rotsee.mds.allOTUs$species)[bacterianumber], cex=0.8)
text(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
     labels=rownames(rotsee.otu.table)[bacterianumber], cex=0.8)



(bacterianumber <- which(rotsee.tax.table@.Data[,"Phylum"]==c("Cyanobacteria")))
points(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
       labels=(rotsee.tax.table@.Data[,"Phylum"]), cex=2, pch=25,col="black", bg=alpha("yellow", alpha=0.4))
text(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
     labels=rownames(rotsee.mds.allOTUs$species)[bacterianumber], cex=0.8)
text(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2]-0.08,
     labels=rownames(rotsee.otu.table)[bacterianumber], cex=0.8)



(bacterianumber <- which(colnames(rotsee.otu.table.t)==("OTU_1686")))                   # or a specific OTU number
points(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
       labels=rownames(rotsee.mds.allOTUs$species)[bacterianumber], cex=2, pch=25,col="black", bg=alpha("yellow", alpha=0.4))
text(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
     labels=rownames(rotsee.mds.allOTUs$species)[bacterianumber], cex=0.8)




#########################################################################################################################
#########################################################################################################################
##  Similarity Percentages (SIMPER) -> i.e., finds the most distinguished OTUs for the three light zones
#########################################################################################################################
#########################################################################################################################

## this one I would run for the different seasons independently...

#(sim<- with(rotsee.env.table, simper(rotsee.vegan.otu.table[,], Light_Zone)))      # annotted with highest Taxon

(sim<- with(rotsee.env.table, simper(t(rotsee.otu.table)[,], interaction(Oxidation_Zone))))    # annotted with OTU number, use this one for plotting points on NMDS
str(sim)
sum(sim$Oxic_Oxycline$average) # total between group variance explained by OTUs
sum(sim$Oxycline_Oxidation_Zone$average)
sum(sim$Oxidation_Zone_Anoxic$average)

## extract most 20 influencal OTUS for the separation of i.e.  the oxidation zone and the oxyclinethe Light Zones
(simOTUs.Oxic_Oxycline = sim$Oxic_Oxycline$species[order(-sim$Oxic_Oxycline$average)][1:20])
(simOTUs.Oxycline_Oxidation_Zone = sim$Oxycline_Oxidation_Zone$species[order(-sim$Oxycline_Oxidation_Zone$average)][1:20]) 
(simOTUs.Oxidation_Zone_Anoxic = sim$Oxidation_Zone_Anoxic$species[order(-sim$Oxidation_Zone_Anoxic$average)][1:20]) 

## plot the most influencal speices on the nmds 
## the most important OTUs according to the simper output are now converted in taxa name index vector
simOTUs<-c(simOTUs.Oxic_Oxycline,
           simOTUs.Oxycline_Oxidation_Zone,
           simOTUs.Oxidation_Zone_Anoxic)   
simOTUsnumbers=c()
for (i in c(1:length(simOTUs))){

  simOTUsnumbers[i]<- which(colnames(t(rotsee.otu.table)) == simOTUs[i])
}

simOTUsnumbers # vector with the position of simOTUs in the vegan table

## mds plot
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
plot(rotsee.mds.allOTUs, type="n", xlim=c(-1.5,1.9))

points(rotsee.mds.allOTUs, display ="species", cex=0.1, pch=21, col="orange", bg=alpha("orange", 0.00001), lwd=0)
#points(rotsee.mds.allOTUs, display ="species", cex=normalized.species.size.vec * 1.3 ,#cex=1.5,
#       pch=21, col=alpha(col.genus[genus.df$id], 0.5), bg=alpha(col.genus[genus.df$id], 0.3), lwd=0)
#legend("topright", inset=c(-0.22,0),c(unique(rotsee.tax.table[,"Phylum"])),pch=c(19),col=alpha(unique(col.genus[genus.df$id]), 0.5), bty="n", cex=0.6)

#points(rotsee.mds.allOTUs, display ="sites", cex=c(goodpoint)*400, alpha= 2, pch=as.numeric(rotsee.env.table$Campaign)+20, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
#       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.5))

points(rotsee.mds.allOTUs, display ="sites", cex=4.5, alpha= 2, pch=as.numeric(as.factor(rotsee.env.table$Oxidation_Zone))+21, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.6))
text(rotsee.mds.allOTUs$points[, 1]+0.05, rotsee.mds.allOTUs$points[, 2]+0.05,
     labels=unlist(dimnames(rotsee.mds.allOTUs$points)[1]),
     col="darkblue",#col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
     cex=0.7)
legend(-1.2,1.4,legend=as.expression(paste("Stress:\n",
                                           round(c(rotsee.mds.allOTUs$stress)*100,2),
                                           "%")),cex=0.8,bty="n")



scaleplot=image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.145,0.33),axis.args = list(cex.axis = 0.7), 
                     zlim=range(-(rotsee.env.table$Depth_m)) )


## plot the selected Simper OTU's on the total community NMDS plot
#-------------------------------------------

points(rotsee.mds.allOTUs$species[simOTUsnumbers,1], rotsee.mds.allOTUs$species[simOTUsnumbers,2], cex=2, pch=21, col="red", bg=alpha("red", 0.1))
## annotate by highest highest taxa
text(rotsee.mds.allOTUs$species[simOTUsnumbers,1]*1.1, rotsee.mds.allOTUs$species[simOTUsnumbers,2]*1.1,
     colnames(rotsee.vegan.otu.table)[simOTUsnumbers], col="darkred", cex=0.7)
## or  annotate by OTUs
text(rotsee.mds.allOTUs$species[simOTUsnumbers,1]*1.1, rotsee.mds.allOTUs$species[simOTUsnumbers,2]*1.1,
     colnames(rotsee.otu.table.t)[simOTUsnumbers], col="darkred")

## plot MOB OTU's on the total community plot
#-------------------------------------------
pmmo_color=c("red", "violet")

(mobmatch = is.na(match(MOBs, rownames(rotsee.otu.table))))
MOBs_present<-MOBs[mobmatch=="FALSE"]

MOBsOTUsnumbers=c()
for (i in c(1:length(MOBs_present))){
  
  MOBsOTUsnumbers[i]<- which(colnames(t(rotsee.otu.table)) == MOBs_present[i])
}

MOBsOTUsnumbers

points(rotsee.mds.allOTUs$species[MOBsOTUsnumbers,1], rotsee.mds.allOTUs$species[MOBsOTUsnumbers,2], cex=2, pch=24, col="orange", bg=alpha(pmmo_color[pMMO_vector], 0.2))
text(rotsee.mds.allOTUs$species[MOBsOTUsnumbers,1], rotsee.mds.allOTUs$species[MOBsOTUsnumbers,2],
     colnames(rotsee.vegan.otu.table)[MOBsOTUsnumbers], cex=0.7, col=pmmo_color[pMMO_vector])
## or assign with OTU

     text(rotsee.mds.allOTUs$species[MOBsOTUsnumbers,1]*1.1, rotsee.mds.allOTUs$species[MOBsOTUsnumbers,2]*1.1,
     colnames(rotsee.otu.table.t)[MOBsOTUsnumbers], col="darkred")

## here you can interactively identify more potential bacterial interactors (deepest phylogenetic identifier or OTU number)
identify(rotsee.mds.allOTUs$species, labels=rownames(rotsee.mds.allOTUs$species), cex=0.6, col="grey" )
identify(rotsee.mds.allOTUs$species, labels=rownames(rotsee.otu.table), cex=0.6, col="grey" )

## add points of species found in the OTU list
(bacterianumber <- which(rownames(rotsee.mds.allOTUs$species)==c("Magnetospirillum")))  # you can identify the species you're looking for 
(bacterianumber <- which(colnames(rotsee.otu.table.t)==("OTU_1686")))                   # or a specific OTU number

points(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
       labels=rownames(rotsee.mds.allOTUs$species)[bacterianumber], cex=2, pch=25,col="black", bg=alpha("yellow", alpha=0.4))
text(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
     labels=rownames(rotsee.mds.allOTUs$species)[bacterianumber], cex=0.8)
#text(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
#     labels=rownames(rotsee.otu.table)[bacterianumber], cex=0.8)


(bacterianumber <- which(rotsee.tax.table@.Data[,"Phylum"]==c("Cyanobacteria")))
points(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
       labels=(rotsee.tax.table@.Data[,"Phylum"]), cex=2, pch=25,col="black", bg=alpha("yellow", alpha=0.4))
text(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
     labels=rownames(rotsee.mds.allOTUs$species)[bacterianumber], cex=0.8)


########################################################################################################################
########################################################################################################################
## vabund --> check this as better alternative for PERMANOVA (ADONiS) and SIMPER, use absolut count data
#########################################################################################################################
########################################################################################################################

## needs a lot of time for calculation.... skip that


# meanvar.plot(rotsee.vegan.otu.table)
# rotsee.glmdata <- mvabund(t(otu_table(bac.rotsee.all)))   # check which one to use , i.e. bac.rotsee
# #(env.glmdata<- rotsee.env.table [,] )  # here it might be good to transform environmental data as performed in the next line
# (env.glmdata<-rotsee.env.table.log)
# 
# (rotsee.glm<- manyglm(rotsee.glmdata ~ Oxidation_Zone , data=env.glmdata, family="negative.binomial",    #choose (multiple) factors for the model,   also poisson family if residual plott shows  linear or curvilinear, or a fan shape
#                      cor.type="I", test="LR",
#                      show.coef=FALSE, show.fitted=TRUE, show.residuals=FALSE ))
# plot(rotsee.glm)
# 
# (rotsee.glm.summary<-summary(rotsee.glm, resamp="residual", nBoot=100, test="LR"))
# #(rotsee.glm.anova<- anova(rotsee.glm, p.uni="none", resamp="montecarlo", test="LR",  #tst= "wald" or       #anova.manyglm   -> are liminc zones different? yes...
# #                        show.time="all", nBoot=10)) # for serious run choose 1000 bootstraps!!
# rotsee.glm.anova.single<- anova(rotsee.glm, p.uni="adjusted", nBoot=10, test="LR") # which  species are likely to be found in which "factor" aka habitat
# 
# summary(rotsee.glm.anova)
# 
# 
# rotsee.glm.anova$uni.p # three spp significant after multiple testing 
# s = sort(rotsee.glm.anova$uni.test[2,],index.return=T,decreasing=T)
# s$x[1:150]/sum(s$x) #the proportions of total test stat due to these 60 OTUs
# sum(s$x[1:150]/sum(s$x) )
# 
# (glm.imp.names <- names(s$x[1:150]))    # get the 10 bacteria which are most distinct beween tested groups/factor
# vabundOTUsnumbers=c()
# for (i in c(1:length(glm.imp.names))){
#   
#   vabundOTUsnumbers[i]<- which(colnames(t(rotsee.otu.table)) == glm.imp.names[i])
# }
# ## plot them on the nMDS
# #-------------------------------------------------------
# par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
# plot(rotsee.mds.allOTUs, type="n", xlim=c(-1.5,1.9))
# 
# #points(rotsee.mds.allOTUs, display ="species", cex=0.1, pch=21, col="orange", bg=alpha("orange", 0.00001), lwd=0)
# points(rotsee.mds.allOTUs, display ="species", cex=normalized.species.size.vec * 1.3 ,#cex=1.5,
#        pch=21, col=alpha(col.genus[genus.df$id], 0.5), bg=alpha(col.genus[genus.df$id], 0.3), lwd=0)
# legend("topright", inset=c(-0.22,0),c(unique(rotsee.tax.table[,"Phylum"])),pch=c(19),col=alpha(unique(col.genus[genus.df$id]), 0.5), bty="n", cex=0.6)
# 
# #points(rotsee.mds.allOTUs, display ="sites", cex=c(goodpoint)*400, alpha= 2, pch=as.numeric(rotsee.env.table$Campaign)+20, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
# #       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.5))
# 
# points(rotsee.mds.allOTUs, display ="sites", cex=4.5, alpha= 2, pch=as.numeric(as.factor(rotsee.env.table$Oxidation_Zone))+21, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
#        bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.6))
# text(rotsee.mds.allOTUs$points[, 1]+0.05, rotsee.mds.allOTUs$points[, 2]+0.05,
#      labels=unlist(dimnames(rotsee.mds.allOTUs$points)[1]),
#      col="darkblue",#col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
#      cex=0.7)
# legend(-1.2,1.4,legend=as.expression(paste("Stress:\n",
#                                            round(c(rotsee.mds.allOTUs$stress)*100,2),
#                                            "%")),cex=0.8,bty="n")
# 
# 
# 
# scaleplot=image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.145,0.33),axis.args = list(cex.axis = 0.7), 
#                      zlim=range(-(rotsee.env.table$Depth_m)) )
# 
# 
# ## Add the vabund OTUs 
# #-------------------------------------------
# points(rotsee.mds.allOTUs$species[vabundOTUsnumbers,1], rotsee.mds.allOTUs$species[vabundOTUsnumbers,2],cex=2, pch=23, col="black", bg=alpha("yellow", 0.2))
# text(rotsee.mds.allOTUs$species[vabundOTUsnumbers,1], rotsee.mds.allOTUs$species[vabundOTUsnumbers,2],
#      colnames(rotsee.vegan.otu.table)[vabundOTUsnumbers], cex=0.7, col="darkgreen")
# #text(rotsee.mds.allOTUs$species[vabundOTUsnumbers,1]*1.1, rotsee.mds.allOTUs$species[vabundOTUsnumbers,2]*1.1,
# #     colnames(rotsee.otu.table.t)[vabundOTUsnumbers], col="darkred")
# 
# 
# ## plot MOB OTU's on the total community plot
# #-------------------------------------------
# pmmo_color=c("red", "violet")
# 
# (mobmatch = is.na(match(MOBs, rownames(rotsee.otu.table))))
# MOBs_present<-MOBs[mobmatch=="FALSE"]
# 
# MOBsOTUsnumbers=c()
# for (i in c(1:length(MOBs_present))){
#   
#   MOBsOTUsnumbers[i]<- which(colnames(t(rotsee.otu.table)) == MOBs_present[i])
# }
# 
# MOBsOTUsnumbers
# 
# points(rotsee.mds.allOTUs$species[MOBsOTUsnumbers,1], rotsee.mds.allOTUs$species[MOBsOTUsnumbers,2], cex=2, pch=24, col="orange", bg=alpha(pmmo_color[pMMO_vector], 0.2))
# text(rotsee.mds.allOTUs$species[MOBsOTUsnumbers,1], rotsee.mds.allOTUs$species[MOBsOTUsnumbers,2],
#      colnames(rotsee.vegan.otu.table)[MOBsOTUsnumbers], cex=0.7, col=pmmo_color[pMMO_vector])
# ## or assign with OTU
# #     text(rotsee.mds.allOTUs$species[MOBsOTUsnumbers,1]*1.1, rotsee.mds.allOTUs$species[MOBsOTUsnumbers,2]*1.1,
# #     colnames(rotsee.otu.table.t)[MOBsOTUsnumbers], col="darkred")
# 
# 
# 
# ## here you can interactively identify potential bacterial interactors
# identify(rotsee.mds.allOTUs$species, labels=rownames(rotsee.mds.allOTUs$species), cex=0.6, col="grey" )
# identify(rotsee.mds.allOTUs$species, labels=rownames(rotsee.otu.table), cex=0.6, col="grey" )
# 
# ## add points of species found in the OTU list
# (bacterianumber <- which(rownames(rotsee.mds.allOTUs$species)=="Magnetospirillum"))
# points(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
#      labels=rownames(rotsee.mds.allOTUs$species)[bacterianumber], cex=2, pch=25,col="black", bg=alpha("yellow", alpha=0.4))
# text(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
#      labels=rownames(rotsee.mds.allOTUs$species)[bacterianumber], cex=0.8)
# 
# 
# (bacterianumber <- which(rotsee.tax.table@.Data[,"Phylum"]==c("Verrucomicrobia")))
# points(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
#        labels=(rotsee.tax.table@.Data[,"Phylum"]), cex=2, pch=25,col="black", bg=alpha("yellow", alpha=0.4))
# text(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
#      labels=rownames(rotsee.mds.allOTUs$species)[bacterianumber], cex=0.8)











########################################################################################################################
########################################################################################################################
## make nmds just with MOBs 
########################################################################################################################
########################################################################################################################
## for downstream analysis you have to choose one of the two standardization methods here:
## choose the seqdepth standardized to all OTUs dataset ... this standardization could be interesting when MOBs 
## are compared to their microbial environment.... anyway...the plot is not much different to the one in the next section...
##----------------------------------------------------------
rotsee.mobs.otu.only.table <- prune_taxa(MOBs, otu_table(bac.rotsee.all))
rotsee.mobs.tax.only.table <- prune_taxa(MOBs,tax_table(bac.rotsee.all))

(rotsee.vegan.MOB.table.1<- make_vegan_otu_table (rotsee.mobs.otu.only.table, rotsee.mobs.tax.only.table))
rotsee.vegan.MOB.table <- rotsee.vegan.MOB.table.1[-c(which(rowSums(rotsee.vegan.MOB.table.1)==0)),]   # if this one gives an error run the next line...otherwise skip the next line
rotsee.vegan.MOB.table <- rotsee.vegan.MOB.table.1


rotsee.mobs.env.table.1 <- rotsee.env.table
rotsee.mobs.env.table <- rotsee.mobs.env.table.1 [-c(which(rowSums(rotsee.vegan.MOB.table.1)==0)),]
rotsee.mobs.env.table <- rotsee.env.table

## run the NMDS

set.seed(77)
rotsee.mds.allMOBs<-metaMDS(rotsee.vegan.MOB.table, distance="bray")

plot(rotsee.mds.allMOBs, type="n")
col.depth<-(colorRampPalette(brewer.pal(5,"Blues"))(1000))
points(rotsee.mds.allMOBs, display ="species", cex=3, pch=25, col="orange", bg=alpha(pmmo_color[pMMO_vector_mobs], 0.2), lwd=0)
points(rotsee.mds.allMOBs, display ="sites", cex=2, alpha= 2, pch=as.numeric(rotsee.mobs.env.table$Campaign)+20,
       col=col.depth[as.numeric(rotsee.mobs.env.table$Depth_m)/max(rotsee.mobs.env.table$Depth_m)*1000],
       bg=alpha(col.depth[as.numeric(rotsee.mobs.env.table$Depth_m)/max(rotsee.mobs.env.table$Depth_m)*1000], 0.5))
text(rotsee.mds.allMOBs$points[, 1], rotsee.mds.allMOBs$points[, 2]+0.1,
     labels=unlist(dimnames(rotsee.mds.allMOBs$points)[1]),
     col="darkblue", cex=0.7)
text(rotsee.mds.allMOBs$species[, 1], rotsee.mds.allMOBs$species[, 2]+0.1,
     labels=unlist(dimnames(rotsee.mds.allMOBs$species)[1]),
     col=pmmo_color[pMMO_vector_mobs], cex=0.7)
legend(-1,0.5,legend=as.expression(paste("Stress:\n",
                                         round(c(rotsee.mds.allMOBs$stress)*100,2),
                                         "%")),cex=0.8,bty="n")
image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.3,0.7), axis.args = list(cex.axis = 0.7),
           zlim=range(-(rotsee.env.table$Depth_m)) )

## or  annotate by OTUs
text(rotsee.mds.allMOBs$species[,1], rotsee.mds.allMOBs$species[,2]-0.1,
     colnames(rotsee.mobs.table.t), col="darkred", cex=0.7)

## environmenal fitting
##----------------------------------------------------------
(mds.allMOBenvfit <- envfit(rotsee.mds.allMOBs, rotsee.mobs.env.table[,-c(1:10,12)], na.rm=TRUE, permutations=1000))
plot(mds.allMOBenvfit, choices = c(1,2), axis = FALSE, p.max = 0.05, col = "grey",  add = TRUE, cex=0.7)

## OTU fitting
(mds.allMOBfit <- envfit(rotsee.mds.allMOBs, rotsee.vegan.MOB.table[,], na.rm=TRUE, permutations=1000))
plot(mds.allMOBfit, choices = c(1,2), axis = FALSE, p.max = 0.05, col = "lightgrey",  add = TRUE, cex=0.8)

## Confidence ellipses
with(rotsee.mobs.env.table, ordiellipse(rotsee.mds.allMOBs, groups= interaction(Light_Zone), display="sites",
                                   draw="polygon", col=c( "red","darkorange", "orange"), border= "lightgrey", alpha=20 ,
                                   label=FALSE))
## ordination hulls
with(rotsee.mobs.env.table, ordihull(rotsee.mds.allMOBs, groups= interaction(Oxidation_Zone, Campaign), display="sites",
                                draw="polygon", col=c("orange"), border= "lightgrey", alpha=20,
                                label=TRUE))

#dev.off()

## add single gam's
##----------------------------------------------------------
gam.PhysicoChemRaw<-ordisurf(as.formula(paste("rotsee.mds.allMOBs ~ ",paste(colnames(rotsee.mobs.env.table[parame<-"Cu_diss_nM"])))),rotsee.mobs.env.table,
                             choices=c(1,2), knots=10,                                         
                             family="gaussian", isotropic= TRUE,
                             col="pink",scaling=3,
                             add = TRUE, display = "sites",
                             nlevels = 10, labcex = 0.8,
                             bubble = TRUE, cex = 1, select = TRUE, method = "REML",
                             gamma = 1, plot = TRUE)
(summary.gam.PhysicoChemRaw <- summary(gam.PhysicoChemRaw) )

legend(-1.2,-0.5, bty="n", legend=as.expression(paste("Explained deviance of\n",parame,"\n",   # change according to what you choose
                                           round(c(summary.gam.PhysicoChemRaw$dev.expl)*100,2),
                                           "%   ","P=",
                                           round(c(summary.gam.PhysicoChemRaw$s.pv[[1]]),3))),cex=0.8)


## or add gam of MOB's
##----------------------------------------------------------

for (i in c(1:ncol(rotsee.vegan.MOB.table))){
  
  
  pdf(paste(rownames(rotsee.mobs.env.table)[i],"_",colnames(rotsee.vegan.MOB.table)[i], " MOB fitting",".pdf", sep=""),useDingbats=FALSE) #, width=4.5, height=4)
  
   plot(rotsee.mds.allMOBs, type="n")
  col.depth<-(colorRampPalette(brewer.pal(5,"Blues"))(1000))
  points(rotsee.mds.allMOBs, display ="species", cex=3, pch=25, col="orange", bg=alpha(pmmo_color[pMMO_vector_mobs], 0.2), lwd=0)
  points(rotsee.mds.allMOBs, display ="sites", cex=2, alpha= 2, pch=as.numeric(as.factor(rotsee.mobs.env.table$Oxidation_Zone))+20, col=col.depth[as.numeric(rotsee.mobs.env.table$Depth_m)/max(rotsee.mobs.env.table$Depth_m)*1000],
         bg=alpha(col.depth[as.numeric(rotsee.mobs.env.table$Depth_m)/max(rotsee.mobs.env.table$Depth_m)*1000], 0.5))
  text(rotsee.mds.allMOBs$points[, 1]+0.1, rotsee.mds.allMOBs$points[, 2]+0.1,
       labels=unlist(dimnames(rotsee.mds.allMOBs$points)[1]),
       col="darkblue", #col.depth[as.numeric(rotsee.env.table.MOB$Depth_m)/max(rotsee.env.table.MOB$Depth_m)*1000], 
       cex=0.7)
  text(rotsee.mds.allMOBs$species[, 1]*1.1, rotsee.mds.allMOBs$species[, 2]+0.1,
       labels=unlist(dimnames(rotsee.mds.allMOBs$species)[1]),
       col=pmmo_color[pMMO_vector_mobs], cex=0.7)
  legend(-1,-1.2,legend=as.expression(paste("Stress:\n",
                                              round(c(rotsee.mds.allMOBs$stress)*100,2),
                                              "%")),cex=0.8,bty="n")
  image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.3,0.7), axis.args = list(cex.axis = 0.7),
             zlim=range(-(rotsee.env.table$Depth_m)) )
  
  ## or  annotate by OTUs
  

  
gam.PhysicoChemRaw<-ordisurf(as.formula(paste("rotsee.mds.allMOBs ~ ",paste(colnames(rotsee.mobs.table.t[,z<-i])))), rotsee.mobs.table.t,
                             choices=c(1,2), knots=10,                                         
                             family="gaussian", isotropic= TRUE,
                             col="pink",scaling=3,
                             add = TRUE, display = "sites",
                             nlevels = 15, labcex = 0.8,
                             bubble = TRUE, cex = 1, select = TRUE, method = "REML",
                             gamma = 1, plot = TRUE)
(summary.gam.PhysicoChemRaw <- summary(gam.PhysicoChemRaw) )
legend(0.7,0.9, bty="n", legend=as.expression(paste("Explained deviance of\n",colnames(rotsee.vegan.MOB.table)[z],"\n", 
                                                      rownames(rotsee.mobs.table.t)[z],"\n",
                                                      round(c(summary.gam.PhysicoChemRaw$dev.expl)*100,2),
                                                      "%   ","P=",
                                                      round(c(summary.gam.PhysicoChemRaw$s.pv[[1]]),3))),cex=0.8)
  
  dev.off()
  
}


## make pdf of environmental variable gams fitted on the MOB structure 
##----------------------------------------------------------

sink(file = "GAMs on 16s based MOBs.txt", append = FALSE, type = c("output"),
     split = FALSE)

for (i in c(gamvector)){
  
  
  pdf(paste(colnames(rotsee.env.table[i]), " MOB GAM fitting",".pdf", sep=""),useDingbats=FALSE) #, width=4.5, height=4)
  
  plot(rotsee.mds.allMOBs, type="n")
  col.depth<-(colorRampPalette(brewer.pal(5,"Blues"))(1000))
  points(rotsee.mds.allMOBs, display ="species", cex=3, pch=25, col="orange", bg=alpha(pmmo_color[pMMO_vector_mobs], 0.2), lwd=0)
  points(rotsee.mds.allMOBs, display ="sites", cex=2, alpha= 2, pch=as.numeric(as.factor(rotsee.mobs.env.table$Oxidation_Zone))+20, col=col.depth[as.numeric(rotsee.mobs.env.table$Depth_m)/max(rotsee.mobs.env.table$Depth_m)*1000],
         bg=alpha(col.depth[as.numeric(rotsee.mobs.env.table$Depth_m)/max(rotsee.mobs.env.table$Depth_m)*1000], 0.5))
  text(rotsee.mds.allMOBs$points[, 1]+0.1, rotsee.mds.allMOBs$points[, 2]+0.1,
       labels=unlist(dimnames(rotsee.mds.allMOBs$points)[1]),
       col="darkblue", #col.depth[as.numeric(rotsee.env.table.MOB$Depth_m)/max(rotsee.env.table.MOB$Depth_m)*1000], 
       cex=0.7)
  text(rotsee.mds.allMOBs$species[, 1]*1.1, rotsee.mds.allMOBs$species[, 2]+0.1,
       labels=unlist(dimnames(rotsee.mds.allMOBs$species)[1]),
       col=pmmo_color[pMMO_vector_mobs], cex=0.7)
  text(rotsee.mds.allMOBs$species[,1], rotsee.mds.allMOBs$species[,2]-0.1,
       colnames(rotsee.mobs.table.t), col="darkred", cex=0.7)
  legend(-1,0.9,legend=as.expression(paste("Stress:\n",
                                              round(c(rotsee.mds.allMOBs$stress)*100,2),
                                              "%")),cex=0.8,bty="n")
  image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.3,0.7), axis.args = list(cex.axis = 0.7),
             zlim=range(-(rotsee.env.table$Depth_m)) )
  
 
  
  gam.PhysicoChemRaw<-ordisurf(as.formula(paste("rotsee.mds.allMOBs ~ ",paste(colnames(rotsee.mobs.env.table[i])))),rotsee.mobs.env.table,
                               choices=c(1,2), knots=10,                                         
                               family="gaussian", isotropic= TRUE,
                               col="darkgrey",scaling=3,
                               add = TRUE, display = "sites",
                               nlevels = 20, labcex = 0.8,
                               bubble = TRUE, cex = 1, select = TRUE, method = "REML",
                               gamma = 1, plot = TRUE)
  (summary.gam.PhysicoChemRaw <- summary(gam.PhysicoChemRaw) )
  
  legend(-1,-0.6,legend=as.expression(paste("Explained deviance of\n ", colnames(rotsee.mobs.env.table[i]),
                                           round(c(summary.gam.PhysicoChemRaw$dev.expl)*100,2),
                                           "%,   ","P=",
                                           round(c(summary.gam.PhysicoChemRaw$s.pv[[1]]),3))),cex=0.8,bty="n")
  
  
  dev.off()
  
}

sink()


########################################################################################################################
## vabund with MOBs only (plot mobs which significantly correlate with a specified model)
#########################################################################################################################

rotsee.glmdata <- mvabund(t(prune_taxa(MOBs, otu_table(bac.rotsee.mobs))))  # check which one to use , i.e. bac.rotsee
(env.glmdata<- rotsee.mobs.env.table [,] )  # here it might be good to transform environmental data as performed in the next line
(env.glmdata<- rotsee.env.table.log)

(rotsee.glm<- manyglm(rotsee.glmdata ~ Oxidation_Zone, data=env.glmdata, family="negative.binomial",    #choose (multiple) factors for the model
                      cor.type="I",
                      show.coef=FALSE, show.fitted=TRUE, show.residuals=FALSE ))
#rotsee.glm.summary<-summary(rotsee.glm, resamp="residual", nBoot=20)
(rotsee.glm.anova<- anova(rotsee.glm, p.uni="none",         #anova.manyglm   -> are liminc zones different? yes...
                          show.time="all", nBoot=200)) # for serious run choose 1000 bootstraps!!


summary(rotsee.glm.anova)
rotsee.glm.anova$uni.p # three spp significant after multiple testing 
s = sort(rotsee.glm.anova$uni.test[2,],index.return=T,decreasing=T)
s$x[1:length(MOBs)]/sum(s$x) #the proportions of total test stat due to these 60 OTUs
sum(s$x[1:length(MOBs)]/sum(s$x) )

(glm.imp.names <- names(s$x[1:3]))    # get the 10 bacteria which are most distinct beween tested groups/factor
vabundOTUsnumbers=c()
for (i in c(1:length(glm.imp.names))){
  
  vabundOTUsnumbers[i]<- which(colnames(t(rotsee.mobs.otu.table)) == glm.imp.names[i])
}
vabundOTUsnumbers
## plot them on the nMDS

plot(rotsee.mds.allMOBs, type="n")
col.depth<-(colorRampPalette(brewer.pal(5,"Blues"))(1000))
points(rotsee.mds.allMOBs, display ="species", cex=3, pch=25, col="orange", bg=alpha(pmmo_color[pMMO_vector_mobs], 0.2), lwd=0)
points(rotsee.mds.allMOBs, display ="sites", cex=2, alpha= 2, pch=as.numeric(rotsee.mobs.env.table$Campaign)+20, col=col.depth[as.numeric(rotsee.mobs.env.table$Depth_m)/max(rotsee.mobs.env.table$Depth_m)*1000],
       bg=alpha(col.depth[as.numeric(rotsee.mobs.env.table$Depth_m)/max(rotsee.mobs.env.table$Depth_m)*1000], 0.5))
text(rotsee.mds.allMOBs$points[, 1]+0.1, rotsee.mds.allMOBs$points[, 2]+0.1,
     labels=unlist(dimnames(rotsee.mds.allMOBs$points)[1]),
     col="darkblue", #col.depth[as.numeric(rotsee.env.table.MOB$Depth_m)/max(rotsee.env.table.MOB$Depth_m)*1000], 
     cex=0.7)
text(rotsee.mds.allMOBs$species[, 1]*1.1, rotsee.mds.allMOBs$species[, 2]+0.1,
     labels=unlist(dimnames(rotsee.mds.allMOBs$species)[1]),
     col=pmmo_color[pMMO_vector_mobs], cex=0.7)
legend(-1,0.9,legend=as.expression(paste("Stress:\n",
                                            round(c(rotsee.mds.allMOBs$stress)*100,2),
                                            "%")),cex=0.8,bty="n")
image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.3,0.7), axis.args = list(cex.axis = 0.7),
           zlim=range(-(rotsee.env.table$Depth_m)) )

## plot vabund mobs
points(rotsee.mds.allMOBs$species[vabundOTUsnumbers,1], rotsee.mds.allMOBs$species[vabundOTUsnumbers,2],cex=3, pch=21, col="darkgreen", bg=alpha("green", 0.1))
#text(rotsee.mds.allMOBs$species[vabundOTUsnumbers,1], rotsee.mds.allMOBs$species[vabundOTUsnumbers,2]+0.1,
#     colnames(rotsee.vegan.MOB.table)[vabundOTUsnumbers], cex=0.9, col="darkgreen")

## assign OTU
text(rotsee.mds.allMOBs$species[vabundOTUsnumbers,1], rotsee.mds.allMOBs$species[vabundOTUsnumbers,2]-0.1,
     colnames(t(rotsee.mobs.otu.table))[vabundOTUsnumbers], cex=0.9, col="darkgreen")




########################################################################################################################
########################################################################################################################
## CCA of all OTUs  (better than RDA here as we want to look at possible ecological niches of MOBs within the water column)
########################################################################################################################
########################################################################################################################

    # data standardization and transformations
    
    ClusterHel           <-decostand(rotsee.vegan.otu.table, "hellinger")     # square root of equally standardized row sums
    ClusterNorm          <-decostand(rotsee.vegan.otu.table, method="total")  # standardized to equal row sum
    rotsee.env.table[,-c(1:10,12)]
    (PhysicoChemicalLog <- rotsee.env.table.log[,envvector[-c(1:5,25:33)]])  # log transformed env data
    (PhysicoChemicalLogsub<-rotsee.env.table.log[,envvector[-c(1:5,25:33)]])
    (PhysicoChemicalLogsub=PhysicoChemicalLogsub[colSums(!is.na(PhysicoChemicalLogsub)) >= nrow(PhysicoChemicalLogsub)] )# remove columns with NA's
 
    
    (PhysicoChemicalRaw <- rotsee.env.table[,envvector[-c(1:5,25:33)]]) # original env data
    (PhysicoChemicalRawsub<-PhysicoChemicalRaw[colSums(!is.na(PhysicoChemicalRaw)) >= nrow(PhysicoChemicalRaw)])
    
    (ClusterHelsub <- ClusterHel[rownames(PhysicoChemicalLogsub),])
    
    
    ## remove collinear variables prior to model selection
    
    #stepwise VIF function with preallocated vectors
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
    ## remove elements that with VIFs larger than 20 (you can also skip this step)
    (col.env.vif<- vif_func(in_frame=PhysicoChemicalLogsub,thresh=20,trace=T))
    PhysicoChemicalLogsub = PhysicoChemicalLogsub[,c(col.env.vif)]
    
    
    # remove=c("Mn_DGT", "Mn_diss_nM", "T")
    # (PhysicoChemicalLogsub = PhysicoChemicalLogsub[ , !names(PhysicoChemicalLogsub) %in% remove ])
    
    
    
    ###################################################################################
    ## CCA analysis model assessments (forward, backward and marginal selection)
    ###################################################################################
    ## the selection is run on a sample subset lacking samples having missing values in the environmental matrix
    
    (Cluster.cca<-cca(ClusterHelsub~ ., PhysicoChemicalLogsub, na.action=na.omit))   # we use the log transformed env variables
    
    (Cluster.cca.forward <- ordistep (cca(ClusterHelsub~1, PhysicoChemicalLogsub, na.action=na.omit), # ordiR2step must inherit from rda...cca not working
                                        scope= formula(Cluster.cca),
                                        direction="forward",pstep=100,nperm=99))   ## choose nperm= 999 or higher for final run of the model
    
    (Cluster.cca.backward <- ordistep (cca(ClusterHelsub~., PhysicoChemicalLogsub, na.action=na.omit),
                                       direction="backward",pstep=100,nperm=99))
    
    (Cluster.cca.both <- ordistep (cca(ClusterHelsub~1, PhysicoChemicalLogsub, na.action=na.omit),
                                     scope= formula(Cluster.cca),
                                     direction="both",pstep=100,nperm=99))
    
    ####################################################################################
    ## the best terms for the final model are assessed
    ####################################################################################
    
    (summarized.selection<- c(rownames(Cluster.cca.forward$CCA$biplot),
                              rownames(Cluster.cca.backward$CCA$biplot),
                              rownames(Cluster.cca.both$CCA$biplot)))
    
    (sortsumsel<- (sort(table(summarized.selection))))
    (constrains.freq.table<- as.data.frame (sortsumsel))

    ## to what ever reason you have to use one of the two FORMULA functions... seems to be a R bug.... just uncomment and run the other if you get an error
    ## the same is true for the removement of single terms further below
    
   (FORMULA <-as.formula(paste("ClusterHel ~ ", paste  (rownames((constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=3))])), # here take again complete data set
                                                        collapse="+"),sep = "")))
   
    (FORMULA <-as.formula(paste("ClusterHel ~ ", paste ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=3))]), # here take again complete data set
                                                       collapse="+"),sep = "")))
  
    ###########################################################################################################
    #
    # Run selected CCA model and tests and apriori fit physico-chemical variables not included in
    # the model
    #
    ###########################################################################################################
    
    (Cluster.cca.forwsel <- cca(FORMULA, PhysicoChemicalLog, na.action=na.exclude))
    summary(Cluster.cca.forwsel)
    
    # test the variation inflation factors and significance of models, axes and constraints
    
    sort(vif.cca(Cluster.cca.forwsel))   # remove variable with vif's >20 from the data set as it is highly colinear with other variables
    
    
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
    
    
    (FORMULA    <-as.formula(paste("ClusterHel ~ ", mgsub(c("\\+Fe_tot_nM", "\\+Fe_diss_nM", "\\+Light_uE..m2s", "\\+Depth_m", "\\+Depth_DGT_m" ), c("","","","",""),
                                                           paste (rownames (constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=3))]), # here take again complete data set
                                                         collapse="+")),sep = "")))
    
    (FORMULA    <-as.formula(paste("ClusterHel ~ ", mgsub(c("\\+Fe_tot_nM", "\\+Fe_diss_nM", "\\+Light_uE..m2s", "\\+Depth_m", "\\+Depth_DGT_m" ), c("","","","",""),
                                                          paste ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=2))]), # here take again complete data set
                                                                       collapse="+")), sep = "")))
    
    
    (Cluster.cca.forwsel <- cca(FORMULA, PhysicoChemicalLog, na.action=na.exclude))
    vif.cca(Cluster.cca.forwsel)
    
    ##################----------------------------------------------------------------------------------------------
    (Cluster.cca.best.forwsel <- vegan::cca(FORMULA, PhysicoChemicalLog, na.action=na.exclude))   # this is the final model
    ##################-----------------------------------------------------------------------------------------------

    # (Cluster.cca.best.forwsel <- cca(ClusterHel ~ Cu_DGT_nM + Condition(Depth_m), PhysicoChemicalLog, na.action=na.exclude))   # alternatively make model with specific variables
    
    
    summary(Cluster.cca.best.forwsel)
    
    anova(Cluster.cca.best.forwsel) # check if total model is singificant
    
    anova(Cluster.cca.best.forwsel, by="axis", permutations = how(nperm=999),  # check for significant axes
          model = c("direct"),
          strata = NULL,
          cutoff = 1)
    
    ## check for significance of environmental terms (all na data has to be removed!)
    # (FORMULAsub    <-as.formula(paste("ClusterHelsub ~ ", mgsub(c("\\+Fe_tot_nM"),
    #                                                             c(""), paste (rownames (constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=3))]), # here take again complete data set
    #                                                                               collapse="+")),sep = "")))
    # 
    # (FORMULAsub    <-as.formula(paste("ClusterHelsub ~ ", mgsub(c("\\+Fe_tot_nM"), c(""),
    #                                                             paste ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=2))]), # here take again complete data set
    #                                                                         collapse="+")),sep = "")))
    # 
    # 
    FORMULAsub <- FORMULA
    
    (Cluster.cca.best.forwsel.sub <- cca(FORMULAsub, PhysicoChemicalLogsub, na.action=na.omit))  
    anova(Cluster.cca.best.forwsel.sub, by="terms", permutations = how(nperm=999),
          model = c("direct"), 
          parallel = getOption("mc.cores"), strata = NULL,
          scope = NULL)
    
    
    ## fit the non-used environmental variables in the RDA model
    (best.model.env.selection <-c(rownames(Cluster.cca.best.forwsel$CCA$biplot)))
    (rest.model.env.selection <-PhysicoChemicalRaw[! c(colnames(PhysicoChemicalRaw)) %in% best.model.env.selection])
    
    (envfit.rest  <- envfit(Cluster.cca.best.forwsel, rest.model.env.selection[,1:length(rest.model.env.selection),drop=FALSE], na.rm=TRUE))

        ## fit the species in the rda model
    #(speciesfit   <- envfit(Cluster.cca.best.forwsel, rotsee.vegan.otu.table))
    
    
    ## and some additional diagnostics
    ##-------------------------------------------
    # finds the so-called "interset correlation" or (weighted) correlation of weighted averages scores and constraints.
    # The defined contrasts are used for factor variables. This is a bad measure since it is a correlation. Further,
    # it focuses on correlations between single contrasts and single axes instead of looking at the multivariate relationship.
    # Fitted vectors (envfit) provide a better alternative. Biplot scores (see scores.cca) are a multivariate alternative for
    # (weighted) correlation between linear combination scores and constraints. 
    intersetcor(Cluster.cca.best.forwsel)
    
    dispvect<-c("sp","wa","cn")
    for (i in c(1:3)){
      print(scores(Cluster.cca.best.forwsel, choices = c(1,2), display =dispvect[i] ,
                   scaling = "symmetric", correlation = TRUE))
    }
    
    # finds the so-called "species - environment correlation" or (weighted) correlation of weighted average scores and
    # linear combination scores. This is a bad measure of goodness of ordination, because it is sensitive to extreme
    # scores (like correlations are), and very sensitive to overfitting or using too many constraints. Better models
    # often have poorer correlations. Function ordispider can show the same graphically.
    spenvcor(Cluster.cca.best.forwsel)
    
    #Calculating the proportion of variance explained for an ordination graph
    goodness(Cluster.cca.best.forwsel, display = c("species"), 
             model = c("CCA", "CA"), statistic = c("explained" ), #"distance"      
             summarize = TRUE, addprevious = FALSE)
    
    # Calculating the variance of each species of the species matrix
    # decomposes the inertia into partial, constrained and unconstrained components
    # for each site or species. Instead of inertia, the function can give the total dispersion
    # or distances from the centroid for each component.                                                                            
    inertcomp(Cluster.cca.best.forwsel, display = c("species"),                     
              unity = FALSE, proportional = FALSE) # "explained"         
    
    
    ##---------------------------------------------------------
    # run the percentage contribution of constrains function 
    ##---------------------------------------------------------
    
    
    ################################################################################################
    #
    #  PerContCon :  Function to assess the importance (% variance explained)
    #                of the environmental parameters in direction of the different canonical axes
    #
    ################################################################################################
    PerContCon<- function(RDAmodel) {                # RDAmodel is an CCA vegan output)     
      
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
      
      writeLines("% explained by environmental variables\n")
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
    
    PerContCon(Cluster.cca.best.forwsel) # here you can see which parameters have higher influence on which canonical axes. 
    # one can plot the specific axes according to structures one whishes to inverstigate
    
    # check R2 of apriori fitted constraints
    (envfit.apriori  <- envfit(Cluster.cca.best.forwsel, PhysicoChemicalRawsub, na.rm=TRUE, choices=c(1,2))) # you can check with the above output
    # which dimensions might be good for a specific model parameter for
    # fitting on and plotting of the CCA constraints
    
    ##########################################################################################################
    # Bi-plot of the CCA  model  (all following plots will be directly saved as PDF's when pdf() is decommented)
    ##########################################################################################################
    my.color1<-c("cyan", "cyan4", "cyan3")
    MOBsOTUsnumbers
    
    #pdf("RDA of Clusters with Env as constraint.pdf",useDingbats=FALSE)
    
    # Basic Plot

    ordiplot<-plot(Cluster.cca.best.forwsel,type="none",choices = c(xx<-1, yy<-2),scaling=3,main="", xlab="",ylab="")     # here you can change the axes to be ploted       
    
    points(ordiplot,"species", pch=21,cex=normalized.species.size.vec,
           col=col.genus[genus.df$id], bg=alpha(col.genus[genus.df$id], 0.99), lwd=0)
    legend("topright",c(unique(rotsee.tax.table[,"Phylum"])),pch=c(21),col=alpha(unique(col.genus[genus.df$id]), 0.8), pt.bg=alpha(unique(col.genus[genus.df$id]),0.5), bty="n", cex=0.5)
    points(ordiplot,"sites", cex=3, alpha= 2, pch=as.numeric(rotsee.env.table$Campaign)+20, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
           bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.5))
    
    text(ordiplot$sites[, 1]*1.3, ordiplot$sites[, 2]*1.3,
         labels=unlist(dimnames(ordiplot$sites)[1]),
         col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],cex=0.7)
    
    image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.3,0.7), axis.args = list(cex.axis = 0.7),
               zlim=range(-(rotsee.env.table$Depth_m)) )
    
    # Confidence ellipses for the Locations (SE +-95%)
    
    (groupz<-sort(unique(rotsee.env.table.log$Oxidation_Zone)))
    with(rotsee.env.table.log,ordiellipse(ordiplot, Oxidation_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[1], show.groups=groupz[1],draw="polygon", alpha=30))
    with(rotsee.env.table.log,ordiellipse(ordiplot, Oxidation_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[1], show.groups=groupz[2],draw="polygon", alpha=30))
    with(rotsee.env.table.log,ordiellipse(ordiplot, Oxidation_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[1], show.groups=groupz[3],draw="polygon", alpha=30))
    with(rotsee.env.table.log,ordiellipse(ordiplot, Oxidation_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[1], show.groups=groupz[4],draw="polygon", alpha=30))
    
    # plot the significant environmental variables from the selected canonical model
    
    (curr.plot.multiplier <- vegan:::ordiArrowMul(scores(Cluster.cca.best.forwsel$CCA$biplot[,c(xx,yy)]),
                                                  at = c(0,0), fill = 1,display="sites", choices = c(1,2)))# multiplier to fit the current plot ratio
    arrows(0, 0, Cluster.cca.best.forwsel$CCA$biplot[, 1]*curr.plot.multiplier*0.7, Cluster.cca.best.forwsel$CCA$biplot[, 2]*curr.plot.multiplier*0.7, lty=1,lwd=1.5,
           length = 0.15, angle = 30,  col="darkgrey")
    
    
    text(Cluster.cca.best.forwsel$CCA$biplot[, 1]*curr.plot.multiplier*0.9, Cluster.cca.best.forwsel$CCA$biplot[, 2]*curr.plot.multiplier*0.9,
         labels=unlist(dimnames(Cluster.cca.best.forwsel$CCA$biplot)[1]), col="darkgrey", cex=0.8)
    
    
    
    # Optionally: plot the rest of the environmental variables if you wish to....
    (envfit.rest  <- envfit(Cluster.cca.best.forwsel, rest.model.env.selection[,1:ncol(rest.model.env.selection)], choices=c(xx,yy),na.rm=TRUE))
    pval.limit <- which(c(envfit.rest$vectors$pvals)<=1) # set max p-value 
    
    (envfit.vectors       <- scores(envfit.rest,"vectors")[c(pval.limit),1:2])  # these are envfit vectors scaled to each one another (i.e. multipied by square root of column r2)
    (env.loading.factor   <- sqrt((envfit.rest$vectors$r)[c(pval.limit)]))      # loading factor for the envfit scores to scale them according to their correlations
    (curr.plot.multiplier <- vegan:::ordiArrowMul(scores(envfit.rest,"vectors")[c(pval.limit),1:2],
                                                  at = c(0,0), fill = 1,display="sites", choices = c(1,2)))# multiplier to fit the current plot ratio
    env.arrows.x  <- as.data.frame(scores(envfit.rest,"vectors"))[c(pval.limit),1]
    env.arrows.y  <- as.data.frame(scores(envfit.rest,"vectors"))[c(pval.limit),2]
    env.names     <-rownames(as.data.frame(scores(envfit.rest,"vectors")))[c(pval.limit)]
    
    arrows(0, 0, env.arrows.x*curr.plot.multiplier*0.55, env.arrows.y*curr.plot.multiplier*0.55, lty=1, col="lightgrey",
           length = 0.15, angle = 30)
    text(env.arrows.x*curr.plot.multiplier*0.72, env.arrows.y*curr.plot.multiplier*0.72,
         labels=env.names, col="lightgrey")    
    
    
    
    # plot MOB centroids

   points(ordiplot$species[MOBsOTUsnumbers,1], ordiplot$species[MOBsOTUsnumbers,2],
          cex=2, pch=25, col="orange", bg=alpha(pmmo_color[pMMO_vector], 0.3))
   
   text(ordiplot$species[MOBsOTUsnumbers,1]+0, ordiplot$species[MOBsOTUsnumbers,2]+0.1,
        attributes(ordiplot$species)[[2]][[1]][MOBsOTUsnumbers],
        col=pmmo_color[pMMO_vector],
        cex=0.7)
   
   text(ordiplot$species[MOBsOTUsnumbers,1]+0, ordiplot$species[MOBsOTUsnumbers,2]-0.1,
        colnames(t(rotsee.otu.table))[MOBsOTUsnumbers],
        col="darkred",
        cex=0.7)
    
    # plot legend and axis with explained variance
    
    legend("topleft",pch=c(21,25),col=alpha(c("blue","darkorange"),0.1),pt.bg=alpha(c("blue","green"),0.5),
           legend=c("Samples","OTUs"),bty="n")
    title(xlab=as.expression(paste("CCA",xx," (",round(c(100/(Cluster.cca.best.forwsel$tot.chi)*(Cluster.cca.best.forwsel$CCA$eig[xx])),2),"%)")),
          ylab=as.expression(paste("CCA",yy," (",round(c(100/(Cluster.cca.best.forwsel$tot.chi)*(Cluster.cca.best.forwsel$CCA$eig[yy])),2),"%)")))
    
    
    #dev.off()
    
    
    ## you might want to check the linearity of important variables in the cca
    
    gam.PhysicoChemRaw<-ordisurf(as.formula(paste("Cluster.cca.best.forwsel ~ ",paste(colnames(PhysicoChemicalRaw[parame<-"Cu_DGT_nM"])))), PhysicoChemicalRaw,
                                 choices=c(1,2), knots=10,                                         
                                 family="gaussian", isotropic= TRUE,
                                 col="darkgrey",scaling=3,
                                 add = TRUE, display = "sites",
                                 nlevels = 20, labcex = 0.8,
                                 bubble = TRUE, cex = 1, select = TRUE, method = "REML",
                                 gamma = 1, plot = TRUE)
    (summary.gam.PhysicoChemRaw <- summary(gam.PhysicoChemRaw) )
    
    legend(1,-2, bty="n", legend=as.expression(paste("Explained deviance of\n",parame,":\n",   # change according to what you choose
                                                    round(c(summary.gam.PhysicoChemRaw$dev.expl)*100,2),
                                                    "%   ","P=",
                                                    round(c(summary.gam.PhysicoChemRaw$s.pv[[1]]),3))),cex=0.8)
    
    
    
    
    
    ## here you can interactively identify potential bacterial interactors
    identify(ordiplot$species, labels=rownames(ordiplot$species), cex=0.6, col="grey" )
    identify(ordiplot$species, labels=rownames(rotsee.otu.table), cex=0.6, col="grey" )
    ## add points of species found in the OTU list
   ( bacterianumber <- c(which(rownames(ordiplot$species)==("Planktothrix"))))
    points(ordiplot$species[bacterianumber,1], ordiplot$species[bacterianumber,2],
           cex=2, pch=25,col="black", bg=alpha("yellow", alpha=0.4))
    text(ordiplot$species[bacterianumber,1], ordiplot$species[bacterianumber,2],
         labels=rownames(ordiplot$species)[bacterianumber], cex=0.8)
    
    
    ## generally you see a horseshoe effect of the CCA
    
    
    ########################################################################################################################
    ###########################################################################################################################
    ## distance-based RDA with all OTUs
    ##########################################################################################################################
    ########################################################################################################################
    
    # data standardization and transformations
    
    ClusterHel           <-decostand(rotsee.vegan.otu.table, "hellinger")     # square root of equally standardized row sums
    ClusterNorm          <-decostand(rotsee.vegan.otu.table, method="total")  # standardized to equal row sum
    
    (PhysicoChemicalLog <- rotsee.env.table.log[,envvector[-c(1:5,25:33)]])  # log transformed env data
    (PhysicoChemicalLogsub<-rotsee.env.table.log[,envvector[-c(1:5,25:33)]])
    (PhysicoChemicalLogsub=PhysicoChemicalLogsub[colSums(!is.na(PhysicoChemicalLogsub)) >= nrow(PhysicoChemicalLogsub)] )# remove columns with NA's
    
    
    (PhysicoChemicalRaw <- rotsee.env.table[,envvector[-c(1:5,25:33)]]) # original env data
    (PhysicoChemicalRawsub<-PhysicoChemicalRaw[colSums(!is.na(PhysicoChemicalRaw)) >= nrow(PhysicoChemicalRaw)])
    
    (col.env.vif<- vif_func(in_frame=PhysicoChemicalLogsub,thresh=20,trace=T))
    PhysicoChemicalLogsub = PhysicoChemicalLogsub[,c(col.env.vif)]
    
    (ClusterHelsub <- ClusterHel[rownames(PhysicoChemicalLogsub),])
    
    
    ###################################################################################
    ## model assessments (forward, backward and marginal selection)
    ###################################################################################
    ## the selection is run on a sample subset lacking samples having missing values in the environmental matrix
    
    (Cluster.capscale<-capscale(ClusterHel~ ., PhysicoChemicalLogsub , distance="bray", metaMDSdist=FALSE, comm=ClusterHel, sqrt.dist = TRUE, na.action=na.omit))   # we use the log transformed env variables
    
    
    (Cluster.capscale.forward <- ordistep (capscale(ClusterHelsub~1, PhysicoChemicalLogsub, distance="bray", metaMDSdist=FALSE, sqrt.dist = TRUE, comm=ClusterHel, na.action=na.omit), # ordiR2step must inherit from rda...cca not working
                                           scope= formula(Cluster.capscale),
                                           direction="forward",pstep=100,nperm=9))
    
    (Cluster.capscale.backward <- ordistep (capscale(ClusterHelsub~., PhysicoChemicalLogsub, distance="bray", metaMDSdist=FALSE, sqrt.dist = TRUE, comm=ClusterHel, na.action=na.omit),
                                            direction="backward",pstep=100,nperm=9))
    
    (Cluster.capscale.both <- ordistep (capscale(ClusterHelsub~1, PhysicoChemicalLogsub, distance="bray", metaMDSdist=FALSE, sqrt.dist = TRUE, comm=ClusterHel, na.action=na.omit),
                                        scope= formula(Cluster.capscale),
                                        direction="both",pstep=100,nperm=9))
    
    ####################################################################################
    ## the best terms for the final model are assessed
    ####################################################################################
    
    (summarized.selection<- c(rownames(Cluster.capscale.forward$CCA$biplot),
                              rownames(Cluster.capscale.backward$CCA$biplot),
                              rownames(Cluster.capscale.both$CCA$biplot)))
    
    (sortsumsel<- (sort(table(summarized.selection))))
    (constrains.freq.table<- as.data.frame (sortsumsel))
    
    ## to what ever reason you have to use one of the two FORMULA functions... seems to be a R bug.... just uncomment and run the other if you get an error
    ## the same is true for the removement of single terms further below
    
    # (FORMULA <-as.formula(paste("ClusterHel ~ ", paste  (rownames((constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=3))])), # here take again complete data set
    #                                                      collapse="+"),sep = "")))
    
    (FORMULA <-as.formula(paste("ClusterHel ~ ", paste ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=2))]), # here take again complete data set
                                                         collapse="+"),sep = "")))
    
    
    ###########################################################################################################
    #
    # Run selected dbRDA model and tests and apriori fit physico-chemical variables not included in
    # the model
    #
    ###########################################################################################################
    
    (Cluster.capscale.forwsel <- capscale(FORMULA, PhysicoChemicalLog, distance="bray", metaMDSdist=FALSE, comm=ClusterHel, sqrt.dist = TRUE, na.action=na.exclude))
    summary(Cluster.capscale.forwsel)
    
    # test the variation inflation factors and significance of models, axes and constraints
    
    sort(vif.cca(Cluster.capscale.forwsel))   # remove variables from the data set that are highly colinear 
    
    # (FORMULA    <-as.formula(paste("ClusterHel ~ ", mgsub(c("\\+Fe_tot_nM"),
    #                                                       c(""), paste  (rownames(constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=2))]), # here take again complete data set
    #                                                                      collapse="+")),sep = "")))
    
    (FORMULA    <-as.formula(paste("ClusterHel ~ ", mgsub(c("\\+Fe_tot_nM", "\\+Fe_diss_nM", "\\+Light_uE..m2s", "\\+Depth_m", "\\+Depth_DGT_m" ), c("","","","",""),
                                                          paste ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=2))]), # here take again complete data set
                                                                 collapse="+")),sep = "")))
    
    
    (Cluster.capscale.best.forwsel <- capscale(FORMULA, PhysicoChemicalLog, distance="bray",comm=ClusterHel, metaMDSdist=FALSE,  sqrt.dist = TRUE, na.action=na.exclude))   # this is the final model
    summary(Cluster.capscale.best.forwsel)
    vif.cca(Cluster.capscale.best.forwsel)
    
    anova(Cluster.capscale.best.forwsel) # check if total model is singificant
    
    anova(Cluster.capscale.best.forwsel, by="axis", permutations = how(nperm=999),  # check for significant axes
          model = c("direct"),
          strata = NULL,
          cutoff = 1)
    
    
    # (FORMULAsub    <-as.formula(paste("ClusterHel ~ ", mgsub(c("\\+Fe_tot_nM"),
    #                                                             c(""), paste  (rownames(constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=1))]), # here take again complete data set
    #                                                                            collapse="+")),sep = "")))
    # 
    # (FORMULA    <-as.formula(paste("ClusterHel ~ ", mgsub(c("\\\\+Fe_tot_nM"),
    #                                                      c(""), paste ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=2))]), # here take again complete data set
    #                                                              collapse="+")),sep = "")))
    
    (Cluster.capscale.best.forwsel.sub <- capscale(FORMULA, PhysicoChemicalLogsub, distance="bray", metaMDSdist=FALSE, comm=ClusterHel, sqrt.dist = TRUE, na.action=na.omit))  # check for significance of environmental terms
    anova(Cluster.capscale.best.forwsel.sub, by="terms", permutations = how(nperm=999),
          model = c("direct"), 
          parallel = getOption("mc.cores"), strata = NULL,
          scope = NULL)
    
    
    # fit the non-used environmental variables in the RDA model
    
    (best.model.env.selection <-c(rownames(Cluster.capscale.best.forwsel$CCA$biplot)))
    (rest.model.env.selection <-PhysicoChemicalRaw[! c(colnames(PhysicoChemicalRaw)) %in% best.model.env.selection])
    
    (envfit.rest  <- envfit(Cluster.capscale.best.forwsel, rest.model.env.selection[,1:length(rest.model.env.selection)], na.rm=TRUE))
    
    
    # and some additional diagnostics
    
    intersetcor(Cluster.capscale.best.forwsel)
    
    dispvect<-c("sp","wa","cn")
    for (i in c(1:3)){
      print(scores(Cluster.capscale.best.forwsel, choices = c(1,2), display =dispvect[i] ,
                   scaling = "symmetric", correlation = TRUE))
    }
    
    screeplot(Cluster.capscale.best.forwsel, bstick = FALSE, type = c("lines"), #"barplots"
              npcs = min(10, if (is.null(Cluster.capscale.best.forwsel$CCA)) Cluster.capscale.best.forwsel$CA$rank else Cluster.capscale.best.forwsel$CCA$rank),
              ptype = "o", bst.col = "red", bst.lty = "solid",
              xlab = "Component", ylab = "Inertia")
    
    spenvcor(Cluster.capscale.best.forwsel)
    
    # goodness(Cluster.capscale.best.forwsel, display = c("sites"), 
    #          model = c("CCA", "CA"), statistic = c("explained", "distance"),
    #          summarize = TRUE, addprevious = FALSE)
    
    inertcomp(Cluster.capscale.best.forwsel, display = c("sites"),
              unity = FALSE, proportional = FALSE)
    
    # run the percentage contribution of constrains function (defined at at the end of the script)
    
    PerContCon(Cluster.capscale.best.forwsel)          # you can plot the axes which show most explained variance by the variable of interest
    
    # check R2 of apriori fitted constraints
    (envfit.apriori  <- envfit(Cluster.capscale.best.forwsel, PhysicoChemicalRawsub, na.rm=TRUE, choices=c(1,2)))  # choose the axes according to the most explained fitting axes (PerContCon output)
    # to check specific variable of interest (i.e. Cu_DGT axes 2,3 r2 is high
    # whereas axes 1,2 r2 is low...)
    
    ##########################################################################################################
    # Bi-plot of the CCA model          (all following plots will be directly saved as PDF's)
    ##########################################################################################################
    my.color1<-c("cyan", "cyan4", "cyan3")
    MOBsOTUsnumbers
    
    #pdf("RDA of Clusters with Env as constraint.pdf",useDingbats=FALSE)
    
    # Basic Plot
    
    ordiplot<-plot(Cluster.capscale.best.forwsel,type="none",choices = c(xx<-1, yy<-2),scaling=3,main="", xlab="",ylab="")     # here you can change the axes to be ploted       
    
    points(ordiplot,"species", pch=21,cex=normalized.species.size.vec,
           col=col.genus[genus.df$id], bg=alpha(col.genus[genus.df$id], 0.99), lwd=0)
    legend("topright", inset=c(-0.1,0),c(unique(rotsee.tax.table[,"Phylum"])),pch=c(19),col=alpha(unique(col.genus[genus.df$id]), 0.5), bty="n", cex=0.6)
    points(ordiplot,"sites", cex=3, alpha= 2, pch=as.numeric(rotsee.env.table$Campaign)+20, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
           bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.5))
    
    text(ordiplot$sites[, 1]*1.3, ordiplot$sites[, 2]*1.3,
         labels=unlist(dimnames(ordiplot$sites)[1]),
         col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],cex=0.7)
    
    image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.3,0.7), axis.args = list(cex.axis = 0.7),
               zlim=range(-(rotsee.env.table$Depth_m)) )
    
    # Confidence ellipses for the Locations (SE +-95%)
    
    (groupz<-sort(unique(rotsee.env.table.log$Oxidation_Zone)))
    with(rotsee.env.table.log,ordiellipse(ordiplot, Oxidation_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[1], show.groups=groupz[1],draw="polygon", alpha=30))
    with(rotsee.env.table.log,ordiellipse(ordiplot, Oxidation_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[1], show.groups=groupz[2],draw="polygon", alpha=30))
    with(rotsee.env.table.log,ordiellipse(ordiplot, Oxidation_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[1], show.groups=groupz[3],draw="polygon", alpha=30))
    with(rotsee.env.table.log,ordiellipse(ordiplot, Oxidation_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[1], show.groups=groupz[4],draw="polygon", alpha=30))
    
    # plot the significant environmental variables from the selected canonical model
    
    (curr.plot.multiplier <- vegan:::ordiArrowMul(scores(Cluster.capscale.best.forwsel$CCA$biplot[,c(xx,yy)]),
                                                  at = c(0,0), fill = 1,display="sites", choices = c(xx,yy)))# multiplier to fit the current plot ratio
    arrows(0, 0, Cluster.capscale.best.forwsel$CCA$biplot[, xx]*curr.plot.multiplier*0.7, Cluster.capscale.best.forwsel$CCA$biplot[, yy]*curr.plot.multiplier*0.7, lty=1,lwd=1.5,
           length = 0.15, angle = 30,  col="darkgrey")
    
    
    text(Cluster.capscale.best.forwsel$CCA$biplot[, xx]*curr.plot.multiplier*0.9, Cluster.capscale.best.forwsel$CCA$biplot[, yy]*curr.plot.multiplier*0.9,
         labels=unlist(dimnames(Cluster.capscale.best.forwsel$CCA$biplot)[1]), col="darkgrey", cex=0.8)
    
    
    
    # Optionally: plot the rest of the environmental variables if you wish to....
    
    (envfit.rest  <- envfit(Cluster.capscale.best.forwsel, rest.model.env.selection[,1:length(rest.model.env.selection)], na.rm=TRUE,  choices=c(xx,yy)))
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
    
    
    
    # plot MOB centroids
    
    points(ordiplot$species[MOBsOTUsnumbers,1], ordiplot$species[MOBsOTUsnumbers,2],
           cex=2, pch=25, col="orange", bg=alpha(pmmo_color[pMMO_vector], 0.3))
    
    text(ordiplot$species[MOBsOTUsnumbers,1]+0, ordiplot$species[MOBsOTUsnumbers,2]+0.1,
         #colnames(rotsee.vegan.otu.table)[MOBsOTUsnumbers],
         attributes(ordiplot$species)[[2]][[1]][MOBsOTUsnumbers],
         col=pmmo_color[pMMO_vector],
         cex=0.7)
    
    text(ordiplot$species[MOBsOTUsnumbers,1]+0, ordiplot$species[MOBsOTUsnumbers,2]-0.1,
         colnames(t(rotsee.otu.table))[MOBsOTUsnumbers],
         col="darkred", cex=0.7)
    
    # plot legend and axis with explained variance
    
    legend("bottomleft",pch=c(21, 25),col=alpha(c("blue","red"),0.5),pt.bg=alpha(c("blue","red"),0.4),
           legend=c("Samples","OTUs"),bty="n")
    title(xlab=as.expression(paste("CAP",xx," (",round(c(100/(Cluster.capscale.best.forwsel$tot.chi)*(Cluster.capscale.best.forwsel$CCA$eig[xx])),2),"%)")),
          ylab=as.expression(paste("CAP",yy," (",round(c(100/(Cluster.capscale.best.forwsel$tot.chi)*(Cluster.capscale.best.forwsel$CCA$eig[yy])),2),"%)")))
    
    
    #dev.off()
    
    
    ## you might want to check the linearity of important variables in the cca
    
    
    gam.PhysicoChemRaw<-ordisurf(as.formula(paste("Cluster.capscale.best.forwsel ~ ",paste(colnames(PhysicoChemicalRaw[parame<-"Cu_DGT_nM"])))), PhysicoChemicalRaw,
                                 choices=c(xx,yy), knots=10,                                         
                                 family="gaussian", isotropic= TRUE,
                                 col="grey",scaling=3,
                                 add = TRUE, display = "sites",
                                 nlevels = 15, labcex = 0.8,
                                 bubble = TRUE, cex = 1, select = TRUE, method = "REML",lwd=0.2,
                                 gamma = 1, plot = TRUE)
    (summary.gam.PhysicoChemRaw <- summary(gam.PhysicoChemRaw) )
    
    legend(0.4,-0.6, bty="n", legend=as.expression(paste("Explained deviance of\n",parame,":\n",   # change according to what you choose
                                                         round(c(summary.gam.PhysicoChemRaw$dev.expl)*100,2),
                                                         "%   ","P=",
                                                         round(c(summary.gam.PhysicoChemRaw$s.pv[[1]]),3))),cex=0.8)
    
    
    
    
    
    
    
    ########################################################################################################################
    ####################################################################################
    ## CCA with the MOBs only data set
    ####################################################################################
    ########################################################################################################################
    # see at the NMDS ... rotsee.vegan.MOB.table can be choosen according to standardization procedure
    

    
    # data standardization and transformations
    
    ClusterHel           <-decostand(rotsee.vegan.mobs.otu.table, "hellinger")     # square root of equally standardized row sums
    ClusterNorm          <-decostand(rotsee.vegan.mobs.otu.table, method="total")  # standardized to equal row sum
    
    (PhysicoChemicalLog <- rotsee.mobs.env.table.log[,envvector[-c(1:5,25:33)]])  # log transformed env data
    (PhysicoChemicalLogsub<-rotsee.mobs.env.table.log[,envvector][-c(1:5,25:33)])
    (PhysicoChemicalLogsub=PhysicoChemicalLogsub[colSums(!is.na(PhysicoChemicalLogsub)) >= nrow(PhysicoChemicalLogsub)] )# remove columns with NA's
    
    
    (PhysicoChemicalRaw <- rotsee.mobs.env.table[,envvector[-c(1:5,25:33)]]) # original env data
    (PhysicoChemicalRawsub<-PhysicoChemicalRaw[colSums(!is.na(PhysicoChemicalRaw)) >= nrow(PhysicoChemicalRaw)])
    
    (col.env.vif<- vif_func(in_frame=PhysicoChemicalLogsub,thresh=20,trace=T))
    PhysicoChemicalLogsub = PhysicoChemicalLogsub[,c(col.env.vif)]
    
    (ClusterHelsub <- ClusterHel[rownames(PhysicoChemicalLogsub),])
    
    
    
    
    ###################################################################################
    ## CCA analysis model assessments (forward, backward and marginal selection)
    ###################################################################################
    ## the selection is run on a sample subset lacking samples having missing values in the environmental matrix
    
    (Cluster.rda<-cca(ClusterHelsub~ ., PhysicoChemicalLogsub, na.action=na.omit))   # we use the log transformed env variables
    
    (Cluster.rda.forward <- ordistep (cca(ClusterHelsub~1, PhysicoChemicalLogsub, na.action=na.omit), # ordiR2step must inherit from rda...cca not working
                                      scope= formula(Cluster.rda),
                                      direction="forward",pstep=100,nperm=99))
    
    (Cluster.rda.backward <- ordistep (cca(ClusterHelsub~., PhysicoChemicalLogsub, na.action=na.omit),
                                       direction="backward",pstep=100,nperm=99))
    
    (Cluster.rda.both <- ordistep (cca(ClusterHelsub~1, PhysicoChemicalLogsub, na.action=na.omit),
                                   scope= formula(Cluster.rda),
                                   direction="both",pstep=100,nperm=99))
    
    ####################################################################################
    ## the best terms for the final model are assessed
    ####################################################################################

    
    
    (summarized.selection<- c(rownames(Cluster.rda.forward$CCA$biplot),
                              rownames(Cluster.rda.backward$CCA$biplot),
                              rownames(Cluster.rda.both$CCA$biplot)))
    
    (sortsumsel<- (sort(table(summarized.selection))))
    (constrains.freq.table<- as.data.frame (sortsumsel))
    
    
    # (FORMULA <-as.formula(paste("ClusterHel ~ ", paste  (rownames((constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=3))])), # here take again complete data set
    #                                                     collapse="+"),sep = "")))
    
    (FORMULA <-as.formula(paste("ClusterHel ~ ", paste  ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=3))]), # here take again complete data set
                                                       collapse="+"),sep = "")))
    
    
    
    
    ###########################################################################################################
    #
    # Run selected CCA model and tests and apriori fit physico-chemical variables not included in
    # the model
    #
    ###########################################################################################################
    
    (Cluster.rda.forwsel <- cca(FORMULA, PhysicoChemicalLog, na.action=na.exclude))
    summary(Cluster.rda.forwsel)
    
    ## test the variation inflation factors and significance of models, axes and constraints
    
    sort(vif.cca(Cluster.rda.forwsel))   # remove Fe_tot_nM from the data set as it is highly colinear with the the dissolved species

    (FORMULA    <-as.formula(paste("ClusterHel ~ ", mgsub(c("\\+Fe_diss_nM" ), c(""),
                                                          paste  ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=3))]), # here take again complete data set
                                                                        collapse="+")),sep = "")))
    
    #(FORMULA    <-as.formula(paste("ClusterHel ~ ", mgsub(c("\\+Light_uE..m2s"), c(""), 
    #                                                      paste  (rownames((constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=3))])), # here take again complete data set
    #                                                      collapse="+")),sep = "")))
    
    (Cluster.rda.best.forwsel <- cca(FORMULA, PhysicoChemicalLog, na.action=na.exclude))   # this is the final model
    summary(Cluster.rda.best.forwsel)
    vif.cca(Cluster.rda.best.forwsel)
    
    anova(Cluster.rda.best.forwsel) # check if total model is singificant
    
    anova(Cluster.rda.best.forwsel, by="axis", permutations = how(nperm=999),  # check for significant axes
          model = c("direct"),
          strata = NULL,
          cutoff = 1)
    

    FORMULAsub <- FORMULA
    # (FORMULAsub    <-as.formula(paste("ClusterHelsub ~ ", mgsub(c("\\+Light_uE..m2s"), c(""), 
    #                                                         paste  ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=3))]), # here take again complete data set
    #                                                                collapse="+")),sep = "")))
    
    #(FORMULAsub    <-as.formula(paste("ClusterHelsub ~ ", mgsub(c("\\+Light_uE..m2s"), c(""), 
    #                                                      paste  (rownames((constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=3))])), # here take again complete data set
    #                                                      collapse="+")),sep = "")))
    
    (Cluster.rda.best.forwsel.sub <- rda(FORMULAsub, PhysicoChemicalLogsub, na.action=na.omit))  # check for significance of environmental terms
    anova(Cluster.rda.best.forwsel.sub, by="terms", permutations = how(nperm=999),
          model = c("direct"), 
          parallel = getOption("mc.cores"), strata = NULL,
          scope = NULL)
    
    
    # fit the non-used environmental variables in the RDA model
    
    (best.model.env.selection <-c(rownames(Cluster.rda.best.forwsel$CCA$biplot)))
    (rest.model.env.selection <-PhysicoChemicalRaw[! c(colnames(PhysicoChemicalRaw)) %in% best.model.env.selection])
    
    (envfit.rest  <- envfit(Cluster.rda.best.forwsel, rest.model.env.selection[,1:length(rest.model.env.selection)], na.rm=TRUE))

    
    # and some additional diagnostics
    
    intersetcor(Cluster.rda.best.forwsel)
    
    spenvcor(Cluster.rda.best.forwsel)
    
    goodness(Cluster.rda.best.forwsel, display = c("species"), 
             model = c("CCA", "CA"), statistic = c("explained", "distance"),
             summarize = TRUE, addprevious = FALSE)
    inertcomp(Cluster.rda.best.forwsel, display = c("species"),
              unity = FALSE, proportional = FALSE)
    
    # run the percentage contribution of constrains function (defined at at the end of the script)... 
    PerContCon(Cluster.rda.best.forwsel)
    # check R2 of apriori fitted constraints
    (envfit.apriori  <- envfit(Cluster.rda.best.forwsel, PhysicoChemicalRawsub, na.rm=TRUE, choices=c(1,2)))
    
    ##########################################################################################################
    # Bi-plot of the RDA model          (all following plots will be directly saved as PDF's)
    ##########################################################################################################
    my.color1<-c("cyan", "cyan4", "cyan3")
    MOBsOTUsnumbers
    
    #pdf("CCA of Clusters with Env as constraint.pdf",useDingbats=FALSE)
    
    # Basic Plot
    ordiplot<-plot(Cluster.rda.best.forwsel,type="none",choices = c(xx<-1, yy<-2),scaling=3,main="", xlab="",ylab="")     # here you can change the axes to be ploted       
    

        points(ordiplot,"species", pch=25,cex=normalized.species.size.vec * 2,
              col=alpha(pmmo_color[pMMO_vector_mobs], 0.8), bg=alpha(pmmo_color[pMMO_vector_mobs],0.3))
        
        points(ordiplot,"sites", cex=3, alpha= 2, pch=as.numeric(as.factor(rotsee.env.table$Oxidation_Zone))+20, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
               bg=alpha(col.depth[as.numeric(rotsee.mobs.env.table$Depth_m)/max(rotsee.mobs.env.table$Depth_m)*1000], 0.5))
        
        text(ordiplot$sites[, 1]*1.3, ordiplot$sites[, 2]*1.3,
             labels=unlist(dimnames(ordiplot$sites)[1]),
             col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],cex=0.7)
        
        image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.3,0.7), axis.args = list(cex.axis = 0.7),
                   zlim=range(-(rotsee.env.table$Depth_m)) )
        
        # Confidence ellipses for the Locations (SE +-95%)
        
        (groupz<-sort(unique(rotsee.env.table.log$Oxidation_Zone)))
        with(rotsee.env.table.log,ordiellipse(ordiplot, Oxidation_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[1], show.groups=groupz[1],draw="polygon", alpha=30))
        with(rotsee.env.table.log,ordiellipse(ordiplot, Oxidation_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[1], show.groups=groupz[2],draw="polygon", alpha=30))
        with(rotsee.env.table.log,ordiellipse(ordiplot, Oxidation_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[1], show.groups=groupz[3],draw="polygon", alpha=30))
        with(rotsee.env.table.log,ordiellipse(ordiplot, Oxidation_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[1], show.groups=groupz[4],draw="polygon", alpha=30))
    # plot the significant environmental variables from the selected canonical model
    
    (curr.plot.multiplier <- vegan:::ordiArrowMul(scores(Cluster.rda.best.forwsel$CCA$biplot[,c(xx,yy)]),
                                                  at = c(0,0), fill = 1,display="sites", choices = c(1,2)))# multiplier to fit the current plot ratio
    arrows(0, 0, Cluster.rda.best.forwsel$CCA$biplot[, xx]*curr.plot.multiplier*0.7, Cluster.rda.best.forwsel$CCA$biplot[, yy]*curr.plot.multiplier*0.7, lty=1,lwd=1.5,
           length = 0.15, angle = 30,  col="darkgrey")
    
    
    text(Cluster.rda.best.forwsel$CCA$biplot[, xx]*curr.plot.multiplier*0.9, Cluster.rda.best.forwsel$CCA$biplot[, yy]*curr.plot.multiplier*0.9,
         labels=unlist(dimnames(Cluster.rda.best.forwsel$CCA$biplot)[1]), col="darkgrey", cex=0.8)
    
    
    
    # Optionally: plot the rest of the environmental variables if you wish to....
    (envfit.rest  <- envfit(Cluster.rda.best.forwsel, rest.model.env.selection[,1:ncol(rest.model.env.selection)], choices=c(xx,yy),na.rm=TRUE))
    pval.limit <- which(c(envfit.rest$vectors$pvals)<=1) # set max p-value 
    
    (envfit.vectors       <- scores(envfit.rest,"vectors")[c(pval.limit),1:2])  # these are envfit vectors scaled to each one another (i.e. multipied by square root of column r2)
    (env.loading.factor   <- sqrt((envfit.rest$vectors$r)[c(pval.limit)]))      # loading factor for the envfit scores to scale them according to their correlations
    (curr.plot.multiplier <- vegan:::ordiArrowMul(scores(envfit.rest,"vectors")[c(pval.limit),1:2],
                                                  at = c(0,0), fill = 1,display="sites", choices = c(1,2)))# multiplier to fit the current plot ratio
    env.arrows.x  <- as.data.frame(scores(envfit.rest,"vectors"))[c(pval.limit),1]
    env.arrows.y  <- as.data.frame(scores(envfit.rest,"vectors"))[c(pval.limit),2]
    env.names     <-rownames(as.data.frame(scores(envfit.rest,"vectors")))[c(pval.limit)]
    
    arrows(0, 0, env.arrows.x*curr.plot.multiplier*0.55, env.arrows.y*curr.plot.multiplier*0.55, lty=1, col="lightgrey",
           length = 0.15, angle = 30)
    text(env.arrows.x*curr.plot.multiplier*0.72, env.arrows.y*curr.plot.multiplier*0.72,
         labels=env.names, col="lightgrey")    
    
    
    
    # plot MOB centroids
    
   # points(ordiplot$species[,1], ordiplot$species[,2],
   #       cex=2, pch=24, col="orange", bg=alpha("red", 0.1))
    
    text(ordiplot$species[,1]+0, ordiplot$species[,2]+0.1,
         attributes(ordiplot$species)[[2]][[1]][],
         col=pmmo_color[pMMO_vector_mobs], cex=0.7)
    
    text(ordiplot$species[,1]+0, ordiplot$species[,2]-0.1,
         colnames(t(rotsee.MOB.otu.table)),
         col=pmmo_color[pMMO_vector_mobs], cex=0.7)
    
    # plot legend and axis with explained variance
    
    legend("topleft",pch=c(21, 25),col=alpha(c("blue","red"),0.5),pt.bg=alpha(c("blue","red"),0.5),
           legend=c("Samples","OTUs"),bty="n")
    title(xlab=as.expression(paste("CCA",xx," (",round(c(100/(Cluster.rda.best.forwsel$tot.chi)*(Cluster.rda.best.forwsel$CCA$eig[xx])),2),"%)")),
          ylab=as.expression(paste("CCA",yy," (",round(c(100/(Cluster.rda.best.forwsel$tot.chi)*(Cluster.rda.best.forwsel$CCA$eig[yy])),2),"%)")))
    
    
    #dev.off()
    
    
    ## you might want to check the linearity of important variables in the cca
    
    
    gam.PhysicoChemRaw<-ordisurf(as.formula(paste("Cluster.rda.best.forwsel ~ ",paste(colnames(PhysicoChemicalRaw[parame<-"Cu_DGT_nM"])))), PhysicoChemicalRaw,
                                 choices=c(1,2), knots=10,                                         
                                 family="gaussian", isotropic= TRUE,
                                 col="grey",scaling=3,
                                 add = TRUE, display = "sites",
                                 nlevels = 15, labcex = 0.8,
                                 bubble = TRUE, cex = 1, select = TRUE, method = "REML",lwd=0.2,
                                 gamma = 1, plot = TRUE)
    (summary.gam.PhysicoChemRaw <- summary(gam.PhysicoChemRaw) )
    
    legend(1,-2, bty="n", legend=as.expression(paste("Explained deviance of\n",parame,":\n",   # change according to what you choose
                                                     round(c(summary.gam.PhysicoChemRaw$dev.expl)*100,2),
                                                     "%   ","P=",
                                                     round(c(summary.gam.PhysicoChemRaw$s.pv[[1]]),3))),cex=0.8)
    
    
    
    
    
#############################################################################################################################################################
##############################################################################################################################################################
## distance-based RDA with MOBs
###########################################################################################################################################################
#########################################################################################################################################################
    # data standardization and transformations
    
    ClusterHel           <-decostand(rotsee.vegan.mobs.otu.table, "hellinger")     # square root of equally standardized row sums
    ClusterNorm          <-decostand(rotsee.vegan.mobs.otu.table, method="total")  # standardized to equal row sum
    
    (PhysicoChemicalLog <- rotsee.mobs.env.table.log[,envvector[-c(1:5,25:33)]])  # log transformed env data
    (PhysicoChemicalLogsub<-rotsee.mobs.env.table.log[,envvector[-c(1:5,25:33)]])
    (PhysicoChemicalLogsub=PhysicoChemicalLogsub[colSums(!is.na(PhysicoChemicalLogsub)) >= nrow(PhysicoChemicalLogsub)] )# remove columns with NA's
    
    (col.env.vif<- vif_func(in_frame=PhysicoChemicalLogsub,thresh=20,trace=T))
    PhysicoChemicalLogsub = PhysicoChemicalLogsub[,c(col.env.vif)]
    
    
    (PhysicoChemicalRaw <- rotsee.mobs.env.table[,envvector[-c(1:5,25:33)]]) # original env data
    (PhysicoChemicalRawsub<-PhysicoChemicalRaw[colSums(!is.na(PhysicoChemicalRaw)) >= nrow(PhysicoChemicalRaw)])
    
    (ClusterHelsub <- ClusterHel[rownames(PhysicoChemicalLogsub),])
    
    ###################################################################################
    ## model assessments (forward, backward and marginal selection)
    ###################################################################################
    ## the selection is run on a sample subset lacking samples having missing values in the environmental matrix
    (Cluster.capscale<-capscale(ClusterHel~ ., PhysicoChemicalLogsub , distance="bray", metaMDSdist=FALSE, comm=ClusterHel, sqrt.dist = TRUE, na.action=na.omit))   # we use the log transformed env variables
    
    
    (Cluster.capscale.forward <- ordistep (capscale(ClusterHelsub~1, PhysicoChemicalLogsub, distance="bray", metaMDSdist=FALSE, sqrt.dist = TRUE, comm=ClusterHel, na.action=na.omit), # ordiR2step must inherit from rda...cca not working
                                           scope= formula(Cluster.capscale),
                                           direction="forward",pstep=100,nperm=99))
    
    (Cluster.capscale.backward <- ordistep (capscale(ClusterHelsub~., PhysicoChemicalLogsub, distance="bray", metaMDSdist=FALSE, sqrt.dist = TRUE, comm=ClusterHel, na.action=na.omit),
                                            direction="backward",pstep=100,nperm=99))
    
    (Cluster.capscale.both <- ordistep (capscale(ClusterHelsub~1, PhysicoChemicalLogsub, distance="bray", metaMDSdist=FALSE, sqrt.dist = TRUE, comm=ClusterHel, na.action=na.omit),
                                        scope= formula(Cluster.capscale),
                                        direction="both",pstep=100,nperm=99))
    
    ####################################################################################
    ## the best terms for the final model are assessed
    ####################################################################################
    (summarized.selection<- c(rownames(Cluster.capscale.forward$CCA$biplot),
                              rownames(Cluster.capscale.backward$CCA$biplot),
                              rownames(Cluster.capscale.both$CCA$biplot)))
    
    (sortsumsel<- (sort(table(summarized.selection))))
    (constrains.freq.table<- as.data.frame (sortsumsel))
    
    ## to what ever reason you have to use one of the two FORMULA functions... seems to be a R bug.... just uncomment and run the other if you get an error
    ## the same is true for the removement of single terms further below
    
    (FORMULA <-as.formula(paste("ClusterHel ~ ", paste  (rownames((constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=3))])), # here take again complete data set
                                                         collapse="+"),sep = "")))
    
    (FORMULA <-as.formula(paste("ClusterHel ~ ", paste ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=2))]), # here take again complete data set
                                                         collapse="+"),sep = "")))
    
    ###########################################################################################################
    # Run selected CCA model and tests and apriori fit physico-chemical variables not included in
    # the model
    ###########################################################################################################
    (Cluster.capscale.forwsel <- capscale(FORMULA, PhysicoChemicalLog, distance="bray", metaMDSdist=FALSE, comm=ClusterHel, sqrt.dist = TRUE, na.action=na.exclude))
    summary(Cluster.capscale.forwsel)
    
    # test the variation inflation factors and significance of models, axes and constraints
    
    sort(vif.cca(Cluster.capscale.forwsel) )  # remove variables from the data set that are highly colinear 
    
    (FORMULA    <-as.formula(paste("ClusterHel ~ ", mgsub(c("\\+Depth_m","\\+Depth_DGT_m"), c("", ""),
                                                          paste  (rownames(constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=2))]), # here take again complete data set
                                                                         collapse="+")),sep = "")))
    
    (FORMULA    <-as.formula(paste("ClusterHel ~ ", mgsub(c("\\+Depth_m","\\+Depth_DGT_m"), c("", ""),
                                                          paste ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=2))]), # here take again complete data set
                                                                 collapse="+")),sep = "")))
    
 
    (Cluster.capscale.best.forwsel <- capscale(FORMULA, PhysicoChemicalLog, distance="bray",comm=rotsee.vegan.mobs.table, metaMDSdist=FALSE,  sqrt.dist = TRUE, na.action=na.exclude))   # this is the final model
    summary(Cluster.capscale.best.forwsel)
    vif.cca(Cluster.capscale.best.forwsel)
    
    anova(Cluster.capscale.best.forwsel) # check if total model is singificant
    
    anova(Cluster.capscale.best.forwsel, by="axis", permutations = how(nperm=999),  # check for significant axes
          model = c("direct"),
          strata = NULL,
          cutoff = 1)
    
    
    #(FORMULAsub    <-as.formula(paste("ClusterHelsub ~ ", mgsub(c("\\+Fe_tot_nM"),
    #                                                            c(""), paste  (rownames(constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=3))]), # here take again complete data set
    #                                                                           collapse="+")),sep = "")))
    
    #(FORMULAsub    <-as.formula(paste("ClusterHel ~ ", mgsub(c("\\+Fe_tot_nM","\\+O2_uM","\\+T_C"),
    #                                                     c("","",""), paste ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=2))]), # here take again complete data set
    #                                                             collapse="+")),sep = "")))
    #FORMULAsub=FORMULA
    
    (Cluster.capscale.best.forwsel.sub <- capscale(FORMULAsub, PhysicoChemicalLogsub, distance="bray", metaMDSdist=FALSE, comm=ClusterHel, sqrt.dist = TRUE, na.action=na.omit))  # check for significance of environmental terms
    anova(Cluster.capscale.best.forwsel.sub, by="terms", permutations = how(nperm=999),
          model = c("direct"), 
          parallel = getOption("mc.cores"), strata = NULL,
          scope = NULL)
    
    
    # fit the non-used environmental variables in the RDA model
    
    (best.model.env.selection <-c(rownames(Cluster.capscale.best.forwsel$CCA$biplot)))
    (rest.model.env.selection <-PhysicoChemicalRaw[! c(colnames(PhysicoChemicalRaw)) %in% best.model.env.selection])
    
    (envfit.rest  <- envfit(Cluster.capscale.best.forwsel, rest.model.env.selection[,1:length(rest.model.env.selection)], na.rm=TRUE))
    
    
    # and some additional diagnostics
    
    intersetcor(Cluster.capscale.best.forwsel)
    
    dispvect<-c("sp","wa","cn")
    for (i in c(1:3)){
      print(scores(Cluster.capscale.best.forwsel, choices = c(1,2), display =dispvect[i] ,
                   scaling = "symmetric", correlation = TRUE))
    }
    
    screeplot(Cluster.capscale.best.forwsel, bstick = FALSE, type = c("lines"), #"barplots"
              npcs = min(10, if (is.null(Cluster.capscale.best.forwsel$CCA)) Cluster.capscale.best.forwsel$CA$rank else Cluster.capscale.best.forwsel$CCA$rank),
              ptype = "o", bst.col = "red", bst.lty = "solid",
              xlab = "Component", ylab = "Inertia")
    
    spenvcor(Cluster.capscale.best.forwsel)
    
    # goodness(Cluster.capscale.best.forwsel, display = c("sites"), 
    #          model = c("CCA", "CA"), statistic = c("explained", "distance"),
    #          summarize = TRUE, addprevious = FALSE)
    
    inertcomp(Cluster.capscale.best.forwsel, display = c("sites"),
              unity = FALSE, proportional = FALSE)
    
    # run the percentage contribution of constrains function (defined at at the end of the script)
    
    PerContCon(Cluster.capscale.best.forwsel)          # you can plot the axes which show most explained variance by the variable of interest
    
    # check R2 of apriori fitted constraints
    (envfit.apriori  <- envfit(Cluster.capscale.best.forwsel, PhysicoChemicalRawsub, na.rm=TRUE, choices=c(1,2)))  # choose the axes according to the most explained fitting axes (PerContCon output)
    # to check specific variable of interest (i.e. Cu_DGT axes 2,3 r2 is high
    # whereas axes 1,2 r2 is low...)
    
    ##########################################################################################################
    # Bi-plot of the CCA model          (all following plots will be directly saved as PDF's)
    ##########################################################################################################
    my.color1<-c("cyan", "cyan4", "cyan3")
    MOBsOTUsnumbers
    
    #pdf("dbRDA of MOB Clusters with Env as constraint.pdf",useDingbats=FALSE)
    
    # Basic Plot
    
    ordiplot<-plot(Cluster.capscale.best.forwsel,type="none",choices = c(xx<-1, yy<-2),scaling=2,main="", xlab="",ylab="")     # here you can change the axes to be ploted       
    
    points(ordiplot,"species", pch=25,cex=normalized.species.size.vec * 2,
           col=alpha(pmmo_color[pMMO_vector_mobs], 0.8), bg=alpha(pmmo_color[pMMO_vector_mobs],0.3))
    
    points(ordiplot,"sites", cex=3, alpha= 2, pch=as.numeric(as.factor(rotsee.env.table$Oxidation_Zone))+20, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
           bg=alpha(col.depth[as.numeric(rotsee.mobs.env.table$Depth_m)/max(rotsee.mobs.env.table$Depth_m)*1000], 0.5))
    
    text(ordiplot$sites[, 1]*1.3, ordiplot$sites[, 2]*1.3,
         labels=unlist(dimnames(ordiplot$sites)[1]),
         col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],cex=0.7)
    
    image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.3,0.7), axis.args = list(cex.axis = 0.7),
               zlim=range(-(rotsee.env.table$Depth_m)) )
    
    # Confidence ellipses for the Locations (SE +-95%)
    (groupz<-sort(unique(rotsee.env.table.log$Oxidation_Zone)))
    with(rotsee.env.table.log,ordiellipse(ordiplot, Oxidation_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[1], show.groups=groupz[1],draw="polygon", alpha=30))
    with(rotsee.env.table.log,ordiellipse(ordiplot, Oxidation_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[1], show.groups=groupz[2],draw="polygon", alpha=30))
    with(rotsee.env.table.log,ordiellipse(ordiplot, Oxidation_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[1], show.groups=groupz[3],draw="polygon", alpha=30))
    with(rotsee.env.table.log,ordiellipse(ordiplot, Oxidation_Zone, kind="se", conf=0.95, label=TRUE, cex=0.8, col=my.color1[1], show.groups=groupz[4],draw="polygon", alpha=30))
    
    # plot the significant environmental variables from the selected canonical model
    
    (curr.plot.multiplier <- vegan:::ordiArrowMul(scores(Cluster.capscale.best.forwsel$CCA$biplot[,c(xx,yy)]),
                                                  at = c(0,0), fill = 1,display="sites", choices = c(xx,yy)))# multiplier to fit the current plot ratio
    arrows(0, 0, Cluster.capscale.best.forwsel$CCA$biplot[, xx]*curr.plot.multiplier*0.7, Cluster.capscale.best.forwsel$CCA$biplot[, yy]*curr.plot.multiplier*0.7, lty=1,lwd=1.5,
           length = 0.15, angle = 30,  col="darkgrey")
    
    
    text(Cluster.capscale.best.forwsel$CCA$biplot[, xx]*curr.plot.multiplier*0.9, Cluster.capscale.best.forwsel$CCA$biplot[, yy]*curr.plot.multiplier*0.9,
         labels=unlist(dimnames(Cluster.capscale.best.forwsel$CCA$biplot)[1]), col="darkgrey", cex=0.8)
    
    
    
    # Optionally: plot the rest of the environmental variables if you wish to....
    
    (envfit.rest  <- envfit(Cluster.capscale.best.forwsel, rest.model.env.selection[,11:length(rest.model.env.selection)], na.rm=TRUE,  choices=c(xx,yy)))
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
    
    # plot MOB centroids
    
    # points(ordiplot$species[,1], ordiplot$species[,2],
    #       cex=2, pch=24, col="orange", bg=alpha("red", 0.1))
    
    text(ordiplot$species[,1]+0, ordiplot$species[,2]+0.1,
         attributes(ordiplot$species)[[2]][[1]][],
         col=pmmo_color[pMMO_vector_mobs], cex=0.7)
    
    text(ordiplot$species[,1]+0, ordiplot$species[,2]-0.1,
         colnames(t(rotsee.mobs.otu.table)),
         col=pmmo_color[pMMO_vector_mobs], cex=0.7)
    
    # plot legend and axis with explained variance
    
    legend("bottomleft",pch=c(21, 25),col=alpha(c("blue","red"),0.5),pt.bg=alpha(c("blue","red"),0.4),
           legend=c("Samples","OTUs"),bty="n")
    title(xlab=as.expression(paste("CAP",xx," (",round(c(100/(Cluster.capscale.best.forwsel$tot.chi)*(Cluster.capscale.best.forwsel$CCA$eig[xx])),2),"%)")),
          ylab=as.expression(paste("CAP",yy," (",round(c(100/(Cluster.capscale.best.forwsel$tot.chi)*(Cluster.capscale.best.forwsel$CCA$eig[yy])),2),"%)")))
    
    
    #dev.off()
    
    
    ## you might want to check the linearity of important variables in the cca
    
    
    gam.PhysicoChemRaw<-ordisurf(as.formula(paste("Cluster.capscale.best.forwsel ~ ",paste(colnames(PhysicoChemicalRaw[parame<-"Cu_DGT_nM"])))), PhysicoChemicalRaw,
                                 choices=c(xx,yy), knots=10,                                         
                                 family="gaussian", isotropic= TRUE,
                                 col="grey",scaling=3,
                                 add = TRUE, display = "sites",
                                 nlevels = 15, labcex = 0.8,
                                 bubble = TRUE, cex = 1, select = TRUE, method = "REML",lwd=0.2,
                                 gamma = 1, plot = TRUE)
    (summary.gam.PhysicoChemRaw <- summary(gam.PhysicoChemRaw) )
    
    legend(0.7,1.5, bty="n", legend=as.expression(paste("Explained deviance of\n",parame,":\n",   # change according to what you choose
                                                        round(c(summary.gam.PhysicoChemRaw$dev.expl)*100,2),
                                                        "%   ","P=",
                                                        round(c(summary.gam.PhysicoChemRaw$s.pv[[1]]),3))),cex=0.8)
    
################################################################################################################################################################    
################################################################################################################################################################
## Bar Plots of subsamples
############################################################################################################################################################### 
###############################################################################################################################################################
    plot_ordered_bar<-function (physeq, x = "Sample", 
                                y = "Abundance", 
                                fill = NULL, 
                                leg_size = 0.5,
                                title = NULL) {
      require(ggplot2)
      require(phyloseq)
      require(plyr)
      require(grid)
      bb <- psmelt(physeq)
      
      samp_names <- aggregate(bb$Abundance, by=list(bb$Sample), FUN=sum)[,1]
      .e <- environment()
      bb[,fill]<- factor(bb[,fill], rev(sort(unique(bb[,fill])))) #fill to genus
      
      
      bb<- bb[order(bb[,fill]),] # genus to fill
      p = ggplot(bb, aes_string(x = x, y = y, 
                                fill = fill), 
                 environment = .e, ordered = FALSE)
      
      
      p = p +geom_bar(stat = "identity", 
                      position = "stack", 
                      color = "black") 
      
      p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
      
      p = p + guides(fill = guide_legend(override.aes = list(colour = NULL), reverse=TRUE)) + 
        theme(legend.key = element_rect(colour = "black")) 
      
      p = p + theme(legend.key.size = unit(leg_size, "cm"))
      
      
      if (!is.null(title)) {
        p <- p + ggtitle(title)
      }
      return(p)
    }
    
    
    
  
    ## barplot raw reads
    MOB_OTU = prune_taxa(MOBs, bac.rotsee.and.MOBs)
    bar <- plot_ordered_bar(MOB_OTU, x="Depth_m", fill="Family", title="MOBs")
    bar + facet_wrap(~OTU)  
    ## barplot standardized to all OTUs sequencing depth
    MOB_OTU = prune_taxa(MOBs, bac.rotsee.and.MOBs.seqdepth)
    bar <- plot_ordered_bar(MOB_OTU, x="Depth_m", fill="Family", title="MOBs")
    bar + facet_wrap(~OTU*Campaign)

    # dev.copy2pdf(file="MOBs by campaign and by Depth.pdf", width = 7, height = 5)
    
    ## write cvs with standardized OTUs 
    # MOB_OTU_table_stseq<- otu_table(MOB_OTU)
    # write.csv(MOB_OTU_table_stseq, file ="rotsee_MOB_otu_table_stseqdepth_AllOTUs.csv")   
    
#########################################################################################################################################
########################################################################################################################################
# Pearson Correlation of MOBs with physico-chemical variables
########################################################################################################################################
#######################################################################################################################################
    ClusterHelPearson   <-decostand(rotsee.vegan.mobs.otu.table, "hellinger")     # square root of equally standardized row sums
    (colnames(ClusterHelPearson) <- paste(colnames(ClusterHelPearson),"_",rownames(rotsee.mobs.otu.table), sep="")) 
    
    PhysicoChemicalLog <- cbind(rotsee.env.table[,c(37),drop=FALSE],log1p(rotsee.env.table[,c(18:28,31,33,34,36,38:40,46)]))  # log transformed env data
    PhysicoChemicalLogsub<-na.omit(PhysicoChemicalLog)
    PhysicoChemicalsub<-rotsee.env.table[,c(18:28,31,33,34,36,37,38:40,46)]
    
    #rotsee.mobs.env.table
    #rotsee.mobs.env.table.log
    names(rotsee.mobs.env.table.log)
    envvector<- c(14,19:47,68:70)
    
    PhysicoChemicalLogsub<-na.omit(rotsee.mobs.env.table.log[envvector])
    
    
    
    #PhysicoChemicalLogsub=rotsee.env.table.log[,-c(1:13, 15:18, 48:67)]
    
    rcorrMatrix<-as.matrix(cbind(ClusterHelPearson[rownames(PhysicoChemicalLogsub),]  ,  PhysicoChemicalLogsub))  
    #rcorrMatrix=rcorrMatrix[,-c()] # remove cols with NAs
    sink(file = "PearsonCorrPhyChemClus.txt", append = FALSE, type = c("output"),
         split = FALSE)
    (CorrTable<-rcorr.adjust(rcorrMatrix[,], type="pearson"))
    sink()
    
    
    
    #(ClusterHelPearson=ClusterHelPearson[,-c()])
    matrixPearson3<-CorrTable$R$r[c(1:ncol(ClusterHelPearson  )),c((ncol(ClusterHelPearson  )+1):nrow(CorrTable$R$r))]
    pvalCellNote<-(CorrTable$R$P[c(1:ncol(ClusterHelPearson  )),c((ncol(ClusterHelPearson  )+1):nrow(CorrTable$R$r))])
    pvalCellNoteStars<- cut(pvalCellNote,  breaks=c(-Inf, 0.001, 0.01, 0.05),  label=c("***", " ** ", "  *  ")) 
    pvalCellNote<-matrix((pvalCellNoteStars), nrow=nrow(matrixPearson3), ncol=ncol(matrixPearson3))
    
    genus.df.mobs.sub=genus.df.mobs#[-c(),]
    #pdf("heatmap physicochemical cluster pearson_16smobs.pdf", height=10, width=20)
    heatmap.2(t(matrixPearson3),
              Rowv=TRUE,
              Colv=TRUE,
              dendrogram= c("column"),
              #distfun = NULL,
              key=TRUE,
              keysize=1, 
              trace="none",
              density.info=c("none"),
              margins=c(10, 8),
              col=colorRampPalette(brewer.pal(9,"Blues"))(100),
              cellnote=t(pvalCellNote),
              notecex=1.0,
              notecol="black",
              na.color=par("bg"),
              scale="none", 
              symkey=TRUE,
              cexRow=1,
              cexCol=1,
              srtRow=35,
              srtCol=35,
              ColSideColors=c(col.genus.mobs[genus.df.mobs.sub$id]),
              key.title="Pearson",
              key.xlab=NA
              #col=redgreen(75),
    )
    
    legend("bottomleft",c(unique(rotsee.mobs.tax.table[,"Order"])),pch=c(21),col=alpha(unique(col.genus.mobs[genus.df.mobs$id]), 1),
           pt.bg=alpha(unique(col.genus.mobs[genus.df.mobs$id]),1), bty="n", cex=0.7)
    
    #dev.off()
    
    # heatmap with zoom function
    d3heatmap(t(matrixPearson3),
              Rowv=TRUE,
              Colv=TRUE,
              dendrogram= c("column"),
              #distfun = NULL,
              key=TRUE,
              keysize=1, 
              trace="none",
              density.info=c("none"),
              margins=c(10, 8),
              col=colorRampPalette(brewer.pal(9,"Blues"))(100),
              cellnote=t(pvalCellNote),
              notecex=1.0,
              notecol="black",
              na.color=par("bg"),
              scale="none", 
              symkey=TRUE,
              cexRow=1,
              cexCol=1,
              srtRow=35,
              srtCol=35,
              #ColSideColors=c(col.genus.mobs[genus.df.mobs$id]),
              key.title="Pearson",
              key.xlab=NA
              #col=redgreen(75),
    )
    
################################################################################################################################################################
################################################################################################################################################################
## correlation plots of MOBs to env split by campaigns
################################################################################################################################################################
################################################################################################################################################################
    ClusterHelPearsonsub = ClusterHelPearson[row.names(PhysicoChemicalLogsub),]
    rcorrMatrixCorr<-as.matrix(cbind( ClusterHelPearsonsub , PhysicoChemicalLogsub))
    #rcorrMatrixCorr<-as.matrix(cbind(ClusterHelPearson  , PhysicoChemicalsub))
    rcorrMatrixCorr<- cbind( rcorrMatrixCorr,na.omit(rotsee.mobs.env.table[,envvector[-c(1:5,25:33)]]))
   

 ## split by campaign 
   pdf(file="Mob correlations September 2014.pdf")  # change name according to campaign

   for (i in c(1:ncol(rotsee.vegan.mobs.otu.table))){
     
               CorMob<-  paste(colnames(rcorrMatrixCorr[i]))
                       for (p in c(1:ncol(PhysicoChemicalLogsub))){
                         
                            #FORMULAcorr = as.formula(paste(CorMob, "~", colnames(PhysicoChemicalLogsub)[p],"| Campaign", sep=""))
                            FORMULAcorr = as.formula(paste(CorMob, "~", colnames(PhysicoChemicalLogsub)[p], sep=""))
                            scatterplot(FORMULAcorr, data=rcorrMatrixCorr)
                            }
   }
   dev.off()

################################################################################################################################################################
################################################################################################################################################################
## Network Analysis (some examples, use the specific scripts for in depth network analysis)
###############################################################################################################################################################
################################################################################################################################################################
ig = make_network(bac.rotsee.mobs, type="samples", max.dist = 0.5
                  , distance="jaccard")
plot_network(ig, bac.rotsee.prune.filt3, color = "Depth_m", shape = "Oxidation_Zone", line_weight = 0.2, 
                           label = "value")


gp100= prune_taxa(names(sort(taxa_sums(bac.rotsee.prune.filt3), TRUE))[1:250], bac.rotsee.prune.filt3)
jg = make_network(gp100, "taxa", "jaccard", 0.7)
plot_network(jg, gp100, "taxa", line_weight = 0.4, label = "value")

################################################################################################################################################################
################################################################################################################################################################
## multiple linear regression with environmental variables vs. mobs.... run this with the complete data set as too many df for env variables compared to samples
################################################################################################################################################################
################################################################################################################################################################
#Some of the MOBs occure just once...remove these as it does not make sense to model these OTUs for single season (rahter state that the on spot environmental
# conditions are important for these specific oTUs)
rotsee.env.table.log[,envvector[-c(1:5,25:33)]]
(rotsee.mlr.env <- na.omit(rotsee.env.table.log[,envvector[-c(1:5,25:33)]]))   # here only take ~half the env variables as samples
(mob.mlr1 <- decostand(as.data.frame(rotsee.vegan.mobs.otu.table[row.names(rotsee.mlr.env),which(colSums(rotsee.vegan.mobs.otu.table)>10)]), method="hellinger") ) # here you can choose the minimal number of reads
(rotsee.mobs.table.mlr1 <- rotsee.mobs.table.t[row.names(rotsee.mlr.env), which(colSums(rotsee.vegan.mobs.otu.table)>10)])

selectvector=vector()
for (i in c(1:ncol(mob.mlr1))){
  if (sum(mob.mlr1[,i]>0)>5){selectvector<-c(selectvector, i)}   # here you can define minimal integration points (aka site occurence) have to be in the model
}
(mob.mlr<- subset(mob.mlr1, select= selectvector))
(mob.mlr.OTU<- subset(rotsee.mobs.table.mlr1, select= selectvector))

## This is the AIC approach...this is currently discussed and it might be better to use the "leap" approach below
##----------------------------------------------------
mylistMLR=list()
sink("Relaimpo_Env_on_MOBs_All_AIC.txt", append=TRUE, split=TRUE)

for (i in c(1:ncol(mob.mlr))){
  
  fit <- lm(mob.mlr[,i] ~ . ,data= rotsee.mlr.env)#
  step <- stepAIC(fit, direction="both")
  (step$anova) # display results
  summary(step)                   # TYPE 2 and Type 3 SS Anova  (because no interaction terms apparent)
  options(contrasts=c("contr.sum", "contr.poly"))
  significant_reg <- anova(step)
  coefficients_reg <- coefficients(step) # model coefficients
  confint(step, level=0.95) # CIs for model parameters
  fitted(step) # predicted values
  residuals(step) # residuals
  vcov(step) # covariance matrix for model parameters
  influence(step) # regression diagnostics
  (relimportance<-calc.relimp(step,type=c("lmg"),rela=TRUE))
  name<-paste("relaimpo:  ",colnames(mob.mlr.OTU)[i], " / ", colnames(mob.mlr)[i],sep="")
  tmp<- list(significance=significant_reg, coefficients=coefficients_reg, relaimpo=relimportance)
  mylistMLR[[name]] <- tmp
  
  
}
mylistMLR
sink()

## for model selection you can also use Leaps...
##-----------------------------------------------------
mob.mlr.leaps<-cbind(mob.mlr, rotsee.mlr.env)
mylistMLR=list()

sink("Relaimpo_Env_on_MOBs_All_Leaps.txt", append=TRUE, split=TRUE)

for (i in c(1:ncol(mob.mlr))){
  (FORMULAleaps = as.formula(paste(colnames(mob.mlr.leaps)[i]," ~ ", paste (colnames(mob.mlr.leaps[,-c(1:ncol(mob.mlr))]), # here take again complete data set
                                                                            collapse="+"),sep = "")))
  regsubsets.out <-  regsubsets( FORMULAleaps,
                                 data = mob.mlr.leaps,
                                 nbest = 10,       # 1 best model for each number of predictors
                                 
                                 nvmax = floor(nrow(mob.mlr1)/9),   
                                 
                                 
                                 
                                 force.in = NULL, force.out = NULL,
                                 method = "exhaustive")
  regsubsets.out
  summary.out <- summary(regsubsets.out)
  as.data.frame(summary.out$outmat)
  #plot(regsubsets.out, scale = "adjr2", main = "Adjusted R^2")
  
  #res.legend <-  subsets(regsubsets.out, statistic="adjr2", legend = FALSE, min.size = 1, main = "Adjusted R^2")
  #res.legend <-  subsets(regsubsets.out, statistic="cp", legend = FALSE, min.size = 1, main = "Mallow Cp")
  #abline(a = 1, b = 1, lty = 2)
  
  #res.legend
  highR2<-which.max(summary.out$adjr2)
  leapsparam_all = summary.out$which[highR2,]
  leapsparam_selected = names(which(leapsparam_all==TRUE))[-1]
  
  (FORMULAleaps = as.formula(paste(colnames(mob.mlr.leaps)[i]," ~ ", paste ((leapsparam_selected), 
                                                                            collapse="+"),sep = "")))
  
  
  fit <- lm(FORMULAleaps ,data= mob.mlr.leaps)#
  summary(fit)                   # TYPE 2 and Type 3 SS Anova
  options(contrasts=c("contr.sum", "contr.poly"))
  significant_reg <- anova(fit)
  coefficients_reg <- coefficients(fit) # model coefficients
  confint(fit, level=0.95) # CIs for model parameters
  fitted(fit) # predicted values
  residuals(fit) # residuals
  vcov(fit) # covariance matrix for model parameters
  influence(fit) # regression diagnostics
  (relimportance<-calc.relimp(fit,type=c("lmg"),rela=TRUE))
  name<-paste("relaimpo:  ",colnames(mob.mlr.OTU)[i], " / ", colnames(mob.mlr)[i],sep="")
  tmp<- list(significance=significant_reg, coefficients=coefficients_reg, relaimpo=relimportance)
  mylistMLR[[name]] <- tmp
}
mylistMLR

sink()

########################################################################################################################
###########################################################################################################################
## Cooccurence Analysis of mobs and other bacteria
########################################################################################################################
###########################################################################################################################
library(corrplot)
library(cooccur)
# The program calculates the observed and expected frequencies of co-occurrence between each pair of species.
# Expected frequency is based on the distribution of each species being random and independent of the other species.
rotsee.mobs.table=rotsee.mobs.otu.table

## mobs and cyanos cooccur
(bacterianumber <- which(rotsee.tax.table@.Data[,"Phylum"]==c("Cyanobacteria")))
cyanonames=names(bacterianumber)
cyanos = prune_taxa(cyanonames, bac.rotsee)
rotsee.cyanobac.table=otu_table(cyanos)
cyano.tax.table=tax_table(cyanos)
rotsee.vegan.cyano.table=make_vegan_otu_table (rotsee.cyanobac.table, cyano.tax.table)
(colnames(rotsee.vegan.cyano.table) <- paste(mgsub(c("\\.", "/","-"),c("","_","_"),(colnames(rotsee.vegan.cyano.table))),"_",rownames(rotsee.cyanobac.table[,]), sep=""))
rotsee.mobs.table.cooccur=rotsee.mobs.table
(rownames(rotsee.mobs.table.cooccur) <- paste(mgsub(c("\\.", "/","-"),c("","_","_"),(colnames(rotsee.vegan.mobs.otu.table))),"_",rownames(rotsee.mobs.table[,]), sep=""))
cooccur_mob_cyano.table= rbind(t(rotsee.vegan.cyano.table), rotsee.mobs.table.cooccur)
OTU_table_cooccur <- decostand(cooccur_mob_cyano.table@.Data, "pa")[,]

## mobs cooccur only
OTU_table_cooccur <- decostand(rotsee.mobs.table@.Data, "pa")[,]
#OTU_table_cooccur   <-rotsee.vegan.otu.table   
(rownames(OTU_table_cooccur) <- paste(mgsub(c("\\.", "/","-"),c("","_","_"),(colnames(rotsee.vegan.mobs.otu.table))),"_",rownames(rotsee.mobs.table[,]), sep=""))# replace characters that  are not parsed propperly

## cooccur analysis 
cooccur.mobs <- cooccur(OTU_table_cooccur[,],
                        type = "spp_site",
                        thresh = FALSE,
                        spp_names = TRUE
                        )
obs.v.exp(cooccur.mobs)
class(cooccur.mobs)
summary(cooccur.mobs)
print(cooccur.mobs)
effect.sizes(cooccur.mobs, standardized = TRUE, matrix = TRUE)
# EXTRACT INFORMATION FOR A FOCAL SPECIES

sppname = paste(rownames(OTU_table_cooccur )[which(rownames(rotsee.mobs.table)=="OTU_129")])  # put in the OTU number
(pair(mod = cooccur.mobs, spp = sppname))

pair.attributes(cooccur.mobs)

# generic cooccur plot
plot(cooccur.mobs)

###########################################################################################################################
###########################################################################################################################
## correlation plot with pearson correlation of mob with i.e. cyanobacteria
###########################################################################################################################
###########################################################################################################################
library(corrplot)

## spearman correlation 
OTU_table_cooccur <- decostand(rotsee.mobs.table@.Data, "pa")
#OTU_table_cooccur   <-rotsee.vegan.otu.table   
(rownames(OTU_table_cooccur) <- paste(mgsub(c("\\.", "/","-"),c("","_","_"),(colnames(rotsee.vegan.mobs.otu.table))),"_",rownames(rotsee.mobs.table[,]), sep=""))# replace characters that  are not parsed propperly

# cyanos and mobs
OTU_table_cooccur <- decostand(cooccur_mob_cyano.table@.Data, "pa")[,]

M<-cor(t(OTU_table_cooccur), method = c("spearman"))

## compute p-values of the correlations
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p.mat <- cor.mtest(t(OTU_table_cooccur))

# diverse plots
corrplot(M, method="circle", type="upper", col=brewer.pal(n=11, name="PuOr"))
corrplot(M, method="color", type="upper", order="hclust")
corrplot(M, method="circle", type="upper", order="hclust", col=brewer.pal(n=11, name="PuOr"), p.mat=p.mat, , sig.level = 0.05)
corrplot(M, method="circle", type="upper", order="hclust", col=brewer.pal(n=11, name="PuOr"), p.mat=p.mat, , sig.level = 0.05, insig="blank")
corrplot(M, method="circle", type="upper",                 col=brewer.pal(n=11, name="PuOr"), p.mat=p.mat, , sig.level = 0.05,
         insig="blank", tl.col=c(rep("green",19),rep("pink",21)))

# or network 
library("qgraph")

Graph_pcor <- qgraph(M,graph = "pcor", layout = "spring",layout = "spring",  minimum=0.2,
                      groups = c(rep("Cyanobacteria",ntaxa(cyanos)),rep("MOB",nrow(rotsee.mobs.table))), 
                      sampleSize = nrow(t(OTU_table_cooccur)),  pastel = TRUE, posCol = "#003399",
                      negCol = "#FF9933", borders = FALSE, vTrans = 200,
                      details = TRUE                      )



###########################################################################################################################
###########################################################################################################################
## compare cooccurence of mobs 16s with the cooccurence of the pmoA sequence data
## import the data and process as in the complete workflow for pmoA
###########################################################################################################################
###########################################################################################################################
otufile     <- "data/P349_ZOTU_pmoA_sintax_clean.tab"  
mapfile     <- "data/Metafile_impuded.txt"  
treefile    <- "data/P349_ZOTU.tre"
refseqfile  <- "data/P349_ZOTU.fa"

pmoA.rotsee <- import_qiime(otufilename = otufile, mapfilename = mapfile ,
                           treefilename = treefile, refseqfilename = refseqfile)
pmoA.rotsee.prune <- pmoA.rotsee
pmoA.rotsee.prune <- prune_taxa(taxa_sums(pmoA.rotsee.prune) > 0, pmoA.rotsee.prune)
ntaxa(pmoA.rotsee.prune) 
filt3 <- genefilter_sample(pmoA.rotsee.prune, filterfun_sample(function(x) x > 50), A=3)
pmoA.rotsee.prune.filt3 <- prune_taxa(filt3, pmoA.rotsee.prune)
ntaxa(pmoA.rotsee.prune.filt3) 

pmoA.rotsee <- pmoA.rotsee.prune.filt3

# pmoA.rotsee <- subset_samples(pmoA.rotsee.prune.filt3, Campaign=="A_June_2013")
# 
# pmoA.rotsee <- subset_samples(pmoA.rotsee.prune.filt3, Campaign=="B_August_2013")
# 
# pmoA.rotsee <- subset_samples(pmoA.rotsee.prune.filt3, Campaign=="C_September_2014")
# 
# pmoA.rotsee <- subset_samples(pmoA.rotsee.prune.filt3, Campaign=="D_September_2015")
# 
# pmoA.rotsee <- subset_samples(pmoA.rotsee.prune.filt3, Campaign!="D_June_2015") # without June

#pmoA.rotsee.taxmerge <- tax_glom(pmoA.rotsee.prune.filt3, taxrank="Genus")
#ntaxa(pmoA.rotsee.taxmerge)

rotsee.pmoA.otu.table=otu_table(pmoA.rotsee)
(rotsee.pmoA.tax.table <- tax_table(pmoA.rotsee))
(rotsee.pmoA.env.table <- get_variable(pmoA.rotsee))

(rotsee.vegan.pmoA.table<- make_vegan_otu_table (rotsee.pmoA.otu.table, rotsee.pmoA.tax.table))

(pmoA_table_cooccur <- decostand(rotsee.pmoA.otu.table@.Data, "pa")[,])
(rownames(pmoA_table_cooccur) <- paste(mgsub(c("\\.", "/","-"),c("","_","_"),(colnames(rotsee.vegan.pmoA.table))),"_",rownames(rotsee.pmoA.otu.table[,]), sep=""))# replace characters that  are not parsed propperly

minusvecsites.16s=c("A01","A02","A03","A05", "D12")
bac.rotsee.16s <- rotsee.mobs.table[,!(colnames(rotsee.mobs.table)%in%minusvecsites.16s)]
rotsee.mobs.table=bac.rotsee.16s

(OTU_table_cooccur <- decostand(rotsee.mobs.table@.Data, "pa"))
#OTU_table_cooccur   <-rotsee.vegan.otu.table   
(rownames(OTU_table_cooccur) <- paste(mgsub(c("\\.", "/","-"),c("","_","_"),(colnames(rotsee.vegan.mobs.otu.table))),"_",rownames(rotsee.mobs.table[,]), sep=""))# replace characters that  are not parsed propperly

(  pmoA_table_cooccur= pmoA_table_cooccur[,match(colnames(pmoA_table_cooccur), colnames(OTU_table_cooccur))]  )

otu_pmoa_table_cooccur = rbind(OTU_table_cooccur, pmoA_table_cooccur)



cooccur.mobs <- cooccur(otu_pmoa_table_cooccur[,],
                        type = "spp_site",
                        thresh = TRUE,
                        spp_names = TRUE)
class(cooccur.mobs)
summary(cooccur.mobs)
prob.table(cooccur.mobs)

# EXTRACT INFORMATION FOR A FOCAL SPECIES

sppname = paste(rownames(OTU_table_cooccur )[which(rownames(rotsee.mobs.table)=="OTU_129")])  # put in the OTU number
(pair(mod = cooccur.mobs, spp = sppname))

# generic cooccur plot
plot(cooccur.mobs)


effect.sizes(cooccur.mobs, standardized = TRUE, matrix = FALSE)

#############################################################################
## correlation plot with pearson correlation of mobs and pmoa
library(corrplot)

##-----------------------------------------------------
### prune data to redcue otu numbers for pmoa
bac.rotsee.pmoA <- tax_glom(pmoA.rotsee.prune.filt3, taxrank="Genus")
(rotsee.pmoA.otu.table <- otu_table(bac.rotsee.pmoA))
(rotsee.tax.table.pmoA <- tax_table(bac.rotsee.pmoA))
rotsee.pmoA.otu.table.vegan= make_vegan_otu_table (rotsee.pmoA.otu.table, rotsee.tax.table.pmoA)


(pmoA_table_cooccur <- decostand(rotsee.pmoA.otu.table@.Data, "hellinger")[,])
(rownames(pmoA_table_cooccur) <- paste(mgsub(c("\\.", "/","-"),c("","_","_"),(colnames(rotsee.pmoA.otu.table.vegan))),"_",rownames(rotsee.pmoA.otu.table[,]), sep=""))# replace characters that  are not parsed propperly

##-----------------------------------------------------
## or use all pmoA otus
bac.rotsee.pmoA=pmoA.rotsee.prune.filt3
(rotsee.pmoA.otu.table <- otu_table(bac.rotsee.pmoA))
(rotsee.tax.table.pmoA <- tax_table(bac.rotsee.pmoA))
rotsee.pmoA.otu.table.vegan= make_vegan_otu_table (rotsee.pmoA.otu.table, rotsee.tax.table.pmoA)

(pmoA_table_cooccur <- decostand(rotsee.pmoA.otu.table@.Data, "hellinger")[,])
(rownames(pmoA_table_cooccur) <- paste(mgsub(c("\\.", "/","-"),c("","_","_"),(colnames(rotsee.vegan.pmoA.table))),"_",rownames(rotsee.pmoA.otu.table[,]), sep=""))# replace characters that  are not parsed propperly

##-----------------------------------------------------
## mobs 16s
minusvecsites.16s=c("A01","A02","A03","A05", "D12")     ## there is no pmoa data for these sites
bac.rotsee.16s <- rotsee.mobs.table[,!(colnames(rotsee.mobs.table)%in%minusvecsites.16s)]
rotsee.mobs.table=bac.rotsee.16s

OTU_table_cooccur <- rotsee.mobs.table@.Data
(OTU_table_cooccur <- decostand(rotsee.mobs.table@.Data, "hellinger"))
#OTU_table_cooccur   <-rotsee.vegan.otu.table   
(rownames(OTU_table_cooccur) <- paste(mgsub(c("\\.", "/","-"),c("","_","_"),(colnames(rotsee.vegan.mobs.otu.table))),"_",rownames(rotsee.mobs.table[,]), sep=""))# replace characters that  are not parsed propperly

## match and merge
(  pmoA_table_cooccur= pmoA_table_cooccur[,match(colnames(pmoA_table_cooccur), colnames(OTU_table_cooccur))]  )

otu_pmoa_table_cooccur = rbind(OTU_table_cooccur, pmoA_table_cooccur)

## run corr analysis
M<-cor(t(otu_pmoa_table_cooccur), method = c("spearman"))


## compute p-values of the correlations
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p.mat <- cor.mtest(t(otu_pmoa_table_cooccur))


# diverse plots
corrplot(M, method="circle", type="upper",tl.cex=0.3, col=brewer.pal(n=11, name="PuOr"))
corrplot(M, method="color", type="upper", order="hclust",tl.cex=0.8, col=brewer.pal(n=11, name="PuOr"),p.mat=p.mat, sig.level = 0.05, , insig="blank")
corrplot(M, method="circle", type="upper",tl.cex=0.3, order="hclust", col=brewer.pal(n=11, name="PuOr"), p.mat=p.mat,  sig.level = 0.01, , insig="blank")
corrplot(M, method="circle", type="upper", order="hclust", col=brewer.pal(n=11, name="PuOr"), p.mat=p.mat,  sig.level = 0.05, insig="blank")
corrplot(M, method="color", type="upper", tl.cex=0.6,  col=brewer.pal(n=11, name="PuOr"), p.mat=p.mat,  sig.level = 0.0001,
         insig="blank", tl.col=c(rep("black",21),rep("pink",121)))

# or network 
library("qgraph")

Graph_pcor <- qgraph(M,graph = "pcor", layout = "spring",layout = "spring",  minimum=0.2,
                     groups = c(rep("Cyanobacteria",ntaxa(cyanos)),rep("MOB",nrow(rotsee.mobs.table))), 
                     sampleSize = nrow(t(OTU_table_cooccur)),  pastel = TRUE, posCol = "#003399",
                     negCol = "#FF9933", borders = FALSE, vTrans = 200,
                     details = TRUE                      )





################################################################
## End of the script
################################################################