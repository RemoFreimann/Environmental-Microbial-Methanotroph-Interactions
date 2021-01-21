## =======================================================================================

## methanotroph and bacteria community Analysis

## =======================================================================================
## R 4.0.0
## ---------------------------------------------------------------------------------------
# unlock lines constraint for working folder if necessary 
## ---------------------------------------------------------------------------------------
.trPaths <- paste(paste(Sys.getenv('APPDATA'), '\\Tinn-R\\tmp\\',
            sep=''), c('', 'search.txt', 'objects.txt', 'file.r', 'selection.r',
            'block.r', 'lines.r'), sep='')

## ---------------------------------------------------------------------------------------
# install potentially missing packages
## ---------------------------------------------------------------------------------------
list.of.packages <- c("vegan","phyloseq", "reshape2", "ggplot2", "tcltk", "fields", "RColorBrewer", "mvabund",
                      "doParallel","foreach", "pheatmap", "grDevices", "relaimpo", "multcomp", "car", "MASS",
                      "leaps","RcmdrMisc", "gplots","cooccur", "ggtree","missMDA", "corrplot", "glmulti")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
source("https://bioconductor.org/biocLite.R")  ## try http:// if https:// URLs are not supported
biocLite("multtest", type="source")
biocLite("phyloseq")

## ---------------------------------------------------------------------------------------
# load packages and set working directory
## ---------------------------------------------------------------------------------------
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
library(glmulti)
library(dendextend)
#install.packages("devtools")
#library(devtools)
#install_github("cpauvert/psadd")
library(psadd)
library("ape")
#setwd(tclvalue(tkchooseDirectory()))
# getwd() # double check working directory

options(max.print=2000000000)
## ---------------------------------------------------------------------------------------
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

sample_data(bac.rotsee)  # sample variables
head(get_variable(bac.rotsee)) # variable overview
rank_names(bac.rotsee) # Available Taxonomic ranks: [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
sort(sample_sums(bac.rotsee)) #Reads per sample: 571 to 349919
taxa_names(bac.rotsee)
ntaxa(bac.rotsee)
tax_table(bac.rotsee)
refseq(bac.rotsee)

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
## Check reads distribution
############################################################################################
sample_sum_df <- data.frame(sum = sample_sums(bac.rotsee))
mean(sample_sum_df[,1])
# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "blue", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

############################################################################################
## Filter OTUs that appear less than x times (i.e. read counts) in at least A sample counts
############################################################################################
filt3 <- genefilter_sample(bac.rotsee.prune, filterfun_sample(function(x) x > 50), A=3)   
## more than 50 reads, in at least 3 samples
bac.rotsee.prune.filt3 <- prune_taxa(filt3, bac.rotsee.prune)
ntaxa(bac.rotsee.prune.filt3) 

#############################################################################################
## Create a rarefyed data set
## can be used for UNIFRAC ordinations, Rarefaction curves, etc. 
############################################################################################
bac.rotsee.prune.filt3.rare <- rarefy_even_depth(bac.rotsee.prune.filt3, sample.size = min(sample_sums(bac.rotsee.prune.filt3)),
                                                rngseed = 33, replace = FALSE, trimOTUs = TRUE)
 bac.rotsee.prune.filt3.rare
 sample_sums(bac.rotsee.prune.filt3.rare)

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


#############################################################################################
##  you can check basic informations and get an overall impression of the data 
#############################################################################################
# print(bac.rotsee) # Summary
# ntaxa(bac.rotsee) # Number of OTUs
# sample_data(bac.rotsee)  # Sample variables
# head(get_variable(bac.rotsee)) # Variable overview
# get_variable(bac.rotsee)
# rank_names(bac.rotsee) # Available Taxonomic ranks: [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
# sort(sample_sums(bac.rotsee)) #Reads per sample
# phy_tree(bac.rotsee) # Tree file
# plot_bar(bac.rotsee) # Number of reads as barplot
# plot_bar(bac.rotsee, x= "X.SampleID", y= "Abundance", fill= "Phylum") 
# plot_heatmap(bac.rotsee) # Plot a first heat map of the data
# plot_heatmap(bac.rotsee, method= "NMDS", distance="bray", taxa.label="Phylum" )
# plot_heatmap(bac.rotsee, sample.order="Cu_diss")
# plot_heatmap(bac.rotsee, sample.order="Limnic_Zone", max.label=1100, taxa.label="Class")
# plot_tree(bac.rotsee, color="Depth_m") # Plot a first tree of the data

##############################################################################################
# here we run the script from Joey and Helmut to get some first visual impression of the data set
##############################################################################################

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
                            title = NULL, facet_grid = NULL) 
{
  require(ggplot2)
  require(phyloseq)
  require(plyr)
  require(grid)
  bb <- psmelt(physeq)
  .e <- environment()
  p = ggplot(bb, aes_string(x = x, y = y, 
                            fill = fill), environment = .e)
    p = p + geom_bar(aes(order=desc(bb[,fill])),
                   stat = "identity", 
                   position = "stack",
                   width = 0.4,
                   color = "black") 
  p = p + theme_HB() # Added by Helmut, adjusts theme
  p = p + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
    p = p + guides(fill = guide_legend(override.aes = list(colour = NULL))) 
  p = p + theme(legend.key = element_rect(colour = "black")) 
    p = p + theme(legend.key.size = unit(leg_size, "cm"))
    if (!is.null(facet_grid)) {
        p <- p + facet_grid(facet_grid)
  }
    if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

## extract top abundant OTUs for plotting
Top30 = names(sort(taxa_sums(bac.rotsee), TRUE)[1:30])
top30OTU = prune_taxa(Top30, bac.rotsee)

plot_ordered_bar(top30OTU, "Class", fill = "Genus")
plot_ordered_bar(top30OTU, "Genus", fill = "Genus")

## some bar plots with rarefied dataset
plot_bar(subset_samples(bac.rotsee.prune.filt3.rare, Campaign!=""), fill="Genus")
plot_bar(subset_samples(bac.rotsee.prune.filt3.rare, Campaign!=""), "Class", fill="Genus", facet_grid=~Campaign)

##############################################################################################
## Barplots:  rarefied, absolute (illumina), standardized absolut and 
## standardized to pmoA qPCR 
##############################################################################################

## load either the first block (rarefied) or the second one (absolute read counts) or the 3th block (standarrdized
## absolute counts) or the 4th block (qPCR corrected pmoA) for producing nice Barplots

#block1
#create campaign subsets with rarefied counts
bac.rotsee.prune.filt3.rare_A <- subset_samples(bac.rotsee.prune.filt3.rare, Campaign=="A_June_2013")
bac.rotsee.prune.filt3.rare_A <- prune_taxa(taxa_sums(bac.rotsee.prune.filt3.rare_A) > 0, bac.rotsee.prune.filt3.rare_A)

bac.rotsee.prune.filt3.rare_B <- subset_samples(bac.rotsee.prune.filt3.rare, Campaign=="B_August_2013")
bac.rotsee.prune.filt3.rare_B <- prune_taxa(taxa_sums(bac.rotsee.prune.filt3.rare_B) > 0, bac.rotsee.prune.filt3.rare_B)

bac.rotsee.prune.filt3.rare_C <- subset_samples(bac.rotsee.prune.filt3.rare, Campaign=="C_September_2014")
bac.rotsee.prune.filt3.rare_C <- prune_taxa(taxa_sums(bac.rotsee.prune.filt3.rare_C) > 0, bac.rotsee.prune.filt3.rare_C)

bac.rotsee.prune.filt3.rare_D <- subset_samples(bac.rotsee.prune.filt3.rare, Campaign=="D_September_2015")
bac.rotsee.prune.filt3.rare_D <- prune_taxa(taxa_sums(bac.rotsee.prune.filt3.rare_D) > 0, bac.rotsee.prune.filt3.rare_D)

#block2
#create campaign subsets with absolut counts 
bac.rotsee.prune.filt3.rare_A <- subset_samples(bac.rotsee.prune.filt3, Campaign=="A_June_2013")
bac.rotsee.prune.filt3.rare_A <- prune_taxa(taxa_sums(bac.rotsee.prune.filt3.rare_A) > 0, bac.rotsee.prune.filt3.rare_A)

bac.rotsee.prune.filt3.rare_B <- subset_samples(bac.rotsee.prune.filt3, Campaign=="B_August_2013")
bac.rotsee.prune.filt3.rare_B <- prune_taxa(taxa_sums(bac.rotsee.prune.filt3.rare_B) > 0, bac.rotsee.prune.filt3.rare_B)

bac.rotsee.prune.filt3.rare_C <- subset_samples(bac.rotsee.prune.filt3, Campaign=="C_September_2014")
bac.rotsee.prune.filt3.rare_C <- prune_taxa(taxa_sums(bac.rotsee.prune.filt3.rare_C) > 0, bac.rotsee.prune.filt3.rare_C)

bac.rotsee.prune.filt3.rare_D <- subset_samples(bac.rotsee.prune.filt3, Campaign=="D_September_2015")
bac.rotsee.prune.filt3.rare_D <- prune_taxa(taxa_sums(bac.rotsee.prune.filt3.rare_D) > 0, bac.rotsee.prune.filt3.rare_D)

#block3 (first run block 2)
# create campaign subsets with relative counts (in %)
standf = function(x) (100 * x / sum(x))
bac.rotsee.prune.filt3.rare_A = transform_sample_counts(bac.rotsee.prune.filt3.rare_A, standf) 
bac.rotsee.prune.filt3.rare_B = transform_sample_counts(bac.rotsee.prune.filt3.rare_B, standf) 
bac.rotsee.prune.filt3.rare_C = transform_sample_counts(bac.rotsee.prune.filt3.rare_C, standf) 
bac.rotsee.prune.filt3.rare_D = transform_sample_counts(bac.rotsee.prune.filt3.rare_D, standf) 

#block 4 (first run block 2)
#create campaign subsets with relative counts and pmoA qPCR data
standf = function(x) (x / sum(x))
bac.rotsee.prune.filt3.rare_A_rel = transform_sample_counts(bac.rotsee.prune.filt3.rare_A, standf)
for(n in 1:nsamples(bac.rotsee.prune.filt3.rare_A_rel)){
  otu_table(bac.rotsee.prune.filt3.rare_A)[,n] <- 
    otu_table(bac.rotsee.prune.filt3.rare_A_rel)[,n]*sample_data(bac.rotsee.prune.filt3.rare_A_rel)$pmoA_copies.L[n]}

bac.rotsee.prune.filt3.rare_B_rel = transform_sample_counts(bac.rotsee.prune.filt3.rare_B, standf)
for(n in 1:nsamples(bac.rotsee.prune.filt3.rare_B_rel)){
  otu_table(bac.rotsee.prune.filt3.rare_B)[,n] <- 
    otu_table(bac.rotsee.prune.filt3.rare_B_rel)[,n]*sample_data(bac.rotsee.prune.filt3.rare_B_rel)$pmoA_copies.L[n]}
bac.rotsee.prune.filt3.rare_C_rel = transform_sample_counts(bac.rotsee.prune.filt3.rare_C, standf) 
for(n in 1:nsamples(bac.rotsee.prune.filt3.rare_C_rel)){
  otu_table(bac.rotsee.prune.filt3.rare_C)[,n] <- 
    otu_table(bac.rotsee.prune.filt3.rare_C_rel)[,n]*sample_data(bac.rotsee.prune.filt3.rare_C_rel)$pmoA_copies.L[n]}
bac.rotsee.prune.filt3.rare_D_rel = transform_sample_counts(bac.rotsee.prune.filt3.rare_D, standf) 
for(n in 1:nsamples(bac.rotsee.prune.filt3.rare_A_rel)){
  otu_table(bac.rotsee.prune.filt3.rare_D)[,n] <- 
    otu_table(bac.rotsee.prune.filt3.rare_D_rel)[,n]*sample_data(bac.rotsee.prune.filt3.rare_D_rel)$pmoA_copies.L[n]}



#Plots for top OTU, by depth for each campaign:
(rotsee.tax.table <- tax_table(bac.rotsee))
(genus.numbers <- length(unique(rotsee.tax.table[,"Genus"])))
col.genus.orig<-colorRampPalette(brewer.pal(11, "PuOr"))(genus.numbers)
col.genus.orig<-brewer.pal(genus.numbers, "PuOr")
#col.genus.orig<-(colorRampPalette(brewer.pal(9,"Reds"))(genus.numbers)) 
#col.genus.orig<-(colorRampPalette(brewer.pal(5,"Greens"))(genus.numbers))
#col.genus.orig<-colorRampPalette(c('red','blue','green'))(genus.numbers)

###########################################################################
##Campaign A
#pdf("June 2013 pmoA profile absolut to qpcr.pdf" ,useDingbats=FALSE)     

Top12 = names(sort(taxa_sums(bac.rotsee.prune.filt3.rare_A), TRUE)[1:121])
top12OTU = prune_taxa(Top12, bac.rotsee.prune.filt3.rare_A)
title="Campaign A - June 2013"
p = plot_ordered_bar(top12OTU, x="Depth_m", fill = "Genus", title=title)
p + coord_flip()  +  
  scale_x_reverse(breaks=c(1:16), limits=c(15,0.8)) +
  #ylim(0,3.6E+7) +
  scale_fill_manual(values=col.genus.orig) +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_line(color="grey"),
        #panel.grid.minor = element_line(color="grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        plot.title=element_text(size=15)) 

dev.off()


#pdf("June 2013 pmoA profile relative_ax.pdf" ,useDingbats=FALSE)     

Top12 = names(sort(taxa_sums(bac.rotsee.prune.filt3.rare_A), TRUE)[1:121])
top12OTU = prune_taxa(Top12, bac.rotsee.prune.filt3.rare_A)
title="Campaign A - June 2013"
p = plot_ordered_bar(top12OTU, x="Depth_m", fill = "Genus", title=title)
p + coord_flip()  +  
  scale_x_reverse(breaks=c(1:16), limits=c(15,0.8)) +
  ylim(0,101) +
  scale_fill_manual(values=col.genus.orig) +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_line(color="grey"),
        #panel.grid.minor = element_line(color="grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        plot.title=element_text(size=15)) 

dev.off()


pdf("August 2013 pmoA profile absolut top6.pdf" , useDingbats=FALSE) 
Top12 = names(sort(taxa_sums(bac.rotsee.prune.filt3.rare_A), TRUE)[1:6])
top12OTU = prune_taxa(Top12, bac.rotsee.prune.filt3.rare_A)
p = plot_ordered_bar(top12OTU, x="Depth_m", fill = "Genus", title=title)
p + facet_wrap(~OTU) +
  scale_fill_manual(values=c("Methylobacter"="#B35806", "Methylocystis" ="#F1A340", "Methylomonas"="#FEE0B6",
                             "Methylosoma"="#D8DAEB", "typeIa"="#998EC3", "typeIb"="#542788"))
dev.off()

p + coord_flip()  + scale_x_reverse() +
  scale_fill_brewer(palette="YlOrRd") 
# p + theme(legend.position="none")


###########################################################################
##Campaign B
#pdf("August 2013 pmoA profile  absolut to qpcr.pdf" ,useDingbats=FALSE)  

Top12 = names(sort(taxa_sums(bac.rotsee.prune.filt3.rare_B), TRUE)[1:121])
top12OTU = prune_taxa(Top12, bac.rotsee.prune.filt3.rare_B)
title="August 2013 pmoA profile absolut"
p = plot_ordered_bar(top12OTU, x="Depth_m", fill = "Genus", title=title)
p + coord_flip()  +  
  scale_x_reverse(breaks=c(1:16), limits=c(15,0.8)) +
 # ylim(0,3.6E+7) +
  scale_fill_manual(values=col.genus.orig) +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_line(color="grey"),
        #panel.grid.minor = element_line(color="grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        plot.title=element_text(size=15)) 

dev.off()


#pdf("August 2013 pmoA profile  relative_ax.pdf" ,useDingbats=FALSE)  

Top12 = names(sort(taxa_sums(bac.rotsee.prune.filt3.rare_B), TRUE)[1:121])
top12OTU = prune_taxa(Top12, bac.rotsee.prune.filt3.rare_B)
title="August 2013 pmoA profile absolut"
p = plot_ordered_bar(top12OTU, x="Depth_m", fill = "Genus", title=title)
p + coord_flip()  +  
  scale_x_reverse(breaks=c(1:16), limits=c(15,0.8)) +
  ylim(0,101) +
  scale_fill_manual(values=col.genus.orig) +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_line(color="grey"),
        #panel.grid.minor = element_line(color="grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        plot.title=element_text(size=15)) 

dev.off()


#pdf("August 2013 pmoA profile absolut top6.pdf" , useDingbats=FALSE) 
Top12 = names(sort(taxa_sums(bac.rotsee.prune.filt3.rare_B), TRUE)[1:6])
top12OTU = prune_taxa(Top12, bac.rotsee.prune.filt3.rare_B)
p = plot_ordered_bar(top12OTU, x="Depth_m", fill = "Genus", title=title)
p + facet_wrap(~OTU) +
  scale_fill_manual(values=c("Methylobacter"="#B35806", "Methylocystis" ="#F1A340", "Methylomonas"="#FEE0B6",
                             "Methylosoma"="#D8DAEB", "typeIa"="#998EC3", "typeIb"="#542788"))
dev.off()

p + coord_flip()  + scale_x_reverse() +
  scale_fill_brewer(palette="YlOrRd") 
# p + theme(legend.position="none")


###########################################################################
# Campaign C
#pdf("September 2014 pmoA profile absolut to qpcr.pdf" ,useDingbats=FALSE)  

Top12 = names(sort(taxa_sums(bac.rotsee.prune.filt3.rare_C), TRUE)[1:121])
top12OTU = prune_taxa(Top12, bac.rotsee.prune.filt3.rare_C)
title="Campaign D - September 2014"
p = plot_ordered_bar(top12OTU, x="Depth_m", fill = "Genus", title=title)
p + coord_flip()  +  
  scale_x_reverse(breaks=c(1:16), limits=c(15,0.8)) +
  #ylim(0,3.6E+7) +
  scale_fill_manual(values=col.genus.orig) +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_line(color="grey"),
        #panel.grid.minor = element_line(color="grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        plot.title=element_text(size=15)) 

dev.off()
 ## extract data for caroles presentation
test_glom_caro =  tax_glom(top12OTU, "Genus")
otu_glom_caro = otu_table(test_glom_caro)
tax_glom_caro = tax_table(test_glom_caro)



#pdf("September 2014 pmoA profile  relativ_ax.pdf" ,useDingbats=FALSE)  

Top12 = names(sort(taxa_sums(bac.rotsee.prune.filt3.rare_C), TRUE)[1:121])
top12OTU = prune_taxa(Top12, bac.rotsee.prune.filt3.rare_C)
title="Campaign D - September 2014"
p = plot_ordered_bar(top12OTU, x="Depth_m", fill = "Genus", title=title)
p + coord_flip()  +  
  scale_x_reverse(breaks=c(1:16), limits=c(15,0.8)) +
  ylim(0,101) +
  scale_fill_manual(values=col.genus.orig) +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_line(color="grey"),
        #panel.grid.minor = element_line(color="grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        plot.title=element_text(size=15)) 

dev.off()


pdf("September 2014 pmoA profile rarefied top6.pdf" , useDingbats=FALSE) 
Top12 = names(sort(taxa_sums(bac.rotsee.prune.filt3.rare_C), TRUE)[1:6])
top12OTU = prune_taxa(Top12, bac.rotsee.prune.filt3.rare_C)
p = plot_ordered_bar(top12OTU, x="Depth_m", fill = "Genus", title=title)
p + facet_wrap(~OTU) +
  scale_fill_manual(values=c("Methylobacter"="#B35806", "Methylocystis" ="#F1A340", "Methylomonas"="#FEE0B6",
                             "Methylosoma"="#D8DAEB", "typeIa"="#998EC3", "typeIb"="#542788"))
dev.off()

p + coord_flip()  + scale_x_reverse() +
  scale_fill_brewer(palette="YlOrRd") 
# p + theme(legend.position="none")


###########################################################################
#Campaign D
#pdf("September 2015 pmoA profile absolut to qpcr.pdf" , useDingbats=FALSE)  

Top12 = names(sort(taxa_sums(bac.rotsee.prune.filt3.rare_D), TRUE)[1:121])
top12OTU = prune_taxa(Top12, bac.rotsee.prune.filt3.rare_D)
title="Campaign D - September 2015"
p = plot_ordered_bar(top12OTU, x="Depth_m", fill = "Genus", title=title)
p + coord_flip()  +  
  scale_x_reverse(breaks=c(1:16), limits=c(15,0.8)) +
  #ylim(0,3.6E+7) +
  scale_fill_manual(values=col.genus.orig) +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_line(color="grey"),
        #panel.grid.minor = element_line(color="grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        plot.title=element_text(size=15)) 

dev.off()



#pdf("September 2015 pmoA profile relativ_ax.pdf" , useDingbats=FALSE)  

Top12 = names(sort(taxa_sums(bac.rotsee.prune.filt3.rare_D), TRUE)[1:121])
top12OTU = prune_taxa(Top12, bac.rotsee.prune.filt3.rare_D)
title="Campaign D - September 2015"
p = plot_ordered_bar(top12OTU, x="Depth_m", fill = "Genus", title=title)
p + coord_flip()  +  
  scale_x_reverse(breaks=c(1:16), limits=c(15,0.8)) +
  ylim(0,101) +
  scale_fill_manual(values=col.genus.orig) +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_line(color="grey"),
        #panel.grid.minor = element_line(color="grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        plot.title=element_text(size=15)) 

dev.off()





#pdf("September 2015 pmoA profile absolut top6.pdf" , useDingbats=FALSE) 
Top12 = names(sort(taxa_sums(bac.rotsee.prune.filt3.rare_D), TRUE)[1:6])
top12OTU = prune_taxa(Top12, bac.rotsee.prune.filt3.rare_D)
p = plot_ordered_bar(top12OTU, x="Depth_m", fill = "Genus", title=title)
p + facet_wrap(~OTU) +
scale_fill_manual(values=c("Methylobacter"="#B35806", "Methylocystis" ="#F1A340", "Methylomonas"="#FEE0B6",
                           "Methylosoma"="#D8DAEB", "typeIa"="#998EC3", "typeIb"="#542788"))
dev.off()

p + coord_flip()  + scale_x_reverse() +
  scale_fill_brewer(palette="YlOrRd") 
# p + theme(legend.position="none")




########################################################################################################################################
# this lines are used to run the rest of the manuscript for different subsets (i.e. Campaign or "Allwithout June as it is different")
########################################################################################################################################

# choose one of the lines

bac.rotsee <- bac.rotsee.prune.filt3.seqdepth
bac.rotsee.pmoA = bac.rotsee                             

bac.rotsee <- subset_samples(bac.rotsee.prune.filt3.seqdepth, Campaign=="A_June_2013")

bac.rotsee <- subset_samples(bac.rotsee.prune.filt3.seqdepth, Campaign=="B_August_2013")

bac.rotsee <- subset_samples(bac.rotsee.prune.filt3.seqdepth, Campaign=="C_September_2014")

bac.rotsee <- subset_samples(bac.rotsee.prune.filt3.seqdepth, Campaign=="D_September_2015")

bac.rotsee <- subset_samples(bac.rotsee.prune.filt3.seqdepth, Campaign!="D_June_2015") # without June

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



#########
## or choose the env table with only the variables acutally measured during a specific campaign
# rotsee.env.table.sort=rotsee.env.table  
# 
# #rotsee.env.table=cbind( rotsee.env.table.sort[,c(1:13,15,16,17,37)],(rotsee.env.table.sort[,c(14,18:35,36:42,44:56)]))    # All
# #rotsee.env.table.log=cbind( rotsee.env.table.sort[,c(1:13,15,16,17,37)],log1p(rotsee.env.table.sort[,c(14,18:35,36:42,44:56)]))
# #envvector<-c(17:56)
# 
# (rotsee.env.table=cbind( rotsee.env.table.sort[,c(1:13,15,16,17,37)],(rotsee.env.table.sort[,c(14,18:36,38:58)])))    # All
# rotsee.env.table.log=cbind( rotsee.env.table.sort[,c(1:13,15,16,17,37)],log1p(rotsee.env.table.sort[,c(14,18:36,38:58)]))
# envvector<-c(11,12,14,15,16,17,18,19,20,21,
#              22,23,24,25,26,27,29,32,34,35,37,38,39,40,
#              48,49:56,57,58)
# head(rotsee.env.table[,envvector])
# #rotsee.mobs.env.table=rotsee.env.table
# #rotsee.mobs.env.table.log=rotsee.env.table.log
# 
# 
# rotsee.env.table=cbind(rotsee.env.table.sort[,c(1:13,15,16,17,37)],(rotsee.env.table.sort[,c(14,18:34,36:40,42,44:58)]))    # without June
# rotsee.env.table.log=cbind( rotsee.env.table.sort[,c(1:13,15,16,17,37)],log1p(rotsee.env.table.sort[,c(14,18:34,36:40,42,44:58)]))
# envvector<-c(17:38,40:43,45:55)
# 
# rotsee.env.table=cbind(rotsee.env.table.sort[,c(1:13,15,16,17,37)],(rotsee.env.table.sort[,c(14,18:26,28,31,33,34,36:41,46,48:58)]))    # June 2013
# rotsee.env.table.log=cbind( rotsee.env.table.sort[,c(1:13,15,16,17,37)],log1p(rotsee.env.table.sort[,c(14,18:26,28,31,33,34,36:41,46,48:58)]))
# envvector<-c(17:38,40:43,45:55)
# 
# rotsee.env.table=cbind( rotsee.env.table.sort[,c(1:13,15,16,17,37)],(rotsee.env.table.sort[,c(14,18:36,38:42,44:56)]))     # August 2013
# rotsee.env.table.log=cbind(rotsee.env.table.sort[,c(1:13,15,16,17,37)],log1p(rotsee.env.table.sort[,c(14,18:36,38:42,44:56)]))
# envvector<-c(17:38,40:43,45:55)
# 
# rotsee.env.table=cbind(rotsee.env.table.sort[,c(1:13,15,16,17,37)],(rotsee.env.table.sort[,c(14,18:40,42:56)]))    # September 2014
# rotsee.env.table.log=cbind( rotsee.env.table.sort[,c(1:13,15,16,17,37)],log1p(rotsee.env.table.sort[c(14,18:40,42:56)]))
# (envvector<-c(17:47))
# 
# rotsee.env.table=cbind(rotsee.env.table.sort[,c(1:13,15,16,17,37)],(rotsee.env.table.sort[,c(14,18:34,36:45,47:56)]))    # September 2015
# rotsee.env.table.log=cbind( rotsee.env.table.sort[,c(1:13,15,16,17,37)],log1p(rotsee.env.table.sort[,c(14,18:36,38:56)]))
# envvector<-c(14:40,42,44:55)



## you can view the tables in a new window:
#            utils::View(rotsee.tax.table)

## you can write the tables into csv files for further analysis in xcel:
#            write.csv(rotsee.otu.table, file ="rotsee_otu_table.csv")



########################################################################################################################################
# plot diversity measures
########################################################################################################################################

richness.pmoA = estimate_richness(bac.rotsee.prune.filt3, split = TRUE, measures = c("Observed",  "chao", "Shannon", "Simpson", "InvSimpson", "Fisher") )
richness.pmoA = cbind(richness.pmoA, rotsee.env.table[,c("Depth_m", "Campaign")])
library(reshape2)
plot_richness(bac.rotsee.prune.filt3, x="Depth_m", measures=c("Observed", "Chao1", "Shannon", "Simpson", "Fisher", "InvSimpson"))  


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

## Bray-Curtis distances between samples
dis <- vegdist(t(otu_table(bac.rotsee.prune.filt3)))
## Calculate multivariate dispersions
envbetadisp_oxi= get_variable(bac.rotsee.prune.filt3)$Oxidation_Zone
envbetadisp_campaign=get_variable(bac.rotsee.prune.filt3)$Campaign

#######################################
## compare betadiversity of oxyc zones between campaigns
envbetadisp_oxi_campaign=interaction(envbetadisp_oxi, envbetadisp_campaign) # oxidation zone * campaign

mod <- betadisper(dis, group=envbetadisp_oxi_campaign,
                  type="median", bias.adjust=TRUE)
mod

## Perform test
anova(mod)
permutest(mod, pairwise = TRUE, permutations = 99)
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)
which(mod.HSD$group[,4]<=0.05)



## with data ellipses instead of hulls
plot(mod, ellipse = TRUE, hull = FALSE) # 1 sd data ellipse
plot(mod, ellipse = TRUE, hull = FALSE, conf = 0.90) # 90% data ellipse

## Draw a boxplot of the distances to centroid for each group
boxplot(mod)

#####################
## additional boxplot
betadistances=mod$distances
betagroups=mod$group
beta_anova=data.frame(betadistances, (betagroups))
beta_anova$betadistances_logit=car::logit(betadistances, percents=FALSE)

## between zones and seasons

barorder=c("Oxycline.A_June_2013" , "Oxidation_Zone.A_June_2013",
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


p <- ggplot(beta_anova, aes(x = X.betagroups., y = betadistances)) + 
  geom_boxplot( aes(colour=X.betagroups.),
                alpha =0.8,
                cex=1,
                width=0.7) +
  
  scale_x_discrete(name="", limits = barorder,
                   labels=c("Oxycline June13", "Oxidation Zone June13",
                            "Oxic August13", "Oxycline August13", "Oxidation Zone August13", "Anoxic August13",
                            "Oxic September14", "Oxycline September14", "Oxidation Zone September14", "Anoxic September14",
                            "Oxic September15", "Oxidation Zone September15", "Anoxic September15")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = c(2.5, 6.5, 10.5)) +
  scale_y_continuous(name="Betadiversity", expand=c(0,0), breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6, 0.7, 0.8, 0.9) , 
                     limits=c(0, 0.95) ) +
  annotate("text", label=paste(plot.posthoc.annot[,1]), y= 0.05 , x= c(1:13), cex=4)
p


#######################################
## only between zones


envbetadisp_oxi_campaign=interaction(envbetadisp_oxi) # oxidation Zone

#######################################
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

## boxplot
p <- ggplot(beta_anova, aes(x = X.betagroups., y = betadistances)) + 
  geom_boxplot( aes(colour=X.betagroups.),
                alpha =0.8,
                cex=1,
                width=0.7) +
  
  scale_x_discrete(name="", limits = barorder,
                   labels=c("Oxic", "Oxycline" , "Oxidation_Zone", "Anoxic")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +
  scale_y_continuous(name="Betadiversity", expand=c(0,0), breaks = c(0.2,0.3,0.4,0.5,0.6, 0.7, 0.7, 0.9) , 
                     limits=c(0.19, 0.95) ) +
  annotate("text", label=paste(plot.posthoc.annot[,1]), y= 0.92 , x= c(1:4), cex=4)
p

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
    
## you can use this function to write the produced OTU table into your working directory
# write.csv(rotsee.vegan.otu.table, "vegan_test.csv")





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
library(Biostrings)
refseq_filt3 = refseq(bac.rotsee.prune.filt3)
refseq_filt3@ranges@NAMES <- paste(mgsub(c("\\.", "/","-"),c("","_","_"),(refseq_filt3@ranges@NAMES)),"_",colnames(rotsee.vegan.otu.table), sep="") # replace characters that  are not parsed propperly
## write the fasta files of the prefiltered OTUs
writeXStringSet(refseq_filt3,file="seq_prune_filt3.fas",format="fasta") 


## this can be used for the phylocom analysis
# comm.pmoA =otu_table(bac.rotsee)
# (newtaxanames <- paste(colnames(rotsee.vegan.otu.table),"_",rownames(comm.pmoA), sep=""))
# bac.rotsee.phytree=bac.rotsee
# taxa_names(bac.rotsee.phytree) <- newtaxanames
# phy.tree.pmoA =phy_tree(bac.rotsee.phytree)
# write.tree(phy.tree.pmoA, file = "testtreepmoA")




########################################################################################################################################
########################################################################################################################################
## "in-depth" NMDS analysis
########################################################################################################################################
########################################################################################################################################


########################################################################################################################################
## NMDS
set.seed(77)
rotsee.mds.allOTUs.pmoa <- rotsee.mds.allOTUs < -metaMDS(rotsee.vegan.otu.table, distance="bray", k=2)
## some analytics of the quality of the NMDS 
stressplot(rotsee.mds.allOTUs)
(goodpoint=goodness(rotsee.mds.allOTUs)) #sum of squared values is equal to squared stress

## you can rotate the MDS according to represent the depth succession on the first axis (or any other paramter...)
#rotsee.mds.allOTUs <- with(rotsee.env.table, MDSrotate(rotsee.mds.allOTUs, -c(Depth_m)))



########################################################################################################################################
## define some colors for the plot
col.depth<-(colorRampPalette(brewer.pal(5,"Blues"))(1000)) # depth scale
col.fact1<-(colorRampPalette(brewer.pal(5,"Greens"))(1000))
col.fact2<-(colorRampPalette(brewer.pal(5,"Reds"))(1000))
YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
col.fact3<-(colorRampPalette(YlOrBr, space="Lab"))


# for less otus (91)
genus.numbers <- length(unique(rotsee.tax.table[,"Genus"]))
#col.genus<-(colorRampPalette(brewer.pal(genus.numbers,"Reds"))(genus.numbers))  # choose this or the next color
col.genus.orig<-brewer.pal(genus.numbers, "PuOr")
genus.df<-transform(rotsee.tax.table, id=match(Genus, unique(Genus)))
test=unique(genus.df$Genus)
match(test, sort(unique(genus.df[,"Genus"]))) # copy the vector to the col.genus.orig
col.genus<-col.genus.orig[c( 2, 6, 3, 1, 5, 4)]
#col.genus<-col.genus.orig[c( 6, 3, 2, 7, 4, 1, 5)]

# for more otus >1000
# col.genus.orig<-colorRampPalette(brewer.pal(11, "PuOr"))(genus.numbers)
# genus.df<-transform(rotsee.tax.table, id=match(Genus, unique(Genus)))
# test=unique(genus.df$Genus)
# match(test, sort(unique(genus.df[,"Genus"]))) # copy the vector to the col.genus.orig
# col.genus<-col.genus.orig[c(  13, 12,  5, 11,  8, 14,  1,  3,  4,  2,  6,  7, 15 , 9 ,10)]

# match to barplots color order from ordered_bar_plots
#target= c("Methylobacter", "Methylocystis", "Methylomonas", "Methylosoma", "typeIa", "typeIb")
#genus.df=genus.df[order(genus.df$Genus),]  



# vector for species point sizes according to raw read counts
bac.rotsee.prune.filt3.otu = otu_table(bac.rotsee.prune.filt3)
min(taxa_sums(bac.rotsee.prune.filt3.otu))
max(taxa_sums(bac.rotsee.prune.filt3.otu))
species.size.vec <- taxa_sums(bac.rotsee.prune.filt3.otu)
x=species.size.vec
normalized.species.size.vec = sqrt((x-min(x))/(max(x)-min(x))*3+0.5)




########################################################################################################################################
## basic nmds plot... will be used for further analysis.... use this code block whenever you want to reset and add stuff from below

#pdf("pmoa NMDS oxzones September 2015.pdf" ,useDingbats=FALSE)      ## if you use this command, there will be a pdf produced of anything which gets printed to device below...until you type "dev.off()"

#par(mar=c(4.1, 4.1, 4.1, 7.1), xpd=TRUE)

# resetPar <- function() {
#  dev.new()
#  op <- par(no.readonly = TRUE)
#  dev.off()
#  op
# }
# par(resetPar())   # you can run this function / line to reset par values to default
plot(rotsee.mds.allOTUs, type="n")
plot(rotsee.mds.allOTUs, type="n", xlim=c(-1.5, 2.5))
plot(rotsee.mds.allOTUs, type="n", xlim=c(-2.2, 0.9), ylim=c(-1.2, 1))


#points(rotsee.mds.allOTUs, display ="species", cex=0.1, pch=21, col="orange", bg=alpha("orange", 0.00001), lwd=0)
points(rotsee.mds.allOTUs, display ="species", cex=normalized.species.size.vec * 2 ,#cex=1.5,
       pch=21, col=col.genus[genus.df$id], bg=alpha(col.genus[genus.df$id], 0.6), lwd=0)
legend("topright",legend=sort(unique(genus.df[,"Genus"])),pch=c(19),col=alpha(col.genus.orig, 0.8), bty="n")

#points(rotsee.mds.allOTUs, display ="sites", cex=c(goodpoint)*400, alpha= 2, pch=as.numeric(rotsee.env.table$Campaign)+20, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
#       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.5))

points(rotsee.mds.allOTUs, display ="sites", cex=3.5, alpha= 2, pch=as.numeric(as.factor(rotsee.env.table$Campaign))+21, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.5))
text(rotsee.mds.allOTUs$points[, 1]+0.05, rotsee.mds.allOTUs$points[, 2]+0.05,
     labels=unlist(dimnames(rotsee.mds.allOTUs$points)[1]),
     col="darkblue",#col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
     cex=0.7)
legend(0,-0.3,legend=as.expression(paste("Stress:\n",
                                          round(c(rotsee.mds.allOTUs$stress)*100,2),
                                          "%")),cex=1,bty="n")


image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.3,0.7), axis.args = list(cex.axis = 0.7),
           zlim=range(-(rotsee.env.table$Depth_m)) )

#identify(rotsee.mds.allOTUs$species, labels=(rotsee.tax.table@.Data[,"Genus"]), cex=0.6, col="grey" )



 dev.off()

########################################################################################################################################
## plot OTU's on the total community plot that are highly abundant


plotabundat<- genefilter_sample(bac.rotsee, filterfun_sample(function(nrotu) nrotu > 1000), A=10) # larger 10000 reads at least in 10 samples
(abundantvect = as.numeric(which(plotabundat)))

text(rotsee.mds.allOTUs$species[abundantvect ,1], rotsee.mds.allOTUs$species[abundantvect ,2],
     colnames(rotsee.vegan.otu.table)[abundantvect ], cex=0.7, col="black")

## or assign with OTU

     text(rotsee.mds.allOTUs$species[abundantvect,1], rotsee.mds.allOTUs$species[abundantvect,2]+0.1,
     colnames(rotsee.otu.table.t)[abundantvect], col="darkgrey", cex = 0.7)


########################################################################################################################################     
## additional graphical elements such as environmental fittin or GAMs


## environmenal fitting (this is a linear fitting on the plot, it is usually better to use the GAM (Generalized Additive Model)
## fitting from further below)
     
(mds.allOTUenvfit <- envfit(rotsee.mds.allOTUs, rotsee.env.table[,-c(1:10,12)], na.rm=TRUE, permutations=1000))
plot(mds.allOTUenvfit, choices = c(1,2), axis = FALSE, p.max = 0.05, col = "darkred",  add = TRUE,
     cex=0.7)
## most of the variables are significantly correlated with the OTU structure. This is due to the "congruent" structuring throughout the depth profile
## also note that half of the samples were kicked out as there are missing values(influencing the varibables that actually have no missing values as well).



## OTU fitting, here you can aposteriori fitt OTUs to check out how well they're representated on the NMDS
#(mds.allOTUfit <- envfit(rotsee.mds.allOTUs, rotsee.vegan.otu.table[,1:10], na.rm=TRUE, permutations=100))
#plot(mds.allOTUfit, choices = c(1,2), axis = FALSE, p.max = 0.01, col = "lightgrey",  add = TRUE, cex=0.8)




## Confidence ellipses for different factors, when they overlap strongly, there might not be a difference between the OTU structure of the relative groupings
with(rotsee.env.table, ordiellipse(rotsee.mds.allOTUs, groups= interaction(Oxidation_Zone), display="sites",
                                draw="polygon", col=c( "darkblue","green","lightblue","orange"), border= "lightgrey", alpha=20 , label=FALSE))
## ordination hulls
with(rotsee.env.table, ordihull(rotsee.mds.allOTUs, groups= interaction(Oxidation_Zone), display="sites",
                                   draw="polygon", col=c("darkblue","green","lightblue","orange"), border= "lightgrey", alpha=20,
                                   label=FALSE))



## GAMS for environmental variables can be plotted on the ordination                                   here you can change the parameter
gam.PhysicoChemRaw<-ordisurf(as.formula(paste("rotsee.mds.allOTUs ~ ",paste(colnames(rotsee.env.table[parame<-"TDN"])))),rotsee.env.table,
                             choices=c(1,2), knots=10,                                         
                             family="gaussian", isotropic= TRUE,
                             col="black",scaling=3,
                             add = TRUE, display = "sites",
                             nlevels = 15, labcex = 0.8,
                             bubble = TRUE, cex = 1, select = TRUE, method = "REML",
                             gamma = 1, plot = TRUE, lwd=0.5)
(summary.gam.PhysicoChemRaw <- summary(gam.PhysicoChemRaw) ) # this summary shows how well a variable correlates with the OTU structure

legend(-0.65,0.9,legend=as.expression(paste("Explained deviance of\n",parame,":\n",
                                           round(c(summary.gam.PhysicoChemRaw$dev.expl)*100,2),
                                           "%,   ","P=",
                                           round(c(summary.gam.PhysicoChemRaw$s.pv[[1]]),3))),cex=0.8,bty="n")



## plot GAMs for specific depth range, this can be interesting, to within which depth range actual correlations occure ( I might program a function for this...)

lowerlimit= 9
upperlimit= 4

# above lower limit [m]
rotsee.mds.allOTUs.8<-scores(rotsee.mds.allOTUs)[which(rotsee.env.table$Depth_m <= lowerlimit ),c(1,2)]
rotsee.env.table.8<-rotsee.env.table[which(rotsee.env.table$Depth_m <= lowerlimit),]

gam.PhysicoChemRaw<-ordisurf(rotsee.mds.allOTUs.8[,c(1,2)],  rotsee.env.table.8[,parame<-"Cu_DGT_nM"],
                             choices=c(1,2), knots=10,                                         
                             family="gaussian", isotropic= TRUE,
                             col="orange",scaling=3,
                             add = TRUE, display = "sites",
                             nlevels = 10, labcex = 0.8,
                             bubble = TRUE, cex = 1, select = TRUE, method = "REML",
                             gamma = 1, plot = TRUE, lwd=0.5)
(summary.gam.PhysicoChemRaw <- summary(gam.PhysicoChemRaw) )

legend(-0.65,0.9,legend=as.expression(paste("Explained deviance of\n",parame,"\n lower water column:\n",
                                           round(c(summary.gam.PhysicoChemRaw$dev.expl)*100,2),
                                           "%,   ","P=",
                                           round(c(summary.gam.PhysicoChemRaw$s.pv[[1]]),3))),cex=0.8,bty="n")



# between upper and lower limit
rotsee.mds.allOTUs.8<-scores(rotsee.mds.allOTUs)[which(rotsee.env.table$Depth_m <= lowerlimit & rotsee.env.table$Depth_m > upperlimit),c(1,2)]
rotsee.env.table.8<-rotsee.env.table[which(rotsee.env.table$Depth_m <= lowerlimit & rotsee.env.table$Depth_m > upperlimit),]

gam.PhysicoChemRaw<-ordisurf(rotsee.mds.allOTUs.8[,c(1,2)],  rotsee.env.table.8[,parame<-"Cu_DGT_nM"],
                             choices=c(1,2), knots=nrow(rotsee.env.table.8)-1,                                         
                             family="gaussian", isotropic= TRUE,
                             col="red",scaling=3,
                             add = TRUE, display = "sites",
                             nlevels = 10, labcex = 0.8,
                             bubble = TRUE, cex = 1, select = TRUE, method = "REML",
                             gamma = 1, plot = TRUE, lwd=0.5)
(summary.gam.PhysicoChemRaw <- summary(gam.PhysicoChemRaw) )

legend(-0.65,0.9,legend=as.expression(paste("Explained deviance of\n",parame,"\n between 3m and 8m:\n",
                                           round(c(summary.gam.PhysicoChemRaw$dev.expl)*100,2),
                                           "%,   ","P=",
                                           round(c(summary.gam.PhysicoChemRaw$s.pv[[1]]),3))),cex=0.8,bty="n")




## plot bar for depth with locater
#position<-locator(n=1)
#position$x= 1.5
#$y= -1.3
#colorbar.plot(position$x, position$y, strip=c(1:1000), strip.width = 0.05, strip.length = 0.1, 
#              zrange = c(1,10), adj.x = 0.5, adj.y = 0.5, col = rev(col.depth), 
#              horizontal = FALSE)
#text(position$x,position$y,"depth")



## make pdf for gams of all environmental variables (will be safed in working directory)
#---------------------------------------------------------------------------------------

# here choose the env variables without NAs for specific campaign 
# normally use the envvector defined for the specific campaign before
head(rotsee.env.table[,-c(1:13,15:17)])

gamvector=c(1:length(rotsee.env.table))
gamvector=gamvector[-c(1:13,15:17, 48:67 )]  # remove factors

sink(file = "GAMs on pmoA.txt", append = FALSE, type = c("output"),
     split = FALSE)

for (i in c(gamvector)){
  
  
  pdf(paste(colnames(rotsee.env.table[i]), " GAM fitting",".pdf", sep=""),useDingbats=FALSE) #, width=4.5, height=4)
  
  
  plot(rotsee.mds.allOTUs, type="n", xlim=c(-1.5, 2.5))
  col.depth<-(colorRampPalette(brewer.pal(5,"Blues"))(1000))
  
  #points(rotsee.mds.allOTUs, display ="species", cex=0.1, pch=21, col="orange", bg=alpha("orange", 0.00001), lwd=0)
  points(rotsee.mds.allOTUs, display ="species", cex=normalized.species.size.vec * 2 ,#cex=1.5,
         pch=21, col=col.genus[genus.df$id], bg=alpha(col.genus[genus.df$id], 0.5), lwd=0)
  legend("topright",legend=sort(unique(genus.df[,"Genus"])),pch=c(19),col=alpha(col.genus.orig, 0.99), bty="n")
  
  #points(rotsee.mds.allOTUs, display ="sites", cex=c(goodpoint)*400, alpha= 2, pch=as.numeric(rotsee.env.table$Campaign)+20, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
  #       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.5))
  
  points(rotsee.mds.allOTUs, display ="sites", cex=3.5, alpha= 2, pch=as.numeric(as.factor(rotsee.env.table$Campaign))+21, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
         bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.7))
  text(rotsee.mds.allOTUs$points[, 1]+0.05, rotsee.mds.allOTUs$points[, 2]+0.05,
       labels=unlist(dimnames(rotsee.mds.allOTUs$points)[1]),
       col="darkblue",#col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
       cex=0.7)
  legend(-1,-1.3,legend=as.expression(paste("Stress:\n",
                                            round(c(rotsee.mds.allOTUs$stress)*100,2),
                                            "%")),cex=1,bty="n")
  
  
  image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.3,0.7), axis.args = list(cex.axis = 0.7),
             zlim=range(-(rotsee.env.table$Depth_m)) )
  
  

  
 
  
  gam.PhysicoChemRaw<-ordisurf(as.formula(paste("rotsee.mds.allOTUs ~ ",paste(colnames(rotsee.env.table[i])))),rotsee.env.table,
                               choices=c(1,2), knots=10,                                         
                               family="gaussian", isotropic= TRUE,
                               col="darkgrey",scaling=3,
                               add = TRUE, display = "sites",
                               nlevels = 15, labcex = 0.8,
                               bubble = TRUE, cex = 1, select = TRUE, method = "REML",
                               gamma = 1, plot = TRUE)
  (summary.gam.PhysicoChemRaw <- summary(gam.PhysicoChemRaw) )
  print(colnames(rotsee.env.table[i]))
  print(summary.gam.PhysicoChemRaw <- summary(gam.PhysicoChemRaw) )
  
  legend(0.65,-1.3,legend=as.expression(paste("Explained deviance of\n ", colnames(rotsee.env.table[i]),
                                            round(c(summary.gam.PhysicoChemRaw$dev.expl)*100,2),
                                            "%,   ","P=",
                                            round(c(summary.gam.PhysicoChemRaw$s.pv[[1]]),3))),cex=0.8,bty="n")

  
  dev.off()

}

sink()

#########################################################################################################################
## here you can interactively identify potential bacterial interactors on the plot


plot(rotsee.mds.allOTUs, type="n")

#points(rotsee.mds.allOTUs, display ="species", cex=0.1, pch=21, col="orange", bg=alpha("orange", 0.00001), lwd=0)
points(rotsee.mds.allOTUs, display ="species", cex=normalized.species.size.vec * 2 ,#cex=1.5,
       pch=21, col=col.genus[genus.df$id], bg=alpha(col.genus[genus.df$id], 0.8), lwd=0)
legend("topright",legend=sort(unique(genus.df[,"Genus"])),pch=c(19),col=alpha(col.genus.orig, 0.8), bty="n")

#points(rotsee.mds.allOTUs, display ="sites", cex=c(goodpoint)*400, alpha= 2, pch=as.numeric(rotsee.env.table$Campaign)+20, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
#       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.5))

points(rotsee.mds.allOTUs, display ="sites", cex=3.5, alpha= 2, pch=as.numeric(as.factor(rotsee.env.table$Campaign))+21, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.5))
text(rotsee.mds.allOTUs$points[, 1]+0.05, rotsee.mds.allOTUs$points[, 2]+0.05,
     labels=unlist(dimnames(rotsee.mds.allOTUs$points)[1]),
     col="darkblue",#col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
     cex=0.7)
legend(-1,-0.7,legend=as.expression(paste("Stress:\n",
                                          round(c(rotsee.mds.allOTUs$stress)*100,2),
                                          "%")),cex=1,bty="n")


image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.3,0.7), axis.args = list(cex.axis = 0.7),
           zlim=range(-(rotsee.env.table$Depth_m)) )




## by hand, select points on the plot and it will annotate them accordingly
identify(rotsee.mds.allOTUs$species, labels=rownames(rotsee.mds.allOTUs$species), cex=0.6, col="grey" ) # by deepest phylogenetic annotation
identify(rotsee.mds.allOTUs$species, labels=rownames(rotsee.otu.table), cex=0.6, col="grey" )           # by OTU number
identify(rotsee.mds.allOTUs$species, labels=(rotsee.tax.table@.Data[,"Genus"]), cex=0.6, col="grey" )


## mark location of species you're looking for 
## phylogenetic search
(bacterianumber <- which(rownames(rotsee.mds.allOTUs$species)==c("deep-sea-1")))  
## phylogenetic search
(bacterianumber <- which((rotsee.tax.table[,"Genus"])==c("typeIb")))
## OTU numbers search
(bacterianumber <- which(colnames(rotsee.otu.table.t)==("ZOTU2055")))  
## mark location
points(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
       labels=rownames(rotsee.mds.allOTUs$species)[bacterianumber], cex=2, pch=25,col="black", bg=alpha("yellow", alpha=0.5))
## annotate them by phylogenetic or OTU annotation
text(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2], 
     labels=rownames(rotsee.mds.allOTUs$species)[bacterianumber], cex=0.8)
text(rotsee.mds.allOTUs$species[bacterianumber,1], rotsee.mds.allOTUs$species[bacterianumber,2],
     labels=rownames(rotsee.otu.table)[bacterianumber], cex=0.8)

#########################################################################################################################
#########################################################################################################################
##  Similarity Percentages (SIMPER) -> i.e., finds the most distinguished OTUs for the three light zones
#########################################################################################################################
#########################################################################################################################

## this one I would run for the different seasons independently...
#(sim<- with(rotsee.env.table, simper(rotsee.vegan.otu.table[,], interaction(Light_Zone))))    # annotted with highest Taxon

(sim<- with(rotsee.env.table, simper(t(rotsee.otu.table)[,], interaction(Oxidation_Zone))))      # annotted with OTU number, use this one for plotting points on NMDS


str(sim)

sum(sim$Oxycline_Oxic$average)# total between group variance explained by OTUs
sum(sim$Oxycline_Oxidation_Zone$average) 
sum(sim$Oxidation_Zone_Anoxic$average)

sum(sim$Oxycline_Anoxic$average)
sum(sim$Oxidation_Zone_Oxic$average)
sum(sim$Oxic_Anoxic$average)

## extract most 6 influencal OTUS for the separation of the specific oxidation zones (i.e. between oxic-oxycline, oxycline-oxidation zone, oxidation zone-anoxic)
(simOTUs.Oxic_Oxycline = sim$Oxycline_Oxic$species[order(-sim$Oxycline_Oxic$average)][1:6])
(simOTUs.Oxycline_Oxidation_Zone = sim$Oxycline_Oxidation_Zone$species[order(-sim$Oxycline_Oxidation_Zone$average)][1:6]) 
(simOTUs.Oxidation_Zone_Anoxic = sim$Oxidation_Zone_Anoxic$species[order(-sim$Oxidation_Zone_Anoxic$average)][1:6])

## plot the most influencal speices on the nmds 
## the most important OTUs according to the simper output are now converted in taxa name index vector
simOTUs<-c(simOTUs.Oxic_Oxycline,
           simOTUs.Oxycline_Oxidation_Zone,
           simOTUs.Oxidation_Zone_Anoxic
           )   
simOTUsnumbers=c()
for (i in c(1:length(simOTUs))){

  simOTUsnumbers[i]<- which(colnames(t(rotsee.otu.table)) == simOTUs[i])
}

simOTUsnumbers # vector with the position of simOTUs in the vegan table


## mds plot
plot(rotsee.mds.allOTUs, type="n")
col.depth<-(colorRampPalette(brewer.pal(5,"Blues"))(1000))

#points(rotsee.mds.allOTUs, display ="species", cex=0.1, pch=21, col="orange", bg=alpha("orange", 0.00001), lwd=0)
points(rotsee.mds.allOTUs, display ="species", cex=normalized.species.size.vec * 2 ,#cex=1.5,
       pch=21, col=col.genus[genus.df$id], bg=alpha(col.genus[genus.df$id], 0.5), lwd=0)
legend("topright",legend=sort(unique(genus.df[,"Genus"])),pch=c(19),col=alpha(col.genus.orig, 0.99), bty="n")

#points(rotsee.mds.allOTUs, display ="sites", cex=c(goodpoint)*400, alpha= 2, pch=as.numeric(rotsee.env.table$Campaign)+20, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
#       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.5))

points(rotsee.mds.allOTUs, display ="sites", cex=3.5, alpha= 2, pch=as.numeric(as.factor(rotsee.env.table$Campaign))+21, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.7))
text(rotsee.mds.allOTUs$points[, 1]+0.05, rotsee.mds.allOTUs$points[, 2]+0.05,
     labels=unlist(dimnames(rotsee.mds.allOTUs$points)[1]),
     col="darkblue",#col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
     cex=0.7)
legend(-1,-0.7,legend=as.expression(paste("Stress:\n",
                                          round(c(rotsee.mds.allOTUs$stress)*100,2),
                                          "%")),cex=1,bty="n")

image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.3,0.7), axis.args = list(cex.axis = 0.7),
           zlim=range(-(rotsee.env.table$Depth_m)) )

## plot the selected Simper OTU's on the total community plot
points(rotsee.mds.allOTUs$species[simOTUsnumbers,1], rotsee.mds.allOTUs$species[simOTUsnumbers,2], cex=2, pch=25, col="yellow", bg=alpha("red", 0.1))
## annotate by highest highest taxa
text(rotsee.mds.allOTUs$species[simOTUsnumbers,1]*1.1, rotsee.mds.allOTUs$species[simOTUsnumbers,2]*1.1,
     colnames(rotsee.vegan.otu.table)[simOTUsnumbers], col="darkred", cex=0.7)
## or  annotate by OTUs
text(rotsee.mds.allOTUs$species[simOTUsnumbers,1]*1.1, rotsee.mds.allOTUs$species[simOTUsnumbers,2]*1.1,
     colnames(rotsee.otu.table.t)[simOTUsnumbers], col="darkred")

#########################################################################################################################
########################################################################################################################
## vabund --> check this as better alternative for PERMANOVA (ADONiS) and SIMPER, use absolut count data
#########################################################################################################################
#########################################################################################################################
xx=(otu_table(bac.rotsee))
xx=t(xx@.Data)
rotsee.glmdata <- mvabund(decostand(xx, "pa"))

meanvar.plot(rotsee.vegan.otu.table)
rotsee.glmdata <- mvabund(t(otu_table(bac.rotsee)))  # check which one to use (absolut count data frame, or presence/absence as above)



(env.glmdata<- rotsee.env.table [,] )  # here it might be good to transform environmental data as performed in the next line
(env.glmdata<- rotsee.env.table.log) 

(rotsee.glm<- manyglm(rotsee.glmdata ~ Oxidation_Zone , data=env.glmdata, family="negative.binomial",    #choose (multiple) factors for the model
                     cor.type="I", test="LR",
                     show.coef=FALSE, show.fitted=TRUE, show.residuals=FALSE ))

#(rotsee.glm.null<- manyglm(rotsee.glmdata ~ 1 , data=env.glmdata, family="negative.binomial",    #null model
#                      cor.type="I", test="LR",
#                      show.coef=FALSE, show.fitted=TRUE, show.residuals=FALSE ))

#anova(rotsee.glm.null, rotsee.glm, nBoot=1)

plot(rotsee.glm)

(rotsee.glm.summary<-summary(rotsee.glm, resamp="residual", nBoot=100, test="LR"))    ## notate Likelihood ratio and values
(rotsee.glm.anova<- anova(rotsee.glm, p.uni="none", test="LR",        #anova.manyglm   -> are liminc zones different? yes...   
                        show.time="all", nBoot=100))        # for serious run choose 1000 bootstraps!!
summary(rotsee.glm.anova)

rotsee.glm.anova.single.sp = anova(rotsee.glm, p.uni="adjusted")  ## for single species


## check out which species are affected mostly
rotsee.glm.anova$uni.p # 
s = sort(rotsee.glm.anova$uni.test[2,],index.return=T,decreasing=T)
s$x[1:10]/sum(s$x)         #the proportions of total test stat due to these 10 OTUs
sum(s$x[1:10]/sum(s$x) )   # this is the % of variation contributed by this 10 OTUs to the total 

(glm.imp.names <- names(s$x[1:10]))    # get the 10 bacteria which are most distinct beween tested groups/factor
vabundOTUsnumbers=c()
for (i in c(1:length(glm.imp.names))){
  
  vabundOTUsnumbers[i]<- which(colnames(t(rotsee.otu.table)) == glm.imp.names[i])
}



## test which variables are most influencual ()
X=env.glmdata[,c("Oxidation_Zone", "Campaign"), drop=FALSE]
X=as.matrix(X, drop=FALSE)
best.r.sq(rotsee.glmdata~ X, data=X, R2="h")

## mds plot
plot(rotsee.mds.allOTUs, type="n")
col.depth<-(colorRampPalette(brewer.pal(5,"Blues"))(1000))

#points(rotsee.mds.allOTUs, display ="species", cex=0.1, pch=21, col="orange", bg=alpha("orange", 0.00001), lwd=0)
points(rotsee.mds.allOTUs, display ="species", cex=normalized.species.size.vec * 2 ,#cex=1.5,
       pch=21, col=col.genus[genus.df$id], bg=alpha(col.genus[genus.df$id], 0.5), lwd=0)
legend("topright",legend=sort(unique(genus.df[,"Genus"])),pch=c(19),col=alpha(col.genus.orig, 0.99), bty="n")

#points(rotsee.mds.allOTUs, display ="sites", cex=c(goodpoint)*400, alpha= 2, pch=as.numeric(rotsee.env.table$Campaign)+20, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
#       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.5))

points(rotsee.mds.allOTUs, display ="sites", cex=3.5, alpha= 2, pch=as.numeric(as.factor(rotsee.env.table$Campaign))+21, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
       bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.7))
text(rotsee.mds.allOTUs$points[, 1]+0.05, rotsee.mds.allOTUs$points[, 2]+0.05,
     labels=unlist(dimnames(rotsee.mds.allOTUs$points)[1]),
     col="darkblue",#col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
     cex=0.7)
legend(-1,-0.7,legend=as.expression(paste("Stress:\n",
                                          round(c(rotsee.mds.allOTUs$stress)*100,2),
                                          "%")),cex=1,bty="n")


image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.3,0.7), axis.args = list(cex.axis = 0.7),
           zlim=range(-(rotsee.env.table$Depth_m)) )

## Add the vabund OTUs 
#-------------------------------------------
points(rotsee.mds.allOTUs$species[vabundOTUsnumbers,1], rotsee.mds.allOTUs$species[vabundOTUsnumbers,2],cex=2, pch=23, col="yellow", bg=alpha("brown", 0.2))
text(rotsee.mds.allOTUs$species[vabundOTUsnumbers,1], rotsee.mds.allOTUs$species[vabundOTUsnumbers,2],
     colnames(rotsee.vegan.otu.table)[vabundOTUsnumbers], cex=0.7, col="darkgreen")
#text(rotsee.mds.allOTUs$species[vabundOTUsnumbers,1]*1.1, rotsee.mds.allOTUs$species[vabundOTUsnumbers,2]*1.1,
#     colnames(rotsee.otu.table.t)[vabundOTUsnumbers], col="darkred")

## short Adonis model to compare (aka PERMANOVA)
variable.to.test="Oxidation_Zone"
env.glmdata.noNAs=na.exclude(env.glmdata[,variable.to.test, drop=FALSE])# Na's for the specific variable will be removed
bac.rotsee.noNAs=t(otu_table(bac.rotsee))
bac.rotsee.noNAs= bac.rotsee.noNAs[c(rownames(env.glmdata.noNAs)),]

adonis.dist.mat<-vegdist((bac.rotsee.noNAs), method="bray")
(adonis.otus=adonis(bac.rotsee.noNAs ~ Oxidation_Zone, data=env.glmdata.noNAs, permutations=999))
coefficients(adonis.otus)

########################################################################################################################
########################################################################################################################
## CCA of all OTUs  (might be better than RDA as possible ecological niches of MOBs within the water column are not linearly distributed)
########################################################################################################################
########################################################################################################################

    # data standardization and transformations
    
    ClusterHel           <-decostand(rotsee.vegan.otu.table, "hellinger")     # square root of equally standardized row sums
    ClusterNorm          <-decostand(rotsee.vegan.otu.table, method="total")  # standardized to equal row sum

    # minus vector for september 2015
   # minusvec=c(1:3,16:18,23,31:40)
   
    
    
    # (PhysicoChemicalLog=(rotsee.env.table.log[,envvector[-c(minusvec)]]))
    # PhysicoChemicalLogsub=na.omit(rotsee.env.table.log[,envvector[-c(minusvec)]])
    # PhysicoChemicalRaw=(rotsee.env.table[,envvector[-c(minusvec)]])
    # PhysicoChemicalRawsub=na.omit(rotsee.env.table[,envvector[-c(minusvec)]])
    # ClusterHelsub <- ClusterHel[rownames(PhysicoChemicalLogsub),]
    
   # (PhysicoChemicalLog=(rotsee.env.table.log[,-c(1:16,19,37,48,49)]))
   #  PhysicoChemicalLogsub=na.omit(rotsee.env.table.log[,-c(1:16,19,37,48,49)])
   #  PhysicoChemicalRaw=(rotsee.env.table[,-c(1:16,19,37,48,49)])
   #  PhysicoChemicalRawsub=na.omit(rotsee.env.table[,-c(1:16,19,37,48,49)])
   #  ClusterHelsub <- ClusterHel[rownames(PhysicoChemicalLogsub),]
    
    
    (PhysicoChemicalLog=(rotsee.env.table.log[,-c(1:13, 15:18, 48:67)]))
    PhysicoChemicalLogsub=na.omit(rotsee.env.table.log[,-c(1:13, 15:18, 48:67)])
    PhysicoChemicalRaw=(rotsee.env.table[,-c(1:13, 15:18, 48:67)])
    PhysicoChemicalRawsub=na.omit(rotsee.env.table[,-c(1:13, 15:18, 48:67)])
    ClusterHelsub <- ClusterHel[rownames(PhysicoChemicalLogsub),]
    
    
    
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
    ## remove elements that with VIFs larger than 5
    (col.env.vif<- vif_func(in_frame=PhysicoChemicalLogsub,thresh=20,trace=T))
    PhysicoChemicalLogsub = PhysicoChemicalLogsub[,c(col.env.vif)]
    
    
    remove=c("Mn_DGT", "Mn_diss_nM", "T")
    (PhysicoChemicalLogsub = PhysicoChemicalLogsub[ , !names(PhysicoChemicalLogsub) %in% remove ])
    
    
    
    ###################################################################################
    ## CCA analysis model assessments (forward, backward and marginal selection)
    ###################################################################################
    ## the selection is run on a sample subset lacking samples having missing values in the environmental matrix
    
    
    
    (Cluster.cca<-cca(ClusterHelsub ~ ., PhysicoChemicalLogsub, na.action=na.omit))   # we use the log transformed env variables
    
    (Cluster.cca.forward <- ordistep (cca(ClusterHelsub~1, PhysicoChemicalLogsub, na.action=na.omit), # ordiR2step must inherit from rda...cca not working
                                        scope= formula(Cluster.cca),
                                        direction="forward",pstep=1000,nperm=999))
    
    (Cluster.cca.backward <- ordistep (cca(ClusterHelsub~., PhysicoChemicalLogsub, na.action=na.omit),
                                       direction="backward",pstep=1000,nperm=999))
    
    (Cluster.cca.both <- ordistep (cca(ClusterHelsub~1, PhysicoChemicalLogsub, na.action=na.omit),
                                     scope= formula(Cluster.cca),
                                     direction="both",pstep=1000,nperm=999))
    
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
    
    # (FORMULA <-as.formula(paste("ClusterHel ~ ", paste  (rownames((constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=1))])), # here take again complete data set
    #                                                      collapse="+"),sep = "")))
   
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
    
    sort(vif.cca(Cluster.cca.forwsel))   # remove colinear variables >20
    
    
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
    
    # (FORMULA    <-as.formula(paste("ClusterHel ~ ", mgsub(c("\\+Mn_tot_nM", "\\+NO3._uM", "\\+Depth_m",  "\\+DIC_mM", "\\+T_C" ),  c("","", "", "", ""),
    #                                                       paste ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=1))][]), # here take again complete data set, change [] to -1 if first parameter should be replaced
    #                                                                     collapse="+")), sep = "")))
    
    #FORMULA= as.formula(ClusterHel ~ CH4_uM + SO42._uM + DIC_mM + Light_uE..m2s)    #use this one if you want to make a custom selection of variables
    
    ## rerun cca with new formula
    (Cluster.cca.forwsel <- cca(FORMULA, PhysicoChemicalLog, na.action=na.exclude))
    sort(vif.cca(Cluster.cca.forwsel))
    
    ##################----------------------------------------------------------------------------------------------
    (Cluster.cca.best.forwsel <- vegan::cca(FORMULA, PhysicoChemicalLog, na.action=na.exclude))   # this is the final model
    ##################-----------------------------------------------------------------------------------------------

    # (Cluster.cca.best.forwsel <- cca(ClusterHel ~ Cu_DGT_nM + Condition(Depth_m), PhysicoChemicalLog, na.action=na.exclude))   # alternatively make model with specific variables
    
    
    summary(Cluster.cca.best.forwsel)
    
    anova(Cluster.cca.best.forwsel) # check if total model is singificant
    

    # check for significance of environmental terms (all na data has to be removed!)
    #(FORMULAsub    <-as.formula(paste("ClusterHelsub ~ ", mgsub(c("\\+Fe_tot_nM", "\\+Depth_DGT_m", "\\+Light_uE..m2s"  ),
    #                                                            c("","",""), paste (rownames (constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=3))]), # here take again complete data set
    #                                                                              collapse="+")),sep = "")))
    
    #(FORMULAsub    <-as.formula(paste("ClusterHelsub ~ ", mgsub(c("\\+pmoA_.", "\\+Depth_DGT_m", "\\+Depth_m"  ), c("","",""),
    #                                                            paste ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=2))][-1]), # here take again complete data set
    #                                                                         collapse="+")),sep = "")))
    
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
    
    
    
    # fit the species in the rda model
    #(speciesfit   <- envfit(Cluster.cca.best.forwsel, rotsee.vegan.otu.table))
    
    ## and some additional diagnostics
    ##-------------------------------------------
    # 
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
              statistic = c("distance"), proportional = FALSE) # "explained"         
    
    
    ##---------------------------------------------------------
    # run the percentage contribution of constrains function 
    ##---------------------------------------------------------
    
    
    ################################################################################################
    #
    #  PerContCon :  Freimannsche Function to assess the importance (% variance explained)
    #                 of the environmental parameters in direction of the different canonical axes
    #                 in constraint multivariate models
    #
    ################################################################################################
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
    
    PerContCon(Cluster.cca.best.forwsel) # here you can see which parameters have higher influence on which canonical axes. 
    # one can plot the specific axes according to structures one whishes to inverstigate
    
    ## check R2 of apriori fitted constraints
    
    (envfit.apriori  <- envfit(Cluster.cca.best.forwsel, PhysicoChemicalRawsub, na.rm=TRUE, choices=c(1,2))) # you can check with the above output
    ## which dimensions might be good for a specific model parameter for
    ## fitting on and plotting of the CCA constraints
    
    ##########################################################################################################
    # Bi-plot of the CCA model         
    ##########################################################################################################
    my.color1<-c("cyan", "cyan4", "cyan3")
   
    
    #pdf("CCA pmoA September 2014.pdf",useDingbats=FALSE)
    
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
    
    
    
    ## Optionally: plot the rest of the environmental variables if you wish to....
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
    text(env.arrows.x*curr.plot.multiplier*0.72, env.arrows.y*curr.plot.multiplier*0.72, cex=0.6,
         labels=env.names, col="lightgrey")    
    
   
    
    ## plot legend and axis with explained variance
    
    
    title(xlab=as.expression(paste("CCA",xx," (",round(c(100/(Cluster.cca.best.forwsel$tot.chi)*(Cluster.cca.best.forwsel$CCA$eig[xx])),2),"%)")),
          ylab=as.expression(paste("CCA",yy," (",round(c(100/(Cluster.cca.best.forwsel$tot.chi)*(Cluster.cca.best.forwsel$CCA$eig[yy])),2),"%)")))
    
    
    dev.off()
    
    
    
    
    ## you might want to check the linearity of important variables in the cca
    
    gam.PhysicoChemRaw<-ordisurf(as.formula(paste("Cluster.cca.best.forwsel ~ ",paste(colnames(PhysicoChemicalRaw[parame<-"Chla_ug.L"])))), PhysicoChemicalRaw,
                                 choices=c(1,2), knots=10,                                         
                                 family="gaussian", isotropic= TRUE,
                                 col="darkgrey",scaling=3,
                                 add = TRUE, display = "sites",
                                 nlevels = 20, labcex = 0.8,
                                 bubble = TRUE, cex = 1, select = TRUE, method = "REML",
                                 gamma = 1, plot = TRUE)
    (summary.gam.PhysicoChemRaw <- summary(gam.PhysicoChemRaw) )
    
    legend(1,1.5, bty="n", legend=as.expression(paste("Explained deviance of\n",parame,":\n",   # change according to what you choose
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
    
    
    ## mark location of species you're looking for 
    ## phylogenetic search
    (bacterianumber <- which(rownames(ordiplot$species)==c("Methylobacter")))  
    ## OTU numbers search
    (bacterianumber <- which(colnames(rotsee.otu.table.t)==("ZOTU239")))  
    ## mark location
    points(ordiplot$species[bacterianumber,1], ordiplot$species[bacterianumber,2],
           labels=rownames(ordiplot$species)[bacterianumber], cex=2, pch=25,col="black", bg=alpha("yellow", alpha=0.4))
    ## annotate them by phylogenetic or OTU annotation
    text(ordiplot$species[bacterianumber,1], ordiplot$species[bacterianumber,2], 
         labels=rownames(ordiplot$species)[bacterianumber], cex=0.8)
    text(ordiplot$species[bacterianumber,1], ordiplot$species[bacterianumber,2],
         labels=rownames(rotsee.otu.table)[bacterianumber], cex=0.8)
    
    ########################################################################################################################
    ###########################################################################################################################
    ## distance-based RDA with all OTUs (can also be done with unifrac distances)
    ##########################################################################################################################
    ########################################################################################################################
    
    # data standardization and transformations
    ClusterHel           <-decostand(rotsee.vegan.otu.table, "hellinger")     # square root of equally standardized row sums
    ClusterNorm          <-decostand(rotsee.vegan.otu.table, method="total")  # standardized to equal row sum
    
    ## minus vector for september 2015
    # minusvec=c(1:3,16:18,23,31:40)
    # 
    # (PhysicoChemicalLog=(rotsee.env.table.log[,envvector[-c(minusvec)]]))
    # PhysicoChemicalLogsub=na.omit(rotsee.env.table.log[,envvector[-c(minusvec)]])
    # PhysicoChemicalRaw=(rotsee.env.table[,envvector[-c(minusvec)]])
    # PhysicoChemicalRawsub=na.omit(rotsee.env.table[,envvector[-c(minusvec)]])
    # ClusterHelsub <- ClusterHel[rownames(PhysicoChemicalLogsub),]
    
    ## use the same PhysicoChemicalLogsub as in the CCA analysis before
    
    
    ###################################################################################
    ## model assessments (forward, backward and marginal selection)
    ###################################################################################
    ## the selection is run on a sample subset lacking samples having missing values in the environmental matrix
    
    (Cluster.capscale<-capscale(ClusterHelsub~ ., PhysicoChemicalLogsub , distance="bray", metaMDSdist=FALSE, comm=ClusterHel, sqrt.dist = TRUE, na.action=na.omit))   # we use the log transformed env variables
    
    
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
    
    # (FORMULA <-as.formula(paste("ClusterHel ~ ", paste  (rownames((constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=3))])), # here take again complete data set
    #                                                      collapse="+"),sep = "")))
    
    (FORMULA <-as.formula(paste("ClusterHel ~ ", paste ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=2))]),
                                                        collapse="+"),sep = "")))
    
    ###########################################################################################################
    #
    # Run selected CCA model and tests and apriori fit physico-chemical variables not included in
    # the model
    #
    ###########################################################################################################
    (Cluster.capscale.forwsel <- capscale(FORMULA, PhysicoChemicalLog, distance="bray", metaMDSdist=FALSE, comm=ClusterHel, sqrt.dist = TRUE, na.action=na.exclude))
    summary(Cluster.capscale.forwsel)
    
    # test the variation inflation factors and significance of models, axes and constraints
    
    sort(vif.cca(Cluster.capscale.forwsel))   # remove variables from the data set that are highly colinear 
    
    
    ## we can introduce a condition to the model to get rid of the depth-correlated varivables
    # (FORMULA    <-as.formula(paste("ClusterHel ~ ", mgsub(c("\\+Light_uE..m2s", "\\+Depth_DGT_m", "\\+Depth_m", "\\+Fe_tot_nM",  "\\+Fe_diss_nM" ),
    #                                                       c("","","","",""), paste  (rownames(constrains.freq.table$sortsumsel[c(which(constrains.freq.table[,1]>=2))]),     
    #                                                                      collapse="+")),"+", paste("Cu_DGT_nM", collapse ="+") ,
    #                                                                                      "+ Condition(", paste("Depth_m", collapse ="+"),")" ,sep = "") ))     
    
    #(FORMULA    <-as.formula(paste("ClusterHel ~ ", mgsub(c("\\+Light_uE..m2s"), c(""),
    #                                                      paste ((constrains.freq.table$summarized.selection[c(which(constrains.freq.table[,2]>=2))]), 
    #                                                             collapse="+")),sep = "")))
    
    
    (Cluster.capscale.best.forwsel <- capscale(FORMULA, PhysicoChemicalLog, distance="bray",comm=ClusterHel, metaMDSdist=FALSE,  sqrt.dist = TRUE, na.action=na.exclude))   # this is the final model
    summary(Cluster.capscale.best.forwsel)
    sort(vif.cca(Cluster.capscale.best.forwsel))
    
    anova(Cluster.capscale.best.forwsel) # check if total model is singificant
    
    anova(Cluster.capscale.best.forwsel, by="axis", permutations = how(nperm=999),  # check for significant axes
          model = c("direct"),
          strata = NULL,
          cutoff = 1)
    
    

    
    
    # fit the non-used environmental variables in the RDA model
    
    (best.model.env.selection <-c(rownames(Cluster.capscale.best.forwsel$CCA$biplot)))
    (rest.model.env.selection <-PhysicoChemicalRaw[! c(colnames(PhysicoChemicalRaw)) %in% best.model.env.selection])
    
    (envfit.rest  <- envfit(Cluster.capscale.best.forwsel, rest.model.env.selection[,11:length(rest.model.env.selection)], na.rm=TRUE))
    
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
    
    #goodness(Cluster.capscale.best.forwsel, display = c("sites"), 
    #        model = c("CCA", "CA"), statistic = c("explained", "distance"),
    #        summarize = TRUE, addprevious = FALSE)
    
    #inertcomp(Cluster.capscale.best.forwsel, display = c("sites"),
    #          statistic = c("distance"), proportional = FALSE)
    
    # run the percentage contribution of constrains function (defined at at the end of the script)
    
    PerContCon(Cluster.capscale.best.forwsel)          # you can plot the axes which show most explained variance by the variable of interest
    
    # check R2 of apriori fitted constraints
    #(envfit.apriori  <- envfit(Cluster.capscale.best.forwsel, PhysicoChemicalRawsub, na.rm=TRUE, choices=c(1,2)))  # choose the axes according to the most explained fitting axes (PerContCon output)
    # to check specific variable of interest (i.e. Cu_DGT axes 2,3 r2 is high
    # whereas axes 1,2 r2 is low...)
    
    ##########################################################################################################
    # Bi-plot of the CCA model          (all following plots will be directly saved as PDF's)
    ##########################################################################################################
    my.color1<-c("cyan", "cyan4", "cyan3")
   
    
    #pdf("RDA of Clusters with Env as constraint.pdf",useDingbats=FALSE)
    
    # Basic Plot
    
    ordiplot<-plot(Cluster.capscale.best.forwsel,type="none",choices = c(xx<-1, yy<-2),scaling=3,main="", xlab="",ylab="")     # here you can change the axes to be ploted       

    points(ordiplot,"species", pch=21,cex=normalized.species.size.vec,
           col=col.genus[genus.df$id], bg=alpha(col.genus[genus.df$id], 0.99), lwd=0)
    legend("topright",c(unique(rotsee.tax.table[,"Genus"])),pch=c(21),col=alpha(unique(col.genus[genus.df$id]), 0.8), pt.bg=alpha(unique(col.genus[genus.df$id]),0.5), bty="n")
    points(ordiplot,"sites", cex=3, alpha= 2, pch=as.numeric(as.factor(rotsee.env.table$Campaign))+20, col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
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
    
    # plot legend and axis with explained variance
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
    ###########################################################################################################################
    # Pearson Correlation of MOBs with physico-chemical variables
    ########################################################################################################################
    ###########################################################################################################################
    # -c(1:10,12)
    ClusterHelPearson   <-decostand(rotsee.vegan.otu.table, "hellinger")     # square root of equally standardized row sums
    (colnames(ClusterHelPearson) <- paste(rownames(rotsee.otu.table),"_",colnames(ClusterHelPearson), sep=""))
    (colnames(ClusterHelPearson) <- paste(colnames(ClusterHelPearson),"_",rownames(rotsee.otu.table), sep=""))
    
    #PhysicoChemicalLog <- cbind(rotsee.env.table[,c(37),drop=FALSE],log1p(rotsee.env.table[,c(14,18:28,31,33,34,36,38:40,46)]))  # log transformed env data
    #PhysicoChemicalLogsub<-na.omit(PhysicoChemicalLog)
    #PhysicoChemicalsub<-rotsee.env.table[,c(14,18:28,31,33,34,36,37,38:40,46)]
    
    rotsee.env.table
    rotsee.env.table.log
    
    # minus vector for september 2015
   # minusvec=c(1:3,16:18,23,31:40)

    
    # c(1:5,8,25,26,27,28,29,30,31,32,33,35,36)
    
   head ( PhysicoChemicalsub<-(rotsee.env.table[-c(1:13, 15:18, 40:59)]))
   head ( PhysicoChemicalLogsub<-(rotsee.env.table.log[-c(1:13, 15:18, 40:59)]))
    
    rcorrMatrix<-as.matrix(cbind(ClusterHelPearson  , PhysicoChemicalsub))  
    
    sink(file = "PearsonCorrPhyChemClus.txt", append = FALSE, type = c("output"),
         split = FALSE)
    (CorrTable<-rcorr.adjust(rcorrMatrix[,], type="pearson"))
    sink()

    
    
    
    
   # ClusterHelPearson=ClusterHelPearson[,-c(2,16,18)]
    matrixPearson3<-CorrTable$R$r[c(1:ncol(ClusterHelPearson  )),c((ncol(ClusterHelPearson  )+1):nrow(CorrTable$R$r))]
    (matrixPearson3=na.omit(matrixPearson3))
    pvalCellNote<-(CorrTable$R$P[c(1:ncol(ClusterHelPearson  )),c((ncol(ClusterHelPearson  )+1):nrow(CorrTable$R$r))])
    ( pvalCellNote=na.omit( pvalCellNote))
    pvalCellNoteStars<- cut(pvalCellNote,  breaks=c(-Inf, 0.001, 0.01, 0.05),  label=c("***", " ** ", "  *  ")) 
    pvalCellNote<-matrix((pvalCellNoteStars), nrow=nrow(matrixPearson3), ncol=ncol(matrixPearson3))
    
    
    
    genus.numbers <- length(unique(rotsee.tax.table[,"Genus"]))
    #col.genus<-(colorRampPalette(brewer.pal(genus.numbers,"Reds"))(genus.numbers))  # choose this or the next color
    col.genus.orig<-brewer.pal(genus.numbers, "PuOr")
    genus.df.pearson<-transform(rotsee.tax.table[-c(as.numeric(attributes(matrixPearson3)$na.action)),], id=match(Genus, unique(Genus))) # or if non OTUs are missing take next line
    genus.df.pearson<-genus.df
    test=unique(genus.df.pearson$Genus)
    match(test, sort(unique(genus.df.pearson[,"Genus"]))) # copy the vector to the col.genus.orig
    col.genus.pearson<-col.genus.orig[c( 2, 6, 3, 1, 5, 4)]
    
    
   #pdf("heatmap pmoA redblue.pdf", height=10, width=20)
   envlinkttree <- heatmap.2(t(matrixPearson3),
              Rowv=TRUE,
              Colv=TRUE,
              dendrogram= c("column"),
              #distfun = NULL,
              key=TRUE,
              keysize=1, 
              trace="none",
              density.info=c("none"),
              margins=c(10, 8),
              #col=colorRampPalette(brewer.pal(9,"Greens"))(100),
              col=colorRampPalette(c("blue", "white", "red"))(100),
              cellnote=t(pvalCellNote),
              notecex=1.0,
              notecol="black",
              na.color=par("bg"),
              scale="none", 
              symkey=TRUE,
              cexRow=1,
              cexCol=0.8,
              srtRow=0,
              srtCol=90,
              ColSideColors=c(col.genus.pearson[genus.df.pearson$id]),
              key.title="Pearson",
              key.xlab=NA
              #col=redgreen(75),
    )
    
   
    
    legend("bottomleft",legend=sort(unique(genus.df[,"Genus"])),pch=c(21),col=alpha(unique(col.genus.orig[genus.df$id]), 1),
           pt.bg=alpha(unique(col.genus.orig[genus.df$id]),1), bty="n")
    
    dev.off()
    
    
    ## extract cluster tree infromation for tanglegram / cophyloplot
    library("ape")
    library("ctc")
    library(phangorn)
    library(phytools)
    
    hc <- as.hclust( envlinkttree$colDendrogram , reorder=TRUE)
    hc.phylo<- as.phylo(hc)
    cPhylo=read.tree(file="maximum Parsimony")   
    #cPhylo=read.tree(file="testtree")
    (association<- cbind(hc.phylo$tip.label, hc.phylo$tip.label))
    # alter edge length if displayed in plot
    hc.phylo$edge.length <- 80 * hc.phylo$edge.length
    cPhylo$edge.length
    
    #pdf("cophyloplot.pdf", height=10, width=20)
    cophyloplot(hc.phylo, cPhylo, assoc = association, rotate=TRUE,
                space=800, length.line = 0, gap=200,  
                use.edge.length = FALSE)
    dev.off()
    
    cor_cophenetic(hc.phylo, cPhylo, method_coef = "pearson")
    
    #### this does not work anymore in R 4.0.0 .... use an older version instead
    # devtools::install_github("rstudio/d3heatmap")
    # library( d3heatmap)
    # # heatmap with zoom function
    # d3heatmap(t(matrixPearson3),
    #           Rowv=TRUE,
    #           Colv=TRUE,
    #           dendrogram= c("column"),
    #           #distfun = NULL,
    #           key=TRUE,
    #           keysize=1, 
    #           trace="none",
    #           density.info=c("none"),
    #           margins=c(10, 8),
    #           col=colorRampPalette(brewer.pal(9,"Blues"))(100),
    #           cellnote=t(pvalCellNote),
    #           notecex=1.0,
    #           notecol="black",
    #           na.color=par("bg"),
    #           scale="none", 
    #           symkey=TRUE,
    #           cexRow=1,
    #           cexCol=1,
    #           srtRow=35,
    #           srtCol=35,
    #           #ColSideColors=c(col.genus.mobs[genus.df.mobs$id]),
    #           key.title="Pearson",
    #           key.xlab=NA
    #           #col=redgreen(75),
    # )
    
    
    
##################################################################################################################################    
## correlation plots of MOBs to env split by campaigns
    
    ClusterHelPearson   <-decostand(rotsee.vegan.otu.table, "hellinger")
    (colnames(ClusterHelPearson) <- paste(mgsub(c("\\.", "/","-"),c("","_","_"),(colnames(ClusterHelPearson))),"_",rownames(rotsee.otu.table), sep="")) # replace characters that  are not parsed propperly
    
 (  rcorrMatrixCorr<-as.matrix(cbind( ClusterHelPearson , PhysicoChemicalLog)))  # log transformed env data
 (  rcorrMatrixCorr<-as.matrix(cbind(ClusterHelPearson  , PhysicoChemicalsub)))  # raw env data
  ( rcorrMatrixCorr<- cbind( rcorrMatrixCorr, rotsee.env.table[,c(1:13, 15,16)]))
    
    # or take raw count data otus
    ClusterHelPearson   <-rotsee.vegan.otu.table
    (colnames(ClusterHelPearson) <- paste(mgsub(c("\\.", "/","-"),c("","_","_"),(colnames(ClusterHelPearson))),"_",rownames(rotsee.otu.table), sep="")) # replace characters that  are not parsed propperly
    
    (  rcorrMatrixCorr<-as.matrix(cbind( ClusterHelPearson , PhysicoChemicalLog)))  # log transformed env data
    (  rcorrMatrixCorr<-as.matrix(cbind(ClusterHelPearson  , PhysicoChemicalsub)))  # raw env data
    ( rcorrMatrixCorr<- cbind( rcorrMatrixCorr, rotsee.env.table[,c(1:13, 15,16)]))
    
   # split by campaign 
   pdf(file="Mob correlations split by campaign.pdf")

   for (i in c(1:ncol(rotsee.otu.table))){
     
               CorMob<-  paste((colnames(rcorrMatrixCorr[i])))
                       for (p in c(1:ncol(PhysicoChemicalLogsub))){
                         
                            FORMULAcorr = as.formula(paste(CorMob, "~", colnames(PhysicoChemicalsub)[p],"| Campaign", sep=""))
                            scatterplot(FORMULAcorr, data=rcorrMatrixCorr)
                            }
   }
   dev.off()
   
   # all together annotated with 
   pdf(file="Mob correlations.pdf")
   for (i in c(1:ncol(rotsee.vegan.otu.table))){
     
     CorMob<-  paste(colnames(rcorrMatrixCorr[i]))
     for (p in c(1:ncol(PhysicoChemicalLog))){
       FORMULAcorr = as.formula(paste(CorMob, "~", colnames(PhysicoChemicalsub)[p],sep=""))
       scatterplot(FORMULAcorr, data=rcorrMatrixCorr)
     }
   }
   dev.off()
   


########################################################################################################################
###########################################################################################################################
## Network Analysis
########################################################################################################################
###########################################################################################################################
ig = make_network(bac.rotsee, type="species", max.dist = 0.5
                  , distance="jaccard")
plot_network(ig, bac.rotsee, color = "Depth_m", shape = "Oxidation_Zone", line_weight = 0.2, 
                           label = "value")


gp100= prune_taxa(names(sort(taxa_sums(bac.rotsee.prune.filt3), TRUE))[1:10], bac.rotsee.prune.filt3)
jg = make_network(gp100, "taxa", "jaccard", 0.7)
plot_network(jg, gp100, "taxa",   line_weight = 0.4, label = "value")

########################################################################################################################
###########################################################################################################################
## multiple linear regression with environmental variables vs. mobs
########################################################################################################################
###########################################################################################################################
#if some MOBs occure just once...remove these as it does not make sense to model these OTUs for single season (rahter state that the on spot environmental
# conditions are important for these specific oTUs)
head(rotsee.mlr.env <- rotsee.env.table.log[-c(1:13,15,16,17,18,   48:56,  57,58,   59:67)])
rotsee.mlr.env<- na.omit(rotsee.mlr.env)
(mob.mlr1 <- decostand(as.data.frame(rotsee.vegan.otu.table[,which(colSums(rotsee.vegan.otu.table)>10)]), method="hellinger") ) # here you can choose the minimal number of reads
mob.mlr1<-mob.mlr1[c(rownames(rotsee.mlr.env)),]
(rotsee.mobs.table.mlr1 <- t(rotsee.otu.table)[, which(colSums(rotsee.vegan.otu.table)>10)])
rotsee.mobs.table.mlr1<-rotsee.mobs.table.mlr1[c(rownames(rotsee.mlr.env)),]

selectvector=vector()
for (i in c(1:ncol(mob.mlr1))){
  if (sum(mob.mlr1[,i]>0)>5){selectvector<-c(selectvector, i)}   # here you can define minimal integration points (aka site occurence) have to be in the model
}
selectvector
(mob.mlr<- subset(mob.mlr1, select= selectvector))
mob.mlr<-mob.mlr[c(rownames(mob.mlr1)),]
(mob.mlr.OTU<- subset(rotsee.mobs.table.mlr1, select= selectvector))
mob.mlr.OTU<-mob.mlr.OTU[c(rownames(mob.mlr1)),]

## This is the AIC approach...this is currently discussed and it might be better to use the "leaps" approach below
##----------------------------------------------------
mylistMLR=list()
sink("Relaimpo_Env_on_MOBs_All_AIC.txt", append=TRUE, split=TRUE)
for (i in c(1:ncol(mob.mlr))){
  
  fit <- lm(mob.mlr[,i] ~ . ,data= rotsee.mlr.env, na.action=na.omit)#
  step <- stepAIC(fit, direction="both", trace=FALSE)
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

## All-subset linear regression using lm() based on AIC  / glmulti (Model selection and multimodel inference made easy)
##-----------------------------------------------------
install.packages(glmulti)
library(glmulti)
mob.mlr.glmulti<-cbind(mob.mlr.OTU, rotsee.mlr.env)
setOldClass("rma.uni")
## set java directory if needed (i.e. if error occurs below)
# Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_271') 
setMethod('getfit', 'rma.uni', function(object, ...) {
  if (object$test=="z") {
    cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=100000)
  } else {
    cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=object$k-object$p)
  }
})

mylistMLR=list()
sink("Relaimpo_Env_on_MOBs_All_glmulti.txt", append=TRUE, split=TRUE)

for (i in c(1:ncol(mob.mlr))){
(FORMULAglmulti = as.formula(paste(colnames(mob.mlr.glmulti)[i]," ~ ", paste (colnames(mob.mlr.glmulti[,-c(1:ncol(mob.mlr))]), # here take again complete data set
                                                                         collapse="+"),sep = "")))

glmulti.lm.out <-  glmulti(FORMULAglmulti, data = mob.mlr.glmulti,
                            level = 1,               # No interaction considered
                            method = "g",            # Exhaustive approach: h , genetic:g
                            crit = "aic",            # AIC as criteria
                            confsetsize = 5,         # Keep 5 best models
                            plotty = F, report = T,  # No plot or interim reports
                            fitfunction = "lm",      # lm function
                            maxsize=10)      
print(glmulti.lm.out)
## Show 5 best models (Use @ instead of $ for an S4 object)
glmulti.lm.out@formulas
## Show result for the best model
summary(glmulti.lm.out@objects[[1]])
plot(glmulti.lm.out, type="s")
round(coef(glmulti.lm.out), 4)

FORMULAglmulti=as.formula(glmulti.lm.out@objects[[1]])

fit <- lm(FORMULAglmulti ,data= mob.mlr.glmulti)#
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

###############################################################
## for model selection you can also use Leaps...
##-----------------------------------------------------
mob.mlr.leaps<-cbind(mob.mlr.OTU, rotsee.mlr.env)
mylistMLR=list()
sink("Relaimpo_Env_on_MOBs_All_Leaps_2.txt", append=TRUE, split=TRUE)

for (i in c(1:ncol(mob.mlr))){
  (FORMULAleaps = as.formula(paste(colnames(mob.mlr.leaps)[i]," ~ ", paste (colnames(mob.mlr.leaps[,-c(1:ncol(mob.mlr))]), # here take again complete data set
                                                                            collapse="+"),sep = "")))
  regsubsets.out <-  regsubsets( FORMULAleaps,
                                 data = mob.mlr.leaps,
                                 really.big=T,
                                 nbest = 20,       # 1 best model for each number of predictors
                                 nvmax = 10,   
                                 force.in = NULL, force.out = NULL,
                                 method = "exhaustiv")
  regsubsets.out
  summary.out <- summary(regsubsets.out)
  as.data.frame(summary.out$outmat)
  plot(regsubsets.out, scale = "adjr2", main = "Adjusted R^2")
  
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
  sumfitposnega=summary(fit)                   # TYPE 2 and Type 3 SS Anova
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
  tmp<- list(posnegacor=sumfitposnega, significance=significant_reg, coefficients=coefficients_reg, relaimpo=relimportance)
  mylistMLR[[name]] <- tmp
}
mylistMLR

sink()

for (i in c(1:121)){ print(mylistMLR[i]$relaimpo$relaimp@lmg)}

###################
## PCA of leaps matrix
###################
## pearson correlation of leaps matrix
mlr.heatmap.data= read.csv("Leaps_MLR_table for heatmap.csv")
row.names(mlr.heatmap.data)=mlr.heatmap.data[,1]
mlr.heatmap.data=mlr.heatmap.data[,-1]
mlr.heatmap.data[is.na(mlr.heatmap.data)] <- 0
## table was put together in excel from the output of the leaps/relaimpo pipeline
HeatMLR   <-mlr.heatmap.data     # square root of equally standardized row sums
(rownames(HeatMLR) <- paste(colnames(rotsee.vegan.otu.table),"_",rownames(mlr.heatmap.data), sep=""))

# matrixPearson3<-CorrTable$R$r[c(1:ncol(ClusterHelPearson  )),c((ncol(ClusterHelPearson  )+1):nrow(CorrTable$R$r))]
# pvalCellNote<-(CorrTable$R$P[c(1:ncol(ClusterHelPearson  )),c((ncol(ClusterHelPearson  )+1):nrow(CorrTable$R$r))])
# pvalCellNoteStars<- cut(pvalCellNote,  breaks=c(-Inf, 0.001, 0.01, 0.05),  label=c("***", " ** ", "  *  ")) 
# pvalCellNote<-matrix((pvalCellNoteStars), nrow=nrow(matrixPearson3), ncol=ncol(matrixPearson3))
HeatMLRmatrix= HeatMLR[,-c(which(colnames(HeatMLR)== "variance.explained" | colnames(HeatMLR)=="OTU.name"))]
#pdf("heatmap all seasons pmoa 2015.pdf", height=10, width=20)
heatmap.2(t(HeatMLRmatrix),
          Rowv=TRUE,
          Colv=TRUE,
          dendrogram= c("column"),
          #distfun = NULL,
          key=TRUE,
          keysize=1, 
          trace="none",
          density.info=c("none"),
          margins=c(10, 8),
          col=colorRampPalette(c("blue", "white", "red"))(100),
          #cellnote=t(pvalCellNote),
          notecex=1.0,
          notecol="black",
          na.color=par("bg"),
          scale="none", 
          symkey=TRUE,
          cexRow=1,
          cexCol=0.8,
          srtRow=35,
          srtCol=35,
          ColSideColors=c(col.genus[genus.df$id]),
          key.title="Pearson",
          key.xlab=NA
          #col=redgreen(75),
)

legend("bottomleft",legend=sort(unique(genus.df[,"Genus"])),pch=c(21),col=alpha(unique(col.genus.orig[genus.df$id]), 1),
       pt.bg=alpha(unique(col.genus.orig[genus.df$id]),1), bty="n")

dev.off()

###################
## PCA of leaps matrix
HeatPCA= HeatMLR[,-c(which(colnames(HeatMLR)== "variance.explained" | colnames(HeatMLR)=="OTU.name"))]
HeatPCA=decostand(HeatPCA+1, "hellinger")
HeatPCA= car::logit(HeatPCA, percents=FALSE)

pca.leaps= rda(HeatPCA)
summary(pca.leaps)

biplot(pca.leaps, type="n", scaling=2)
points(pca.leaps, display ="sites", cex=(HeatMLR$variance.explained)^2 *3+1 ,#cex=1.5,
       pch=21, col=col.genus[genus.df$id], bg=alpha(col.genus[genus.df$id], 0.8), lwd=0)
#ordipointlabel(pca.leaps, display ="sites", cex=0.5, col=col.genus[genus.df$id])
legend("topright",legend=sort(unique(genus.df[,"Genus"])),pch=c(19),col=alpha(col.genus.orig, 0.99), bty="n")
points(pca.leaps, display ="species", cex=1.5, alpha= 2, pch=22, col="red",
       bg=alpha("red", 0.4), scaling=2)
text(pca.leaps, display ="species", cex=0.7, col="red", scaling=2)
ordihull(pca.leaps, groups=genus.df[,"Genus"], col=unique(col.genus.orig[genus.df$id]))
text(pca.leaps$CA$u[, 1]*1.99, pca.leaps$CA$u[, 2]*1.99,
     labels=unlist(dimnames(pca.leaps$CA$u)[1]),
     col="darkblue",#col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
     cex=0.3)
legend(-1,-0.7,legend=as.expression(paste("Stress:\n",
                                          round(c(rotsee.mds.allOTUs$stress)*100,2),
                                          "%")),cex=1,bty="n")

biplot(pca.leaps, scaling=2)

########################################################################################################################
###########################################################################################################################
## Cooccurence Analysis
########################################################################################################################
###########################################################################################################################
# The program calculates the observed and expected frequencies of co-occurrence between each pair of species.
# Expected frequency is based on the distribution of each species being random and independent of the other species.
rotsee.otu.table=otu_table(bac.rotsee)
OTU_table_cooccur <- decostand(rotsee.otu.table@.Data, "pa")[,]
#OTU_table_cooccur   <-rotsee.vegan.otu.table   
(rownames(OTU_table_cooccur) <- paste(mgsub(c("\\.", "/","-"),c("","_","_"),(colnames(rotsee.vegan.otu.table))),"_",rownames(rotsee.otu.table[,]), sep=""))# replace characters that  are not parsed propperly
 
cooccur.mobs <- cooccur(OTU_table_cooccur[c(1:90),],
                              type = "spp_site",
                              thresh = TRUE,
                              spp_names = TRUE)
class(cooccur.mobs)
summary(cooccur.mobs)
prob.table(cooccur.mobs)

# EXTRACT INFORMATION FOR A FOCAL SPECIES
sppname = paste(rownames(OTU_table_cooccur )[which(rownames(rotsee.otu.table)=="ZOTU66")])  # put in the OTU number
(pair(mod = cooccur.mobs, spp = sppname))

# generic cooccur plot
plot(cooccur.mobs)


## general stuff...

################################################################################################
#
# Merging samples
# 
# merge_samples can be very useful if you would like to see what happens to an analysis
# if you remove the indivual effects between replicates or between samples from a particular
# explanatory variable. With the merge_samples function, the abundance values of merged
# samples are summed, so make sure to do any preprocessing to account for differences
# in sequencing effort before merging or you will achieve a sequencing-effort-weighted
# average (which you may want, but keep in mind).
#
# see       http://joey711.github.io/phyloseq/merge.html      for details on merging
################################################################################################

################################################################################################
#
# Merging taxa
#
# One of the sources of noise in a microbial census dataset is a fine-scaled definition
# for OTU that blurs a pattern that might be otherwise evident were we to consider a higher
# taxonomic rank. The best way to deal with this is to use the agglomeration functions tip_glom
# or tax_glom that merge similar OTUs based on a phylogenetic or taxonomic threshold, respectively.
# However, for either function to work, it must be capable of merging two or more OTUs that have
# been deemed asequivalents.
#
################################################################################################

## merges the first 300 OTUs
merge.rotsee<- merge_taxa(bac.rotsee.prune.filt3, taxa_names(bac.rotsee.prune.filt3)[1:1000])
plot_tree(merge.rotsee)

## All tips of the tree separated by a cophenetic distance smaller than h will be agglomerated into one taxa using merge_taxa
plot_tree(bac.rotsee, label.tips="taxa_names", size="abundance", title="Before tip_glom()")
plot_tree(tip_glom(bac.rotsee, h=0.2), label.tips="taxa_names", size="abundance", title="After tip_glom()") 

# This method merges species that have the same taxonomy at a certain taxaonomic rank. Its approach is analogous
# to tip_glom, but uses categorical data instead of a tree. In principal, other categorical data known for all
# taxa could also be used in place of taxonomy, but for the moment, this must be stored in the taxonomyTable of the data.
# Also, columns/ranks to the right of the rank chosen to use for agglomeration will be replaced with NA, because
# they should be meaningless following agglomeration.
colnames(tax_table(bac.rotsee.prune.filt3))
bac.rotsee.taxmerge <- tax_glom(bac.rotsee.prune.filt3, taxrank="Phylum")
## How many taxa before/after agglomeration?
ntaxa(bac.rotsee.prune.filt3)
ntaxa(bac.rotsee.taxmerge)
plot_tree(bac.rotsee.taxmerge, label.tips="Phylum", size="abundance", sizebase=, base.spacing = 0.01, color="Cu_part_nM")

plot_tree(d, nodelabf = nodeplotboot(60, 60, 3), color = , size= "abundance",
          shape = "Phylum", ladderize = "left") + coord_polar(theta = "y")

gpac <- subset_taxa(d, Class == "Alphaproteobacteria")
plot_tree(gpac, nodelabf = nodeplotboot(60, 60, 3), color = , size= "abundance",
          shape = "Phylum", ladderize = "left") + coord_polar(theta = "y")

plot_tree(gpac, color = , shape = "Family", label.tips = "Genus", 
          size = "abundance", plot.margin = 0.6)

myTaxa = names(sort(taxa_sums(d), decreasing = TRUE)[1:20])
ex1 = prune_taxa(myTaxa, d)
plot(phy_tree(ex1), show.node.label = TRUE)
plot_tree(ex1, color = , label.tips = "Phylum", ladderize = , justify = "left" , size = "Abundance")


################################################################
## End of the script
################################################################