################################################################
# load the dataset and packages and quickly check for data integrity
################################################################

install.packages("mixOmics")
###########################################################################################
library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center', 
                      fig.show=TRUE, fig.keep = 'all', out.width = '50%') 

library(mixOmics)
library(phyloseq)

#########################################################################################
# Produce data set for 16s and pmoA
#########################################################################################
options(max.print = 2000000000)

## ---------------------------------------------------------------------------------------
# Create pmoA 
## ---------------------------------------------------------------------------------------
otufile     <- "data/P349_ZOTU_pmoA_sintax.tab" 
mapfile     <- "data/Metafile_imputed.txt" 
treefile    <- "data/P349_ZOTU.tre"
refseqfile  <- "data/P349_ZOTU.fa"

## ---------------------------------------------------------------------------------------
# Create primary S4 data structure from the quime output files
## ---------------------------------------------------------------------------------------
bac.rotsee.pmoA <- import_qiime(otufilename = otufile, mapfilename = mapfile,
                                treefilename = treefile, refseqfilename = refseqfile)

(rotsee.otu.sort.pmoA <-otu_table(bac.rotsee.pmoA))
(rotsee.tax.sort.pmoA <- tax_table(bac.rotsee.pmoA))
(rotsee.tree.sort.pmoA <- phy_tree(bac.rotsee.pmoA))
(rotsee.env.sort.pmoA <- sample_data(bac.rotsee.pmoA))
(rotsee.refseq.sort.pmoA <- refseq(bac.rotsee.pmoA))

## sort pmoA according to neighbour tree
library(ggtree)
d=fortify(rotsee.tree.sort.pmoA)
dd = subset(d, isTip)
phylogenetic_order_vector.pmoA=dd$label[order(dd$y, decreasing=TRUE)]

(rotsee.otu.sort.pmoA <- rotsee.otu.sort.pmoA[match(phylogenetic_order_vector.pmoA,rownames(rotsee.otu.sort.pmoA)), ])
(rotsee.tax.sort.pmoA <- rotsee.tax.sort.pmoA[match(phylogenetic_order_vector.pmoA,rownames(rotsee.tax.sort.pmoA)), ])
(bac.rotsee.pmoA = merge_phyloseq(rotsee.otu.sort.pmoA, rotsee.tax.sort.pmoA, rotsee.env.sort.pmoA, rotsee.tree.sort.pmoA , rotsee.refseq.sort.pmoA))

(rotsee.otu.sort.pmoA <-otu_table(bac.rotsee.pmoA))
(rotsee.tax.sort.pmoA <- tax_table(bac.rotsee.pmoA))
(rotsee.tree.sort.pmoA <- phy_tree(bac.rotsee.pmoA))
(rotsee.env.sort.pmoA <- sample_data(bac.rotsee.pmoA))

## ---------------------------------------------------------------------------------------
## Check and remove OTUs with ZERO reads
## ---------------------------------------------------------------------------------------
bac.rotsee.prune.pmoA <- bac.rotsee.pmoA
bac.rotsee.prune.pmoA <- prune_taxa(taxa_sums(bac.rotsee.prune.pmoA) > 0, bac.rotsee.prune.pmoA)
ntaxa(bac.rotsee.prune.pmoA) 
sample_sums(bac.rotsee.prune.pmoA)

## ---------------------------------------------------------------------------------------
## Filter OTUs that appear less than x times in at least A sample counts
## ---------------------------------------------------------------------------------------
filt3 <- genefilter_sample(bac.rotsee.prune.pmoA, filterfun_sample(function(x) x > 50), A=3)   
bac.rotsee.prune.filt3.pmoA <- prune_taxa(filt3, bac.rotsee.prune.pmoA)
ntaxa(bac.rotsee.prune.filt3.pmoA)
sample_sums(bac.rotsee.prune.filt3.pmoA)

## ---------------------------------------------------------------------------------------
##   standardize samples to the median sequencing depth 
## ---------------------------------------------------------------------------------------
total = median(sample_sums(bac.rotsee.prune.filt3.pmoA))
standf = function(x, t=total) (t * (x / sum(x)))
bac.rotsee.prune.filt3.seqdepth.pmoA = transform_sample_counts(bac.rotsee.prune.filt3.pmoA, standf)
sample_sums(bac.rotsee.prune.filt3.seqdepth.pmoA)

## ---------------------------------------------------------------------------------------
##   this is the data set for pmoA used for down stream analysis
## ---------------------------------------------------------------------------------------
#bac.rotsee.pmoA <-  bac.rotsee.prune.filt3.pmoA
bac.rotsee.pmoA <-  bac.rotsee.prune.filt3.seqdepth.pmoA
ntaxa(bac.rotsee.pmoA) 

## ---------------------------------------------------------------------------------------
##   Prune for speific taxonomic ranks if you wish to (i.e. for 16s analysis)
## ---------------------------------------------------------------------------------------
# bac.rotsee.pmoA <- tax_glom(bac.rotsee.pmoA, taxrank="Genus")

## ---------------------------------------------------------------------------------------
##   produce the tables for downstream analysis
## ---------------------------------------------------------------------------------------
(rotsee.otu.table.pmoA <- otu_table(bac.rotsee.pmoA))
## use this for the key OTUs
#(rotsee.otu.table.pmoA <- rotsee.otu.table.pmoA / colSums(rotsee.otu.table.pmoA))
(rotsee.tax.table.pmoA <- tax_table(bac.rotsee.pmoA))

(rotsee.env.table <- get_variable(bac.rotsee.pmoA)) #mapfile 

# set negative values in the metafile to 0, except for 13C-CH4 isotopes
rotsee.env.table[,-which(names(rotsee.env.table)=="X13CH4")][rotsee.env.table[,-which(names(rotsee.env.table)=="X13CH4")] < 0] <- 0 
rotsee.env.table <- rotsee.env.table[,-c(1)]

## shapiro wilk test to assess normality of data (p>0.05 is normaly distributed)
names(rotsee.env.table)
rotsee.shapiro.env= rotsee.env.table[,-c(1:13, 15:18, 48:67)]

lengthvec=c()
for (i in c(1:length(rotsee.shapiro.env))){
  
        testshapiro= shapiro.test(rotsee.shapiro.env[,i])
        if (testshapiro$p.value <= 0.05){
            lengthvec[i]=print(colnames(rotsee.shapiro.env[i]))
                    }
  }
print(na.omit(lengthvec))

# produce log transformed env table
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
names(rotsee.env.table.log)
head(rotsee.env.table.log)


####################
## remove impuded variables that have no data for one or more seasons
bad.impudation = c("NH4", "PO4", "DIC", "NO2", "NO3", "X13CH4", "Chla", "SO4", "Cr_tot")
rotsee.env.table      =  rotsee.env.table[,!names(rotsee.env.table) %in% bad.impudation]
rotsee.env.table.log  =  rotsee.env.table.log[,!names(rotsee.env.table.log) %in% bad.impudation]

########################################################################################################################################
# function to produce data sets that can be used with the package
########################################################################################################################################
rotsee.otu.table.t.pmoA<-t(rotsee.otu.table.pmoA)     
rotsee.vegan.otu.table.pmoA <- rotsee.otu.table.t.pmoA

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


rotsee.vegan.otu.table.pmoA<- make_vegan_otu_table (rotsee.otu.table.pmoA, rotsee.tax.table.pmoA)

## run this if you did not tax_glom 
(colnames(rotsee.vegan.otu.table.pmoA) <- paste(colnames(rotsee.vegan.otu.table.pmoA),"_",rownames(rotsee.otu.table.pmoA), sep=""))

## ---------------------------------------------------------------------------------------
# this is the X matrix for DAIBLO
## ---------------------------------------------------------------------------------------
(pmoAdata = rotsee.vegan.otu.table.pmoA)

## ---------------------------------------------------------------------------------------
# Create 16s data set
## ---------------------------------------------------------------------------------------
otufile     <- "data/P349_OTU_ab2_utax.tab"       
mapfile     <- "data/Metafile_imputed.txt"     
treefile    <- "data/P349_OTU_ab2.tre"
refseqfile  <- "data/P349_OTU_ab2.fa"
## ---------------------------------------------------------------------------------------
# Create primary S4 data structure from the quime output files
## ---------------------------------------------------------------------------------------
bac.rotsee.16s <- import_qiime(otufilename = otufile, mapfilename = mapfile,
                               treefilename = treefile, refseqfilename = refseqfile)

##minusvecsites.16s=c("A01","A02","A03","A04","A05","D12")
#bac.rotsee.16s <- subset_samples(bac.rotsee.16s, !X.SampleID == minusvecsites.16s) 
bac.rotsee.16s <- subset_samples(bac.rotsee.16s, !X.SampleID == "A01")
bac.rotsee.16s <- subset_samples(bac.rotsee.16s, !X.SampleID == "A02")
bac.rotsee.16s <- subset_samples(bac.rotsee.16s, !X.SampleID == "A03")
bac.rotsee.16s <- subset_samples(bac.rotsee.16s, !X.SampleID == "A04")
bac.rotsee.16s <- subset_samples(bac.rotsee.16s, !X.SampleID == "A05")
bac.rotsee.16s <- subset_samples(bac.rotsee.16s, !X.SampleID == "D12") 

head(rotsee.otu.sort.16s <-otu_table(bac.rotsee.16s))
(rotsee.tax.sort.16s <- tax_table(bac.rotsee.16s))
(rotsee.tree.sort.16s <- phy_tree(bac.rotsee.16s))
(rotsee.env.sort.16s <- sample_data(bac.rotsee.16s))
(rotsee.refseq.sort.16s <- refseq(bac.rotsee.16s))

sort(rownames(rotsee.env.sort.16s))

## ---------------------------------------------------------------------------------------
## remove or filter algea befor cynobacteria selection etc
## ---------------------------------------------------------------------------------------
chlpst <- subset_taxa(bac.rotsee.16s, Class %in% "Chloroplast")
ntaxa(chlpst) # [1] 7

## Prune chloroplasts
bac.rotsee.prune <- subset_taxa(bac.rotsee.16s,!(Class %in% c("Chloroplast"))) #, "Euglenozoa", "Chlorophyta", "Streptophyta","Cryptophyta", "Haptophyceae")))

ntaxa(bac.rotsee.prune)
bac.rotsee.16s = bac.rotsee.prune

## ---------------------------------------------------------------------------------------
## remove 16sMOBs if wished
## ---------------------------------------------------------------------------------------
#keepTaxa <- setdiff(taxa_names(bac.rotsee.16s), MOBs )
#bac.rotsee.16s <- prune_taxa(keepTaxa, bac.rotsee.16s)

## ---------------------------------------------------------------------------------------
## Check and remove OTUs with ZERO reads
## ---------------------------------------------------------------------------------------
bac.rotsee.prune.16s <- bac.rotsee.16s
bac.rotsee.prune.16s <- prune_taxa(taxa_sums(bac.rotsee.prune.16s) > 0, bac.rotsee.prune.16s)
ntaxa(bac.rotsee.prune.16s) 

## ---------------------------------------------------------------------------------------
## Filter OTUs that appear less than x times in at least A sample counts
## ---------------------------------------------------------------------------------------
filt3 <- genefilter_sample(bac.rotsee.prune.16s, filterfun_sample(function(x) x >= 3), A=2)   # A=2 as sd of 0 not allowed in diablo
bac.rotsee.prune.filt3.16s <- prune_taxa(filt3, bac.rotsee.prune.16s)
ntaxa(bac.rotsee.prune.filt3.16s)
sample_sums(bac.rotsee.prune.filt3.16s)

## ---------------------------------------------------------------------------------------
##   standardize samples to the median sequencing depth 
## ---------------------------------------------------------------------------------------
total = median(sample_sums(bac.rotsee.prune.filt3.16s))
standf = function(x, t=total) (t * (x / sum(x)))
bac.rotsee.prune.filt3.seqdepth.16s = transform_sample_counts(bac.rotsee.prune.filt3.16s, standf)
sample_sums(bac.rotsee.prune.filt3.seqdepth.16s)

## additional subfiltering
(bac.rotsee.prune.filt3.seqdepth.16s = filter_taxa(bac.rotsee.prune.filt3.seqdepth.16s,
                                                   function(x) sum(x > sample_sums(bac.rotsee.prune.filt3.seqdepth.16s)/100*0.05) > (0.1*length(x)), TRUE)) # 5% of all reads and at least in 5% of the samples
ntaxa(bac.rotsee.prune.filt3.seqdepth.16s)
sample_sums(bac.rotsee.prune.filt3.seqdepth.16s)

## ---------------------------------------------------------------------------------------
##   this is the data set for 16s used for down stream analysis
## ---------------------------------------------------------------------------------------
#bac.rotsee.16s <-  bac.rotsee.prune.filt3.16s
bac.rotsee.16s <-  bac.rotsee.prune.filt3.seqdepth.16s

## ---------------------------------------------------------------------------------------
##   Smooth data on a speific taxonomic rank if you wish to
## ---------------------------------------------------------------------------------------
# bac.rotsee.16s <- tax_glom(bac.rotsee.16s, taxrank="Phylum")
# bac.rotsee.16s.class <- tax_glom(bac.rotsee.16s, taxrank="Class")
# bac.rotsee.16s.order <- tax_glom(bac.rotsee.16s, taxrank="Order")

## ---------------------------------------------------------------------------------------
#  sort for cyanobacteria
# (bac.rotsee.16s.phylum.cyanos= subset_taxa(bac.rotsee.16s,  Phylum==c("Cyanobacteria")))
# filt3 <- genefilter_sample(bac.rotsee.16s.phylum.cyanos, filterfun_sample(function(x) x >= 5), A=2)   # A=2 as sd of 0 not allowed in diablo
# bac.rotsee.16s.phylum.cyanos <- prune_taxa(filt3, bac.rotsee.16s.phylum.cyanos)
# ntaxa(bac.rotsee.16s.phylum.cyanos)
# sample_sums(bac.rotsee.16s.phylum.cyanos)

## ---------------------------------------------------------------------------------------

## ---------------------------------------------------------------------------------------
# merge specific class and order for cobalamine  (either use the tax_glom smoothed data from above or include all OTUs)
## ---------------------------------------------------------------------------------------

## the first block is for the tax_glom data set the second one for all OTUs

# # bacterioidetes 
# bac.rotsee.16s.class.cobalamine.Flavo= subset_taxa(bac.rotsee.16s.class,  Class==c("Flavobacteriia"))
# bac.rotsee.16s.class.cobalamine.Sphingo= subset_taxa(bac.rotsee.16s.class,  Class==c( "Sphingobacteriia"))
# bac.rotsee.16s.class.cobalamine.Cyto= subset_taxa(bac.rotsee.16s.class,  Class==c("Cytophagia"))

# (bac.rotsee.16s.class.cobalamine.Flavo= subset_taxa(bac.rotsee.16s,  Class==c("Flavobacteriia")))
# (bac.rotsee.16s.class.cobalamine.Sphingo= subset_taxa(bac.rotsee.16s,  Class==c( "Sphingobacteriia")))
# (bac.rotsee.16s.class.cobalamine.Cyto= subset_taxa(bac.rotsee.16s,  Class==c("Cytophagia")))

 # cyanobacteria
# bac.rotsee.16s.order.Nostocales= subset_taxa(bac.rotsee.16s.order,  Order==c("Nostocales"))
# bac.rotsee.16s.order.Oscillatoriales= subset_taxa(bac.rotsee.16s.order,  Order==c("Oscillatoriales"))
# bac.rotsee.16s.order.Chroococcales= subset_taxa(bac.rotsee.16s.order,  Order==c("Chroococcales"))

# (bac.rotsee.16s.order.Nostocales= subset_taxa(bac.rotsee.16s,  Order==c("Nostocales")))
# (bac.rotsee.16s.order.Oscillatoriales= subset_taxa(bac.rotsee.16s,  Order==c("Oscillatoriales")))
# (bac.rotsee.16s.order.Chroococcales= subset_taxa(bac.rotsee.16s,  Order==c("Chroococcales")))

#alphabacteria
# bac.rotsee.16s.order.Sphingomonadales = subset_taxa(bac.rotsee.16s.order,  Order==c("Sphingomonadales"))
# bac.rotsee.16s.order.Caulobacterales= subset_taxa(bac.rotsee.16s.order,  Order==c("Caulobacterales"))
# bac.rotsee.16s.order.Rhodobacterales= subset_taxa(bac.rotsee.16s.order,  Order==c("Rhodobacterales"))
# bac.rotsee.16s.order.Rhizobiales= subset_taxa(bac.rotsee.16s.order,  Order==c("Rhizobiales"))

# (bac.rotsee.16s.order.Sphingomonadales = subset_taxa(bac.rotsee.16s,  Order==c("Sphingomonadales")))
# (bac.rotsee.16s.order.Caulobacterales= subset_taxa(bac.rotsee.16s,  Order==c("Caulobacterales")))
# (bac.rotsee.16s.order.Rhodobacterales= subset_taxa(bac.rotsee.16s,  Order==c("Rhodobacterales")))
# (bac.rotsee.16s.order.Rhizobiales= subset_taxa(bac.rotsee.16s,  Order==c("Rhizobiales")))

#gammabacteria
# bac.rotsee.16s.order.Pseudomonadales= subset_taxa(bac.rotsee.16s.order,  Order==c("Pseudomonadales"))
# bac.rotsee.16s.order.Chromatiales= subset_taxa(bac.rotsee.16s.order,  Order==c("Chromatiales"))
# bac.rotsee.16s.order.Alteromonadales= subset_taxa(bac.rotsee.16s.order,  Order==c("Alteromonadales"))
# bac.rotsee.16s.order.Xanthomonadales= subset_taxa(bac.rotsee.16s.order,  Order==c("Xanthomonadales"))

# (bac.rotsee.16s.order.Pseudomonadales= subset_taxa(bac.rotsee.16s,  Order==c("Pseudomonadales")))
# (bac.rotsee.16s.order.Chromatiales= subset_taxa(bac.rotsee.16s,  Order==c("Chromatiales")))
# (bac.rotsee.16s.order.Alteromonadales= subset_taxa(bac.rotsee.16s,  Order==c("Alteromonadales")))
# (bac.rotsee.16s.order.Xanthomonadales= subset_taxa(bac.rotsee.16s,  Order==c("Xanthomonadales")))

## tax_glom data set
# cobal.16s.merged = merge_phyloseq(bac.rotsee.16s.class.cobalamine.Flavo,
#                                   bac.rotsee.16s.class.cobalamine.Sphingo,
#                                   bac.rotsee.16s.class.cobalamine.Cyto,
#                                   bac.rotsee.16s.order.Nostocales,
#                                   bac.rotsee.16s.order.Oscillatoriales,
#                                   bac.rotsee.16s.order.Chroococcales,
#                                   bac.rotsee.16s.order.Sphingomonadales, 
#                                   bac.rotsee.16s.order.Caulobacterales,
#                                   bac.rotsee.16s.order.Rhodobacterales,
#                                   bac.rotsee.16s.order.Rhizobiales,
#                                   bac.rotsee.16s.order.Pseudomonadales,
#                                   bac.rotsee.16s.order.Chromatiales,
#                                   bac.rotsee.16s.order.Alteromonadales,
#                                   bac.rotsee.16s.order.Xanthomonadales)

## all otu data set
# cobal.16s.merged = merge_phyloseq(otu_table(bac.rotsee.16s.class.cobalamine.Flavo),
#                                   tax_table(bac.rotsee.16s.class.cobalamine.Flavo),
# 
#                                   otu_table(bac.rotsee.16s.class.cobalamine.Sphingo),
#                                   tax_table(bac.rotsee.16s.class.cobalamine.Sphingo),
#                                   
#                                   otu_table(bac.rotsee.16s.class.cobalamine.Cyto),
#                                   tax_table(bac.rotsee.16s.class.cobalamine.Cyto),
#                                                                     
#                                   otu_table(bac.rotsee.16s.order.Nostocales),
#                                   tax_table(bac.rotsee.16s.order.Nostocales),
#                                                                     
#                                   otu_table(bac.rotsee.16s.order.Oscillatoriales),
#                                   tax_table(bac.rotsee.16s.order.Oscillatoriales),
#                                                                     
#                                   otu_table(bac.rotsee.16s.order.Chroococcales),
#                                   tax_table (bac.rotsee.16s.order.Chroococcales),
#                                                                     
#                                   otu_table(bac.rotsee.16s.order.Sphingomonadales), 
#                                   tax_table(bac.rotsee.16s.order.Sphingomonadales),
#                                                                     
#                                   otu_table(bac.rotsee.16s.order.Caulobacterales),
#                                   tax_table(bac.rotsee.16s.order.Caulobacterales),
#                                                                     
#                                   otu_table(bac.rotsee.16s.order.Rhodobacterales),
#                                   tax_table(bac.rotsee.16s.order.Rhodobacterales),
#                                                                     
#                                   otu_table(bac.rotsee.16s.order.Rhizobiales),
#                                   tax_table(bac.rotsee.16s.order.Rhizobiales),
#                                                                     
#                                   otu_table(bac.rotsee.16s.order.Chromatiales),
#                                   tax_table(bac.rotsee.16s.order.Chromatiales),
#                                                                     
#                                   otu_table(bac.rotsee.16s.order.Alteromonadales),
#                                   tax_table(bac.rotsee.16s.order.Alteromonadales),
#                                                                     
#                                   otu_table(bac.rotsee.16s.order.Xanthomonadales),
#                                   tax_table(bac.rotsee.16s.order.Xanthomonadales))


## define one of above as bac.rotsee.16s

## ---------------------------------------------------------------------------------------
##   produce the tables for downstream analysis
## ---------------------------------------------------------------------------------------

## this is for the normal tax glom at phylum level or for all otus
head(rotsee.otu.table.16s <- otu_table(bac.rotsee.16s))
(rotsee.tax.table.16s <- tax_table(bac.rotsee.16s))

## for subset of most abundant key 16s OTUs
# total = median(sample_sums(bac.rotsee.16s))
# standf = function(x, t=total) (t * (x / sum(x)))
# bac.rotsee.prune.filt3.seqdepth.16s = transform_sample_counts(bac.rotsee.16s, standf)
# sample_sums(bac.rotsee.prune.filt3.seqdepth.16s)
# (rotsee.otu.table.16s <- otu_table(bac.rotsee.prune.filt3.seqdepth.16s))
# (rotsee.tax.table.16s <- tax_table(bac.rotsee.prune.filt3.seqdepth.16s))

## this is for the cyanobacteria
# (rotsee.otu.table.16s <- otu_table(bac.rotsee.16s.phylum.cyanos))
# (rotsee.tax.table.16s <- tax_table(bac.rotsee.16s.phylum.cyanos))

# this is for the otus potentially synthetising cobalamine
# (rotsee.otu.table.16s <- otu_table(cobal.16s.merged))
# (rotsee.tax.table.16s <- tax_table(cobal.16s.merged))

########################################################################################################################################
# function to produce data sets that can be used with the packages below
########################################################################################################################################
rotsee.otu.table.t.16s<-t(rotsee.otu.table.16s)     
#rotsee.vegan.otu.table.16s <- rotsee.otu.table.t.16s
head(rotsee.vegan.otu.table.16s<- make_vegan_otu_table (rotsee.otu.table.16s, rotsee.tax.table.16s))

# run this if you did not tax_glom
(colnames(rotsee.vegan.otu.table.16s) <- paste(colnames(rotsee.vegan.otu.table.16s),"_",rownames(rotsee.otu.table.16s), sep=""))
colnames(rotsee.vegan.otu.table.16s)

## ---------------------------------------------------------------------------------------
# this is the Y matrix 
## ---------------------------------------------------------------------------------------
(bacdata = rotsee.vegan.otu.table.16s)

## ---------------------------------------------------------------------------------------
# reorder data frames to match sampling order
## ---------------------------------------------------------------------------------------
## shuffle so that sites are in correct order and remove categorical variables from the env table
pmoAdata=pmoAdata[match(rownames(bacdata), rownames(pmoAdata)), ]
names(rotsee.env.table.log)

# rund one of the lines below
#(envdata=rotsee.env.table.log[,-c(1:13, 15:18, 40:59)]) # this is with all env parameters
(envdata=rotsee.env.table.log[,-c(1:13, 14, 15:18,22, 26, 30, 40:59)]) # remove potential biologically irrelevant or trending variables

(envdata=envdata[match(rownames(bacdata), rownames(envdata)), ])

# check that the subjects are matched in the two data sets
(cbind(rownames(pmoAdata), rownames(bacdata), rownames(envdata)))

#################################################
# Run the Analysis (supervised model)
#################################################
# extract training data
data = list(pmoA = pmoAdata, 
            rRNA = bacdata, 
            env = envdata)
str(data)
# check dimension
lapply(data, dim)

# outcome
(Trdata=rotsee.env.table.log[match(rownames(bacdata), rownames(rotsee.env.table.log)), ])
Tr = Trdata$Oxidation_Zone
summary(Tr)


# block connections
design = matrix(c(0.1), ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0
design 

## ---------------------------------------------------------------------------------------
## Tuning the number of components
## ---------------------------------------------------------------------------------------
#  fit DIABLO model without variable selection
sgccda.res = mixOmics:::block.splsda(X = data, Y = Tr, ncomp = 5, design = design)

perf.diablo = mixOmics:::perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10, progresBar=TRUE, cpus=16)   # change number of cores according to your computers cpu

# check the number of components and distance mehtods giving the best accuracy
plot(perf.diablo) 
perf.diablo$choice.ncomp$WeightedVote

# choose the best 
(ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"])
#ncomp=2
## ---------------------------------------------------------------------------------------
# define keep X components 
## ---------------------------------------------------------------------------------------
lapply(data, dim)

test.keepX = list (pmoA = c(seq(10, 20, 5), seq(30, 50, 5), seq(60,120,10)),    # for all pmoA and all 16s OTUs
                   rRNA = c(seq(10, 100, 5),seq(101,301,10), seq(302,414,14)),        # seq(302,502,25)),
                   env = c(seq(4, 9, 1),seq(10, 20, 2)))

# test.keepX = list (pmoA = c(seq(10, 50, 10), seq(60,120,30)),   # for all pmoA and most important key 16s OTUs
#                    rRNA = c(seq(10, 50, 10), seq(60,550,50)),
#                    env = c(seq(10, 18, 2), seq(20,24,4)))

# test.keepX = list (pmoA = c(seq(10, 50, 10), seq(60,120,30)),   # for all pmoA and merged 16s
#                    rRNA = c(seq(10, 18, 2), seq(20,40,5)),
#                    env = c(seq(10, 18, 2), seq(20,24,4)))

# test.keepX = list (pmoA = c(seq(4,6,1)),                       # for pmoa genus merged run and merged 16s
#                    rRNA = c(seq(10, 18, 2), seq(20,50,5)),
#                    env = c(seq(10, 18, 2), seq(20,24,4))) 

# test.keepX = list (pmoA = c(seq(10, 50, 10), seq(60,120,30)),    # for cyanos all pmoAs
#                    rRNA = c(seq(6, 18, 2)),
#                    env = c(seq(10, 18, 2), seq(20,24,4)))

# test.keepX = list (pmoA = c(seq(10, 50, 10), seq(60,120,30)),    # for cobalamine
#                    rRNA = c(seq(20, 90, 10)),
#                    env = c(seq(6, 12, 2),seq(14,18,2),seq(20,24,2)))


## a testrun to quicky check if the script works 
# test.keepX = list (pmoA = c(seq(10, 20, 5)),    # for all pmoA and all 16s OTUs
#                    rRNA = c(seq(10, 100, 30)),        
#                    env = c(seq(4, 9, 1)))
# test.keepX ## no numbers should occure twice
# tune.TCGA = mixOmics:::tune.block.splsda(X = data, Y = Tr, ncomp = ncomp, 
#                                          test.keepX = test.keepX, design = design,
#                                          validation = 'Mfold', folds = 2, nrepeat = 1,
#                                          cpus = 8, dist = "mahalanobis.dist", progressBar = TRUE)   # change distance measure if needed


## run the real model
test.keepX ## no numbers should occure twice
tune.TCGA = mixOmics:::tune.block.splsda(X = data, Y = Tr, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 3,
                              cpus = 8, dist = "mahalanobis.dist", progressBar = TRUE)   # change distance measure if needed

tune.TCGA$choice.keepX
list.keepX = tune.TCGA$choice.keepX
list.keepX

## ---------------------------------------------------------------------------------------
## Final model 
## ---------------------------------------------------------------------------------------

sgccda.res = mixOmics:::block.splsda(X = data, Y = Tr, ncomp = ncomp, keepX = list.keepX, design = design)
sgccda.res$design

## The selected variables (use them also for the variation partitioning)
test1=mixOmics:::selectVar(sgccda.res, block = 'pmoA', comp = 1)$pmoA$name 
test2=mixOmics:::selectVar(sgccda.res, block = 'pmoA', comp = 2)$pmoA$name
test3=mixOmics:::selectVar(sgccda.res, block = 'pmoA', comp = 3)$pmoA$name
#test4=mixOmics:::selectVar(sgccda.res, block = 'pmoA', comp = 4)$pmoA$name
(pmoAnames=unique(c(test1, test2, test3)))

test1=mixOmics:::selectVar(sgccda.res, block = 'rRNA', comp = 1)$rRNA$name 
test2=mixOmics:::selectVar(sgccda.res, block = 'rRNA', comp = 2)$rRNA$name 
test3=mixOmics:::selectVar(sgccda.res, block = 'rRNA', comp = 3)$rRNA$name 
(baciesnames=unique(c(test1, test2, test3)))

test1=mixOmics:::selectVar(sgccda.res, block = 'env', comp = 1)$env$name 
test2=mixOmics:::selectVar(sgccda.res, block = 'env', comp = 2)$env$name 
test3=mixOmics:::selectVar(sgccda.res, block = 'env', comp = 3)$env$name 
(envnames=unique(c(test1, test2, test3)))

mixOmics:::plotDiablo(sgccda.res, ncomp = 1)
mixOmics:::plotDiablo(sgccda.res, ncomp =2)
mixOmics:::plotDiablo(sgccda.res, ncomp =3)

## ---------------------------------------------------------------------------------------
## Plots
## ---------------------------------------------------------------------------------------

# The sample plot with the plotIndiv function projects each sample into the space spanned by the components of each block.
mixOmics:::plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO', ellipse=TRUE)

# In the arrow plot below, the start of the arrow indicates the centroid between all data sets for a given sample and
# the tips of the arrows the location of that sample in each block. 
mixOmics:::plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')

#The correlation circle plot highlights the contribution of each selected variable to each component
#pdf("Correlationplot_121pmoA_55916s_nonImputed.pdf" ,useDingbats=FALSE)  
mixOmics:::plotVar(sgccda.res, var.names = TRUE, style = 'graphics', legend = TRUE, 
        pch = c(16, 17, 15), cex = c(0.65,0.65,0.65), col = c('darkorchid', 'brown1', 'darkgreen'))
dev.off()
# The circos plot represents the correlations between variables of different types, represented on the side quadrants

#pdf("Circosplot_rel_121pmoA_55916s_nonImputed.pdf" ,useDingbats=FALSE)  
mixOmics:::circosPlot(sgccda.res, cutoff = 0.6, line = TRUE, showIntraLinks=FALSE, color.Y = c("darkblue","green","lightblue","orange"),
           color.blocks= c('darkorchid', 'brown1', 'lightgreen'),
           color.cor = c("green","red"), size.labels = 1.5)

dev.off()


# relevance network
mixOmics:::network(sgccda.res, blocks = c(1,2,3),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.5)


## make gephi network
library(igraph)
my.network = mixOmics:::network(sgccda.res, blocks = c(1,2,3),
                     color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.0) # take 0 cut-off so the additional labeling is in accordance

(groupname.pmoA = dimnames(my.network$M_pmoA_rRNA)[[1]])
(pmoA.name.index = na.omit(match(groupname.pmoA, colnames(rotsee.vegan.otu.table.pmoA))))
(namevector.gephipmoA <-as.data.frame((rotsee.tax.table.pmoA[,"Genus"]))[pmoA.name.index,])

(groupname.16s = dimnames(my.network$M_rRNA_env)[[1]])
(bac16s.name.index = na.omit(match(groupname.16s, colnames(rotsee.vegan.otu.table.16s))))
(namevector.gephi16s <-as.data.frame((rotsee.tax.table.16s[,"Phylum"]))[bac16s.name.index,])

(groupname.Env = dimnames(my.network$M_rRNA_env)[[2]])
lengthEnv = length(groupname.Env)
namevector.gephiEnv=c(rep("Environment", lengthEnv))

V(my.network$gR)$Groupname = c(as.character(namevector.gephipmoA), as.character(namevector.gephi16s),as.character(namevector.gephiEnv)) 
V(my.network$gR)$GroupSize = c(rep(50, (length(namevector.gephipmoA))), rep(20, length(namevector.gephi16s)), rep(15, lengthEnv))

write.graph(my.network$gR, file = "Network_121pmoA_4146s_Imputed.gml", format = "gml")


#pdf("Loadings_121pmoA_18cyanos_24env_negrelwithheavyC13CH4.pdf" ,useDingbats=FALSE) 
mixOmics:::plotLoadings(sgccda.res, comp = 1, contrib = 'max', method = 'median')
mixOmics:::plotLoadings(sgccda.res, comp = 2, contrib = 'max', method = 'median')
mixOmics:::plotLoadings(sgccda.res, comp = 3, contrib = 'max', method = 'median')
#mixOmics:::plotLoadings(sgccda.res, comp = 4, contrib = 'max', method = 'median')
dev.off()

rowcol=c("white", "white","white","white")
mixOmics:::cimDiablo(sgccda.res,  margins = c(15, 15), color.Y=rowcol)

#pdf("CIM_rel_121pmoA_18cyanos_24env_negrelwithheavyC13CH4.pdf" ,useDingbats=FALSE) 
mixOmics:::cimDiablo(sgccda.res,  margins = c(15, 15))
dev.off()

# Performance of the model
perf.diablo = mixOmics:::perf(sgccda.res, validation = 'Mfold', M = 10, nrepeat = 10, 
                   dist = 'centroids.dist')
perf.diablo  ## lists the different outputs
# Performance with Majority vote
perf.diablo$MajorityVote.error.rate
# Performance with Weighted prediction
perf.diablo$WeightedVote.error.rate
auc.splsda = mixOmics:::auroc(sgccda.res, roc.block = "miRNA", roc.comp = 2)



################################################################
## End of the script
################################################################
