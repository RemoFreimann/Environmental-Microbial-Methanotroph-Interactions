install.packages("picante")

library(picante)
library(vegan)
library(ggplot2)
library(car)
library(multcomp)
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
##----------------------------------------------

# The function mpd will calculate the mean pairwise distance between all
# species or individuals in each community. Similarly, the mntd function
# calculates the mean nearest taxon distance, the average distance separating
# each speciesor individual in the community from its closest heterospecific relative. 

# Measures of 'standardized effect size' of phylogenetic community structure
# can be calculated for MPD and MNTD by compared observed phylogenetic relatedness
# to the pattern expected under some null model of phylogeny or community randomization.
# Standardized effect sizes describe the difference between average phylogenetic distances
# in the observed communities versus null communities generated with some
# randomization method,standardized by the standard deviation of phylogenetic
# distances in the null data
# 
# Positive SES values (mpd.obs.z > 0) and high quantiles (mpd.obs.p > 0.95) 
# indicate phylogenetic evenness, while negative SES values and low quantiles
# (mpd.obs.p < 0.05) indicate phylogenetic clustering, relative to the null model.
# MPD is generally thought to be more sensitive to tree-wide patterns of phylogenetic
# clustering and eveness, while MNTD is moresensitive to patterns of evenness and
# clustering closer to the tips of the phylogeny


#################################################################################################
##---------------------------------------------------------------------------------------------
## prepare 16s data 
##---------------------------------------------------------------------------------------------
#################################################################################################
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

#save chlorplast data, just in case we need it later
chlpst <- phyloseq::subset_taxa(bac.rotsee, Class %in% "Chloroplast")
ntaxa(chlpst) 

## Remove chloroplasts from Bacterial OTU table
bac.rotsee.prune <- phyloseq::subset_taxa(bac.rotsee.prune, !(Class %in% "Chloroplast"))
ntaxa(bac.rotsee.prune)

## Check and remove OTUs with ZERO reads
bac.rotsee.prune <- phyloseq::prune_taxa(taxa_sums(bac.rotsee.prune) > 0, bac.rotsee.prune)
ntaxa(bac.rotsee.prune) 

## ---------------------------------------------------------------------------------------
#   FILTERING (should be done according to your hypthesis, i.e. if you're looking for rare
#              species also or just want to process the main players)
## ---------------------------------------------------------------------------------------
## Filter OTUs that appear at least in three different samples
filt3 <- genefilter_sample(bac.rotsee.prune, filterfun_sample(function(x) x > 3), A=1)
bac.rotsee.prune.filt3 <- prune_taxa(filt3, bac.rotsee.prune)
ntaxa(bac.rotsee.prune.filt3) 

## ---------------------------------------------------------------------------------------
## Standardization of sample reads to the median sequencing depth applied to each sample  (readcounds, 100% and 1)
## ---------------------------------------------------------------------------------------
total = median(sample_sums(bac.rotsee.prune.filt3))
standf = function(x, t=total) (t * (x / sum(x)))
bac.rotsee.prune.filt3.seqdepth = transform_sample_counts(bac.rotsee.prune.filt3, standf)

(rotsee.otu.sort <-otu_table(bac.rotsee.prune.filt3.seqdepth))
(rotsee.tax.sort <- tax_table(bac.rotsee.prune.filt3.seqdepth))
(rotsee.tree.sort <- phy_tree(bac.rotsee.prune.filt3.seqdepth))
(rotsee.env.sort <- sample_data(bac.rotsee.prune.filt3.seqdepth))
(rotsee.refseq.sort <- refseq(bac.rotsee.prune.filt3.seqdepth))

comm.16s=rotsee.otu.sort
phy.tree.16s=rotsee.tree.sort
metadata=rotsee.env.sort

##----------------------------------------------
#######################################
## run for 16s
#######################################
comm16s = decostand(t(comm.16s), "total")
phy16s  = phy.tree.16s

# convert phylogenety to a distance matrix
phy.dist.16s <- cophenetic(phy16s)

# Standardized effect size of mean pairwise distances in communities
# When used with a phylogenetic distance matrix, equivalent
# to -1 times the Nearest Relative Index (NRI).
# ses.mpd is negative  for clustered communities and positive or evenly 
# spread communties
# calculate ses.mpd (equals -1x NRI)
comm16s.sesmpd <- ses.mpd(comm16s, phy.dist.16s, null.model= "richness", abundance.weighted = TRUE, 
                       runs = 999, iterations = 1000)
head(comm16s.sesmpd)
# compare ses.mpd between habitats
plot(comm16s.sesmpd$mpd.obs.z ~ as.factor(metadata$Oxidation_Zone)    , xlab = "Zone", ylab = "SES(MPD)")
abline(h = 0, col = "gray")
#t.test(comm.sesmpd$mpd.obs.z ~ metadata$Oxidation_Zone)

## test between zones
mpds=comm16s.sesmpd$mpd.obs.z
oxigroups=as.factor(metadata$Oxidation_Zone)
oxi_anova=data.frame(mpds, as.factor(oxigroups))
colnames(oxi_anova)[2]<-"oxigroup"
barorder=c("Oxic", "Oxycline" , "Oxidation_Zone", "Anoxic")
#Anova of groups and posthoctests
options(contrasts=c(factor="contr.sum",ordered="contr.poly"))      ## SPSS type III ... marginal 
(Anova(Tukeylm<-lm(mpds ~ oxigroups, data=(oxi_anova)),
       type=3, singular.ok=TRUE))
(TukeComp<-summary(glht(Tukeylm, linfct= mcp(oxigroups="Tukey"), p.adjust="holmes")))
mod.cld <- cld(TukeComp)

(plot.posthoc.annot <-(mod.cld$mcletters$Letters))
plot.posthoc.annot= data.frame(plot.posthoc.annot)
plot.posthoc.annot$id=rownames(plot.posthoc.annot)
plot.posthoc.annot=transform(plot.posthoc.annot, id=match(id, barorder))
plot.posthoc.annot=arrange(plot.posthoc.annot, id)

## boxplot with specific order 
#pdf("ses of mpd 16s.pdf" ,useDingbats=FALSE) 
p <- ggplot(oxi_anova, aes(x = oxigroups, y = mpds)) + 
  geom_boxplot( aes(colour=oxigroups),
                alpha =0.8,
                cex=1,
                width=0.7) +
  theme_bw() +
  scale_x_discrete(name="", limits = barorder,
                   labels=c("Oxic", "Oxycline" , "Oxidation_Zone", "Anoxic")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +
  geom_hline(yintercept = 0, lty=2) +
  scale_y_continuous(name="ses of mpd", expand=c(0,0), 
                     breaks = c(seq(-6, 1, 1)) , 
                     limits=c(-6, 1) ) +
  annotate("text", label=paste(plot.posthoc.annot[,1]), y= 0.92 , x= c(1:4), cex=4)
p
dev.off()

# plot species present in a fescue community
plot(phy, show.tip.label = FALSE, main = "Community A06")
tiplabels(tip = which(phy$tip.label %in% colnames(comm)[comm["A06", ] > 0.01]), pch = 20, col="green" )


# Standardized effect size of mean nearest taxon distances in communities
# When used with a phylogenetic distance matrix, equivalent
# to -1 times the Nearest Taxon Index (NTI).
#  negative values of ses.mntd indicate that species co-occur with 
# more closely related species than expected, and positive values 
# indicate that closely related species do not co-ccur
# calculate ses.mntd
comm16s.sesmntd <- ses.mntd(comm16s, phy.dist.16s, null.model = "richness", abundance.weighted = FALSE, 
                         runs = 999, iterations=1000)
head(comm16s.sesmntd)
# compare ses.mntd between habitats
plot(comm16s.sesmntd$mntd.obs.z ~ as.factor(metadata$Oxidation_Zone), xlab = "Habitat", ylab = "SES(MNTD)")
abline(h = 0, col = "gray")
#t.test(comm.sesmntd$mntd.obs.z ~ metadata$Oxidation_Zone)
## test between zones only between zones
mntd=comm16s.sesmntd$mntd.obs.z
oxigroups=as.factor(metadata$Oxidation_Zone)
oxi_anova=data.frame(mntd, as.factor(oxigroups))
colnames(oxi_anova)[2]<-"oxigroup"
barorder=c("Oxic", "Oxycline" , "Oxidation_Zone", "Anoxic")
#Anova of groups and posthoctests
options(contrasts=c(factor="contr.sum",ordered="contr.poly"))      ## SPSS type III ... marginal 
(Anova(Tukeylm<-lm(mntd ~ oxigroups, data=(oxi_anova)),
       type=3, singular.ok=TRUE))
(TukeComp<-summary(glht(Tukeylm, linfct= mcp(oxigroups="Tukey"), p.adjust="holmes")))
mod.cld <- cld(TukeComp)

(plot.posthoc.annot <-(mod.cld$mcletters$Letters))
plot.posthoc.annot= data.frame(plot.posthoc.annot)
plot.posthoc.annot$id=rownames(plot.posthoc.annot)
plot.posthoc.annot=transform(plot.posthoc.annot, id=match(id, barorder))
plot.posthoc.annot=arrange(plot.posthoc.annot, id)

## boxplot with specific order 
#pdf("ses of mntd 16s.pdf" ,useDingbats=FALSE) 
p <- ggplot(oxi_anova, aes(x = oxigroups, y = mntd)) + 
  geom_boxplot( aes(colour=oxigroups),
                alpha =0.8,
                cex=1,
                width=0.7) +
  theme_bw() +
  scale_x_discrete(name="", limits = barorder,
                   labels=c("Oxic", "Oxycline" , "Oxidation_Zone", "Anoxic")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +
  geom_hline(yintercept = 0, lty=2) +
  scale_y_continuous(name="ses of mntd", expand=c(0,0), 
                     breaks = c(seq(-11, 1, 1)) , 
                     limits=c(-11, 1) ) +
  annotate("text", label=paste(plot.posthoc.annot[,1]), y= 0.92 , x= c(1:4), cex=4)
p
dev.off()


# plot species present in a fescue community
plot(phy, show.tip.label = FALSE, main = "Community A01")
tiplabels(tip = which(phy$tip.label %in% colnames(comm)[comm["A01", ] > 0.01]), pch = 20, col="green" )




#################################################################################################
##---------------------------------------------------------------------------------------------
## prepare MOB data
##---------------------------------------------------------------------------------------------
#################################################################################################
otufile     <- "data/P349_ZOTU_pmoA_sintax_clean.tab"  
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

## use this if you dont want to standardize to seq depth
#bac.rotsee.prune.filt3.seqdepth <- bac.rotsee.prune.filt3


(rotsee.otu.sort <-otu_table(bac.rotsee.prune.filt3.seqdepth))
(rotsee.tax.sort <- tax_table(bac.rotsee.prune.filt3.seqdepth))
(rotsee.tree.sort <- phy_tree(bac.rotsee.prune.filt3.seqdepth))
(rotsee.env.sort <- sample_data(bac.rotsee.prune.filt3.seqdepth))
(rotsee.refseq.sort <- refseq(bac.rotsee.prune.filt3.seqdepth))

comm.pmoA=rotsee.otu.sort
phy.tree.pmoA=rotsee.tree.sort
metadatapmoA=rotsee.env.sort

##############################################################
## run for pmoA
##############################################################
comm = decostand(t(comm.pmoA), "total")
phy  = phy.tree.pmoA

phy.dist <- cophenetic(phy)

ix= which(rownames(metadatapmoA)%in% c("A01", "A02" ,"A03", "A05", "D12"))
metadatapmoA<-metadata[-ix,]

# calculate ses.mpd (equals -1x NRI)
comm.pmoA.sesmpd <- ses.mpd(comm, phy.dist, null.model= "richness", abundance.weighted = TRUE, 
                       runs = 999, iterations = 1000)
head(comm.pmoA.sesmpd)
# compare ses.mpd between habitats
plot(comm.pmoA.sesmpd$mpd.obs.z ~ as.factor(metadatapmoA$Oxidation_Zone)    , xlab = "Zone", ylab = "SES(MPD)")
abline(h = 0, col = "gray")
#t.test(comm.pmoA.sesmpd$mpd.obs.z ~ metadata$Oxidation_Zone)
## test between zones only between zones
mpds=comm.pmoA.sesmpd$mpd.obs.z
oxigroups=as.factor(metadatapmoA$Oxidation_Zone)
oxi_anova=data.frame(mpds, oxigroups)
barorder=c("Oxic", "Oxycline" , "Oxidation_Zone", "Anoxic")
#Anova of groups and posthoctests
options(contrasts=c(factor="contr.sum",ordered="contr.poly"))      ## SPSS type III ... marginal 
(Anova(Tukeylm<-lm(mpds ~ oxigroups, data=(oxi_anova)),
                          type=3, singular.ok=TRUE))
(TukeComp<-summary(glht(Tukeylm, linfct= mcp(oxigroups="Tukey"), p.adjust="holmes")))
mod.cld <- cld(TukeComp)

(plot.posthoc.annot <-(mod.cld$mcletters$Letters))
plot.posthoc.annot= data.frame(plot.posthoc.annot)
plot.posthoc.annot$id=rownames(plot.posthoc.annot)
plot.posthoc.annot=transform(plot.posthoc.annot, id=match(id, barorder))
plot.posthoc.annot=arrange(plot.posthoc.annot, id)

## boxplot with specific order 
#pdf("ses of mpd pmoA.pdf" ,useDingbats=FALSE)    
p <- ggplot(oxi_anova, aes(x = oxigroups, y = mpds)) + 
  geom_boxplot( aes(colour=oxigroups),
                alpha =0.8,
                cex=1,
                width=0.7) +
  theme_bw() +
  scale_x_discrete(name="", limits = barorder,
                   labels=c("Oxic", "Oxycline" , "Oxidation_Zone", "Anoxic")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +
  geom_hline(yintercept = 0, lty=2) +
  scale_y_continuous(name="ses of mpd", expand=c(0,0), 
                     breaks = c(seq(-5.5, 2, 1)) , 
                     limits=c(-5.5, 2) ) +
  annotate("text", label=paste(plot.posthoc.annot[,1]), y= 0.92 , x= c(1:4), cex=4)
p
dev.off()
# plot species present in a fescue community
plot(phy, show.tip.label = FALSE, main = "Community A20")
tiplabels(tip = which(phy$tip.label %in% colnames(comm)[comm["A20", ] > 0.01]), pch = 20, col="green" )


# calculate ses.mntd
comm.pmoA.sesmntd <- ses.mntd(comm, phy.dist, null.model = "richness", abundance.weighted = FALSE, 
                         runs = 999, iterations=1000)
head(comm.pmoA.sesmntd)
# compare ses.mntd between habitats
plot(comm.pmoA.sesmntd$mntd.obs.z ~ as.factor(metadatapmoA$Oxidation_Zone), xlab = "Habitat", ylab = "SES(MNTD)")
abline(h = 0, col = "gray")
##t.test(comm.pmoA.sesmntd$mntd.obs.z ~ metadata$Oxidation_Zone)
## test between zones only between zones
mntd=comm.pmoA.sesmntd$mntd.obs.z
oxigroups=as.factor(metadatapmoA$Oxidation_Zone)
oxi_anova=data.frame(mntd, oxigroups)
barorder=c("Oxic", "Oxycline" , "Oxidation_Zone", "Anoxic")
#Anova of groups and posthoctests
options(contrasts=c(factor="contr.sum",ordered="contr.poly"))      ## SPSS type III ... marginal 
(Anova(Tukeylm<-lm(mntd ~ oxigroups, data=(oxi_anova)),
                   type=3, singular.ok=TRUE))
(TukeComp<-summary(glht(Tukeylm, linfct= mcp(oxigroups="Tukey"), p.adjust="holmes")))
mod.cld <- cld(TukeComp)

(plot.posthoc.annot <-(mod.cld$mcletters$Letters))
plot.posthoc.annot= data.frame(plot.posthoc.annot)
plot.posthoc.annot$id=rownames(plot.posthoc.annot)
plot.posthoc.annot=transform(plot.posthoc.annot, id=match(id, barorder))
plot.posthoc.annot=arrange(plot.posthoc.annot, id)

## boxplot with specific order 
#pdf("ses of mntd pmoA.pdf" ,useDingbats=FALSE) 
p <- ggplot(oxi_anova, aes(x = oxigroups, y = mntd)) + 
  geom_boxplot( aes(colour=oxigroups),
                alpha =0.8,
                cex=1,
                width=0.7) +
  theme_bw() +
  scale_x_discrete(name="", limits = barorder,
                   labels=c("Oxic", "Oxycline" , "Oxidation_Zone", "Anoxic")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +
  geom_hline(yintercept = 0, lty=2) +
  scale_y_continuous(name="ses of mntd", expand=c(0,0), 
                     breaks = c(seq(-5, 1.5, 1)) , 
                     limits=c(-5, 1.8) ) +
  annotate("text", label=paste(plot.posthoc.annot[,1]), y= 0.92 , x= c(1:4), cex=4)
p
dev.off()

# plot species present in a fescue community
plot(phy, show.tip.label = FALSE, main = "Community B03")
tiplabels(tip = which(phy$tip.label %in% colnames(comm)[comm["B03", ] > 0.01]), pch = 20, col="green" )


################################################################
## End of the script
################################################################