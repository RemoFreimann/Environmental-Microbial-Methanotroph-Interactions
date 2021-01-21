library(ggplot2)
library(fileds)
library(vegan)

##-----------------------------------------------------------------------------
## nmds of 16s with colors according to the modularity of pearson correlation network as assessedd in gephy
##-----------------------------------------------------------------------------
## take data from the Spearman script
(bac.rotsee.16s.mds.gephi= rotsee.vegan.otu.table.16s)
#(colnames(bac.rotsee.16s.mds.gephi) <- paste(colnames(rotsee.vegan.otu.table.16s),"_",rownames(rotsee.otu.table.16s), sep=""))
(ix= which(rownames(bac.rotsee.16s.mds.gephi)%in% c("A01", "A02" ,"A03", "A05", "D12")))
(bac.rotsee.16s.mds.gephi<-bac.rotsee.16s.mds.gephi[-ix,])
(bac.rotsee.16s.mds.gephi)
rotsee.env.table.modul<- rotsee.env.table[-ix,]
#(ix= which(colnames(bac.rotsee.16s.mds.gephi)%in% c(MOBs)))
#(bac.rotsee.16s.mds.gephi<-bac.rotsee.16s.mds.gephi[,-ix])
#(bac.rotsee.16s.mds.gephi)

## colors for the pearson gephi network 
## calculate the modularitie classes in gephi and reimport them here
gephi.modularity.colors.data= read.csv("modularity class.csv")
head(gephi.modularity.colors.data)
(iz = which(colnames(bac.rotsee.16s.mds.gephi) %in% (gephi.modularity.colors.data$Label)))
(bac.rotsee.16s.mds.gephi.filtered <- bac.rotsee.16s.mds.gephi[,iz])
(is = which(gephi.modularity.colors.data$Label %in% colnames(bac.rotsee.16s.mds.gephi.filtered)))
gephi.modularity.colors = gephi.modularity.colors.data[is, ]
#gephi.modularity.colors$Label=factor(gephi.modularity.colors$Label) 

bac.rotsee.16s.mds.gephi.filtered<-bac.rotsee.16s.mds.gephi.filtered[ , as.character(gephi.modularity.colors$Label)]
(gephi.modularity.colorvector <- gephi.modularity.colors$modularity_class )
sort(unique(gephi.modularity.colorvector))

test=(which(gephi.modularity.colors$modularity_class==13))
colnames(bac.rotsee.16s.mds.gephi.filtered)[test]
gephi.modularity.colors$Label[test] 

## define some colors
col.depth<-(colorRampPalette(brewer.pal(5,"Blues"))(1000))

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pal=sample(col_vector, 15)
#
colors.gephi.modularity <- brewer.pal(11, "RdYlGn")
pal <- colorRampPalette(colors.gephi.modularity)
tol21rainbow=pal(15)
#
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA",
                "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744",
                "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77",
                "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455",
                "#DD7788")
# here copy hex colors from gephi
tol21rainbow=  c("#8A633D",
                 "#DE9C58",
                 "#FFD081",
                 "#FFBB3C",
                 "#FFB94B",
                 "#FAB538",
                 "#945E26",
                 "#BA8B59",
                 "#D79429",
                 "#FFC88C",
                 "#FFD261",
                 "#FFBF6D",
                 "#F1893E",
                 "#B17427",
                 "#FFD531")   # earth tone colors in Gephi


##  run NMDS
set.seed(77)
rotsee.mds.16s.gephi<-metaMDS(bac.rotsee.16s.mds.gephi.filtered, distance="bray", trace=TRUE, plot=TRUE, k=2)
rotsee.mds.16s.gephi <- with(rotsee.env.table.modul, MDSrotate(rotsee.mds.16s.gephi, Depth_m))

plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
points(rotsee.mds.16s.gephi, display ="species", cex=0.1, pch=21, col="orange", bg=alpha("orange", 0.00001), lwd=0)
points(rotsee.mds.16s.gephi, display ="sites", cex=1.5  ,#cex=1.5,
       pch=21, col="green",
       bg="black", lwd=0)
dev.off()
## make subplots of gephi modules
#pdf("modularity groups earthtones_point sites.pdf", useDingbats = FALSE)

par.zero= par()
par=par.zero
par(mfrow=c(6,4))

hist(gephi.modularity.colorvector, xlim=c(-1,21), breaks=c(-1:21))

plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
points(rotsee.mds.16s.gephi, display ="species", cex=0.5, alpha= 3, pch=21,
       col=alpha("grey", 0.4),
       bg=alpha("lightgrey", 0.2))
points(rotsee.mds.16s.gephi, display ="sites", cex=0.6, alpha= 2, pch=as.numeric(as.factor(rotsee.env.table.modul$Campaign))+21,
       col=col.depth[as.numeric(rotsee.env.table.modul$Depth_m)/max(rotsee.env.table.modul$Depth_m)*1000],
       bg=alpha(col.depth[as.numeric(rotsee.env.table.modul$Depth_m)/max(rotsee.env.table.modul$Depth_m)*1000], 0.5))
# text(rotsee.mds.16s.gephi$points[, 1]+0.05, rotsee.mds.16s.gephi$points[, 2]+0.05,
#      labels=unlist(dimnames(rotsee.mds.16s.gephi$points)[1]),
#      col="darkblue",#col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
#      cex=0.4)
legend(-1,-1,legend=as.expression(paste("Stress:\n",
                                         round(c(rotsee.mds.16s.gephi$stress)*100,2),
                                         "%")),cex=1,bty="n")
image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.3,0.7), axis.args = list(cex.axis = 0.7),
           zlim=range(-(rotsee.env.table.modul$Depth_m)) )
with(rotsee.env.table.modul, ordiellipse(rotsee.mds.16s.gephi, groups= interaction(Oxidation_Zone), display="sites",
                                   draw="polygon", col=c( "darkblue","green","lightblue","orange"), border= "lightgrey", alpha=20 , label=FALSE))



plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 0
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 0.5), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 1
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 0.5), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 2
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 1), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 3
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 1), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 4
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 0.5), lwd=0)
 plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 5
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=1.5,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 1), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 6
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=1.5,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 1), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 7
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 1), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 8
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 0.5), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 9
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 0.5), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 10
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 0.5), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 11
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 0.5), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 12
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 0.5), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 13
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 0.5), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 14
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 0.5), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 15
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 1), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 16
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 1), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 17
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 1), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 18
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 1), lwd=0)
plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 19
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 1), lwd=0)plot(rotsee.mds.16s.gephi, type="n", xlim=c(-1.7,1.6), ylim=c(-1.2, 1.6))
modchoose= 20
points(rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 1],
       rotsee.mds.16s.gephi$species[which(gephi.modularity.colorvector== modchoose) , 2],
       cex=0.7,  pch=21,
       col=alpha(tol21rainbow[modchoose+1], 1),
       bg=alpha(tol21rainbow[modchoose+1], 1), lwd=0)

dev.off()




