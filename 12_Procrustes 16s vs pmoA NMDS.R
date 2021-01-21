########################################################################################################################################
#  run firts the script for 16s rRNA and safe the respective ordination plot under a specific name
#  repeat this with the pmoA script
########################################################################################################################################
library(vegan)
library(RColorBrewer)
library(scales) 


rotsee.mds.allOTUs.pmoa   #= rotsee.mds.allOTUs # this is the mds from the pmoA script after data cleaning etc.
rotsee.mds.allOTUs.16s    #= rotsee.mds.allOTUs # this is the mds from the 16s script after data cleaning etc.

plot(rotsee.mds.allOTUs.pmoa)
plot(rotsee.mds.allOTUs.16s)

## remove some sites from 16s that have no pmoA data
procrustepmoa = rotsee.mds.allOTUs.pmoa$points

ix= which(rownames(rotsee.mds.allOTUs.16s$points)%in% c("A01", "A02" ,"A03", "A05", "D12"))
procruste16s =rotsee.mds.allOTUs.16s$points[-ix,]
  
(rotsee.env.table.procruts =rotsee.env.table[-ix,])

## run the procrustes rotation
procru= procrustes(X = procrustepmoa, Y = procruste16s, scale =TRUE, symmetric= TRUE)
summary(procru)
plot(procru)

#pdf("Procrustes pmoa 16s.pdf" ,useDingbats=FALSE)    
procruplot= plot(procru, type="none")
points(procruplot$points, cex=2.5, alpha= 2, pch=as.numeric(rotsee.env.table.procrust$Campaign)+21,
       col="red",
       bg=alpha("red", 0.5))
points(procruplot$heads, cex=2.5, alpha= 2, pch=as.numeric(rotsee.env.table.procruts$Campaign)+21,
       col="green",
       bg=alpha("green", 0.5))
arrows(procruplot$points[,1], procruplot$points[,2], procruplot$heads[,1], procruplot$heads[,2],
       col="grey", length = 0.1, angle = 30)
#text(procruplot$points, labels=rownames(procruplot$heads), cex=0.8)
text(procruplot$heads+0.01, labels=rownames(procruplot$heads), cex=0.8)

(procruprotest= protest(X = procrustepmoa, Y = procruste16s, scores = "sites", permutations = 9999))
legend(-0.1,-0.1,legend=as.expression(paste("Correlation",":\n",
                                            round(c(procruprotest$scale)*100,2),
                                            "%,   ","P=",
                                            round(procruprotest$signif  ,3))),cex=0.8,bty="n")
#image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.3,0.7), axis.args = list(cex.axis = 0.7),
#           zlim=range(-(rotsee.env.table$Depth_m)) )


plot(procru, kind=2, xaxt="n")
axis(side=1, at=c(1:70), labels=c( unlist(attributes(procru$Yrot)$dimnames[1])), las=2, cex=0.3)

residuals(procru)


## alternatively you can run a mantel test
xdis=vegdist(decostand(rotsee.vegan.otu.table.pmoA, "hellinger"))
ydis=vegdist(decostand(rotsee.vegan.otu.table.16s[-ix,], "hellinger"))

mantel(xdis, ydis, method="pearson", permutations=9999, strata = NULL,
       na.rm = FALSE, parallel = getOption("mc.cores"))
