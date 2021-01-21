################################################################
# load the dataset and packages and quickly check for data integrity
################################################################
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(fields)
library(mvabund)

rotsee.env.table = read.csv("data/Metafile_imputed.csv" , header = TRUE)
head(rotsee.env.table)

## make the sample IDs the column names
rownames(rotsee.env.table) = rotsee.env.table[,1]
rotsee.env.table=rotsee.env.table[,-1]
names(rotsee.env.table)

###################################################################
## here you can transform your data accordingly  (i.e. check normal distribution, heteroscedasticity etc.).
## here we produce log transformed environment data within the table
## if needed.
###################################################################
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

###################################################################
## normal line plots of variable profiles you might want to check
## (choose your variable of interes)
###################################################################
p= ggplot(rotsee.env.table, aes(Depth_m,X13CH4)) + geom_point() + geom_smooth() +
   coord_flip()  + scale_x_reverse() + facet_wrap(~Campaign) +
  scale_x_reverse(breaks=c(1:16), limits=c(15,1))
p

################################################################
## make a subset of the environmental variables and remove 
## imputed variables missing complete profiles:
## deltaC13, DIC, NO3, NO2, SO4, PO4, NH4 und Chl-a. 
##
## you skip this sectionif you want to analyse the complete
## data set with the imputed profiles.
################################################################
removeImputed = c("NO2", "NO3", "X13CH4", "SO4", "Chla", "PO4", "NH4", "DIC")

rotsee.env.table = rotsee.env.table[,!names(rotsee.env.table) %in% removeImputed]
names(rotsee.env.table)

rotsee.env.table.log = rotsee.env.table.log[,!names(rotsee.env.table.log) %in% removeImputed]
names(rotsee.env.table.log)

################################################################
## pca of environmental variables to get an overview
################################################################
PhysicoChemicalLogsub = rotsee.env.table.log[,-c(1:18,40:59)]
col.zone= c( "darkblue","green","lightblue","orange")
col.depth<-(colorRampPalette(brewer.pal(5,"Blues"))(1000)) # depth scale

pca.env= rda(PhysicoChemicalLogsub)

## use line below to procduce a pdf (also uncomment dev.off() further below)
#pdf("PCA_environmental variables_not imputed.pdf",useDingbats=FALSE)

scrs <- scores(pca.env, display = c("sites", "species"), scaling = 2)
ordiplot=plot(pca.env, type="n", scaling=2, main="", xlab="",ylab="")
points(pca.env, display ="species", cex=2 ,#cex=1.5,
       pch=21, col="lightgrey", bg=alpha("lightgrey", 0.7), lwd=0)
points(pca.env, display ="sites", cex=1.5, alpha= 2, pch=as.numeric(as.factor(rotsee.env.table$Campaign))+21, 
       col=col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
       #bg=alpha(col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000], 0.5))
       bg=alpha(col.zone[as.numeric(as.factor(rotsee.env.table$Oxidation_Zone))], 0.5))

#text(scrs$sites[, 1]+0.05, scrs$sites[, 2]+0.05,
#     labels=unlist(dimnames(scrs$sites)[1]),
#     col="darkblue",#col.depth[as.numeric(rotsee.env.table$Depth_m)/max(rotsee.env.table$Depth_m)*1000],
#     cex=0.7)

image.plot(legend.only=TRUE,col = rev(col.depth), smallplot= c(.92,.94,0.3,0.7), axis.args = list(cex.axis = 0.7),
           zlim=range(-(rotsee.env.table$Depth_m)) )
with(rotsee.env.table, ordiellipse(ordiplot, groups= interaction(Oxidation_Zone), display="sites",
                                   draw="polygon", col=c( "darkblue","green","lightblue","orange"), border= "darkgrey", alpha=20 , label=FALSE))


arrows(0, 0, scrs$species[, 1]*1, scrs$species[, 2]*1, lty=1,lwd=0.5,
       length = 0.0, angle = 15,  col=alpha("lightgrey", 0.5))

text(scrs$species[, 1]+0.05, scrs$species[, 2]+0.05,
     labels=unlist(dimnames(scrs$species)[1]),
     col="darkblue", cex=0.7)
#text(scrs$species[, 1], scrs$species[, 2],
#     labels=unlist(dimnames(scrs$species)[1]), col="black", cex=0.8)

title(xlab=as.expression(paste("PC1"," (",round(c(pca.env$CA$eig[1]),2),"%)")),
      ylab=as.expression(paste("PC2 "," (",round(c(pca.env$CA$eig[2]),2),"%)")))


dev.off()

################################################################
## test for physicom-chemically difference between oxidation zones (mvabund)
################################################################

(env.glmdata<- mvabund(rotsee.env.table.log[,c(19:39,60:62)]) )

(rotsee.glm<- manyglm(env.glmdata ~ Oxidation_Zone , data=rotsee.env.table, family="negative.binomial",    
                      cor.type="I", test="LR",
                      show.coef=FALSE, show.fitted=TRUE, show.residuals=FALSE ))

plot(rotsee.glm)

(rotsee.glm.summary<-summary(rotsee.glm, resamp="residual", nBoot=1000, test="LR"))    ## note Likelihood ratio and values
(rotsee.glm.anova<- anova(rotsee.glm, p.uni="none", test="LR",                         ##anova.manyglm   -> are liminc zones different?   
                          show.time="all", nBoot=1000))        
summary(rotsee.glm.anova)

(rotsee.glm.anova.single.sp = anova(rotsee.glm, p.uni="adjusted"))                     ## check importance of single variables

################################################################
## End of the script
################################################################