############################################################################
##
## Variation Partitioning 
##
############################################################################

################################################################
# load the dataset and packages and quickly check for data integrity
################################################################
library(vegan)
library(caret)
library(igraph)
## use the data from the DIABLO script accordingly
pmoAdata
bacdata
envdata

pmoAnames
baciesnames
envnames

##########################################################################################################
## Variation Partitioning of Env and 16s

##--------------------------------------------------------------------------------------------------------------------
# take the variables that were choosen within the DIABLO blocks
pmoAdata.diablo = pmoAdata[,pmoAnames]
bacdata.diablo = bacdata[,baciesnames]
envdata.diablo = envdata[,envnames]

##--------------------------------------------------------------------------------------------------------------------
## filter for chosen elements within the DIABLO blocks according to the cut-off used in the Gephi/Cytoscape graph (here 0.5)
cuttoffvalue = 0.4991
choose.network = mixOmics:::network(sgccda.res, blocks = c(1,2,3), cutoff = cuttoffvalue) 

(index.pmoA.16s= which(choose.network$M_pmoA_rRNA > cuttoffvalue, arr.ind=T))
(index.16s.env= which(choose.network$M_rRNA_env > cuttoffvalue, arr.ind=T))
(index.pmoA.env= which(choose.network$M_pmoA_env > cuttoffvalue, arr.ind=T))

(pmoA.index=unique(c(rownames(index.pmoA.16s), rownames(index.pmoA.env))))
(bac.index =unique( c(colnames(choose.network$M_pmoA_rRNA)[index.pmoA.16s[,2]], rownames(index.16s.env))))
(env.index =unique( c(colnames(choose.network$M_pmoA_env)[index.pmoA.env[,2]], colnames(choose.network$M_rRNA_env)[index.16s.env[,2]])))

##-------
pmoAdata.diablo =  pmoAdata.diablo[,pmoA.index] 
# (remove1 <- readClipboard())
# (pmoAdata.diablo = pmoAdata.diablo[ , !colnames(pmoAdata.diablo) %in% remove1 ])

##-------
bacdata.diablo =  bacdata.diablo[, bac.index[bac.index %in% colnames(bacdata.diablo)]]
# (remove2 <- readClipboard())
# (bacdata.diablo = bacdata.diablo[ , !colnames(bacdata.diablo) %in% remove2 ])

##-------
envdata.diablo = envdata.diablo[,env.index]
# (remove3 <- readClipboard())
# (envdata.diablo = envdata.diablo[ , !colnames(envdata.diablo) %in% remove3 ])

##--------------------------------------------------------------------------------------------------------------------
pmoAdata.diablo = decostand(pmoAdata.diablo, "hellinger") # transform pmoA

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
df1=scale(envdata.diablo)
(col.env.vif<- vif_func(in_frame=df1,thresh=5,trace=T))
reduced_Data1 = df1[,c(col.env.vif)]
## or remove elements that correlate (pearson >0.5)
#df2=cor(df1)
#hc = findCorrelation(df2, cutoff=0.5) # putt any value as a "cutoff" 
#hc = sort(hc)
#reduced_Data1 = df1[,-c(hc)]
head(reduced_Data1)  
(reduced_Data1= scale(reduced_Data1, center = TRUE, scale = TRUE))  # z-scores of data

## Bacterial data prep
## remove elements that with VIFs larger than 3
 (df3=bacdata.diablo)
 ## change . to - etc
colnames(df3) <- gsub("\\.", "-", colnames(df3))
colnames(df3) <- gsub("\\[", "", colnames(df3))
colnames(df3) <- gsub("\\]", "", colnames(df3))
colnames(df3) <- gsub("\\X.", "", colnames(df3))
colnames(df3) 

(col.16s.vif_3<- vif_func(in_frame=df3,thresh=3,trace=T))
## change . to - etc 
(col.16s.vif_3= gsub("\\.", "-", col.16s.vif_3))
(col.16s.vif_3= gsub("\\[", "", col.16s.vif_3))
(col.16s.vif_3= gsub("\\]", "", col.16s.vif_3))
(col.16s.vif_3= gsub("\\X.", "", col.16s.vif_3))
reduced_Data2 = df3[,c(col.16s.vif_3)]
## or remove elements that correlate (pearson >0.5)
# df3=bacdata.diablo
# df4=cor(df3)
# hc2 = findCorrelation(df4, cutoff=0.5) # putt any value as a "cutoff" 
# hc2 = sort(hc2)
# reduced_Data2 = df3[,-c(hc2)]

head(reduced_Data2)
colnames(reduced_Data2)=gsub("\\-", "_", colnames(reduced_Data2) )
reduced_Data2 = decostand(reduced_Data2 , "hellinger") # transform bacies

reduced_Data =as.data.frame(cbind(reduced_Data1, reduced_Data2))
colnames(reduced_Data)


capscale(pmoAdata.diablo ~ ., as.data.frame(reduced_Data), dist="bray")

xx=c(colnames(reduced_Data))
as.formula(paste("pmoAdata.diablo", paste( xx, collapse=" + "), sep=" ~ "))   # paste the result below
paste("pmoAdata.diablo", paste( xx, collapse=" + "), sep=" ~ ")
## paste the variables in the varpart and the capscale functions accordingly below

##--------------------------------------------------------------------------------------------------------------------
## run variation partitioning followed by testing of the fractions 
pmoAdata.diablo.bray= vegdist(pmoAdata.diablo, "bray")
varp <- varpart (pmoAdata.diablo.bray, ~ Turb_NTU + Chla_ug.L + NO2._uM + Zn_diss_nM,
                 
                                       ~  RFN20_OTU_445 + PSB_M_3_OTU_167 + TM7_1_OTU_100 + Fluviicola_OTU_247 + 
                                          Chitinophagaceae_OTU_240 + Fluviicola_OTU_1464 + Fluviicola_OTU_216 + 
                                          WCHB1_41_OTU_484 + Sulfuritalea_OTU_496 + Bacteroidales_OTU_655 + 
                                          Limnohabitans_OTU_349 + Gemmataceae_OTU_325 + C111_OTU_519 + 
                                          Novosphingobium_OTU_573 + Ruminococcaceae_OTU_225 + Sphingobacteriaceae_OTU_101 + 
                                          Sphingomonadales_OTU_90 + Sphingobacteriales_OTU_202 + Phycisphaerales_OTU_536,
                                          data = reduced_Data)

varp
plot(varp)

## check partitioning and if fractions are significant
# fractions [a+b+c]:
(pmoAdata.hel.all <- capscale(pmoAdata.diablo ~ Turb_NTU + Chla_ug.L + NO2._uM + Zn_diss_nM + 
                                RFN20_OTU_445 + PSB_M_3_OTU_167 + TM7_1_OTU_100 + Fluviicola_OTU_247 + 
                                Chitinophagaceae_OTU_240 + Fluviicola_OTU_1464 + Fluviicola_OTU_216 + 
                                WCHB1_41_OTU_484 + Sulfuritalea_OTU_496 + Bacteroidales_OTU_655 + 
                                Limnohabitans_OTU_349 + Gemmataceae_OTU_325 + C111_OTU_519 + 
                                Novosphingobium_OTU_573 + Ruminococcaceae_OTU_225 + Sphingobacteriaceae_OTU_101 + 
                                Sphingomonadales_OTU_90 + Sphingobacteriales_OTU_202 + Phycisphaerales_OTU_536, reduced_Data, dist="bray" ))

RsquareAdj (pmoAdata.hel.all)
anova(pmoAdata.hel.all)

# pmoAdata.hel.all.forw.0 = capscale(pmoAdata.diablo ~ 1, data=reduced_Data)
# pmoAdata.hel.all.forw.all = capscale(pmoAdata.diablo ~ ., data=reduced_Data)
# sel.os <- ordistep (pmoAdata.hel.all.forw.0, scope = formula (pmoAdata.hel.all.forw.all), direction = 'forward')
# sel.os


# fraction [a] env only:
(pmoAdata.hel.partial.env <- vegan::capscale(pmoAdata.diablo ~ Turb_NTU + Chla_ug.L + NO2._uM + Zn_diss_nM + 
                                        Condition(RFN20_OTU_445 + PSB_M_3_OTU_167 + TM7_1_OTU_100 + Fluviicola_OTU_247 + 
                                                    Chitinophagaceae_OTU_240 + Fluviicola_OTU_1464 + Fluviicola_OTU_216 + 
                                                    WCHB1_41_OTU_484 + Sulfuritalea_OTU_496 + Bacteroidales_OTU_655 + 
                                                    Limnohabitans_OTU_349 + Gemmataceae_OTU_325 + C111_OTU_519 + 
                                                    Novosphingobium_OTU_573 + Ruminococcaceae_OTU_225 + Sphingobacteriaceae_OTU_101 + 
                                                    Sphingomonadales_OTU_90 + Sphingobacteriales_OTU_202 + Phycisphaerales_OTU_536
                                                  ), reduced_Data, dist="bray" ))
RsquareAdj (pmoAdata.hel.partial.env)
anova(pmoAdata.hel.partial.env)

# fraction [c] baccies only:
(pmoAdata.hel.partial.bacies <- capscale(pmoAdata.diablo ~   RFN20_OTU_445 + PSB_M_3_OTU_167 + TM7_1_OTU_100 + Fluviicola_OTU_247 + 
                                           Chitinophagaceae_OTU_240 + Fluviicola_OTU_1464 + Fluviicola_OTU_216 + 
                                           WCHB1_41_OTU_484 + Sulfuritalea_OTU_496 + Bacteroidales_OTU_655 + 
                                           Limnohabitans_OTU_349 + Gemmataceae_OTU_325 + C111_OTU_519 + 
                                           Novosphingobium_OTU_573 + Ruminococcaceae_OTU_225 + Sphingobacteriaceae_OTU_101 + 
                                           Sphingomonadales_OTU_90 + Sphingobacteriales_OTU_202 + Phycisphaerales_OTU_536 +
                                     Condition(Turb_NTU + Chla_ug.L + NO2._uM + Zn_diss_nM ), reduced_Data, dist="bray"))
RsquareAdj (pmoAdata.hel.partial.bacies)
anova(pmoAdata.hel.partial.bacies)

# fraction [a+b]:
(pmoAdata.hel.partial.envonly <- capscale(pmoAdata.diablo ~   Turb_NTU + Chla_ug.L + NO2._uM + Zn_diss_nM , reduced_Data, dist="bray" ))
RsquareAdj (pmoAdata.hel.partial.envonly)
anova(pmoAdata.hel.partial.envonly)

# fractions [b+c]:
(pmoAdata.hel.partial.baconly <- capscale(pmoAdata.diablo ~  RFN20_OTU_445 + PSB_M_3_OTU_167 + TM7_1_OTU_100 + Fluviicola_OTU_247 + 
                                            Chitinophagaceae_OTU_240 + Fluviicola_OTU_1464 + Fluviicola_OTU_216 + 
                                            WCHB1_41_OTU_484 + Sulfuritalea_OTU_496 + Bacteroidales_OTU_655 + 
                                            Limnohabitans_OTU_349 + Gemmataceae_OTU_325 + C111_OTU_519 + 
                                            Novosphingobium_OTU_573 + Ruminococcaceae_OTU_225 + Sphingobacteriaceae_OTU_101 + 
                                            Sphingomonadales_OTU_90 + Sphingobacteriales_OTU_202 + Phycisphaerales_OTU_536  , reduced_Data, dist="bray"))
RsquareAdj (pmoAdata.hel.partial.baconly)
anova(pmoAdata.hel.partial.baconly)


## added some additional env variables central within the network modules
# or ad some additional env parameters within the specific modules to the reduced_Data1, i.e. Cu_DGT and CH4 
## this is just to play around .... not used for analysis
reduced_Data3= scale(envdata.diablo, center = TRUE, scale = TRUE)
names(envdata)   
reduced_Data4 =as.data.frame(cbind(reduced_Data3, reduced_Data2))
colnames(reduced_Data4)

xxx=c(colnames(reduced_Data4))
as.formula(paste("pmoAdata.diablo", paste(xxx, collapse=" + "), sep=" ~ "))   # paste the result below

pmoAdata.diablo.bray= vegdist(pmoAdata.diablo, "bray")
varp <- varpart (pmoAdata.diablo.bray, ~ Chla_ug.L + CH4_uM + Fe_DGT_nM + Cu_DGT_nM + Cu_diss_nM + Fe_diss_nM,
                 ~  Bacteroidales_OTU_276 + Solirubrobacterales_OTU_71 + Methylophilaceae_OTU_4 + 
                   Victivallaceae_OTU_388 + Planctomyces_OTU_423 + Chitinophagaceae_OTU_93, data = reduced_Data4)
varp
plot(varp)



################################################################
## End of the script
################################################################

