## The original table loaded (macExpDat) can be found in dropbox or box
#Important data structurs
#macExpDat: dataframe containing all expression values for all genes in dataset (and all tissues in dataset)
#macMax: dataframe containing max expression values for all genes across tissues
#macGenes : dataframe with ens gId and expression value for all top quartile genes

#LOAD DATA
datdir <- "~/Dropbox/WillseyLab/CPPs/data"
macExpDat <- read.delim(file.path(datdir, "macExpEnsembl.txt"))
dim(macExpDat)
macExpDat <- data.frame(gId = macExpDat$ENSEMBL, macExpDat[,c(4:13)])

#CREATE DATAFRAME OF MAXES
macMax <- data.frame(gId=macExpDat$gId, max=apply(macExpDat[,-1], 1, max))
head(macMax)

#BASIC ANALYSIS: FIND QUARTILES
summary(macMax$max)

#FILTER TO ONLY TOP QUARTILE
macGenes <- filter(macMax, max >= 338.2)
dim(macGenes)

write.table(macGenes, file=file.path(datdir, "topQuartMacGenes.txt"))

############ ENRICHMENT ANALYSIS ##############
# PIPELINE1

#original list:
datdir <- "~/Dropbox/WillseyLab/CPPs/Pipeline1/data"
wrkdir <- "~/Dropbox/WillseyLab/CPPs/Pipeline1/working"
outdir <- "~/Dropbox/WillseyLab/CPPs/Pipeline1/output"

load(file.path(datdir, "mmAAFinal_buildDF1.RData"))
mmAA <- mmAAtoSave
dim(mmAA) #52659 x 13
sum(!duplicated(mmAA$gId)) #22216
gId <- unique(mmAA$gId) #22216
macgId <- unique(macGenes$gId) #5856 unique gIds in top quart
sum(unique(macExpDat$gId) %in% gId) #18028


#########Pipeline1: mac top versus mac total

dim(mmAA) #52659 x 13
sum(!duplicated(mmAA$gId)) #22216
gId <- unique(mmAA$gId)
topmacgId <- unique(macGenes$gId) #5856
sum(topmacgId %in% gId) #5770
sum(unique(macExpDat$gId) %in% gId) #18028
allmacgId <- unique(macExpDat$gId) #18421 
sum(allmacgId %in% gId) #18028
topmacgId <- topmacgId[topmacgId %in% gId]
allmacgId <- allmacgId[allmacgId %in% gId]


#Enrichment
sum(allmacgId %in% gId) #18028
sum(topmacgId %in% gId) #5770

phyper(5770-1, 5770, 18028-5770, 18028)

#####################################################
#######filtered by R content in sliding window of 15; mmAAR
load(file.path(wrkdir, "mmAAfilteredbyR.RData"))
#Enrichment
gId <- unique(mmAAR$gId) #12721
length(gId)
sum(allmacgId %in% gId) #10887
sum(topmacgId %in% gId) #3601
phyper(3601-1, 5770, 18028-5770, 10887)
######################################################

######################################################
#filtered by SRP; mmAASP
load(file.path(wrkdir, "mmAAfilteredbyR_SRP.RData"))
#Enrichment
gId <- unique(mmAASP$gId) 
length(gId) #1935
sum(allmacgId %in% gId) #1708
sum(topmacgId %in% gId) #387
phyper(387-1, 5770, 18028-5770, 1708)
######################################################


#####################################################
#######filtered by presence of at least 1 H in prediction; mmAASS
load(file.path(wrkdir, "mmAAfilteredbyR_SRP_H.RData"))
#Enrichment
gId <- unique(mmAASS$gId) 
length(gId) #1382
sum(topmacgId %in% gId) #295
sum(allmacgId %in% gId) #1214

phyper(295-1, 5770, 18028-5770, 1214)
######################################################


#####################################################
#######only those I have predictions for; mmAApred
load(file.path(wrkdir, "mmAAfilteredbyR_SRP_H_pred.RData"))
#Enrichment
gId <- unique(mmAApred$gId) 
length(gId) #1381
sum(allmacgId %in% gId) #1214
sum(topmacgId %in% gId) #295
phyper(295-1, 5770, 18028-5770, 1214)
######################################################


#####################################################
#######Filtered by arg content on helixes
load(file.path(wrkdir, "mmAAfilteredbyR_SRP_H_RonH.RData"))
#Enrichment
gId <- unique(mmAAfiltered$gId) 
length(gId) #935
sum(allmacgId %in% gId) 
sum(topmacgId %in% gId)
phyper(201-1, 5770, 18028-5770, 839)
######################################################

# PIPELINE2
#load dataframes