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
#write.table(macExpDat, file=file.path(datdir, "allQuartMacGenes.txt"), sep="\t", quote=F)

#################### TAKING MAXES
#LOAD DATA
datdir <- "~/Dropbox/WillseyLab/CPPs/data"
macExpDat <- read.delim(file.path(datdir, "macExpEnsembl.txt"))
dim(macExpDat)
macExpDat <- data.frame(gId = macExpDat$ENSEMBL, macExpDat[,c(4:13)])
#write.table(macExpDat, file=file.path(datdir, "allQuartMacGenes.txt"), sep="\t", quote=F)

P1 <- aggregate(Peritoneum_untreated_rep1 ~ ENSEMBL, macExpDat, max)
P2 <- aggregate(Peritoneum_untreated_rep2 ~ ENSEMBL, macExpDat, max)
Lu1 <- aggregate(Lung_untreated_rep1 ~ ENSEMBL, macExpDat, max)
Li1 <- aggregate(Liver_untreated_rep1 ~ ENSEMBL, macExpDat, max)
S1 <- aggregate(Spleen_untreated_rep1 ~ ENSEMBL, macExpDat, max)
I1 <- aggregate(Intestine_untreated_rep1 ~ ENSEMBL, macExpDat, max)
I2 <- aggregate(Intestine_untreated_rep2 ~ ENSEMBL, macExpDat, max)
A1 <- aggregate(Adipose_untreated_rep1 ~ ENSEMBL, macExpDat, max)
A2 <- aggregate(Adipose_untreated_rep2 ~ ENSEMBL, macExpDat, max)
B1 <- aggregate(BMDM_untreated_rep1 ~ ENSEMBL, macExpDat, max)

d1 <- full_join(P1, P2)
d2 <- full_join(d1, Lu1)
d3 <- full_join(d2, Li1)
d4 <- full_join(d3, S1)
d5 <- full_join(d4, I1)
d6 <- full_join(d5, I2)
d7 <- full_join(d6, A1)
d8 <- full_join(d7, A2)
d9 <- full_join(d8, B1)
dim(d9)
mbyg <- d9

write.table(mbyg, file.path(datdir, "macExpMaxesPerTissue.txt"), sep="\t", quote=F)

head(mbyg)
mbyg.max <- apply(mbyg[,-1], 1, max)
mbyg.max <- data.frame(mbyg.max)
head(mbyg.max)
row.names(mbyg.max) <- mbyg$ENSEMBL
head(mbyg.max)
summary(mbyg.max$mbyg.max)
mbyg.max$gIds <- rownames(mbyg.max)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 100.9   125.5   155.5   926.0   530.1 42860.0 
topQ <- filter(mbyg.max, mbyg.max >= 530.1)
dim(mbyg.max)
dim(topQ)
head(topQ)

write.table(topQ, file=file.path(datdir, "topQuartMacGenes.txt"), sep="\t", quote=F)
write.table(mbyg.max, file=file.path(datdir, "allMacGeneMaxes.txt"), sep="\t", quote=F)

topQgId <- topQ$gId #are all unique Ids
allQgId <- mbyg.max$gId #also all unique Ids

############ ENRICHMENT ANALYSIS ##############
# PIPELINE1

#original list:
datdir <- "~/Dropbox/WillseyLab/CPPs/Pipeline1/data"
wrkdir <- "~/Dropbox/WillseyLab/CPPs/Pipeline1/working"
outdir <- "~/Dropbox/WillseyLab/CPPs/Pipeline1/output"

load(file.path(datdir, "mmAAFinal_buildDF1.RData"))
mmAA <- mmAAtoSave
dim(mmAA) #61440 12
sum(!duplicated(mmAA$gId)) #22732
gId <- unique(mmAA$gId) #22732

top <- topQgId[topQgId %in% gId]
all <- allQgId[allQgId %in% gId]
length(top) #4592
length(all) #18188


#########Pipeline1: mac top versus mac total


#####################################################
#######filtered by R content in sliding window of 15; mmAAR
load(file.path(wrkdir, "mmAAfilteredbyR.RData"))
#Enrichment
gIdR <- unique(mmAAR$gId) 
length(gIdR) #12994
sum(all %in% gIdR) #10990
sum(top %in% gIdR) #2825
t <- 4592
o <- 18188-4592
h <- sum(topmacgId %in% gIdR)
s <- sum(allmacgId %in% gIdR) #10990
phyper(2825-1, 4592, 18188-4592, 10990, lower.tail=FALSE)
######################################################

######################################################
#filtered by SRP; mmAASP
load(file.path(wrkdir, "mmAAfilteredbyR_SRP.RData"))
#Enrichment
gIdSP <- unique(mmAASP$gId) 
length(gIdSP) #1958
sum(all %in% gIdSP) #1715
sum(top %in% gIdSP) #312
phyper(312-1, 4592, 18188-4592, 1715, lower.tail=FALSE)
######################################################


#####################################################
#######filtered by presence of at least 1 H in prediction; mmAASS
#load(file.path(wrkdir, "mmAAfilteredbyR_SRP_H.RData"))
#Enrichment
#gId <- unique(mmAASS$gId) 
#length(gId) #1379
#sum(topmacgId %in% gId) #293
#sum(allmacgId %in% gId) #1211

#phyper(293-1, 5831, 18188-5831, 1211)
######################################################


#####################################################
#######only those I have predictions for; mmAApred
load(file.path(wrkdir, "mmAAfilteredbyR_SRP_H_pred.RData"))
#Enrichment
gIdpred <- unique(mmAApred$gId) 
length(gIdpred) #1378
sum(all %in% gIdpred) #1211
sum(top %in% gIdpred) #293
phyper(238-1, 4592, 18188-4592, 1211, lower.tail=FALSE)
######################################################


#####################################################
#######Filtered by arg content on helixes
load(file.path(wrkdir, "mmAAfilteredbyR_SRP_H_RonH.RData"))
#Enrichment
gIdRH <- unique(mmAAfiltered$gId) 
length(gIdRH) #934
sum(all %in% gIdRH) 
sum(top %in% gIdRH)
phyper(163-1, 4592, 18188-4592, 838, lower.tail=FALSE)
######################################################

# PIPELINE2
#load dataframes
allQuart <- read.delim(file.path(datdir, "allQuartMacGenes.txt"))
topQuart <- read.delim(file.path(datdir, "allMacGeneMaxes.txt"))
allgId <- unique(allQuart$gId)
topgId <- unique(topQuart$gId)

#DF1: mmAA
mmAA <- read.delim("~/Dropbox/WillseyLab/CPPs/MmAAtable.txt")
gId1 <- unique(mmAA$gId)
length(gId1)
sum(allgId %in% gId1) #18188
sum(topgId %in% gId1) #5831


datdir2 <- '~/Dropbox/WillseyLab/CPPs/Pipeline2/working'
#DF2: mmSP
load(file.path(datdir2, "mmSP.RData"))
gId2 <- unique(mmSP$gId)
length(gId2)
sum(all %in% gId2)
a <- sum(all %in% gId2)
sum(top %in% gId2)
t <- sum(top %in% gId2)

phyper(t-1, 4592, 18188-4592, a, lower.tail=FALSE)


#DF3 mmAAlongFULL - all with sublength >= 100 (no SP filtering yet)
load(file.path(datdir2, "mmAAlongFULL.RData"))
gId3 <- unique(mmAAlongFULL$gId)
length(gId3)
sum(all %in% gId3)
a <- sum(all%in% gId3)
sum(top %in% gId3)
t <- sum(top %in% gId3)


#DF4 - SP in mmaalongfull
gId4 <- gId3[gId3 %in% gId2]
length(gId4)
sum(all %in% gId4)
sum(top %in% gId4)
phyper(601, 4459, 17819-4459, 3154, lower.tail=FALSE)


#DF5 mmAAlong for background, df=RKdf
load(file.path(datdir2, "RKdftId2gId.RData"))
gId5 <- unique(RKtIds$gId)
length(gId5)
sum(all %in% gId5)
a <- sum(allgId %in% gId5)
sum(top %in% gId5)
t <- sum(top %in% gId5)
phyper(562-1, 4459, 17819-4459, 2943, lower.tail=FALSE)

#All Gids with SS predictions and mac exp data for long reads
#skip for now, removing this part of the analysis
load(file.path(datdir2, "mmAAWithStructPredsOnly.RData"))
load(file.path(datdir2, "struct2ryAll.RData"))
load(file.path(datdir2, "tId2gId.RData"))
SSgIds <- left_join(struct2ryAll, tId2gId, by=c("fName"="tId"))
dim(SSgIds) #8029 3
SSgIds <- unique(SSgIds$gId)
SSgIds <- SSgIds[SSgIds %in% mmAAlongFULL$gId]
length(SSgIds) #2986
sum(allgId %in% SSgIds) #2437
sum(topgId %in% SSgIds)

#DF5: CPP.pred
load(file.path(datdir2, "mmAAWithStructPredsOnly.RData"))
gId6 <- unique(CPP.pred$gId)
length(gId6)
sum(all %in% gId6)
a <- sum(all %in% gId6)
sum(top %in% gId6)
t <- sum(top %in% gId6)
#phyper(t-1, 600, 2409-60, a)
phyper(t-1, 4459, 17819-4459, a, lower.tail=FALSE)

#DF7 mmAAshortFULL - all with sublength >= 100 (no SP filtering yet)
load(file.path(datdir2, "mmAAshortFULL.RData"))
gId7 <- unique(mmAAshortFULL$gId)
length(gId7)
sum(all %in% gId7)
a <- sum(allgId %in% gId7)
sum(top%in% gId7)
t <- sum(top %in% gId7)

#DF8 mmAAshort and SPs
gId8 <- gId7[gId7 %in% gId2] #gId2: SRPs
length(gId8)
sum(all %in% gId8)
sum(top %in% gId8)
phyper(185-1, 1496, 4466-1496, 714, lower.tail=FALSE)

#DF9 RKdfShort
load(file.path(datdir2, "RKdfShort.RData"))
head(RKdfShort)
gId9 <- unique(RKdfShort$gId)
length(gId9) #296
sum(all %in% gId9)
a <- sum(all %in% gId9)
sum(top %in% gId9)
t <- sum(top %in% gId9)
phyper(52-1, 1496, 4466-1496, 227, lower.tail=FALSE)


#filter to short ones with mac exp and sec strt pred
#load(file.path(datdir2, "struct2ryAll.RData"))
#load(file.path(datdir2, "tId2gId.RData"))
#SSgIds <- left_join(struct2ryAll, tId2gId, by=c("fName"="tId"))
#SSgIds <- unique(SSgIds$gId)
#length(SSgIds)
#SSgIdShort <- SSgIds[SSgIds %in% gId6] #601
#sum(SSgIdShort %in% allgId)
#sum(SSgIdShort %in% topgId)



#DF10:CPPshort.pred - equivalent of CPP.pred
load(file=file.path(datdir2, "CPPshort_pred.RData"))
load(file.path(datdir2, "tId2gId.RData"))
CPPshort.pred <- left_join(CPPshort.pred, tId2gId)
gId10 <- unique(CPPshort.pred$gId)
length(gId10) #98
sum(gId10 %in% all) #70
sum(gId10 %in% top) #19
#phyper(19-1, 166, 493-166, 70)
phyper(18-1, 1496, 4466-1496, 70, lower.tail=FALSE)

#allPepsRK
load(file.path(wrkdir, "allPepsRK.RData"))
#load(file.path(datdir2, "tId2gId.RData"))
#allPepsRK<- left_join(allPepsRK, tId2gId)
gId11 <- unique(allPepsRK$gId)
length(gId11)
sum(gId11 %in% all) 
sum(gId11 %in% top) 
phyper(251-1, 4459, 18188-4459, 1209, lower.tail=FALSE)




