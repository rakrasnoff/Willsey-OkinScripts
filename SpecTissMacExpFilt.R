#############LOOKING FOR GENES EXPRESSED IN PARTICULAR TISSUES##########
datdir <- "~/Dropbox/WillseyLab/CPPs/data"
macMaxDat <- read.delim(file.path(datdir, "macExpMaxesPerTissue.txt"))
dim(macMaxDat)
head(macExpDat)
mac <- data.frame(gId = macMaxDat$ENSEMBL, perit1 = as.numeric(macMaxDat$Peritoneum_untreated_rep1), perit2=as.numeric(macMaxDat$Peritoneum_untreated_rep2),
                  lung = as.numeric(macMaxDat$Lung_untreated_rep1), liver=as.numeric(macMaxDat$Liver_untreated_rep1), spleen = as.numeric(macMaxDat$Spleen_untreated_rep1),
                  intest1=as.numeric(macMaxDat$Intestine_untreated_rep1), intest2=as.numeric(macMaxDat$Intestine_untreated_rep2), adipose1=as.numeric(macMaxDat$Adipose_untreated_rep1),
                  adipose2=as.numeric(macMaxDat$Adipose_untreated_rep2), bmdm=as.numeric(macMaxDat$BMDM_untreated_rep1))
head(mac)
dim(mac)

mact <- data.frame(t(mac[,-1]))
head(mact[,c(1:10)])
colnames(mact) <- mac$gId
head(mact[,c(1:10)])


#correlations
lungV <- c("0", "0", "1", "0", "0", "0", "0", "0", "0", "0")
lungV <- as.numeric(lungV)
peritV <- as.numeric(c("1", "1", "0", "0", "0", "0", "0", "0", "0", "0"))
liverV <- as.numeric(c("0", "0", "0", "1", "0", "0", "0", "0", "0", "0"))
spleenV <- as.numeric(c("0", "0", "0", "0", "1", "0", "0", "0", "0", "0"))
intestV <- as.numeric(c("0", "0", "0", "0", "0", "1", "1", "0", "0", "0"))
adiposeV <- as.numeric(c("0", "0", "0", "0", "0", "0", "0", "1", "1", "0"))
bmdmV <- as.numeric(c("0", "0", "0", "0", "0", "0", "0", "0", "0", "1"))

lungcor <- apply(mact, 2, function(x) cor(lungV, x, method="spearman"))
peritcor <- apply(mact, 2, function(x) cor(peritV, x, method="spearman"))
livercor <- apply(mact, 2, function(x) cor(liverV, x, method="spearman"))
spleencor <- apply(mact, 2, function(x) cor(spleenV, x, method="spearman"))
intestcor <- apply(mact, 2, function(x) cor(intestV, x, method="spearman"))
adiposecor <- apply(mact, 2, function(x) cor(adiposeV, x, method="spearman"))
bmdmcor <- apply(mact, 2, function(x) cor(bmdmV, x, method="spearman"))

lungcor <- data.frame(gId=names(lungcor), cor=lungcor)
peritcor <- data.frame(gId=names(peritcor), cor=peritcor)
livercor <- data.frame(gId=names(livercor), cor=livercor)
spleencor <- data.frame(gId=names(spleencor), cor=spleencor)
intestcor <- data.frame(gId=names(intestcor), cor=intestcor)
adiposecor <- data.frame(gId=names(adiposecor), cor=adiposecor)
bmdmcor <- data.frame(gId=names(bmdmcor), cor=bmdmcor)

lungGenes <- filter(lungcor, cor >= .6)
peritGenes <- filter(peritcor, cor >= .6)
liverGenes <- filter(livercor, cor >= .6)
spleenGenes <- filter(spleencor, cor >= .6)
intestGenes <- filter(intestcor, cor >= .6)
adiposeGenes <- filter(adiposecor, cor >= .6)
bmdmGenes <- filter(bmdmcor, cor >= .6)

specGenes <- full_join(lungGenes, peritGenes, by=c("gId"="gId"))
specGenes <- full_join(specGenes, liverGenes, by=c("gId"="gId"))
specGenes <- full_join(specGenes, spleenGenes, by=c("gId"="gId"))
specGenes <- full_join(specGenes, intestGenes, by=c("gId"="gId"))
specGenes <- full_join(specGenes, adiposeGenes, by=c("gId"="gId"))
specGenes <- full_join(specGenes, bmdmGenes, by=c("gId"="gId"))
dim(specGenes)
sum(duplicated(specGenes$gId))

sum(spleenGenes$gId %in% adiposeGenes$gId) #checked several, seems to be little to no overlap between lists

write.table(specGenes, file.path(datdir, "specExpdMacGenes.txt"), sep="\t", quote=F)

#############ENRICHMENT ANALYSIS################
# PIPELINE1

#original list:
basedatdir <- "~/Dropbox/WillseyLab/CPPs/data"
datdir <- "~/Dropbox/WillseyLab/CPPs/Pipeline1/data"
wrkdir <- "~/Dropbox/WillseyLab/CPPs/Pipeline1/working"
outdir <- "~/Dropbox/WillseyLab/CPPs/Pipeline1/output"

#First, load my full dataframe and see how many specifically expressed genes overlap
load(file.path(datdir, "mmAAFinal_buildDF1.RData"))
mmAA <- mmAAtoSave
allgIds <- unique(mmAA$gId) #22732
sum(specGenes$gId %in% allgIds) #4312 of 4358
specGenes <- specGenes$gId

macGenes <- read.delim(file.path(basedatdir, "macExpMaxesPerTissue.txt")) #already unique gIds
macGenes <- macGenes$ENSEMBL
sum(macGenes %in% allgIds) #18188, consistent with previous filtering

#####################################################
#######filtered by R content in sliding window of 15; mmAAR
load(file.path(wrkdir, "mmAAfilteredbyR.RData"))
#Enrichment
sum(unique(mmAAR$gId) %in% macGenes) #10990
sum(unique(mmAAR$gId) %in% specGenes) #2715
#total = 18188
#targets = 4312
#hits = 2715
#sample size = 10990

phyper(2715-1, 4312, 18188-4312, 10990, lower.tail=FALSE)






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
