#### This script begins by loading macrophage expression data. The second part uses this data for filtering. 


#Filter by macrophage expression
library(Biobase)
library(GEOquery)
options(stringsAsFactors = TRUE)

datdir <-"~/Documents/Willsey-OkinProject/MacExpDat"
wrkdir <-"~/Documents/Willsey-OkinProject/MacExpDat"
outdir <-"~/Documents/Willsey-OkinProject/MacExpDat"

info <- getGEO("GSE56711", GSEMatrix=FALSE)
Meta(info)

gse <- getGEO("GSE56711", GSEMatrix=TRUE)
show(gse)
#annotation(gse)
#colnames(fData(gse))
expSet <- exprs(gse[[1]])
show(pData(phenoData(gse[[1]]))[1:5,c(1,6,8)])
head(expSet)
nrow(expSet)

identical(colnames(expSet), rownames(pData(phenoData(gse[[1]]))))
colnames(expSet) <- pData(phenoData(gse[[1]]))$title
head(expSet)
expDat <- as.data.frame(expSet)
head(expDat)

#GET INFO ABOUT ILMN IDS

geneInfo <- read.delim(file.path(datdir, "illIDforGSE56711.txt"))
head(geneInfo)
colnames(geneInfo)
nrow(geneInfo)
geneInfo <- geneInfo[,c(1,9,12,13)]
colnames(geneInfo)
head(geneInfo)

#CONVERT ILMN IDS TO ENTREZ IDS USING geneInfo
expDat$ID <- rownames(expDat)
macExpDat <- left_join(geneInfo, expDat)
head(macExpDat)
macExpTemp <- subset(macExpDat, macExpDat$Entrez_Gene_ID != "NA")
macExpTemp$Avg <- apply(macExpTemp[,-c(1:4)], 1, mean)
head(macExpTemp)
dim(macExpTemp)
macExpDatAvg <- macExpTemp[,c(2,3,15)]
macExpDat <- macExpTemp
head(macExpDat)

write.table(macExpDatAvg, file=file.path(outdir, "macExpDatAvg.txt"), sep = "\t", quote = F)
write.table(macExpDat, file = file.path(outdir, "macExpDat.txt"), sep = "\t", quote = F)

################### Filtering ##################
####################
#Look at expression in macrophages
datdir <-"~/Dropbox/WillseyLab/CPPs/MacExpDat"
wrkdir <-"~/Documents/Willsey-OkinProject/MacExpDat"
outdir <-"~/Dropbox/WillseyLab/CPPs/MacExpDat"

#first, need to translate macrophage expression data to ensembl IDs
macExpDat <- read.delim(file.path(datdir, "macExpDat.txt"))
macExpDat$Entrez_Gene_ID <- as.character(macExpDat$Entrez_Gene_ID)

mkeys <- (keys(org.Mm.eg.db, keytype="ENTREZID"))
cols <- c("ENSEMBL")
key <- select(org.Mm.eg.db, keys=mkeys, columns = cols, keytype="ENTREZID")
head(key)
dim(key)

macExpEns <- left_join(macExpDat, key, by=c("Entrez_Gene_ID" = "ENTREZID"))
head(macExpEns)
tail(macExpEns)

write.table(macExpEns, file.path(outdir, "macExpEnsembl.txt"), sep = "\t", quote=F)

# #now check expression
# macExpEns <- macExpEns[,c(2:14,16)]
# macExpList <- list(macExpEns[,c(4:13)])
# summary(unlist(macExpList)) #getting quartiles
# 
# ##now just looking at maxes
# maxMacExp <- apply(macExpEns[,c(4:13)], 1, max)
# summary(maxMacExp)
# ###
# 
# checkMacExp <- left_join(AAbyRc30Trim, macExpEns, by=c("gId" = "ENSEMBL"))
# head(checkMacExp)
# dim(checkMacExp)
# checkMacExp <- checkMacExp[!duplicated(checkMacExp$seq),]
# head(checkMacExp)
# dim(checkMacExp)
# DF <- data.frame(matrix(unlist(checkMacExp), nrow=nrow(checkMacExp)),stringsAsFactors=FALSE)
# names(DF) <- names(checkMacExp)
# head(DF)
# checkMacExp <- DF
# checkMacExp <- checkMacExp[complete.cases(checkMacExp),]
# dim(checkMacExp)
# maxMacExp <- apply(checkMacExp[,c(10:20)], 1, max)
# filtMacExp <- filter(checkMacExp, maxMacExp >= 122.1)
# dim(filtMacExp)
# 
# ###check enrichment
# sum(mouseCPPs$mouseEnsembl %in% checkMacExp$gId)
# /nrow(macExpEns)
# sum(mouseCPPs$mouseEnsembl %in% filtMacExp$gId)/nrow(filtMacExp)
# filtMacExp[filtMacExp$gId %in% mouseCPPs$mouseEnsembl,]
# 
