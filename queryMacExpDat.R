#source("http://www.bioconductor.org/biocLite.R")
#biocLite("GEOquery")
library(Biobase)
library(GEOquery)
options(stringsAsFactors = TRUE)

datdir <-"~/Documents/Willsey-OkinProject/MacExpDat"
wrkdir <-"~/Documents/Willsey-OkinProject/MacExpDat"
outdir <-"~/Documents/Willsey-OkinProject/MacExpDat"

info <- getGEO(filename = file.path(datdir, "GSE21764_family.soft.gz"), GSEMatrix=FALSE)
Meta(info)

gse <- getGEO(filename = file.path(datdir, "GSE21764_family.soft.gz"), GSEMatrix=TRUE)
show(gse)
#annotation(gse)
#colnames(fData(gse))
show(pData(phenoData(expSet[[1]]))[1:5,c(1,6,8)])
expDat <- (exprs(eset[[1]]))
head(expDat)
nrow(expDat)

#GET INFO ABOUT ILMN IDS

geneInfo <- read.delim(file.path(datdir, "GPL6885-11608.txt"))
head(geneInfo)
colnames(geneInfo)
nrow(geneInfo)
geneInfo <- geneInfo[,c(1,9,12,13)]
colnames(geneInfo)
head(geneInfo)

#CONVERT ILMN IDS TO ENTREZ IDS USING geneInfo
ilmnIDs <- as.vector(rownames(expDat))
expDatRN <- data.frame(expDat, ilmnIDs)
macExpDat <- merge(geneInfo, expDatRN, by.x = "ID", by.y = "ilmnIDs")
macExpTemp <- subset(macExpDat, macExpDat$Entrez_Gene_ID != "NA")
macExpTemp <- macExpTemp[,-c(seq(6, 29, 2))] #get rid of treatment conditions
macExpTemp$Avg <- apply(macExpTemp[,-c(1:4)], 1, mean)
head(macExpTemp)
macExpDat <- macExpTemp[,c(2,3,17)]
head(macExpDat)

#load list of mouse peptidases, with gene ids available
mousePeps <- read.delim(file.path(datdir, "mousePeptidases.txt"))
head(mousePeps)
#Convert to ENTREZIDs for expression analysis
library(org.Mm.eg.db)
#keytypes(org.Mm.eg.db)
#create vector of converted ids
convertedIds<- select(org.Mm.eg.db, as.vector(mousePeps$Gene), "ENTREZID", "SYMBOL")
#append vector to dataframe if applicable
mousePeps$convertedIds <- convertedIds$ENTREZID[match(mousePeps$Gene, convertedIds$SYMBOL)]
head(mousePeps)
mousePeps <- subset(mousePeps, mousePeps$convertedIds != "<NA>")
head(mousePeps,50)




