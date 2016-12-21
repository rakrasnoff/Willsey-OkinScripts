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

























################################################################################
#load in files from website
#download control files to folder, read in all files in folder, merge and average
getwd()
wrkdir <- "/Users/rebeccakrasnoff/Documents/Willsey-OkinProject/GSE21764_family.xml"
setwd(wrkdir)
#get files with .txt extension

#file_list <- Sys.glob("*.txt")
#files <- lapply(file_list, function(x) read.delim(x, header=FALSE))
#controlInd <- seq(2, 26, 2)
#controls <- files[controlInd]
#controls[1]

wrkdir <- "/Users/rebeccakrasnoff/Documents/Willsey-OkinProject/GSE21764_family.xml"
find_controls = function(mypath){
  filenames=Sys.glob("*.txt")
  files <- lapply(filenames, function(x) read.delim(x, header=FALSE))
  controlInd <- seq(2, 26, 2)
  controls <- files[controlInd] }

control_files = find_controls(wkdir)
control_files
tempdat <- left_join(as.data.frame(control_files[1]), as.data.frame(control_files[2]), by = "V1")
head(tempdat)

tempdat <- left_join(as.data.frame(tempdat), as.data.frame(control_files[3]), by = "V1")
tempdat <- left_join(as.data.frame(tempdat), as.data.frame(control_files[4]), by = "V1")
tempdat <- left_join(as.data.frame(tempdat), as.data.frame(control_files[5]), by = "V1")
tempdat <- left_join(as.data.frame(tempdat), as.data.frame(control_files[6]), by = "V1")
tempdat <- left_join(as.data.frame(tempdat), as.data.frame(control_files[7]), by = "V1")
tempdat <- left_join(as.data.frame(tempdat), as.data.frame(control_files[8]), by = "V1")
tempdat <- left_join(as.data.frame(tempdat), as.data.frame(control_files[9]), by = "V1")
tempdat <- left_join(as.data.frame(tempdat), as.data.frame(control_files[10]), by = "V1")
tempdat <- left_join(as.data.frame(tempdat), as.data.frame(control_files[11]), by = "V1")
tempdat <- left_join(as.data.frame(tempdat), as.data.frame(control_files[12]), by = "V1")

# find_exp = function(mypath){
#   filenames=Sys.glob("*.txt")
#   files <- lapply(filenames, function(x) read.delim(x, header=FALSE))
#   expInd <- seq(1, 25, 2)
#   exps <- files[expInd] }
# 
# exp_files = find_exp(wkdir)
# summary(exp_files)
# head(tempdat)
# tempdat <- left_join(as.data.frame(tempdat), as.data.frame(exp_files[1]), by = "V1")
# tempdat <- left_join(as.data.frame(tempdat), as.data.frame(exp_files[2]), by = "V1")
# tempdat <- left_join(as.data.frame(tempdat), as.data.frame(exp_files[3]), by = "V1")
# tempdat <- left_join(as.data.frame(tempdat), as.data.frame(exp_files[4]), by = "V1")
# tempdat <- left_join(as.data.frame(tempdat), as.data.frame(exp_files[5]), by = "V1")
# tempdat <- left_join(as.data.frame(tempdat), as.data.frame(exp_files[6]), by = "V1")
# tempdat <- left_join(as.data.frame(tempdat), as.data.frame(exp_files[7]), by = "V1")
# tempdat <- left_join(as.data.frame(tempdat), as.data.frame(exp_files[8]), by = "V1")
# tempdat <- left_join(as.data.frame(tempdat), as.data.frame(exp_files[9]), by = "V1")
# tempdat <- left_join(as.data.frame(tempdat), as.data.frame(exp_files[10]), by = "V1")
# tempdat <- left_join(as.data.frame(tempdat), as.data.frame(exp_files[11]), by = "V1")
# tempdat <- left_join(as.data.frame(tempdat), as.data.frame(exp_files[12]), by = "V1")
# 
# head(tempdat)
# names(tempdat) <- c("ids", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8", "e9", "e10", "e11", "e12")
head(tempdat)
dim(tempdat)

names(tempdat) <- c("ids", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12")
expDat <- tempdat
head(expDat)
#GET INFO ABOUT ILMN IDS

geneInfo <- read.delim(file.path(datdir, "GPL6885-11608.txt"))
head(geneInfo)
colnames(geneInfo)
nrow(geneInfo)
geneInfo <- geneInfo[,c(1,9,12,13)]
colnames(geneInfo)
head(geneInfo)

#CONVERT ILMN IDS TO ENTREZ IDS USING geneInfo
ilmnIDs <- as.vector(expDat$ids)
#expDatRN <- data.frame(expDat, ilmnIDs)
macExpDat <- left_join(geneInfo, expDat, by = c("ID" = "ids"))
head(macExpDat)
macExpTemp <- subset(macExpDat, macExpDat$Entrez_Gene_ID != "NA")
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
