#Script for Filtering Full Data Frame
#Important data structures
#mmAA: contains unfiltered list of genes with corresponding information (SRP, seq, etc)
#mmAAR: mmAA, filtered by R content
##mmAAME: mmAA, filtered by R content and macrophage expression
##mmAASP: mmAA, filtered by R content, macrophage expression, and SRP presence
##missingFasta: all tIds in mmAASP that do not have secondary structure predictions; find here:
    #datdir, missingPreds953.Fa
##mmAASS: mmAASP minus any peptides with no helixes; brings me to 2860
##mmAApred: mmAASS, minus any peptides missing secondary structure predictions, 1907 x 14
##hexSplitList: output of function hexSplit in list form; contains individual helixes from each 
  ##transcript in mmAApred
##hexSplitdf: dataframe form of hexSplitList
##mmAA2cut: with all cut pieces of helixes and sequences, then filtered by R content on hexes; final DF

library(seqRFLP)

#############################################################################
#1. LOAD DATA FRAME
datdir <- "~/Dropbox/WillseyLab/CPPs/Pipeline1/data"
wrkdir <- "~/Dropbox/WillseyLab/CPPs/Pipeline1/working"
outdir <- "~/Dropbox/WillseyLab/CPPs/Pipeline1/output"

#Load macrophage expression data for enrichment analysis
maxMac <- read.delim(file.path(datdir, "topQuartMacGenes.txt"))

#mmAA <- read.delim(file.path(datdir, "mmAAFinal_buildDF1.txt"))
load(file.path(datdir, "mmAAFinal_buildDF1.RData"))
mmAA <- mmAAtoSave
head(mmAA)
names(mmAA)
dim(mmAA) #52659 x 13

#############################################################################
#2. Filter by R content (sliding window)
#mmAA$maxRper15 <- as.list(mmAA$maxRper15)
head(mmAA)
names(mmAAR) <- c("tId", "gId", "seq", "subseq", "length", "Rc", "Rp", "sublength", "maxRper8", "macExp", "sigPep", "Prediction", "Confidence")
mmAAR <- mmAA[mmAA$maxRper15 > 0, ]
dim(mmAA)
dim(mmAAR) #51875 x 13
head(mmAAR)
mmAAR$keep <- 0
mmAAR$keep[mmAAR$length < 15] <- 1
mmAAR$keep[mmAAR$length >= 15 & mmAAR$maxRper15 >= 4] <- 1
mmAAR <- subset(mmAAR, keep == 1)
dim(mmAAR) #25525 x 14
head(mmAAR)

save(mmAAR, file=file.path(wrkdir, "mmAAfilteredbyR.RData"))
#############################################################################
#3.Filter by Mac Exp
#mmAAME <- mmAA #to filter original list
mmAAME <- mmAAR # to continue filtering previously filtered list
mmAAME$keep <- 0
mmAAME$keep[is.na(mmAAME$macExp)] <- 1
mmAAME$keep[!is.na(mmAAME$macExp) & mmAAME$macExp >= 122.5] <- 1
sum(mmAAME$keep)
#20684 (keeping those with NAs)

############4. Filter by Signal Peptides
#mmAASP <- mmAAME
mmAASP <- mmAAR
mmAASP <- subset(mmAASP, sigPep == "Y")
dim(mmAASP)
#2938

save(mmAASP, file=file.path(wrkdir, "mmAAfilteredbyR_SRP.RData"))
#saved version: no macrophage expression filter
#which of these do I have secondary structure predictions for
head(mmAASP)
sum(is.na(mmAASP$Prediction) == T) #953 are missing
missingPreds <- mmAASP[is.na(mmAASP$Prediction),]
dim(missingPreds)
missingFasta = dataframe2fas(missingPreds[,c(1,4)], file=file.path(datdir, "missingPreds953.Fa"))

#############################################################################
#3.Filter by Secondary Structure
head(mmAASP)
mmAASS <- mmAASP
mmAASS$helix <- 0
maxRper15 <- data.frame(mmAASS$tId, mmAASS$maxRper15)
mmAASS <- mmAASS[,-9] #removing maxes, saved to their own dataframe
mmAASS <- mutate(mmAASS, helix = str_count(mmAASS$Prediction, "H"))
head(mmAASS)

mmAASS <- mmAASS[mmAASS$helix != 0,]
dim(mmAASS) #this brings me to 2860; however, 953 are still NA

save(mmAASS, file=file.path(wrkdir, "mmAAfilteredbyR_SRP_H.RData"))

#####Keeping only those with SS predictions now
mmAApred <- filter(mmAASS, is.na(Prediction) == F)
save(mmAApred, file=file.path(wrkdir, "mmAAfilteredbyR_SRP_H_pred.RData"))

hexSplit <- function(Prediction) {
df <- strsplit(Prediction, "")
df <- data.frame(df)
df$pos <- seq(1:nrow(df))
names(df) <- c("hex", "pos")
df$hex[df$hex!="H"] <- "-"
test <- array(0, nrow(df))
test <- data.frame(start=test, end=test)
if (df$hex[1] == "H") {
  test$start[1] = 1
}
if(df$hex[nrow(df)] == "H") {
  test$end[nrow(df)]= nrow(df)
}
for (i in 1:(nrow(df)-1)) {
  if (df$hex[i] == "-" & df$hex[i+1] == "H") 
  {test$start[i] <- df$pos[i+1]}
  if (df$hex[i] == "H" & df$hex[i+1] == "-")
  {test$end[i] <- df$pos[i]}
}
test <- data.frame(test)
testhex<-data.frame(start=test$start[test$start!="0"])
testhex$end <- test$end[test$end != "0"]
testhex
}

hexSplitList <- lapply(1:nrow(mmAApred), function(x) hexSplit(mmAApred$Prediction[x]))
names(hexSplitList) <- mmAApred$tId[1:nrow(mmAApred)]
hexSplitdf <- do.call(rbind, hexSplitList)
hexSplitdf$tId <- do.call(rbind, strsplit(rownames(hexSplitdf), "[.]"))[,1]
head(hexSplitdf)
dim(hexSplitdf)

save(hexSplitdf, file=file.path(wrkdir, "hexSplitdf.RData"))

#merge info with sequences
mmAA2merge <- data.frame(tId = mmAApred$tId, subseq = mmAApred$subseq)
mmAA2cut <- full_join(hexSplitdf, mmAA2merge)
dim(mmAA2cut) #12319
dim(trialdf)
mmAA2cut <- mutate(mmAA2cut, hseq = substring(subseq, start, end))
head(mmAA2cut)
mmAA2cut <- mutate(mmAA2cut, Rs = str_count(hseq, "R"))
mmAA2cut <- mutate(mmAA2cut, length = str_count(hseq))
mmAA2cut <- mutate(mmAA2cut, percent = (Rs/length)*100)
mmAA2cut <- filter(mmAA2cut, percent > 0)
dim(mmAA2cut) #5928 x 8
length(unique(mmAA2cut$tId)) #1589 

countRs <- function(mmAA) {
  aastring <- AAString(mmAA)
  nsR <- letterFrequencyInSlidingView(aastring, 8, "R") }

mmAA2cut$maxRper8 <- 0
mmAA2cut$maxRper8[mmAA2cut$length < 8] <- mmAA2cut$Rs[mmAA2cut$length < 8]
mmAA2cut$maxRper8[mmAA2cut$length >= 8] <- lapply(lapply(mmAA2cut$hseq[mmAA2cut$length >= 8], countRs), max)
head(mmAA2cut)

mmAA2cut$keep <- 0
mmAA2cut$keep[mmAA2cut$length < 8] <- 1
mmAA2cut$keep[mmAA2cut$length >= 8 & mmAA2cut$maxRper8 >= 2] <- 1
mmAA2cut <- subset(mmAA2cut, keep == 1)
dim(mmAA2cut) #2730
length(unique(mmAA2cut$tId)) #1256
#merge with other data
maxRper8 <- data.frame(tId2=mmAA2cut$tId, maxRper8=do.call(rbind, mmAA2cut$maxRper8))
mmAA2cut <- mmAA2cut[,-9]

mmAAfiltered <- left_join(mmAA2cut, mmAApred, by=c("tId"="tId"))
dim(mmAA2cut)
dim(mmAAfiltered)
mmAAfiltered <- cbind(mmAAfiltered, maxRper8)
sum(mmAAfiltered$tId == mmAAfiltered$tId2)
dim(mmAAfiltered)
head(mmAAfiltered)
mmAAfiltered <- data.frame(gId=mmAAfiltered$gId, tId=mmAAfiltered$tId, subseq=mmAAfiltered$subseq.x, sublength=mmAAfiltered$sublength, sigPep=mmAAfiltered$sigPep, macExp=mmAAfiltered$macExp, Prediction=mmAAfiltered$Prediction, helix=mmAAfiltered$helix, hseq=mmAAfiltered$hseq, hlength=mmAAfiltered$length.x, hexRs=mmAAfiltered$Rs, percentRhex=mmAAfiltered$percent, maxRper8hex=mmAAfiltered$maxRper8)
dim(mmAAfiltered)
sum(!duplicated(mmAAfiltered$tId))
sum(!duplicated(mmAAfiltered$gId))

tIds <- unique(mmAAfiltered$tId) 

save(mmAAfiltered, file=file.path(wrkdir, "mmAAfilteredbyR_SRP_H_RonH.Rdata"))

####
##Save results
save(mmAAfiltered, file=file.path(outdir, "CPPpredictionsHelixSlidingWindowOf8.RData"))
load(file=file.path(outdir, "CPPpredictionsHelixSlidingWindowOf8.RData"))
#mmAA2cut$maxRper8 <- do.call(rbind, mmAA2cut$maxRper8)
write.table(mmAAfiltered, file.path(outdir, "Pipeline1Output.txt"), sep = "\t", quote = F)




###check enrichment
mouseCPPs <- read.delim(file.path(datdir, "mouseCPPs.txt"))
sum(CPPs$gId %in% mouseCPPs$mouseEnsembl)
sum(CPPs$gId %in% unique(mouseCPPs$mouseEnsembl)) #same, still 5, and my list has unique tIds, but not gIds. 
CPPs_unique <- CPPs[!duplicated(CPPs$gId),] #935 gene Ids
mouseEns_unique <- unique(mouseCPPs$mouseEnsembl)

q <- 4 #number of hits in my list -1
m <- 18 #number of known CPPs: nrow(mouseCPPs)
n <- #61440 unique tIds, 22732 unique gIds
k <- 935 #length of my list; nrow(allPeps1 or allPeps2)


phyper(4, 18, 22732-18, 935) #.999


18/22732 #.000792
5/935 #.00535

sum(mouseEns_unique %in% missingPreds$gId)

