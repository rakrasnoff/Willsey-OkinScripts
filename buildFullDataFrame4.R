require(Biostrings)
options(stringsAsFactors = FALSE)
library(stringr)
library(dplyr)
library(plyr)
library(zoo)
library(seqRFLP)
library(org.Mm.eg.db)
library(tidyr)

####################################################################################################
### Set directories
####################################################################################################
SPdatdir <- "~/Dropbox/WillseyLab/CPPs/sigPepFiles" #signal peptide files
MEdatdir <- "~/Dropbox/WillseyLab/CPPs/MacExpDat" #macrophage expression data
outdir <- "~/Dropbox/WillseyLab/CPPs"


####################################################################################################
### Load dataframe
####################################################################################################
#Start by loading dataframe with mouse peptides. This list comes from GEO database. List initially 
#loaded in CPP_load_DNA script
mmAA <- read.delim("~/Dropbox/WillseyLab/CPPs/MmAAtable.txt")
mmAA_original <- mmAA #for use later on
head(mmAA)
dim(mmAA) #61440 by 3


####################################################################################################
### Remove any sequence that occurs after a stop
####################################################################################################

mmAA <- mutate(mmAA, length = str_count(seq)) #this creates a column with peptide length
mmAA <- separate(mmAA, seq, c("subseq"), remove=F, sep = "[*]", extra="drop")
head(mmAA)
tail(mmAA)

####################################################################################################
###3. Filter by presence of SRP
####################################################################################################
#Because I already have SRP data for the full peptides, I will first remove any peptides with no SRPs
files <- dir(SPdatdir)

#Read in files
SPs <- lapply(files, function(x) { 
  SPlist <- read.delim(file.path("~/Dropbox/WillseyLab/CPPs/sigPepFiles", x))
})
head(SPs[[1]])
head(SPs[[2]])

#Rename, reformat
for (i in 2:31){
  names(SPs[[i]]) <- names(SPs[[1]])
  SPs[[i]] <- SPs[[i]][,c(1:12)]
}
summary(SPs)

SPdf <- rbind.fill(SPs)
dim(SPdf)
SPdf <- SPdf[,c(1,10)]
head(SPdf)
mmAASP <- left_join(mmAA, SPdf, by=c("tId"= "name"))
dim(mmAASP) 
head(mmAASP)
names(mmAASP) <- c("tId", "gId", "seq", "subseq", "length", "sigPep")

#Now, filtering
mmSP <- filter(mmAASP, sigPep != "N")
dim(mmSP) #gives me 7771, none are duplicated
#Now, each of these 7771 sequences needs to be run through the prediction server
#save(mmSP, file=file.path(outdir, "mmSP.RData"))
load(file.path(outdir, "mmSP.RData"))
mmAA <- mmSP #7771
#Filtering by SRP round 2- I ran all of the subsequences, and none of the y/n predictions changed

####################################################################################################
###4. Filter by R and K presence  FROM HERE DOWN, ONLY HAVE INFO ON LONG TRANSCRIPTS (>100)
####################################################################################################
#mmAA <- mutate(mmAA, Mlength = str_count(mmAA$Mseq)) #this creates a column with peptide length
mmAA <- mutate(mmAA, sublength = str_count(subseq))
mmAA <- mutate(mmAA, Rc = str_count(mmAA$subseq, "R")) #number of Rs in the entire peptide
mmAA <- mutate(mmAA, Kc = str_count(mmAA$subseq, "K"))
mmAA <- mutate(mmAA, Rp = Rc/length * 100) #percentage of Rs in the peptide
mmAA <- mutate(mmAA, Kp = Kc/length *100)
head(mmAA)
dim(mmAA) #7771 x 11

#This function will count the number of Args and Lys in a window of defined size. To make it simpler to use in an apply function, 
#to change the size of the window, just alter it in the function below (rather than setting as a variable).
countRKs <- function(seq, name) {
  #seqInd <- match("seq", colnames(mmAA))
  #seq <- mmAA[,seqInd]
  #nameInd <- match("tId", colnames(mmAA))
  #name <- mmAA[,nameInd]
  n <- str_count(seq)
  l <- 100 
  seqtable <- data.frame(substring(seq, 1:(n-l+1), l:n))
  names(seqtable) <- c("seq")
  AAseqlist <- lapply(seqtable$seq, AAString)
  seqtable$pos <- rownames(seqtable)
  seqtable$tId <- name
  seqtable <- dplyr::mutate(seqtable, rnames = paste(tId,".",pos, sep=""))
  rownames(seqtable) <- seqtable$rnames
  nsR <- lapply(AAseqlist, letterFrequencyInSlidingView, 16, "RK") 
  seqtable$maxRK <- data.frame(maxRK=unlist(lapply(nsR, max))) #may work if i get rid of the data.frame here
  seqtable <- seqtable[,c(1,2,3,5)]
  seqtable
}

mmAAlong <- mmAA[mmAA$sublength >= 100,] #6916 rows

#########RKwindow <- mapply(countRKs, mmAAlong$subseq, mmAAlong$tId, SIMPLIFY=FALSE)######
#save(RKwindow, file=file.path(outdir, "RKwindowOutput.RData"))
######USE LOAD INSTEAD OF RERUNNING###############

load(file.path(outdir, "RKwindowOutput.RData"))
names(RKwindow) <- mmAAlong$tId
#names(RKtest) <- mmAAlong$tId[1:3]
sum(duplicated(mmAAlong$tId)) #none
maxes <- lapply(1:7002, function(x) max(RKwindow[[x]]$maxRK))
RKwindow2 <- RKwindow[maxes>=4]
length(names(RKwindow2)) #trims to 6295
length(names(RKwindow)) #from 7002
RKwindow2 <- lapply(RKwindow2, function(x) subset(x, maxRK >= 4))
RKwindow2 <- lapply(RKwindow2, function(x) mutate(x, maxRK = unlist(as.list(maxRK))))
#RKtest <- lapply(RKtest, function(x) mutate(x, maxRK=unlist(as.list(maxRK))))

#row.names(RKwindow2) <- names(RKwindow2)
RKdf <- rbind.fill(RKwindow2)
#dftest <- rbind.fill(RKtest)
head(RKdf)
names(RKdf) <- c("seq100", "pos", "tId", "maxRK")
dim(RKdf) #1675892       4

min(RKdf$maxRK) #checking that filter worked; it did
sum(duplicated(RKdf$tId))
sum(!duplicated(RKdf$tId)) #6295

##########Saving RKdf and mmAAlong#########
save(RKdf, mmAAlong, file=file.path(outdir, "RKdf_mmAAlong.RData"))
load(file.path(outdir, "RKdf_mmAAlong.RData"))
####################################################################################################
###5. Add alpha helix content
####################################################################################################

#read in secondary structure results
struct2ry <- read.delim("~/Dropbox/WillseyLab/CPPs/2ndryStruct/merged_jpred_results_2-1-17.txt")
struct2ry2 <- read.delim("~/Dropbox/WillseyLab/CPPs/2ndryStruct/merged_jpred_results_batch2_2-16-17.txt")
struct2ry3 <- read.delim("~/Dropbox/WillseyLab/CPPs/2ndryStruct/merged_partial_results_03-16-17.txt") #6011
head(struct2ry)
head(struct2ry2)
head(struct2ry3)
class(struct2ry$Prediction)
struct2ry$Prediction <- gsub(",","", struct2ry$Prediction)
struct2ry$Confidence <- gsub(",","", struct2ry$Confidence)
struct2ry2$Prediction <- gsub(",","", struct2ry2$Prediction)
struct2ry2$Confidence <- gsub(",","", struct2ry2$Confidence)
names(struct2ry)
sum(struct2ry$fName %in% struct2ry2$Name) #7
struct2ry2 <- struct2ry2[!struct2ry2$Name %in% struct2ry$fName,]
dim(struct2ry2) #1919, 3
names(struct2ry2) <- c("fName", "Prediction", "Confidence")
struct2ryAll <- rbind(struct2ry, struct2ry2) #have 3885 rows it seems - is this how many I had before?
dim(struct2ryAll) #3885
sum(struct2ry3$Name %in% struct2ryAll$fName) #1867 overlap; will keep earlier ones
struct2ry3 <- struct2ry3[!(struct2ry3$Name %in% struct2ryAll$fName),] #4144 left over
dim(struct2ry3) #4144, 3
names(struct2ry3) <- c("fName", "Prediction", "Confidence")
struct2ryAll <- rbind(struct2ryAll, struct2ry3)
dim(struct2ryAll) #8029, 3

#combine with mmAA dataframe
CPPdf <- left_join(RKdf, struct2ryAll, by = c("tId"="fName"))
dim(CPPdf) #1675892, 6
CPPdf <- left_join(CPPdf, mmAAlong, by=c("tId"="tId"))
dim(CPPdf) #1675892, 16 
head(CPPdf)
CPP <- data.frame(tId=CPPdf$tId, gId=CPPdf$gId, fullseq=CPPdf$seq, fulllength=CPPdf$length, subseq=CPPdf$subseq, sublength=CPPdf$sublength, 
                  sigPep=CPPdf$sigPep, Rc=CPPdf$Rc, Kc=CPPdf$Kc, seq100=CPPdf$seq100, 
                  pos=CPPdf$pos, maxRK=CPPdf$maxRK, Prediction=CPPdf$Prediction, Confidence=CPPdf$Confidence)
dim(CPP) #1675892, 14

head(CPP, 100)
####################################################################################################
###6. Filter by alpha helixes
####################################################################################################

head(CPP)
CPP.pred <- CPP[is.na(CPP$Prediction)==F,] #this filters out remaining peptides with no structure predictions
dim(CPP)
dim(CPP.pred) #777889
sum(!duplicated(CPP$tId)) #6295
sum(!duplicated(CPP.pred$tId)) #4881 - meaning I have predictions for 4881 peptides here. Above, I have predictions for 8029 peptides; however, 
#only 4883 of these are for peptides >100 AA. Thus, I am missing predictions for >1000 peptides that contain a signal peptide, and are at least 100 AA long. 

######IMPORTANT: THESE PEPTIDES DO NOT HAVE SECONDARY STRUCTURE PREDICTIONS###########
missing <- CPP[!(CPP$tId %in% CPP.pred$tId),]
sum(!duplicated(missing$tId))
write.table(missing, file.path(outdir, "longPepMissing2ndryPred.txt"), sep = "\t", quote=F)
###################################################################################

#Now, I will extract all unique sequences:
CPP.u <- CPP.pred[!(duplicated(CPP.pred$tId)),]
dim(CPP.u) #4881 x 14
#Get rid of seq100 - here, it is just giving me a random seq100, based on whatever was left over
CPP.u <- CPP.u[,-c(10,11,12)] #also getting rid of position and maxRK for seq100
dim(CPP.u) #4881 x 11


#Now, I will run hexsplit on CPP.u
hexSplit <- function(Prediction) {
  pred.df <- strsplit(Prediction, "")
  pred.df <- data.frame(pred.df)
  names(pred.df) <- "hex"
  pred.df <- data.frame(pred.df[(pred.df$hex != ","),])
  pred.df$pos <- seq(1:nrow(pred.df))
  names(pred.df) <- c("hex", "pos")
  pred.df$hex[pred.df$hex!="H"] <- "-"
  zeros <- array(0, nrow(pred.df))
  pos.df <- data.frame(start=zeros, end=zeros)
  if (pred.df$hex[1] == "H") {
    pos.df$start[1] = 1
  }
  if(pred.df$hex[nrow(pred.df)] == "H") {
    pos.df$end[nrow(pred.df)]= nrow(pred.df)
  }
  for (i in 1:(nrow(pred.df)-1)) {
    if (pred.df$hex[i] == "-" & pred.df$hex[i+1] == "H") 
    {pos.df$start[i] <- pred.df$pos[i+1]}
    if (pred.df$hex[i] == "H" & pred.df$hex[i+1] == "-")
    {pos.df$end[i] <- pred.df$pos[i]}
  }
  pos.df <- data.frame(pos.df)
  pos.hex<-data.frame(start=pos.df$start[pos.df$start!="0"])
  pos.hex$end <- pos.df$end[pos.df$end != "0"]
  pos.hex
}

#splitHexes <- lapply(CPP.u$Prediction, hexSplit)
#save(splitHexes, file="~/Dropbox/WillseyLab/CPPs/hexsplitoutput.txt")
load("~/Dropbox/WillseyLab/CPPs/hexsplitoutput.txt")

names(splitHexes) <- CPP.u$tId
splitHexesDf <- do.call(rbind, splitHexes)
save(splitHexesDf, file=file.path(outdir, "splitHexesDf.RData"))
head(splitHexesDf)
splitHexesDf$tId <- do.call(rbind, strsplit(rownames(splitHexesDf), "[.]"))[,1]
head(splitHexesDf)
dim(splitHexesDf) #27316      3

#merge info with sequences
mmAA2merge <- data.frame(tId = CPP.u$tId, subseq = CPP.u$subseq, sublength=CPP.u$sublength)
mmAA2cut <- full_join(splitHexesDf, mmAA2merge)
dim(mmAA2cut)
dim(splitHexesDf)
mmAA2cut <- filter(mmAA2cut, start != "NA")
dim(mmAA2cut)
mmAA2cut <- mutate(mmAA2cut, start50 = start-50)
mmAA2cut <- mutate(mmAA2cut, end50 = end+50)
mmAA2cut$cut1[mmAA2cut$start50 <= 0] <- 1
mmAA2cut$cut1[mmAA2cut$start50 > 0] <- mmAA2cut$start50[mmAA2cut$start50 > 0]
mmAA2cut$cut2[mmAA2cut$end50 > mmAA2cut$sublength] <- mmAA2cut$sublength[mmAA2cut$end50>mmAA2cut$sublength]
mmAA2cut$cut2[mmAA2cut$end50 <= mmAA2cut$sublength] <- mmAA2cut$end50[mmAA2cut$end50 <= mmAA2cut$sublength]
mmAA2cut <- mutate(mmAA2cut, hseqFull = substring(subseq, cut1, cut2)) #contains surrounding areas
mmAA2cut <- mutate(mmAA2cut, hseq = substring(subseq, start, end)) #Only has sequence that overlaps directly with helix
head(mmAA2cut)
mmAA2cut <- mutate(mmAA2cut, hexKsFull = str_count(hseqFull, "K"))
mmAA2cut <- mutate(mmAA2cut, hexKs = str_count(hseq, "K"))
mmAA2cut <- mutate(mmAA2cut, hexRsFull = str_count(hseqFull, "R"))
mmAA2cut <- mutate(mmAA2cut, hexRs = str_count(hseq, "R"))

mmAA2cut <- mutate(mmAA2cut, hexlengthFull=str_count(hseqFull))
mmAA2cut <- mutate(mmAA2cut, hexlength = str_count(hseq))

mmAA2cut <- mutate(mmAA2cut, RKpercentFull=((hexRsFull+hexKsFull)/hexlengthFull)*100)
mmAA2cut <- mutate(mmAA2cut, RKpercent = ((hexRs+hexKs)/hexlength)*100)

CPP.hex <- filter(mmAA2cut, RKpercent > 0)
dim(mmAA2cut)
dim(CPP.hex)
head(CPP.hex)

length(unique(mmAA2cut$tId)) #4662 unique ones

#Not sure what else I want to do to filter, so I'm going to move on to shorter ones now

####################################################################################################
###7. Filtering short transcripts
####################################################################################################
mmAA <- mutate(mmAA, sublength=str_count(subseq))
mmAAshort <- filter(mmAA, sublength < 100)
dim(mmAAshort) #769, 7

countRKshort <- function(mmAA) {
  aastring <- AAString(mmAA)
  nsR <- letterFrequencyInSlidingView(aastring, 16, "RK") }

RKwindowsShort <- lapply(mmAAshort$subseq, countRKshort)
names(RKwindowsShort) <- mmAAshort$tId
mmAAshort$maxRKper16<- lapply(lapply(mmAAshort$subseq, countRKshort), max)
mmAAshort <- filter(mmAAshort, maxRKper16 >= 4)
dim(mmAAshort) #leaves me with 373 

CPPshort <- inner_join(mmAAshort, struct2ryAll, by = c("tId"="fName"))
dim(CPPshort) #105 x 14 : missing 268 of them

######IMPORTANT: THESE PEPTIDES DO NOT HAVE SECONDARY STRUCTURE PREDICTIONS###########
missingshort <- mmAAshort[!(mmAAshort$tId %in% struct2ryAll$fName),]
sum(!duplicated(missingshort$tId))
write.table(missingshort, file.path(outdir, "shortPepMissing2ndryPred.txt"), sep = "\t", quote=F)
###################################################################################


hexSplit <- function(Prediction) {
  pred.df <- strsplit(Prediction, "")
  pred.df <- data.frame(pred.df)
  names(pred.df) <- "hex"
  pred.df <- data.frame(pred.df[(pred.df$hex != ","),])
  pred.df$pos <- seq(1:nrow(pred.df))
  names(pred.df) <- c("hex", "pos")
  pred.df$hex[pred.df$hex!="H"] <- "-"
  zeros <- array(0, nrow(pred.df))
  pos.df <- data.frame(start=zeros, end=zeros)
  if (pred.df$hex[1] == "H") {
    pos.df$start[1] = 1
  }
  if(pred.df$hex[nrow(pred.df)] == "H") {
    pos.df$end[nrow(pred.df)]= nrow(pred.df)
  }
  for (i in 1:(nrow(pred.df)-1)) {
    if (pred.df$hex[i] == "-" & pred.df$hex[i+1] == "H") 
    {pos.df$start[i] <- pred.df$pos[i+1]}
    if (pred.df$hex[i] == "H" & pred.df$hex[i+1] == "-")
    {pos.df$end[i] <- pred.df$pos[i]}
  }
  pos.df <- data.frame(pos.df)
  pos.hex<-data.frame(start=pos.df$start[pos.df$start!="0"])
  
  pos.hex$end <- pos.df$end[pos.df$end != "0"]
  pos.hex
}


splitHexesShort <- lapply(CPPshort$Prediction, hexSplit)
save(splitHexesShort, file="~/Dropbox/WillseyLab/CPPs/hexsplitoutputShort.txt")
load("~/Dropbox/WillseyLab/CPPs/hexsplitoutput.txt")

names(splitHexesShort) <- CPPshort$tId
splitHexesDfshort <- do.call(rbind, splitHexesShort)
splitHexesDfshort$tId <- do.call(rbind, strsplit(rownames(splitHexesDfshort), "[.]"))[,1]
save(splitHexesDfshort, file=file.path(outdir, "splitHexesDfshort.RData"))
dim(splitHexesDfshort) #218     3

#merge info with sequences
CPPshort.pred <- data.frame(tId = CPPshort$tId, subseq = CPPshort$subseq, sublength=CPPshort$sublength)
CPPshort.pred <- full_join(splitHexesDfshort, CPPshort.pred)
dim(CPPshort.pred)
CPPshort.pred <- filter(CPPshort.pred, start != "NA")
dim(CPPshort.pred)
CPPshort.pred <- mutate(CPPshort.pred, start50 = start-50)
CPPshort.pred <- mutate(CPPshort.pred, end50 = end+50)
CPPshort.pred$cut1[CPPshort.pred$start50 <= 0] <- 1
CPPshort.pred$cut1[CPPshort.pred$start50 > 0] <- CPPshort.pred$start50[CPPshort.pred$start50 > 0]
CPPshort.pred$cut2[CPPshort.pred$end50 > CPPshort.pred$sublength] <- CPPshort.pred$sublength[CPPshort.pred$end50>CPPshort.pred$sublength]
CPPshort.pred$cut2[CPPshort.pred$end50 <= CPPshort.pred$sublength] <- CPPshort.pred$end50[CPPshort.pred$end50 <= CPPshort.pred$sublength]
CPPshort.pred <- mutate(CPPshort.pred, hseqFull = substring(subseq, cut1, cut2)) #contains surrounding areas
CPPshort.pred <- mutate(CPPshort.pred, hseq = substring(subseq, start, end)) #Only has sequence that overlaps directly with helix
head(CPPshort.pred)
CPPshort.pred <- mutate(CPPshort.pred, hexKsFull = str_count(hseqFull, "K"))
CPPshort.pred <- mutate(CPPshort.pred, hexKs = str_count(hseq, "K"))
CPPshort.pred <- mutate(CPPshort.pred, hexRsFull = str_count(hseqFull, "R"))
CPPshort.pred <- mutate(CPPshort.pred, hexRs = str_count(hseq, "R"))

CPPshort.pred <- mutate(CPPshort.pred, hexlengthFull=str_count(hseqFull))
CPPshort.pred <- mutate(CPPshort.pred, hexlength = str_count(hseq))

CPPshort.pred <- mutate(CPPshort.pred, RKpercentFull=((hexRsFull+hexKsFull)/hexlengthFull)*100)
CPPshort.pred <- mutate(CPPshort.pred, RKpercent = ((hexRs+hexKs)/hexlength)*100)

CPPshort.hex <- filter(CPPshort.pred, RKpercent > 0)
dim(CPPshort.pred)
dim(CPPshort.hex) #112 remaining
sum(!duplicated(CPPshort.hex$tId)) #69 unique
head(CPPshort.hex)
#CPPshort.hex <- filter(CPPshort.pred, RKpercent > 25)
dim(CPPshort.hex) #35 remaining
sum(!duplicated(CPPshort.hex$tId)) #30 unique

##One possible change is to reduce number of RKs required on helix, just require them in surrounding area, for now, am not filtering based on this, will filter 
#when all reads are recombined

########################################################################
#Next: rejoin with longer reads
#mmAA2cut
#CPPshort.hex

#renaming for simplicity
longPeps <- filter(mmAA2cut, (hexRs+hexKs) > 0)
longPeps$RKpercent[is.na(longPeps$RKpercent)==T] <- 0 #0/0 giving some NAs - removed and made 0
longPeps$RKpercentFull[is.na(longPeps$RKpercent)==T] <- 0

shortPeps <- CPPshort.hex

colnames(longPeps)
colnames(shortPeps)


longPeps <- data.frame(tId=longPeps$tId, subseq=longPeps$subseq, sublength=longPeps$sublength, start=longPeps$start, end=longPeps$end,
                       start50=longPeps$start50, end50=longPeps$end50, cut1=longPeps$cut1, cut2=longPeps$cut2, hexlength=longPeps$hexlength, 
                       hexlengthFull=longPeps$hexlengthFull, hseq=longPeps$hseq, hseqFull=longPeps$hseqFull, hexRs=longPeps$hexRs, 
                       hexRsFull=longPeps$hexRsFull, hexKs = longPeps$hexKs, hexKsfull = longPeps$hexKsFull, RKpercent=longPeps$RKpercent, 
                       RKpercentFull=longPeps$RKpercentFull)
shortPeps <- data.frame(tId=shortPeps$tId, subseq=shortPeps$subseq, sublength=shortPeps$sublength, start=shortPeps$start, end=shortPeps$end,
                       start50=shortPeps$start50, end50=shortPeps$end50, cut1=shortPeps$cut1, cut2=shortPeps$cut2, hexlength=shortPeps$hexlength, 
                       hexlengthFull=shortPeps$hexlengthFull, hseq=shortPeps$hseq, hseqFull=shortPeps$hseqFull, hexRs=shortPeps$hexRs, 
                       hexRsFull=shortPeps$hexRsFull, hexKs = shortPeps$hexKs, hexKsfull = shortPeps$hexKsFull, RKpercent=shortPeps$RKpercent, 
                       RKpercentFull=shortPeps$RKpercentFull)



allPeps <- rbind(longPeps, shortPeps)
dim(allPeps)
min(allPeps$RKpercent) #1.5625
min(allPeps$RKpercentFull) #1.162791

allPeps <- filter(allPeps, RKpercentFull >= 25)
dim(allPeps) #leaves me with 43
min(allPeps$RKpercent) #9.090909
allPeps
allPeps1 <- allPeps
sum(!duplicated(allPeps$tId)) #35 are unique
allPeps <- filter(allPeps, RKpercent >= 25)
dim(allPeps) #gives me 26 at 25%
allPeps2 <- allPeps
##############Enrichment analysis
# ###check enrichment
datdir <- "~/Dropbox/WillseyLab/CPPs"
mouseCPPs <- read.delim(file.path(datdir, "mouseCPPs.txt"))
head(mouseCPPs)

#need to merge gene names with tIds
tId2gId <- read.delim("~/Dropbox/WillseyLab/CPPs/MmAAtable.txt")
tId2gId <- tId2gId[,c(1:2)]
allPepsgId1 <- left_join(allPeps1, tId2gId)
allPepsgId2 <- left_join(allPeps2, tId2gId)

mouseCPPs <- mouseCPPs[!duplicated(mouseCPPs$mouseEnsembl),] #here, removing all duplicate gIds
#there are duplicate gIds because some peptides have been found to be CPPs in a) more than 1 experiment, or
#b) there are multiple regions that act as CPPs. For now, to simplify test, am removing all duplicates from both and just
#searching on the gId level for genes that seem to have CPP properties.... can easily change this if need be
#to change this, get rid of this step, and then be sure to keep duplicate gIds below

sum(mouseCPPs$mouseEnsembl %in% tId2gId$gId)/nrow(tId2gId) #.0002929687
sum(mouseCPPs$mouseEnsembl %in% allPepsgId1$gId)/nrow(allPepsgId1) #.02325581 (don't keep any more this way)
sum(mouseCPPs$mouseEnsembl %in% allPepsgId2$gId)/nrow(allPepsgId2) #.03846154
sum(mouseCPPs$mouseEnsembl %in% allPepsgId1$gId) #1
sum(mouseCPPs$mouseEnsembl %in% allPepsgId2$gId) #1

# /nrow(macExpEns)
# sum(mouseCPPs$mouseEnsembl %in% filtMacExp$gId)/nrow(filtMacExp)
# filtMacExp[filtMacExp$gId %in% mouseCPPs$mouseEnsembl,]


#This isn't the issue
sum(missing$gId %in% mouseCPPs$mouseEnsembl)
sum(missingshort$gId %in% mouseCPPs$mouseEnsembl)
sum(mouseCPPs$mouseEnsembl %in% mmAA$gId)


#####Enrichment analysis
q <- 0 #number of hits in my list -1
m <- 18 #number of known CPPs: nrow(mouseCPPs)
n <- nrow(mmAA_original) - nrow(mouseCPPs) #total - CPPs
k <- 26 #length of my list; nrow(allPeps1 or allPeps2)


phyper(q, m, n, k)
dhyper(q, m, n, k)



#####
#Test
mmAA2 <- mutate(mmAA_original, length = str_count(seq)) #this creates a column with peptide length
mmAA2 <- separate(mmAA2, seq, c("subseq"), remove=F, sep = "[*]", extra="drop")
head(mmAA2)
mmAA2 <- mutate(mmAA2, sublength = str_count(subseq))
mmAA2 <- mutate(mmAA2, Rc = str_count(subseq, "R")) #number of Rs in the entire peptide
mmAA2 <- mutate(mmAA2, Kc = str_count(subseq, "K"))
mmAA2 <- mutate(mmAA2, totalRp = Rc/length * 100) #percentage of Rs in the peptide
mmAA2 <- mutate(mmAA2, totalKp = Kc/length *100)
mmAA2 <- mutate(mmAA2, totalRKc = Rc + Kc)
head(mmAA2)
dim(mmAA2) #61440, 10
mmAA2 <- filter(mmAA2, totalRp+totalKp>0)
dim(mmAA2) #61162, 10
#This function will count the number of Args in a window of defined size. To make it simpler to use in an apply function, 
#to change the size of the window, just alter it in the function below (rather than setting as a variable).
counteR <- function(mmAA) {
  aastring <- AAString(mmAA)
  nsR <- letterFrequencyInSlidingView(aastring, 32, "RK") }



mmAA2$maxRKper32 <- 0
mmAA2$maxRKper32[mmAA2$length < 32] <- mmAA2$totalRKc[mmAA2$length < 32]
mmAA2$maxRper32[mmAA2$length >= 32] <- lapply(lapply(mmAA2$seq[mmAA2$length >= 32], counteR), max)
mmAA2$maxRKper32perc <- 0
mmAA2$maxRKper32perc[mmAA2$length < 32] <- mmAA2$totalRKc[mmAA2$length < 32]/mmAA2$sublength







#with origina mmAA, with no changes made to it
head(mmAA)
AA <- separate(mmAA, seq, c("subseq"), remove=F, sep = "[*]", extra="drop")
AA <- dplyr::select(AA, -seq)
AA <- dplyr::rename(AA, seq=subseq)
AA <- dplyr::mutate(AA, length = str_count(seq))
#AA <- filter(AA, length <= 32)
AA <- dplyr::mutate(AA, Rs = str_count(seq, "R"))
AA <- dplyr::mutate(AA, Ks = str_count(seq, "K"))
AA <- filter(AA, (Rs+Ks) >= 1)
dim(AA) #61162

countRKs <- function(seq, name) {
  #seqInd <- match("seq", colnames(mmAA))
  #seq <- mmAA[,seqInd]
  #nameInd <- match("tId", colnames(mmAA))
  #name <- mmAA[,nameInd]
  n <- str_count(seq)
  l <- 32
  seqtable <- data.frame(substring(seq, 1:(n-l+1), l:n))
  names(seqtable) <- c("seq")
  AAseqlist <- lapply(seqtable$seq, AAString)
  seqtable$pos <- rownames(seqtable)
  seqtable$tId <- name
  seqtable <- dplyr::mutate(seqtable, rnames = paste(tId,".",pos, sep=""))
  rownames(seqtable) <- seqtable$rnames
  nsR <- lapply(seqtable$seq, str_count, "R|K") 
  seqtable$maxRK <- data.frame(maxRK=unlist(lapply(nsR, max))) #may work if i get rid of the data.frame here
  #seqtable <- seqtable[,c(1,2,3,5)]
  seqtable
}

AAshort <- mutate(AA[AA$length<32,], RKper32=Rs+Ks)
dim(AAshort) #755

AA <- filter(AA, length>=32)
dim(AA) #60407

RKs <- mapply(countRKs, AA$seq, AA$tId)





names(RKs) <- AA$tId
RKmaxes <- lapply(1:length(names(RKs)), function(x) max(RKs[[x]]))


RKs <- RKwindow[RKmaxes>=4]
RKsdf <- do.call(rbind, RKs)
length(names(RKwindow2)) #trims to 6295
length(names(RKwindow)) #from 7002

RKwindow2 <- lapply(RKwindow2, function(x) subset(x, maxRK >= 4))
RKwindow2 <- lapply(RKwindow2, function(x) mutate(x, maxRK = unlist(as.list(maxRK))))
#RKtest <- lapply(RKtest, function(x) mutate(x, maxRK=unlist(as.list(maxRK))))

#row.names(RKwindow2) <- names(RKwindow2)
RKdf <- rbind.fill(RKwindow2)
#dftest <- rbind.fill(RKtest)
head(RKdf)
names(RKdf) <- c("seq100", "pos", "tId", "maxRK")
dim(RKdf) #1675892       4

min(RKdf$maxRK) #checking that filter worked; it did
sum(duplicated(RKdf$tId))
sum(!duplicated(RKdf$tId)) #6295





