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
#load(file.path(outdir, "mmSP.RData"))
mmAA <- mmSP
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

mmAAshort <- filter(mmAA, sublength < 100)
dim(mmAAshort)

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
missingshort <- CPP[!(CPPshort$tId %in% struct2ryAll$fName),]
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
CPPshort.hex <- filter(CPPshort.pred, RKpercent > 25)
dim(CPPshort.hex) #35 remaining
sum(!duplicated(CPPshort.hex$tId)) #30 unique

########################################################################
#Next: rejoin with longer reads

























#filter to remove unnecessary rows, make df smaller
mmAAR <- mmAA[mmAA$maxRper15 > 0, ]
mmAAR$keep <- 0
mmAAR$keep[mmAAR$length < 15] <- 1
mmAAR$keep[mmAAR$length >= 15 & mmAAR$maxRper15 >= 4] <- 1
mmAAR <- subset(mmAAR, keep == 1)
dim(mmAAR)

mmAAlong <- mmAAR[mmAAR$length >= 15,]
Rwindows <- lapply(1:29629, function(x) countRswithNames(mmAAlong[x,]))
RwindowsList <- Rwindows
Rwindows <- lapply(Rwindows, function(x) mutate(x, pos.num=seq(1:nrow(x))))
Rwindows <- lapply(Rwindows, subset, R>=4)
#stopped here; need to add pos.num thing to short ones before combining
Rwindows <- do.call(rbind, Rwindows)


#make dataframe with only tIds aand seqs, so I can add sequences here
seqs <- data.frame(tId=mmAAlong$tId, seq=mmAAlong$seq)
names(Rwindows) <- c("R", "tId")

Rdflong <- left_join(Rwindows, seqs)


dfshort <- mmAA[mmAA$length < 15,]
dfshort <- filter(dfshort, Rc > 0)
Rdfshort <- data.frame(R=dfshort$Rc, tId=dfshort$tId, seq=dfshort$seq)
Rdfshort$pos.num = 1 #need to make sure columns are in the correct order
dftest <- rbind(Rdfshort, Rdflong)



######



mmAA1 <- mmAA
#################################################################################################
#################################################################################################


###################################### MACROPHAGE EXP DATA ######################################
#################################################################################################
#Add macrophage expression values
macExpDat <- read.delim(file.path(MEdatdir, "macExpEnsembl.txt"))
#now check expression
macExpList <- list(macExpDat[,c(4:13)])
summary(unlist(macExpList)) #getting quartiles
head(macExpList[[1]])
## I will take these maxes and add them to my dataframe. I am only taking maxes because as long as the protein
## is expressed at a minimum of 110 in one sample, it will stay for further analysis. However, this data can be revisited
## if we decide that particular tissues are most relevant, or decide to filter based on average, etc.
macExpDat$max <- apply(macExpDat[,c(4:13)], 1, max)
maxMacExp <- macExpDat[,c(14:15)]
names(maxMacExp) <- c("ensId", "macExp")
mmAA <- left_join(mmAA, maxMacExp, by=c("gId" = "ensId"))
head(mmAA)
tail(mmAA)
mmAA <- mmAA[!duplicated(mmAA$seq),]
head(mmAA)
dim(mmAA)

###################################### SIGNAL PEPTIDE DATA ######################################
#################################################################################################


###################################### ALPHA HELIXES ############################################
#################################################################################################
#Write to file for running through prediction server
AAFasta <- read.delim(file.path(datdir, "filteredListPossibleCPPs.txt"))
fastaDF <- AAFasta[,c(1,3)]
#fastaDF$seq <- gsub("[*]","",fastaDF$seq) #now fixed up above
fast1 <- dataframe2fas(fastaDF, file = "~/Dropbox/WillseyLab/CPPs/filteredFAs/AAs1.fa")


#for rerunning where stars messed stuff up
#test <- data.frame(mmAASS$tId, mmAASS$seq, mmAASS$length, mmAASS$Prediction)
#names(test) <- c("tId", "seq", "length", "Prediction")
##test$seq <- gsub("$[*]", "", test$seq)
#test$predLength <- data.frame(str_count(test[,4]))
#test <- mutate(test, lengthLessOne = length - 1)
#test$eq <- array(0, 2938)
#head(test)
#test$eq[is.na(test$predLength)] <- 0
#test$eq[test$length == test$predLength] <- 0
#test$eq[test$length != test$predLength] <- 1
#test$eq[test$lengthLessOne == test$predLength] <- 0
#test2 <- test[test$eq == 1,]
#test2$seq <- gsub("[*]", "", test2$seq)
#test2Fas <- test2[,c(1,2)]
#fasta2 <- dataframe2fas(test2Fas, file="~/Dropbox/WillseyLab/CPPs/filteredFAs/AAsrerun.fa")


#read in secondary structure results
struct2ry <- read.delim("~/Dropbox/WillseyLab/CPPs/2ndryStruct/merged_jpred_results_2-1-17.txt")
struct2ry2 <- read.delim("~/Dropbox/WillseyLab/CPPs/2ndryStruct/merged_jpred_results_batch2_2-16-17.txt")
head(struct2ry)
head(struct2ry2)
class(struct2ry$Prediction)
struct2ry$Prediction <- gsub(",","", struct2ry$Prediction)
struct2ry$Confidence <- gsub(",","", struct2ry$Confidence)
struct2ry2$Prediction <- gsub(",","", struct2ry2$Prediction)
struct2ry2$Confidence <- gsub(",","", struct2ry2$Confidence)
names(struct2ry)
struct2ry <- struct2ry[!struct2ry$fName %in% struct2ry2$Name,]
names(struct2ry2) <- c("fName", "Prediction", "Confidence")
struct2ryAll <- rbind(struct2ry, struct2ry2)


#combine with mmAA dataframe
mmAA2ndry <- left_join(mmAA, struct2ryAll, by = c("tId"="fName"))


#################################################################################################
#################################################################################################
#Save dataframe
mmAAtoSave <- mmAA2ndry
mmAAtoSave$maxRper15 <- do.call(rbind, mmAAtoSave$maxRper15)
mmAAtoSave$maxRper15 <- data.frame(mmAAtoSave$maxRper15)
names(mmAAtoSave) <- c("tId", "gId", "seq", "length", "Rc", "Rp", "maxRper15", "macExp", "sigPep", "Prediction", "Confidence")
#write.table(mmAAtoSave, file.path(outdir, "finalAADataframe.txt"), sep = "\t", quote=F, row.name=F)
save(mmAAtoSave, file=file.path(outdir, "mmAAFinal.RData"))



