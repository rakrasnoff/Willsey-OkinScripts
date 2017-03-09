require(Biostrings)
options(stringsAsFactors = FALSE)
library(stringr)
library(dplyr)
library(plyr)
library(zoo)
library(seqRFLP)
library(org.Mm.eg.db)

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
dim(mmAA)


####################################################################################################
### Remove any sequence that occurs after a stop
####################################################################################################
#
mmAA <- mutate(mmAA, length = str_count(mmAA$seq)) #this creates a column with peptide length
#seqsplit <- do.call(rbind, lapply(mmAA$seq, function(x) str_split(x, "[*]")))[,1]
mmAA <- mutate(mmAA, seq2 = as.list(str_split(seq,"[*]"))) 
mmAA <- mutate(mmAA, seq2 = gsub("[,]", , seq2))
head(mmAA)
tail(mmAA)

#This function creates a dataframe with each piece of AA sequence, separated by stop codons
stopSplit <- function(seq, tId) {
  splits <- data.frame(seq=strsplit(seq, "[*]"))
  splits$tId <- tId
  names(splits) <- c("seq", "tId")
  splits[1,]
}

mmAAstop <- mapply(stopSplit, mmAA$seq, mmAA$tId, SIMPLIFY=FALSE)
names(mmAAstop) <- mmAA$tId
mmAAstopdf <- rbind.fill(mmAAstop)
dim(mmAAstopdf) #60343 rows

#This segment is for looking at starts - don't worry about it now#
#Now, filter by Met presence; must have at least 1 met to keep
#mmAAstopdf <- mutate(mmAAstopdf, Mc = str_count(mmAAstopdf$seq, "M"))
#mmAAstopdf <- filter(mmAAstopdf, Mc >= 1)
#dim(mmAAstopdf) #This gets rid of a few hundred, bringing me to 59104
#Now, need to cut from first M in sequence to last AA of sequence
#pos = regexpr('pattern', x) # Returns position of 1st match in a string
#mmAAstopdf <- mutate(mmAAstopdf, Mpos = regexpr("M", seq))
#sum(mmAAstopdf$Mpos > 1) #M does not start the sequence for 5535 sequences
#mmAAstopdf <- mutate(mmAAstopdf, length = str_count(seq)) #this gives me the last position, to make it easy to cut
#mmAAstopdf <- mutate(mmAAstopdf, Mseq = substring(seq, Mpos, length))
#Internal check: the length should differ for 5535 sequences
#mmAAstopdf <- mutate(mmAAstopdf, length2 = str_count(Mseq))
#sum(mmAAstopdf$length != mmAAstopdf$length2) #Equals 5535, seems to have worked properly
#mmAAstopdf <- mutate(mmAAstopdf, Mpos2 = regexpr("M", Mseq)) #make sure it cut correctly
#sum(mmAAstopdf$Mpos2 > 1) #gives 0 - definitely cut correctly

#I will now filter by SRP in 2 ways
#head(mmAAstopdf)

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
mmAAstopdf <- left_join(mmAAstopdf, SPdf, by=c("tId"= "name"))
dim(mmAAstopdf) 
head(mmAAstopdf)
names(mmAAstopdf) <- c("seq", "tId", "Mc", "Mpos", "length", "Mseq", "length2", "Mpos2", "sigPep")

mmSP <- filter(mmAAstopdf, sigPep != "N")
dim(mmSP) #gives me 7642
#Now, each of these 7642 sequences needs to be run through the prediction server
#mmAA <- mmSP

#Filtering by SRP round 2
mmAA4Fasta <- data.frame(tId = make.unique(mmSP$tId, sep="_"), Mseq = mmSP$Mseq) #using make unique here, because otherwise server only runs once
#fast1 <- dataframe2fas(mmAA4Fasta[c(1:2000),], file = "~/Dropbox/WillseyLab/CPPs/SRPpred2/input/AAs1.fa")
#fast2 <- dataframe2fas(mmAA4Fasta[c(2001:4000),], file = "~/Dropbox/WillseyLab/CPPs/SRPpred2/input/AAs2.fa")
#fast3 <- dataframe2fas(mmAA4Fasta[c(4001:6000),], file = "~/Dropbox/WillseyLab/CPPs/SRPpred2/input/AAs3.fa")
#fast4 <- dataframe2fas(mmAA4Fasta[c(6001:7642),], file = "~/Dropbox/WillseyLab/CPPs/SRPpred2/input/AAs4.fa")

SPdatdir2 <- "~/Dropbox/WillseyLab/CPPs/SRPpred2/predictions"
files <- dir(SPdatdir2)

#Read in files
SPs <- lapply(files, function(x) { 
  SPlist <- read.delim(file.path("~/Dropbox/WillseyLab/CPPs/SRPpred2/predictions", x))
})
head(SPs[[1]])
head(SPs[[2]])

#Rename, reformat
for (i in 2:4){
  names(SPs[[i]]) <- names(SPs[[1]])
  SPs[[i]] <- SPs[[i]][,c(1:12)]
}
summary(SPs)

SPdf <- rbind.fill(SPs)
dim(SPdf) #correct number of rows
SPdf <- SPdf[,c(1,10)]
head(SPdf)

sum(duplicated(SPdf$name))

mmSP <- left_join(mmSP, SPdf, by=c("tId"= "name"))
dim(mmSP) 
head(mmSP)
sum(duplicated(mmSP$tId))
names(mmSP) <- c("seq", "tId", "Mc", "Mpos", "length", "Mseq", "length2", "Mpos2", "sigPep", "sigPep2")

#dummy <- filter(mmSP, sigPep2 !="Y")
mmSP <- filter(mmSP, sigPep2 != "N")

dim(mmSP) #gives me 7540 (doesn't change by much)
mmAA <- mmSP

####################################################################################################
###4. Filter by R and K presence
####################################################################################################
#mmAA <- mutate(mmAA, Mlength = str_count(mmAA$Mseq)) #this creates a column with peptide length
mmAA <- mutate(mmAA, Rc = str_count(mmAA$Mseq, "R")) #number of Rs in the entire peptide
mmAA <- mutate(mmAA, Kc = str_count(mmAA$Mseq, "K"))
mmAA <- mutate(mmAA, Rp = Rc/length * 100) #percentage of Rs in the peptide
mmAA <- mutate(mmAA, Kp = Kc/length *100)
head(mmAA)
dim(mmAA)

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

mmAAlong <- mmAA[mmAA$length2 >= 100,] #6916 rows
RKwindow <- mapply(countRKs, mmAAlong$Mseq, mmAAlong$tId, SIMPLIFY=FALSE)
#save(RKwindow, file=file.path(outdir, "RKwindowOutput.RData"))
#load(file.path(outdir, "RKwindowOutput.RData"))
names(RKwindow) <- mmAAlong$tId
sum(duplicated(mmAAlong$tId))
mmAAlong$tId[duplicated(mmAAlong$tId)]
maxes <- lapply(1:6916, function(x) max(RKwindow[[x]]$maxRK))
RKwindow2 <- RKwindow[maxes>=4]
length(names(RKwindow2)) #trims to 6228
length(names(RKwindow)) #from 6916
RKwindow2 <- lapply(RKwindow2, function(x) subset(x, maxRK >= 4))
RKwindow2 <- lapply(RKwindow2, function(x) mutate(x, maxRK = unlist(as.list(maxRK))))

#row.names(RKwindow2) <- names(RKwindow2)
RKdf <- rbind.fill(RKwindow2)
dim(RKdf) #1668297       4

min(RKdf$maxRK) #checking that filter worked; it did
sum(duplicated(RKdf$tId))
sum(!duplicated(RKdf$tId)) #6224 #why isn't this exactly 6228, or more different?

####################################################################################################
###5. Add alpha helix content
####################################################################################################

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
sum(struct2ry$fName %in% struct2ry2$Name)
struct2ry2 <- struct2ry2[!struct2ry2$Name %in% struct2ry$fName,]
names(struct2ry2) <- c("fName", "Prediction", "Confidence")
struct2ryAll <- rbind(struct2ry, struct2ry2) #have 3885 rows it seems - is this how many I had before?


#combine with mmAA dataframe
CPPdf <- left_join(RKdf, struct2ryAll, by = c("tId"="fName"))
CPPdf <- left_join(CPPdf, mmAAlong, by=c("tId"="tId"))
dim(CPPdf) #1669988      19
head(CPPdf)

head(CPPdf, 100)
#However, I have NAs for many sequences, because we are filtering differently than before. I need to get a list of IDs that I need predictions for. Even better, list of Ids with subsequences, and my rownames
emptySeq <- CPPdf[is.na(CPPdf$Prediction)==T,]
dim(emptySeq) #it's over a million sequences...
length(unique(emptySeq$tId)) #but only 4423 peptides #now down to 4371 peptides, which is good news - I did cut out a few by splitting things up 
#Even though it's not enough to make a difference, it is reassuring that the effect was in the same direction as predicted
head(emptySeq)
#I need the original peptide lengths - will need to cut and merge with mmAAlong
max(emptySeq$length)
#Need to write these sequences to fasta
empty4Fasta <- data.frame(tId=emptySeq$tId, seqAll=emptySeq$seqAll)
empty4Fasta <- empty4Fasta[!duplicated(empty4Fasta$tId), ]
dim(empty4Fasta) #4371 long, the correct length
#fast1 <- dataframe2fas(empty4Fasta, file = "~/Dropbox/WillseyLab/CPPs/RKEnrichedPeptides.fa") commented out so I don't overwrite my original list (for comparison purposes)...
fast2 <- dataframe2fas(empty4Fasta, file = "~/Dropbox/WillseyLab/CPPs/RKEnrichedPeptides2.fa")
#Are all 4371 within the 4423 I previously identified?Need to figure out how many Ids overlap, but also how many times sequence length for Ids overlap

testFast <- readAAStringSet("~/Dropbox/WillseyLab/CPPs/RKEnrichedPeptides.fa") 
namesTestFast <- paste(names(testFast))
length(namesTestFast)
namesTF <- as.data.frame(namesTestFast)
class(namesTF[,1])
head(namesTF)
tIdTF <- do.call("rbind", strsplit(namesTF[,1], "[. :]"))[,1]
#gIdTF <- do.call("rbind", strsplit(namesTF[,1], "[. :]"))[,11]
length(tIdTF)
seqTF = as.data.frame(paste(testFast))
dim(seqTF)
TF2 <- data.frame(tIdTF, seqTF)
names(TF2) <- c("tId", "seq")
head(TF2)

sum(unique(emptySeq$tId) %in% TF2$tId) #they are all there (all 4371 of them, so that's good)

####################################################################################################
###6. Filter by alpha helixes
####################################################################################################

head(CPPdf)
CPPdf.C <- CPPdf[is.na(CPPdf$Prediction)==F,]
dim(CPPdf)
dim(CPPdf.C)


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

#trialList <- lapply(1:nrow(CPPdf.C), function(x) hexSplit(CPPdf.C$Prediction[x]))
trialList <- lapply(CPPdf.C$Prediction, hexSplit)
names(trialList) <- mmAAshort$tId[1:nrow(mmAAshort)]
trialdf <- do.call(rbind, trialList)
trialdf$tId <- do.call(rbind, strsplit(rownames(trialdf), "[.]"))[,1]
head(trialdf)

#merge info with sequences
mmAA2merge <- data.frame(tId = mmAAshort$tId, seq = mmAAshort$seq)
mmAA2cut <- full_join(trialdf, mmAA2merge)
dim(mmAA2cut)
dim(trialdf)
mmAA2cut <- mutate(mmAA2cut, hseq = substring(seq, start, end))
head(mmAA2cut)
mmAA2cut <- mutate(mmAA2cut, Rs = str_count(hseq, "R"))
mmAA2cut <- mutate(mmAA2cut, length = str_count(hseq))
mmAA2cut <- mutate(mmAA2cut, percent = (Rs/length)*100)
mmAA2cut <- filter(mmAA2cut, percent > 0)
dim(mmAA2cut)
length(unique(mmAA2cut$tId)) 

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
dim(mmAA2cut)
length(unique(mmAA2cut$tId)) #1259
save(mmAA2cut, file=file.path(outdir, "CPPpredictionsHelixSlidingWindowOf8.RData"))

tIds <- unique(mmAA2cut$tId) 





























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



