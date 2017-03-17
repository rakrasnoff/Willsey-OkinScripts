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
RKwindow <- mapply(countRKs, mmAAlong$subseq, mmAAlong$tId, SIMPLIFY=FALSE)
#RKtest <- mapply(countRKs, mmAAlong$subseq[1:3], mmAAlong$tId[1:3], SIMPLIFY=FALSE)
#save(RKwindow, file=file.path(outdir, "RKwindowOutput.RData"))
#load(file.path(outdir, "RKwindowOutput.RData"))
names(RKwindow) <- mmAAlong$tId
names(RKtest) <- mmAAlong$tId[1:3]
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
CPP <- data.frame(tId=CPPdf$tId, gId=CPPdf$gId, fullseq=CPPdf$seq, fulllength=CPPdf$length, subseq=CPPdf$subseq, sublength=CPPdf$sublength, 
                  sigPep=CPPdf$sigPep, Rc=CPPdf$Rc, Kc=CPPdf$Kc, seq100=CPPdf$seq100, 
                  pos=CPPdf$pos, maxRK=CPPdf$maxRK, Prediction=CPPdf$Prediction, Confidence=CPPdf$Confidence)

head(CPP, 100)
####################################################################################################
###6. Filter by alpha helixes
####################################################################################################

head(CPPdf)
CPP.C <- CPP[is.na(CPP$Prediction)==F,]
dim(CPP)
dim(CPP.C)

#First, trim predictions to correct windows of 100
class(CPP.C$pos) <- "numeric"
CPP.C <- mutate(CPP.C, pred100=(substring(Prediction, pos, pos+99)))
head(CPP.C)


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
trialList <- lapply(CPP.C$pred100, hexSplit)
#needs updates after this point
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



