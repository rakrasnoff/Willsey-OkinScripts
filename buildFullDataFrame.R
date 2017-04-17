require(Biostrings)
options(stringsAsFactors = FALSE)
library(stringr)
library(dplyr)
library(plyr)
library(zoo)
library(seqRFLP)
library(org.Mm.eg.db)

##Create dataframe with all information necessary for filtering

#Start by loading dataframe with mouse peptides. This list comes from GEO database. List initially 
#loaded in CPP_load_DNA script
mmAA <- read.delim("~/Dropbox/WillseyLab/CPPs/MmAAtable.txt")

#Now set directories
SPdatdir <- "~/Dropbox/WillseyLab/CPPs/sigPepFiles" #signal peptide files
MEdatdir <- "~/Dropbox/WillseyLab/CPPs/MacExpDat" #macrophage expression data
outdir <- "~/Dropbox/WillseyLab/CPPs"

#Adding R content 
mmAA$seq <- gsub("[*]", "", mmAA$seq) #added this for second run; will fix issues with length and astrixes
mmAA <- mutate(mmAA, length = str_count(mmAA$seq)) #this creates a column with peptide length
mmAA <- mutate(mmAA, Rc = str_count(mmAA$seq, "R")) #number of Rs in the entire peptide
mmAA <- mutate(mmAA, Rp = Rc/length * 100) #percentage of Rs in the peptide
head(mmAA)

#This function will count the number of Args in a window of defined size. To make it simpler to use in an apply function, 
#to change the size of the window, just alter it in the function below (rather than setting as a variable).
countRs <- function(mmAA) {
  aastring <- AAString(mmAA)
  nsR <- letterFrequencyInSlidingView(aastring, 15, "R") }

mmAA$maxRper15 <- 0
mmAA$maxRper15[mmAA$length < 15] <- mmAA$Rc[mmAA$length < 15]
mmAA$maxRper15[mmAA$length >= 15] <- lapply(lapply(mmAA$seq[mmAA$length >= 15], countRs), max)


#RrichWindows <- RrichWindows[lapply(RrichWindows, max) >= 4]




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
mmAA <- left_join(mmAA, SPdf, by=c("tId"= "name"))

dim(mmAA)
head(mmAA)
names(mmAA) <- c("tId", "gId", "seq", "length", "Rc", "Rp", "maxRper15", "macExp", "SigPep")

###################################### ALPHA HELIXES ############################################
#################################################################################################
#Write to file for running through prediction server
#AAFasta <- read.delim(file.path(datdir, "filteredListPossibleCPPs.txt"))
#fastaDF <- AAFasta[,c(1,3)]
#fastaDF$seq <- gsub("[*]","",fastaDF$seq) #now fixed up above
#fast1 <- dataframe2fas(fastaDF, file = "~/Dropbox/WillseyLab/CPPs/filteredFAs/AAs1.fa")


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


