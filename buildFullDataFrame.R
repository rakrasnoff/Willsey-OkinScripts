require(Biostrings)
options(stringsAsFactors = FALSE)
library(stringr)
library(dplyr)
library(plyr)
library(zoo)
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
mmAA$maxRper15 <- do.call(rbind, mmAA$maxRper15)
dim(mmAA)
head(mmAA)
names(mmAA) <- c("tId", "gId", "seq", "length", "Rc", "Rp", "maxRper15", "macExp", "SigPep")

###################################### ALPHA HELIXES ############################################
#################################################################################################


#################################################################################################
#################################################################################################
#Save dataframe
write.table(mmAA, file.path(outdir, "fullAADataFrame.txt"), sep = "\t", quote=F)



