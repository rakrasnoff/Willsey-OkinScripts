library(seqRFLP)

#Filter dataframe and check enrichment
datdir <- "~/Dropbox/WillseyLab/CPPs"
outdir <- "~/Dropbox/WillseyLab/CPPs/output"

mmAA <- read.delim(file.path(datdir, "finalAADataframe.txt"))
load(file.path(outdir, "mmAAFinal.RData"))
mmAA <- mmAAtoSave
head(mmAA)
names(mmAA)

#Filter by R content (sliding window)

mmAA$maxRper15 <- as.list(mmAA$maxRper15)
head(mmAA)
mmAAR <- mmAA[mmAA$maxRper15 > 0, ]
#mmAAR$length <- as.numeric(mmAAR$length)
dim(mmAAR)
head(mmAAR)
mmAAR$keep <- 0
mmAAR$keep[mmAAR$length < 15] <- 1
mmAAR$keep[mmAAR$length >= 15 & mmAAR$maxRper15 >= 4] <- 1
mmAAR <- subset(mmAAR, keep == 1)
dim(mmAAR)
head(mmAAR)

#Filter by Mac Exp
mmAAME <- mmAA #to filter original list
mmAAME <- mmAAR # to continue filtering previously filtered list
mmAAME$keep <- 0
mmAAME$keep[is.na(mmAAME$macExp)] <- 1
mmAAME$keep[!is.na(mmAAME$macExp) & mmAAME$macExp >= 122.5] <- 1
sum(mmAAME$keep)
#20677 (keeping those with NAs)

#Filter by Signal Peptides
mmAASP <- mmAAME
mmAASP <- subset(mmAASP, SigPep == "Y")
dim(mmAASP)
#2938



#Filter by Secondary Structure
head(mmAASP)
mmAASS <- mmAASP
mmAASS$helix <- 0
maxRper15 <- data.frame(mmAASS$tId, mmAASS$maxRper15)
mmAASS <- mmAASS[,-7]
mmAASS <- mutate(mmAASS, helix = str_count(mmAASS$Prediction, "H"))
head(mmAASS)

mmAASS <- mmAASS[mmAASS$helix != 0,]
dim(mmAASS) #this brings me to 2862; however, 972 are still NA

#peptide: use mmAASS
#argHelix <- function(peptide) {
 # pred <- strsplit(mmAASS$Prediction[mmAASS$tId == peptide], "")
  #conf <- strsplit(mmAASS$Confidence[mmAASS$tId == peptide], "")
#  seq <- strsplit(mmAASS$seq[mmAASS$tId == peptide], "")
 # argHelixDF <- data.frame(pred, conf, seq)
 # names(argHelixDF) <- c("pred", "conf", "seq")
#  argHelixDF
#}


#argHelixTest <- function(peptide) {
  #struct <- peptide$pred[peptide$seq == "R"]
  #argPerHelix <- sum(str_count(struct, "H"))
#}


#argHelix(mmAASS$tId[1])

#mmAASS$seq[9] <- substring(mmAASS$seq[9], 1, 184) #for whatever reason, didn't rerun this one
#mmAASS$seq <- sub("[*]", "", mmAASS$seq)
#mmAASS$seq <- sub("$[*]", "", mmAASS$seq)
#mmAASS$seq[18] <- sub("[*]", "", mmAASS$seq[18]) #seems like I have to run it once for each star I am removing...
#mmAASS <- mutate(mmAASS, length = str_count(mmAASS$seq))
#head(mmAASS)
#mmAASS <- mutate(mmAASS, predLength = str_count(mmAASS$Prediction))
#head(mmAASS, 20)
#mmAASS$length == mmAASS$predLength
mmAASS <- mutate(mmAASS, seq2=substring(seq, 1, predLength))

hexPos <- function(pred, seq) {
pred <- strsplit(pred, "")
#conf <- strsplit(mmAASS$Confidence[1], "")
seq <- strsplit(seq,"")
argHelixDF <- data.frame(pred, seq)
names(argHelixDF) <- c("pred", "seq")
argHelixDF$pos <- seq(1:nrow(argHelixDF))
argHelixDF$pos[argHelixDF$pred == "H"]
#argHelix <- subset(argHelixDF, pred == "H")
#helixLength <- nrow(argHelix)
#sum(argHelix$seq == "R")
}


mmAAshort <- mmAASS[1:2862,]
mmAAshort <- mutate(mmAAshort, hexPost = lapply(1:2862, function(x) hexPos(mmAAshort$Prediction[x], mmAAshort$seq2[x])))
head(mmAAshort)
#seq <- mmAAshort$seq[1]
#pos <- mmAAshort$hexPost[1]

hexSeq <- function(seq, pos) {
seq <- strsplit(seq, "")
df <- data.frame(seq)
hexSeq <- df[pos[[1]],]
}

mmAAshort <- mutate(mmAAshort, hexSeq=lapply(1:2862, function(x) hexSeq(mmAAshort$seq2[x], mmAAshort$hexPost[x])))

mmAAshort <- mutate(mmAAshort, hexR = str_count(mmAAshort$hexSeq, "R"))

mmAAshort <- filter(mmAAshort, hexR > 0)
dim(mmAAshort)
head(mmAAshort)

mmAAshort <- mutate(mmAAshort, Rhexperc = hexR/helix * 100)
#sum(mouseCPPs$mouseEnsembl %in% mmAAshort$gId)
mmAAtrimmed <- filter(mmAAshort, Rhexperc >= 25)
sum(mouseCPPs$mouseEnsembl %in% mmAAtrimmed$gId)
dim(mmAAtrimmed)


save(mmAAtrimmed, file=file.path(outdir, "filtered10R25percTotalHelix.RData")) #difficult to save beause of lists


#####

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

trialList <- lapply(1:nrow(mmAAshort), function(x) hexSplit(mmAAshort$Prediction[x]))
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

CPPs <- subset(mmAAshort, mmAAshort$tId %in% tIds)
dim(CPPs)
head(CPPs)
sum(CPPs$gId %in% mouseCPPs$mouseEnsembl)


####

#write table to file
mmAA2cut$maxRper8 <- do.call(rbind, mmAA2cut$maxRper8)
#write.table(mmAASS, file.path(datdir, "filteredListPossibleCPPs.txt"), sep = "\t", quote = F)
#mmAA <- read.delim(file.path(datdir, "filteredListPossibleCPPs.txt"))


# ###check enrichment
mouseCPPs <- read.delim(file.path(datdir, "mouseCPPs.txt"))
head(mouseCPPs)
sum(mouseCPPs$mouseEnsembl %in% mmAASP$gId)/nrow(mmAASP)
sum(mouseCPPs$mouseEnsembl %in% mmAA$gId)/nrow(mmAA)
# /nrow(macExpEns)
# sum(mouseCPPs$mouseEnsembl %in% filtMacExp$gId)/nrow(filtMacExp)
# filtMacExp[filtMacExp$gId %in% mouseCPPs$mouseEnsembl,]



