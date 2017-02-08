library(seqRFLP)

#Filter dataframe and check enrichment
datdir <- "~/Dropbox/WillseyLab/CPPs"
outdir <- "~/Dropbox/WillseyLab/CPPs/output"

mmAA <- read.delim(file.path(datdir, "finalAADataframe.txt"))
head(mmAA)

#Filter by R content (sliding window)

mmAA$maxRper15 <- as.list(mmAA$maxRper15)
head(mmAA)
mmAAR <- filter(mmAA, maxRper15 > 0)
mmAAR$length <- as.numeric(mmAAR$length)
dim(mmAAR)
head(mmAAR)
mmAAR$keep <- 0
mmAAR$keep[mmAAR$length < 15] <- 1
mmAAR$keep[mmAAR$length >= 15 & mmAAR$maxRper15 >= 4] <- 1
mmAAR <- filter(mmAAR, keep == 1)
dim(mmAAR)
head(mmAAR)

#Filter by Mac Exp
#mmAAME <- mmAA #to filter original list
mmAAME <- mmAAR # to continue filtering previously filtered list
mmAAME$keep <- 0
mmAAME$keep[is.na(mmAAME$macExp)] <- 1
mmAAME$keep[!is.na(mmAAME$macExp) & mmAAME$macExp >= 122.5] <- 1
sum(mmAAME$keep)
#20684 (keeping those with NAs)

#Filter by Signal Peptides
mmAASP <- mmAAME
mmAASP <- filter(mmAASP, SigPep == "Y")
dim(mmAASP)
#2938



#Filter by Secondary Structure
head(mmAASP)
mmAASS <- mmAASP
mmAASS$helix <- 0
mmAASS <- mutate(mmAASS, helix = str_count(mmAASS$Prediction, "H"))
head(mmAASS)

argHelix <- function(peptide) {
  pred <- strsplit(mmAASS$Prediction[mmAASS$tId == peptide], "")
  conf <- strsplit(mmAASS$Confidence[mmAASS$tId == peptide], "")
  seq <- strsplit(mmAASS$seq[mmAASS$tId == peptide], "")
  argHelixDF <- data.frame(pred, conf, seq)
  names(argHelixDF) <- c("pred", "conf", "seq")
  argHelixDF
}


argHelixTest <- function(peptide) {
  struct <- peptide$pred[peptide$seq == "R"]
  argPerHelix <- sum(str_count(struct, "H"))
}

#for rerunning 
test <- data.frame(mmAASS$tId, mmAASS$seq, mmAASS$length, mmAASS$Prediction)
names(test) <- c("tId", "seq", "length", "Prediction")
#test$seq <- gsub("$[*]", "", test$seq)
test$predLength <- data.frame(str_count(test[,4]))
test <- mutate(test, lengthLessOne = length - 1)
test$eq <- array(0, 2938)
head(test)
test$eq[is.na(test$predLength)] <- 0
test$eq[test$length == test$predLength] <- 0
test$eq[test$length != test$predLength] <- 1
test$eq[test$lengthLessOne == test$predLength] <- 0

test2 <- test[test$eq == 1,]
test2$seq <- gsub("[*]", "", test2$seq)
test2Fas <- test2[,c(1,2)]
fasta2 <- dataframe2fas(test2Fas, )




#write table to file
mmAASS$maxRper15 <- do.call(rbind, mmAASP$maxRper15)
write.table(mmAASS, file.path(datdir, "filteredListPossibleCPPs.txt"), sep = "\t", quote = F)
mmAA <- read.delim(file.path(datdir, "filteredListPossibleCPPs.txt"))



# ###check enrichment
mouseCPPs <- read.delim(file.path(datdir, "mouseCPPs.txt"))
head(mouseCPPs)
sum(mouseCPPs$mouseEnsembl %in% mmAASP$gId)/nrow(mmAASP)
sum(mouseCPPs$mouseEnsembl %in% mmAA$gId)/nrow(mmAA)
# /nrow(macExpEns)
# sum(mouseCPPs$mouseEnsembl %in% filtMacExp$gId)/nrow(filtMacExp)
# filtMacExp[filtMacExp$gId %in% mouseCPPs$mouseEnsembl,]



