library(seqRFLP)

#Filter dataframe and check enrichment
datdir <- "~/Dropbox/WillseyLab/CPPs"
outdir <- "~/Dropbox/WillseyLab/CPPs/output"

mmAA <- read.delim(file.path(datdir, "fullAADataFrame.txt"))
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

mmAASP$maxRper15 <- do.call(rbind, mmAASP$maxRper15)
write.table(mmAASP, file.path(datdir, "filteredListPossibleCPPs.txt"), sep = "\t", quote = F)
mmAA <- read.delim(file.path(datdir, "filteredListPossibleCPPs.txt"))

# ###check enrichment
mouseCPPs <- read.delim(file.path(datdir, "mouseCPPs.txt"))
head(mouseCPPs)
sum(mouseCPPs$mouseEnsembl %in% mmAASP$gId)/nrow(mmAASP)
sum(mouseCPPs$mouseEnsembl %in% mmAA$gId)/nrow(mmAA)
# /nrow(macExpEns)
# sum(mouseCPPs$mouseEnsembl %in% filtMacExp$gId)/nrow(filtMacExp)
# filtMacExp[filtMacExp$gId %in% mouseCPPs$mouseEnsembl,]

#Write to file for 2ry structure prediction
mmAA <- read.delim(file.path(datdir, "filteredListPossibleCPPs.txt"))
fastaDF <- mmAA[,c(1,3)]
fastaDF$seq <- gsub("[*].*$","",fastaDF$seq)
fast1 <- dataframe2fas(fastaDF[c(1:2938),], file = "~/Dropbox/WillseyLab/CPPs/filteredFAs/AAs1.fa")
#fast2 <- dataframe2fas(fastaDF[c(501:1000),], file = "~/Dropbox/WillseyLab/CPPs/filteredFAs/AAs2.fa")
#fast3 <- dataframe2fas(fastaDF[c(1001:1500),], file = "~/Dropbox/WillseyLab/CPPs/filteredFAs/AAs3.fa")
#fast4 <- dataframe2fas(fastaDF[c(1501:2000),], file = "~/Dropbox/WillseyLab/CPPs/filteredFAs/AAs4.fa")
#fast5 <- dataframe2fas(fastaDF[c(2001:2500),], file = "~/Dropbox/WillseyLab/CPPs/filteredFAs/AAs5.fa")
#fast6 <- dataframe2fas(fastaDF[c(2501:2938),], file = "~/Dropbox/WillseyLab/CPPs/filteredFAs/AAs6.fa")


