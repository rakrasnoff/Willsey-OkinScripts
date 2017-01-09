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

# ###check enrichment
mouseCPPs <- read.delim(file.path(datdir, "mouseCPPs.txt"))
head(mouseCPPs)
sum(mouseCPPs$mouseEnsembl %in% mmAASP$gId)/nrow(mmAASP)
sum(mouseCPPs$mouseEnsembl %in% mmAA$gId)/nrow(mmAA)
# /nrow(macExpEns)
# sum(mouseCPPs$mouseEnsembl %in% filtMacExp$gId)/nrow(filtMacExp)
# filtMacExp[filtMacExp$gId %in% mouseCPPs$mouseEnsembl,]