#Filter dataframe and check enrichment
datdir <- "~/Dropbox/WillseyLab/CPPs"
outdir <- "~/Dropbox/WillseyLab/CPPs/output"

mmAA <- read.delim(file.path(datdir, "fullAADataFrame.txt"))
head(mmAA)

#Filter by R content (sliding window)

mmAA$maxRper15 <- as.list(mmAA$maxRper15)
head(mmAA)
mmAA2 <- filter(mmAA, maxRper15 > 0)
mmAA2$length <- as.numeric(mmAA2$length)
dim(mmAA2)
head(mmAA2)
mmAA2$keep <- 0
mmAA2$keep[mmAA2$length < 15] <- 1
mmAA2$keep[mmAA2$length >= 15 & mmAA2$maxRper15 >= 4] <- 1
mmAA3 <- filter(mmAA2, keep == 1)
dim(mmAA3)
head(mmAA3)

#Filter by Mac Exp
mmAA4 <- mmAA
mmAA4$keep <- 0
mmAA4$keep[is.na(mmAA$macExp)] <- 1
mmAA4$keep[!is.na(mmAA$macExp) & mmAA$macExp >= 122.5] <- 1
sum(mmAA4$keep)


# ###check enrichment
# sum(mouseCPPs$mouseEnsembl %in% checkMacExp$gId)
# /nrow(macExpEns)
# sum(mouseCPPs$mouseEnsembl %in% filtMacExp$gId)/nrow(filtMacExp)
# filtMacExp[filtMacExp$gId %in% mouseCPPs$mouseEnsembl,]