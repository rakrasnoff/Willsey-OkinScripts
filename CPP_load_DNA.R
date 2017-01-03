#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#browseVignettes("Biostrings")
require(Biostrings)
options(stringsAsFactors = FALSE)
library(stringr)
library(dplyr)
library(zoo)

######Set directories
datdir <- "~/Dropbox/WillseyLab/CPPs/mouse_genome"
outdir <- "~/Dropbox/WillseyLab/CPPs"
#######

######Load date using the biostrings function readDNAStringSet
s = readDNAStringSet(file.path(datdir, "Mus_musculus.GRCm38.cds.all.fa"))
######

######Translate nucleotides to amino acids
translated <- translate(s, genetic.code=GENETIC_CODE, if.fuzzy.codon="solve")
translated[5]
AA <- translated #[c(1:10)] #subset optional
namesAA <- paste(names(AA))
length(namesAA)
namesAA <- as.data.frame(namesAA)
class(namesAA[,1])
head(namesAA)
tId <- do.call("rbind", strsplit(namesAA[,1], "[. :]"))[,1]
gId <- do.call("rbind", strsplit(namesAA[,1], "[. :]"))[,11]
length(tId)
length(gId)
seq = as.data.frame(paste(AA))
dim(seq)
AA2 <- data.frame(tId, gId, seq)
names(AA2) <- c("tId", "gId", "seq")
head(AA2)
#########

######Write results to file
write.table(AA2, "~/Dropbox/WillseyLab/CPPs/MmAAtable.txt", sep = "\t")
#######


# #### Adding R content 
# AAr <- mutate(AA2, length = str_count(AA2[,3]))
# AAr <- mutate(AAr, Rc = str_count(AA2[,3], "R"))
# AAr <- mutate(AAr, Rp = Rc/length * 100)
# head(AAr)
# 
# #####Write results to file
# write.table(AAr, "~/Dropbox/WillseyLab/CPPs/MmAArtable.txt", sep = "\t")
# 
# 
# ###Writing Fasta
# AA4Fasta <- AA2[,c(2,3)]
# AA4Fasta$seq <- gsub("[*].*$","",AA4Fasta$seq)
# head(AA4Fasta)
# AA4Fasta <- AA4Fasta[!duplicated(AA4Fasta[,1]),]
# fasTestG <- dataframe2fas(AA4Fasta[c(1:2000),], file = "~/Dropbox/WillseyLab/CPPs/allAAsGids.fa")
# 
# 
# AA4Fasta <- AA2[,c(1,3)]
# AA4Fasta$seq <- gsub("[*].*$","",AA4Fasta$seq)
# head(AA4Fasta)
# AA4Fasta <- AA4Fasta[,c(1:2)]
# AA4Fasta <- AA4Fasta[!duplicated(AA4Fasta[,1]),]
# fasTestT <- dataframe2fas(AA4Fasta[c(1:2000),], file = "~/Dropbox/WillseyLab/CPPs/allAAsTids.fa")
# #########################
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #########################Begin filtering########################
# AAbyRp <- filter(AA2, Rp >= 20) #list number 1
# dim(AAbyRp)
# AAbyRp
# #seqs for secondary structure prediction
# write.table(AAbyRp[,c(2:3)], file = file.path(datdir, "AAbyRpSeqs.txt"), sep = "\t", quote = F)
# 
# AAbyRp <- arrange(AAbyRp, desc(length))
# AAbyRp$length
# AAbyRp$gId
# 
# write.table(AAbyRp, file = file.path(datdir, "~/Desktop/possCPPsbyRper.txt"), sep = "\t", quote = F)
# 
# ###########
# AAbyRc <- filter(AA2, Rc >=5) #selecting only those with at least 5 Rs, but there are too many, so need to filter further
# dim(AAbyRc)
# 
# countRs <- function(aa) {
# aastring <- AAString(aa)
# nsR <- letterFrequencyInSlidingView(aastring, 15, "R") }
# 
# # nsRs <- lapply(AAbyRp$seq[c(1:10)], countRs)
# # maxes <- lapply(nsRs, max)
# # maxesdf <- as.data.frame(maxes)
# 
# #AAsample <- AAbyRp[c(1:10),]
# AAbyRcShort <- filter(AAbyRc, length < 15) ##List number 2: those under 30 with at least 5 Rs
# AAbyRc15 <- filter(AAbyRc, length >= 15)
# dim(AAbyRc15)
# AAbyRc15 <- mutate(AAbyRc15, maxper15 = lapply(lapply(AAbyRc15[,3], countRs), max))
# head(AAbyRc15)
# AAbyRc15Trim <- filter(AAbyRc15, maxper15 >= 4) #those with lengths of 30 where at least 20% are Rs
# dim(AAbyRc15Trim)
# 
# #################### Create distribution of R counts per window of 30
# AA15 <- filter(AA2, length >= 15) #can only use those of length 30 or longer
# dim(AA15)
# Rcontent <- lapply(AA15$seq, countRs)
# summary(Rcontent[1])
# Rcontent2 <- unlist(Rcontent)
# summary(Rcontent2)
# plot(Rcontent2)
# ####################
# 
# #Check for enrichment of known CPPs in AAbyRc30Trim versus my unfiltered list of >=30 length
# mouseCPPs <- read.delim("~/Documents/Willsey-OkinProject/mouseCPPs.txt")
# head(mouseCPPs)
# 
# sum(mouseCPPs$mouseEnsembl %in% AA15$gId)/nrow(AA15)
# sum(mouseCPPs$mouseEnsembl %in% AAbyRc15Trim$gId)/nrow(AAbyRc15Trim)
# 
# ####################
# #Look at expression in macrophages
# datdir <-"~/Documents/Willsey-OkinProject/MacExpDat"
# wrkdir <-"~/Documents/Willsey-OkinProject/MacExpDat"
# outdir <-"~/Documents/Willsey-OkinProject/MacExpDat"
# 
# #first, need to translate macrophage expression data to ensembl IDs
# macExpDat <- read.delim(file.path(outdir, "macExpDat.txt"))
# macExpDat$Entrez_Gene_ID <- as.character(macExpDat$Entrez_Gene_ID)
# 
# mkeys <- (keys(org.Mm.eg.db, keytype="ENTREZID"))
# cols <- c("ENSEMBL")
# key <- select(org.Mm.eg.db, keys=mkeys, columns = cols, keytype="ENTREZID")
# head(key)
# dim(key)
# 
# macExpEns <- left_join(macExpDat, key, by=c("Entrez_Gene_ID" = "ENTREZID"))
# head(macExpEns)
# tail(macExpEns)
# 
# #now check expression
# macExpEns <- macExpEns[,c(2:14,16)]
# macExpList <- list(macExpEns[,c(4:13)])
# summary(unlist(macExpList)) #getting quartiles
# 
# ##now just looking at maxes
# maxMacExp <- apply(macExpEns[,c(4:13)], 1, max)
# summary(maxMacExp)
# ###
# 
# checkMacExp <- left_join(AAbyRc30Trim, macExpEns, by=c("gId" = "ENSEMBL"))
# head(checkMacExp)
# dim(checkMacExp)
# checkMacExp <- checkMacExp[!duplicated(checkMacExp$seq),]
# head(checkMacExp)
# dim(checkMacExp)
# DF <- data.frame(matrix(unlist(checkMacExp), nrow=nrow(checkMacExp)),stringsAsFactors=FALSE)
# names(DF) <- names(checkMacExp)
# head(DF)
# checkMacExp <- DF
# checkMacExp <- checkMacExp[complete.cases(checkMacExp),]
# dim(checkMacExp)
# maxMacExp <- apply(checkMacExp[,c(10:20)], 1, max)
# filtMacExp <- filter(checkMacExp, maxMacExp >= 122.1)
# dim(filtMacExp)
# 
# ###check enrichment
# sum(mouseCPPs$mouseEnsembl %in% checkMacExp$gId)
# /nrow(macExpEns)
# sum(mouseCPPs$mouseEnsembl %in% filtMacExp$gId)/nrow(filtMacExp)
# filtMacExp[filtMacExp$gId %in% mouseCPPs$mouseEnsembl,]
