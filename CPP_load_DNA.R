source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
browseVignettes("Biostrings")
require(Biostrings)
install.packages("seqinr")
options(stringsAsFactors = FALSE)
library(stringr)
library(dplyr)
library(zoo)


s = readDNAStringSet("~/Downloads/Mus_musculus.GRCm38.cdna.all.fa")

translated <- translate(s, genetic.code=GENETIC_CODE, if.fuzzy.codon="solve")

translated[5]
AA <- translated #[c(1:10)] #subset optional

namesAA <- paste(names(AA))
length(namesAA)
namesAA <- as.data.frame(namesAA)
class(namesAA[,1])
tId <- do.call("rbind", strsplit(namesAA[,1], "[. :]"))[,1]
gId <- do.call("rbind", strsplit(namesAA[,1], "[. :]"))[,11]
length(tId)
length(gId)
seq = as.data.frame(paste(AA))
dim(seq)

AA2 <- data.frame(tId, gId, seq)
names(AA2) <- c("tId", "gId", "seq")
head(AA2)
head(testAA)

AA2 <- mutate(AA2, length = str_count(AA2[,3]))
AA2 <- mutate(AA2, Rc = str_count(AA2[,3], "R"))
AA2 <- mutate(AA2, Rp = Rc/length * 100)
head(AA2)

AAbyRp <- filter(AA2, Rp >= 20) #list number 1
dim(AAbyRp)
AAbyRp

AAbyRp <- arrange(AAbyRp, desc(length))
AAbyRp$length
AAbyRp$gId

write.table(AAbyRp, file = "~/Desktop/possCPPs.txt", sep = "\t", quote = F)

######
AAbyRc <- filter(AA2, Rc >=5) #selecting only those with at least 5 Rs, but there are too many, so need to filter further
dim(AAbyRc)

countRs <- function(aa) {
aastring <- AAString(aa)
nsR <- letterFrequencyInSlidingView(aastring, 30, "R") }

# nsRs <- lapply(AAbyRp$seq[c(1:10)], countRs)
# maxes <- lapply(nsRs, max)
# maxesdf <- as.data.frame(maxes)

#AAsample <- AAbyRp[c(1:10),]
AAbyRcShort <- filter(AAbyRc, length < 30) ##List number 2: those under 30 with at least 5 Rs
AAbyRc30 <- filter(AAbyRc, length >= 30)
dim(AAbyRc30)
AAbyRc30 <- mutate(AAbyRc30, maxper30 = lapply(lapply(AAbyRc30[,3], countRs), max))
head(AAbyRc30)
AAbyRc30Trim <- filter(AAbyRc30, maxper30 >= 6) #those with lengths of 30 where at least 20% are Rs
dim(AAbyRc30Trim)


