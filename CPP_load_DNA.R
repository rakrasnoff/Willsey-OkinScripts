source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
browseVignettes("Biostrings")

test <- read.delim("~/Downloads/Mus_musculus.GRCm38.cdna.all.fa", sep = ".")
head(test)
names(test)

s = readDNAStringSet("~/Downloads/Mus_musculus.GRCm38.cdna.all.fa")
subseq(s, start=c(1, 2, 3), end=c(3, 6, 5))
consensusMatrix(s, baseOnly=TRUE)

df <- subseq(s)
head(df)

ambiguityMap <- (c("AGTC"))
names(ambiguityMap) <- "N"

testAA <- translate(s, genetic.code=GENETIC_CODE, if.fuzzy.codon="solve")
head(testAA)
