source("https://bioconductor.org/biocLite.R")
biocLite("BioSeqClass")
library(BioSeqClass)

datdir <- "~/Documents/Willsey-OkinProject/mouse_genome"
AAbyRp <- read.delim("~/Desktop/possCPPsbyRper.txt")
head(AAbyRp)

test1 <- predictPROTEUS(AAbyRp$seq[1])


test4 <- getDSSP(AAbyRp$seq[1])
test3 <- predictPROTEUS(AAbyRp$seq[1],proteus2.organism="euk")

test2 <- featureSSC(AAbyRp$seq[1])


###tutorial
file = file.path(path.package("BioSeqClass"), "example", "acetylation_K.fasta")  
tmp = readAAStringSet(file) 
proteinSeq = as.character(tmp)
PROTEUS = predictPROTEUS(proteinSeq[1:2],proteus2.organism="euk") 

######## protr demo
install.packages("protr")
require(protr)

extracell = readFASTA(system.file('protseq/extracell.fasta'
                                  , package = 'protr'))
head(extracell)
length(extracell)
extracell <- extracell[(sapply(extracell, protcheck))]
length(extracell)
x1 = t(sapply(extracell, extractAPAAC))
labels = as.factor(rep(0, length(extracell)))


test00 <- readFASTA('~/Desktop/testAAfasta.fa')
