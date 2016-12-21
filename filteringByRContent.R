#Filtering by R content

#########################Begin filtering########################
AAbyRp <- filter(AA2, Rp >= 20) #list number 1
dim(AAbyRp)
AAbyRp
#seqs for secondary structure prediction
write.table(AAbyRp[,c(2:3)], file = file.path(datdir, "AAbyRpSeqs.txt"), sep = "\t", quote = F)

AAbyRp <- arrange(AAbyRp, desc(length))
AAbyRp$length
AAbyRp$gId

write.table(AAbyRp, file = file.path(datdir, "~/Desktop/possCPPsbyRper.txt"), sep = "\t", quote = F)

###########
AAbyRc <- filter(AA2, Rc >=5) #selecting only those with at least 5 Rs, but there are too many, so need to filter further
dim(AAbyRc)

countRs <- function(aa) {
  aastring <- AAString(aa)
  nsR <- letterFrequencyInSlidingView(aastring, 15, "R") }

# nsRs <- lapply(AAbyRp$seq[c(1:10)], countRs)
# maxes <- lapply(nsRs, max)
# maxesdf <- as.data.frame(maxes)

#AAsample <- AAbyRp[c(1:10),]
AAbyRcShort <- filter(AAbyRc, length < 15) ##List number 2: those under 30 with at least 5 Rs
AAbyRc15 <- filter(AAbyRc, length >= 15)
dim(AAbyRc15)
AAbyRc15 <- mutate(AAbyRc15, maxper15 = lapply(lapply(AAbyRc15[,3], countRs), max))
head(AAbyRc15)
AAbyRc15Trim <- filter(AAbyRc15, maxper15 >= 4) #those with lengths of 30 where at least 20% are Rs
dim(AAbyRc15Trim)

#################### Create distribution of R counts per window of 30
AA15 <- filter(AA2, length >= 15) #can only use those of length 30 or longer
dim(AA15)
Rcontent <- lapply(AA15$seq, countRs)
summary(Rcontent[1])
Rcontent2 <- unlist(Rcontent)
summary(Rcontent2)
plot(Rcontent2)
####################

#Check for enrichment of known CPPs in AAbyRc30Trim versus my unfiltered list of >=30 length
mouseCPPs <- read.delim("~/Documents/Willsey-OkinProject/mouseCPPs.txt")
head(mouseCPPs)

sum(mouseCPPs$mouseEnsembl %in% AA15$gId)/nrow(AA15)
sum(mouseCPPs$mouseEnsembl %in% AAbyRc15Trim$gId)/nrow(AAbyRc15Trim)


