install.packages("stringr")
library(stringr)
dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_66.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
library(xlsx)

CPPtable <- read.csv("~/Documents/Willsey-OkinProject/humanMouseNaturalCpps.csv", header = TRUE, stringsAsFactors = FALSE)

CPPtable$Rs <- str_count(CPPtable$Sequence, "R")
CPPtable$Plength <- nchar(CPPtable$Sequence)
n <- nrow(CPPtable)
for (i in 1:29) {
  CPPtable$Rpercent[i] <- (CPPtable$Rs[i]/CPPtable$Plength[i])*100
}
head(CPPtable)

meanpercent <- mean(CPPtable$Rpercent)
meanlength <- mean(CPPtable$Plength)
meanRs <- mean(CPPtable$Rs)
SDRs <- sd(CPPtable$Rs)
SDlength <- sd(CPPtable$Plength)
median(CPPtable$Rs)

CPPstats <- data.frame("means" = c(meanRs, meanlength, meanpercent), "SDs" = c(SDRs, SDlength, sd(CPPtable$Rpercent)),  row.names = c("Rs", "Length", "Percent"))
CPPstats$medians <- c(median(CPPtable$Rs), median(CPPtable$Plength), median(CPPtable$Rpercent))
CPPstats$min <- c(min(CPPtable$Rs), min(CPPtable$Plength), min(CPPtable$Rpercent))
CPPstats$max <- c(max(CPPtable$Rs), max(CPPtable$Plength), max(CPPtable$Rpercent))


CPPtable[,c(1,2,7)]

write.table(CPPtable, "humanMouseNaturalCppsArgDat.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(CPPstats, "humanMouseNaturalCppStats.txt", sep = "\t", row.names = FALSE, quote = FALSE)
