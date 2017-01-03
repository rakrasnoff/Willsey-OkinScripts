#Get signal peptide data using SignalP4 server

### Set directories
datdir <- "~/Dropbox/WillseyLab/CPPs"
outdir <- "~/Dropbox/WillseyLab/CPPs/fastas"

### Load data
AA2 <- read.delim(file.path(datdir, "MmAAtable.txt"))
head(AA2)

###Writing Fasta with transcript ids
AA4Fasta <- AA2[,c(1,3)]
AA4Fasta$seq <- gsub("[*].*$","",AA4Fasta$seq)
head(AA4Fasta)
AA4Fasta <- AA4Fasta[!duplicated(AA4Fasta[,1]),]
long <- 9000 < str_count(AA4Fasta$seq)
sum(long)
dim(AA4Fasta)
AA4Fasta <- AA4Fasta[!(6000 < str_count(AA4Fasta$seq)),] 
dim(AA4Fasta)

##
fast1 <- dataframe2fas(AA4Fasta[c(1:2000),], file = "~/Dropbox/WillseyLab/CPPs/AAs1.fa")
fast2 <- dataframe2fas(AA4Fasta[c(2001:4000),], file = "~/Dropbox/WillseyLab/CPPs/AAs2.fa")
fast3 <- dataframe2fas(AA4Fasta[c(4001:6000),], file = "~/Dropbox/WillseyLab/CPPs/AAs3.fa")
fast4 <- dataframe2fas(AA4Fasta[c(6001:8000),], file = "~/Dropbox/WillseyLab/CPPs/AAs4.fa")
fast5 <- dataframe2fas(AA4Fasta[c(8001:10000),], file = "~/Dropbox/WillseyLab/CPPs/AAs5.fa")
fast6 <- dataframe2fas(AA4Fasta[c(10001:12000),], file = "~/Dropbox/WillseyLab/CPPs/AAs6.fa")
fast7 <- dataframe2fas(AA4Fasta[c(12001:14000),], file = "~/Dropbox/WillseyLab/CPPs/AAs7.fa")
fast8 <- dataframe2fas(AA4Fasta[c(14001:16000),], file = "~/Dropbox/WillseyLab/CPPs/AAs8.fa")
fast9 <- dataframe2fas(AA4Fasta[c(16001:18000),], file = "~/Dropbox/WillseyLab/CPPs/AAs9.fa")
fast10 <- dataframe2fas(AA4Fasta[c(18001:20000),], file = "~/Dropbox/WillseyLab/CPPs/AAs10.fa")
fast11 <- dataframe2fas(AA4Fasta[c(20001:22000),], file = "~/Dropbox/WillseyLab/CPPs/AAs11.fa")
fast12 <- dataframe2fas(AA4Fasta[c(22001:24000),], file = "~/Dropbox/WillseyLab/CPPs/AAs12.fa")
fast13 <- dataframe2fas(AA4Fasta[c(24001:26000),], file = "~/Dropbox/WillseyLab/CPPs/AAs13.fa")
fast14 <- dataframe2fas(AA4Fasta[c(26001:28000),], file = "~/Dropbox/WillseyLab/CPPs/AAs14.fa")
fast15 <- dataframe2fas(AA4Fasta[c(28001:30000),], file = "~/Dropbox/WillseyLab/CPPs/AAs15.fa")
fast16 <- dataframe2fas(AA4Fasta[c(30001:32000),], file = "~/Dropbox/WillseyLab/CPPs/AAs16.fa")
fast17 <- dataframe2fas(AA4Fasta[c(32001:34000),], file = "~/Dropbox/WillseyLab/CPPs/AAs17.fa")
fast18 <- dataframe2fas(AA4Fasta[c(34001:36000),], file = "~/Dropbox/WillseyLab/CPPs/AAs18.fa")
fast19 <- dataframe2fas(AA4Fasta[c(36001:38000),], file = "~/Dropbox/WillseyLab/CPPs/AAs19.fa")
fast20 <- dataframe2fas(AA4Fasta[c(38001:40000),], file = "~/Dropbox/WillseyLab/CPPs/AAs20.fa")
fast21 <- dataframe2fas(AA4Fasta[c(40001:42000),], file = "~/Dropbox/WillseyLab/CPPs/AAs21.fa")
fast22 <- dataframe2fas(AA4Fasta[c(42001:44000),], file = "~/Dropbox/WillseyLab/CPPs/AAs22.fa")
fast23 <- dataframe2fas(AA4Fasta[c(44001:46000),], file = "~/Dropbox/WillseyLab/CPPs/AAs23.fa")
fast24 <- dataframe2fas(AA4Fasta[c(46001:48000),], file = "~/Dropbox/WillseyLab/CPPs/AAs24.fa")
fast25 <- dataframe2fas(AA4Fasta[c(48001:50000),], file = "~/Dropbox/WillseyLab/CPPs/AAs25.fa")
fast26 <- dataframe2fas(AA4Fasta[c(50001:52000),], file = "~/Dropbox/WillseyLab/CPPs/AAs26.fa")
fast27 <- dataframe2fas(AA4Fasta[c(52001:54000),], file = "~/Dropbox/WillseyLab/CPPs/AAs27.fa")
fast28 <- dataframe2fas(AA4Fasta[c(54001:56000),], file = "~/Dropbox/WillseyLab/CPPs/AAs28.fa")
fast29 <- dataframe2fas(AA4Fasta[c(56001:58000),], file = "~/Dropbox/WillseyLab/CPPs/AAs29.fa")
fast30 <- dataframe2fas(AA4Fasta[c(58001:60000),], file = "~/Dropbox/WillseyLab/CPPs/AAs30.fa")
fast31 <- dataframe2fas(AA4Fasta[c(60001:61418),], file = "~/Dropbox/WillseyLab/CPPs/AAs31.fa")
##

# ####
# datdir <- "~/Dropbox/WillseyLab/CPPs/sigPepFiles"
# files <- dir(datdir)
# 
# SPs <- lapply(files, function(x) { 
#   SPlist <- read.delim(file.path("~/Dropbox/WillseyLab/CPPs/sigPepFiles", x[i]))
# })
# head(SPs[[1]])
# head(SPs[[2]])
# 
# head(AA2)
# 
# for (i in 2:31){
#   names(SPs[[i]]) <- names(SPs[[1]])
#   SPs[[i]] <- SPs[[i]][,c(1:12)]
# }
# summary(SPs)
# 
# AAsp <- rbind.fill(SPs)
# dim(AAsp)
# AAsp2 <- left_join(AA2, AAsp, by=c("tId"= "name"))
# dim(AAsp2)
# head(AAsp2)
# names(AAsp2) <- c("tId", "gId", "seq", "Cmax", "pos", "Ymax", "pos.1", "Smax", "pos.2", "Smean", "D", "SP", "Dmaxcut", "Networks.used")
# AAbySP <- AAsp2[AAsp2$SP == "Y",]
# dim(AAbySP)
# 
# 
