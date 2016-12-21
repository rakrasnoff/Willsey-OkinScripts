M2H <- read.delim("~/Documents/General/Mouse_to_human_gene_conversionKey.txt")
known <- read.delim("~/Documents/Willsey-OkinProject/humanMouseNaturalCPPs.txt")
head(M2H)
head(known)
names(known)
toTranslate <- known$GeneName[known$Species == "human"]
toKeep <-  known$GeneName[known$Species == "mouse"]
toKeep <- as.data.frame(toKeep)
names(toKeep) <- c("mouseGeneName")
toTranslate <- as.data.frame(toTranslate)
names(toTranslate) <- c("humanGeneName")

joined <- left_join(toTranslate, M2H)
head(joined)
dim(joined)

known_mouse <- joined$mouseGeneName
known_mouse_ens <- joined$mouseEnsembl

toKeep2 <- left_join(toKeep, M2H)
toKeep2 <- toKeep2[,c(1,3)]

known_mouse_Id_ens <- joined[,c(2:3)]
head(known_mouse_Id_ens)
mouseCPPs <- rbind(known_mouse_Id_ens, toKeep2)
mouseCPPs
mouseCPPs <- mouseCPPs[complete.cases(mouseCPPs),]

write.table(mouseCPPs, "~/Documents/Willsey-OkinProject/mouseCPPs.txt", sep = "\t", quote = F)
