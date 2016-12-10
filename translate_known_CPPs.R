M2H <- read.delim("~/Documents/General/Mouse_to_human_gene_conversionKey.txt")
known <- read.delim("~/Documents/Willsey-OkinProject/humanMouseNaturalCPPs.txt")
head(M2H)
head(known)
names(known)
toTranslate <- known$GeneName[known$Species == "human"]
toTranslate <- as.data.frame(toTranslate)
names(toTranslate) <- c("humanGeneName")

joined <- left_join(toTranslate, M2H)
head(joined)
dim(joined)

known_mouse <- joined$mouseGeneName
known_mouse_ens <- joined$mouseEnsembl
