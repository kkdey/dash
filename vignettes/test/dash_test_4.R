
#####################  dash on metagenomics data   ####################

library(metagenomeSeq)
data("lungData")
lungData
pheno <- pData(lungData)
featureData(lungData)
ff <- fData(lungData)
MRcounts(lungData)
dd <- MRcounts(mouseData)
