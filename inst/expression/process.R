
library("affy")
library("mouse4302cdf")

AB <- ReadAffy()

barreraExpressionX <- mas5(AB)

barreraExpressionX$Tissue <- sapply(strsplit(sampleNames(barreraExpressionX),split="\\."),"[",3)

save(barreraExpressionX, file="barreraExpressionX_mas5.RData")

