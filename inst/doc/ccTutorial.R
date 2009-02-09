###################################################
### chunk number 1: prepare
###################################################
options(length=55, digits=3)
options(SweaveHooks=list(along=function() par(mar=c(2.5,4.2,4,1.5), font.lab=2),
                         boxplot=function() par(mar=c(5,5,1,1), font.lab=2)))
set.seed(1)


###################################################
### chunk number 2: loadpackage
###################################################
library("Ringo")
library("biomaRt")
library("topGO")
library("xtable")
library("ccTutorial")
library("lattice")


###################################################
### chunk number 3: locateData
###################################################
pairDir <- system.file("PairData",package="ccTutorial") 
list.files(pairDir, pattern="pair$")


###################################################
### chunk number 4: remark1 eval=FALSE
###################################################
## # the following chunk 'readNimblegen' requires at least 1GB of RAM
## # and takes about 10 minutes. If time and memory are issues, you can
## # skip this step, see chunk 'remark2' below.


###################################################
### chunk number 5: readNimblegen
###################################################
RGs <- lapply(sprintf("files_array%d.txt",1:4),
  readNimblegen, "spottypes.txt", path=pairDir)


###################################################
### chunk number 6: loadProbeAnno
###################################################
data("probeAnno")
allChrs <- chromosomeNames(probeAnno)


###################################################
### chunk number 7: showmm9genes
###################################################
data("mm9genes")
mm9genes[sample(nrow(mm9genes), 4), 
   c("name", "chr", "strand", "start", "end", "symbol")]


###################################################
### chunk number 8: loadMm9.gene2GO
###################################################
data("mm9.gene2GO")


###################################################
### chunk number 9: loadGenesGOAnnotation
###################################################
data("mm9.g2p")


###################################################
### chunk number 10: arrayGenes
###################################################
arrayGenes <- names(mm9.g2p)[listLen(mm9.g2p)>=5]
arrayGenesWithGO <- intersect(arrayGenes, names(mm9.gene2GO))


###################################################
### chunk number 11: remark2 eval=FALSE
###################################################
## # the following chunk 'preprocess' requires at least 1GB of RAM
## # and takes about 5 minutes. If time and memory are issues, 
## # instead of running that chunk, you can load the result 'X', the
## # ExpressionSet holding the fold changes after preprocessing, by
## data("X")


###################################################
### chunk number 12: preprocess
###################################################
MAs <- lapply(RGs, function(thisRG)
  preprocess(thisRG[thisRG$genes$Status=="Probe",], 
             method="nimblegen", returnMAList=TRUE))
MA <- do.call(rbind, MAs)
X  <- asExprSet(MA)
sampleNames(X) <- paste(X$Cy5, X$Tissue, sep=".")


###################################################
### chunk number 13: chipAlongChromActc1
###################################################
chipAlongChrom(X, chrom="2", xlim=c(113.8725e6,113.8835e6), ylim=c(-3,5),
               probeAnno=probeAnno, gff=mm9genes, paletteName='Set2')


###################################################
### chunk number 14: smoothing
###################################################
smoothX <- computeRunningMedians(X, probeAnno=probeAnno, 
  modColumn="Tissue", allChr=allChrs, winHalfSize=450, min.probes=5)
sampleNames(smoothX) <- paste(sampleNames(X),"smoothed",sep=".")


###################################################
### chunk number 15: smoothAlongChromActc1
###################################################
chipAlongChrom(X, chrom="2", xlim=c(113.8725e6,113.8835e6), ylim=c(-3,5),
               probeAnno=probeAnno, gff=mm9genes, paletteName='Set2')
chipAlongChrom(smoothX, chrom="2", xlim=c(113.8725e6,113.8835e6), ilwd=4,
               probeAnno=probeAnno, paletteName='Dark2', add=TRUE)


###################################################
### chunk number 16: computeY0
###################################################
y0 <- apply(exprs(smoothX), 2, upperBoundNull, prob=0.99)


###################################################
### chunk number 17: histogramsSmoothed
###################################################
myPanelHistogram <- function(x, ...){
  panel.histogram(x, col=brewer.pal(8,"Dark2")[panel.number()], ...)
  panel.abline(v=y0[panel.number()], col="red")
}

h = histogram( ~ y | z, 
      data = data.frame(
        y = as.vector(exprs(smoothX)), 
        z = rep(X$Tissue, each = nrow(smoothX))), 
      layout = c(1,2), nint = 50, 
      xlab = "smoothed reporter level [log2]",
      panel = myPanelHistogram)

print(h)


###################################################
### chunk number 18: computeY0Echo eval=FALSE
###################################################
## y0 <- apply(exprs(smoothX), 2, upperBoundNull, prob=0.99)


###################################################
### chunk number 19: cherFinding
###################################################
chersX <- findChersOnSmoothed(smoothX, 
   probeAnno = probeAnno, 
   thresholds = y0, 
   allChr = allChrs, 
   distCutOff = 450, 
   minProbesInRow = 5, 
   cellType = X$Tissue)


###################################################
### chunk number 20: relateChers
###################################################
chersX <- relateChers(chersX, mm9genes, upstream=5000)


###################################################
### chunk number 21: loadCherFinding eval=FALSE
###################################################
## # since especially the call to relateChers takes some time, we load the
## ## pre-saved image here:
## data("chersX")


###################################################
### chunk number 22: showChers
###################################################
chersXD <- as.data.frame(chersX)
head(chersXD[
  order(chersXD$maxLevel, decreasing=TRUE), 
  c("chr", "start", "end", "cellType", "features", "maxLevel", "score")])


###################################################
### chunk number 23: plotCher1
###################################################
plot(chersX[[which.max(chersXD$maxLevel)]], smoothX, probeAnno=probeAnno, 
     gff=mm9genes, paletteName="Dark2", ylim=c(-1,6))


###################################################
### chunk number 24: showCellType
###################################################
table(chersXD$cellType)


###################################################
### chunk number 25: getGenesEnrichedPerTissue
###################################################
brainGenes <- getFeats(chersX[sapply(chersX, cellType)=="brain"])
heartGenes <- getFeats(chersX[sapply(chersX, cellType)=="heart"])
brainOnlyGenes <- setdiff(brainGenes, heartGenes)
heartOnlyGenes <- setdiff(heartGenes, brainGenes)


###################################################
### chunk number 26: useTopGO
###################################################
sigGOTable <- function(selGenes, GOgenes=arrayGenesWithGO, 
 gene2GO=mm9.gene2GO[arrayGenesWithGO], ontology="BP", maxP=0.001)
{
  inGenes <- factor(as.integer(GOgenes %in% selGenes))
  names(inGenes) <- GOgenes
  GOdata <- new("topGOdata", ontology=ontology, allGenes=inGenes, 
                annot=annFUN.gene2GO, gene2GO=gene2GO)
  myTestStat <- new("elimCount", testStatistic=GOFisherTest, 
                    name="Fisher test", cutOff=maxP)
  mySigGroups <- getSigGroups(GOdata, myTestStat)
  sTab <- GenTable(GOdata, mySigGroups, topNodes=length(usedGO(GOdata)))
  names(sTab)[length(sTab)] <- "p.value"
  sTab <- subset(sTab, as.numeric(p.value) < maxP)
  sTab$Term <- sapply(mget(sTab$GO.ID, env=GOTERM), Term)
  return(sTab)
}

brainRes <- sigGOTable(brainOnlyGenes)
print(brainRes)


###################################################
### chunk number 27: useTopGOHeart
###################################################
heartRes <- sigGOTable(heartOnlyGenes)
print(heartRes)


###################################################
### chunk number 28: loadExpressionData
###################################################
data("barreraExpressionX")


###################################################
### chunk number 29: loadArrayGenesToProbeSets
###################################################
data("arrayGenesToProbeSets")


###################################################
### chunk number 30: compareChIPAndExpression
###################################################
bX <- exprs(barreraExpressionX)
allH3K4me3Genes  <- union(brainGenes, heartGenes)
allH3K4ProbeSets <- unlist(arrayGenesToProbeSets[allH3K4me3Genes])
noH3K4ProbeSets  <- setdiff(rownames(bX), allH3K4ProbeSets)
brainH3K4ExclProbeSets <- unlist(arrayGenesToProbeSets[brainOnlyGenes])
heartH3K4ExclProbeSets <- unlist(arrayGenesToProbeSets[heartOnlyGenes])

brainIdx <- barreraExpressionX$Tissue=="Brain"

brainExpression <- list(
  H3K4me3BrainNoHeartNo  = bX[noH3K4ProbeSets, brainIdx],
  H3K4me3BrainYes        = bX[allH3K4ProbeSets, brainIdx],
  H3K4me3BrainYesHeartNo = bX[brainH3K4ExclProbeSets, brainIdx],
  H3K4me3BrainNoHeartYes = bX[heartH3K4ExclProbeSets, brainIdx]
)


###################################################
### chunk number 31: H3K4me3VsExpression
###################################################
boxplot(brainExpression, col=c("#666666","#999966","#669966","#996666"), 
        names=NA, varwidth=TRUE, log="y", 
        ylab='gene expression level in brain cells')
mtext(side=1, at=1:length(brainExpression), padj=1, font=2, 
      text=rep("H3K4me3",4), line=1)
mtext(side=1, at=c(0.2, 1:length(brainExpression)), padj=1, font=2, 
      text=c("brain/heart","-/-","+/+","+/-","-/+"), line=2)


###################################################
### chunk number 32: testExpressionGreater
###################################################
with(brainExpression, 
     wilcox.test(H3K4me3BrainYesHeartNo, H3K4me3BrainNoHeartNo, 
                 alternative="greater"))


###################################################
### chunk number 33: sessionInfo
###################################################
toLatex(sessionInfo())


###################################################
### chunk number 34: printMm9Genes
###################################################
print(xtable(mm9genes[sample(nrow(mm9genes), 4), 
   c("name", "chr", "strand", "start", "end", "symbol")],
label="tab-mm9genes",
caption="\\sl An excerpt of the object 'mm9genes'."),
type="latex", table.placement="h!t", size="scriptsize",
include.rownames=FALSE)


###################################################
### chunk number 35: printChersXD
###################################################
print(xtable(head(chersXD[order(chersXD$maxLevel, decreasing=TRUE), 
c("chr", "start", "end", "cellType", "features", "maxLevel", "score")]),
label="tab-chersXD",
caption="\\sl The first six lines of object 'chersXD'."),
type="latex", table.placement="h!t", size="scriptsize",
include.rownames=FALSE)


###################################################
### chunk number 36: printBrainRes
###################################################
## for having prettier tables in the PDF, we use 'xtable' here:
print(xtable(brainRes, label="tab-brainResGO", caption="\\sl GO terms that are significantly over-represented among genes showing H3K4me3 enrichment specifically in brain cells"), type="latex", table.placement="h!t", size="scriptsize", include.rownames=FALSE)


###################################################
### chunk number 37: printHeartRes
###################################################
## for having prettier tables in the PDF, we use 'xtable' here:
print(xtable(heartRes, label="tab-heartResGO", caption="\\sl GO terms that are significantly over-represented among genes showing H3K4me3 enrichment specifically in heart cells"), type="latex", table.placement="h!b", size="scriptsize", include.rownames=FALSE)


