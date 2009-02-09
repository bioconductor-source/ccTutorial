###################################################
### chunk number 1: installPackages eval=FALSE
###################################################
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("Ringo", "biomaRt", "topGO", "ccTutorial"))


###################################################
### chunk number 2: prepare
###################################################
options(length=60, "stringsAsFactors"=FALSE)
set.seed(123)
options(SweaveHooks=list(along=function() par(mar=c(2.5,4.2,4,1.5), font.lab=2),
                         boxplot=function() par(mar=c(5,5,1,1), font.lab=4),
                         dens=function() par(mar=c(4.1, 4.1, 0.1, 0.1), font.lab=2)))


###################################################
### chunk number 3: loadpackage
###################################################
library("Ringo")
library("biomaRt")
library("topGO")
library("xtable")
library("ccTutorial")


###################################################
### chunk number 4: locateData
###################################################
pairDir <- system.file("PairData",package="ccTutorial") 
list.files(pairDir, pattern="pair$")


###################################################
### chunk number 5: exampleFilesTxt
###################################################
read.delim(file.path(pairDir,"files_array1.txt"), header=TRUE)


###################################################
### chunk number 6: readNimblegen
###################################################
RGs <- lapply(sprintf("files_array%d.txt",1:4),
  readNimblegen, "spottypes.txt", path=pairDir)


###################################################
### chunk number 7: showRG
###################################################
head(RGs[[1]]$R)
head(RGs[[1]]$G)
tail(RGs[[1]]$genes)


###################################################
### chunk number 8: showProbeStatus
###################################################
table(RGs[[1]]$genes$Status)


###################################################
### chunk number 9: imageRG
###################################################
RG1breaks <- c(0,quantile(RGs[[1]]$G, probs=seq(0,1,by=0.1)),2^16)
jpeg("ccTutorialArrayImages.jpg", quality=100,
     height=1074*1.5, width=768*1.5)
par(mar=c(0.01,0.01,2.2,0.01))
layout(matrix(c(1,2,5,6,3,4,7,8,9,10,13,14,11,12,15,16),
       ncol=4,byrow=TRUE))
for (this.set in 1:4){
  thisRG <- RGs[[this.set]]
  for (this.channel in c("green","red")){
    my.colors <- colorRampPalette(c("black",paste(this.channel,c(4,1),
                                    sep="")))(length(RG1breaks)-1)
    for (arrayno in 1:2){
      image(thisRG, arrayno, channel=this.channel, 
            mybreaks=RG1breaks, mycols=my.colors, frame.plot=FALSE)
      mtext(side=3, line=0.2, font=2, text=colnames(thisRG[[toupper(
            substr(this.channel,1,1))]])[arrayno])
}}}
dev.off()


###################################################
### chunk number 10: corPlotRG2G
###################################################
corPlot(log2(RGs[[2]]$G))


###################################################
### chunk number 11: corPlotRG2R
###################################################
corPlot(log2(RGs[[2]]$R))


###################################################
### chunk number 12: posToProbeAnno
###################################################
probeAnno <- posToProbeAnno(file.path(system.file("exonerateData",
  package="ccTutorial"), "allChromExonerateOut.txt"))
allChrs <- chromosomeNames(probeAnno)


###################################################
### chunk number 13: posToProbeAnnoExtent
###################################################
genome(probeAnno) <- "M. musculus (mm9)"
arrayName(probeAnno) <- "2005-06-17_Ren_MM5Tiling"


###################################################
### chunk number 14: showProbeAnno
###################################################
show(probeAnno)
ls(probeAnno)


###################################################
### chunk number 15: showProbeAnno2
###################################################
table(probeAnno["9.unique"])


###################################################
### chunk number 16: reporterSpacing
###################################################
startDiffByChr <- lapply(as.list(allChrs), function(chr){
  chrsta <- probeAnno[paste(chr,"start",sep=".")]
  chruni <- probeAnno[paste(chr,"unique",sep=".")]
  ## get start positions of unique reporter match positions
  return(diff(sort(chrsta[chruni=="0"])))})
startDiff <- unlist(startDiffByChr, use.names=FALSE)
table(cut(startDiff, breaks=c(0,50,99,100,200,1000,max(startDiff))))


###################################################
### chunk number 17: makeGffWithBiomaRt
###################################################
ensembl  <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
gene.ids <- unique(unlist(lapply(as.list(c(1:19,"X","Y")), 
   function(this.chr) 
    getBM(attributes="ensembl_gene_id", filters="chromosome_name", 
          values=this.chr, mart=ensembl)[,1]), use.names=FALSE))
    sel.attributes=c("ensembl_gene_id", "mgi_symbol", "chromosome_name",
          "strand", "start_position","end_position", "description")
mm9genes <- getBM(attributes=sel.attributes, filters="ensembl_gene_id",
                  value=gene.ids, mart=ensembl)


###################################################
### chunk number 18: replaceMm9GenesNames
###################################################
mm9genes$name    <- mm9genes$"ensembl_gene_id"
mm9genes$gene    <- mm9genes$"ensembl_gene_id"
mm9genes$chr     <- mm9genes$"chromosome_name"
mm9genes$symbol  <- mm9genes$"mgi_symbol"
mm9genes$start   <- mm9genes$"start_position"
mm9genes$end     <- mm9genes$"end_position"
mm9genes$feature <- rep("gene",nrow(mm9genes))


###################################################
### chunk number 19: removeDuplicatedSymbols
###################################################
if (any(duplicated(mm9genes$name))){
  dupl <- unique(mm9genes$name[duplicated(mm9genes$name)])
  G <- lapply(as.list(dupl), function(this.gene){
    this.gff <- subset(mm9genes,name == this.gene)
    if (nrow(unique(this.gff[,c("name","chr","start","end",
        "description")]))>1) return(this.gff[1,,drop=FALSE])
    non.zero.gff <- subset(this.gff, nchar(symbol)>0)
    this.other.sym <- NULL
    if (nrow(non.zero.gff)> 0){
      shortest <- which.min(nchar(non.zero.gff$symbol))
      this.new.sym <- non.zero.gff$symbol[shortest]
      if (nrow(non.zero.gff)>1)
        this.other.sym <- paste("Synonyms", 
           paste(non.zero.gff$symbol[-shortest],collapse=","),sep=":")
    } else { this.new.sym <- "" }
    this.gff$symbol[1] <- this.new.sym
    if (!is.null(this.other.sym))
      this.gff$description[1] <- paste(this.gff$description[1],
                                       this.other.sym,sep=";")
    return(this.gff[1,,drop=FALSE])
  })
  mm9genes <- rbind(mm9genes[-which(mm9genes$name %in% dupl),],
                    do.call("rbind",G))
}


###################################################
### chunk number 20: reorderMm9
###################################################
mm9genes <- mm9genes[order(mm9genes$chr, mm9genes$start),
  c("name","chr","strand","start","end","symbol","description","feature")]
rownames(mm9genes) <- NULL


###################################################
### chunk number 21: loadMM9Genes eval=FALSE
###################################################
## data(mm9genes)


###################################################
### chunk number 22: setSeed1
###################################################
set.seed(1)


###################################################
### chunk number 23: showmm9genes
###################################################
mm9genes[sample(seq(nrow(mm9genes)),4), 
  c("name", "chr", "strand", "start", "end", "symbol")]


###################################################
### chunk number 24: getGenesGOAnnotation
###################################################
ensembl  <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
ontoGOs <- lapply(as.list(c("biological_process","cellular_component", 
                  "molecular_function")), function(onto){
  ontoBM <- getBM(mart=ensembl, attributes=c("ensembl_gene_id", 
                   paste("go",onto,"id", sep="_"), 
                   paste("go",onto,"linkage_type", sep="_")), 
                  filters="ensembl_gene_id", value=mm9genes$name)
  names(ontoBM) <- c("ensembl_gene_id","go","evidence_code")
  ontoBM <- subset(ontoBM,!( evidence_code %in% c("","IEA","NAS","ND")))
})
mm9GO <- do.call("rbind", ontoGOs)


###################################################
### chunk number 25: getGenesGOAnnotation eval=FALSE
###################################################
## mm9GO <- getBM(attributes=c("ensembl_gene_id","go","evidence_code"), 
##           filters="ensembl_gene_id", value=mm9genes$name, mart=ensembl)
## mm9GO <- subset(mm9GO, !( evidence_code %in% c("", "IEA","NAS","ND")))


###################################################
### chunk number 26: getGOperGene
###################################################
mm9.gene2GO <- with(mm9GO, split(go, ensembl_gene_id))


###################################################
### chunk number 27: loadGenesGOAnnotation eval=FALSE
###################################################
## data(mm9.gene2GO)
## data(mm9.g2p)


###################################################
### chunk number 28: mappingGenesToProbes
###################################################
mm9.g2p <- features2Probes(gff=mm9genes, probeAnno=probeAnno)


###################################################
### chunk number 29: showMm9g2p
###################################################
table(cut(listLen(mm9.g2p),breaks=c(-1,0,10,50,100,500,1200)))


###################################################
### chunk number 30: arrayGenes
###################################################
arrayGenes <- names(mm9.g2p)[listLen(mm9.g2p)>=5]
arrayGenesWithGO <- intersect(arrayGenes, names(mm9.gene2GO))


###################################################
### chunk number 31: preprocess
###################################################
MAs <- lapply(RGs, function(thisRG)
  preprocess(thisRG[thisRG$genes$Status=="Probe",], 
             method="nimblegen", returnMAList=TRUE))
MA <- do.call("rbind",MAs)
X  <- asExprSet(MA)
sampleNames(X) <- paste(X$Cy5, X$Tissue, sep=".")


###################################################
### chunk number 32: showX
###################################################
show(X)


###################################################
### chunk number 33: chipAlongChromCrmp1
###################################################
chipAlongChrom(X, chrom="5", xlim=c(37.63e6,37.64e6), probeAnno=probeAnno,
               gff=mm9genes, ylim=c(-3,5), paletteName="Set2")


###################################################
### chunk number 34: smoothing
###################################################
smoothX <- computeRunningMedians(X, probeAnno=probeAnno, 
   modColumn="Tissue", allChr=allChrs, winHalfSize=450, min.probes=5)
sampleNames(smoothX) <- paste(sampleNames(X),"smoothed",sep=".")


###################################################
### chunk number 35: smoothAlongChromCrmp1
###################################################
chipAlongChrom(X, chrom="5", xlim=c(37.63e6,37.64e6), probeAnno=probeAnno,
               gff=mm9genes, paletteName="Set2", ylim=c(-3,5))
chipAlongChrom(smoothX, chrom="5", xlim=c(37.63e6,37.64e6), add=TRUE,
               probeAnno=probeAnno,ilwd=4, paletteName="Dark2")


###################################################
### chunk number 36: getChersXD
###################################################
data("chersX")
chersXD <- as.data.frame(chersX)


###################################################
### chunk number 37: tryCmarrt
###################################################
library("CMARRT")
cmarrtDat <- do.call("rbind", lapply(as.list(allChrs), function(chr){
  areUni <- probeAnno[paste(chr,"unique",sep=".")]==0
  chrDat <- data.frame("chr"=rep(chr, sum(areUni)),
     "start"=probeAnno[paste(chr,"start",sep=".")][areUni],
     "stop"=probeAnno[paste(chr,"end",sep=".")][areUni],
     "logR"=exprs(X)[probeAnno[paste(chr,"index",sep=".")][areUni],1])
}))

cmarrtRes <- cmarrt.ma(cmarrtDat, M=0.5, frag.length=900,
                       window.opt = "fixed.gen.dist")

cmarrtReg <- cmarrt.peak(cmarrtRes, alpha=0.05, method="BY", minrun=4) 
cmarrtRegDf <- lapply(cmarrtReg, as.data.frame)$cmarrt.bound
names(cmarrtRegDf)[1:3] <- c("chr","start","end")


###################################################
### chunk number 38: showCmarrt
###################################################
head(cmarrtRegDf)


###################################################
### chunk number 39: getRingoAndCmarrt
###################################################
ringoChersChr9 <- subset(chersXD, chr=="9" & cellType=="brain")
cmarrtChersChr9 <- subset(cmarrtRegDf, chr=="9")
dim(ringoChersChr9)
dim(cmarrtChersChr9)


###################################################
### chunk number 40: overlapRingoVsCmarrt
###################################################
chersChr9Overlap <- as.matrix(
   regionOverlap(ringoChersChr9, cmarrtChersChr9))
minRegChr9Len <- outer(with(ringoChersChr9, end-start+1), 
                       with(cmarrtChersChr9, end-start+1), pmin)
fracChr9Overlap <- chersChr9Overlap /minRegChr9Len


###################################################
### chunk number 41: showOverlapRingoVsCmarrt
###################################################
summary(apply(fracChr9Overlap, 1, max))


###################################################
### chunk number 42: loadCherFinding
###################################################
data(chersX)
chersXD <- as.data.frame(chersX)


###################################################
### chunk number 43: getGenesEnrichedPerTissue
###################################################
brainGenes <- getFeats(chersX[sapply(chersX, cellType)=="brain"])
heartGenes <- getFeats(chersX[sapply(chersX, cellType)=="heart"])
brainOnlyGenes <- setdiff(brainGenes, heartGenes)
heartOnlyGenes <- setdiff(heartGenes, brainGenes)


###################################################
### chunk number 44: useTopGO
###################################################
sigGOTable <- function(selGenes, GOgenes=arrayGenesWithGO, 
                       gene2GO=mm9.gene2GO[arrayGenesWithGO],
                       ontology="BP", maxP=0.001)
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
heartRes <- sigGOTable(heartOnlyGenes)


###################################################
### chunk number 45: computeRegionsOverlap
###################################################
brainRegions <- subset(chersXD, cellType=="brain")
heartRegions <- subset(chersXD, cellType=="heart")
chersOBL <- as.matrix(regionOverlap(brainRegions, heartRegions))
minRegLen <- outer(with(brainRegions, end-start+1),
                   with(heartRegions, end-start+1), pmin)
fracOverlap <- chersOBL/minRegLen


###################################################
### chunk number 46: compTissueSpecificRegions
###################################################
brainSpecReg <-  brainRegions[rowMax(fracOverlap)<0.7,]
heartSpecReg <-  heartRegions[rowMax(t(fracOverlap))<0.7,]
mean(is.element(unlist(strsplit(brainSpecReg$features, 
   split="[[:space:]]"), use.names=FALSE), brainOnlyGenes))
selGenes <- intersect(unlist(strsplit(brainSpecReg$features, 
   split="[[:space:]]"), use.names=FALSE), heartGenes)


###################################################
### chunk number 47: targetPos
###################################################
targetPos <- seq(-5000, 10000, by=250)


###################################################
### chunk number 48: compareH3K4me3InSepRegionGenes
###################################################
selGenesVal <- vector("list", length(selGenes))
for (i in 1:length(selGenes)){
  if (i %% 1000 == 0) cat(i,"... ")
  ix <- mm9.g2p[[selGenes[i]]]
  if (length(ix)==0) next
  iy <- exprs(smoothX)[names(ix),,drop=FALSE]
  if (sum(!is.na(iy[,1]))<2) next
  iy <- apply(iy, 2, function(a)
              approx(x=ix, y=a, xout=targetPos)$y)
  selGenesVal[[i]] <- iy
}
names(selGenesVal) <- selGenes


###################################################
### chunk number 49: loadSelGenesVal eval=FALSE
###################################################
## load(file="selGenesVal.RData")


###################################################
### chunk number 50: getDataQuantilesOfSelGenesVal
###################################################
quants <- c(0.1,0.3,0.5,0.7, 0.9)
valsBySample <- lapply(list("H3K4me3.brain.smoothed",
   "H3K4me3.heart.smoothed"), function(thisSample){
      do.call("cbind", lapply(selGenesVal, function(theseValues){
      return(theseValues[,thisSample,drop=FALSE])}))})
names(valsBySample) <- c("H3K4me3.brain.smoothed",
                         "H3K4me3.heart.smoothed")
quantsBySample <- lapply(valsBySample, function(theseValues)
  apply(theseValues, 1, quantile, probs=quants, na.rm=TRUE))


###################################################
### chunk number 51: getProbePosDensity
###################################################
probePosDens <- density(unlist(mm9.g2p, use.names=FALSE))
probePosY <- approx(x=probePosDens$x, y=probePosDens$y, 
                    xout=targetPos)$y
probePosY <- probePosY/max(probePosY)


###################################################
### chunk number 52: sepRegGenesSmoothedQuantiles
###################################################
plot(x=0,y=0, ylim=c(-0.2,3.5), xlim=range(targetPos),
     xaxt="n", type="n", xlab="Distance to gene start [bp]",
     ylab="smoothed fold-change")
axis(side=1, at=seq(-5000, 10000, by=1000),
     labels=paste(-5:10,"kb",sep=""))
mycolors <- c("green","orange")
for (i in 1:2){
  lines(x=targetPos,y=quantsBySample[[i]][3,]*probePosY, 
        lwd=2, col=mycolors[i])
  lines(x=targetPos,y=quantsBySample[[i]][5,]*probePosY, 
        lwd=2, lty=2, col=mycolors[i])
}
legend(x="topleft", fill=c("green","orange"), bty="n", 
       legend=c("H3K4me3 brain", "H3K4me3 heart"))
legend(x="topright", lty=c(2,1), lwd=2, bty="n",
       legend=c("90% quantile", "median"))


###################################################
### chunk number 53: topGOsepRegGenes
###################################################
sepRegRes <- sigGOTable(selGenes=selGenes)
print(sepRegRes)


###################################################
### chunk number 54: printSepRegRes
###################################################
print(xtable(sepRegRes, label="tab-sepRegResGO", caption="\\sl GO terms that are significantly over-represented among genes that show different H3K4me3 regions in heart and brain cells"), type="latex", table.placement="htb", size="scriptsize", include.rownames=FALSE)


###################################################
### chunk number 55: processBarerraExpressionData
###################################################
library("affy")
library("mouse4302cdf")
AB <- ReadAffy(celfile.path=system.file("expression", 
                                        package="ccTutorial"))
barreraExpressionX <- mas5(AB)
barreraExpressionX$Tissue <- sapply(
   strsplit(sampleNames(barreraExpressionX),split="\\."),"[",3)


###################################################
### chunk number 56: loadExpressionData eval=FALSE
###################################################
## ## to speed up creation of the vignette, we load the saved result of the above code chunks
## data("barreraExpressionX")


###################################################
### chunk number 57: mapEnsToAffy
###################################################
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
bmRes <- getBM(attributes=c("ensembl_gene_id","affy_mouse430_2"), 
               filters="ensembl_gene_id", value=arrayGenes, 
               mart=ensembl)
bmRes <- subset(bmRes, nchar(affy_mouse430_2)>0)
arrayGenesToProbeSets <- split(bmRes[["affy_mouse430_2"]], 
                               bmRes[["ensembl_gene_id"]])


###################################################
### chunk number 58: loadArrayGenes eval=FALSE
###################################################
## data("arrayGenesToProbeSets")


###################################################
### chunk number 59: lookArrayGenesToProbeSets
###################################################
table(listLen(arrayGenesToProbeSets))


###################################################
### chunk number 60: sessionInfo
###################################################
toLatex(sessionInfo())


