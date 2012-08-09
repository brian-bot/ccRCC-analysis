## PROGRAM TO PULL IN RAW CEL FILES FROM SYNAPSE AND PERFORM A SUPERVISED NORMALIZATION
#####

options(stringsAsFactors=F)

require(affy)
require(snm)
require(mGenomics)
require(mg.hgu133plus2.db)
require(synapseClient)
# synapseLogin() # if not specified in .Rprofile


## DOWNLOAD THE SMOKING CELFILES
smokeEnt <- downloadEntity("syn316208")
smokeCels <- file.path(smokeEnt$cacheDir, smokeEnt$files)

## DOWNLOAD THE OBESITY CELFILES
obeseEnt <- downloadEntity("syn316372")
obeseCels <- file.path(obeseEnt$cacheDir, obeseEnt$files)

## PULL IN CLINICAL DATA FILES FOR MAPPING
clinEnt <- loadEntity("syn333113")
clinDat <- clinEnt$objects$clin
clinDat <- clinDat[, c("patient", "cel.files", "batch", "tumor") ]

outcomeEnt <- loadEntity("syn333100")
outcomeDat <- outcomeEnt$objects$smob.outc
outcomeDat$batch <- NULL

clinAll <- merge(clinDat, outcomeDat, by="patient", all.x=T, all.y=F)
rownames(clinAll) <- clinAll$cel.files
clinAll$cel.files <- NULL

## KICK OUT FILE THAT SHOULD NOT HAVE GOTTEN TRANSFERED
allCels <- c(smokeCels, obeseCels)
allCels <- allCels[ basename(allCels) %in% rownames(clinAll) ]
clinAll <- clinAll[ basename(allCels), ]
clinAll$fullCelname <- allCels

#####
## JOINT NORMALIZATION OF ALL DATASETS
#####
rawExprs <- ReadAffy(filenames=clinAll$fullCelname)
clinAll$scanDate <- rawExprs@protocolData@data$ScanDate
clinAll$scanDate <- sapply(strsplit(clinAll$scanDate, " "), "[[", 1)
clinAll$scanDate <- as.Date(clinAll$scanDate, format="%m/%d/%y")
thisOrder <- order(clinAll$scanDate)
clinAll <- clinAll[thisOrder, ]
rawExprs <- rawExprs[, thisOrder]
exprMat <- log2(pm(rawExprs))

## OUTLIER PRESENT -- AND POSSIBLY POOR PERFORMING BATCH
png("~/rawBoxplotScanDate.png")
boxplot(rawExprs, col=as.factor(clinAll$scanDate), xaxt="n", xlab="celfiles colored by scan date", ylab = "log2(pm probe expression)")
dev.off()
png("~/rawBoxplotCohort.png")
boxplot(rawExprs, col=as.factor(clinAll$batch), xaxt="n", xlab="celfiles colored by obesity/smoking cohort", ylab = "log2(pm probe expression)")
dev.off()

p75 <- apply(exprMat, 2, quantile, p=0.75)
outlier <- which(p75 == min(p75))
exprMat <- exprMat[, -outlier]
clinAll <- clinAll[-outlier, ]

## LOOK AT INTENSITY DEPENDENT EFFECTS
intFit <- snm(exprMat,
              int.var=data.frame(array=factor(1:ncol(exprMat))))
intSvd <- fs(intFit$norm.dat)
plot(intSvd$d)
xyplot(intSvd$v[,2] ~ intSvd$v[,1], groups=clinAll$scanDate, xlab="eigengene 1", ylab="eigengene 2", main="svd colored by scan date")

png("~/normIntensityonlySvd1.png")
plot(intSvd$v[,1], col=factor(clinAll$scanDate), xaxt="n", xlab="sample ordered by scan date", ylab = "eigengene 1", main="colored by scan date")
dev.off()
png("~/normIntensityonlySvd2.png")
plot(intSvd$v[,2], col=factor(clinAll$tumor), xaxt="n", xlab="sample ordered by scan date", ylab = "eigengene 2", main="colored by tumor/normal")
dev.off()

## ONE BATCH THAT IS OUT TO LUNCH -- SAME AS IN BOXPLOTS


## TRY FITTING A BATCH-SPECIFIC INTENSITY EFFECT
# intBatchFit <- snm(exprMat,
#                    int.var=data.frame(batch=factor(clinAll$scanDate)))
# intBatchSvd <- fs(intBatchFit$norm.dat)
# plot(intBatchSvd$d)
# xyplot(intBatchSvd$v[,2] ~ intBatchSvd$v[,1], groups=clinAll$scanDate=="2006-07-30")
## EVEN FITTING INTENSITY DEPENDENT EFFECT WE STILL SEE BATCH EFFECTS -- REMOVE THAT BATCH (OLDEST)

outliers <- which(clinAll$scanDate == "2006-07-30")
exprMat <- exprMat[, -outliers]
clinAll <- clinAll[-outliers, ]

## FIT JUST INTENSITY DEPENDENT WITH OUTLIERS REMOVED AND CHECK SIGNAL
intNewFit <- snm(exprMat,
                 int.var=data.frame(array=factor(1:ncol(exprMat))))
intNewSvd <- fs(intNewFit$norm.dat)
plot(intNewSvd$d)
xyplot(intNewSvd$v[,2] ~ intNewSvd$v[,1], groups=clinAll$scanDate)

png("~/normIntensityFilterSvd1.png")
plot(intNewSvd$v[,1], col=factor(clinAll$tumor), xaxt="n", xlab="sample ordered by scan date", ylab = "eigengene 1", main="colored by tumor/normal")
dev.off()
png("~/normIntensityFilterSvd2.png")
plot(intNewSvd$v[,2], col=factor(clinAll$batch), xaxt="n", xlab="sample ordered by scan date", ylab = "eigengene 2", main="colored by obesity/smoking cohort")
dev.off()
# s <- sample(nrow(exprMat), 20000)
# plot(rowMeans(exprMat[s,]), intNewSvd$u[s,2])
## THE OTHER TWO SCAN DATES FROM OBESITY COHORT STICKING OUT FROM SMOKING COHORT

## FIT WITH ADJUSTMENT VARIABLE FOR BATCH AND INCLUDE TUMOR AND PATIENT AS BIOLOGICAL VARIABLES
normFit <- snm(raw.dat = exprMat,
             adj.var = model.matrix(~as.factor(clinAll$batch)),
             bio.var = model.matrix(~as.factor(clinAll$tumor)),
             int.var = data.frame(int=as.factor(1:nrow(clinAll))),
             rm.adj = TRUE,
             diagnose = FALSE)
evalThis <- fs(normFit$norm.dat)
plot(evalThis$d)
xyplot(evalThis$v[,2] ~ evalThis$v[,1], groups=clinAll$scanDate)

png("~/normFullmodelSvd1.png")
plot(evalThis$v[,1], col=factor(clinAll$tumor), xaxt="n", xlab="sample ordered by scan date", ylab = "eigengene 1", main="colored by tumor/normal")
dev.off()
png("~/normFullmodelSvd2.png")
plot(evalThis$v[,2], col=factor(clinAll$scanDate), xaxt="n", xlab="sample ordered by scan date", ylab = "eigengene 2", main="colored by scan date")
dev.off()
## THIS LOOKS REALLY GOOD -- AT LEAST FOR TUMOR/NORMAL QUESTIONS -- CALLING IT SUFFICIENTLY 'NORMALIZED'
## WOULD LIKELY NORMALIZE SEPARATELY FOR SMOKING OR OBESITY RELATED ANALYSES


## NOW PREPARE DATA TO RUN THROUGH ENTIRE mGenomics AFFY WORKFLOW
myDir <- tempfile()
dir.create(myDir)
for( i in clinAll$fullCelname ){
  system(paste("ln -s ", i, " ", myDir, sep=""), ignore.stdout=T, ignore.stderr=T)
}
## ENSURE ORDERING WITH CLINICAL INFORMATION IS CORRECT
celFiles <- list.celfiles(myDir, full.names = T, recursive = T)
clinAll <- clinAll[basename(celFiles), ]
adjVar <- model.matrix(~as.factor(clinAll$batch))
bioVar <- model.matrix(~as.factor(clinAll$tumor))
intVar <- data.frame(int=as.factor(1:nrow(clinAll)))


finalFit <- affyWorkflow.SNM(rawDataDir=myDir, bio.var=bioVar, adj.var=adjVar, int.var=intVar, rm.adj=TRUE)

## CHECK ONE MORE TIME BEFORE PUSHING UP
eset <- finalFit$hgu133plus2$eset
statistics <- finalFit$hgu133plus2$statistics
thisOrder <- order(clinAll$scanDate)
clinAll <- clinAll[thisOrder, ]
eset <- eset[, thisOrder]

finalSvd <- fs(exprs(eset))
plot(finalSvd$d)
xyplot(finalSvd$v[,2] ~ finalSvd$v[,1], groups=clinAll$scanDate)

png("~/normFinalSummarizedSvd1.png")
plot(finalSvd$v[,1], col=factor(clinAll$tumor), xaxt="n", xlab="sample ordered by scan date", ylab = "eigengene 1", main="colored by tumor/normal")
dev.off()
png("~/normFinalSummarizedSvd2.png")
plot(finalSvd$v[,2], col=factor(clinAll$scanDate), xaxt="n", xlab="sample ordered by scan date", ylab = "eigengene 2", main="colored by scan date")
dev.off()


## OK -- NOW PUSH UP TO SYNAPSE (syn424079)
normEnt <- createEntity(Data(list(name="post normalization objects", 
                                  description="r objects to be used for analysis post normalization and qc including an eset containing matrix of normalized and summarized expression values as per Mecham (Bioinformatics 2010, PMID: 20363728), a statistics object containing singular values and probe weights, and a clinical data object containing sample-specific information.",
                                  parentId=smokeEnt$properties$parentId)))
normEnt <- addObject(normEnt, eset)
normEnt <- addObject(normEnt, statistics)
normEnt <- addObject(normEnt, clinAll)
normEnt <- storeEntity(normEnt)
