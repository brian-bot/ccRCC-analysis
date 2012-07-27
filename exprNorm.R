## PROGRAM TO PULL IN RAW CEL FILES FROM SYNAPSE AND PERFORM A SUPERVISED NORMALIZATION
#####

options(stringsAsFactors=F)

require(affy)
require(snm)
require(mGenomics)
#require(preprocessCore)
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
exprMat <- log2(pm(rawExprs))
scanDate <- rawExprs@protocolData@data$ScanDate
clinAll$scanDate <- sapply(strsplit(scanDate, " "), "[[", 1)

## OUTLIER PRESENT -- AND POSSIBLY POOR PERFORMING BATCH
boxplot(rawExprs, col=as.factor(clinAll$scanDate))

p75 <- apply(exprMat, 2, quantile, p=0.75)
outlier <- which(p75 == min(p75))
rawExprs <- rawExprs[, -outlier]
exprMat <- exprMat[, -outlier]
clinAll <- clinAll[-outlier, ]

## LOOK AT INTENSITY DEPENDENT EFFECTS
intFit <- snm(exprMat,
              int.var=data.frame(array=factor(1:ncol(exprMat))))
intSvd <- fs(intFit$norm.dat)
plot(intSvd$d)
xyplot(intSvd$v[,2] ~ intSvd$v[,1], groups=clinAll$scanDate)
## ONE BATCH THAT IS OUT TO LUNCH -- SAME AS IN BOXPLOTS


## TRY FITTING A BATCH-SPECIFIC INTENSITY EFFECT 
intBatchFit <- snm(exprMat,
                   int.var=data.frame(batch=factor(clinAll$scanDate)))
intBatchSvd <- fs(intBatchFit$norm.dat)
plot(intBatchSvd$d)
xyplot(intBatchSvd$v[,2] ~ intBatchSvd$v[,1], groups=clinAll$scanDate=="07/30/06")
## EVEN FITTING INTENSITY DEPENDENT EFFECT WE STILL SEE BATCH EFFECTS -- REMOVE THAT BATCH (OLDEST)

outliers <- which(clinAll$scanDate == "07/30/06")
rawExprs <- rawExprs[, -outliers]
exprMat <- exprMat[, -outliers]
clinAll <- clinAll[-outliers, ]

## FIT JUST INTENSITY DEPENDENT WITH OUTLIERS REMOVED AND CHECK SIGNAL
intNewFit <- snm(exprMat,
                 int.var=data.frame(array=factor(1:ncol(exprMat))))
intNewSvd <- fs(intNewFit$norm.dat)
plot(intNewSvd$d)
xyplot(intNewSvd$v[,2] ~ intNewSvd$v[,1], groups=clinAll$scanDate)
plot(intNewSvd$v[,1], col=factor(clinAll$batch))
plot(intNewSvd$v[,2], col=factor(clinAll$batch))
s <- sample(nrow(exprMat), 20000)
plot(rowMeans(exprMat[s,]), intNewSvd$u[s,2])
## THE OTHER TWO SCAN DATES FROM OBESITY COHORT STICKING OUT FROM SMOKING COHORT

## FIT WITH ADJUSTMENT VARIABLE FOR BATCH AND INCLUDE TUMOR AND PATIENT AS BIOLOGICAL VARIABLES
eval1 <- snm(raw.dat = exprMat,
             adj.var = model.matrix(~as.factor(clinAll$batch)),
             bio.var = model.matrix(~as.factor(clinAll$tumor)),
             int.var = data.frame(int=as.factor(1:nrow(clinAll))),
             rm.adj = TRUE,
             diagnose = FALSE)
evalThis <- fs(eval1$norm.dat)
plot(evalThis$d)
xyplot(evalThis$v[,2] ~ evalThis$v[,1], groups=clinAll$scanDate)
plot(evalThis$v[,1], col=factor(clinAll$batch))
plot(evalThis$v[,2], col=factor(clinAll$batch))
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
normSvd <- fs(exprs(eset))
plot(normSvd$d)
xyplot(normSvd$v[,2] ~ normSvd$v[,1], groups=clinAll$scanDate)
plot(normSvd$v[,1])



## OK -- NOW PUSH UP TO SYNAPSE
normEnt <- createEntity(Data(list(name="snmNormalizedExprs", parentId=smokeEnt$properties$parentId)))
normEnt <- addObject(normEnt, eset)
normEnt <- addObject(normEnt, statistics)
normEnt <- storeEntity(normEnt)