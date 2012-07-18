## PROGRAM TO PULL IN RAW CEL FILES FROM SYNAPSE AND PERFORM A SUPERVISED NORMALIZATION
#####

options(stringsAsFactors=F)

require(affy)
require(snm)
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


#####
## JOINT NORMALIZATION OF ALL DATASETS
#####
rawExprs <- ReadAffy(filenames=allCels)
exprMat <- exprs(rawExprs)

exprNorm1 <- snm(raw.dat=exprMat, 
                 bio.var=model.matrix(~as.factor(clinAll[colnames(exprMat), "tumor"])),
                 int.var = data.frame(patient=as.factor(clinAll[colnames(exprMat), "patient"])))


