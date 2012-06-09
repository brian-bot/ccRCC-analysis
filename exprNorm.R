## PROGRAM TO PULL IN RAW CEL FILES FROM SYNAPSE AND PERFORM AN UNSUPERVISED NORMALIZATION
#####

options(stringsAsFactors=F)

require(affy)
require(frma)
require(synapseClient)
# synapseLogin() # if not specified in .Rprofile


## DOWNLOAD THE SMOKING CELFILES
smokeEnt <- downloadEntity("syn316208")
smokeCels <- file.path(smokeEnt$cacheDir, smokeEnt$files)

## DOWNLOAD THE OBESITY CELFILES
obeseEnt <- downloadEntity("syn316372")
obeseCels <- file.path(obeseEnt$cacheDir, obeseEnt$files)

## KICK OUT FILE THAT SHOULD NOT HAVE GOTTEN TRANSFERED
allCels <- c(smokeCels, obeseCels)
allCels <- allCels[-grep("design.basenorm", allCels, fixed=T)]

#####
## JUST A DEMONSTRATION OF WHAT WE CAN DO
#####
## WE REALLY NEED THE TUMOR/NORMAL DISTINCTION AND PATIENT IDS (KI NUMBERS, NOT MAYO IDS) TO DO THIS CORRECTLY
## THERE ARE A FEW PATIENTS THAT HAD SAMPLES RUN IN BOTH COHORTS
rawExprs <- ReadAffy(filenames=allCels)
frmaExprs <- frma(rawExprs, background="none", normalize="quantile")

frmaEnt <- createEntity(ExpressionData(list(name="frma normalized data (combined)", parentId = propertyValue(smokeEnt, "parentId"))))
frmaEnt <- addObject(frmaEnt, frmaExprs)
frmaEnt <- storeEntity(frmaEnt)
## frmaEnt ID = syn317285
