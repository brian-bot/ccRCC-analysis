options(stringsAsFactors=F)

require(synapseClient)
require(Biobase)

## GRAB THE NORMALIZED EXPRESSION DATA
normEnt <- loadEntity('syn424079')
normMat <- exprs(normEnt$objects$eset)
clinAll <- normEnt$objects$clinAll

## LOOK AT ONLY THOSE PATIENTS THAT HAVE MATCH TUMOR/NORMAL
idx <- names(table(clinAll$patient)[table(clinAll$patient) == 2])
idx <- clinAll$patient %in% idx
clinAll <- clinAll[idx, ]
normMat <- normMat[, idx]

## GET PLATFORM ANNOTATIONS
qry <- synapseQuery('select id, name from entity where entity.parentId=="syn308579"')
ent <- loadEntity(qry$entity.id[ qry$entity.name == annotation(normEnt$objects$eset) ])
annotations <-  ent$objects$hgu133plus2_annotations
rownames(annotations) <- annotations[, "Probe Set"]
rownames(normMat) <- annotations[rownames(normMat), "Symbol"]

#####
## TAKE TUMOR/NORMAL RATIOS (DIFFERENCES)
#####
normRatios <- sapply(as.list(unique(clinAll$patient)), function(x){
  idt <- which(clinAll$patient == x & clinAll$tumor == "yes")
  idn <- which(clinAll$patient == x & clinAll$tumor == "no")
  normMat[, idt] - normMat[, idn]
})
colnames(normRatios) <- as.character(unique(clinAll$patient))


#####
## LOOK AT TUMOR/NORMAL TEST (IF DIFFERENCE IS SIGNIFICANTLY DIFFERENT FROM ZERO)
#####
resTumorNormal <- apply(normRatios, 1, function(x){
  fit1 <- t.test(x)
  c(fit1$estimate, fit1$p.value)
})
resTumorNormal <- t(resTumorNormal)
colnames(resTumorNormal) <- c("meanTumorNormal", "pValue")
resTumorNormal <- resTumorNormal[order(resTumorNormal[, "pValue"]), ]

hist(resTumorNormal[, "pValue"])
hist(p.adjust(resTumorNormal[, "pValue"], method="BH"))
plot(resTumorNormal[,1], -1*log10(resTumorNormal[, 2]))
abline(h=-1*log10(0.01))
abline(v=c(log2(3/2), log2(2/3)))


## LOOK AT DIFFERENCE BETWEEN SIGN SCORE CATEGORIES (DIFFERENCE OF DIFFERENCES)
resSigncat <- sapply(as.list(rownames(normRatios)), function(x){
  tmp <- unique(clinAll[, c("patient", "signcat")])
  tmp$exprRatio <- as.numeric(normRatios[x, as.character(tmp$patient)])
  fit1 <- t.test(exprRatio ~ signcat, data=tmp)
  c(fit1$estimate, fit1$p.value)
})
resSigncat <- t(resSigncat)
colnames(resSigncat) <- c("meanLowSigncat", "meanHighSigncat", "pValue")
rownames(resSigncat) <- rownames(normMat)
resSigncat <- resSigncat[order(resSigncat[, "pValue"]), ]

hist(resSigncat[, "pValue"])
hist(p.adjust(resSigncat[, "pValue"], method="BH"))
plot(resSigncat[,2]-resSigncat[,1], -1*log10(resSigncat[, 3]))
abline(h=-1*log10(0.01))
abline(v=c(log2(3/2), log2(2/3)))


## LOOK AT DIFFERENCE BETWEEN SIGN SCORE CATEGORIES (DIFFERENCE OF DIFFERENCES)
resNecrosis <- sapply(as.list(rownames(normRatios)), function(x){
  tmp <- unique(clinAll[, c("patient", "necrosis")])
  tmp$exprRatio <- as.numeric(normRatios[x, as.character(tmp$patient)])
  fit1 <- t.test(exprRatio ~ necrosis, data=tmp)
  c(fit1$estimate, fit1$p.value)
})
resNecrosis <- t(resNecrosis)
colnames(resNecrosis) <- c("meanNoNecrosis", "meanYesNecrosis", "pValue")
rownames(resNecrosis) <- rownames(normMat)
resNecrosis <- resNecrosis[order(resNecrosis[, "pValue"]), ]


hist(resNecrosis[, "pValue"])
hist(p.adjust(resNecrosis[, "pValue"], method="BH"))
plot(resNecrosis[,2]-resNecrosis[,1], -1*log10(resNecrosis[, 3]))
abline(h=-1*log10(0.01))
abline(v=c(log2(3/2), log2(2/3)))

